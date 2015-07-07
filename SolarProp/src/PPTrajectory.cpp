#include "PPTrajectory.h"

PPTrajectory::PPTrajectory(double ri, double thi, double phii, double ei, const std::string& paramFileName,
        const OutputFormat& oFormat) : traj(1, PPPoint(ri, thi, phii, ei, 0)), status(BoundaryHit::None),
    b(0, 0, 0), curlBoverB(0, 0, 0), kTensor(0, 0, 0, 0, 0, 0, 0, 0, 0), kpar(0), kperp(0), vdrift(0, 0, 0),
    updatedB(false), updatedVdrift(false), updatedKTensor(false),
	generator(std::chrono::system_clock::now().time_since_epoch().count()), ndistro(0.0, 1.0),
	params(paramFileName), outputFormat(oFormat) {
	// Check whether the particle is starting in the heliosphere...
	if (traj.back().getR() >= params.getRHP()) {
		status = BoundaryHit::Heliopause;
	} else if (traj.back().getR() <= params.getRSun()) {
		// If the particle enters the sun
		status = BoundaryHit::Sun;
	}

	// Update all parameters
	updateB();
	updateKTensor();
	updateVdrift();
}

double PPTrajectory::vel() const {
	return cAUs * sqrt(1 - pow(params.getM() / (traj.back().getE() + params.getM()), 2));
}

// Returns rigidity in T*AU
double PPTrajectory::rigidity() const {
	// Conversion factor: 1 GV/c = 2.23*10^-11 T*AU
	return 2.23 * pow(10, -11) * sqrt(traj.back().getE() * (traj.back().getE() + 2 * params.getM()));
}

double PPTrajectory::lambdaPar() {
	if (rigidity() >= params.getRig0()) {
		return params.getLambda0() * (rigidity() / params.getRig0()) 
			* (1 + traj.back().getR() / params.getR0());
	} else {
		return params.getLambda0() * (1 + traj.back().getR() / params.getR0());
	}

    /*
    // Get field at earth
    updateB();
    double be = params.getB0() * sqrt(1 + pow(params.getOmega() * params.getR0() / params.getVsw(), 2));

    return params.getLambda0() * (rigidity() / params.getRig0()) * (be / b.magnitude());
    */
}

double PPTrajectory::fs() const {
	double wts2 = pow(rigidity() / params.getRig0(), 2);
	return wts2 / (1 + wts2);
}

double PPTrajectory::gamma() const {
	return (traj.back().getE() + 2 * params.getE0()) / (traj.back().getE() + params.getE0());
}

PPTrajectory& PPTrajectory::updateB() {
	if (!updatedB) {
		double r = traj.back().getR();
		double sinth = sin(traj.back().getTh());
		double h = heaviside(traj.back().getTh(), getThp());

		// Compute B
		b.r = params.getAc() * params.getB0() * pow(params.getR0() / r, 2) * h;
		b.phi = -b.r * params.getOmega() * r * sinth / params.getVsw();

		updatedB = true;
	}

	return *this;
}

PPTrajectory& PPTrajectory::updateKTensor() {
	if (!updatedKTensor) {
		// Make sure the magnetic field has been updated!
		updateB();

		kpar = (vel() / 3) * lambdaPar();
		kperp = params.getDiffFact() * kpar;

		// This component is independent of angle
		kTensor.rr = (b.r*b.r * kpar + b.phi*b.phi * kperp) / (b.r*b.r + b.phi*b.phi);
		kTensor.thth = kperp;
		kTensor.phiphi = (b.phi*b.phi * kpar + b.r*b.r * kperp) / (b.r*b.r + b.phi*b.phi);
		kTensor.rphi = b.r * b.phi * (kpar - kperp) / (b.r*b.r + b.phi*b.phi);
        kTensor.phir = kTensor.rphi;

		// These should always be zero!!!
		kTensor.rth = 0;
        kTensor.thr = 0;
		kTensor.phith = 0; 
		kTensor.thphi = 0; 

		updatedKTensor = true;
	}

	return *this;
}

PPTrajectory& PPTrajectory::updateVdrift() {
	if (!updatedVdrift) {
		// Make sure the magnetic field has been updated!
		updateB();

		double r = traj.back().getR();
		double sinth = sin(traj.back().getTh());
		double costh = cos(traj.back().getTh());
		double rig = rigidity();
		// Larmor radius
		double rL = rig / b.magnitude();
		// Convenience parameter
		double gammaGC = r * params.getOmega() * sinth / params.getVsw();

		// Compute v_d^HCS
		// For convenience, initialize the parts of the r and phi components of vdHCS which are the same
		double vdHCSr = 0;
		double vdHCSphi = 0;
		// The formula for HCS drifts only applies for |d| <= 2 rL!
		if (fabs(r * costh) <= 2 * rL) {
			vdHCSr = (0.457 - 0.412 * fabs(r * costh / rL) + 0.0915 * pow(r * costh / rL, 2))
				* params.getAc() * params.getQSign() * vel();

			// Put in the factors that are different between the components
			vdHCSphi = vdHCSr / sqrt(1.0 + gammaGC * gammaGC);
			vdHCSr *= gammaGC / sqrt(1.0 + gammaGC * gammaGC);
		}

		// Compute v_d^gc factor
		double vdgcFact = 2 * rig * vel() * r * params.getQSign() * params.getAc() 
			/ (3 * params.getB0() * pow(params.getR0(), 2) * pow(1 + gammaGC * gammaGC, 2));
		vdgcFact *= heaviside(traj.back().getTh(), getThp());

		// Incorporate the drift reduction factor and combine all components
		vdrift.r = fs() * (vdHCSr + vdgcFact * (-gammaGC * costh / sinth));
		vdrift.th = fs() * (vdgcFact * (2 + gammaGC * gammaGC) * gammaGC);
		vdrift.phi = fs() * (vdHCSphi + vdgcFact * gammaGC * gammaGC * costh / sinth);

		updatedVdrift = true;
	}

	return *this;
}

double PPTrajectory::getDr(const double& dWr, const double& dWphi) const {
    // Convenience variables
	double r = traj.back().getR();
	double gammaGC = r * params.getOmega() * sin(traj.back().getTh()) / params.getVsw();

    // dr/ds
	double dsCdr = (kpar * (3 * params.getR0() + r * (3 + gammaGC * gammaGC)) + kperp * gammaGC*gammaGC
        * (2*params.getR0() * (2 + gammaGC*gammaGC) + r * (5 + 3 * gammaGC*gammaGC)))
        / (r * (params.getR0() + r) * pow(1 + gammaGC*gammaGC, 2))
        - params.getVsw()
        - vdrift.r;

    // dr/dWr
	double dWrCdr = sqrt(2*kTensor.rr - 2*pow(kTensor.rphi, 2) / kTensor.phiphi);

    // dr/dWphi
	double dWphiCdr = sqrt(2 / kTensor.phiphi) * kTensor.rphi;
    
    return dsCdr * params.getDs() + dWrCdr * dWr + dWphiCdr * dWphi;
}

double PPTrajectory::getDth(const double& dWth) const {
    // Convenience variable
	double r = traj.back().getR();

	// dth/ds
	double dsCdth = kperp * cos(traj.back().getTh()) / (r * r * sin(traj.back().getTh())) - vdrift.th / r;

	// dth/dWth
	double dWthCdth = sqrt(2 * kTensor.thth) / r;
    
    return dsCdth * params.getDs() + dWthCdth * dWth;
}

double PPTrajectory::getDphi(const double& dWphi) const {
    // Convenience variables
	double r = traj.back().getR();
	double sinth = sin(traj.back().getTh());
	double gammaGC = r * params.getOmega() * sinth / params.getVsw();

    // dphi/ds
	double dsCdphi = ((kperp - kpar) * params.getOmega() * ((gammaGC*gammaGC * + 3) * r + 2 * params.getR0()))
        / ((params.getR0() + r) * r * params.getVsw() * pow(1 + gammaGC*gammaGC, 2))
        - vdrift.phi / (r * sinth);

    // dphi/dWphi
	double dWphiCdphi = sqrt(2 * kTensor.phiphi) / (r * sinth);
    
    return dsCdphi * params.getDs() + dWphiCdphi * dWphi;
}

double PPTrajectory::getDe() const {
	return params.getDs() * 2 * traj.back().getE() * params.getVsw() * gamma() / (3 * traj.back().getR());
}

const BoundaryHit& PPTrajectory::step() {
	// Update all parameters
	updateB();
	updateKTensor();
	updateVdrift();
	
	// Generate Wiener terms
	double dWr = ndistro(generator) * sqrt(params.getDs());
	double dWth = ndistro(generator) * sqrt(params.getDs()); 
	double dWphi = ndistro(generator) * sqrt(params.getDs());

    // Make particle at the next point
	traj.push_back(PPPoint(traj.back().getR() + getDr(dWr, dWphi),
                traj.back().getTh() + getDth(dWth),
                traj.back().getPhi() + getDphi(dWphi),
                traj.back().getE() + getDe(),
                traj.back().getS() + params.getDs()));
	
    // Check if the particle hit the heliopause or sun
    return checkStatus();
}

const BoundaryHit& PPTrajectory::checkStatus() {
	if (traj.back().getR() >= params.getRHP()) {
		status = BoundaryHit::Heliopause;
        return status;
	} else { // Particle may have gone through the sun.  Check this.
        // Get previous point
        const double& r = traj.end()[-2].getR();
        const double& th = traj.end()[-2].getTh();
        const double& phi = traj.end()[-2].getPhi();

        // Check if the particle has entered or jumped through the sun
        const double& rNew = traj.back().getR();
        const double& thNew = traj.back().getTh();
        const double& phiNew = traj.back().getPhi();

        // Get shortest distance from the line between the new and old points and the sun
        double dToOrigin = sqrt(pow(r * rNew, 2) * (pow(cos(thNew) * sin(th), 2)
                    - 1/2 * cos(phi - phiNew) * sin(2*th) * sin(2*thNew)
                    + pow(sin(thNew), 2) * (pow(cos(th), 2) + pow(sin(th) * sin(phi - phiNew), 2))))
            / sqrt(r*r + rNew*rNew
                    - 2*r*rNew * (cos(th) * cos(thNew) + cos(phi - phiNew) * sin(th) * sin(thNew)));

        // If the particle entered the sun, change the status.  Otherwise, continue.
        if (dToOrigin <= params.getRSun()) {
            status = BoundaryHit::Sun;
            // Put the particle at the origin!
            // TODO: this is a bit kludgy.  Should write the status to XML!                  2015-05-05 14:22
            traj.back().setR(0.0);

            return status;
        } else {
            // Made a new point, so everything needs to be updated again
            updatedB = false;
            updatedVdrift = false;
            updatedKTensor = false;

            return status;
        }
    }
}

const BoundaryHit& PPTrajectory::integrate(int n) {
	// For positive n, integrate n times, or until a boundary is hit.
	if (n > 0) {
		for (int i = 0; i < n; i++) {
			// If a boundary was hit, stop integrating
			if (step() != BoundaryHit::None) {
				break;
			}
		}
	} else {
		// For non-positive n, integrate until a boundary is hit
		int i = 0;
		while (status == BoundaryHit::None) {
			step();
			i++;
		}
	}
	
	return status;
}

PPTrajectory& PPTrajectory::writeToXML(const std::string& fname) {
	// Open the file
	std::ofstream xmlWriter(fname);

	if (xmlWriter.is_open()) {
		// Generate the XML
		xmlWriter << toXML();

		// Flush the buffer and close the file
		xmlWriter.close();
	} // Else, unable to open file.  Not sure what to do here.

	return *this;
}

std::string PPTrajectory::toXML() const {
	std::string xml;

	// Write simulation parameters
	xml += "<trajectory>\n" + params.toXML(1);
    xml += "\t<boundary>";
    if (status == BoundaryHit::Sun) {
        xml += "sun</boundary>\n";
    } else {
        xml += "heliopause</boundary>\n";
    }
	xml += "\t<points>\n";

	// Full precision version of to_string
	std::stringstream converter;
	converter << std::setprecision(std::numeric_limits<double>::digits10);

    // Get backwards time at exit
    double sMax = traj.back().getS();

    if (outputFormat == OutputFormat::All) {
        // Write all points
        // Create a reverse iterator
        auto point = traj.rbegin();
	
        // Put all the points in the stream
        for (; point != traj.rend(); point++) {
            xml += point->toXML(2, sMax);
        }
    } else if (outputFormat == OutputFormat::FirstLast) {
        // Only record first and last points
        xml += traj.back().toXML(2, sMax);
        xml += traj.front().toXML(2, sMax);
    }

	xml += "\t</points>\n</trajectory>\n\n";

	return xml;
}



