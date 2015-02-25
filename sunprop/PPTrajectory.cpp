#include "PPTrajectory.h"

PPTrajectory::PPTrajectory(double ri, double thi, double phii, double ei, const std::string& paramFileName)
	: traj(1, PPPoint(ri, thi, phii, ei, 0)), status(BoundaryHit::None),
	generator(std::chrono::system_clock::now().time_since_epoch().count()), ndistro(0.0, 1.0),
	b(0, 0, 0), curlBoverB(0, 0, 0), kTensor(0, 0, 0, 0, 0, 0, 0, 0, 0), vdrift(0, 0, 0), updatedB(false),
	updatedVdrift(false), updatedKTensor(false), params(paramFileName) { }

double PPTrajectory::vel() const {
	return cAUs * sqrt(1 - pow(params.getM() / traj.back().getE(), 2));
}

// Returns rigidity in T*AU
double PPTrajectory::rigidity() const {
	// Conversion factor: 1 GV/c = 2.23*10^-11 T*AU
	return 2.23 * pow(10, -11) * sqrt(pow(traj.back().getE(), 2) - pow(params.getM(), 2));
}

double PPTrajectory::lambdaPar() const {
	double r = traj.back().getR();

	if (rigidity() >= params.getRig0()) {
		return params.getLambda0() * (rigidity() / params.getRig0()) * (1 + r / params.getR0());
	} else {
		return params.getLambda0() * (1 + r / params.getR0());
	}	
}

// Factor in drift velocity
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
		double costh = cos(traj.back().getTh());
		double h = heaviside(traj.back().getTh(), getThp());

		// Compute B
		// TODO: might never need to know the components of B, just curl B/|B| and |B|!!!		 2015-02-23 14:04
		b.r = params.getAc() * params.getB0() * pow(params.getR0() / r, 2) * h;
		b.phi = -params.getAc() * params.getB0() * pow(params.getR0() / r, 2) * h * params.getOmega()
			* r * sinth / params.getVsw();

		// TODO: what happens to the curl at the HCS???											 2015-02-23 14:04
		//		-Assuming that H' vanishes everywhere.											 2015-02-25 11:30
		curlBoverB.r = curlBoverB.th = curlBoverB.phi = params.getAc() * params.getOmega() * h
			/ pow(pow(params.getVsw(), 2) + pow(r * params.getOmega() * sinth, 2), 3.0/2.0);
		curlBoverB.r *= -costh * (2 * pow(params.getVsw(), 2) + pow(r * params.getOmega() * sinth, 2));
		curlBoverB.th *= (2 * pow(params.getVsw(), 2) * params.getOmega() * sinth + r*r
			* pow(params.getOmega() * sinth, 3)) / params.getOmega();
		curlBoverB.phi *= r * params.getVsw() * params.getOmega() * sinth * costh;

		updatedB = true;
	}

	return *this;
}

PPTrajectory& PPTrajectory::updateKTensor() {
	if (!updatedKTensor) {
		// Make sure the magnetic field has been updated!
		updateB();

		double kpar = (vel() / 3) * lambdaPar();
		double kperp = 0.01 * kpar;

		// This component is independent of angle
		kTensor.thth = kperp;
		kTensor.rr = (b.r*b.r * kpar + b.phi*b.phi * kperp) / (b.r*b.r + b.phi*b.phi);
		kTensor.phiphi = (b.phi*b.phi * kpar + b.r*b.r * kperp) / (b.r*b.r + b.phi*b.phi);
		kTensor.rphi = kTensor.phir = b.r * b.phi * (kpar - kperp) / (b.r*b.r + b.phi*b.phi);

		// These should always be zero!!!
		//kTensor.rth = kTensor.thr = 0;
		//kTensor.phith = kTensor.thphi = 0;

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
		double phi = traj.back().getPhi();
		double rig = rigidity();
		// Larmor radius
		double rL = rig / b.magnitude();

		// Compute v_d^HCS
		// For speed, initialize the parts of the r and phi components of vdHCS which are the same
		double vdHCSr = (0.457 - 0.412 * abs(r * costh) / rL + 0.0915 * r*r * costh*costh / (rL*rL))
			* params.getAc() * params.getQSign() * vel() / sqrt(1.0 + pow(r * sinth * params.getOmega()
			/ params.getVsw(), 2));
		double vdHCSphi = vdHCSr;

		// Put in the factors that are different between the components
		vdHCSr *= r * params.getOmega() * sinth / params.getVsw();

		// Compute v_d^gc factor
		double vdgcFact = params.getQSign() * rig * vel() / (3 * b.magnitude());

		// TODO: not sure how to handle the fact that |B| = 0 here...

		// Incorporate the drift reduction factor and combine all components
		vdrift.r = fs() * (vdHCSr + vdgcFact * curlBoverB.r);
		vdrift.th = fs() * vdgcFact * curlBoverB.th;
		vdrift.phi = fs() * (vdHCSphi + vdgcFact * curlBoverB.phi);

		updatedVdrift = true;
	}

	return *this;
}

PPTrajectory& PPTrajectory::step() {
	// Update all parameters
	updateB();
	updateKTensor();
	updateVdrift();
	
	// Make stuff clearer and maybe slightly faster
	double r = traj.back().getR();
	double sinth = sin(traj.back().getTh());
	double costh = cos(traj.back().getTh());
	double phi = traj.back().getPhi();

	// Coefficients in dr
	double dsCdr = (2 / r + 2 * pow(params.getVsw(), 2) * (1 / (r * pow(params.getVsw(), 2) + r*r*r
		* pow(params.getOmega() * sinth, 2)) - 1 / (r * pow(params.getVsw(), 2) + params.getDiffFact() * r*r*r
		* pow(params.getOmega() * sinth, 2)) - 1 / (r + params.getR0()))) * kTensor.rr - vdrift.r
		- params.getVsw();
	double dWrCdr = sqrt(2 * kTensor.rr - 2 * pow(kTensor.rphi, 2) / kTensor.phiphi);
	double dWphiCdr = sqrt(2 / kTensor.phiphi) * kTensor.rphi;
	
	// Coefficients in dth
	double dsCdth = -vdrift.th / r + kTensor.thth * costh / (r*r * sinth);
	double dWthCdth = sqrt(2 * kTensor.thth) / r;
	
	// Coefficients in dphi
	double dsCdphi = kTensor.rphi / sinth * ((r + 2 * params.getR0()) / (r*r * (r + params.getR0()))
		- 2 * pow(params.getOmega(), 2) / (pow(r * params.getOmega(), 2) + pow(params.getVsw() / sinth, 2)))
		- vdrift.phi / (r * sinth);
	double dWphiCdphi = sqrt(2 * kTensor.phiphi) / (r * sinth);
	
	// Coefficients in dE
	double dsCdE = 2 * traj.back().getE() * params.getVsw() * gamma() / (3 * r);
	
	// Generate Wiener terms
	double dWr = ndistro(generator) * sqrt(params.getDs());
	double dWth = ndistro(generator) * sqrt(params.getDs()); 
	double dWphi = ndistro(generator) * sqrt(params.getDs());
	
	// Move the particle to the new point
	traj.push_back(PPPoint(r + (dsCdr * params.getDs() + dWrCdr * dWr + dWphiCdr * dWphi),
		traj.back().getTh() + (dsCdth * params.getDs() + dWthCdth * dWth),
		phi + (dsCdphi * params.getDs() + dWphiCdphi * dWphi),
		traj.back().getE() + dsCdE * params.getDs(),
		traj.back().getS() + params.getDs()));

	// If the particle has reached or passed the heliopause
	if (traj.back().getR() >= params.getRHP()) {
		status = BoundaryHit::Heliopause;
	} else if (traj.back().getR() <= params.getRSun()) {
		// If the particle enters the sun
		status = BoundaryHit::Sun;
	}

	// Made a new point, so everything needs to be updated again
	updatedB = false;
	updatedVdrift = false;
	updatedKTensor = false;

	return *this;
}

BoundaryHit PPTrajectory::integrate(int n) {
	// For positive n, integrate n times, or until a boundary is hit.
	if (n > 0) {
		for (int i = 0; i < n; i++) {
			step();
			// If a boundary is hit, stop integrating
			if (status != BoundaryHit::None) {
				break;
			}
		}
	} else {
		// For negative n, integrate until a boundary is hit
		while (status == BoundaryHit::None) {
			step();
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
	xml += "<trajectory>\n" + params.toXML();
	xml += "\t<points>\n";

	// Full precision version of to_string
	std::stringstream converter;
	converter << std::setprecision(std::numeric_limits<double>::digits10);

	// Create a reverse iterator
	auto point = traj.rbegin();
	double sMax = point->getS();
	
	// Put all the points in the stream
	for (; point != traj.rend(); point++) {
		converter.str("");
		converter << point->getR();
		xml += "\t\t<point>\n\t\t\t<r>" + converter.str() + "</r>\n";

		converter.str("");
		converter << point->getTh();
		xml += "\t\t\t<th>" + converter.str() + "</th>\n";

		converter.str("");
		converter << point->getPhi();
		xml += "\t\t\t<phi>" + converter.str() + "</phi>\n";

		converter.str("");
		converter << point->getE();
		xml += "\t\t\t<e>" + converter.str() + "</e>\n";

		converter.str("");
		converter << sMax - point->getS();
		xml += "\t\t\t<t>" + converter.str() + "</t>\n\t\t</point>\n";
	}

	xml += "\t</points>\n</trajectory>\n\n";

	return xml;
}



