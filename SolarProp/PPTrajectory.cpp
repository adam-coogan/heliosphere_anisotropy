#include "PPTrajectory.h"

PPTrajectory::PPTrajectory(double ri, double thi, double phii, double ei, const std::string& paramFileName)
	: traj(1, PPPoint(ri, thi, phii, ei, 0)), status(BoundaryHit::None),
	generator(std::chrono::system_clock::now().time_since_epoch().count()), ndistro(0.0, 1.0),
	b(0, 0, 0), curlBoverB(0, 0, 0), kTensor(0, 0, 0, 0, 0, 0, 0, 0, 0), kpar(0), kperp(0),
   	vdrift(0, 0, 0), updatedB(false), updatedVdrift(false), updatedKTensor(false), params(paramFileName) { }

// Correct																					 2015-04-02 11:47
double PPTrajectory::vel() const {
	return cAUs * sqrt(1 - pow(params.getM() / traj.back().getE(), 2));
}

// Returns rigidity in T*AU
// Correct.  Though maybe the sign should change for electrons vs positrons?				 2015-04-02 11:42
double PPTrajectory::rigidity() const {
	// Conversion factor: 1 GV/c = 2.23*10^-11 T*AU
	return 2.23 * pow(10, -11) * sqrt(pow(traj.back().getE(), 2) - pow(params.getM(), 2));
}

// Correct																					 2015-04-02 11:44
double PPTrajectory::lambdaPar() const {
	if (rigidity() >= params.getRig0()) {
		return params.getLambda0() * (rigidity() / params.getRig0()) 
			* (1 + traj.back().getR() / params.getR0());
	} else {
		return params.getLambda0() * (1 + traj.back().getR() / params.getR0());
	}	

	// lambdaPar for the wavy paper
	/*
	return pow(params.getR0() / traj.back().getR(), 2) * pow((params.getVsw() * params.getVsw()
		+ pow(traj.back().getR() * params.getOmega() * sin(traj.back().getTh()), 2))
		/ (params.getVsw() * params.getVsw() + pow(params.getR0() * params.getOmega(), 2)), 1.0 / 2.0);
	*/
}

// Correct																					 2015-04-02 11:44
double PPTrajectory::fs() const {
	double wts2 = pow(rigidity() / params.getRig0(), 2);
	return wts2 / (1 + wts2);
}

// Correct																					 2015-04-02 11:48
double PPTrajectory::gamma() const {
	return (traj.back().getE() + 2 * params.getE0()) / (traj.back().getE() + params.getE0());
}

// Correct.  Trying to calculate vgc as in the wavy paper, rather than using curlBoverB directly. 2015-04-02 14:33
PPTrajectory& PPTrajectory::updateB() {
	if (!updatedB) {
		double r = traj.back().getR();
		double sinth = sin(traj.back().getTh());
		double costh = cos(traj.back().getTh());
		double h = heaviside(traj.back().getTh(), getThp());

		// Compute B
		b.r = params.getAc() * params.getB0() * pow(params.getR0() / r, 2) * h;
		b.phi = -params.getAc() * params.getB0() * pow(params.getR0() / r, 2) * h * params.getOmega()
			* r * sinth / params.getVsw();

		updatedB = true;
	}

	return *this;
}

// Correct																					 2015-04-02 11:38
PPTrajectory& PPTrajectory::updateKTensor() {
	if (!updatedKTensor) {
		// Make sure the magnetic field has been updated!
		updateB();

		kpar = (vel() / 3) * lambdaPar();
		kperp = 0.01 * kpar;

		// This component is independent of angle
		kTensor.rr = (b.r*b.r * kpar + b.phi*b.phi * kperp) / (b.r*b.r + b.phi*b.phi);
		kTensor.thth = kperp;
		kTensor.phiphi = (b.phi*b.phi * kpar + b.r*b.r * kperp) / (b.r*b.r + b.phi*b.phi);
		kTensor.rphi = kTensor.phir = b.r * b.phi * (kpar - kperp) / (b.r*b.r + b.phi*b.phi);

		// These should always be zero!!!
		kTensor.rth = kTensor.thr = 0;
		kTensor.phith = kTensor.thphi = 0;

		updatedKTensor = true;
	}

	return *this;
}

// Correct.  Uses wavy paper's expression for vgc.											 2015-04-02 14:47
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
		// Convenience parameter
		double gammaGC = r * params.getOmega() * sinth / params.getVsw();

		// Compute v_d^HCS
		// For convenience, initialize the parts of the r and phi components of vdHCS which are the same
		double vdHCSr = 0;
		double vdHCSphi = 0;
		// The formula for HCS drifts only applies for |d| <= 2 rL!
		if (fabs(r * costh) <= 2 * rL) {
			vdHCSr = (0.457 - 0.412 * fabs(r * costh) / rL + 0.0915 * pow(r * costh / rL, 2))
				* params.getAc() * params.getQSign() * vel() / sqrt(1.0 + gammaGC * gammaGC);
			vdHCSphi = vdHCSr;

			// Put in the factors that are different between the components
			vdHCSr *= gammaGC;
		}

		// Compute v_d^gc factor
		// I THINK the units work here, now...
		double vdgcFact = 2 * params.getAc() * vel() * rig * r / (3 * params.getB0() * params.getR0()
			* params.getQSign() * pow(1 + gammaGC * gammaGC, 2));

		// Incorporate the drift reduction factor and combine all components
		vdrift.r = fs() * (vdHCSr + vdgcFact * (-gammaGC * costh / sinth));
		vdrift.th = fs() * vdgcFact * (2 + gammaGC * gammaGC) * gammaGC;
		vdrift.phi = fs() * (vdHCSphi + vdgcFact * gammaGC * gammaGC * costh / sinth);

		std::cout << "vHCS = " << vdHCSr << ",\tvGCr = " << vdgcFact * (-gammaGC * costh / sinth) << std::endl;

		updatedVdrift = true;
	}

	return *this;
}

// Correct																					 2015-04-02 16:47
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
	// Convenience parameter
	double gammaGC = r * params.getOmega() * sinth / params.getVsw();

	// Coefficients in dr
	double dsCdr = (kpar * (2 * params.getR0() + r * (3 + gammaGC * gammaGC)) + kperp * gammaGC * gammaGC
		* (2 * params.getR0() * (2 + gammaGC * gammaGC) + r * (5 + 3 * gammaGC * gammaGC)))
		/ (r * (r + params.getR0()) * pow(1 + gammaGC * gammaGC, 2))
		- params.getVsw()
		- vdrift.r;

	//std::cout << vdrift.r * params.getDs() << std::endl;
	//std::cout << (kpar * (2 * params.getR0() + r * (3 + gammaGC * gammaGC)) + kperp * gammaGC * gammaGC
	//	* (2 * params.getR0() * (2 + gammaGC * gammaGC) + r * (5 + 3 * gammaGC * gammaGC)))
	//	/ (r * (r + params.getR0()) * pow(1 + gammaGC * gammaGC, 2)) << std::endl;

	double dWrCdr = sqrt(2 * kTensor.rr - 2 * pow(kTensor.rphi, 2) / kTensor.phiphi);

	double dWphiCdr = sqrt(2 / kTensor.phiphi) * kTensor.rphi;
	
	// Coefficients in dth
	double dsCdth = kperp * costh / (r * r * sinth) - vdrift.th / r;

	double dWthCdth = sqrt(2 * kTensor.thth) / r;
	
	// Coefficients in dphi
	double dsCdphi = (kperp - kpar) * gammaGC * (2 * params.getR0() + r * (3 + gammaGC * gammaGC))
		/ (r * r * (r + params.getR0()) * pow(1 + gammaGC * gammaGC, 2) * sinth)
		- vdrift.phi / (r * sinth);

	double dWphiCdphi = sqrt(2 * kTensor.phiphi) / (r * sinth);
	//std::cout << "dWphiCdphi = " << dWphiCdphi << std::endl;
	
	// Coefficients in dE
	double dsCdE = 2 * traj.back().getE() * params.getVsw() * gamma() / (3 * r);
	
	// Generate Wiener terms
	double dWr = ndistro(generator) * sqrt(params.getDs());
	double dWth = ndistro(generator) * sqrt(params.getDs()); 
	double dWphi = ndistro(generator) * sqrt(params.getDs());

	//std::cout << "dWphiCdphi = " << dWphiCdphi << std::endl;
	//std::cout << "|dr|\t= " << fabs(dsCdr * params.getDs() + dWrCdr * dWr + dWphiCdr * dWphi) << std::endl;
	//std::cout << "\t= " << dsCdr * params.getDs() << "\t" << dWrCdr * dWr << "\t" <<  dWphiCdr * dWphi << std::endl;
	//std::cout << "|dth|  = " << fabs(dsCdth * params.getDs() + dWthCdth * dWth) << std::endl;
	//std::cout << "|dphi| = " << fabs(dsCdphi * params.getDs() + dWphiCdphi * dWphi) << std::endl;
	
	// Move the particle to the new point
	traj.push_back(PPPoint(r + dsCdr * params.getDs() + dWrCdr * dWr + dWphiCdr * dWphi,
		traj.back().getTh() + dsCdth * params.getDs() + dWthCdth * dWth,
		phi + dsCdphi * params.getDs() + dWphiCdphi * dWphi,
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
		int i = 0;
		while (status == BoundaryHit::None) {
			step();
			i++;
		}
		std::cout << "t_exit = " << i * params.getDs() / 86400 << " days" << std::endl;
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



