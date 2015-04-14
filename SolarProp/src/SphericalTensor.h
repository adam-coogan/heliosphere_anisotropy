#ifndef SPHERICALTENSOR_H
#define SPHERICALTENSOR_H

// Tensor in spherical coordinates
struct SphericalTensor {
	SphericalTensor(double rr0, double phiphi0, double thth0, double rphi0, double phir0, double rth0,
			double thr0, double thphi0, double phith0)
		: rr(rr0), phiphi(phiphi0), thth(thth0), rphi(rphi0), phir(phir0), rth(rth0), thr(thr0),
		thphi(thphi0), phith(phith0) { };

	double rr, phiphi, thth, rphi, phir, rth, thr, thphi, phith;
};

#endif


