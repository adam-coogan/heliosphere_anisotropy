#ifndef SPHERICALVECTOR_H
#define SPHERICALVECTOR_H

#include <cmath>
#include <iostream>

// Vector in spherical coordinates
struct SphericalVector {
	SphericalVector(double r0, double th0, double phi0) : r(r0), th(th0), phi(phi0) { };

	double r, th, phi;

	friend std::ostream& operator<<(std::ostream& os, const SphericalVector& sv) {
		os << "(v_r, v_th, v_phi) = (" << sv.r << ", " << sv.th << ", " << sv.phi << ")";
		return os;
	};

	SphericalVector& operator+=(const SphericalVector& sv) {
		r += sv.r;
		th += sv.th;
		phi += sv.phi;

		return *this;
	};

	SphericalVector& operator-=(const SphericalVector& sv) {
		r -= sv.r;
		th -= sv.th;
		phi -= sv.phi;

		return *this;
	};

	SphericalVector& operator=(const SphericalVector& sv) {
		// Check for self-assignment
		if (this != &sv) {
			r = sv.r;
			th = sv.th;
			phi = sv.phi;
		}

		return *this;
	};

	SphericalVector& operator*=(const double a) {
		r *= a;
		th *= a;
		phi *= a;

		return *this;
	};

	const SphericalVector operator+(const SphericalVector& sv) const {
		return SphericalVector(*this) += sv;
	};

	const SphericalVector operator-(const SphericalVector& sv) const {
		return SphericalVector(*this) -= sv;
	};

	template<typename T>
	const SphericalVector operator*(const T a) const {
		return SphericalVector(*this) *= a;
	};

	double magnitude2() const { return r * r + th * th + phi * phi; };

	double magnitude() const { return sqrt(magnitude2()); };
};

template<typename T>
const SphericalVector operator*(const T a, const SphericalVector& sv) {
	return sv * a;
};

#endif


