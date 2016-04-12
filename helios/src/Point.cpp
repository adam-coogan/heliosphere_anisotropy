#include "Point.h"

// Constructors

// Default point is the origin
Point::Point() { }

Point::Point(const Point& pt) {
    *this = pt;
}

Point::Point(double rn, double thn, double phn) {
    set(rn, thn, phn);
}

// Getters

double Point::getR() const {
    return r;
}

double Point::getTh() const {
    return th;
}

double Point::getPh() const {
    return ph;
}

// Setters

void Point::setR(double rn) {
    if (rn < 0) {
        r = std::abs(rn);
        setTh(M_PI + th);
    } else {
        r = rn;
    }
}

void Point::setTh(double thn) {
    th = thn;

    // Renormalize th
    while (th > M_PI) {
        th = 2 * M_PI - th;
        ph = ph - M_PI;
    }

    while (th < 0) {
        th = -th;
        ph += M_PI;
    }

    // Renormalize ph
    renormalizePh();
}

void Point::setPh(double phn) {
    ph = phn;

    // Renormalize ph
    renormalizePh();
}

void Point::set(double rn, double thn, double phn) {
    setR(rn);
    ph = phn; // Ok since setTh() calls renormalizePh()
    setTh(thn);
}

// Incrementers

void Point::incrR(double deltaR) {
    setR(r + deltaR);
}

void Point::incrTh(double deltaTh) {
    setTh(r + deltaTh);
}

void Point::incrPh(double deltaPh) {
    setPh(r + deltaPh);
}

// Function to renormalize ph
void Point::renormalizePh() {
    while (ph > 2 * M_PI) {
        ph -= 2 * M_PI;
    }

    while (ph < 0) {
        ph += 2 * M_PI;
    }
}

// Assignment operators
Point& Point::operator=(const Point& pt) {
    // Check for self assignment, though it doesn't matter for this simple class
    if (this != &pt) {
        set(pt.getR(), pt.getPh(), pt.getPh());
    }

    return *this;
}

Point& Point::operator+=(const Point& pt) {
    set(r + pt.getR(), th + pt.getTh(), ph + pt.getPh());

    return *this;
}

Point& Point::operator-=(const Point& pt) {
    set(r - pt.getR(), th - pt.getTh(), ph - pt.getPh());

    return *this;
}

template<typename T>
Point& Point::operator*=(const T& scalar) {
    set(scalar * r, scalar * th, scalar * ph);

    return *this;
}

Point Point::operator+(const Point& pt) const {
    return Point(*this) += pt;
}

Point Point::operator-(const Point& pt) const {
    return Point(*this) -= pt;
}

template<typename T>
Point Point::operator*(const T& scalar) const {
    return Point(*this) *= scalar;
}

bool Point::operator==(const Point& pt) const {
    return (r == pt.getR()) && (th == pt.getTh()) && (ph == pt.getPh());
}

bool Point::operator!=(const Point& pt) const {
    return !(*this == pt);
}


