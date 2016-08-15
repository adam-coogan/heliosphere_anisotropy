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

std::tuple<double, double, double> Point::getXYZ() const {
    return std::make_tuple(r * sin(th) * cos(ph), r * sin(th) * sin(ph), r * cos(th));
}

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
        ph = ph + M_PI; // Ok since setTh() calls renormalizePh()
        setTh(M_PI - th);
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
    ph = phn; // Ok since setTh() calls renormalizePh()
    setTh(thn);
    setR(rn);
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
        set(pt.getR(), pt.getTh(), pt.getPh());
    }

    return *this;
}

Point& Point::operator+=(const Point& pt) {
    // Do this in cartesian coordinates
    double newX = r * sin(th) * cos(ph) + pt.getR() * sin(pt.getTh()) * cos(pt.getPh());
    double newY = r * sin(th) * sin(ph) + pt.getR() * sin(pt.getTh()) * sin(pt.getPh());
    double newZ = r * cos(th) + pt.getR() * cos(pt.getTh());

    r = sqrt(newX*newX + newY*newY + newZ*newZ);
    th = acos(newZ / r);
    ph = atan2(newY, newX); // The image of atan2 is [-pi, pi], so must renormalize ph
    renormalizePh();

    /*
    // NM addition
    set(r + pt.getR(), th + pt.getTh(), ph + pt.getPh());
    */

    return *this;
}

Point& Point::operator-=(const Point& pt) {
    return (*this) += (-1 * pt);

    /*
    // NM addition
    set(r - pt.getR(), th - pt.getTh(), ph - pt.getPh());

    return *this;
    */
}

Point Point::operator+(const Point& pt) const {
    return Point(*this) += pt;
}

Point Point::operator-(const Point& pt) const {
    return Point(*this) -= pt;
}

bool Point::operator==(const Point& pt) const {
    return (r == pt.getR()) && (th == pt.getTh()) && (ph == pt.getPh());
}

bool Point::operator!=(const Point& pt) const {
    return !(*this == pt);
}

std::ostream& operator<<(std::ostream& os, const Point& pt) {
    return os << "(" << pt.getR() << ", " << pt.getTh() << ", " << pt.getPh() << ")";
}

double Point::dist(const Point& pt1, const Point& pt2) {
    return sqrt(pow(pt1.getR(), 2) + pow(pt2.getR(), 2) - 2 * pt1.getR() * pt2.getR()
            * (sin(pt1.getTh()) * sin(pt2.getTh()) * cos(pt1.getPh() - pt2.getPh())
                + cos(pt1.getTh()) * cos(pt2.getTh())));
}


