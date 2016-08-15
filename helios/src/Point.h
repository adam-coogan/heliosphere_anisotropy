#ifndef POINT
#define POINT 

#include <cmath>
#include <ostream>
#include <stdexcept>

class Point {
    public:
        //! Default constructor
        Point();

        //! Copy constructor
        Point(const Point& pt);

        Point(double rn, double thn, double phn);

        std::tuple<double, double, double> getXYZ() const;

        // Getters
        double getR() const;
        double getTh() const;
        double getPh() const;

        /*!
         * \arg rn new radial coordinate for the point
         * \return this point with r set to rn.  If rn < 0, th is set to pi + th.
         */
        void setR(double rn);

        /*!
         * \arg thn new polar angle for the point
         * \return this point with th set to thn.  Ensures th and ph lie in [0, pi] and [0, 2 pi].
         */
        void setTh(double thn);

        /*!
         * \arg phn new azimuthal angle for the point
         * \return this point with ph set to phn.  Ensures ph lies in [0, 2 pi].
         */
        void setPh(double phn);

        //! Resets all coordinates, renormalizing if necessary
        // TODO: replace with =
        void set(double rn, double thn, double phn);

        //! Incrementers
        void incrR(double deltaR);
        void incrTh(double deltaTh);
        void incrPh(double deltaPh);
        
        //! Assignment operator
        Point& operator=(const Point& pt);
        //! Compound assignment operators
        Point& operator+=(const Point& pt);
        Point& operator-=(const Point& pt);
        template<typename T>
        Point& operator*=(const T& scalar);
        template<typename T>
        Point& operator/=(const T& scalar);
        //! Binary arithmetic operators
        Point operator+(const Point& pt) const;
        // TODO: implement this.  It'll (very slightly) clean up the increments in the simulation.
        //Point operator+(const std::tuple<double, double, double>& pt) const;
        Point operator-(const Point& pt) const;
        template<typename T>
        Point operator*(const T& scalar) const;
        template<typename T>
        Point operator/(const T& scalar) const;
        //! Comparison operators
        bool operator==(const Point& pt) const;
        bool operator!=(const Point& pt) const;

        //! Gets distance between two points
        static double dist(const Point& pt1, const Point& pt2);

    private:
        double r;
        double th;
        double ph;

        /*!
         * \return this point with ph renormalized to lie within [0, 2 pi]
         */
        void renormalizePh();
};

//! Function for printing a point
std::ostream& operator<<(std::ostream& os, const Point& pt);

template<typename T>
Point& Point::operator*=(const T& scalar) {
    setR(scalar * r);
    
    /*
    // NM multiplication
    set(scalar * r, scalar * th, scalar * ph);
    */

    return *this;
}

template<typename T>
Point& Point::operator/=(const T& scalar) {
    if (scalar != 0) {
        return *this *= 1.0 / scalar;
    } else {
        throw std::domain_error("Point::operator/=: cannot divide by zero.");
    }
}

template<typename T>
Point Point::operator*(const T& scalar) const {
    return Point(*this) *= scalar;
}

template<typename T>
Point Point::operator/(const T& scalar) const {
    return Point(*this) /= scalar;
}

// Can't be a member function since the scalar is on the LHS
template<typename T>
Point operator*(const T& scalar, const Point& pt) {
    return Point(pt) *= scalar;
}

#endif


