#ifndef WAVY3D
#define WAVY3D 

#include "Basic3D.h"
#include "Point.h"
#include <algorithm>
#include <array>
#include <tuple>

/*!
 * Adds a wavy HCS into the basic heliosphere model.  See doi:10.1007/s10509-012-1003-z (Strauss et al 2012)
 * for details.
 * \param Par the parameter container type.  Should extend Basic3D and add tilt angle parameter.
 */
template<class Par>
class Wavy3D : public Basic3D<Par> {
    public:
        /*!
         * Constructor is the same as for Basic3D.  The parameter class handles loading the tilt angle.
         */
        Wavy3D(std::string paramFileName) : Basic3D<Par>(paramFileName) { };

    protected:
        // Needs to be totally rewritten
        virtual void updateVdr();

        /*!
         * This method is virtual since there are difference ways to calculate this.
         * \return angular extent of HCS at given radius and azimuthal angle.
         * TODO: should this take a Point instead?  Probably not.
         */
        virtual double hcsExtent(double rS, double phS) const;

        /*!
         * Gets distance to closest point on HCS.
         * \return tuple containing the point in the HCS closest to the particle and the distance between the
         * particle and point in AU.
         */
        std::tuple<Point, double> getLHCS() const;

        /*! 
         * \return xi, the angle between tangent to HCS and radial line.  The correct sign is also determined.
         */
        virtual double getXi(const Point& hcsPt) const;

        // Members are hidden during inheritance from a template class: these are all nondependent, and C++
        // won't look for them in the dependent base class.  Since these variables are used a lot, this seems
        // better than using this->... or TrajectoryBase<P>::.
        using Basic3D<Par>::params;
        using Basic3D<Par>::pos;
        using Basic3D<Par>::P;
        using Basic3D<Par>::beta;
        using Basic3D<Par>::sinPsi;
        using Basic3D<Par>::cosPsi;
        using Basic3D<Par>::speedOfLight;
        using Basic3D<Par>::gamma;
        using Basic3D<Par>::vd;
        using Basic3D<Par>::bMag;

    private:
        //! Default Nelder-Mead parameters from Wikipedia.  The NM algorithm is how I implement getLHCS().
        const double alphaNM = 1; //!Used for reflected point
        const double gammaNM = 2; //!Used for expanded point
        const double rhoNM = 1/2; //!Used for contracted point
        const double sigmaNM = 1/2; //!Used for reducing the simplex

        /*! Convenience method use to put a point in the HCS
         * \arg pt the point to put in the HCS
         * \return pt with the same r and ph, but with th updated to be the polar angle for the HCS
         */
        Point& setHCSTheta(Point& pt) const;
};

template<class Par>
double Wavy3D<Par>::hcsExtent(double rS, double phS) const {
    return M_PI / 2 + asin(sin(params.getAlpha()) * sin(phS - params.getPhPhase()
                + params.getOmega() * rS / params.getVsw()));
}

template<class Par>
Point& Wavy3D<Par>::setHCSTheta(Point& pt) const {
    pt.setTh(hcsExtent(pt.getR(), pt.getPh()));
    return pt;
}

// TODO: write this
// TODO: return closest point's coordinates
// by TrajectoryBase.
template<class Par>
std::tuple<Point, double> Wavy3D<Par>::getLHCS() const {
    // TODO: figure out how to make reasonable initial guesses!!!  Initialize points with default constructor
    // for now.
    // Using this tuple saves time on calculating distances and arguably makes more sense.  Too bad C++'s
    // tuple syntax is kind of poopy.
    std::array<std::tuple<Point, double>, 3> pts;

    // TODO: make this a member variable.  Should be precision goal rather than number of iterations
    int nIter = 10;

    // Loop until the simplex has converged to a satisfactory level of precision
    for (int i = 0; i < nIter; ++i) {
        // Find distances from each point to the particle
        for (auto& p : pts) {
            std::get<1>(p) = Point::dist(pos, std::get<0>(p));
        }

        // 1. Use a lambda expression to sort pts
        std::sort(pts.begin(), pts.end(), [this](const std::tuple<Point, double>& a,
                    const std::tuple<Point, double>& b){ return std::get<1>(a) < std::get<1>(b); });

        // 2. Compute centroid (origin) of points 1 and 2.  Note that its theta value does not matter!
        Point ptO = (std::get<0>(pts[0]) + std::get<0>(pts[1])) / 2;

        // 3. Compute reflected point, BEING CAREFUL TO CALCULATE ITS THETA VALUE!
        Point ptR = ptO + alphaNM * (ptO - std::get<0>(pts[2]));
        setHCSTheta(ptR); // put point in the HCS.  CRUCIAL!
        double ptRDist = Point::dist(pos, ptR);

        if (ptRDist >= std::get<1>(pts[0]) && ptRDist < std::get<1>(pts[1])) {
            pts[2] = std::tuple<Point, double>(ptR, ptRDist);
        } else {
            if (ptRDist < std::get<1>(pts[0])) {
                // Try expanding if reflected point is the best so far
                Point ptE = ptO + gammaNM * (ptR - ptO);
                setHCSTheta(ptE); // put point in the HCS.  CRUCIAL!
                double ptEDist = Point::dist(pos, ptE);

                if (ptEDist < ptRDist) {
                    pts[2] = std::tuple<Point, double>(ptE, ptEDist);
                } else {
                    pts[2] = std::tuple<Point, double>(ptR, ptRDist);
                }
            } else {
                // Compute contracted point
                Point ptC = ptO + rhoNM * (std::get<0>(pts[2]) - ptO);
                setHCSTheta(ptC); // put point in the HCS.  CRUCIAL!
                double ptCDist = Point::dist(pos, ptC);

                // If contracted point is better, replace worst one with it
                if (ptCDist < std::get<1>(pts[2])) {
                    pts[2] = std::tuple<Point, double>(ptC, ptCDist);
                } else {
                    // Nothing worked very well.  Reduce the simplex.
                    for (auto& p : {pts[1], pts[2]}) { // ok to hard code this!
                        std::get<0>(p) = std::get<0>(pts[0]) + sigmaNM * (std::get<0>(p)
                                - std::get<0>(pts[0]));
                        setHCSTheta(std::get<0>(p)); // put point in the HCS.  CRUCIAL!
                    }
                }
            }
        }
    }

    return std::tuple<Point, double>(Point(), 0); // placeholder
}

template<class Par>
double Wavy3D<Par>::getXi(const Point& hcsPt) const {
    double tanXi = params.Omega() * hcsPt.getR() / (params.getVsw() * sinPsi * sin(hcsPt.getTh()))
        * sqrt(pow(sin(params.getAlpha()), 2) - pow(cos(hcsPt.getTh()), 2));
    double sgnXi = sgn(cos(hcsPt.getPh() - params.getPhPhase()
                + params.getOmega() * hcsPt.getR() / params.getVsw()));

    return sgnXi * atan(tanXi);
}

template<class Par>
void Wavy3D<Par>::updateVdr() {
    // Sign of drift velocity.  Changes as the particle passes through the HCS.
    double vdSign = 1;

    // Heaviside function comes in here
    if (pos.getTh() > hcsExtent(pos.getR(), pos.getPh())) {
        vdSign = -1;
    } else if (pos.getTh() == hcsExtent(pos.getR(), pos.getPh())) {
        vdSign = 0;
    }

    // Gradient and curvature part of drift velocity
    double vdCoeff = 4.468e-5 * (2 * P * beta * pos.getR()) // (1 GV) / (1 nT * 1 au) = 4.468e-5 au / s
        / (3 * pow(1 + gamma*gamma, 2) * sgn(params.getCharge()) * params.getAc() * params.getB0()
                * pow(params.getR0(), 2));
    vd.r = vdSign * vdCoeff * (-gamma / tan(pos.getTh()));
    vd.th = vdSign * vdCoeff * (2 + gamma*gamma) * gamma;
    vd.ph = vdSign * vdCoeff * gamma*gamma / tan(pos.getTh());

    // Larmor radius
    double rL = P / bMag * 0.0223; // (1 GV/c) / (1 nT) = 0.0223 au

    // Check whether it's even possible for L < 2 r_L to be true
    if (abs(pos.getTh() - M_PI / 2) <= params.getAlpha()
            || (pos.getTh() < M_PI / 2 - params.getAlpha()
                && pos.getR() * cos(params.getAlpha() + pos.getTh()) <= 2 * rL)
            || (pos.getTh() > M_PI / 2 + params.getAlpha()
                && -pos.getR() * cos(pos.getTh() - params.getAlpha()) <= 2 * rL)) {
        // If it is possible, actually compute L
        Point hcsPoint;
        double hcsL;
        std::tie(hcsPoint, hcsL) = getLHCS(); // unpacks the tuple returned by getLHCS()

        // If L is really less than 2 r_L, compute HCS drift velocity
        if (hcsL <= 2 * rL) {
            // Standard approximation
            double hcsDriftFact = params.getAc() * sgn(params.getCharge())
                * (0.457 - 0.412 * hcsL / rL + 0.0915 * pow(hcsL / rL, 2)) * speedOfLight*beta;
            
            // Get the angle xi
            double xi = getXi(hcsPoint);

            // This is different with a wavy HCS.
            vd.r += hcsDriftFact * sinPsi * cos(xi);
            vd.th += hcsDriftFact * sin(xi);
            vd.ph += hcsDriftFact * cosPsi * cos(xi);
#if DEBUG
            std::cout << "hcsDriftFact = " << hcsDriftFact << std::endl;
#endif
        }
    }

    // Drift reduction factor
    double fs = 10 * pow(P / params.getRig0(), 2) / (1 + 10 * pow(P / params.getRig0(), 2));

    vd.r *= fs;
    vd.th *= fs;
    vd.ph *= fs;

#if DEBUG
    std::cout << "vdCoeff = " << vdCoeff
        << "\nvdr = " << vd.r
        << "\nvdth = " << vd.th
        << "\nvdph = " << vd.ph << std::endl;
#endif
}

#endif


