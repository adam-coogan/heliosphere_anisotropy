#ifndef WAVY3D
#define WAVY3D 

#include "Basic3D.h"
#include "Point.h"
#include <algorithm>
#include <array>
#include <fstream>
#include <iterator>
#include <tuple>

#define DEBUGNM true

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

        /*!
         * Gets distance to closest point on HCS.
         * \return tuple containing the point in the HCS closest to the particle and the distance between the
         * particle and point in AU.
         * TODO: MAKE PROTECTED AFTER DONE TESTING!
         */
        std::tuple<Point, double> getLHCS() const;

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
        int itersNM = 10; // TODO: maybe set a precision goal instead?

        /*!
         * Sets up simplex used in getLHCS().
         * \return three distinct points in the HCS near the particle.  Note these are NOT sorted by distance
         * from the particle's current location.
         */
        std::array<std::tuple<Point, double>, 3> initializeSimplex() const;

        /*!
         * Helper method for getLHCS().
         * \arg pts a list of three points.  The points MUST be in the HCS (ie, they MUST have already had
         * setHCSTheta() called on them)!  pts is modified: each point's distance from the particle is
         * calculated and recorder.  These distances are then used to sort pts in ascending order.
         */
        void sortSimplexPts(std::array<std::tuple<Point, double>, 3>& pts) const;

        /*! Convenience method use to put a point in the HCS
         * \arg pt the point to put in the HCS
         * \return pt with the same r and ph, but with th updated to be the polar angle for the HCS
         */
        Point& setHCSTheta(Point& pt) const;

#if DEBUGNM
        /*!
         * Used to debug NM algorithm.  Produces a data file containing simplex at each step.  Format is:
         *  r1a,th1a,ph1a,d1a
         *  r1b,th1b,ph1b,d1b
         *  r1c,th1c,ph1c,d1c
         * where the number indicates the timestep and d is the distance to the HCS
         * \arg outPath data file name
         */
        void writeNM(const std::string& outPath, const std::string& data) const;
#endif
};

#if DEBUGNM
template<class Par>
void Wavy3D<Par>::writeNM(const std::string& outPath, const std::string& data) const {
    // Write run data to a CSV
    std::ofstream writer(outPath);

    if (writer.is_open()) {
        writer << data;
        writer.close();
    }
}
#endif

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

template<class Par>
std::array<std::tuple<Point, double>, 3> Wavy3D<Par>::initializeSimplex() const {
    // Using this tuple saves time on calculating distances and arguably makes more sense.  Too bad C++'s
    // tuple syntax is kind of poopy.  Note that clang's warning here is a bug.
    // TODO: test these guesses
    std::array<std::tuple<Point, double>, 3> pts = {
        std::make_tuple(Point(std::max(pos.getR() - 2 * params.getOmega()/params.getVsw(), params.getRSun()),
                    pos.getTh(), pos.getPh() - M_PI / 12), 0),
        std::make_tuple(Point(std::min(pos.getR() + 2 * params.getOmega()/params.getVsw(), params.getRHP()),
                    pos.getTh(), pos.getPh() - M_PI / 12), 0),
        std::make_tuple(Point(pos.getR(), pos.getTh(), pos.getPh() + M_PI / 6), 0)
    };
    
    // Put points in the HCS.  CRUCIAL!
    for (auto& p : pts) {
        setHCSTheta(std::get<0>(p));
    }

    return pts;
}

// TODO: make initial point parameters member variables.  Make method to intialize the simplex.
template<class Par>
std::tuple<Point, double> Wavy3D<Par>::getLHCS() const {
    /*if (params.getAlpha() != 0) {*/
    // Generate initial guesses
    std::array<std::tuple<Point, double>, 3> pts = initializeSimplex();

#if DEBUGNM
    std::string nmStr;
#endif

    // Loop until the simplex has converged to a satisfactory level of precision
    for (int i = 0; i < itersNM; i++) {
        // 1. Sort the simplex points
        sortSimplexPts(pts);

        // DEBUG
#if DEBUGNM
        for (auto pt : pts) {
            const Point& p = std::get<0>(pt);
            nmStr.append(std::to_string(p.getR()) + "," + std::to_string(p.getTh()) + ","
                    + std::to_string(p.getPh()) + "," + std::to_string(std::get<1>(pt)) + "\n");
        }
#endif

        // 2. Compute centroid (origin) of points 1 and 2.  Note that its theta value does not matter!
        Point ptO = (std::get<0>(pts[0]) + std::get<0>(pts[1])) / 2;

        // 3. Compute reflected point, BEING CAREFUL TO CALCULATE ITS THETA VALUE!
        Point ptR = ptO + alphaNM * (ptO - std::get<0>(pts[2]));
        setHCSTheta(ptR); // put point in the HCS.  CRUCIAL!
        double ptRDist = Point::dist(pos, ptR);

        if (ptRDist >= std::get<1>(pts[0]) && ptRDist < std::get<1>(pts[1])) {
            //std::cout << "Exit @ step 3" << std::endl;
            pts[2] = std::tuple<Point, double>(ptR, ptRDist);
        } else {
            if (ptRDist < std::get<1>(pts[0])) {
                // Try expanding if reflected point is the best so far
                Point ptE = ptO + gammaNM * (ptR - ptO);
                setHCSTheta(ptE); // put point in the HCS.  CRUCIAL!
                double ptEDist = Point::dist(pos, ptE);

                //std::cout << "ptE distance = " << ptEDist << std::endl;

                if (ptEDist < ptRDist) {
                    //std::cout << "Exit @ step 4a" << std::endl;
                    pts[2] = std::tuple<Point, double>(ptE, ptEDist);
                } else {
                    //std::cout << "Exit @ step 4b" << std::endl;
                    pts[2] = std::tuple<Point, double>(ptR, ptRDist);
                }
            } else {
                // Compute contracted point
                Point ptC = ptO + rhoNM * (std::get<0>(pts[2]) - ptO);
                setHCSTheta(ptC); // put point in the HCS.  CRUCIAL!
                double ptCDist = Point::dist(pos, ptC);

                //std::cout << "ptC distance = " << ptCDist << std::endl;

                // If contracted point is better, replace worst one with it
                if (ptCDist < std::get<1>(pts[2])) {
                    //std::cout << "Exit @ step 5" << std::endl;
                    pts[2] = std::tuple<Point, double>(ptC, ptCDist);
                } else {
                    //std::cout << "Exit @ step 6" << std::endl;
                    // Nothing worked very well.  Reduce the simplex.
                    for (int i = 1; i <= 2; i++) { // ok to hard code this!
                        std::get<0>(pts[i]) = std::get<0>(pts[0]) + sigmaNM * (std::get<0>(pts[i])
                                - std::get<0>(pts[0]));
                        setHCSTheta(std::get<0>(pts[i])); // put point in the HCS.  CRUCIAL!
                    }
                }
            }
        }
        //std::cout << std::endl;
    }

#if DEBUGNM
    writeNM("/Users/acoogan/Dropbox/heliosphere_anisotropy/nmdata/nmdata.csv", nmStr);
#endif

    // Return the best point
    return pts[0];
    /*
    } else {
        return std::make_tuple(Point(pos.getR() * sin(pos.getTh()), M_PI / 2, pos.getPh()),
                std::abs(pos.getR() * cos(pos.getTh())));
    }
    */
}

template<class Par>
void Wavy3D<Par>::sortSimplexPts(std::array<std::tuple<Point, double>, 3>& pts) const {
    // Find distances from each point to the particle
    for (auto& p : pts) {
        std::get<1>(p) = Point::dist(pos, std::get<0>(p));
    }

    // 1. Use a lambda expression to sort pts
    std::sort(std::begin(pts), std::end(pts), [](const std::tuple<Point, double>& a,
                const std::tuple<Point, double>& b){ return std::get<1>(a) < std::get<1>(b); });
}

template<class Par>
double Wavy3D<Par>::getXi(const Point& hcsPt) const {
    double tanXi = params.getOmega() * hcsPt.getR() / (params.getVsw() * sinPsi * sin(hcsPt.getTh()))
        * sqrt(pow(sin(params.getAlpha()), 2) - pow(cos(hcsPt.getTh()), 2));
    double sgnXi = sgn(cos(hcsPt.getPh() - params.getPhPhase()
                + params.getOmega() * hcsPt.getR() / params.getVsw()));

    return sgnXi * atan(tanXi);
}

// TODO: break into helper for computing HCS drift velocity?
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
    if (std::abs(pos.getTh() - M_PI / 2) <= params.getAlpha()
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


