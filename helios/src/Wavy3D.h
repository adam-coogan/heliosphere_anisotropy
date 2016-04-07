#ifndef WAVY3D
#define WAVY3D 

#include "Basic3D.h"

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
         */
        virtual double hcsExtent(double rS, double phS) const;

        // Gets distance to closest point on HCS
        virtual double getL() const;

        // Members are hidden during inheritance from a template class: these are all nondependent, and C++
        // won't look for them in the dependent base class.  Since these variables are used a lot, this seems
        // better than using this->... or TrajectoryBase<P>::.
        using Basic3D<Par>::params;
        using Basic3D<Par>::r;
        using Basic3D<Par>::th;
        using Basic3D<Par>::ph;
        using Basic3D<Par>::P;
        using Basic3D<Par>::beta;
        using Basic3D<Par>::sinPsi;
        using Basic3D<Par>::cosPsi;
        using Basic3D<Par>::speedOfLight;
        using Basic3D<Par>::gamma;
        using Basic3D<Par>::vd;
        using Basic3D<Par>::bMag;
};

template<class Par>
double Wavy3D<Par>::hcsExtent(double rS, double phS) const {
    return M_PI / 2 + asin(sin(params.getAlpha()) * sin(phS - params.phPhase
                + params.Omega() * rS / params.Vsw()));
}

template<class Par>
void Wavy3D<Par>::updateVdr() {
    // Sign of drift velocity.  Changes as the particle passes through the HCS.
    double vdSign = 1;

    // Heaviside function comes in here
    if (th > hcsExtent(r, ph)) {
        vdSign = -1;
    } else if (th == hcsExtent(r, ph)) {
        vdSign = 0;
    }

    // Gradient and curvature part of drift velocity
    double vdCoeff = 4.468e-5 * (2 * P * beta * r) // (1 GV) / (1 nT * 1 au) = 4.468e-5 au / s
        / (3 * pow(1 + gamma*gamma, 2) * sgn(params.getCharge()) * params.getAc() * params.getB0()
                * pow(params.getR0(), 2));
    vd.r = vdSign * vdCoeff * (-gamma / tan(th));
    vd.th = vdSign * vdCoeff * (2 + gamma*gamma) * gamma;
    vd.ph = vdSign * vdCoeff * gamma*gamma / tan(th);

    // Compute HCS part of drifts
    double d = fabs(r * cos(th));
    double rL = P / bMag * 0.0223; // (1 GV/c) / (1 nT) = 0.0223 au

    if (d <= 2 * rL) {
        double hcsDriftFact = params.getAc() * sgn(params.getCharge())
            * (0.457 - 0.412 * d / rL + 0.0915 * pow(d / rL, 2)) * speedOfLight*beta;
        vd.r += hcsDriftFact * sinPsi;
        vd.ph += hcsDriftFact * cosPsi;
#if DEBUG
        std::cout << "hcsDriftFact = " << hcsDriftFact << std::endl;
#endif
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


