#ifndef BASIC3D
#define BASIC3D

//#include "Basic3DParams.h"
#include "DiffusionTensor.h"
#include "DriftVelocity.h"
#include "MagneticField.h"
#include "TrajectoryBase.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <ctime>
#include <random>

/*!
 * Gets the sign of a number.
 * \param val value to get sign of.
 * \return sign of val, or 0 if val is 0.
 */
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// Makes the random number generator syntax cleaner
typedef boost::mt19937 MT19937;
typedef boost::normal_distribution<double> NDistribution;
typedef boost::variate_generator<MT19937, NDistribution> VariateGenerator;

/*!
 * Implements the basic 3D model found in doi:10.1088/0004-637X/735/2/83 (Strauss et al).
 * \param Par the parameter container type.  A container of type Basic3D should be used with this class.
 */
template<class Par>
class Basic3D : public TrajectoryBase<Par> {
    public:
        //! Override base class constructor to set up random number generator.
        Basic3D(std::string paramFileName);

        /*!
         * Advances simulation by one timestep.
         */
        virtual Status step();

        virtual Status getStatus();

    protected:
        //! Jupiter's angular position
        double phJup;

        //! Magnitude of magnetic field
        double bMag;
        //! Diffusion tensor.
        DiffusionTensor k;
        //! Drift velocity.
        DriftVelocity vd;

        //! Mersenne twister PRNG
        MT19937 engine;
        //! Normal distribution
        NDistribution distribution;
        //! Uses \a engine and \a distribution to generate Gaussian random numbers with mean 0 and standard
        //! deviation 1.
        VariateGenerator generator;

        /*!
         * Updates magnitude of magnetic field.
         */
        virtual void updateB();

        /*!
         * Updates diffusion tensor.
         */
        virtual void updateK();

        /*!
         * Updates drift velocity.
         */
        virtual void updateVdr();

        /*!
         * Updates convenience variables, magnetic field and diffusion tensor, in that order.
         */
        void updateVars();

        /*!
         * Updates angular variable associated with magnetic field.
         */
        void updatePsi();

        /*!
         * Updates variable associated with magnetic field.
         */
        void updateGamma();

        // Members are hidden during inheritance from a template class: these are all nondependent, and C++
        // won't look for them in the dependent base class.  Since these variables are used a lot, this seems
        // better than using this->... or TrajectoryBase<P>::.
        using TrajectoryBase<Par>::status;
        using TrajectoryBase<Par>::params;
        using TrajectoryBase<Par>::s;
        using TrajectoryBase<Par>::r;
        using TrajectoryBase<Par>::th;
        using TrajectoryBase<Par>::ph;
        using TrajectoryBase<Par>::ek;
        using TrajectoryBase<Par>::P;
        using TrajectoryBase<Par>::beta;
        using TrajectoryBase<Par>::tanPsi;
        using TrajectoryBase<Par>::sinPsi;
        using TrajectoryBase<Par>::cosPsi;
        using TrajectoryBase<Par>::updatePsi;
        using TrajectoryBase<Par>::speedOfLight;
        using TrajectoryBase<Par>::gamma;

    private:
        void updateKHelper();
};

template<class Par>
Basic3D<Par>::Basic3D(std::string paramFileName) : TrajectoryBase<Par>(paramFileName),
        phJup(params.getPh0Jup()), engine(static_cast<unsigned int>(std::random_device{}())),
        distribution(0, 1), generator(engine, distribution) { }

template<class Par>
void Basic3D<Par>::updateVars() {
    TrajectoryBase<Par>::updateVars();

    updatePsi();
    updateGamma();
    // Only update diffusion tensor, etc after updating all other variables.  
    updateB();
    updateK();
    updateVdr();
}

template<class T>
void Basic3D<T>::updatePsi() {
    tanPsi = params.getOmega() * (r - params.getRSun()) * sin(th) / params.getVsw();
    cosPsi = 1 / sqrt(1 + tanPsi*tanPsi);
    sinPsi = sqrt(1 - cosPsi*cosPsi);
}

template<class T>
void Basic3D<T>::updateGamma() {
    gamma = r * params.getOmega() * sin(th) / params.getVsw();
}

template<class Par>
Status Basic3D<Par>::getStatus() {
    TrajectoryBase<Par>::getStatus();

    // Check whether particle ended in Jupiter's magnetosphere
    if (r > params.getRBeginJup() && r < params.getREndJup()
            && ph > phJup - params.getDphJup() && ph < phJup + params.getDphJup()
            && th > params.getThJup() - params.getDthJup() && th < params.getThJup() + params.getDthJup()) {
        status = Status::Jupiter;
    }

    return status;
}

template<class Par>
Status Basic3D<Par>::step() {
    // Don't both stepping if simulation has finished
    if (status == Status::Running) {
        // Updates everything
        updateVars();

        // Compute all the differentials
        double dr_ds = 2 / r * k.rr + k.rr_dr - (params.getVsw() + vd.r);
        double dr_dWr = sqrt(2 * k.rr - 2 * k.rph*k.rph / k.phph);
        double dr_dWph = k.rph * sqrt(2 / k.phph);

        double dth_ds = 1 / (r*r * tan(th)) * k.thth - vd.th / r;
        double dth_dWth = sqrt(2 * k.thth) / r;

        double dph_ds = 1 / (r*r * sin(th)) * k.rph + 1 / (r * sin(th)) * k.rph_dr - vd.ph / (r * sin(th));
        double dph_dWph = sqrt(2 * k.phph) / (r * sin(th));

        double dek_ds = 2 * params.getVsw() / (3 * r) * beta*beta * (ek + params.getM());

        // Generate the Wiener terms
        double dWr = generator() * sqrt(params.getDs());
        double dWth = generator() * sqrt(params.getDs());
        double dWph = generator() * sqrt(params.getDs());

#if DEBUG
        std::cout << "dr_ds = " << dr_ds << std::endl;
        std::cout << "dr_dWr = " << dr_dWr << std::endl;
        std::cout << "dr_dWph = " << dr_dWph << std::endl;
        std::cout << "dth_ds = " << dth_ds << std::endl;
        std::cout << "dth_dWth = " << dth_dWth << std::endl;
        std::cout << "dph_ds = " << dph_ds << std::endl;
        std::cout << "dph_dWph = " << dph_dWph << std::endl;
        std::cout << "dE_ds = " << dek_ds << std::endl;
#endif

        // Move the particle
        r += dr_ds * params.getDs() + dr_dWr * dWr + dr_dWph * dWph;
        th += dth_ds * params.getDs() + dth_dWth * dWth;
        ph += dph_ds * params.getDs() + dph_dWph * dWph;
        ek += dek_ds * params.getDs();
        // Move Jupiter
        phJup -= params.getOmegaJup() * params.getDs();
        // Advance the backwards clock!
        s += params.getDs();

        // Renormalize the angles
        TrajectoryBase<Par>::renormalize();

        // Check the status
        getStatus();
    }

    return status;
}

template<class Par>
void Basic3D<Par>::updateB() {
    // Magnitude of magnetic field
    bMag = params.getB0() / (sqrt(1 + pow(1 - params.getRSun(), 2)) * r*r * cosPsi);
}

template<class Par>
void Basic3D<Par>::updateK() {
    // Temporarily increment r to numerically differentiate k
    r += params.getDeltar();
    // Update convenience variables that depend on r
    updatePsi();
    updateKHelper();
    double krr_r = k.rr;
    double krph_r = k.rph;

    // Compute tensor elements at actual point
    r -= params.getDeltar();
    updatePsi();
    updateKHelper();

    // Compute r derivatives
    k.rr_dr = (krr_r - k.rr) / params.getDeltar();
    k.rph_dr = (krph_r - k.rph) / params.getDeltar();

#if DEBUG
    std::cout << "Krr = " << k.rr
        << "\n\tKrr2 = " << krr_r
        << "\nKphph = " << k.phph
        << "\nKrph = " << k.rph
        << "\n\tKrph2 = " << krph_r
        << "\nKthth = " << k.thth
        << "\ndKrr_dr = " << k.rr_dr
        << "\n\tdelta_r = " << params.getDeltar()
        << "\ndKrph_dr = " << k.rph_dr << std::endl;
#endif
}

template<class Par>
void Basic3D<Par>::updateKHelper() {
    // Parallel and perpendicular diffusion components
    double kpar = (speedOfLight*beta / 3) * params.getLambda0() * (1 + r / params.getRRefLambda())
        * (P > params.getRig0() ? P / params.getRig0() : 1);
    double kperp = params.getKperp_kpar() * kpar;
#if DEBUG
    std::cout << "kpar = " << kpar << std::endl;
    std::cout << "kperp = " << kperp << std::endl;
    std::cout << "cosPsi = " << cosPsi << std::endl;
    std::cout << "sinPsi = " << sinPsi << std::endl;
#endif

    k.rr = kpar * cosPsi*cosPsi + kperp * sinPsi*sinPsi; // TODO: make an abstact base class for 3D models...
    k.phph = kpar * sinPsi*sinPsi + kperp * cosPsi*cosPsi; // kperp = kperp,r
    k.rph = (kperp - kpar) * cosPsi * sinPsi; // kperp = kperp,r
    k.thth = kperp; // kperp = kperp,th
}

template<class Par>
void Basic3D<Par>::updateVdr() {
    // Sign of drift velocity.  Changes as the particle passes through the HCS.
    double vdSign = 1;

    // Heaviside function comes in here
    if (th > M_PI / 2) {
        vdSign = -1;
    } else if (th == M_PI / 2) {
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


