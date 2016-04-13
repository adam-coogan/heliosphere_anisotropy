#ifndef TRAJECTORYBASE
#define TRAJECTORYBASE

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include "Point.h"

/*!
 * Describes a trajectory's status.  Contains Jupiter since there's no way to extend enums.
 */
enum class Status {Sun, Heliopause, Running, Jupiter};

/*!
 * Contains data (coordinates, energy, convenience variables, diffusion tensor, magnetic field) and functions
 * relevant for computing a trajectory.
 * T is the parameter container type.
 */
template<class T>
class TrajectoryBase {
    public:
        TrajectoryBase(const std::string& paramFileName);

        /*!
         * Reinitializes the simulation using the same parameters as it was constructed with.
         * \return the status of the reinitialized simulation.
         */
        Status initialize();

        /*!
         * Reinitializes the simulation using new parameters.
         * \param paramFileName path to new file to load parameters from.
         * \return the status of the reinitialized simulation.
         */
        Status initialize(const std::string& paramFileName);

        /*!
         * Reinitializes the simulation.
         * \param r0 new starting radius
         * \param th0 new starting polar angle
         * \param ph0 new starting azimuthal angle
         * \param ek0 new starting energy
         * \return the status of the reinitialized simulation.
         */
        Status initialize(double r0, double th0, double ph0, double ek0);

        /*!
         * Returns a string containing the trajectory's current coordinates.
         * \return string containing "r,th,ph,ek,s".
         */
        std::string stateToString() const;

        //! Steps the trajectory forward.  The implementation details are completely model-dependent.
        // TODO: are there any model-independent things I can abstract out?  The Parker equation (as in, the
        // calculation of dr, dth, dph, dek), for example?
        virtual Status step() = 0;

        /*!
         * Recalculates and returns thetrajectory's status.
         * \return Sun if r < rSun, Heliopause if r > rHP, Running otherwise.
         */
        virtual Status getStatus();

        /*!
         * Gets the simulation's parameters.
         * \return the simulation's parameter object.
         */
        const T& getParams() const { return params; };

        /*!
         * \sa Parameters::paramString()
         */
        std::string paramString() const { return params.paramString(); };

        /*!
         * Gets path to parameter file.
         */
        const std::string& getParamFileName() const { return params.getParamFileName(); };

    protected:
        //! Simulation's current status.
        Status status;

        //! Contains heliosphere (Vsw, Omega, etc) and particle (charge, mass) parameters.
        T params;

        //! Speed of light in au / s
        static constexpr double speedOfLight = 0.002004;

        //! Backwards time (seconds).
        double s;
        //! Particle's current position
        Point pos;
        //! Kinetic energy (GeV).
        double ek;

        // Useful convenience variables.
        //! Rigidity (GV).
        double P;
        //! Velocity divided by speed of light.
        double beta;
        //! \a tanPsi = \a Omega (\a r - \a rSun) * sin(\a th) / \a Vsw.
        double tanPsi, sinPsi, cosPsi;
        //! \a gamma = \a r * \a Omega * sin(\a th) / \a Vsw.
        double gamma;

        /*!
         * Updates \a P, \a beta, \a tanPsi and \a gamma.
         */
        virtual void updateVars();

        // Functions for updating convenience variables.
        void updateP();
        void updateBeta();
};

template<class T>
TrajectoryBase<T>::TrajectoryBase(const std::string& paramFileName) : params() {
    initialize(paramFileName);
}

template<class T>
Status TrajectoryBase<T>::initialize(double r0, double th0, double ph0, double ek0) {
    // Set initial phase space coordinates
    pos.set(r0, th0, ph0);
    ek = ek0;
    s = 0;

    updateVars();

    // Check status
    return getStatus();
}

template<class T>
Status TrajectoryBase<T>::initialize() {
    // Set initial phase space coordinates
    return initialize(params.getR0(), params.getTh0(), params.getPh0(), params.getEk0());
}

template<class T>
Status TrajectoryBase<T>::initialize(const std::string& paramFileName) {
    // Load new parameters
    params.loadParameters(paramFileName);

    return initialize();
}

template<class T>
Status TrajectoryBase<T>::getStatus() {
    // Check status
    if (pos.getR() < params.getRSun()) {
        status = Status::Sun;
    } else if (pos.getR() > params.getRHP()) {
        status = Status::Heliopause;
    } else {
        status = Status::Running;
    }

    return status;
}

template<class T>
std::string TrajectoryBase<T>::stateToString() const {
    // Use high precision
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::digits10);

    out << static_cast<long double>(pos.getR()) << ","
        << static_cast<long double>(pos.getTh()) << ","
        << static_cast<long double>(pos.getPh()) << ","
        << static_cast<long double>(ek) << ","
        << static_cast<long double>(s);

    return out.str();
}

template<class T>
void TrajectoryBase<T>::updateP() {
    P = sqrt(ek*ek + 2 * ek * params.getM()) / std::abs(params.getCharge());
}

template<class T>
void TrajectoryBase<T>::updateBeta() {
    beta = sqrt(ek*ek + 2 * ek * params.getM()) / (ek + params.getM());
}

template<class T>
void TrajectoryBase<T>::updateVars() {
    updateP();
    updateBeta();
}

#endif


