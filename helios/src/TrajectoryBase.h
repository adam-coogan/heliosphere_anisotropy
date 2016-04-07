#ifndef TRAJECTORYBASE
#define TRAJECTORYBASE

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

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
        //! Distance from center of sun (au).
        double r;
        //! Polar angle (radians).
        double th;
        //! Azimuthal angle (radians).
        double ph;
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

        /*!
         * Recalculates angles to make sure th is in [0, pi] and ph is in [0, 2pi].
         */
        void renormalize();

        // Functions for updating convenience variables.
        void updateP();
        void updateBeta();
};

template<class T>
TrajectoryBase<T>::TrajectoryBase(const std::string& paramFileName) : params() {
    initialize(paramFileName);
}

template<class T>
Status TrajectoryBase<T>::initialize() {
    // Set initial phase space coordinates
    r = params.getR0();
    th = params.getTh0();
    ph = params.getPh0();
    ek = params.getEk0();
    s = 0;

    updateVars();

    // Check status
    return getStatus();
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
    if (r < params.getRSun()) {
        status = Status::Sun;
    } else if (r > params.getRHP()) {
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

    out << static_cast<long double>(r) << ","
        << static_cast<long double>(th) << ","
        << static_cast<long double>(ph) << ","
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

template<class T>
void TrajectoryBase<T>::renormalize() {
    while (th > M_PI) {
        th = 2 * M_PI - th;
        ph = ph - M_PI;
    }

    while (th < 0) {
        th = -th;
        ph += M_PI;
    }

    while (ph > 2 * M_PI) {
        ph -= 2 * M_PI;
    }

    while (ph < 0) {
        ph += 2 * M_PI;
    }
}

#endif


