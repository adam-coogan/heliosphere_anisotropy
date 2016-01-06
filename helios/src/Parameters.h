#ifndef PARAMETERS
#define PARAMETERS 

#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace po = boost::program_options;

/*!
 * Stores parameters relevant to a Helios simulation, like Omega, B0, m, chargeSign, etc.
 */
class Parameters {
    public:
        /*!
         * Initializes parameter descriptions.  A subclass must add new parameter descriptions for each new
         * member variable in its constructor.  Parameter descriptions must have a notifier function of the
         * form
         *      [this](const double& varPar) { this->var = varPar; },
         * where var is the parameter's name.
         */
        Parameters();

        /*!
         * Loads (or reloads) simulation parameters.  Since several parameters have no sensible default, this
         * object is not valid until this method is called.
         * \param paramFile name of file containing parameters relative to the directory containing the
         * Helios binary.
         */
        void loadParameters(const std::string& paramFileName);

        // TODO: check if valid before returning!  Throw exception if not.
        // Should not be overwriteable.
        double getM() const { return m; };
        int getCharge() const { return charge; };
        double getR0() const { return r0; };
        double getTh0() const { return th0; };
        double getPh0() const { return ph0; };
        double getEk0() const { return ek0; };
        double getRHP() const { return rHP; };
        double getRSun() const { return rSun; };
        double getOmega() const { return omega; };
        double getVsw() const { return kmsToaus * Vsw; };

        /*!
         * Creates a string representation of the simulation's parameters.
         * \return a string containing lines of the form "paramName=paramValue".
         */
        virtual std::string paramString() const;

        /*!
         * Gets path to parameter file.
         */
        const std::string& getParamFileName() const { return paramFileName; };

    protected:
        //! Description of parameters to read.
        // Should be a member variable so parameter descriptions can be accessed after the values are read.
        po::options_description desc;

        //! Parameter file
        std::string paramFileName;

        //! Used to assign parameter values to internal variables
        template <typename T>
        void assignParam(T& var, const T& val) { var = val; };

        // Store parameters as member variables to save time on lookup.
        double m;
        int charge;
        double r0;
        double th0;
        double ph0;
        double ek0;
        double rHP;
        double rSun;
        double omega;
        double Vsw;

        //! True if the parameters were successfully read from a file.
        // TODO: use this
        bool valid;

        //! 1 km / s = 6.685e-9 au / s.
        static constexpr double kmsToaus = 6.685e-9;
};

#endif


