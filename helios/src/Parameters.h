#ifndef PARAMETERS
#define PARAMETERS 

#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
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
         * Initializes parameter descriptions.  Subclasses must add new parameter descriptions in their 
         * constructor.
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
        double m() const { return params["m"].as<double>(); };
        int charge() const { return params["charge"].as<int>(); };
        double r0() const { return params["r0"].as<double>(); };
        double th0() const { return params["th0"].as<double>(); };
        double ph0() const { return params["ph0"].as<double>(); };
        double ek0() const { return params["ek0"].as<double>(); };
        double rHP() const { return params["rHP"].as<double>(); };
        double rSun() const { return params["rSun"].as<double>(); };
        double omega() const { return params["omega"].as<double>(); };
        double Vsw() const { return kmsToaus * params["Vsw"].as<double>(); };

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
        //! Stores program parameters.
        po::variables_map params;

        //! Description of parameters to read.
        // Should be a member variable so parameter descriptions can be accessed after the values are read.
        po::options_description desc;

        //! True if the parameters were successfully read from a file.
        // TODO: use this
        bool valid;

        //! Parameter file
        std::string paramFileName;

        //! 1 km / s = 6.685e-9 au / s.
        static constexpr double kmsToaus = 6.685e-9;
};

#endif


