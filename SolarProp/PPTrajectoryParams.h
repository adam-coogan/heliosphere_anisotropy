/*
 * Container for all constant numerical PPTrajectory parameters.  Reads values from a CSV file.
 */
#ifndef PPTRAJECTORYPARAMS_H
#define PPTRAJECTORYPARAMS_H 

#include <fstream>
#include <iomanip>
#include <limits> // For converting doubles to strings
#include <map>
#include <sstream>
#include <string>

// TODO: put m, e0 in parameter csv															 2015-02-24 18:19

class PPTrajectoryParams {
	public:
		// No default constructor.
		PPTrajectoryParams(const std::string& paramFileName0);
		~PPTrajectoryParams() { };

		// Only want getters
		double getDs() const { return ds; };
		double getLambda0() const { return lambda0; };
		double getDiffFact() const { return diffFact; };
		double getRig0() const { return rig0; };
		double getB0() const { return b0; };
		double getAc() const { return ac; };
		double getR0() const { return r0; };
		double getVsw() const { return vsw; };
		double getOmega() const { return omega; };
		double getRHP() const { return rHP; };
		double getRSun() const { return rSun; };
		double getM() const { return m; };
		double getQSign() const { return qSign; };
		double getE0() const { return e0; };

		const std::string& getParamFileName() const { return paramFileName; };

		std::string toXML() const;

	private:
		PPTrajectoryParams& assignParams(const std::map<std::string, double>& params);
		// Backwards time step size (s)
		double ds;

		// Solar system parameters
		double lambda0; // Constant in parallel mean free path/kpar (AU)
		double diffFact; // Factor relating parallel and perpendicular diffusion coefficients
		double rig0; // Rigidity constant in kpar (T*AU)
		double b0; // Reference value for B (T)
		double ac; // Polarity of solar cycle (+/- 1)
		double r0; // Distance of Earth from the Sun (AU)
		double vsw; // Solar wind velocity (AU/s)
		double omega; // Sun's equatorial angular velocity (rad/s)
		double rHP; // Distance to heliopause (AU)
		double rSun; // Radius of sun (AU)

		// Particle parameters.  Only allowing electrons and positrons for now.
		double qSign; // Sign of particle charge (+/- 1)
		// TODO: check, since in the paper this is rest energy/nucleon.						 2015-02-23 14:46
		const double m; // Electron mass (GeV)
		const double e0; // Electron rest energy (GeV)

		// Name of parameter file
		const std::string paramFileName;
};

#endif


