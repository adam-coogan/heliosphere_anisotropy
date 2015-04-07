/*
 * Lightweight container for run parameters.
 */

#ifndef RUNCONFIG_H
#define RUNCONFIG_H

#include "SphericalVector.h"
#include <exception>
#include <iostream>
#include <string>
#include <vector>

class RunConfig {
	public:
		/*
		 * Reads run configuration from a csv file.
		 * Arguments:
		 *	paramFileName0: the csv file containing the run configuration parameters.  If any are missing or
		 *		of the wrong form, a RunConfigException will be thrown.
		 */
		RunConfig(const std::string& runConfigFileName0);

		// Getters
		double getR() const { return start.r; };
		double getTh() const { return start.th; };
		double getPhi() const { return start.phi; };
		const std::vector<double>& getEs() const { return eis; };
		int getThreadsPerEnergy() const { return threadsPerEnergy; };
		int getRunsPerThread() const { return runsPerThread; };

		const std::string& getRunConfigFileName() const { return runConfigFileName; };

		std::string toXML() const;

	private:
		// Internal helper function for defining coordinates
		void assignDouble(double& var, const std::string& name, const std::string& val, bool defined) {
			if (!defined) {
				try {
					var = std::stod(val);
				} catch (const std::invalid_argument& e) {
					throw InvalidRunConfigException(runConfigFileName, name);
				}
			} else {
				throw RunConfigMultipleDefs(runConfigFileName, name);
			}
		};

		void assignInt(int& var, const std::string& name, const std::string& val, bool defined) {
			if (!defined) {
				try {
					var = std::stod(val);
				} catch (const std::invalid_argument& e) {
					throw InvalidRunConfigException(runConfigFileName, name);
				}
			} else {
				throw RunConfigMultipleDefs(runConfigFileName, name);
			}
		};

		// TODO: distinguish between path and file name!!!									 2015-04-06 12:18
		const std::string configFileName;

		// Point at which to begin SDE integration
		double ri, thi, phii;

		// Initial energies to simulate
		std::vector<double> eis;

		// Number of threads to run per energy and number of simulations to carry out in each thread
		int threadsPerEnergy, runsPerThread;

		// Used to check whether required parameters are all specified in the run config file
		const static std::vector<std::string> requiredConfigs;
};

class RunConfigException : public std::runtime_error {
	public:
		RunConfigException(const std::string& runConfigFileName)
			: runtime_error("Run configuration file " + runConfigFileName + ": ") { };

	protected:
		// Useful to have to avoid returning pointers to stack objects...
		static std::ostringstream ss;
};

#endif


