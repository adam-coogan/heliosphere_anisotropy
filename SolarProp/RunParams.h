/*
 * Lightweight container for run parameters.
 */

#ifndef RUNPARAMS_H
#define RUNPARAMS_H

#include "Params.h"
#include <string>
#include <vector>

class RunConfig {
	public:
		/*
		 * Reads run configuration from a csv file.
		 * Arguments:
		 *	rcFileName: the csv file containing the run configuration parameters.  If any are missing or
		 *		of the wrong form, a RunConfigException will be thrown.  
		 *		Note that ri can be < 0, thi can lie outside of [0, pi) and phii can lie outside of [0, 2pi).
		 *		Such coordinates will be fixed when the trajectory is constructed.
		 */
		RunConfig(const std::string& rcFileName) : startPoint(rcFileName, {"ri", "thi", "phii"}),
				intParams(rcFileName, {"threads_per_energy", "runs_per_thread"}),
				flags(rcFileName, {"output_format", "output_dir"}) {
			// Must have at least one thread per energy
			if (intParams["threads_per_energy"] <= 0) {
				throw ParamInvalidException(rcFileName, "threads_per_energy");
			}

			// Must run at least one MC per thread
			if (intParams["runs_per_thread"] <= 0) {
				throw ParamInvalidException(rcFileName, "runs_per_thread");
			}

			// The output format must specify whether to record all points or just the first and last
			if (flags["output_format"] != "all" && flags["output_format"] != "firstlast") {
				throw ParamInvalidException(rcFileName, "runs_per_thread");
			}

			// If all of these parameters are valid, read the energies
		};

		// Getters
		double getR() const { return startPoint["ri"]; };
		double getTh() const { return startPoint["thi"]; };
		double getPhi() const { return startPoint["phii"]; };
		int getThreadsPerEnergy() const { return intParams["threads_per_energy"]; };
		int getRunsPerThread() const { return intParams["runs_per_thread"]; };
		const std::vector<double>& getEs() const { return eis; };

		// All parameters come from the same file
		const std::string& getParamFileName() const { return startPoint.getParamFileName(); };

		// TODO: write ei to XML!
		std::string toXML() const {
			return "\t<run_configs>\n" + startPoint.toXML() + intParams.toXML() + flags.toXML()
				+ "\t</params>\n";
		};

	private:
		// Initial energies to simulate
		std::vector<double> eis;

		// Contains coordinates at which to begin SDE integration
		Params<double> startPoint;

		// Contains number of threads to run per energy and number of simulations to carry out in each thread
		Params<int> intParams;

		// Contains output format and directory
		Params<std::string> flags;
};

#endif


