/*
 * Container for all constant numerical PPTrajectory parameters.  Reads values from a CSV file.  Exceptions
 * that can be thrown when loading parameters are also defined in this header.
 */
#ifndef HELIOPARAMS_H
#define HELIOPARAMS_H

#include "Params.h"
#include <string>
#include <vector>

class HelioParams {
	public:
		/*
		 * Reads heliospheric parameters from a file.
		 *	paramFileName0: the csv file containing the required parameters.  If any are missing or of the
		 *		wrong form, a ParamException will be thrown.
		 *	TODO: important: catch missing parameter exceptions and rethrow!!!
		 */
		HelioParams(const std::string& pFileName) : params(pFileName, {"ds", "lambda0", "rig0", "b0", "r0",
				"vsw", "rHP", "omega", "rSun", "diffFact", "m", "e0", "qSign", "ac"}) {
			// Check whether charge sign is +/- 1
			if (params["qSign"] != 1 && params["qSign"] != -1) {
				throw ParamInvalidException(pFileName, "qSign");
			}

			// Check whether solar polarity is +/- 1
			if (params["ac"] != 1 && params["ac"] != -1) {
				throw ParamInvalidException(pFileName, "ac");
			}
		};

		// Only want getters.  Need to use at since [] doesn't have a const version.
		double getDs() const { return params["ds"]; };
		double getLambda0() const { return params["lambda0"]; };
		double getDiffFact() const { return params["diffFact"]; };
		double getRig0() const { return params["rig0"]; };
		double getB0() const { return params["b0"]; };
		double getR0() const { return params["r0"]; };
		double getVsw() const { return params["vsw"]; };
		double getOmega() const { return params["omega"]; };
		double getRHP() const { return params["rHP"]; };
		double getRSun() const { return params["rSun"]; };
		double getM() const { return params["m"]; };
		double getE0() const { return params["e0"]; };
		double getAc() const { return params["ac"]; };
		double getQSign() const { return params["qSign"]; };

		const std::string& getParamFileName() const { return params.getParamFileName(); };

		std::string toXML(int indents) const { return "\t<params>\n" + params.toXML(indents + 1)
            + "\t</params>\n"; };

	private:
		// Stores all numerical parameters
		Params<double> params;
};

#endif


