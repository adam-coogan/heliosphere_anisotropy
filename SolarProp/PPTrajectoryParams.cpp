#include "PPTrajectoryParams.h"
#include <iostream>

const std::vector<std::string> PPTrajectoryParams::requiredParams = {"ds", "lambda0", "rig0", "b0", "r0", "vsw", "rHP", "omega",
	"rSun", "diffFact", "qSign", "m", "e0", "ac"};

std::ostringstream ParamException::ss;

PPTrajectoryParams::PPTrajectoryParams(const std::string& paramFileName0) : paramFileName(paramFileName0) {
	// Open the file
	std::ifstream paramFile(paramFileName);

	// Record any extra parameters
	std::vector<std::string> extraParams;

	// Make sure file was opened!
	if (paramFile.is_open()) {
		// Stores line
		std::string line;
		// String stream for line
		std::istringstream lineSS;
		// Current parameter name
		std::string name;
		// Current parameter value
		std::string val;

		// TODO: throw exception if parameters are multiply defined
		// Read lines from file
		while(std::getline(paramFile, line)) {
			// Check if line is a comment
			if (line.size() > 0 && line[0] != '#' && line[0] != '\n' && line[0] != '\r') { 
				// Failbit gets set when lineSS encounters a blank line, so it needs to be cleared.
				lineSS.clear();
				// Put next line in string stream
				lineSS.str(line);

				// Get parameter name and value.  Ignore anything else in the line.
				std::getline(lineSS, name, ',');
				std::getline(lineSS, val, ',');

				// See if we encountered one of the boolean parameters
				if (name == "posCharge") {
					if (val == "true") {
						params["qSign"] = 1;
					} else if (val == "false") {
						params["qSign"] = -1;
					} else {
						// Else invalid CSV!  Throw exception.
						throw InvalidParamException(paramFileName, name);
					}
				} else if (name == "agt0") {
					if (val == "true") {
						params["ac"] = 1;
					} else if (val == "false") {
						params["ac"] = -1;
					} else {
						// Else invalid CSV!  Throw exception.
						throw InvalidParamException(paramFileName, name);
					}
				} else if (std::find(requiredParams.begin(), requiredParams.end(), name)
						!= requiredParams.end() && name != "ac" && name != "qSign") {
					// Check whether name is a required non-boolean parameter.
					try {
						params[name] = std::stod(val);
					} catch (const std::invalid_argument& e) {
						// Thrown if stod fails
						throw InvalidParamException(paramFileName, name);
					}
				} else {
					// Unnecessary parameter encountered.  Add to the list.
					extraParams.push_back(name);
				}
			}
		}

		// Close file when done
		paramFile.close();

		// Make sure all parameters were read.  There are somewhat cleaner ways to do this...
		std::vector<std::string> missingParams;

		for (const std::string& requiredParam : requiredParams) {
			if (params.count(requiredParam) == 0) {
				// The parameter was missing.  Add it to the list.  In the case of the boolean parameters,
				// refer to the name to be used in the parameter file rather than the numerical value.
				if (requiredParam == "qSign") {
					missingParams.push_back("posCharge");
				} else if (requiredParam == "ac") {
					missingParams.push_back("agt0");
				} else {
					missingParams.push_back(requiredParam);
				}
			}
		}

		// If parameters are missing or extra ones are present, throw an exception
		if (missingParams.size() != 0 || extraParams.size() != 0) {
			throw ParamNotFoundException(paramFileName, missingParams, extraParams);
		}
	} else {
		// Else the parameter file could not be opened.  Throw an exception!
		throw ParamFileNotFoundException(paramFileName);
	}
}

std::string PPTrajectoryParams::toXML() const {
	// Full precision version of to_string
	std::stringstream converter;
	converter << std::setprecision(std::numeric_limits<double>::digits10);

	std::string xml = "\t<params>\n";
	
	for (const auto& param : params) {
		converter.str("");
		converter << param.second;
		xml += "\t\t<" + param.first + ">" + converter.str() + "</" + param.first + ">\n";
	}

	xml += "\t</params>\n";

	return xml;
}


