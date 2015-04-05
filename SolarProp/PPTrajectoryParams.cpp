#include "PPTrajectoryParams.h"
#include <iostream>

std::ostringstream ParamNotFoundException::ss;

PPTrajectoryParams::PPTrajectoryParams(const std::string& paramFileName0) : paramFileName(paramFileName0) {
	// Open the file
	std::ifstream paramFile(paramFileName);

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

		// Read lines from file
		while(std::getline(paramFile, line)) {
			// Put line in string stream
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
					throw InvalidParamException(name);
				}
			} else if (name == "agt0") {
				if (val == "true") {
					params["ac"] = 1;
				} else if (val == "false") {
					params["ac"] = -1;
				} else {
					// Else invalid CSV!  Throw exception.
					throw InvalidParamException(name);
				}
			} else {
				// Store in parameter map
				try {
					params[name] = std::stod(val);
				} catch (const std::invalid_argument& e) {
					// Thrown if stod fails
					throw InvalidParamException(name);
				}
			}
		}

		// Close file when done
		paramFile.close();

		// Make sure all parameters were read
		std::vector<std::string> missingParams;

		for (const std::string& requiredParam : getRequiredParams()) {
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

		// If parameters are missing, throw an exception
		if (missingParams.size() != 0) {
			throw ParamNotFoundException(paramFileName, missingParams);
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


