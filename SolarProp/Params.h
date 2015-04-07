#ifndef PARAMS_H
#define PARAMS_H 

#include "ParamExceptions.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits> // For converting doubles to strings
#include <map>
#include <sstream>
#include <string>
#include <vector>

template <class T>
class Params {
	public:
		/*
		 * Reads parameters from a csv file.  Blank lines and lines beginning with # are ignored, as are
		 * lines containing definitions of parameters not in this Params' requiredParams.
		 *	Arguments:
		 *		paramFileName: the csv file containing the required parameters.  If any are missing or of the
		 *			wrong form, a ParamException will be thrown.
		 *		requiredPs: a vector containing the names of the parameters that are required to be defined
		 *			in paramFileName.
		 */
		Params(const std::string& pFileName, const std::vector<std::string>& requiredPs)
				: paramFileName(pFileName), requiredParams(requiredPs) {
			// Open the file
			std::ifstream paramFile(paramFileName);

			// Make sure file was opened!
			if (paramFile.is_open()) {
				// Stores line
				std::string line;
				// String stream for line
				std::stringstream ss;
				// Current parameter name
				std::string name;
				// Current parameter value
				std::string val;

				// Read lines from file
				while(std::getline(paramFile, line)) {
					// Check if line is a comment or newline of some sort
					if (line.size() > 0 && line[0] != '#' && line[0] != '\n' && line[0] != '\r') { 
						// Failbit gets set when ss encounters a blank line, so it needs to be cleared.
						ss.clear();
						// Put next line in string stream
						ss.str(line);

						// Get parameter name and value.  Ignore anything else in the line.
						std::getline(ss, name, ',');
						std::getline(ss, val, ',');

						// Clear eofbit, which gets set after extracting val, the last thing in ss
						ss.clear();
						// Stick val into ss so it can be extracted in the correct format
						ss.str(val);

						// See if we encountered a required parameter.  Everything else in the file is
						// ignored.
						if (std::find(requiredParams.begin(), requiredParams.end(), name)
								!= requiredParams.end()) {
							// Check whether it's in the map already
							if (params.count(name) == 0) {
								// Make sure parameter was converted correctly.
								if (!(ss >> params[name])) {
									// If it was not, it must be of the wrong type
									throw ParamInvalidException(paramFileName, name);
								}
							} else {
								// If it's in the map, throw an exception
								throw ParamMultipleDefsException(paramFileName, name);
							}
						}
					}
				}

				// Close file when done
				paramFile.close();

				// Check whether all required parameters are in the parameter map
				std::vector<std::string> missingParams;

				for (const auto& requiredParam : requiredParams) {
					if (params.count(requiredParam) == 0) {
						// The parameter was missing.  Add it to the list. 
						missingParams.push_back(requiredParam);
					}
				}

				// If parameters are missing or extra ones are present, throw an exception
				if (missingParams.size() != 0) {
					throw ParamsNotFoundException(paramFileName, missingParams);
				}
			} else {
				// Else the parameter file could not be opened.  Throw an exception!
				throw ParamFileNotFoundException(paramFileName);
			}
		};

		/*
		 * Accesses elements in the parameter set.
		 *	Arguments:
		 *		name: the name of a parameter to look up the value of.
		 *	Returns: a reference to the parameter with the given name.
		 */
		T& operator[](const std::string& name) { return params[name]; };

		/*
		 * Accesses elements in the parameter set.
		 *	Arguments:
		 *		name: the name of a parameter to look up the value of.
		 *	Returns: a const reference to the parameter with the given name.
		 */
		const T& operator[](const std::string& name) const { return params.at(name); };

		const std::string& getParamFileName() const { return paramFileName; };

		/*
		 * Gives a string listing the values of the parameters in this object.
		 * Returns:
		 *	A string containing lines of the form <name>value</name> for each parameter.
		 */
		std::string toXML() const {
			// Full precision version of to_string
			std::stringstream converter;
			converter << std::setprecision(std::numeric_limits<double>::digits10);

			std::string xml;
			
			for (const auto& param : params) {
				converter.str("");
				converter << param.second;
				xml += "\t\t<" + param.first + ">" + converter.str() + "</" + param.first + ">\n";
			}

			return xml;
		};

	private:
		const std::string paramFileName;

		std::vector<std::string> requiredParams;

		std::map<std::string, T> params;
};

#endif


