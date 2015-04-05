/*
 * Container for all constant numerical PPTrajectory parameters.  Reads values from a CSV file.  Exceptions
 * that can be thrown when loading parameters are also defined in this header.
 */
#ifndef PPTRAJECTORYPARAMS_H
#define PPTRAJECTORYPARAMS_H 

#include <fstream>
#include <iomanip>
#include <limits> // For converting doubles to strings
#include <map>
#include <sstream>
#include <string>
#include <vector>

// TODO: make the exceptions in this file subclass some virtual class.						 2015-04-05 13:48

class PPTrajectoryParams {
	public:
		// No default constructor.
		PPTrajectoryParams(const std::string& paramFileName0);
		~PPTrajectoryParams() { };

		// These parameters must be specified in the parameter file.
		static std::vector<std::string> getRequiredParams() {
			return {"ds", "lambda0", "rig0", "b0", "r0", "vsw", "rHP", "omega", "rSun", "diffFact", "qSign",
				"m", "e0", "ac"};
		};

		// Only want getters.  Need to use at since [] doesn't have a const version.
		double getDs() const { return params.at("ds"); };
		double getLambda0() const { return params.at("lambda0"); };
		double getDiffFact() const { return params.at("diffFact"); };
		double getRig0() const { return params.at("rig0"); };
		double getB0() const { return params.at("b0"); };
		double getAc() const { return params.at("ac"); };
		double getR0() const { return params.at("r0"); };
		double getVsw() const { return params.at("vsw"); };
		double getOmega() const { return params.at("omega"); };
		double getRHP() const { return params.at("rHP"); };
		double getRSun() const { return params.at("rSun"); };
		double getM() const { return params.at("m"); };
		double getQSign() const { return params.at("qSign"); };
		double getE0() const { return params.at("e0"); };

		const std::string& getParamFileName() const { return paramFileName; };

		std::string toXML() const;

	private:
		std::map<std::string, double> params;

		// Name of parameter file
		const std::string paramFileName;
};

class ParamFileNotFoundException : public std::runtime_error {
	public:
		ParamFileNotFoundException(const std::string& paramFile) : runtime_error("SolarProp: parameter file "
			 + paramFile + " could not be opened.") { };
};

class ParamNotFoundException : public std::runtime_error {
	public:
		ParamNotFoundException(const std::string& paramFile, const std::vector<std::string>& missingPars)
			: runtime_error("SolarProp: parameter file " + paramFile + " is missing a parameter: "),
			missingParams(missingPars){ };

		const char* what() const noexcept {
			// Clear ss
			ss.str("");

			// Put message and missing parameters into ss
			ss << runtime_error::what();

			for (auto mp = missingParams.begin(); mp != missingParams.end() - 1; mp++) {
				ss << *mp << ", ";
			}

			ss << *(missingParams.end() - 1);

			// Avoid returning a pointer to data inside this exception as it lives on the stack
			return ss.str().c_str();
		};
	
	private:
		const std::vector<std::string>& missingParams;

		static std::ostringstream ss;
};

class InvalidParamException : public std::runtime_error {
	public:
		InvalidParamException(const std::string& param) : runtime_error("SolarProp: parameter " + param
			+ " is invalid.") { };
};

#endif


