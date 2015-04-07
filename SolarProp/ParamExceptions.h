#ifndef PARAMEXCEPTIONS
#define PARAMEXCEPTIONS 

#include <exception>
#include <sstream>
#include <string>
#include <vector>

class ParamException : public std::runtime_error {
	public:
	   	ParamException(const std::string& paramFileName) : runtime_error("Parameter file " + paramFileName
				+ ": ") { };

	protected:
		// Useful to have to avoid returning pointers to stack objects...
		static std::ostringstream ss;
};

class ParamFileNotFoundException : public ParamException {
	public:
		ParamFileNotFoundException(const std::string& paramFileName) : ParamException(paramFileName) { };

		const char* what() const noexcept {
			ParamException::ss.str("");
			ParamException::ss << ParamException::what() << "could not be opened";

			return ss.str().c_str();
		};
};

class ParamMultipleDefsException : public ParamException {
	public:
		ParamMultipleDefsException(const std::string& paramFileName, const std::string& varName)
		   	: ParamException(paramFileName), name(varName) { };

		const char* what() const noexcept {
			ss.str("");
			ss.clear();

			ss << ParamException::what() << "parameter " << name << " multiply defined";

			return ss.str().c_str();
		};

	private:
		const std::string name;
};

class ParamInvalidException : public ParamException {
	public:
		ParamInvalidException(const std::string& paramFileName, const std::string& varName)
			: ParamException(paramFileName), name(varName) { };

		const char* what() const noexcept {
			ss.str("");
			ss << ParamException::what() << "parameter " << name << " is invalid";

			return ParamException::ss.str().c_str();
		};
	
	private:
		const std::string name;
};

class ParamsNotFoundException : public ParamException {
	public:
		ParamsNotFoundException(const std::string& paramFileName, const std::vector<std::string>& missingPars)
			: ParamException(paramFileName), missingParams(missingPars) { };

		const char* what() const noexcept {
			// Clear ss
			ParamException::ss.str("");
			ParamException::ss << ParamException::what();

			if (missingParams.size() != 0) {
				// Put message and missing parameters into ss
				ParamException::ss << "missing parameter(s): ";

				for (auto mp = missingParams.begin(); mp != missingParams.end() - 1; mp++) {
					ParamException::ss << *mp << ", ";
				}

				ParamException::ss << *(missingParams.end() - 1) << ".  ";
			}


			// Avoid returning a pointer to data inside this exception as it lives on the stack
			return ParamException::ss.str().c_str();
		};

		// Should I return a reference?
		std::vector<std::string> getMissingParams() const { return missingParams; };
	
	private:
		// Any required parameters that were omitted from the parameter file
		const std::vector<std::string> missingParams;
};

#endif


