#include "PPTrajectoryParams.h"
#include <iostream>

// TODO: exceptions																			 2015-02-24 14:28
// TODO: put in default parameter values													 2015-02-24 14:30
PPTrajectoryParams::PPTrajectoryParams(const std::string& paramFileName0) : m(0.000511), e0(m),
	paramFileName(paramFileName0) {
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
		// Map for storing results
		std::map<std::string, double> params;

		// Read lines from file
		while(std::getline(paramFile, line)) {
			// Put line in string stream
			lineSS.str(line);

			// Get parameter name and value.  Ignore anything else in the line.
			std::getline(lineSS, name, ',');
			std::getline(lineSS, val, ',');

			// See if we encountered one of the boolean parameters
			if (name == "electron") {
				if (val == "true") {
					qSign = -1;
				} else if (val == "false") {
					qSign = 1;
				} // Else invalid CSV!  Throw exception.
			} else if (name == "agt0") {
				if (val == "true") {
					ac = +1;
				} else if (val == "false") {
					ac = -1;
				} // Else invalid CSV!  Throw exception.
			} else {
				// Store in parameter map
				// TODO: catch this exception!												 2015-02-24 17:23
				params[name] = std::stod(val);
			}
		}

		// Close file when done
		paramFile.close();

		// Assign parameters
		assignParams(params);
	} // Else the parameter file could not be opened.  
	// TODO: throw an exception and invalidate the PPTrajectory!							 2015-02-24 12:04
}

// TODO: throw an exception if a parameter cannot be found!									 2015-02-24 12:05
PPTrajectoryParams& PPTrajectoryParams::assignParams(const std::map<std::string, double>& params) {
	// Look up values in the map
	ds = params.at("ds");

	lambda0 = params.at("lambda0");
	diffFact = params.at("diffFact");
	rig0 = params.at("rig0");
	b0 = params.at("b0");
	r0 = params.at("r0");
	vsw = params.at("vsw");
	omega = params.at("omega");
	rHP = params.at("rHP");
	rSun = params.at("rSun");

	return *this;
}

std::string PPTrajectoryParams::toXML() const {
	// Full precision version of to_string
	std::stringstream converter;
	converter << std::setprecision(std::numeric_limits<double>::digits10);

	// This is fucking stupid.
	std::string xml = "\t<params>\n";
	converter << ds;
	xml += "\t\t<ds>" + converter.str() + "</ds>\n";

	converter.str("");
	converter << lambda0;
	xml += "\t\t<lambda0>" + converter.str() + "</lambda0>\n";

	converter.str("");
	converter << diffFact;
	xml += "\t\t<diffFact>" + converter.str() + "</diffFact>\n";

	converter.str("");
	converter << rig0;
	xml += "\t\t<rig0>" + converter.str() + "</rig0>\n";

	converter.str("");
	converter << b0;
	xml += "\t\t<b0>" + converter.str() + "</b0>\n";

	converter.str("");
	converter << ac;
	xml += "\t\t<ac>" + converter.str() + "</ac>\n";

	converter.str("");
	converter << r0;
	xml += "\t\t<r0>" + converter.str() + "</r0>\n";

	converter.str("");
	converter << vsw;
	xml += "\t\t<vsw>" + converter.str() + "</vsw>\n";

	converter.str("");
	converter << omega;
	xml += "\t\t<omega>" + converter.str() + "</omega>\n";

	converter.str("");
	converter << rHP;
	xml += "\t\t<rHP>" + converter.str() + "</rHP>\n";

	converter.str("");
	converter << rSun;
	xml += "\t\t<rSun>" + converter.str() + "</rSun>\n";

	converter.str("");
	converter << m;
	xml += "\t\t<m>" + converter.str() + "</m>\n";

	converter.str("");
	converter << e0;
	xml += "\t\t<e0>" + converter.str() + "</e0>\n";

	converter.str("");
	converter << qSign;
	xml += "\t\t<qSign>" + converter.str() + "</qSign>\n\t</params>\n";

	return xml;
}


