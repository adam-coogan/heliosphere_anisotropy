#include "RunConfig.h"
#include <iostream>

// TODO (BIG): make RunConfig and PPTrajectoryParams subclass a Parameter class, with generic methods
// for storing numerical configs in a map.  The requiredNumericalParams vector will specify what to look
// for in the csv.  Subclasses can extend it with a handful of non-numerical configs.

const std::vector<std::string> RunConfig::requiredConfigs = {"ri", "thi", "phii", "ei", "threads_per_energy",
	"runs_per_thread", "output_format", "output_dir"};

std::ostringstream RunConfigException::ss;

RunConfig::RunConfig(const std::string& configFileName0) : configFileName(configFileName0) 
		: ri(-1), thi(-1), phii(-1), threadsPerEnergy(-1), runsPerThread(-1) {
	// Open the file
	std::ifstream configFile(configFileName);

	// Record any extra run configs
	std::vector<std::string> extraConfigs;

	// Make sure file was opened!
	if (configFile.is_open()) {
		// Stores line
		std::string line;
		// String stream for line
		std::stringstream ss;
		// Current config name
		std::string name;
		// Current config value
		std::string val;

		// TODO: throw exception if run configs are multiply defined
		// Read lines from file
		while(std::getline(configFile, line)) {
			// Check if line is a comment
			if (line.size() > 0 && line[0] != '#' && line[0] != '\n' && line[0] != '\r') { 
				// Failbit gets set when ss encounters a blank line, so it needs to be cleared.
				ss.clear();
				// Put next line in string stream
				ss.str(line);

				// Get config name and value.  Ignore anything else in the line.
				std::getline(ss, name, ',');
				std::getline(ss, val, ',');
				// Clear eofbit, which gets set after extracting val, the last thing in ss
				ss.clear();
				// Load val into ss to enable formatted output
				ss.str(val)

				if (std::find(requiredConfigs.begin(), requiredConfigs.end(), name)
						!= requiredConfigs.end()) {
					if (name == "ri") {
						// Check whether ri was already set.
						if (ri == -1) {
							// Check if conversion was successful.
							if (!(ss >> ri)) {
								throw InvalidConfigException(configFileName, name);
							}
						} else {
							throw "multiple defs";
						}
					} else if (name == "thi") {
						if (thi == -1) {
							if (!(ss >> thi)) {
								throw InvalidConfigException(configFileName, name);
							}
						} else {
							throw "multiple defs";
						}
					} else if (name == "phii") {
						if (phii == -1) {
							if (!(ss >> phii)) {
								throw InvalidConfigException(configFileName, name);
							}
						} else {
							throw "multiple defs";
						}
					} else if (name == "threads_per_energy") {
						if (threadsPerEnergy == -1) {
							if (!(ss >> threadsPerEnergy) || threadsPerEnergy <= 0) {
								throw InvalidConfigException(configFileName, name);
							}
						} else {
							throw "multiple defs";
						}
					} else if (name == "runs_per_thread") {
						if (runsPerThread == -1) {
							if (!(ss >> runsPerThread) || runsPerThread <= 0) {
								throw InvalidConfigException(configFileName, name);
							}
						} else {
							throw "multiple defs";
						}
					} else if (name == "output_format") {
						if (outputFormat == "") {
							if (val == "firstlast" || val == "all") {
								outputFormat = val;
							} else {
								throw InvalidConfigException(configFileName, name);
							}
						} else {
							throw "multiple defs";
						}
					} else if (name == "output_dir") {
						if (outputDir == "") {
							if (val == "firstlast" || val == "all") {
								outputFormat = val;
							} else {
								throw InvalidConfigException(configFileName, name);
							}
						} else {
							throw "multiple defs";
						}
					} else if (name == "ei") {
					}
				} else {
					// Unnecessary run config encountered.  Add to the list.
					extraConfigs.push_back(name);
				}
			}
		}

		// Close file when done
		configFile.close();

		// Make sure all configs were read.  There are somewhat cleaner ways to do this...
		std::vector<std::string> missingConfigs;

		for (const std::string& requiredConfig : requiredConfigs) {
			if (configs.count(requiredConfig) == 0) {
				// The config was missing.  Add it to the list.
				if (requiredConfig == "qSign") {
					missingConfigs.push_back("posCharge");
				} else if (requiredConfig == "ac") {
					missingConfigs.push_back("agt0");
				} else {
					missingConfigs.push_back(requiredConfig);
				}
			}
		}

		// If configs are missing or extra ones are present, throw an exception
		if (missingConfigs.size() != 0 || extraConfigs.size() != 0) {
			throw ConfigNotFoundException(configFileName, missingConfigs, extraConfigs);
		}
	} else {
		// Else the config file could not be opened.  Throw an exception!
		throw ConfigFileNotFoundException(configFileName);
	}
}

std::string PPTrajectoryParams::toXML() const {
	// Full precision version of to_string
	std::stringstream converter;
	converter << std::setprecision(std::numeric_limits<double>::digits10);

	std::string xml = "\t<configs>\n";
	
	for (const auto& config : configs) {
		converter.str("");
		converter << config.second;
		xml += "\t\t<" + config.first + ">" + converter.str() + "</" + config.first + ">\n";
	}

	xml += "\t</configs>\n";

	return xml;
}


