#include "SolarPropMC.h"
#include <functional>
#include <iomanip> // TODO: use this to format file names nicely, truncating trailing 0s!!!
#include <iostream>
#include <mutex>
#include <sstream>
#include <thread>
#include <tuple>
#include <vector>

/*
 * Main MC method.  Runs MC samples in separate threads with specified configurations.
 * Arguments:
 *	argv[1]: configuration file for the run.  Specifies number of threads per energy, number of runs per
 *		energy, observed particle energies, output type ("all" vs "exit") and integrator.
 *	argv[2] (optional): file specifying values of the heliospheric parameters.  Default is default_params.csv.
 */
int main(int argc, char* argv[]) {
	if (argc > 1 && argc <= 3) {
		// Assume default parameter
		std::string paramFileName = "default_params.csv";

		// If a non-default parameter file was specified, use it instead
		if (argc == 3) {
			paramFileName = argv[2];
		}

		// Load run configurations
		RunConfig rcs(argv[1]);

		// Run the MCs!
		SolarPropMC::runMCs(rcs, paramFileName);

		return 0;
	} else {
		std::cout << "SolarProp: must supply an argument specifying the pseudoparticle's initial "
			"energy in GeV." << std::endl;

		return 1;
	}
}


