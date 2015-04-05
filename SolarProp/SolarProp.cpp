#include "PPTrajectory.h"

#include <iostream>
#include <sstream>
#include <iomanip>

// TODO: make an executable for running MCs.  Should take number of threads and number of samples as
//	arguments.  Should be able to control what information is written to XML: just initial and exit points?
//	The whole trajectory?																	 2015-04-03 23:30
// TODO: make parameter file an argument, but supply a default.								 2015-04-03 23:32
// TODO: threading stuff.  For now, just take energy as first argument.						 2015-04-05 12:46
// TODO: sample from points around Earth!													 2015-04-05 12:49

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
		std::string paramFile = "default_params.csv";

		// If a non-default parameter file was specified, use it instead
		if (argc == 3) {
			paramFile = argv[2];
		}

		// Try constructing the trajectory.  Many problems can occur when loading parameters.
		try {
			// Particles start at Earth
			PPTrajectory traj(1, PPPoint::pi / 2.0, 0, std::stod(argv[1]), paramFile);

			std::cout << "Generating trajectory for particle with measured energy " << traj.getE0()
				<< " GeV..." << std::endl;

			// Run the full simulation
			traj.integrate(0);

			std::cout << "Done." << std::endl;

			// Write results to XML
			traj.writeToXML("runs/" + std::string(argv[1]) + "_GeV.xml");

			return 0;
		} catch (std::runtime_error& e) {
			// Catch exceptions that occur when parameters are loaded
			std::cout << e.what() << std::endl;
		}
	} else {
		std::cout << "SolarProp error: must supply an argument specifying the pseudoparticle's initial "
			"energy in GeV." << std::endl;

		return 1;
	}
}


