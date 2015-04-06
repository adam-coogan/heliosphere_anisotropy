#include "PPTrajectory.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

// TODO: make an executable for running MCs.  Should take number of threads and number of samples as
//	arguments.  Should be able to control what information is written to XML: just initial and exit points?
//	The whole trajectory?																	 2015-04-03 23:30
// TODO: make parameter file an argument, but supply a default.								 2015-04-03 23:32
// TODO: threading stuff.  For now, just take energy as first argument.						 2015-04-05 12:46
//	-Test this!

/*
 * Runs a set of samples.
 * Arguments:
 *	numSims: number of MC samples to run.
 *	paramFile: heliosphere parameter configuration file.
 *	TODO: pass run configuration file!  Should figure out a good way of naming the XML files or something.
 */
void runMCs(int numSims, const std::string& paramFile, double ei) {
	try {
		std::cout << "Generating " << numSims << " trajectories for particle with measured energy "
			<< ei << " GeV..." << std::endl;

		for (int i = 0; i < numSims; i++) {
			// Particles start at Earth
			// TODO: Earth's position should probably not be hardcoded in like this.  Should have a special
			// set of constructors for starting trajectories at and around Earth!			 2015-04-05 15:22
			PPTrajectory traj(1, PPPoint::pi / 2.0, 0, ei, paramFile);

			// Run the full simulation
			traj.integrate(0);

			// Write results to XML
			traj.writeToXML("runs/" + std::to_string(ei) + "_GeV_run_" + std::to_string(i) + ".xml");
		}

		std::cout << "Done." << std::endl;
	} catch (ParamException& e) {
		// Catch exceptions that occur when parameters are loaded
		std::cout << e.what() << std::endl;

		// Rethrow e
		//throw;
	}
}

void testFun(){
	std::cout << "testFun" << std::endl;
}

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

		// Create a vector of threads, adding the appropriate number per energy
		std::vector<double> energies = {std::stod(argv[1])};
		std::vector<std::vector<std::thread> > threadLists(energies.size());

		int runsPerThread = 1;

		// Start the threads
		for (auto& threadList : threadLists) {
			for (int i = 0; i < runsPerThread; i++) {
				threadList.push_back(std::thread(runMCs, 1, paramFile, 100));
			}
		}

		// Block until all threads finish
		for (auto& threadList : threadLists) {
			for (auto& thread : threadList) {
				thread.join();
			}
		}

		return 0;
	} else {
		std::cout << "SolarProp: must supply an argument specifying the pseudoparticle's initial "
			"energy in GeV." << std::endl;

		return 1;
	}
}


