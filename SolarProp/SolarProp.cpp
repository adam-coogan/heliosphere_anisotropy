#include "PPTrajectory.h"
#include "RunConfig.h"

#include <functional>
#include <iomanip> // TODO: use this to format file names nicely, truncating trailing 0s!!!
#include <iostream>
#include <mutex>
#include <sstream>
#include <thread>
#include <tuple>
#include <vector>

// TODO: make an executable for running MCs.  Should take number of threads and number of samples as
//	arguments.  Should be able to control what information is written to XML: just initial and exit points?
//	The whole trajectory?																	 2015-04-03 23:30
// TODO: make this object oriented, or maybe put the functions inside a namespace???  (Can hide runEnergyMCs
// using an unnamed namespace!)
/*
 * Runs a set of MCs for a thread for the given initial energy and start point.
 * Arguments:
 *	ei: the energy for the run.
 *	rcs: the run configuration.  This contains the number of MCs to run in this thread.
 *	paramFileName: heliosphere parameter configuration file.
 */
void runEnergyMCs(const double& ei, const RunConfig& rcs, const std::string& paramFileName, int& runCount,
		std::mutex& mtx) {
	try {
		for (int i = 0; i < rcs.getRunsPerThread(); i++) {
			// Create trajectory starting at point specified in run configuration file
			PPTrajectory traj(rcs.getR(), rcs.getTh(), rcs.getPhi(), ei, paramFileName);

			// Protect access to std::cout
			mtx.lock();
			std::cout << "Running an integration..." << std::endl;
			mtx.unlock();

			// Run the full simulation
			traj.integrate(0);

			// Protect access to runCount with the mutex
			mtx.lock();
			// Write results to XML and increment the run counter
			traj.writeToXML(rcs.getOutputDir() + std::to_string(ei) + "_GeV_run_" + std::to_string(runCount++)
					+ ".xml");
			mtx.unlock();
		}
	} catch (ParamException& e) {
		// Catch exceptions that occur when parameters are loaded
		std::cout << e.what() << std::endl;
	}
}

/*
 * Runs MCs for all initial energies.
 * Arguments:
 *	rcs: the run configurations.
 */
// TODO: how do I handle exceptions in threads???
void runMCs(const RunConfig& rcs, const std::string& paramFileName) {
	// Create a vector of pairs of energies, the number of runs that have been completed for that energy,
	// and the energies' associated threads
	std::vector<std::tuple<double, int, std::vector<std::thread>>> threadLists;
	// Used to protect run counters
	std::mutex mtx;

	for (const double& ei : rcs.getEis()) {
		// Create the collection of threads and counters
		threadLists.push_back(std::make_tuple(ei, 0, std::vector<std::thread>()));

		std::cout << "Generating " << rcs.getThreadsPerEnergy() << " threads each with "
			<< rcs.getRunsPerThread() << " trajectories for particle with measured energy " << ei
			<< " GeV..." << std::endl;

		// Add the correct numer of threads for each energy
		for (int i = 0; i < rcs.getThreadsPerEnergy(); i++) {
			// Add a new thread executing runEnergyMCs to the back.  Pass it a reference to the run counter
			// for the energy.  References must be wrapped in std::ref.
			//mtx.lock();
			std::cout << "Added thread " << i << " for energy " << ei << std::endl;
			//mtx.unlock();

			std::get<2>(threadLists.back()).push_back(std::thread(runEnergyMCs, std::ref(ei), std::ref(rcs),
						std::ref(paramFileName), std::ref(std::get<1>(threadLists.back())), std::ref(mtx)));
		}
	}

	// Block until all threads finish.  Loop over the pairs of energies and their threads...
	for (auto& threadList : threadLists) {
		// ...and then over the threads!
		for (auto& thread : std::get<2>(threadList)) {
			thread.join();
		}
	}

	std::cout << "All simulations complete!" << std::endl;
}

/*
 * Main MC method.  Runs MC samples in separate threads with specified configurations.
 * Arguments:
 *	argv[1]: configuration file for the run.  Specifies number of threads per energy, number of runs per
 *		energy, observed particle energies, output type ("all" vs "exit") and integrator.
 *	argv[2] (optional): file specifying values of the heliospheric parameters.  Default is default_params.csv.
 */
int main(int argc, char* argv[]) {
	// TODO: by default, start at Earth and take an energy as a command line arg.
	// TODO: important: use boost::program_options!
	if (argc > 1 && argc <= 3) {
		// Assume default parameter
		std::string paramFileName = "default_params.csv";

		// If a non-default parameter file was specified, use it instead
		if (argc == 3) {
			paramFileName = argv[2];
		}

		// Load run configurations
		RunConfig rcs(argv[1]);

		// FAILS: run the MCs!
		runMCs(rcs, paramFileName);

		return 0;
	} else {
		std::cout << "SolarProp: must supply an argument specifying the pseudoparticle's initial "
			"energy in GeV." << std::endl;

		return 1;
	}
}


