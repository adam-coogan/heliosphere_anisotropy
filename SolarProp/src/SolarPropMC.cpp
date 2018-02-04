#include "SolarPropMC.h"

namespace {
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
            int integrationID;

			for (int i = 0; i < rcs.getRunsPerThread(); i++) {
				// Create trajectory starting at point specified in run configuration file
				PPTrajectory traj(rcs.getR(), rcs.getTh(), rcs.getPhi(), ei, paramFileName,
                        rcs.getOutputFormat());

				// Run the full simulation
				traj.integrate(0);

				// Protect access to runCount
                {
                    std::lock_guard<std::mutex> lock(mtx);

                    // Get the intergration number
                    integrationID = runCount++;
                }

                // Write results to XML and increment the run counter
                // TODO: catch exceptions...
                traj.writeToXML(rcs.getOutputDir() + "run_" + std::to_string(integrationID) + ".xml");

				// Protect access to cout
				{
					std::lock_guard<std::mutex> lock(mtx);
					std::cout << "Finished and wrote integration " << integrationID << std::endl;
				}
			}
		} catch (ParamException& e) {
			// Catch exceptions that occur when parameters are loaded
			std::cout << e.what() << std::endl;
		}
	}
}

/*
 * Runs MCs for all initial energies.
 * Arguments:
 *	rcs: the run configurations.
 */
// TODO: how do I handle exceptions in threads???
void SolarPropMC::runMCs(const RunConfig& rcs, const std::string& paramFileName) {
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
			std::get<2>(threadLists.back()).push_back(std::thread(runEnergyMCs, std::ref(ei), std::ref(rcs),
						std::ref(paramFileName), std::ref(std::get<1>(threadLists.back())), std::ref(mtx)));

            // Protect access to std::cout
            {
                std::lock_guard<std::mutex> lock(mtx);
                std::cout << "Added thread " << i << " for energy " << ei << " GeV" << std::endl;
            }
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

