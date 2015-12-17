#ifndef SIMULATION
#define SIMULATION

#include "TrajectoryBase.h"
#include <fstream>
#include <iostream>

/*!
 * Manages a Helios simulation.
 * \param T the trajectory class to use for the simulation.
 */
template<class T>
class Simulation {
    public:
        /*!
         * Constructor.
         * \param configPath path to simulation configuration file.
         * \param numRuns number of trajectories to generate.
         * \param verbose if true, status messages will be printed throughout the simulation.  True by
         * default.
         */
        Simulation(const std::string& configPath, int numRuns, bool verbose = true);

        //void reset(const std::string& configPath, int numRuns, bool verbose = true);

        /*!
         * Writes data generated during simulation.  Data is written in csv format.
         * \param file in which to write simulation data.  Should have .csv extension.
         */
        void write(const std::string& outPath) const;

        /*!
         * Carries out a simulation of \a numRuns trajectories from Earth to the heliopause.
         */
        void run();

        /*!
         * \return number of trajectories that terminate in Jupiter's magnetosphere.
         */
        int getJupiterCount() const { return jupiterCount; };

        /*!
         * \return number of trajectories that terminate in the sun.
         */
        int getSunCount() const { return sunCount; };

    private:
        T traj;
        int numRuns;
        bool verbose;

        //! String containing header data and phase space endpoints for each run in the simulation.
        std::string runsString;
        int jupiterCount, sunCount;
        //! Time between beginning and end of the simulation.
        double duration;
};

template<class T>
Simulation<T>::Simulation(const std::string& configPath, int numRuns, bool verbose) : traj(configPath),
        numRuns(numRuns), verbose(verbose) {
    // String containing run data
    runsString = "# Run exit points.  Columns are r (AU), th (rad), ph (rad), ek (GeV), s (s).\n# First line"
        " contains initial point of trajectory.  Generated with simulation settings from "
        + traj.getParamFileName() + "\n" + traj.stateToString();
}

template<class T>
void Simulation<T>::run() {
    if (verbose) {
        std::cout << "Tracing " << numRuns << " particles detected with energy "
            << traj.getParams().getParamFileName() << " GeV at Earth back to the heliopause..." << std::endl;
    }

    // Variable for tracking what percent of the simulation has finished
    int percentDone = 0;
    int runsTo1Percent = numRuns / 100;

    // Measure how long the simulation takes.  Store time since program started.
    std::clock_t start = std::clock();

    // Variable for storing current simulation's status
    Status stepStatus;

    // Generate the runs
    for (int successes = 0; successes < numRuns; ) {
        while (true) {
#if DEBUG
            std::cout << traj.stateToString() << std::endl;
#endif
            // Step the simulation
            stepStatus = traj.step();

#if DEBUG
            std::cout << std::endl;
#endif

            if (stepStatus == Status::Jupiter) {
                // Trajectory ended in Jupiter.  Run another trajectory and don't increment i.
                jupiterCount++;
                break;
            } else if (stepStatus == Status::Sun) {
                // Trajectory ended in the sun.  Run another trajectory and don't increment i.
                sunCount++;
                break;
            } else if (stepStatus == Status::Heliopause) {
                // Particle ended at the heliopause!  Increment the event counter, add the end point to
                // the list and run another trajectory.
                successes++;
                runsString += "\n" + traj.stateToString();

                // Print percent indicator
                runsTo1Percent--;

                if (runsTo1Percent <= 0) {
                    percentDone += 1;
                    runsTo1Percent = numRuns / 100;

                    if (verbose) {
                        std::cout << percentDone << "% of runs complete" << std::endl;
                    }
                }

                break;
            }
        }

        // Reinitialize simulation variables
        traj.initialize();
    }

    // Measure new time since program started, subtract and convert to seconds or minutes
    double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    if (verbose) {
        if (duration < 60) {
            std::cout << "Simulation completed in " << duration << " seconds" << std::endl;
        } else {
            std::cout << "Simulation completed in " << duration / 60 << " minutes" << std::endl;
        }
    }
}

template<class T>
void Simulation<T>::write(const std::string& outPath) const {
    // Write run data to a CSV
    std::ofstream writer(outPath);

    if (writer.is_open()) {
        writer << runsString;
        writer.close();
    }

    if (verbose) {
        std::cout << "Wrote exit point data to " << outPath << std::endl;
        // Beep!
        std::cout << '\a';
    }
}

#endif


