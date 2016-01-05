#include "Basic3D.h"
#include "Simulation.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <limits>
#include <string>

namespace po = boost::program_options;

po::variables_map parseCLArgs(int argc, char *argv[]) {
    // Description of command line options
    po::options_description clDesc;
    clDesc.add_options()
        ("configpath,c", po::value<std::string>()->default_value("config/params.config", "configuration file"))
        ("outpath,o", po::value<std::string>(), "output file path (including csv extension)")
        ("numruns,n", po::value<int>()->default_value(1), "number of runs")
    ;

    // Parse command line arguments
    po::variables_map clArgs;
    po::store(po::parse_command_line(argc, argv, clDesc), clArgs);
    po::notify(clArgs);

    return clArgs;
}

int main(int argc, char *argv[]) {
    std::cout.precision(std::numeric_limits<double>::digits10); // Should be max_digits10, but icpc complains

    if (argc >= 2) {
        std::cout << "Parsing command line arguments..." << std::endl;
        // Parse command line arguments...
        po::variables_map clArgs = parseCLArgs(argc, argv);
        // ...and unpack them into variables.
        const int numRuns = clArgs["numruns"].as<int>();
        const std::string configPath = clArgs["configpath"].as<std::string>();
        const std::string outPath = clArgs["outpath"].as<std::string>();

        std::cout << "Building simulation..." << std::endl;

        // Run the simulation
        Simulation<Basic3D> sim(configPath, numRuns);

        sim.run();

        // Write results
        sim.write(outPath);
    } else {
        std::cout << "Must provide output file name." << std::endl;
    }
}


