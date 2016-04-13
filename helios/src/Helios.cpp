#include "Basic3D.h"
#include "Basic3DParams.h"
#include "HeliosCLParser.h"
#include "Simulation.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <limits>
#include <string>

int main(int argc, char *argv[]) {
    std::cout.precision(std::numeric_limits<double>::digits10); // Should be max_digits10, but icpc complains

    if (argc >= 2) {
        std::cout << "Parsing command line arguments..." << std::endl;
        // Parse command line arguments...
        po::variables_map clArgs = HeliosCLParser::parseCLArgs(argc, argv);
        // ...and unpack them into variables.
        const int numRuns = clArgs["numruns"].as<int>();
        const std::string configPath = clArgs["configpath"].as<std::string>();
        const std::string outPath = clArgs["outpath"].as<std::string>();

        std::cout << "Building simulation..." << std::endl;

        // Run the simulation
        Simulation<Basic3D<Basic3DParams>> sim(configPath, numRuns);

        sim.run();

        // Write results
        sim.write(outPath);
    } else {
        std::cout << "Must provide output file name." << std::endl;
    }
}


