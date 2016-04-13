#ifndef HELIOSCLPARSER
#define HELIOSCLPARSER 

#include <boost/program_options.hpp>
#include <string>

namespace HeliosCLParser {
    namespace po = boost::program_options;

    /*!
     * Parses command line arguments.
     * \return an object containing the number of runs, configuration path and output file path for a simulation.
     */
    boost::program_options::variables_map parseCLArgs(int argc, char *argv[]);
}

#endif

