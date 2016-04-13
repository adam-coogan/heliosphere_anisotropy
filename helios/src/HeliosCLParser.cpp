#include "HeliosCLParser.h"

namespace HeliosCLParser {
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
}

