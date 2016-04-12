#include "Basic3D.h"
#include "Basic3DParams.h"
#include "Wavy3D.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <limits>
#include <string>

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    std::cout.precision(std::numeric_limits<double>::digits10); // Should be max_digits10, but icpc complains

    // Description of command line options
    po::options_description clDesc;
    clDesc.add_options()
        ("configpath,c", po::value<std::string>()->default_value("/Users/acoogan/Dropbox/"
             "heliosphere_anisotropy/helios/config/params.config", "configuration file"))
        ("outpath,o", po::value<std::string>(), "output file") // TODO: make this
    ;

    // Parse command line arguments
    po::variables_map clParams;
    po::store(po::parse_command_line(argc, argv, clDesc), clParams);
    po::notify(clParams);

    // Make a test trajectory
    Basic3D<Basic3DParams> traj(clParams["configpath"].as<std::string>());

    std::cout << "Built a trajectory." << std::endl;
    std::cout << "Here's its state (r, th, ph, ek, s): " << traj.stateToString() << std::endl;
    std::cout << "Here are its parameters:\n" << traj.paramString() << std::endl;

    std::cout << "Stepping: " << (traj.step() == Status::Running? "running" : "terminated") << std::endl;
    std::cout << "Here's its state (r, th, ph, ek, s): " << traj.stateToString() << std::endl;
}



