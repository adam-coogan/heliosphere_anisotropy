#include "Parameters.h"

Parameters::Parameters() : desc("Simulation parameters") {
    // Invalid until parameters are loaded
    valid = false;

    // Add parameter descriptions
    // TODO: some of these should be moved to Basic3DParams.  Only keep the ones needed for the functions in
    // TrajectoryBase!!!
    desc.add_options()
        // Type of particle
        ("m", po::value<double>(), "particle's mass (GeV)")
        ("charge", po::value<int>(), "particle's charge (+/- 1)")
        // Initial position
        ("r0", po::value<double>()->default_value(1.0), "initial distance from center of sun (au)")
        ("th0", po::value<double>()->default_value(M_PI / 2), "initial polar angle (rad)")
        ("ph0", po::value<double>()->default_value(0.0), "initial azimuthal angle (rad)")
        ("ek0", po::value<double>()->default_value(0.0), "particle's initial kinetic energy (GeV)")
        // Minimum set of solar parameters
        ("rHP", po::value<double>()->default_value(140), "distance from center of sun to heliopause (au)")
        ("rSun", po::value<double>()->default_value(0.005), "sun's radius (au)") // TODO: rInner = rSun!
        ("omega", po::value<double>()->default_value(2 * M_PI / (25.4 * 24 * 3600)),
                                                     "sun's angular velocity (rad / s)")
        ("Vsw", po::value<double>()->default_value(400), "solar wind velocity (km / s)")
    ;
}

void Parameters::loadParameters(const std::string& paramFileNameArg) {
    try {
        paramFileName = paramFileNameArg;
        // Open parameter file
        std::ifstream paramFile;
        paramFile.open(paramFileName);

        // Read file into variable map
        po::store(po::parse_config_file(paramFile, desc), params);
        po::notify(params);
        // TODO: validate parameters

        valid = true;
    } catch (std::ifstream::failure e) {
        std::cout << "Parameter file " << paramFileName << " could not be opened." << std::endl;
        valid = false;
    }
}

std::string Parameters::paramString() const {
    if (valid) {
        // Use high precision
        std::ostringstream out;
        out << std::setprecision(std::numeric_limits<double>::digits10);

        out << "r0=" << r0() << " # au\n"
            << "th0=" << th0() << " # rad\n"
            << "ph0=" << ph0() << " # rad\n"
            << "ek0=" << ek0() << " # GeV\n"
            << "rHP=" << rHP() << " # au\n"
            << "rSun=" << rSun() << " # au\n"
            << "m=" << m() << " # GeV\n"
            << "charge=" << charge() << "\n"
            << "omega=" << omega() << " # rad / s\n"
            << "Vsw=" << Vsw() << " # km / s\n";

        return out.str();
    } else {
        // TODO: throw exception
        return "Parameters have not been loaded.";
    }
}


