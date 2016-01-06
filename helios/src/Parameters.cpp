#include "Parameters.h"

Parameters::Parameters() : desc("Simulation parameters") {
    // Invalid until parameters are loaded
    valid = false;

    // Add parameter descriptions
    // TODO: some of these should be moved to Basic3DParams.  Only keep the ones needed for the functions in
    // TrajectoryBase!!!
    desc.add_options()
        // Type of particle
        ("m", po::value<double>()
                ->notifier([this](const double& mPar) { this->m = mPar; }),
                "particle's mass (GeV)")
        ("charge", po::value<int>()
                ->notifier([this](const double& chargePar) { this->charge = chargePar; }),
                "particle's charge (+/- 1)")
        // Initial position
        ("r0", po::value<double>()->default_value(1.0)
                ->notifier([this](const double& r0Par) { this->r0 = r0Par; }),
                "initial distance from center of sun (au)")
        ("th0", po::value<double>()->default_value(M_PI / 2)
                ->notifier([this](const double& th0Par) { this->th0 = th0Par; }),
                "initial polar angle (rad)")
        ("ph0", po::value<double>()->default_value(0.0)
                ->notifier([this](const double& ph0Par) { this->ph0 = ph0Par; }),
                "initial azimuthal angle (rad)")
        ("ek0", po::value<double>()->default_value(0.0)
                ->notifier([this](const double& ek0Par) { this->ek0 = ek0Par; }),
                "particle's initial kinetic energy (GeV)")
        // Minimum set of solar parameters
        ("rHP", po::value<double>()->default_value(140)
                ->notifier([this](const double& rHPPar) { this->rHP = rHPPar; }),
                "distance from center of sun to heliopause (au)")
        ("rSun", po::value<double>()->default_value(0.005)
                ->notifier([this](const double& rSunPar) { this->rSun = rSunPar; }),
                "sun's radius (au)") // TODO: rInner = rSun!
        ("omega", po::value<double>()->default_value(2 * M_PI / (25.4 * 24 * 3600))
                ->notifier([this](const double& omegaPar) { this->omega = omegaPar; }),
                "sun's angular velocity (rad / s)")
        ("Vsw", po::value<double>()->default_value(400)
                ->notifier([this](const double& VswPar) { this->Vsw = VswPar; }),
                "solar wind velocity (km / s)")
    ;
}

void Parameters::loadParameters(const std::string& paramFileNameArg) {
    try {
        paramFileName = paramFileNameArg;
        // Open parameter file
        std::ifstream paramFile;
        paramFile.open(paramFileName);

        // Temporarily read file into variable map
        po::variables_map params;
        po::store(po::parse_config_file(paramFile, desc), params);
        // Transfer parameter values into local variables
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

        out << "r0=" << getR0() << " # au\n"
            << "th0=" << getTh0() << " # rad\n"
            << "ph0=" << getPh0() << " # rad\n"
            << "ek0=" << getEk0() << " # GeV\n"
            << "rHP=" << getRHP() << " # au\n"
            << "rSun=" << getRSun() << " # au\n"
            << "m=" << getM() << " # GeV\n"
            << "charge=" << getCharge() << "\n"
            << "omega=" << getOmega() << " # rad / s\n"
            << "Vsw=" << getVsw() << " # km / s\n";

        return out.str();
    } else {
        // TODO: throw exception
        return "Parameters have not been loaded.";
    }
}


