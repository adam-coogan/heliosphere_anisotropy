#include "Wavy3DParams.h"

Wavy3DParams::Wavy3DParams() : Basic3DParams() {
    desc.add_options()
        ("alpha", po::value<double>()->default_value(0)
                ->notifier([this](const double& alphaPar) { this->alpha = alphaPar; }),
                "heliospheric current sheet tilt angle")
        ("phPhase", po::value<double>()->default_value(0)
                ->notifier([this](const double& phPhasePar) { this->phPhase = phPhasePar; }),
                "azimuthal phase in function for computing heliospheric current sheet's angular extent")
    ;
}

std::string Wavy3DParams::paramString() const {
    if (valid) {
        // Use high precision
        std::ostringstream out;
        out << std::setprecision(std::numeric_limits<double>::digits10);

        out << Parameters::paramString()
            << "alpha=" << getAlpha() << " # rad\n"
            << "phPhase=" << getPhPhase() << " # rad\n";
        
        return out.str();
    } else {
        // TODO: throw exception
        return "Parameters have not been loaded.";
    }
}


