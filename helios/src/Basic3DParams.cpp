#include "Basic3DParams.h"


Basic3DParams::Basic3DParams() : Parameters() {
    desc.add_options()
        // Strauss' model's parameters
        ("deltar", po::value<double>()->default_value(0.01)
                ->notifier([this](const double& deltarPar) { this->deltar = deltarPar; }),
                "used for numerical radial derivatives (au)")
        ("lambda0", po::value<double>()->default_value(0.15)
                ->notifier([this](const double& lambda0Par) { this->lambda0 = lambda0Par; }),
                "parallel mean free path constant (au)")
        ("rRefLambda", po::value<double>()->default_value(1.0)
                ->notifier([this](const double& rRefLambdaPar) { this->rRefLambda = rRefLambdaPar; }),
                "reference distance in mean free path (au)")
        ("kperp_kpar", po::value<double>()->default_value(0.01)
                ->notifier([this](const double& kperp_kparPar) { this->kperp_kpar = kperp_kparPar; }),
                "ratio of k_perp to k_parallel")
        ("rig0", po::value<double>()->default_value(1.0)
                ->notifier([this](const double& rig0Par) { this->rig0 = rig0Par; }),
                "reference rigidity in mean free path (GV)")
        ("ds", po::value<double>()->default_value(350)
                ->notifier([this](const double& dsPar) { this->ds = dsPar; }),
                "timestep (s)")
        ("b0", po::value<double>()->default_value(5.0)
                ->notifier([this](const double& b0Par) { this->b0 = b0Par; }),
                "reference magnetic field strength (nT)")
        ("Ac", po::value<int>()
                ->notifier([this](const double& AcPar) { this->Ac = AcPar; }),
                "HMF polarity (+/- 1)") // TODO: should be read in from a file
        //("alpha", po::value<double>(), "HCS tilt angle (rad)") // TODO: implement and read from file
        // Jupiter parameters
        ("ph0Jup", po::value<double>()->default_value(M_PI)
                ->notifier([this](const double& ph0JupPar) { this->ph0Jup = ph0JupPar; }),
                "Jupiter's initial azimuthal position")
        ("omegaJup", po::value<double>()->default_value(2 * M_PI / (4333 * 3600 * 24))
                ->notifier([this](const double& omegaJupPar) { this->omegaJup = omegaJupPar; }),
                "Jupiter's orbital period")
        ("rBeginJup", po::value<double>()->default_value(5.2 - 0.0477 * 2)
                ->notifier([this](const double& rBeginJupPar) { this->rBeginJup = rBeginJupPar; }),
                "inner boundary of Jupiter's volume")
        ("rEndJup", po::value<double>()->default_value(5.2 + 0.095 * 2)
                ->notifier([this](const double& rEndJupPar) { this->rEndJup = rEndJupPar; }),
                "outer boundary of Jupiter's volume")
        ("dthJup", po::value<double>()->default_value(0.009 * 2)
                ->notifier([this](const double& dthJupPar) { this->dthJup = dthJupPar; }),
                "polar size of Jupiter's volume")
        ("dphJup", po::value<double>()->default_value(0.009 * 2)
                ->notifier([this](const double& dphJupPar) { this->dphJup = dphJupPar; }),
                "azimuthal size of Jupiter's volume")
        ("thJup", po::value<double>()->default_value(M_PI / 2)
                ->notifier([this](const double& thJupPar) { this->thJup = thJupPar; }),
                "Jupiter's polar angle")
    ;
}

std::string Basic3DParams::paramString() const {
    if (valid) {
        // Use high precision
        std::ostringstream out;
        out << std::setprecision(std::numeric_limits<double>::digits10);

        out << Parameters::paramString()
            << "deltar=" << getDeltar() << " # au\n"
            << "lambda0=" << getLambda0() << " # au\n"
            << "rRefLambda=" << getRRefLambda() << " # au\n"
            << "rig0=" << getRig0() << " # GV\n"
            << "b0=" << getB0() << " # nT\n"
            << "kperp_kpar=" << getKperp_kpar() << "\n"
            << "Ac=" << getAc() << "\n"
            << "ds=" << getDs() << " # s\n"
            << "ph0Jup=" << getPh0Jup() << " # rad\n"
            << "omegaJup=" << getOmegaJup() << " # rad / s\n"
            << "rBeginJup=" << getRBeginJup() << " # au\n"
            << "rEndJup=" << getREndJup() << " # au\n"
            << "dthJup=" << getDthJup() << " # rad\n"
            << "dphJup=" << getDphJup() << " # rad\n"
            << "thJup=" << getThJup() << " # rad\n";

        return out.str();
    } else {
        // TODO: throw exception
        return "Parameters have not been loaded.";
    }
}


