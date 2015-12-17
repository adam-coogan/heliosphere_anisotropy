#include "Basic3DParams.h"


Basic3DParams::Basic3DParams() : Parameters() {
    desc.add_options()
        // Strauss' model's parameters
        ("deltar", po::value<double>()->default_value(0.01), "used for numerical radial derivatives (au)")
        ("lambda0", po::value<double>()->default_value(0.15), "parallel mean free path constant (au)")
        ("rRefLambda", po::value<double>()->default_value(1.0), "reference distance in mean free path (au)")
        ("kperp_kpar", po::value<double>()->default_value(0.01), "ratio of k_perp to k_parallel")
        ("rig0", po::value<double>()->default_value(1.0), "reference rigidity in mean free path (GV)")
        ("ds", po::value<double>()->default_value(350), "timestep (s)")
        ("b0", po::value<double>()->default_value(5.0), "reference magnetic field strength (nT)")
        ("Ac", po::value<int>(), "HMF polarity (+/- 1)") // TODO: should be read in from a file
        //("alpha", po::value<double>(), "HCS tilt angle (rad)") // TODO: implement and read from file
        // Jupiter parameters
        ("ph0Jup", po::value<double>()->default_value(M_PI), "Jupiter's initial azimuthal position")
        ("omegaJup", po::value<double>()->default_value(2 * M_PI / (4333 * 3600 * 24)), "Jupter's orbital "
                                                                                                    "period")
        ("rBeginJup", po::value<double>()->default_value(5.2 - 0.0477 * 2), "inner boundary of Jupiter's "
                                                                                                    "volume")
        ("rEndJup", po::value<double>()->default_value(5.2 + 0.095 * 2), "outer boundary of Jupiter's volume")
        ("dthJup", po::value<double>()->default_value(0.009 * 2), "polar size of Jupiter's volume")
        ("dphJup", po::value<double>()->default_value(0.009 * 2), "azimuthal size of Jupiter's volume")
        ("thJup", po::value<double>()->default_value(M_PI / 2), "Jupiter's polar angle")
    ;
}

std::string Basic3DParams::paramString() const {
    if (valid) {
        // Use high precision
        std::ostringstream out;
        out << std::setprecision(std::numeric_limits<double>::digits10);

        out << Parameters::paramString()
            << "deltar=" << deltar() << " # au\n"
            << "lambda0=" << lambda0() << " # au\n"
            << "rRefLambda=" << rRefLambda() << " # au\n"
            << "rig0=" << rig0() << " # GV\n"
            << "b0=" << b0() << " # nT\n"
            << "kperp_kpar=" << kperp_kpar() << "\n"
            << "Ac=" << Ac() << "\n"
            << "ds=" << ds() << " # s\n"
            << "ph0Jup=" << ph0Jup() << " # rad\n"
            << "omegaJup=" << omegaJup() << " # rad / s\n"
            << "rBeginJup=" << rBeginJup() << " # au\n"
            << "rEndJup=" << rEndJup() << " # au\n"
            << "dthJup=" << dthJup() << " # rad\n"
            << "dphJup=" << dphJup() << " # rad\n"
            << "thJup=" << thJup() << " # rad\n";

        return out.str();
    } else {
        // TODO: throw exception
        return "Parameters have not been loaded.";
    }
}


