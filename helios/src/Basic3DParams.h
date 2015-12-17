#ifndef BASIC3DPARAMS
#define BASIC3DPARAMS 

#include "Parameters.h"
#include <iomanip>
#include <iostream>
#include <limits>

class Basic3DParams : public Parameters {
    public:
        Basic3DParams();

        double deltar() const { return params["deltar"].as<double>(); };
        double lambda0() const { return params["lambda0"].as<double>(); };
        double rRefLambda() const { return params["rRefLambda"].as<double>(); };
        double rig0() const { return params["rig0"].as<double>(); };
        double ds() const { return params["ds"].as<double>(); };
        double kperp_kpar() const { return params["kperp_kpar"].as<double>(); };
        int Ac() const { return params["Ac"].as<int>(); };
        double b0() const { return params["b0"].as<double>(); };
        double ph0Jup() const { return params["ph0Jup"].as<double>(); };
        double omegaJup() const { return params["omegaJup"].as<double>(); };
        double rBeginJup() const { return params["rBeginJup"].as<double>(); };
        double rEndJup() const { return params["rEndJup"].as<double>(); };
        double dthJup() const { return params["dthJup"].as<double>(); };
        double dphJup() const { return params["dphJup"].as<double>(); };
        double thJup() const { return params["thJup"].as<double>(); };

        /*!
         * Adds Jupiter parameters to the parameter string.
         * \sa Parameters::paramString()
         */
        virtual std::string paramString() const;
};

#endif


