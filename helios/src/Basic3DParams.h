#ifndef BASIC3DPARAMS
#define BASIC3DPARAMS 

#include "Parameters.h"
#include <iomanip>
#include <iostream>
#include <limits>

class Basic3DParams : public Parameters {
    public:
        Basic3DParams();

        double getDeltar() const { return deltar; };
        double getLambda0() const { return lambda0; };
        double getRRefLambda() const { return rRefLambda; };
        double getRig0() const { return rig0; };
        double getDs() const { return ds; };
        double getKperp_kpar() const { return kperp_kpar; };
        int getAc() const { return Ac; };
        double getB0() const { return b0; };
        double getPh0Jup() const { return ph0Jup; };
        double getOmegaJup() const { return omegaJup; };
        double getRBeginJup() const { return rBeginJup; };
        double getREndJup() const { return rEndJup; };
        double getDthJup() const { return dthJup; };
        double getDphJup() const { return dphJup; };
        double getThJup() const { return thJup; };

        /*!
         * Adds Jupiter parameters to the parameter string.
         * \sa Parameters::paramString()
         */
        virtual std::string paramString() const;

    protected:
        double deltar;
        double lambda0;
        double rRefLambda;
        double rig0;
        double ds;
        double kperp_kpar;
        int Ac;
        double b0;
        double ph0Jup;
        double omegaJup;
        double rBeginJup;
        double rEndJup;
        double dthJup;
        double dphJup;
        double thJup;
};

#endif


