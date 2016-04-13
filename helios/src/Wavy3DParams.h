#ifndef WAVY3DPARAMS
#define WAVY3DPARAMS 

/*!
 * The most basic implementation of a wavy HCS model.  The tilt angle is hard coded, as is the phase in the
 * function to calculate the angular extent of the HCS.  Future versions will provide data to calculate these
 * as a function of time.
 */

#include "Basic3DParams.h"
#include <string>

class Wavy3DParams : public Basic3DParams {
    public:
        Wavy3DParams();

        double getAlpha() const { return alpha; };
        double getPhPhase() const { return phPhase; };

        /*!
         * Adds wavy sheet parameters to the parameter string.
         */
        virtual std::string paramString() const;

    protected:
        double alpha;
        double phPhase;
};

#endif


