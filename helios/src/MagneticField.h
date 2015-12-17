#ifndef MAGNETICFIELD
#define MAGNETICFIELD 

#include <cmath>

/*!
 * Magnetic field container.
 */
struct MagneticField {
    double r, th, ph;

    /*!
     * \return current magnitude of magnetic field.
     */
    double magnitude() {
        return sqrt(r*r + th*th + ph*ph);
    }
};
#endif


