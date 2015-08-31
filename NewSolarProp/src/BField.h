#include <cmath>

class BField {
    public:
        // Components
        double br, bph;

        BField() : br(0), bph(0) { };

        BField(double r, double th) {
            updateB(r, th);
        };

        double magnitude() const { return sqrt(br*br + bph*bph); };

        BField& updateB(double rc, double thc) {
            double bFact = Ac * B0 * r0*r0 / (rc*rc);// * (1 - 2 * heaviside(thc - M_PI / 2));

            br = bFact;
            bph = -bFact * (rc - rSun) * Omega / Vsw * sin(thc);

            return *this;
        };
};
