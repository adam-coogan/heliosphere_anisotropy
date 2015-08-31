#include <cmath>

class KTensor {
    public:
        double rr, phph, rph, thth;
        double dKrr_dr, dKphph_dr, dKrph_dr, dKthth_dr;
        double dKrr_dth, dKphph_dth, dKrph_dth, dKthth_dth;

        KTensor() : rr(0), phph(0), rph(0), thth(0) { };

        KTensor(double r, double th, double v, double P) {
            updateK(r, th, v, P);
        }

        KTensor& updateK(double rc, double thc, double vc, double Pc) {
            // Compute elements at shifted points and store elements
            updateKElements(rc + deltar, thc, vc, Pc);

            // TODO: CHANGE DERIVATIVES BACK!!!
            double krr_r = rr;
            //double kphph_r = phph;
            double krph_r = rph;
            //double kthth_r = thth;

            // Make sure theta is in [0, pi].  No need to worry about phi here.
            if (thc + deltath > M_PI) {
                updateKElements(rc, thc - deltath, vc, Pc);
            } else {
                updateKElements(rc, thc + deltath, vc, Pc);
            }

            //double krr_th = rr;
            //double kphph_th = phph;
            //double krph_th = rph;
            //double kthth_th = thth; 

            // Compute elements at actual point
            updateKElements(rc, thc, vc, Pc);

            // Compute r derivatives
            dKrr_dr = (krr_r - rr) / deltar;
            ///dKphph_dr = (kphph_r - phph) / deltar;
            dKrph_dr = (krph_r - rph) / deltar;
            //dKthth_dr = (kthth_r - thth) / deltar;

            // Compute th derivatives
            //dKrr_dth = (krr_th - rr) / deltath;
            //dKphph_dth = (kphph_th - phph) / deltath;
            //dKrph_dth = (krph_th - rph) / deltath;
            //dKthth_dth = (kthth_th - thth) / deltath;
            dKthth_dth = 0; // TODO: CHANGE BACK!!!

            return *this;
        };

    private:
        KTensor& updateKElements(double rc, double thc, double vc, double Pc) {
            double kpar = (vc / 3) * lambda0 * (1 + rc / r0) * (Pc >= P0? Pc / P0 : 1);

            double kperp = kperp_kpar * kpar;

            // Convert to spherical coordinates
            double tanPsi = Omega * (rc - rSun) * sin(thc) / Vsw;
            double cosPsi = 1 / sqrt(1 + tanPsi*tanPsi);
            double sinPsi = sqrt(1 - cosPsi*cosPsi);

            rr = kpar * cosPsi*cosPsi + kperp * sinPsi*sinPsi;
            phph = kpar * sinPsi*sinPsi + kperp * cosPsi*cosPsi;
            rph = (kperp - kpar) * cosPsi * sinPsi;
            thth = kperp;

            return *this;
        };
}
