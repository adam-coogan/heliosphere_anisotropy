#ifndef POINT
#define POINT 

#include <cmath>
#include <stdexcept>

class Point {
    public:
        //! Default constructor
        Point() { set(0, 0, 0); };

        Point(double rn, double thn, double phn) {
            set(rn, thn, phn);
        };

        double getR() const { return r; };
        double getTh() const { return th; };
        double getPh() const { return ph; };

        /*!
         * \arg rn new radial coordinate for the point
         * \return this point with r set to rn.  If rn < 0, th is set to pi + th.
         */
        void setR(double rn) {
            if (rn < 0) {
                r = std::abs(rn);
                setTh(M_PI + th);
            } else {
                r = rn;
            }
        };

        //! Increments r
        void incrR(double deltaR) { setR(r + deltaR); };

        /*!
         * \arg thn new polar angle for the point
         * \return this point with th set to thn
         */
        void setTh(double thn) {
            th = thn;
        };

        //! Increments th
        void incrTh(double deltaTh) { setTh(r + deltaTh); };

        /*!
         * \arg phn new azimuthal angle for the point
         * \return this point with ph set to phn
         */
        void setPh(double phn) {
            ph = phn;
        };

        //! Increments ph
        void incrPh(double deltaPh) { setPh(r + deltaPh); };

        //! Resets all coordinates, renormalizing if necessary
        void set(double rn, double thn, double phn) {
            setR(rn);
            ph = phn; // Ok since setTh() calls renormalizePh
            setTh(thn);
        };

    private:
        double r;
        double th;
        double ph;

        /*!
         * \return this point with th and ph renormalized to lie within [0, pi] and [0, 2 pi]
         */
        void renormalizeTh() {
            while (th > M_PI) {
                th = 2 * M_PI - th;
                ph = ph - M_PI;
            }

            while (th < 0) {
                th = -th;
                ph += M_PI;
            }
        };

        /*!
         * \return this point with ph renormalized to lie within [0, 2 pi]
         */
        void renormalizePh() {
            while (ph > 2 * M_PI) {
                ph -= 2 * M_PI;
            }

            while (ph < 0) {
                ph += 2 * M_PI;
            }
        };
};

#endif


