#ifndef POINT
#define POINT 

#include <cmath>
#include <stdexcept>

class Point {
    public:
        Point(double rn, double thn, double phn) {
            set(rn, thn, phn);
        };

        /*!
         * \arg rn new radial coordinate for the point
         * \return this point with r set to rn.  If rn < 0, th is set to pi + th.
         */
        Point& setR(double rn) {
            if (rn < 0) {
                r = std::abs(rn);
                setTh(M_PI + th);
            } else {
                r = rn;
            }

            return *this;
        };

        /*!
         * \arg thn new polar angle for the point
         * \return this point with th set to thn
         */
        Point& setTh(double thn) {
            th = thn;
            return renormalizeTh();
        };

        /*!
         * \arg phn new azimuthal angle for the point
         * \return this point with ph set to phn
         */
        Point& setPh(double phn) {
            ph = phn;
            return renormalizePh();
        };

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
        Point& renormalizeTh() {
            while (th > M_PI) {
                th = 2 * M_PI - th;
                ph = ph - M_PI;
            }

            while (th < 0) {
                th = -th;
                ph += M_PI;
            }

            return renormalizePh();
        };

        /*!
         * \return this point with ph renormalized to lie within [0, 2 pi]
         */
        Point& renormalizePh() {
            while (ph > 2 * M_PI) {
                ph -= 2 * M_PI;
            }

            while (ph < 0) {
                ph += 2 * M_PI;
            }

            return *this;
        };
};

#endif


