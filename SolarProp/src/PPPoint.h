#ifndef	PPPOINT_H
#define PPPOINT_H

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

// TODO: make a .cpp file

// Represents a point in a pseudoparticle's trajectory.  s is the backwards time -t.
class PPPoint {
	public:
		PPPoint(double r0, double th0, double phi0, double ei, double s0)
				: e(ei), s(s0) {
			setR(r0);
			setTh(th0);
			setPhi(phi0);
		};

		// Getters
		double getR() const { return r; };
		double getTh() const { return th; };
		double getPhi() const { return phi; };
		double getE() const { return e; };
		double getS() const { return s; };

		PPPoint& setR(double r0) {
            r = r0;

            // Make sure r is non-negative
			if (r < 0) {
				setTh(M_PI - th);
				setPhi(phi + M_PI);
				r = -r0;
			}

			return *this;
		};

		PPPoint& setTh(double th0) {
            th = th0;

            // Check whether 0 <= th <= pi
			if (th < 0) {
                while (th < 0) {
                    th += 2 * M_PI;
                }
				setTh(th);
			} else if (th > M_PI && th <= 2 * M_PI) {
                th = 2 * M_PI - th;
                setPhi(phi + M_PI);
			} else if (th > 2 * M_PI) {
                setTh(th - 2 * M_PI);
            }

			return *this;
		};

		PPPoint& setPhi(double phi0) {
            phi = phi0;

            // Check whether 0 <= phi <= 2 pi
			if (phi < 0) {
				while (phi < 0) {
					phi += 2 * M_PI;
				}
			} else if (phi > 2 * M_PI) {
				while (phi > 2 * M_PI) {
					phi -= 2 * M_PI;
				}
			}

			return *this;
		};

		PPPoint& setE(double ei) {
			e = ei;
			return *this;
		};

		PPPoint& setS(double s0) {
			s = s0;
			return *this;
		};

        /*
         * Generates a string containing XML representing this PPPoint.
         * Arguments:
         *  indents: number of tabs to put in front of each line.
         *  sMax: maximum value of backwards time.  sMax - s = t.
         * Returns:
         *  string with XML fields for each of this PPPoint's coordinates.
         */
        std::string toXML(int indents, const double& sMax = 0.0) const {
            std::string tabs;

            for (int i = 0; i < indents; i++) {
                tabs += "\t";
            }

            // Full precision version of to_string.  Clear the stringstream, put each variable in, extract
            // and concatenate with xml.
            std::stringstream converter;
            converter << std::setprecision(std::numeric_limits<double>::digits10);
            // TODO: write a helper function...
            converter << r;
            std::string xml = tabs + "<point>\n" + tabs + "\t<r>" + converter.str() + "</r>\n";

            converter.str("");
            converter << th;
            xml += tabs + "\t<th>" + converter.str() + "</th>\n";

            converter.str("");
            converter << phi;
            xml += tabs + "\t<phi>" + converter.str() + "</phi>\n";

            converter.str("");
            converter << e;
            xml += tabs + "\t<e>" + converter.str() + "</e>\n";

            converter.str("");
            converter << sMax - s;
            xml += tabs + "\t<t>" + converter.str() + "</t>\n" + tabs + "</point>\n";

            return xml;
        }

	private:
		double r, th, phi, e, s;
};

#endif


