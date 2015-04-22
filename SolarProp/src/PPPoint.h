#ifndef	PPPOINT_H
#define PPPOINT_H

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
			if (r < 0) {
				setTh(PPPoint::pi - th);
				setPhi(phi + PPPoint::pi);
				r = -r0;
			} else {
				r = r0;
			}

			return *this;
		};

		PPPoint& setTh(double th0) {
			if (th < 0) {
				th = -th;
				setPhi(phi + PPPoint::pi);
			} else if (th > PPPoint::pi) {
				while(th  > PPPoint::pi) {
					th = 2 * PPPoint::pi - th;
				}
				setPhi(phi + PPPoint::pi);
			} else {
				th = th0;
			}

			return *this;
		};

		PPPoint& setPhi(double phi0) {
			if (phi < 0) {
				while (phi < 0) {
					phi += 2 * PPPoint::pi;
				}
			} else if (phi > 2 * PPPoint::pi) {
				while (phi > 2 * PPPoint::pi) {
					phi -= 2 * PPPoint::pi;
				}
			} else {
				phi = phi0;
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
	
		constexpr static double pi = 3.14159265359;

	private:
		double r, th, phi, e, s;
};

#endif


