#ifndef	PPPOINT_H
#define PPPOINT_H

// Represents a point in a pseudoparticle's trajectory.  s is the backwards time -t.
class PPPoint {
	public:
		PPPoint(double r0, double th0, double phi0, double ei, double s0)
			: e(ei), s(s0) {
			setR(r0);
			setTh(th0);
			setPhi(phi0);
		};
		~PPPoint() { };

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
	
		constexpr static double pi = 3.14159265359;

	private:
		double r, th, phi, e, s;
};

#endif


