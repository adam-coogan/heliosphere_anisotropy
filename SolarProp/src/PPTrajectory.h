/*
 * This code implements the stochastic differential equation model for transport of electrons/positrons in the
 * heliosphere that is described in "Modeling the Modulation of Galactic and Jovian Electrons by Stochastic
 * Processes," by Strauss, Potgieter, Busching and Kopp, APJ 735:32 (13pp), 2011 July 10.
 */

// TODO: invalidate if parameter file could not be found!!!									 2015-02-24 11:40
// TODO: should be able to give angle at which the particle was detected at Earth.			 2015-02-25 01:11
// TODO: make an integrator class that determines what method the trajectory uses.  Should be able to give it
//		 the pseudoparticle integrator, the deterministic one I've implemented in Mathematica, or one that
//		 calculates single particle trajectories including the random components, which is what I need to
//		 find a way to implement.															 2015-03-06 16:45

// Prints name and value of a variable p to the given stream
#define DEBUG_PRINT_PARAM(stream, p) \
	// (stream) << #p << " = " << (p) << std::endl;
#define DEBUG_PRINT(stream, str) \
	// (stream) << (str)
#define DEBUG_FLUSH(stream) \
	// (stream) << std::endl

#ifndef PPTRAJECTORY_H
#define PPTRAJECTORY_H

#include "PPPoint.h"
#include "HelioParams.h"
#include "SphericalVector.h"
#include "SphericalTensor.h"

#include <chrono> // For initializing the random number generator
#include <cmath>
#include <fstream>
#include <limits>
#include <random> // For the Wiener terms
#include <sstream>
#include <string> // For writing to XML
#include <vector>

// Enumerates possible statuses of PPTrajectory's integrator
enum class BoundaryHit {Sun, Heliopause, None};

// Represents a pseudoparticle trajectory.  The pseudoparticle has the ability to move itself forward in time.
class PPTrajectory {
	public:
		/*
		 * Constructor.
		 * Arguments:
		 *	ri, thi, phii: "initial" (in backwards time) coordinates relative to sun and ecliptic plane.
		 *	ei: energy of particle at time of detection.
		 *	electron: true if the particle is an electron, false if it is a positron.
		 *	gti: true if A_c (the sign of the solar cycle) is +1, false if it is -1.
		 *	paramFileName (optional): file containing simulation parameters.
		 */	
		PPTrajectory(double ri, double thi, double phii, double ei, const std::string& paramFileName);
		~PPTrajectory() { };

		/*
		 * Steps the integration forward.  
		 * Arguments:
		 *	n: number of steps to take backward.  A negative or 0 value means the integration will continue
		 *	   until the particle reaches the heliopause.  This is the default.
		 * Returns:
		 *	A BoundaryHit specifying whether the particle has hit the heliopause, entered the sun or whether
		 *	it has not reached a boundary yet.
		 */
		BoundaryHit integrate(int n = -1);

		/* Queries the integration status.
		 * Returns:
		 *	A BoundaryHit specifying whether the particle has hit the heliopause, entered the sun or whether
		 *	it has not reached a boundary yet.
		 */
		BoundaryHit getStatus() const { return status; };

		/*
		 * Returns current trajectory to a file as XML.  XML format:
		 * <trajectory datetime=...>
		 *		<params>
		 *			<ds> ds </ds>
		 *			<lambda0> lambda0 </lambda0>
		 *			...
		 *		</params>
		 *		<points>
		 *			<point>
		 *				<r> r </r>
		 *				<th> th </th>
		 *				...
		 *				<t> sMax-s </t>
		 *			</point>
		 *			...
		 *		</points>
		 * </trajectory>
		 * Returns:
		 *	*this.
		 */
		std::string toXML() const;

		/*
		 * Convenience function that writes the XML returned from toXML to a file.
		 * Arguments:
		 *	fname: the name of the file in which to write the XML.
		 * Returns:
		 *	*this.
		 */
		PPTrajectory& writeToXML(const std::string& fname);

		/*
		 * Getters for the "initial conditions" (in backwards time) of the trajectory.
		 */
		double getR0() const { return traj.front().getR(); };
		double getTh0() const { return traj.front().getTh(); };
		double getPhi0() const { return traj.front().getPhi(); };
		double getE0() const { return traj.front().getE(); };

	private:
		// Vector of points in the pseudoparticle trajectory.
		std::vector<PPPoint> traj;
		//Indicates whether the boundary has been reached
		BoundaryHit status;

		// Magnetic field
		SphericalVector b;
		// curl(B/|B|)
		SphericalVector curlBoverB;
		// Diffusion tensor, k.  Depends on b.
		SphericalTensor kTensor;
		double kpar, kperp;
		// Drift velocity
		SphericalVector vdrift;
		// True if relevant quantiy has been updated, false after moving to a new point until step is called
		bool updatedB, updatedVdrift, updatedKTensor; 

		// Used for generating random variables from N(0, 1) for the Wiener terms in the SDE
		std::default_random_engine generator;
		std::normal_distribution<double> ndistro;

		// Parameter container class
		const HelioParams params;

		// Unchanging constants: speed of light and pi
		constexpr static double cAUs = 0.0020039888; // c in AU/s
		constexpr static double pi = 3.14159265359;

		/* 
		 * Steps the integration forward once.  This can be called even after the particle has hit the sun or
		 * heliopause.  Called by integrate.  Updates status if the particle reaches a boundary.
		 * Returns:
		 *	*this.
		 */
		PPTrajectory& step();

		/*
		 * Computes velocity.
		 * Returns:
		 *	Particle velocity in AU/s.
		 */
		double vel() const;

		/*
		 * Computes rigidity.
		 * Returns:
		 *	Particle rigidity in T*AU.
		 */
		double rigidity() const;

		/*
		 * Gives lambda_parallel (eq 4).
		 * Returns:
		 *	The component of the diffusion tensor k parallel to the magnetic field in HMF (heliospheric
		 *	magnetic field) aligned coordinates in AU.
		 */
		double lambdaPar() const;

		/*
		 * Drift reduction factor at current point (eq 10).
		 * Returns:
		 *	Drift velocity reduction factor.  No units.
		 */
		double fs() const;

		/*
		 * Gamma parameter at current point (eq 3).  This factor appears in the adiabatic cooling term in the
		 * transport equation.
		 * Returns:
		 *	Gamma.  No units.
		 */
		double gamma() const;

		// TODO: implement wavy HCS.  For now keep it flat.									 2015-02-23 14:15
		/*
		 * Gets theta'.  This is the angular extent of the current sheet at a given point.
		 * Returns:
		 *	Theta'.
		 */
		double getThp() const { return pi / 2; };

		/*
		 * Heaviside step function.  Used to separate the HMF into regions of opposite polarity.
		 * Returns:
		 *	1 if th <= thp
		 *	//0 if th = thp
		 *	-1 if th > thp.
		 */
		double heaviside(double th, double thp) const {
			if (th <= thp) {
				return 1;
			//} else if (th == thp) {
			//	return 0;
			} else {
				return -1;
			}
		};

		/*
		 * Gives magnetic field at current point (eq 5).  Recomputes b (in T) and curlBoverB at the current
		 * point.
		 * Returns:
		 *	*this.
		 */
		PPTrajectory& updateB(); // I think this is where I need to use move...

		/*
		 * Generate the diffusion tensor.  Calls updateB if it has not yet been called at the current point.
		 * Components of k have units of AU^2/s.
		 * Returns:
		 *	*this.
		 */
		PPTrajectory& updateKTensor();

		/*
		 * Computes pitch-angle averaged guiding center drift velocity in spherical coordinates.  Updates
		 * vdrift, which has units of AU/s.  This takes into account the gradient and curvature components as
		 * well as those due to drift in the HCS.
		 * Returns:
		 *	*this.
		 */
		PPTrajectory& updateVdrift();

		/*
		 * Creates a string containing the values of the parameters used in the current model.
		 * Returns:
		 *	A reference to a string containing this run's parameter values, formatted as attributes of an XML
		 *	tag.
		 */
		std::string paramsToString() const;
};

#endif


