#include "PPTrajectory.h"

#include <iostream>
#include <sstream>
#include <iomanip>

int main(int argc, char* argv[]) {
	if (argc == 2) {
		PPTrajectory traj(1, 3.14159265359 / 2 * 1.01, 0, std::stod(argv[1]));

		std::cout << "Generating trajectory for particle with measured energy " << traj.getE0() << "..."
			<< std::endl;

		// Run the full simulation
		traj.integrate(0);

		std::cout << "Simulation complete!" << std::endl;

		// Write results to XML
		traj.writeToXML("runs/debug.xml");
	} else {
		// Should this print to std::cerr?
		std::cout << "SolarProp error: must supply an argument specifying the pseudoparticle's initial energy"
			" in GeV." << std::endl;

		return 1;
	}
}


