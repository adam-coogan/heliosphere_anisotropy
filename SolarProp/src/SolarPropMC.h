#ifndef SOLARPROPMC_H
#define SOLARPROPMC_H

#include "PPTrajectory.h"
#include "RunConfig.h"
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace SolarPropMC {
	void runMCs(const RunConfig& rcs, const std::string& paramFileName);
}

#endif


