#include "Wavy3D.h"
#include "Wavy3DParams.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <limits>
#include <string>

namespace po = boost::program_options;

int main() {
    std::cout.precision(std::numeric_limits<double>::digits10); // Should be max_digits10, but icpc complains

    /******************************* Point tests ****************************/

    std::cout << "Testing Point class" << std::endl;

    Point p1 = Point();
    Point p2 = Point(1, M_PI / 4, M_PI * 3);
    Point p3 = Point(-3, M_PI / 3, M_PI / 10);
    Point p4 = Point(4, 5 * M_PI / 3, M_PI / 10);
    Point p5 = Point(-5, 5 * M_PI / 3, M_PI / 10);
    Point p6 = Point(2, M_PI / 3, - 3 * M_PI / 10);

    std::cout << "Constructor tests" << std::endl;
    std::cout << "(0, 0, 0)? : " << p1 << std::endl;
    std::cout << "(1, " << M_PI / 4 << ", " << M_PI / 3 << ")? : " << p2 << std::endl;
    std::cout << "(3, " << M_PI * 2/3 << ", " << M_PI * 11/10 << ")? : " << p3 << std::endl;
    std::cout << "(4, " << M_PI * 1/3 << ", " << M_PI * 11 / 10 << ")? : " << p4 << std::endl;
    std::cout << "(5, " << M_PI * 2/3 << ", " << M_PI / 10 << ")? : " << p5 << std::endl;
    std::cout << "(2, " << M_PI * 1/3 << ", " << M_PI * 17 / 10 << ")? : " << p6 << std::endl;
    std::cout << "\tAll tests passed!" << std::endl;

    std::cout << "\nOperator overload tests" << std::endl;
    std::cout << p1 << "? : " << p1 + p2 << std::endl;
    std::cout << p1 << "? : " << p1 - p2 << std::endl;
    std::cout << "(4, " << M_PI * 2/3 << ", " << M_PI * 14 / 10 << ")? : " << 2 * p6 << std::endl;
    std::cout << "(4, " << M_PI * 2/3 << ", " << M_PI * 14 / 10 << ")? : " << p6 * 2 << std::endl;
    std::cout << "true? : " << (p6 == p6) << std::endl;
    std::cout << "false? : " << (3 * p6 == p6 * 2) << std::endl;
    std::cout << "(2, " << M_PI * 1/6 << ", " << M_PI * 11 / 20 << ")? : " << p4 / 2 << std::endl;
    std::cout << p1 << "? : " << p4 + p5 << std::endl;
    std::cout << p1 << "? : " << p4 - p3 << std::endl;

    /******************************* Nelder-Mead tests ****************************/

    std::cout << "\nNM test" << std::endl;

    // Construct a wavy trajectory
    Wavy3D<Wavy3DParams> wavy("config/wavyparams.config");

    // Reset the particle's position (energy doesn't matter)
    wavy.initialize(1, M_PI / 4, M_PI / 4, 1);

    // Get distance to HCS
    std::tuple<Point, double> hcsPoint = wavy.getLHCS();
    std::cout << "Point in HCS: " << std::get<0>(hcsPoint) << " au" << std::endl;
    std::cout << "Distance to HCS: " << std::get<1>(hcsPoint) << " au" << std::endl;

    return 0;
}


