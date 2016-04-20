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
    std::cout << "\t(0, 0, 0)? : " << p1 << std::endl;
    std::cout << "\t(1, " << M_PI / 4 << ", " << M_PI << ")? : " << p2 << std::endl;
    std::cout << "\t(3, " << M_PI * 2/3 << ", " << M_PI * 11/10 << ")? : " << p3 << std::endl;
    std::cout << "\t(4, " << M_PI * 1/3 << ", " << M_PI * 11 / 10 << ")? : " << p4 << std::endl;
    std::cout << "\t(5, " << M_PI * 2/3 << ", " << M_PI / 10 << ")? : " << p5 << std::endl;
    std::cout << "\t(2, " << M_PI * 1/3 << ", " << M_PI * 17 / 10 << ")? : " << p6 << std::endl;
    std::cout << "All tests passed!" << std::endl;

    std::cout << "\nOperator overload tests" << std::endl;
    std::cout << "\t" << p2 << "? : " << p1 + p2 << std::endl;
    std::cout << "\t(1, " << 3.0/4 * M_PI << ", " << 2*M_PI << ")? : " << p1 - p2 << std::endl;
    std::cout << "\t(4, " << M_PI * 1/3 << ", " << M_PI * 17 / 10 << ")? : " << 2 * p6 << std::endl;
    std::cout << "\t(4, " << M_PI * 1/3 << ", " << M_PI * 17 / 10 << ")? : " << p6 * 2 << std::endl;
    std::cout << "\t(2, " << M_PI * 1/3 << ", " << M_PI * 11 / 10 << ")? : " << p4 / 2 << std::endl;

    std::cout << "Addition" << std::endl;
    double newX4p5 = 4 * sin(M_PI / 3) * cos(11 * M_PI / 10) + 5 * sin(M_PI * 2/3) * cos(M_PI / 10);
    double newY4p5 = 4 * sin(M_PI / 3) * sin(11 * M_PI / 10) + 5 * sin(M_PI * 2/3) * sin(M_PI / 10);
    double newZ4p5 = 4 * cos(M_PI / 3) + 5 * cos(M_PI * 2/3);
    double newR4p5 = sqrt(newX4p5*newX4p5 + newY4p5*newY4p5 + newZ4p5*newZ4p5);
    std::cout << "\t(" << newR4p5 << ", " << acos(newZ4p5 / newR4p5) << ", " << atan(newY4p5 / newX4p5)
                << ")? : " << p4 + p5 << std::endl;

    std::cout << "Subtraction" << std::endl;
    double newX4m3 = 4 * sin(M_PI / 3) * cos(11 * M_PI / 10) - 3 * sin(M_PI * 2/3) * cos(M_PI * 11 / 10);
    double newY4m3 = 4 * sin(M_PI / 3) * sin(11 * M_PI / 10) - 3 * sin(M_PI * 2/3) * sin(M_PI * 11 / 10);
    double newZ4m3 = 4 * cos(M_PI / 3) - 3 * cos(M_PI * 2/3);
    double newR4m3 = sqrt(newX4m3*newX4m3 + newY4m3*newY4m3 + newZ4m3*newZ4m3);
    std::cout << "\t(" << newR4m3 << ", " << acos(newZ4m3 / newR4m3) << ", "
        << 2*M_PI + atan2(newY4m3, newX4m3) << ")? : " << p4 - p3 << std::endl;

    std::cout << "Equality" << std::endl;
    std::cout << "\ttrue (1)? : " << (p6 == p6) << std::endl;
    std::cout << "\tfalse (0)? : " << (3 * p6 == p6 * 2) << std::endl;

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


