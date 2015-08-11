#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>

////////////////////////////////////////

struct KTensor {
    double rr, phph, rph, thth;

    KTensor() : rr(0), phph(0), rph(0), thth(0) { };

    KTensor(double rr_, double phph_, double rph_, double thth_)
        : rr(rr_), phph(phph_), rph(rph_), thth(thth_) { };
};

struct BField {
    double r, th, ph;

    BField() : r(0), th(0), ph(0) { };

    BField(double r_, double ph_) : r(r_), ph(ph_) { };
};

// Used to indicate the simulation's status after each timestep
enum class Status {Sun, Heliopause, Running};

// Constants
double ds = 350.0; // Timestep, s
double lambda0 = 0.15; // AU
double P0 = 1.0; // GV
double r0 = 1.0; // AU
double B0 = 0.00000000354 * 4.485 * std::pow(10, 10); // GV / AU
double Be = 0.00000000517 * 4.485 * std::pow(10, 10); // GV / AU
double vsw = 400 * 6.685 * std::pow(10, -9); // AU / s
int qSign = -1; // Charge sign: +/- 1
int A = 1; // Solar cycle sign
double m = 5.11 * std::pow(10, -4); // Electron mass, GeV
double omega = 2 * M_PI / (25.4 * 24 * 3600); // rad / s
double rHP = 140; // AU
double rSun = 0.005; // AU
double kperp_kpar = 0.01; // kperp / kpar

// Coordinates
double s;
double r;
double th;
double ph;
// Kinetic energy
double ek;
// Rigidity
double P;
// Gamma constant
double gammaD;
// Velocity
double v;

// Magnetic field
//double Br;
//double Bth;
//double Bph;
BField B;

// K tensor
//double Krr, Krph, Kthth, Kphph;
KTensor K;
// Partial derivatives of K
double dKrr_dr, dKrph_dr;

// Drift velocity
double vdr, vdth, vdph;

// Gaussian random number generator with mean 0 and standard deviation 1
std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 1.0);

////////////////////////////////////////

// Heaviside step function
int heaviside(double th) {
    if (th < M_PI / 2.0) {
        return 1;
    } else if (th > M_PI / 2.0) {
        return -1;
    } else {
        return 0;
    }
}

// Computes rigidity for particle mass m (defined above)
void updateRig() {
    P = std::sqrt(ek * (ek + 2 * m));
}

// Computes the gamma parameter used in the drift velocity
double updateGammaD(double r_, double th_) {
    return omega * r_ * std::sin(th_) / vsw;
}

// Gets particle velocity
void updateV() {
    // c = 0.002004 AU / s
    v = 0.002004 * std::sqrt(ek * (ek + 2 * m)) / (ek + m);
}

// Updates the magnetic field
BField updateB(double r_, double th_) {
    double bFact = B0 * pow(r0 / r_, 2);

    return BField(bFact, -bFact * omega * (r_ - rSun) * std::sin(th_) / vsw);
}

// Updates the diffusion tensor K and its relevant partial derivatives (dKrr/dr, dKrph/dr)
KTensor updateK(double r_, double th_) {
    double kpar = v / 3 * (1 + r_ / r0) * lambda0;
    if (P >= P0) {
        kpar *= (P / P0);
    }
    double kperp = kperp_kpar * kpar;

    /*
    double gD = updateGammaD(r_, th_);

    dKrr_dr = (kpar * std::pow(r_, 4) - std::pow(r_, 2) * (r_ - rSun) * (kperp * (rSun - 2 * r0 - 3 * r_)
                + kpar * (r_ + 2 * r0 + rSun)) * std::pow(gD, 2) + kperp * std::pow((r_ - rSun) * gD, 4))
        / ((r_ + r0) * std::pow(std::pow(r_, 2) + std::pow((r_ - rSun) * gD, 2), 2));
    dKrph_dr = ((kpar - kperp) * r_ * gD * (std::pow(r_, 2) * (rSun - r0 - 2 * r_)
                + std::pow((r_ - rSun) * gD, 2) * (r0 + rSun)))
        / ((r_ + r0) * std::pow(std::pow(r_, 2) + std::pow((r_ - rSun) * gD, 2), 2));
    */

    // Angles used in transformation from HMF-aligned to spherical coordinates
    double tanPsi = omega * (r_ - rSun) * std::sin(th_) / vsw;
    double cosPsi = 1.0 / std::sqrt(1 + std::pow(tanPsi, 2));
    double sinPsi = std::sqrt(1 - std::pow(cosPsi, 2));

    return KTensor(kpar * std::pow(cosPsi, 2) + kperp * std::pow(sinPsi, 2),
            kpar * std::pow(sinPsi, 2) + kperp * std::pow(cosPsi, 2),
            (kperp - kpar) * cosPsi * sinPsi,
            kperp);
}

void updateDK() {
    double deltar = 0.01; // AU

    KTensor K1 = updateK(r + deltar, th);

    // Compute the derivatives numerically, as Strauss does
    dKrr_dr = (K1.rr - K.rr) / deltar;
    dKrph_dr = (K1.rph - K.rph) / deltar;
}

// Updates drift velocity
void updateVd() {
    // Gradient and curvature drifts
    double vdCoeff = 2 * P * v * r / (3 * A * qSign * B0 * std::pow(r0 * (1 + gammaD*gammaD), 2))
        * heaviside(th);
    
    vdr = vdCoeff * (-gammaD / std::tan(th));
    vdth = vdCoeff * (2 * gammaD + std::pow(gammaD, 3));
    vdph = vdCoeff * (std::pow(gammaD, 2) / std::tan(th));

    // Account for HCS drifts when the particle is close to the current sheet
    double d = std::abs(r * std::cos(th));
    // This is always finite since the call to heaviside() is not in updateB()
    double rL = P / std::sqrt(B.r*B.r + B.th*B.th + B.ph*B.ph);
    if (d <= 2 * rL) {
        double tanPsi = omega * (r - rSun) * std::sin(th) / vsw;
        double cosPsi = 1 / std::sqrt(1 + std::pow(tanPsi, 2));
        double sinPsi = std::sqrt(1 - std::pow(cosPsi, 2));

        vdr += A * qSign * (0.457 - 0.421 * d / rL + 0.0915 * std::pow(d / rL, 2)) * v * sinPsi;
        vdph += A * qSign * (0.457 - 0.421 * d / rL + 0.0915 * std::pow(d / rL, 2)) * v * cosPsi;
    }

    // Drift reduction factor
    double fs = std::pow(P / P0, 2) / (1 + std::pow(P / P0, 2));

    vdr *= fs;
    vdth *= fs;
    vdph *= fs;
}

// Updates convenience variables, magnetic field, K and K's partials
void updateVars() {
    updateRig();
    gammaD = updateGammaD(r, th);
    updateV();
    B = updateB(r, th);
    K = updateK(r, th);
    updateDK();
    updateVd();
}

// Sets up the initial conditions for the simulation.  The particles starts at Earth, which has coordinates
// (r, th, ph) = (1, pi/2, 0).
void initialize(double ek_i) {
    s = 0.0;
    r = 1.0;
    th = M_PI / 2;
    ph = 0;
    // Use initial energy of 100 MeV for now
    ek = ek_i;

    // Set up initial variables
    updateVars();
}

// Steps the simulation forward (ie, backwards in time) by ds
Status step() {
    updateVars();

    // Partials for incrementing the coordinates and energy
    double dr_ds = 2 / r * K.rr + dKrr_dr - vsw - vdr;
    double dr_dWr = std::sqrt(2 * K.rr - 2 * std::pow(K.rph, 2) / K.phph);
    double dr_dWph = std::sqrt(2 / K.phph) * K.rph;

    double dth_ds = 1 / (r*r * std::tan(th)) * K.thth - vdth / r;
    double dth_dWth = std::sqrt(2 * K.thth) / r;

    double dph_ds = 1 / (r*r * std::sin(th)) * K.rph + 1 / (r * std::sin(th)) * dKrph_dr
        - vdph / (r * std::sin(th));
    double dph_dWph = std::sqrt(2 * K.phph) / (r * std::sin(th));

    double dek_ds = 2 * vsw * ek / (3 * r) * ((ek + 2 * m) / (ek + m)); // Last term is Gamma

    // Generate the Wiener terms
    double dWr = distribution(generator) * std::sqrt(ds);
    double dWth = distribution(generator) * std::sqrt(ds);
    double dWph = distribution(generator) * std::sqrt(ds);

    // Increment the coordinates!
    r += dr_ds * ds + dr_dWr * dWr + dr_dWph * dWph;
    th += dth_ds * ds + dth_dWth * dWth;
    ph += dph_ds * ds + dph_dWph * dWph;
    ek += dek_ds * ds;
    s += ds;

    // Renormalize the coordinates.  TODO: this is wrong!  It will fail for big jumps.
    if (r < 0) {
        r = -r;
        th = M_PI / 2 - th;
        ph = -ph;
    }
    while (th > M_PI) {
        th = 2 * M_PI - th;
        ph = ph - M_PI;
    }
    while (th < 0) {
        th = -th;
        ph = ph + M_PI;
    }
    while (ph > 2 * M_PI) {
        ph = ph - 2 * M_PI;
    }
    while (ph < 0) {
        ph = ph + 2 * M_PI;
    }

    if (r >= rHP) {
        return Status::Heliopause;
    } else if (r <= rSun) {
        return Status::Sun;
    } else {
        return Status::Running;
    }
}

// Returns a CSV string containing the particle's current coordinates, kinetic energy and elapsed time
std::string stateToString() {
    return std::to_string(r) + "," + std::to_string(th) + "," + std::to_string(ph) + "," + std::to_string(ek)
        + "," + std::to_string(s);
}

////////////////////

// argv[1]: initial energy (GeV)
// argv[2]: number of runs ending at the heliopause to simulate
// argv[3]: output file name.  The file will be a csv located in the rundata directory.
int main(int argc, char *argv[]) {
    // Indicates that only the final point needs to be stored
    //bool endOnly = true;

    // Make sure both command line arguments were provided
    if (argc == 4) {
        // Initial energy (GeV)
        double ek_i = std::stod(argv[1]);
        // Number of runs to perform
        double runs = std::stod(argv[2]);
        // Output file name
        std::string fName(argv[3]);
        // Directory to which run data will be written
        const std::string runDir("rundata");

        std::cout << "Tracing " << runs << " particles detected with energy " << ek_i << " GeV at Earth back"
            " to the heliopause..." << std::endl;

        // String containing run data.  This will be written to a CSV file.
        std::string runsString("# Run exit points.  Columns are r (AU), th (rad), ph (rad), ek (GeV), s "
                "(s).");

        // Measure how long the simulation takes.  Store time since program started.
        std::clock_t start = std::clock();

        // Variable for tracking what percent of the simulation has finished
        int percentDone = 0;
        int runsTo10Percent = runs / 10;

        // Generate the runs
        for (int successes = 0; successes < runs; ) {
            initialize(ek_i);

            while (true) {
                if (step() == Status::Sun) {
                    // Trajectory ended in the sun.  Run another trajectory and don't increment i.
                    break;
                } else if (step() == Status::Heliopause) {
                    // Particle ended at the heliopause!  Increment the event counter, add the end point to
                    // the list and run another trajectory.
                    successes++;
                    runsString += "\n" + stateToString();

                    // Print percent indicator
                    runsTo10Percent--;

                    if (runsTo10Percent <= 0) {
                        percentDone += 10;
                        std::cout << percentDone << "% of runs complete" << std::endl;
                        runsTo10Percent = runs / 10;
                    }

                    break;
                }
            }
        }

        // Measure new time since program started, subtract and convert to seconds or minutes
        double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
        if (duration < 60) {
            std::cout << "Simulation completed in " << duration << " seconds" << std::endl;
        } else {
            std::cout << "Simulation completed in " << duration / 60 << " minutes" << std::endl;
        }

        // Write run data to a CSV
        std::ofstream writer(runDir + "/" + fName + ".csv");

        if (writer.is_open()) {
            writer << runsString;
            writer.close();
        }

        std::cout << "Wrote exit point data to " << runDir + "/" + fName + ".csv" << std::endl;

        // Beep!
        std::cout << '\a';

        return 0;
    } else if (argc == 2) {
        // Run some tests
        if (std::string(argv[1]) == "B") {
            r = 1;
            th = M_PI / 2;
            ph = 0;
            BField Bf = updateB(r, th);
            std::cout << "(r, th, ph) = (" << r << ", " << th << ", " << ph << "): Br = " << Bf.r
                << ", Bph = " << Bf.ph << std::endl;

            r = 10;
            th = M_PI / 3;
            ph = M_PI / 3;
            Bf = updateB(r, th);
            std::cout << "(r, th, ph) = (" << r << ", " << th << ", " << ph << "): Br = " << Bf.r
                << ", Bph = " << Bf.ph << std::endl;

            r = 10;
            th = -M_PI / 3;
            ph = 3 * M_PI / 2;
            Bf = updateB(r, th);
            std::cout << "(r, th, ph) = (" << r << ", " << th << ", " << ph << "): Br = " << Bf.r
                << ", Bph = " << Bf.ph << std::endl;
        }
    } else {
        std::cout << "Need three command line arguments: the initial energy (in GeV), number of runs to "
            "perform and output file name." << std::endl;
        return 1;
    }
}
