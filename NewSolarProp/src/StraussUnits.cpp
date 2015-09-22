#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#define DEBUG false
#define DEBUG_RANDOMS false

#define PI 3.141592653589793

//// Constants

// Initial coordinates
const double r0 = 1.0; // AU
const double th0 = PI / 2.0; // rad
const double ph0 =  0.0; // rad
// Jupiter initial phi position
const double ph0Jup = PI;
// Used for numerically differentiating K
const double deltar = 0.01; // AU
// Parallel mean free path constant
const double lambda0 = 0.15; // AU
// Reference distance in mean free path
const double rRefLambda = 1.0; // AU
// k_perp / k_parallel
const double kperp_kpar = 0.01;
// Reference rigidity
const double P0 = 1.0; // GV
// Polarity of HMF
const double Ac = -1;
// Distance from sun to heliopause
const double rHP = 140; // AU
// Sun's radius
const double rSun = 0.005; // AU
// Inner boundary used in Strauss' code
const double rInner = rSun; // AU
// Particle mass
const double m = 0.000511;
// Sign of particle's charge
const int qSign = -1;

// Program time
const double protime = 1.496e8 / 400;
// Timestep
const double ds = 0.001;
// Angular velocity of sun
const double Omega1 = 2 * PI / (25.4 * 24 * 3600);
// Solar wind velocity
const double Vsw = 1;
// Reference field strength
const double B0 =  5 * 0.06;

////// Jupiter parameters
const double omegaJup = 2 * PI / (4333.0 * 3600.0 * 24.0) * protime;
const double rJup = 5.2; // AU
const double thJup = PI / 2.0;
const double dphJup = 0.009 * 2.0;
const double dthJup = 0.009 * 2.0;
const double rBeginJup = rJup - 0.0477 * 2.0;
const double rEndJup = rJup + 0.095 * 2.0;

////////////////////////////////////////

/*
int heaviside(double x) {
    if (x > 0) {
        return 1;
    } else if (x < 0) {
        return -1;
    } else {
        return 0;
    }
}
*/

class KTensor {
    public:
        double rr, phph, rph, thth;
        double dKrr_dr, dKphph_dr, dKrph_dr, dKthth_dr;
        double dKrr_dth, dKphph_dth, dKrph_dth, dKthth_dth;

        KTensor() : rr(0), phph(0), rph(0), thth(0) { };

        KTensor(double r, double th, double P) {
            updateK(r, th, P);
        }

        KTensor& updateK(double rc, double thc, double Pc) {
            // Compute elements at shifted points and store elements
            updateKElements(rc + deltar, thc, Pc);

            // Record relevant elements
            double krr_r = rr;
            double krph_r = rph;

            // Compute elements at actual point
            updateKElements(rc, thc, Pc);

            // Compute r derivatives
            dKrr_dr = (krr_r - rr) / deltar;
            dKrph_dr = (krph_r - rph) / deltar;

            // Compute th derivatives
            dKthth_dth = 0;

#if DEBUG
            std::cout << "Krr = " << rr
                << "\n\tKrr2 = " << krr_r
                << "\nKphph = " << phph
                << "\nKrph = " << rph
                << "\n\tKrph2 = " << krph_r
                << "\nKthth = " << thth
                << "\ndKrr_dr = " << dKrr_dr
                << "\n\tdelta_r = " << deltar
                << "\ndKrph_dr = " << dKrph_dr << std::endl;
#endif
            
            return *this;
        };

    private:
        KTensor& updateKElements(double rc, double thc, double Pc) {
            double beta = Pc / sqrt(Pc*Pc + m*m);
            double kpar = (beta * 750 / 3) * lambda0 * (1 + rc / rRefLambda) * (Pc >= P0? Pc / P0 : 1);

            double kperp = kperp_kpar * kpar;

            // Convert to spherical coordinates
            double Omega = Omega1 * protime;
            double tanPsi = Omega * (rc - rSun) * sin(thc) / Vsw;
            double cosPsi = 1 / sqrt(1 + tanPsi*tanPsi);
            double sinPsi = sqrt(1 - cosPsi*cosPsi);

            rr = kpar * cosPsi*cosPsi + kperp * sinPsi*sinPsi;
            phph = kpar * sinPsi*sinPsi + kperp * cosPsi*cosPsi;
            rph = (kperp - kpar) * cosPsi * sinPsi;
            thth = kperp;

            return *this;
        };
};

// Used to indicate the simulation's status after each timestep
enum class Status {Sun, Heliopause, Jupiter, Running};

//// Coordinates, fields and tensors

// Coordinates
double s;
double r;
double th;
double ph;
// Jupiter
double phJup;
// Kinetic energy
double ek; // GeV
// Rigidity
double P; // GV

// Magnetic field
double Bmag;
// Diffusion tensor in spherical coordinates
KTensor K;
// Partial derivatives of K.  TODO: delete!!!
//double dKrr_dr, dKrph_dr;
// Drift velocity components
double vdr, vdth, vdph;

// Gaussian random number generator with mean 0 and standard deviation 1
double seed1 = 975635;
double rand1, rand2;
const double a = pow(7.0, 5);
const double M = pow(2.0, 101) - 1;
const double c = 0.0;

void straussRand() {
    double randHelper1 = fmod(a * seed1 + c, M);
    double randHelper2 = fmod(a * randHelper1 + c, M);

    seed1 = randHelper2;

    double sDev1 = randHelper1 / M;
    double sDev2 = randHelper2 / M;

    rand1 = sqrt(-2.0 * log(sDev1)) * cos(2.0 * PI * sDev2);
    rand2 = sqrt(-2.0 * log(sDev1)) * sin(2.0 * PI * sDev2);

#if DEBUG_RANDOMS
    std::cout << "seed1 = " << seed1 << std::endl;
    std::cout << "r1, r2 = " << rand1 << ", " << rand2 << std::endl;
#endif
}

////////////////////////////////////////

double getP(double e) {
    return sqrt(e * (e + 2 * m));
}

void updateVd() {
    double Omega = Omega1 * protime;

    double gamma = r * Omega * sin(th) / Vsw;
    double beta = P / sqrt(P*P + m*m);
    double vdCoeff = 2.0/(3.0 * Ac * qSign * B0) * P * beta * r / pow(1 + gamma*gamma, 2);
    double vdSign = 1;

    // Heaviside function comes in here
    if (th > PI / 2) {
        vdSign = -1;
    } else if (th == PI / 2) {
        vdSign = 0;
    }

    vdr = vdSign * vdCoeff * (-gamma / tan(th));
    vdth = vdSign * vdCoeff * (2 + gamma*gamma) * gamma;
    vdph = vdSign * vdCoeff * gamma*gamma / tan(th);

    double d = fabs(r * cos(th));
    Bmag = B0 / (sqrt(1 + pow(1 - rSun, 2)) * r*r * 1/sqrt(1 + pow(Omega*(r-rSun)*sin(th)/Vsw, 2)));
    double rL = P/750 * 1 / Bmag;

    if (d <= 2 * rL) {
        vdr += Ac * qSign * (0.457 - 0.412 * d / rL + 0.0915 * pow(d / rL, 2)) * beta * 750;
    }

    // Drift reduction factor
    double fs = 10 * pow(P / P0, 2) / (1 + 10 * pow(P / P0, 2));

    vdr *= fs;
    vdth *= fs;
    vdph *= fs;

#if DEBUG
    std::cout << "vdCoeff = " << vdCoeff
        << "\nvdr = " << vdr
        << "\nvdth = " << vdth
        << "\nvdph = " << vdph << std::endl;
#endif
}

void updateVars() {
    P = getP(ek);
    K.updateK(r, th, P);
    updateVd();
}

void initialize(double ek0) {
    s = 0;
    r = r0;
    th = th0;
    ph = ph0;
    ek = ek0;
    // Reset Jupiter
    phJup = ph0Jup;
}

Status step() {
    updateVars();

    double dr_ds = 2 / r * K.rr + K.dKrr_dr - (Vsw + vdr);
    double dr_dWr = sqrt(2 * K.rr - 2 * K.rph*K.rph / K.phph);
    double dr_dWph = K.rph * sqrt(2 / K.phph);

    double dth_ds = 1 / (r*r) * K.dKthth_dth + 1 / (r*r * tan(th)) * K.thth - vdth / r;
    double dth_dWth = sqrt(2 * K.thth) / r;

    double dph_ds = 1 / (r*r * sin(th)) * K.rph + 1 / (r * sin(th)) * K.dKrph_dr - vdph / (r * sin(th));
    double dph_dWph = sqrt(2 * K.phph) / (r * sin(th));

    double Gamma = (ek + 2 * m) / (ek + m);
    double dek_ds = 2 * Vsw / (3 * r) * Gamma * ek;

    // Generate the Wiener terms
    /*
    double dWr = generator() * sqrt(ds);
    double dWth = generator() * sqrt(ds);
    double dWph = generator() * sqrt(ds);
    */
    straussRand();
    double dWr = sqrt(ds) * rand1;
    double dWph = sqrt(ds) * rand2;
    straussRand();
    double dWth = sqrt(ds) * rand2;

#if DEBUG
    std::cout << "dr_ds = " << dr_ds << std::endl;
    std::cout << "dr_dWr = " << dr_dWr << std::endl;
    std::cout << "dr_dWph = " << dr_dWph << std::endl;
    std::cout << "dth_ds = " << dth_ds << std::endl;
    std::cout << "dth_dWth = " << dth_dWth << std::endl;
    std::cout << "dph_ds = " << dph_ds << std::endl;
    std::cout << "dph_dWph = " << dph_dWph << std::endl;
    std::cout << "dE_ds = " << dek_ds << std::endl;
#endif

    // Move the particle
    r += dr_ds * ds + dr_dWr * dWr + dr_dWph * dWph;
    th += dth_ds * ds + dth_dWth * dWth;
    ph += dph_ds * ds + dph_dWph * dWph;
    ek += dek_ds * ds;
    // Move Jupiter
    phJup -= omegaJup * ds;
    // Advance the backwards clock!
    s += ds;

    // Renormalize coordinates
    while (th > PI) {
        th = 2 * PI - th;
        ph = ph - PI;
    }

    while (th < 0) {
        th = -th;
        ph += PI;
    }

    while (ph > 2 * PI) {
        ph -= 2 * PI;
    }

    while (ph < 0) {
        ph += 2 * PI;
    }


    if (r > rHP) {
        return Status::Heliopause;
    } else if (r < rInner) { // Just copying Strauss here...
        return Status::Sun;
    } else if (r > rBeginJup && r < rEndJup
            && ph > phJup - dphJup && ph < phJup + dphJup
            && th > thJup - dthJup && th < thJup + dthJup) {
        return Status::Jupiter;
    } else {
        return Status::Running;
    }
}

// Returns a CSV string containing the particle's current coordinates, kinetic energy and elapsed time
std::string stateToString() {
    return std::to_string(r) + "," + std::to_string(th) + "," + std::to_string(ph) + "," + std::to_string(ek)
        + "," + std::to_string(s * 4.3287 * 86400);
}

void printState() {
    std::cout << "s = " << s * 4.3287 * 86400 << "\nr, th, ph, E = " << r << "," << th << ","
        << ph << "," << ek << std::endl;
}

////////////////////

// argv[1]: initial energy (GeV)
// argv[2]: number of runs ending at the heliopause to simulate
// argv[3]: output file name.  The file will be a csv located in the rundata directory.
int main(int argc, char *argv[]) {
    std::cout.precision(std::numeric_limits<double>::max_digits10);

    // Make sure both command line arguments were provided
    if (argc == 4) {
        // Initial energy (GeV)
        double ek0 = std::stod(argv[1]);
        // Number of runs to perform
        double runs = std::stod(argv[2]);
        // Output file name
        std::string fName(argv[3]);
        // Directory to which run data will be written
        const std::string runDir("rundata");

        std::cout << "Tracing " << runs << " particles detected with energy " << ek0 << " GeV at Earth back"
            " to the heliopause..." << std::endl;

        // String containing run data.  This will be written to a CSV file.
        std::string runsString("# Run exit points.  Columns are r (AU), th (rad), ph (rad), ek (GeV), s "
                "(s).");

        // Measure how long the simulation takes.  Store time since program started.
        std::clock_t start = std::clock();

        // Count particles that hit Jupiter and the sun
        int jupiterCount = 0, sunCount = 0;

        // Variable for tracking what percent of the simulation has finished
        int percentDone = 0;
        int runsTo1Percent = runs / 100;

        // Variable for storing current simulation's status
        Status stepStatus;

        // Generate the runs
        for (int successes = 0; successes < runs; ) {
            // Reinitialize simulation variables
            initialize(ek0);

            while (true) {
                // Step the simulation
#if DEBUG
                printState();
#endif
                stepStatus = step();
#if DEBUG
                std::cout << std::endl;
#endif

                if (stepStatus == Status::Jupiter) {
                    // Trajectory ended in Jupiter.  Run another trajectory and don't increment i.
                    jupiterCount++;
                    break;
                } else if (stepStatus == Status::Sun) {
                    // Trajectory ended in the sun.  Run another trajectory and don't increment i.
                    sunCount++;
                    break;
                } else if (stepStatus == Status::Heliopause) {
                    // Particle ended at the heliopause!  Increment the event counter, add the end point to
                    // the list and run another trajectory.
                    successes++;
                    runsString += "\n" + stateToString();

                    // Print percent indicator
                    runsTo1Percent--;

                    if (runsTo1Percent <= 0) {
                        percentDone += 1;
                        std::cout << percentDone << "% of runs complete" << std::endl;
                        runsTo1Percent = runs / 100;
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

        // Output sun and Jupiter counts
        std::cout << sunCount << " trajectories ended in the sun, " << jupiterCount << " in Jupiter"
            << std::endl;

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
    } else {
        std::cout << "Need three command line arguments: the initial energy (in GeV), number of runs to "
            "perform and output file name." << std::endl;
        return 1;
    }
}


