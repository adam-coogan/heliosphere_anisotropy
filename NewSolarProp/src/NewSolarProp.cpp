/*
 * Rewriting everything again, because I could have some ill-thought-out code in SolarProp.cpp...
 */
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

//// Constants

// Initial coordinates
const double r0 = 1.0; // AU
const double th0 = M_PI / 2; // rad
const double ph0 = 0; // rad
// Used for numerically differentiating K
const double deltar = 0.01; // AU
// Used for numerically differentiating K
const double deltath = M_PI / 180; // rad
// Timestep
const double ds = 345.6; // s
// Parallel mean free path constant
//const double lambda0 = 0.15; // AU, a
//const double lambda0 = 0.15; // AU, b
//const double lambda0 = 0.05; // AU, c
//const double lambda0 = 0.10; // AU, d
// k_perp / k_parallel
const double kperp_kpar = 0.01; // x
//const double kperp_kpar = 0.02; // y
//const double kperp_kpar = 0.10; // z
// Reference rigidity
const double P0 = 1.0; // GV
// Polarity of HMF
const double Ac = 1;
// Reference field strength
const double B0 = 3.535 * pow(10, -9) * 4.485 * std::pow(10, 10); // GV / AU
// Field strength at Earth
const double Bearth = 5.180 * pow(10, -9) * 4.485 * pow(10, 10); // GV / AU
// Angular velocity of sun
const double Omega = 2 * M_PI / (25.4 * 24 * 3600); // rad / s
// Speed of light
const double speedOfLight = 0.002004; // AU / s
// Distance from sun to heliopause
const double rHP = 120; // AU
// Inner boundary used in Strauss' code
const double rInner = 0.01; // AU
// Sun's radius
const double rSun = 0.005; // AU
// Solar wind velocity
const double Vsw = 400 * 6.685 * std::pow(10, -9); // AU / s
// Particle mass
const double m = 0.000511;
// Sign of particle's charge
const int qSign = -1;

////////////////////////////////////////

int heaviside(double x) {
    if (x > 0) {
        return 1;
    } else if (x < 0) {
        return -1;
    } else {
        return 0;
    }
}

class BField {
    public:
        // Components
        double br, bph;

        BField() : br(0), bph(0) { };

        double magnitude() const { return sqrt(br*br + bph*bph); };

        BField& updateB(double rc, double thc) {
            double bFact = Ac * B0 * r0*r0 / (rc*rc);// * (1 - 2 * heaviside(thc - M_PI / 2));

            br = bFact;
            bph = -bFact * (rc - rSun) * Omega / Vsw * sin(thc);

            return *this;
        };
};

class KTensor {
    public:
        double rr, phph, rph, thth;
        double dKrr_dr, dKphph_dr, dKrph_dr, dKthth_dr;
        double dKrr_dth, dKphph_dth, dKrph_dth, dKthth_dth;

        KTensor() : rr(0), phph(0), rph(0), thth(0) { };

        KTensor& updateK(double rc, double thc, double vc, double Pc) {
            // Compute elements at shifted points and store elements
            updateKElements(rc + deltar, thc, vc, Pc);
            double krr_r = rr;
            double kphph_r = phph;
            double krph_r = rph;
            double kthth_r = thth;

            // Make sure theta is in [0, pi].  No need to worry about phi here.
            if (thc + deltath > M_PI) {
                updateKElements(rc, thc - deltath, vc, Pc);
            } else {
                updateKElements(rc, thc + deltath, vc, Pc);
            }

            double krr_th = rr;
            double kphph_th = phph;
            double krph_th = rph;
            double kthth_th = thth;

            // Compute elements at actual point
            updateKElements(rc, thc, vc, Pc);

            // Compute r derivatives
            dKrr_dr = (krr_r - rr) / deltar;
            dKphph_dr = (kphph_r - phph) / deltar;
            dKrph_dr = (krph_r - rph) / deltar;
            dKthth_dr = (kthth_r - thth) / deltar;

            // Compute th derivatives
            dKrr_dth = (krr_th - rr) / deltath;
            dKphph_dth = (kphph_th - phph) / deltath;
            dKrph_dth = (krph_th - rph) / deltath;
            dKthth_dth = (kthth_th - thth) / deltath;

            return *this;
        };

    private:
        KTensor& updateKElements(double rc, double thc, double vc, double Pc) {
            // a
            double kpar = (vc / 3) * 0.15 * (1 + rc / r0) * (Pc >= P0? Pc / P0 : 1);

            // b
            //double kpar = (speedOfLight / 3) * 0.10 * (1 + rc / r0) * (Pc >= P0? Pc / P0 : 1);

            /* // c
            double bMagnitude = B0 / Vsw * pow(r0 / rc, 2)
                * sqrt(pow(Omega * (rc - rSun) * sin(thc), 2) + Vsw*Vsw);
            double kpar = (vc / 3) * 0.05 * Bearth / bMagnitude * (Pc >= P0? Pc / P0 : 1);
            */

            /* // d
            double bMagnitude = B0 / Vsw * pow(r0 / rc, 2)
                * sqrt(pow(Omega * (rc - rSun) * sin(thc), 2) + Vsw*Vsw);
            double kpar = (vc / 3) * 0.15 * Bearth / bMagnitude * Pc / P0;
            */

            double kperp = kperp_kpar * kpar;

            // Convert to spherical coordinates
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
enum class Status {Sun, Heliopause, Running};

//// Coordinates, fields and tensors

// Coordinates
double s;
double r;
double th;
double ph;
// Kinetic energy
double ek; // GeV
// Rigidity
double P; // GV
// Velocity
double v; // AU / s

// Magnetic field
BField B;
// Diffusion tensor in spherical coordinates
KTensor K;
// Partial derivatives of K
double dKrr_dr, dKrph_dr;
// Drift velocity components
double vdr, vdth, vdph;

// Gaussian random number generator with mean 0 and standard deviation 1
/* std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 1.0); */
typedef boost::mt19937 MT19937;
typedef boost::normal_distribution<double> NDistribution;
MT19937 engine(static_cast<unsigned int>(std::time(0)));
NDistribution distribution(0, 1);
boost::variate_generator<MT19937, NDistribution> generator(engine, distribution);

////////////////////////////////////////

double getV(double e) {
    return speedOfLight * sqrt(e * (e + 2 * m)) / (e + m);
}

double getP(double e) {
    return sqrt(e * (e + 2 * m));
}

void updateVd() {
    double gamma = r * Omega * sin(th) / Vsw;
    double vdCoeff = 2 * P * v * r / (3 * qSign * B0 * r0*r0 * Ac * pow(1 + gamma * gamma, 2));

    // Heaviside function comes in here
    if (th > M_PI / 2) {
        vdCoeff *= -1;
    } else if (th == M_PI / 2) {
        vdCoeff = 0;
    }

    vdr = vdCoeff * (-gamma / tan(th));
    vdth = vdCoeff * (2 + gamma*gamma) * gamma;
    vdph = gamma*gamma / tan(th);

    double d = fabs(r * cos(th));
    double rL = P / B.magnitude();

    if (d <= 2 * rL) {
        double tanPsi = Omega * (r - rSun) * std::sin(th) / Vsw;
        double cosPsi = 1 / std::sqrt(1 + std::pow(tanPsi, 2));
        double sinPsi = std::sqrt(1 - std::pow(cosPsi, 2));

        double vdhcs = Ac * qSign * (0.457 - 0.412 * d / rL + 0.0915 * pow(d / rL, 2)) * v;
        vdr += sinPsi * vdhcs;
        vdph += cosPsi * vdhcs;
    }

    //double fs = 1; // Lots of people don't put in a drift reduction factor
    double fs = 10 * pow(P / P0, 2) / (1 + 10 * pow(P / P0, 2)); // THIS IS CORRECT!

    vdr *= fs;
    vdth *= fs;
    vdph *= fs;
}

void updateVars() {
    v = getV(ek);
    P = getP(ek);

    B.updateB(r, th);
    K.updateK(r, th, v, P);

    updateVd();
}

void initialize(double ek_i) {
    s = 0;
    r = r0;
    th = th0;
    ph = ph0;
    ek = ek_i;

    updateVars();
}

Status step() {
    updateVars();

    // TODO: try different sign combinations!!!
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
    double dWr = generator() * sqrt(ds);
    double dWth = generator() * sqrt(ds);
    double dWph = generator() * sqrt(ds);

    r += dr_ds * ds + dr_dWr * dWr + dr_dWph * dWph;
    th += dth_ds * ds + dth_dWth * dWth;
    ph += dph_ds * ds + dph_dWph * dWph;
    ek += dek_ds * ds;
    s += ds;

    // Renormalize coordinates
    while (th > M_PI) {
        th = 2 * M_PI - th;
        ph = ph - M_PI;
    }

    while (th < 0) {
        th = -th;
        ph += M_PI;
    }

    while (ph > 2 * M_PI) {
        ph -= 2 * M_PI;
    }

    while (ph < 0) {
        ph += 2 * M_PI;
    }


    if (r > rHP) {
        return Status::Heliopause;
    } else if (r < rInner) { // Just copying Strauss here...
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
        int runsTo1Percent = runs / 100;

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


