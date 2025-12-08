#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <filesystem>

#include "Aircraft.h"
#include "Dynamics.h"
#include "Controllers.h"
#include "Utils.h"

using namespace Eigen;

void runSimulation(bool rudderedConfig, const VectorXd& xtrim, const VectorXd& utrim, 
                   double dt, double Tsim)
{
    int N = int(Tsim / dt);

    // Aircraft Parameters
    AircraftParams params;
    params.ruddered = rudderedConfig;
    params.inertiaMat <<
        197.2026, 0.0,       146.0218,
        0.0,      1808.7634, 0.0,
        146.0218, 0.0,       1611.5609;

    Aircraft ac(params);
    Dynamics dyn(ac);

    // Compute linear model A,B
    MatrixXd A, B;
    dyn.linearizeAt(xtrim, utrim, A, B, 1e-6);

    // Compute LQR gain K
    MatrixXd Q = MatrixXd::Zero(12,12);
    Q(9,9) = 0.1; Q(10,10) = 0.1; Q(11,11) = 0.1;

    Eigen::Vector4d diag_vals;
    diag_vals << 5.25, 4.19, 10.48, 4.0;
    MatrixXd R = diag_vals.asDiagonal() * 55.0;

    MatrixXd K = Controllers::computeLQRGainFromContinuous(A, B, Q, R, dt);

    // Save Linear Model
    if (rudderedConfig) {
        Utils::writeCSV("output/A_ruddered.csv", A);
        Utils::writeCSV("output/B_ruddered.csv", B);
        Utils::writeCSV("output/K_ruddered.csv", K);
    } else {
        Utils::writeCSV("output/A_rudderless.csv", A);
        Utils::writeCSV("output/B_rudderless.csv", B);
        Utils::writeCSV("output/K_rudderless.csv", K);
    }

    // Run Nonlinear Simulation
    MatrixXd hist(N+1, 16);
    VectorXd x = xtrim;
    VectorXd u = utrim;
    hist.row(0) << x.transpose(), u.transpose();

    for (int k = 0; k < N; ++k) {
        VectorXd ufb = -K * (x - xtrim) + utrim;

        // Input limits
        ufb(0) = Utils::clamp(ufb(0), -20*M_PI/180, 20*M_PI/180);
        ufb(1) = Utils::clamp(ufb(1), -28*M_PI/180, 23*M_PI/180);

        if (!rudderedConfig)
            ufb(2) = Utils::clamp(ufb(2), -90*M_PI/180, 90*M_PI/180);
        else
            ufb(2) = Utils::clamp(ufb(2), -17.7*M_PI/180, 17.7*M_PI/180);

        ufb(3) = Utils::clamp(ufb(3), 0.0, 1.0);

        x = dyn.rk4Step(x, ufb, dt);
        hist.row(k+1) << x.transpose(), ufb.transpose();
    }

    if (rudderedConfig)
        Utils::writeCSV("output/ruddered_sim_history.csv", hist);
    else
        Utils::writeCSV("output/rudderless_sim_history.csv", hist);

    std::cout << (rudderedConfig ? 
        "Ruddered simulation complete.\n" :
        "Rudderless simulation complete.\n");
}

int main() {
    std::filesystem::create_directories("output");
    double dt = 0.01;
    double Tsim = 100.0;

    // Rudderless trim
    VectorXd xtrim_rl = VectorXd::Zero(12);
    VectorXd utrim_rl = VectorXd::Zero(4);
    xtrim_rl(0) = 68.0;
    xtrim_rl(2) = 0.0224;
    xtrim_rl(7) = 0.0003;
    utrim_rl(1) = -0.046;
    utrim_rl(3) = 0.2229;

    // Ruddered trim
    VectorXd xtrim_r = xtrim_rl;
    xtrim_r(2) = 0.0197;
    VectorXd utrim_r = VectorXd::Zero(4);
    utrim_r(1) = -0.027;
    utrim_r(3) = 0.2333;

    runSimulation(false, xtrim_rl, utrim_rl, dt, Tsim);
    runSimulation(true,  xtrim_r,  utrim_r,  dt, Tsim);

    std::cout << "Both simulations complete.\n";
    return 0;
}
