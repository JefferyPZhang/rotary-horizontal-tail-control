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
                   double dt, double Tsim){
    
    // Num iterations
    int N = int(Tsim / dt);

    // Aircraft params
    AircraftParams params;
    params.ruddered = rudderedConfig;
    params.inertiaMat << 
        197.2026, 0.0,     146.0218,
        0.0,     1808.7634, 0.0,
        146.0218, 0.0,     1611.5609;

    // Create aircraft + dynamics
    Aircraft ac(params);
    Dynamics dyn(ac);

    // Linearize
    MatrixXd A, B;
    dyn.linearizeAt(xtrim, utrim, A, B, 1e-6);

    // LQR gain
    MatrixXd Q = MatrixXd::Zero(12,12);
    Q(9,9) = 0.1;
    Q(10,10) = 0.1;
    Q(11,11) = 0.1;
  
    Eigen::Vector4d diag_vals;
    diag_vals << 5.25, 4.19, 10.48, 4.0;
    Eigen::Matrix4d R = diag_vals.asDiagonal() * 55.0;
    MatrixXd K = Controllers::computeLQRGainFromContinuous(A, B, Q, R, dt);

    // Simulation Loop
    MatrixXd hist(N+1, 16);
    VectorXd x = xtrim;
    VectorXd u = utrim;

    hist.row(0) << x.transpose(), u.transpose();

    for(int k=0; k < N; ++k){
        VectorXd ufb = -K * (x - xtrim) + utrim;

        // Input constraints
        ufb(0) = Utils::clamp(ufb(0), -20.0*M_PI/180.0, 20.0*M_PI/180.0); // aileron
        ufb(1) = Utils::clamp(ufb(1), -28.0*M_PI/180.0, 23.0*M_PI/180.0); // elevator
        if (!rudderedConfig)
            ufb(2) = Utils::clamp(ufb(2), -90.0*M_PI/180.0, 90.0*M_PI/180.0); // tail
        else
            ufb(2) = Utils::clamp(ufb(2), -17.7*M_PI/180.0, 17.7*M_PI/180.0); // rudder
        ufb(3) = Utils::clamp(ufb(3), 0.0, 1.0); // throttle

        x = dyn.rk4Step(x, ufb, dt);
        hist.row(k+1) << x.transpose(), ufb.transpose();
    }
    // 6. Save results
    if (rudderedConfig)
        Utils::writeCSV("output/ruddered_sim_history.csv", hist);
    else
        Utils::writeCSV("output/rudderless_sim_history.csv", hist);

    std::cout << (rudderedConfig ? 
                 "Ruddered simulation complete.\n" : 
                 "Rudderless simulation complete.\n");
}

int main(){
    std::filesystem::create_directories("output");
    double dt = 0.01;
    double Tsim = 100.0;

    // Set rudderless trim state
    VectorXd xtrim_rl = VectorXd::Zero(12);
    VectorXd utrim_rl = VectorXd::Zero(4);
    xtrim_rl(0) = 68.0;      // u
    xtrim_rl(2) = 0.0224;    // w
    xtrim_rl(7) = 0.0003;    // pitch angle
    utrim_rl(1) = -0.046;    // elevator
    utrim_rl(3) = 0.2229;    // throttle

    // Set ruddered trim state
    VectorXd xtrim_r = xtrim_rl;
    xtrim_r(2) = 0.0197;
    VectorXd utrim_r = VectorXd::Zero(4);
    utrim_r(1) = -0.027;     // elevator
    utrim_r(3) = 0.2333;     // throttle

    runSimulation(false, xtrim_rl, utrim_rl, dt, Tsim);  // rudderless
    runSimulation(true,  xtrim_r,  utrim_r,  dt, Tsim);  // ruddered
    std::cout << "Both simulations complete.\n";
    return 0;
}
