#pragma once
#include <Eigen/Dense>

struct AircraftParams {
  double m = 1111.0; // kg (Cessna 172R Skyhawk R)
  double S = 17.08; // m^2 ()
  double cref = 1.57; // m
  double bref = 11.0; // m
  double rho = 1.225; // kg/m^3
  Eigen::Matrix3d inertiaMat;
  double g = 9.80665;
  bool ruddered;
};

class Aircraft {
public:
  Aircraft(const AircraftParams &p);
  // Compute aerodynamic coeff.
  void aeroCoeffs(double alpha, double beta,
                  double delta_a, double delta_e, double delta_r, double delta_tail,
                  double &CL, double &CD, double &CY) const;

  // Compute forces, moments in body frame
  void computeForcesMoments(const Eigen::VectorXd &x, const Eigen::VectorXd &u,
                            Eigen::Vector3d &F_body, Eigen::Vector3d &M_body) const;

  const AircraftParams params;
};
