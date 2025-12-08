#define _USE_MATH_DEFINES
#include "Aircraft.h"
#include <cmath>
#include <algorithm>

Aircraft::Aircraft(const AircraftParams &p) : params(p) {}

void Aircraft::aeroCoeffs(double alpha, double beta,
                         double delta_a, double delta_e, double delta_r, double delta_tail,
                         double &CL, double &CD, double &CY) const {
  // Polynomial fits from the paper (Singhasenee)
  if (params.ruddered) {
    // Ruddered (Eqns 51-53)
    CL = -1.979 * alpha*alpha*alpha - 0.1339 * alpha*alpha + 5.1882 * alpha + 0.1514;
    CD = 1.1827 * alpha*alpha + 0.0582 * alpha + 0.0089;
    CY = -0.17299 * beta + 0.1293 * delta_r;
  } else {
    // Rudderless (Eqns 54-56)
    CL = -1.9779 * alpha*alpha*alpha - 0.1339 * alpha*alpha + 5.1846 * alpha + 0.1512;
    CD = 1.1821 * alpha*alpha + 0.0581 * alpha + 0.0085;
    CY = (0.3578 * std::cos(2.0*delta_tail) - 0.3578) * beta + 0.0096 * delta_e * std::sin(delta_tail);
  }
}

void Aircraft::computeForcesMoments(const Eigen::VectorXd &x, const Eigen::VectorXd &u,
                                    Eigen::Vector3d &F_body, Eigen::Vector3d &M_body) const {
  // x: [u v w p q r phi theta psi PN PE PD]
  // u: [delta_a, delta_e, delta_r_or_tail, delta_t]
  const double u_b = x(0), v_b = x(1), w_b = x(2);
  const double p = x(3), q = x(4), r = x(5);
  const double phi = x(6), theta = x(7), psi = x(8);

  // Magnitude of airspeed
  const double V = std::sqrt(u_b*u_b + v_b*v_b + w_b*w_b);
  double alpha = 0.0, beta = 0.0;
  if (V > 1e-8) {
    alpha = std::atan2(w_b, u_b); // small-angle
    beta  = std::asin(std::clamp(v_b / V, -1.0, 1.0));
  }

  // Controls
  const double delta_a = u(0);
  const double delta_e = u(1);
  const double delta3   = u(2); // either delta_r or delta_tail
  const double delta_t  = u(3);

  const double delta_r = params.ruddered ? delta3 : 0.0;
  const double delta_tail = params.ruddered ? 0.0 : delta3;

  // Aerodynamic coefficients in stability axes
  double CL, CD, CY;
  aeroCoeffs(alpha, beta, delta_a, delta_e, delta_r, delta_tail, CL, CD, CY);

  // Dynamic pressure and reference area
  const double qbar = 0.5 * params.rho * V * V;
  const double qS   = qbar * params.S;

  // Forces in stability axes
  Eigen::Vector3d F_stab;
  F_stab(0) = -CD * qS;
  F_stab(1) =  CY * qS;
  F_stab(2) = -CL * qS;

  // Rotate stability -> body using rotation about body-y by -alpha:
  const double ca = std::cos(alpha);
  const double sa = std::sin(alpha);
  Eigen::Matrix3d R_s2b;
  R_s2b <<  ca,  0.0, -sa,
            0.0, 1.0,  0.0,
            sa,  0.0,  ca;

  Eigen::Vector3d F_b = R_s2b * F_stab;

  // Gravity in body frame
  double cphi = std::cos(phi), sphi = std::sin(phi);
  double cth  = std::cos(theta), sth = std::sin(theta);
  double cpsi = std::cos(psi), spsi = std::sin(psi);

  // Rotation from body
  Eigen::Matrix3d R_b2i;
  R_b2i << cth*cpsi, sphi*sth*cpsi - cphi*spsi, cphi*sth*cpsi + sphi*spsi,
           cth*spsi, sphi*sth*spsi + cphi*cpsi, cphi*sth*spsi - sphi*cpsi,
           -sth,     sphi*cth,                  cphi*cth;
  Eigen::Matrix3d R_i2b = R_b2i.transpose();

  Eigen::Vector3d grav_i;
  grav_i << 0.0, 0.0, params.m * params.g;
  Eigen::Vector3d Fg_b = R_i2b * grav_i;

  // Add gravity to body forces
  F_b += Fg_b;

  // Add propulsive thrust along body x-axis (simple linear mapping)
  double Fthrust = 0.25 * params.m * params.g * delta_t;
  F_b(0) += Fthrust;

  // Final body forces output
  F_body = F_b;

  // Moments
  Eigen::Vector3d etaBar;
  Eigen::Matrix3d dC_domega = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d dC_du = Eigen::Matrix3d::Zero();

  if (params.ruddered) {
    // Ruddered coefficients
    etaBar(0) = -0.1422 * (beta*beta*beta) - 0.0112 * beta;
    etaBar(1) =  1.3246 * alpha*alpha*alpha - 1.9218 * alpha*alpha - 2.069 * alpha - 0.0461;
    etaBar(2) =  0.1193 * beta;

    dC_domega << -0.54597,     0.0,       0.031232,
                   0.0,       -14.841,     0.0,
                  -0.006825,   0.0,      -0.11954;

    dC_du <<  0.009614,   0.000262,  0.000183,
               0.0,      -0.029862,  0.0,
              -0.000401,  0.0,     -0.001080;
    dC_du *= (180.0 / M_PI);
  } else {
    // Rudderless coefficients
    etaBar(0) = (0.0102 * std::cos(2.0*delta_tail) - 0.0148) * beta;
    etaBar(1) = (-0.6548 * std::cos(2.0*delta_tail) - 1.3045) * alpha - 0.0787;
    etaBar(2) = (-0.15115 * std::cos(2.0*delta_tail) + 0.15115) * beta - 0.0987 * alpha * std::sin(2.0*delta_tail);

    // dC_domega (rudderless)
    dC_domega(0,0) = -0.481;
    dC_domega(0,1) = 0.065 * std::sin(2.0*delta_tail);
    dC_domega(0,2) = 0.052 + 0.010 * std::sin(2.0*delta_tail - M_PI/2.0);

    dC_domega(1,0) = 0.124 * std::sin(2.0*delta_tail);
    dC_domega(1,1) = -8.251 + 6.728 * std::sin(2.0*delta_tail - M_PI/2.0);
    dC_domega(1,2) = -1.048 * std::sin(2.0*delta_tail);

    dC_domega(2,0) = 0.004 + 0.015 * std::sin(2.0*delta_tail - M_PI/2.0);
    dC_domega(2,1) = -0.968 * std::sin(2.0*delta_tail);
    dC_domega(2,2) = -0.150 - 0.149 * std::sin(2.0*delta_tail - M_PI/2.0);

    // dC/du rudderless
    dC_du(0,0) =  0.009614;
    dC_du(0,1) =  0.000262 * std::sin(delta_tail);
    dC_du(0,2) =  0.0;

    dC_du(1,0) =  0.0;
    dC_du(1,1) = -0.0299 * std::cos(delta_tail);
    dC_du(1,2) =  0.0;

    dC_du(2,0) = -0.000401;
    dC_du(2,1) = -0.00427 * std::sin(delta_tail);
    dC_du(2,2) =  0.0;

    dC_du *= (180.0 / M_PI);
  }

  // Dynamic term
  Eigen::Vector3d omega_vec(p, q, r);
  Eigen::Vector3d dynTerm = Eigen::Vector3d::Zero();
  if (V > 1e-8) dynTerm = dC_domega * (omega_vec * (params.cref / V));

  // Control term
  Eigen::Vector3d ctrlTerm;
  ctrlTerm(0) = dC_du(0,0) * delta_a + dC_du(0,1) * delta_e + dC_du(0,2) * delta3;
  ctrlTerm(1) = dC_du(1,0) * delta_a + dC_du(1,1) * delta_e + dC_du(1,2) * delta3;
  ctrlTerm(2) = dC_du(2,0) * delta_a + dC_du(2,1) * delta_e + dC_du(2,2) * delta3;

  Eigen::Vector3d Cvec = etaBar + dynTerm + ctrlTerm; // [Cl, Cm, Cn] in stability axes

  double Cl = Cvec(0);
  double Cm = Cvec(1);
  double Cn = Cvec(2);

  // Moments in STABILITY axes
  Eigen::Vector3d M_stab;
  M_stab(0) = Cl * qS * params.bref;
  M_stab(1) = Cm * qS * params.cref;
  M_stab(2) = Cn * qS * params.bref;

  // Rotate moments from stability -> body
  Eigen::Vector3d M_b = R_s2b * M_stab;

  M_body = M_b;
}
