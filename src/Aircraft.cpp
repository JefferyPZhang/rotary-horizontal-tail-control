#define _USE_MATH_DEFINES
#include "Aircraft.h"
#include <cmath>
#include <algorithm>

Aircraft::Aircraft(const AircraftParams &p) : params(p) {
  if (params.inertiaMat.isZero(0)) {
    Eigen::Matrix3d I; // Van Dyke, Eq 3.1
    I << 197.2026, 0.0, 146.0218,
         0.0, 1808.7634, 0.0,
         146.0218, 0.0, 1611.5609;
    const_cast<AircraftParams&>(params).inertiaMat = I;
  }
}

void Aircraft::aeroCoeffs(double alpha, double beta,
                         double delta_a, double delta_e, double delta_r, double delta_tail,
                         double &CL, double &CD, double &CY) const {
  if (params.ruddered) {
    // Ruddered (Singhasenee paper Eq 51-53)
    CL = -1.979 * alpha*alpha*alpha - 0.1339 * alpha*alpha + 5.1882 * alpha + 0.1514;
    CD = 1.1827 * alpha*alpha + 0.0582 * alpha + 0.0089;
    CY = -0.17299 * beta + 0.1293 * delta_r;
  } else {
    // Rudderless (Singhasenee paper Eq 54-56)
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
  const double V = std::sqrt(u_b*u_b + v_b*v_b + w_b*w_b);
  double alpha = 0.0, beta = 0.0;
  if (V > 1e-6) {
    alpha = std::atan2(w_b, u_b);
    beta  = std::asin(std::clamp(v_b / V, -1.0, 1.0));
  }
  double delta_a = u(0), delta_e = u(1), delta3 = u(2), delta_t = u(3);
  double delta_r = params.ruddered ? delta3 : 0.0;
  double delta_tail = params.ruddered ? 0.0 : delta3;
  double CL, CD, CY;
  aeroCoeffs(alpha, beta, delta_a, delta_e, delta_r, delta_tail, CL, CD, CY);
  const double qbar = 0.5 * params.rho * V * V;

  // Forces via lift/drag projection
  const double Lift = CL * qbar * params.S;
  const double Drag = CD * qbar * params.S;
  const double Side = CY * qbar * params.S;
  const double X_aero = -Drag * std::cos(alpha) + Lift * std::sin(alpha);
  const double Z_aero = -Drag * std::sin(alpha) - Lift * std::cos(alpha);
  const double Y_aero = Side;
  const double Fthrust = 0.25 * params.m * params.g * delta_t;

  // Total Body Forces
  F_body = Eigen::Vector3d(X_aero + Fthrust, Y_aero, Z_aero);

  // Moments
  Eigen::Vector3d etaBar;
  Eigen::Matrix3d dC_domega;
  Eigen::Matrix3d dC_du;

  if (params.ruddered) {
    // Ruddered static terms (Singhasenee paper Eq 3.4-3.6)
    etaBar(0) = -0.1422 * (beta*beta*beta) - 0.0112 * beta;
    etaBar(1) =  1.3246 * alpha*alpha*alpha - 1.9218 * alpha*alpha - 2.069 * alpha - 0.0461;
    etaBar(2) =  0.1193 * beta;

    // dC/domega for ruddered (Singhasenee paper Eq 46)
    dC_domega << -0.54597, 0.0,      0.031232,
                   0.0,     -14.841,  0.0,
                  -0.006825, 0.0,   -0.11954;

    // dC/du for ruddered (Singhasenee paper Eq 46)
    dC_du <<  0.009614,      0.000262,  0.000183,
               0.0,          -0.029862,  0.0,
              -0.000401,     0.0,      -0.001080;

    dC_du *= (180.0 / M_PI); // to deg
  } else {
    // Rudderless static eta (Singhasenee paper Eq 48)
    etaBar(0) = (0.0102 * std::cos(2.0*delta_tail) - 0.0148) * beta;
    etaBar(1) = (-0.6548 * std::cos(2.0*delta_tail) - 1.3045) * alpha - 0.0787;
    etaBar(2) = (-0.15115 * std::cos(2.0*delta_tail) + 0.15115) * beta - 0.0987 * alpha * std::sin(2.0*delta_tail);

    // dC/domega (rudderless) (Singhasenee paper Eq 49)
    dC_domega(0,0) = -0.481;
    dC_domega(0,1) = 0.065 * std::sin(2.0*delta_tail);
    dC_domega(0,2) = 0.052 + 0.010 * std::sin(2.0*delta_tail - M_PI/2.0);
    dC_domega(1,0) = 0.124 * std::sin(2.0*delta_tail);
    dC_domega(1,1) = -8.251 + 6.728 * std::sin(2.0*delta_tail - M_PI/2.0);
    dC_domega(1,2) = -1.048 * std::sin(2.0*delta_tail);
    dC_domega(2,0) = 0.004 + 0.015 * std::sin(2.0*delta_tail - M_PI/2.0);
    dC_domega(2,1) = -0.968 * std::sin(2.0*delta_tail);
    dC_domega(2,2) = -0.150 - 0.149 * std::sin(2.0*delta_tail - M_PI/2.0);

    // dC/du for rudderless (Singhasenee paper Eq 50)
    dC_du(0,0) =  0.009614;
    dC_du(0,1) =  0.000262 * std::sin(delta_tail);
    dC_du(0,2) =  0.0;
    dC_du(1,0) =  0.0;
    dC_du(1,1) = -0.0299 * std::cos(delta_tail);
    dC_du(1,2) =  0.0;
    dC_du(2,0) = -0.000401;
    dC_du(2,1) = -0.00427 * std::sin(delta_tail);
    dC_du(2,2) =  0.0;

    dC_du *= (180.0 / M_PI);  // to deg
  }

  // Dynamic term scaling
  Eigen::Vector3d omega(p, q, r);
  Eigen::Vector3d dynTerm = dC_domega * (omega * (params.cref / std::max(V, 1e-6)));

  // Control term
  Eigen::Vector3d ctrlTerm;

  // Ordered control vector
  ctrlTerm(0) = dC_du(0,0) * delta_a + dC_du(0,1) * delta_e + dC_du(0,2) * delta3;
  ctrlTerm(1) = dC_du(1,0) * delta_a + dC_du(1,1) * delta_e + dC_du(1,2) * delta3;
  ctrlTerm(2) = dC_du(2,0) * delta_a + dC_du(2,1) * delta_e + dC_du(2,2) * delta3;

  Eigen::Vector3d Cvec = etaBar + dynTerm + ctrlTerm;

  double Cl = Cvec(0), Cm = Cvec(1), Cn = Cvec(2);
  double qbar_local = qbar;
  const double Mx = Cl * qbar_local * params.S * params.bref;
  const double My = Cm * qbar_local * params.S * params.cref;
  const double Mz = Cn * qbar_local * params.S * params.bref;

  M_body = Eigen::Vector3d(Mx, My, Mz);
}
