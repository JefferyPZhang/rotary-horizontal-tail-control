#include "Dynamics.h"

Dynamics::Dynamics(const Aircraft &ac_) : ac(ac_) {}

// Nonlinear dynamics
Eigen::VectorXd Dynamics::f(const Eigen::VectorXd &x, const Eigen::VectorXd &u) const {
  Eigen::VectorXd xdot = Eigen::VectorXd::Zero(12);

  // Forces & moments in body frame
  Eigen::Vector3d F_body, M_body;
  ac.computeForcesMoments(x, u, F_body, M_body);

  // Mass & inertia
  double m = ac.params.m;
  Eigen::Matrix3d I = ac.params.inertiaMat;

  // Body velocities and angular rates
  Eigen::Vector3d vel(x(0), x(1), x(2));
  Eigen::Vector3d omega(x(3), x(4), x(5));

  // Translational acceleration (body frame):
  Eigen::Vector3d acc = F_body / m - omega.cross(vel);
  xdot(0) = acc(0);
  xdot(1) = acc(1);
  xdot(2) = acc(2);

  // Rotational dynamics
  Eigen::Vector3d omega_dot = I.inverse() * (M_body - omega.cross(I * omega));
  xdot(3) = omega_dot(0);
  xdot(4) = omega_dot(1);
  xdot(5) = omega_dot(2);

  // Euler angle kinematics
  double phi = x(6), theta = x(7), psi = x(8);
  double cphi = std::cos(phi), sphi = std::sin(phi);
  double cth  = std::cos(theta), sth = std::sin(theta);

  Eigen::Matrix3d T;
  T << 1.0, sphi * std::tan(theta), cphi * std::tan(theta),
       0.0, cphi,                  -sphi,
       0.0, sphi / cth,            cphi / cth;

  Eigen::Vector3d euler_dot = T * omega;
  xdot(6) = euler_dot(0);
  xdot(7) = euler_dot(1);
  xdot(8) = euler_dot(2);

  // Navigation: body -> inertial
  double cpsi = std::cos(psi), spsi = std::sin(psi);
  Eigen::Matrix3d R_1v;
  Eigen::Matrix3d R_21;
  Eigen::Matrix3d R_b2;
  R_1v  << cpsi,  spsi,  0,
           -spsi, cpsi,  0,
           0,     0,     1;
  R_21  << cth,   0,     -sth,
           0,     1,     0,
           sth,   0,     cth;
  R_b2  << 1,     0,     0,
           0,     cphi,  sphi,
           0,     -sphi, cphi;
  Eigen::Matrix3d R_b2i = R_b2 * R_21 * R_1v;
  Eigen::Vector3d vel_inertial = R_b2i * vel;
  xdot(9)  = vel_inertial(0);
  xdot(10) = vel_inertial(1);
  xdot(11) = vel_inertial(2);

  return xdot;
}

// RK4 integrator
Eigen::VectorXd Dynamics::rk4Step(const Eigen::VectorXd &x, const Eigen::VectorXd &u, double dt) const {
  Eigen::VectorXd k1 = f(x, u);
  Eigen::VectorXd k2 = f(x + 0.5*dt*k1, u);
  Eigen::VectorXd k3 = f(x + 0.5*dt*k2, u);
  Eigen::VectorXd k4 = f(x + dt*k3, u);
  return x + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
}

// Central-difference Jacobians
void Dynamics::linearizeAt(const Eigen::VectorXd &x0, const Eigen::VectorXd &u0,
                           Eigen::MatrixXd &Aout, Eigen::MatrixXd &Bout, double eps) const {
  int n = x0.size();
  int m = u0.size();
  Aout = Eigen::MatrixXd::Zero(n, n);
  Bout = Eigen::MatrixXd::Zero(n, m);

  // central difference for A
  for (int i = 0; i < n; ++i) {
    Eigen::VectorXd xp = x0; xp(i) += eps;
    Eigen::VectorXd xm = x0; xm(i) -= eps;
    Eigen::VectorXd fp = f(xp, u0);
    Eigen::VectorXd fm = f(xm, u0);
    Aout.col(i) = (fp - fm) / (2.0 * eps);
  }

  // central difference for B
  for (int j = 0; j < m; ++j) {
    Eigen::VectorXd up = u0; up(j) += eps;
    Eigen::VectorXd um = u0; um(j) -= eps;
    Eigen::VectorXd fp = f(x0, up);
    Eigen::VectorXd fm = f(x0, um);
    Bout.col(j) = (fp - fm) / (2.0 * eps);
  }
}
