#include "Dynamics.h"
#include <cassert>

Dynamics::Dynamics(const Aircraft &ac_) : ac(ac_) {}

Eigen::VectorXd Dynamics::f(const Eigen::VectorXd &x, const Eigen::VectorXd &u) const {
  assert(x.size() == 12 && u.size() == 4);
  Eigen::VectorXd xdot = Eigen::VectorXd::Zero(12);

  Eigen::Vector3d F_body, M_body;
  ac.computeForcesMoments(x, u, F_body, M_body);

  double m = ac.params.m;
  Eigen::Vector3d vel(x(0), x(1), x(2));
  Eigen::Vector3d omega(x(3), x(4), x(5));

  Eigen::Vector3d acc = F_body / m - omega.cross(vel);
  xdot(0)=acc(0); xdot(1)=acc(1); xdot(2)=acc(2);

  // Rotational dynamics
  Eigen::Matrix3d I = ac.params.inertiaMat;
  Eigen::Vector3d Iomega = I * omega;
  Eigen::Vector3d omega_dot = I.inverse() * (M_body - omega.cross(Iomega));
  xdot(3)=omega_dot(0); xdot(4)=omega_dot(1); xdot(5)=omega_dot(2);

  // Euler kinematics
  double phi = x(6), theta = x(7), psi = x(8);
  Eigen::Matrix3d T;
  T << 1, std::sin(phi)*std::tan(theta), std::cos(phi)*std::tan(theta),
       0, std::cos(phi), -std::sin(phi),
       0, std::sin(phi)/std::cos(theta), std::cos(phi)/std::cos(theta);
  Eigen::Vector3d euler_dot = T * omega;
  xdot(6)=euler_dot(0); xdot(7)=euler_dot(1); xdot(8)=euler_dot(2);

  // Navigation: body -> inertial rotation
  double cphi = std::cos(phi), sphi = std::sin(phi);
  double cth  = std::cos(theta), sth = std::sin(theta);
  double cpsi = std::cos(psi), spsi = std::sin(psi);
  Eigen::Matrix3d R;
  R << cth*cpsi, sphi*sth*cpsi - cphi*spsi, cphi*sth*cpsi + sphi*spsi,
       cth*spsi, sphi*sth*spsi + cphi*cpsi, cphi*sth*spsi - sphi*cpsi,
       -sth,     sphi*cth,                  cphi*cth;
  Eigen::Vector3d vel_inertial = R * vel;
  xdot(9)=vel_inertial(0); xdot(10)=vel_inertial(1); xdot(11)=vel_inertial(2);

  return xdot;
}

// 4th-order Runge-Kutta
Eigen::VectorXd Dynamics::rk4Step(const Eigen::VectorXd &x, const Eigen::VectorXd &u, double dt) const {
  Eigen::VectorXd k1 = f(x,u);
  Eigen::VectorXd k2 = f(x + 0.5*dt*k1, u);
  Eigen::VectorXd k3 = f(x + 0.5*dt*k2, u);
  Eigen::VectorXd k4 = f(x + dt*k3, u);
  return x + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
}

// Compute Jacobians
void Dynamics::linearizeAt(const Eigen::VectorXd &x0, const Eigen::VectorXd &u0,
Eigen::MatrixXd &Aout, Eigen::MatrixXd &Bout, double eps) const {
  int n = x0.size();
  int m = u0.size();

  Aout = Eigen::MatrixXd::Zero(n, n);
  Bout = Eigen::MatrixXd::Zero(n, m);

  // Central difference for A
  for (int i = 0; i < n; ++i){
    Eigen::VectorXd xp = x0;
    Eigen::VectorXd xm = x0;
    xp(i) += eps;
    xm(i) -= eps;
    Eigen::VectorXd fp = f(xp, u0);
    Eigen::VectorXd fm = f(xm, u0);
    Aout.col(i) = (fp - fm) / (2.0 * eps);
  }

  // Central difference for B
  for (int j = 0; j < m; ++j){
    Eigen::VectorXd up = u0;
    Eigen::VectorXd um = u0;
    up(j) += eps;
    um(j) -= eps;
    Eigen::VectorXd fp = f(x0, up);
    Eigen::VectorXd fm = f(x0, um);
    Bout.col(j) = (fp - fm) / (2.0 * eps);
  }
}

