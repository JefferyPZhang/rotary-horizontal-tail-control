#pragma once
#include <Eigen/Dense>
#include "Aircraft.h"

class Dynamics {
public:
  Dynamics(const Aircraft &ac);
  Eigen::VectorXd f(const Eigen::VectorXd &x, const Eigen::VectorXd &u) const;
  void linearizeAt(const Eigen::VectorXd &x0, const Eigen::VectorXd &u0,
                   Eigen::MatrixXd &Aout, Eigen::MatrixXd &Bout, double eps = 1e-6) const;
  Eigen::VectorXd rk4Step(const Eigen::VectorXd &x, const Eigen::VectorXd &u, double dt) const;
private:
  const Aircraft &ac;
};
