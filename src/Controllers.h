#pragma once
#include <Eigen/Dense>

namespace Controllers {
  Eigen::MatrixXd solveDiscreteRiccati(const Eigen::MatrixXd &Ad, const Eigen::MatrixXd &Bd,
                                       const Eigen::MatrixXd &Q, const Eigen::MatrixXd &R,
                                       int maxIters=1000, double tol=1e-9);
  Eigen::MatrixXd computeLQRGainFromContinuous(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                                               const Eigen::MatrixXd &Q, const Eigen::MatrixXd &R,
                                               double dt=0.01);
}
