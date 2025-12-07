#include "Controllers.h"
#include <Eigen/Dense>

namespace Controllers {

Eigen::MatrixXd solveDiscreteRiccati(const Eigen::MatrixXd &Ad, const Eigen::MatrixXd &Bd,
                               const Eigen::MatrixXd &Q, const Eigen::MatrixXd &R,
                               int maxIters, double tol) {
  int n = Ad.rows();
  Eigen::MatrixXd P = Q;
  for (int k=0;k<maxIters;++k){
    Eigen::MatrixXd BtP = Bd.transpose() * P;
    Eigen::MatrixXd S = R + BtP * Bd;
    Eigen::MatrixXd Sinv = S.inverse();
    Eigen::MatrixXd Ktemp = Sinv * BtP * Ad;
    Eigen::MatrixXd Pnext = Ad.transpose() * P * Ad - Ad.transpose() * P * Bd * Ktemp + Q;
    if ((Pnext - P).norm() < tol) { P = Pnext; break; }
    P = Pnext;
  }
  return P;
}

Eigen::MatrixXd computeLQRGainFromContinuous(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
                                             const Eigen::MatrixXd &Q, const Eigen::MatrixXd &R,
                                             double dt) {
  // Bilinear Discretization (Tustin Transform)
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(A.rows(), A.cols());
  Eigen::MatrixXd Ad = (I - 0.5*dt*A).inverse() * (I + 0.5*dt*A);
  Eigen::MatrixXd Bd = (I - 0.5*dt*A).inverse() * (dt * B);
  Eigen::MatrixXd P = solveDiscreteRiccati(Ad,Bd,Q,R);
  Eigen::MatrixXd S = R + Bd.transpose() * P * Bd;
  Eigen::MatrixXd K = S.inverse() * (Bd.transpose() * P * Ad);
  return K;
}

} // namespace
