#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace Utils {
  using Matrix = Eigen::MatrixXd;
  using Vector = Eigen::VectorXd;

  double clamp(double x, double a, double b);
  void writeCSV(const std::string &filename, const Eigen::MatrixXd &M, const std::vector<std::string> &colnames = {});
}
