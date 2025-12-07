#include "Utils.h"
#include <fstream>
#include <iomanip>
#include <cmath>

namespace Utils {
  double clamp(double x, double a, double b) {
    if (x < a) return a;
    if (x > b) return b;
    return x;
  }

  // ChatGPT-generated .csv writer
  void writeCSV(const std::string &filename, const Eigen::MatrixXd &M, const std::vector<std::string> &colnames) {
    std::ofstream f(filename);
    f<<std::fixed<<std::setprecision(10);
    if (!colnames.empty()) {
      for (size_t i=0;i<colnames.size();++i){
        f<<colnames[i]<<(i+1<colnames.size() ? "," : "\n");
      }
    }
    for (int r=0;r<M.rows();++r){
      for (int c=0;c<M.cols();++c){
        f<<M(r,c)<<(c+1<M.cols()? ",":"");
      }
      f<<"\n";
    }
    f.close();
  }
}
