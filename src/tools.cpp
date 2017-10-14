#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd error(estimations[0].size());
  error.fill(0.);
  int n = estimations.size();
  VectorXd diff(error.size());
  for (int i = 0; i < n; ++i) {
    diff = estimations[i] - ground_truth[i];
    diff = diff.array()*diff.array();
    error += diff;
  }
  error /= n;
  error = error.array().sqrt();
  return error;
}
