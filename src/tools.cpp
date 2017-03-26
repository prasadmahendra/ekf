#include <iostream>
#include <math.h>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth) {
  // Calculate the RMSE
  
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if ( estimations.size() == 0 || estimations.size() != ground_truth.size() ) {
    return rmse;
  }
  
  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    VectorXd est = estimations[i];
    VectorXd actual = ground_truth[i];
    VectorXd delta = est - actual;
    VectorXd delta_sq = delta.array() * delta.array();
    rmse = rmse + delta_sq;
  }
  
  //calculate the mean
  rmse = rmse / (double)estimations.size();
  
  //calculate the squared root
  rmse = rmse.array().sqrt();
  
  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  // Calculate a Jacobian.
  MatrixXd Hj(3,4);
  
  float px2 = pow(px, 2);
  float py2 = pow(py, 2);
  float px2_py2 = px2 + py2;
  float px2_py2_sqrt = sqrt(px2_py2);
  
  //check division by zero
  if (px2_py2 == 0) {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    Hj << 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0;
    return Hj;
  }
  
  //compute the Jacobian matrix
  Hj << px / px2_py2_sqrt, py / px2_py2_sqrt, 0, 0,
        (-1 * py) / px2_py2, px / px2_py2, 0, 0,
        (py * (vx * py - vy * px)) / pow(px2_py2, 3/2), (px * (vy * px - vx * py)) / pow(px2_py2, 3/2), px / px2_py2_sqrt, py / px2_py2_sqrt;
  
  return Hj;
}
