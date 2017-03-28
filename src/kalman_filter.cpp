#include "kalman_filter.h"
#include <stdlib.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in,
                        MatrixXd &P_in,
                        VectorXd &u_in,
                        MatrixXd &F_in,
                        MatrixXd &H_in,
                        MatrixXd &R_in,
                        MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  u_ = u_in;
  F_ = F_in;
  H_ = H_in;
  Q_ = Q_in;
  debug = false;
}

void KalmanFilter::Predict() {
  // predict the state
  x_ = F_ * x_ + u_;                      // new predicted state (x) and covariance (P)
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;                 // where Q is the process covariance (process noise due to uncertainty in motion/velocity)
}

void KalmanFilter::UpdateEKF(MeasurementPackage::SensorType sensorType, const VectorXd &z, const MatrixXd &R_lidar, const MatrixXd &R_radar) {
  // update the state by using Extended Kalman Filter equations
  
  MatrixXd H;           // H here is sensor dependent (see below)
  MatrixXd R;           // R here is sensor dependent (see below)
  VectorXd z_pred;      // x_ mapped back to sensor space (eg: polar or cartesian)
  
  if (sensorType == MeasurementPackage::SensorType::RADAR ) {
    H = tools.CalculateJacobian(x_);  // Linear approximation of the non-linear measurement function matrix H
    R = R_radar;                      // Measurement covariance matrix - radar (sensor measurement uncertainty)

    /*
      z is in polar coordinates. x_ is in cartesian. convert x to polar coords first so that we
      can calculate z_pred for our error calculation (y)
     
      Recall that:
     
      Measurement x (Radar):
      ======================
     
      z_pred = (​ρ )  <-- x_ = (px)
               (ϕ )           (py)
               (ρ˙)           (vx)
                              (vy)
     
      H(radar) = ( √(pow(px, 2) + pow(py, 2))                           )
                 ( arctan(px / py)                                      )
                 ( ((px * vx) + (py + vy)) / √(pow(px, 2) + pow(py, 2)) )
     
      where:
          z_pred = Radar measurement prediction in polar coordinates
          x_     = Radar measurement prediction in cartesian coordinates
          ​ρ      = rho; radian distance from origin
          ϕ      = phi; bearing
          ρ˙     = rho dot; Radial velocity
     
     
      Additional Notes:
     
      In C++, atan2() returns values between -pi and pi. When calculating phi in y = z - h(x) for radar measurements, 
      the resulting angle phi in the y vector should be adjusted so that it is between -pi and pi. The Kalman filter is expecting small angle values between the range -pi and pi.
    */

    double px = x_(0);    // pos x
    double py = x_(1);    // pos y
    double vx = x_(2);    // velocity x
    double vy = x_(3);    // velocity y
    
    double rho = sqrt(pow(px, 2) + pow(py, 2));                 // Radian distance from origin
    double phi = px != 0 || py != 0 ? atan2(py, px) : 0;        // bearing; tan(ϕ) = opp/adj = phi; arctan(phi) = phi
    double phi_dot = rho > 0 ? (px*vx + py*vy) / rho : 0;       // Radial velocity
    
    z_pred = VectorXd(3);
    z_pred << rho, phi, phi_dot;
  } else {
    H = H_;                           // Measurement function matrix H (Lidar)
    R = R_lidar;                      // Measurement covariance matrix - lidar (sensor measurement uncertainty)
    z_pred = H * x_;                  // map predicted measurement to object state
  }
  
  VectorXd y = z - z_pred;            // error calculation (y) given new object state (z)
  MatrixXd Ht = H.transpose();        // H matrix transposed
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H * PHt + R;           // S matrix
  MatrixXd Si = S.inverse();          // Inverse of the S matrix
  MatrixXd K = PHt * Si;              // The Kalman gain
  
  // new estimate
  x_ = x_ + (K * y);                  // new updated state (x) and state covariance (P)
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;              // updated state covariance matrix

  if ( debug ) {
    DumpMatrix("H", H);
    DumpMatrix("z", z);
    DumpMatrix("z_pred", z_pred);
    DumpMatrix("x_", x_);
    DumpMatrix("P_", P_);
  }
}

void KalmanFilter::DumpMatrix(std::string name, const MatrixXd &x) {
  std::cout << "Matrix dumper: " << std::endl;
  std::cout << name << " = " << std::endl << x << std::endl;
  std::cout << std::endl;
}
