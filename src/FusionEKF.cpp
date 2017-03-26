#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


/*  EKF

    Prediction step:
    ===============
 
    x′ = F∗x + noise
 
    where:
      x       = state of the object being tracked
      F       = state transition matrix (position of x after time Δt); F is a matrix that, when multiplied with x, predicts where the object will be after time Δt.
      x′      = predicted state of x
      noise   = process uncertainty
 
    Update step: 
    ============
 
    For the update step, we use the measurement function to map the state vector into the measurement space of the sensor. To give a concrete example, lidar only 
    measures an object's position. But the extended Kalman filter models an object's position and velocity. So multiplying by the measurement function H matrix will 
    drop the velocity information from the state vector x. Then the lidar measurement position and our belief about the object's position can be compared.

    z = H∗x + w
    
    where:
      z       = z is the measurement vector. For a lidar sensor, the z vector contains the position−x and position−y measurements.
      w       = measurement uncertainty (sensor measurement noise).
      H       = H measurement function matrix; maps state back to measurement space (ex: cartesian or polar coordinates. polar for radar, cartesian for lidar).
                H is the matrix that projects your belief about the object's current state into the measurement space of the sensor. For lidar, this is a fancy way
                of saying that we discard velocity information from the state variable since the lidar sensor only measures position: The state vector x contains
                information about [p​x, p​y, v​x, v​y] whereas the z vector will only contain [px,py].
 

    State transition matdix F:
    ==========================
    F(4,4) = (1,  0,  Δt, 0 )
             (0,  1,  0,  Δt)
             (0,  0,  1,  0 )
             (0,  0,  0,  1 )
 
    Process covariance matrix Q:
    ===========================
 
    P​′ ​= FP(F​t) + Q
 
​    Q(4,4)  = (pow(Δt,4)/4​​ * pow(σ(​ax), 2),     0,                              pow(Δt,3)/2​​ * pow(σ(​ax), 2),       0                           )
              (0,                               pow(Δt,4)/4​​ * pow(σ(​ay), 2),    0,                                 pow(Δt,3)/2​​ * pow(σ(​ay), 2) )
​              (pow(Δt,3)/2​​ * pow(σ(​ax), 2),     0,                              pow(Δt,2)​​ * pow(σ(​ax), 2),         0                           )
              (0,                               pow(Δt,3)/2​​ * pow(σ(​ay), 2),    0,                                 pow(Δt,2)​​ * pow(σ(​ay), 2)   )
 ​ 
 
    where:
      pow(σ(​ax), 2)    = sigma (ax) squared; acceleration covariance in the x direction squared
      pow(σ(​ay), 2)    = sigma (ay) squared; acceleration covariance in the y direction squared
 
    
    Measurement x (Lidar):
    ======================
        x = (px)
            (py)

    Measurement x (Radar):
    ======================
        x = (​ρ )
            (ϕ )
            (ρ˙)
 
        where:
              ​ρ      = rho; radian distance from origin
              ϕ      = phi; bearing
              ρ˙     = rho dot; Radial velocity
​​  
 
    Measurement function matrix H (Lidar):
    ======================================
      2d observation space (px) <- 4d state (px)
                           (py)             (py)
                                            (vx)
                                            (vy)
 
      H(lidar) = (1, 0, 0, 0)
                 (0, 1, 0, 0)

 
    Measurement function matrix H (Radar):
    ======================================
 ​     x = (​ρ ) <- 4d state (px)
          (ϕ )             (py)
          (ρ˙)             (vx)
                           (vy)
 
 
      H(radar) = ( √(pow(px, 2) + pow(py, 2))                           )
                 ( arctan(px / py)                                      )
                 ( ((px * vx) + (py + vy)) / √(pow(px, 2) + pow(py, 2)) )
 
 
      Linear approximation of H(radar) using multi-dimensional Taylor series expansion:
 
      H​j =  ( px / √(pow(px,2) + pow(py, 2)),                                     py / (√(pow(px,2) + pow(py, 2)),                                      0,                                   0                            )
            ( -(py / (pow(px,2) + pow(py, 2)),                                    px / (pow(px,2) + pow(py, 2)),                                        0,                                   0                            )
            ( py * (vx * py - vy * px) / pow((pow(px,2) + pow(py, 2)), 3/2),      px * (vy * px - vx * py) / pow((pow(px,2) + pow(py, 2)), 3/2),        px / √(pow(px,2) + pow(py, 2)),      py / √(pow(px,2) + pow(py, 2))
 ​​ 

    Measurement Covariance (noise) Matrix R (Lidar):
    ===============================================
    
    Represents the uncertainty in the position measurements we receive from the laser sensor; The dimensions of the R matrix is square and each side of its matrix 
    is the same length as the number of measurements parameters.
 
    R = E[ωωT] = (​pow(σ(​px), 2),   0            )
                 (0,               ​pow(σ(​py), 2))
 
 
    where:
        ω       = measurement noise
        σ(​px)   = measurement covariance in the x direction
        σ(​py)   = measurement covariance in the y direction
 
    Measurement Covariance (noise) Matrix R (Radar):
    ================================================
 
    R = (pow(σ(​ρ)​, 2),     0,                 0              )
        (0,                pow(σ(ϕ)​, 2),      0              )
        (0,                0,                 pow(σ(ρ˙)​, 2)  )


    Additional notes:
 
    1. Generally, the parameters for the random noise measurement matrix will be provided by the sensor manufacturer
    2. for measurement updates with lidar, we can use the H matrix for calculating y, S, K and P.
       for radar, H​j is used to calculate S, K and P.

*/

FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // measurement covariance matrix - laser (sensor measurement uncertainty)
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement covariance matrix - radar (sensor measurement uncertainty)
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // measurement function (laser)
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  ekf_.H_ = H_laser_;
  
  // Object state
  ekf_.x_ = VectorXd(4);
  
  // state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
  
  // the initial state transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
  
  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << 0,0,0,0,
             0,0,0,0,
             0,0,0,0,
             0,0,0,0;
  
  ekf_.u_ = VectorXd(4);
  ekf_.u_ << 0,0,0,0;
}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    //ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      VectorXd radar  = measurement_pack.raw_measurements_;
      
      // where:
      // radar[0]      = rho (ρ); radian distance from origin
      // radar[1]      = phi (ϕ); bearing
      // radar[2]      = rho_dot (ρ˙); Radial velocity
      
      double rho      = radar[0];
      double phi      = radar[1];
      double rho_dot  = radar[2];
      
      double x        = cos(phi) * rho;       // cos(phi) = adj/hyp = x / rho
      double y        = sin(phi) * rho;       // sin(phi) = opp/hyp = y / rho
      double vx       = cos(phi) * rho_dot;   // cos(phi) = adj/hyp = vx / rho_dot
      double vy       = sin(phi) * rho_dot;   // sin(phi) = opp/hyp = vy / rho_dot
      
      ekf_.x_ << x, y, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      VectorXd lidar = measurement_pack.raw_measurements_;
      // where:
      // lidar[0] = x position
      // lidar[1] = y position
      
      double x  = lidar[0];
      double y  = lidar[1];
      double vx = 0;        // lidar measurements do not contain velocity
      double vy = 0;        // lidar measurements do not contain velocity
      ekf_.x_ << x, y, vx, vy;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

 
  // Compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // Update the state transition matrix F according to the new elapsed time.
  // Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  // Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
  float noise_ax = 9;
  float noise_ay = 9;
  
  // Update the process noise covariance matrix Q.
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  ekf_.Q_ <<  dt_4/4*noise_ax,    0,                dt_3/2*noise_ax,  0,
              0,                  dt_4/4*noise_ay,  0,                dt_3/2*noise_ay,
              dt_3/2*noise_ax,    0,                dt_2*noise_ax,    0,
              0,                  dt_3/2*noise_ay,  0,                dt_2*noise_ay;
  
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  ekf_.UpdateEKF(measurement_pack.sensor_type_, measurement_pack.raw_measurements_, R_laser_, R_radar_);

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
