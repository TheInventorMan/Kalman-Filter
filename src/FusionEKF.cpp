#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // Initializing matrices
  MatrixXd R_laser_(2, 2);
  MatrixXd R_radar_(3, 3);
  MatrixXd H_laser_(2, 4);

  // Measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // Measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

   // Initialize state matrices
   VectorXd x_(4);
   MatrixXd F_(4,4);
   MatrixXd P_(4,4);
   MatrixXd Q_(4,4);

   x_ << 0.0, 0.0, 1.0, 1.0;

   F_ << 1.0, 0.0, 1.0, 0.0,
         0.0, 1.0, 0.0, 1.0,
         0.0, 0.0, 1.0, 0.0,
         0.0, 0.0, 0.0, 1.0;

   // Don't have any information about noise and covariance yet
   P_ << 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0;

   Q_ << 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0;

   ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    // First measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Extract radar data from pack
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];

      // Convert to cartesian and set as initial state
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);

      ekf_.x_ << px, py, vx, vy;

      // Fairly high certainty about both position and velocity
      ekf_.P_ << 1.0, 0.0, 0.0, 0.0,
                 0.0, 1.0, 0.0, 0.0,
                 0.0, 0.0, 1.0, 0.0,
                 0.0, 0.0, 0.0, 1.0;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Extract laser data
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];

      ekf_.x_ << px, py, 0.0, 0.0;

      // No information about velocity, so arbitrarily high variance.
      ekf_.P_ << 1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 100.0, 0.0,
                0.0, 0.0, 0.0, 100.0;
    }
    // Done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /**
   * Prediction
   */
   // Get elapsed time since last measurement
   float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
   previous_timestamp_ = measurement_pack.timestamp_;

   // Update state transition matrix with time element.
   ekf_.F_(0, 2) = dt;
   ekf_.F_(1, 3) = dt;

   // Compute noise matrix Q
   double dt2 = dt*dt;
   double dt3 = dt2*dt/2;
   double dt4 = dt2*dt2/4;

   double nx = 9.0; //noise_ax
   double ny = 9.0; //noise_ay

   ekf_.Q_ << dt4*nx, 0.0, dt3*nx, 0.0,
              0.0, dt4*ny, 0.0, dt3*ny,
              dt3*nx, 0.0, dt2*nx, 0.0,
              0.0, dt3*ny, 0.0, dt2*ny;

  // Make initial prediction
  ekf_.Predict();

  /**
   * Update
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Set radar measurement covariance matrix and update with EKF equations.
    MatrixXd R_radar_(3,3);
    R_radar_ << 0.09, 0, 0,
            0, 0.0009, 0,
            0, 0, 0.09;
    ekf_.R_ = R_radar_;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Set laser measurement covariance matrix and update with regular Kalman equations.
    MatrixXd R_laser_(2,2);
    R_laser_ << 0.0225, 0,
            0, 0.0225;
    ekf_.R_ = R_laser_;

    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // Print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
