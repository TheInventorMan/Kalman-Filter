#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
   x_ = F_ * x_;
   MatrixXd Ft = F_.transpose();
   P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
   VectorXd y = z - H_ * x_;

   MatrixXd Ht = H_.transpose();
   MatrixXd S = H_ * P_ * Ht + R_;

   MatrixXd Si = S.inverse();
   MatrixXd K =  P_ * Ht * Si;

   MatrixXd I_ = MatrixXd::Identity(x_.size(), x_.size());

   // new state
   x_ = x_ + (K * y);
   P_ = (I_ - (K * H_)) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
   //h(x) function:
   float px = x_(0);
   float py = x_(1);
   float vx = x_(2);
   float vy = x_(3);

   float rho = sqrt(px*px + py*py);
   float phi = atan2(py,px);

   if (fabs(rho)<0.00001){
     rho = 0.00001;
   }

   float rho_dot = (px*vx + py*vy) / rho;

   VectorXd hx(3);
   hx << rho, phi, rho_dot;

   VectorXd y = z - hx;

   for (int i=0; y(1) > M_PI;  y(1) -= 2*M_PI) {

   }
   for (int i=0; y(1) < -M_PI; y(1) += 2*M_PI) {

   }
   MatrixXd Hj = tools.CalculateJacobian(x_);
   MatrixXd Hjt = Hj.transpose();
   MatrixXd S =  Hj * P_ * Hjt + R_;
   MatrixXd Si = S.inverse();
   MatrixXd K =  P_ * Hjt * Si;

   MatrixXd I_ = MatrixXd::Identity(x_.size(), x_.size());
   // new state
   x_ = x_ + (K * y);
   P_ = (I_ - (K * Hj)) * P_;

}
