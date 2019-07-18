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
   cout << "Begin prediction" << endl;
   x_ = F_ * x_;
   MatrixXd Ft = F_.transpose();
   P_ = F_ * P_ * Ft + Q_;
   cout << "End prediction" << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
   cout << "Begin kalman update" << endl;
   VectorXd y = z - H_ * x_;

   MatrixXd Ht = H_.transpose();
   MatrixXd S = H_ * P_ * Ht + R_;

   MatrixXd Si = S.inverse();
   MatrixXd K =  P_ * Ht * Si;

   MatrixXd I_ = MatrixXd::Identity(x_.size(), x_.size());

   // new state
   x_ = x_ + (K * y);
   P_ = (I_ - (K * H_)) * P_;
   cout << "End kalman update" << endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
   //h(x) function:
   cout << "begin ekf update" << endl;
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
   cout << "begin jacobian" << endl;
   MatrixXd Hj = tools.CalculateJacobian(x_);
   MatrixXd Hjt = Hj.transpose();
   cout << "begin S calc" << endl;
   MatrixXd S1 = Hj * P_;
   cout << "check 1" << endl;
   MatrixXd S2 = S1 * Hjt;
   cout << "check 2" << endl;
   cout << R_ << endl;
   MatrixXd S = S2 + R_;
   cout << "check 3" << endl;


   MatrixXd Si = S.inverse();
   cout << "begin K calc" << endl;
   MatrixXd K =  P_ * Hjt * Si;
   cout << "begin I calc" << endl;
   MatrixXd I_ = MatrixXd::Identity(x_.size(), x_.size());
   // new state
   cout << "begin state update" << endl;
   x_ = x_ + (K * y);
   P_ = (I_ - (K * Hj)) * P_;
   cout << "end ekf update" << endl;
}
