#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  TODO:
    * predict the state
  */
  x_  = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
  MatrixXd s = H_ * P_ * H_.transpose() + R_;
  MatrixXd k = P_ * H_.transpose() * s.inverse();
  x_ = x_ + k * y;

  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  P_ = (I - k * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  // convert x object state to rho, theta, rho_dot coordinates
  double rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  double theta = atan(x_(1) / x_(0));
  double rho_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;

  //calculate h(x_) instead of h * x_
  VectorXd h = VectorXd(3);
  h << rho, theta, rho_dot;

  VectorXd y = z - h;
  MatrixXd s = H_ * P_ * H_.transpose() + R_;
  MatrixXd k = P_ * H_.transpose() * s.inverse();
  x_ = x_ + k * y;

  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  P_ = (I - k * H_) * P_;

}
