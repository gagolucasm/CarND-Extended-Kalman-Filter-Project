#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

 x_ = F_*x_;
 MatrixXd Ft = F_.transpose();
 P_ = F_*P_*Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {


  MatrixXd PHt = P_ * H_.transpose();
  MatrixXd S = H_ * PHt + R_;
  MatrixXd K = PHt * S.inverse();
  

  x_ = x_ + K *(z - H_ * x_);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  double range = sqrt( pow(x_[0],2) + pow(x_[1],2) );
  double bearing;
  double range_rate;
  if (fabs(range > 0.001)) {
    bearing = atan(x_[1] / x_[0]);
    range_rate = ((x_[0] * x_[2] + x_[1] * x_[3]) / range);
  } else {
    // Too low, use 0  
    bearing = 0;
    range_rate = 0;
  }


  MatrixXd z_pred(3, 1);
  z_pred << range, bearing, range_rate;  

  MatrixXd PHt = P_ * H_.transpose();
  MatrixXd S = H_ * PHt + R_;
  MatrixXd K = PHt * S.inverse();
  
  //new estimate
  VectorXd y = z - z_pred;
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
