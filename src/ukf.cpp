#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5; // state dim
  n_aug_ = 7; // aug dim
  n_sigma_ = 1+2*n_aug_; // sigma dim
  // initial state vector
  x_ = VectorXd(n_x_);
  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .5;
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;
  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;
  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;
  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;
  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  // Sigma points spreading
  lambda_ = 3. - n_aug_;
  // Weights of sigma points
  weights_ = VectorXd(n_sigma_);
  weights_.fill(0.5/(lambda_+n_aug_));
  weights_(0) = lambda_/(lambda_+n_aug_);
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (is_initialized_ == false) {
    x_ = VectorXd(n_x_);
    P_ = MatrixXd::Identity(n_x_, n_x_);
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0];
      double theta = meas_package.raw_measurements_[1];
      x_ << rho*std::cos(theta), rho*std::sin(theta), 0., 0., 0.;
    } else {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0., 0., 0.;
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  double delta_t = (meas_package.timestamp_ - time_us_)/1e6;
  time_us_ = meas_package.timestamp_; 
  Prediction(delta_t);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  }
}

MatrixXd UKF::GenerateSigmaPoints(const VectorXd& x_aug, const MatrixXd& P_aug) {
  double coef = std::sqrt(n_aug_+lambda_);
  MatrixXd sqrtP = P_aug.llt().matrixL();
  MatrixXd sigma_pts(n_aug_, n_sigma_);
  sigma_pts.col(0) = x_aug;
  for (int i = 0; i < n_aug_; ++i) {
    sigma_pts.col(i+1) = x_aug + coef*sqrtP.col(i);
    sigma_pts.col(i+1+n_aug_) = x_aug - coef*sqrtP.col(i);
  }
  return sigma_pts;
}

void UKF::PredictSigmaPoints(double delta_t) {
  VectorXd x_aug(n_aug_);
  x_aug.fill(0.);
  x_aug.head(n_x_) = x_;

  MatrixXd P_aug(n_aug_, n_aug_);
  P_aug.fill(0.);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;
  auto sigma_pts = GenerateSigmaPoints(x_aug, P_aug);
  VectorXd delta_x(n_x_);
  delta_x.fill(0.);
  VectorXd delta_noise(n_x_);
  delta_noise.fill(0.);
  for (int i = 0; i < n_sigma_; ++i) {
    double px = sigma_pts(0, i);
    double py = sigma_pts(1, i);
    double v  = sigma_pts(2, i);
    double psi  = sigma_pts(3, i);
    double psid = sigma_pts(4, i);
    double nu_a = sigma_pts(5, i);
    double nu_psidd = sigma_pts(6, i);
    if (std::fabs(psid) < 1e-4) {
      delta_x(0) = v*std::cos(psi)*delta_t;
      delta_x(1) = v*std::sin(psi)*delta_t;
    }
    else {
      delta_x(0) = v*(std::sin(psi+delta_t*psid)-std::sin(psi))/psid;
      delta_x(1) = v*(std::cos(psi)-std::cos(psi+delta_t*psid))/psid;
      delta_x(3) = psid*delta_t;
    }
    delta_noise(0) = 0.5*delta_t*delta_t*std::cos(psi)*nu_a;
    delta_noise(1) = 0.5*delta_t*delta_t*std::sin(psi)*nu_a;
    delta_noise(2) = delta_t*nu_a;
    delta_noise(3) = 0.5*delta_t*delta_t*nu_psidd;
    delta_noise(4) = delta_t*nu_psidd;
    Xsig_pred_.col(i) = sigma_pts.col(i).head(n_x_) + delta_x + delta_noise;
    Xsig_pred_(3,i) = std::fmod(Xsig_pred_(3,i), PI2);
  } 
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  PredictSigmaPoints(delta_t);
  x_.fill(0.);
  for (int i = 0; i < n_sigma_; ++i) {
    x_ += weights_(i)*Xsig_pred_.col(i);
  }
  P_.fill(0.);
  for (int i = 0; i < n_sigma_; ++i) {
    VectorXd diff = Xsig_pred_.col(i) - x_;
    diff(3) = std::fmod(diff(3), PI2);
    P_ += weights_(i)*diff*diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& meas_package) {
  // Copy px and py directly
  MatrixXd Zsig = Xsig_pred_.topLeftCorner(2, n_sigma_);
  // Copy the mean of px and py
  VectorXd z_pred = x_.head(2);
  // Compute covariance matrix
  MatrixXd S(2, 2);
  S.fill(0.);
  VectorXd diff(2); 
  for (int i = 0; i < n_sigma_; ++i) {
    diff = Zsig.col(i) - z_pred;
    S += weights_(i)*diff*diff.transpose();
  }
  S(0, 0) += std_laspx_*std_laspx_;
  S(1, 1) += std_laspy_*std_laspy_;
  // Finish the update and calculate the normalized innovation squared
  double NIS = Update(meas_package, z_pred, S, Zsig);
  ++n_lidar_;
  if (NIS > chi95_lidar_) ++n_over95_lidar_;
  std::cout << "Percentage of LIDAR NIS observations over 95 chi^2 percentile: " << double(n_over95_lidar_)*100./n_lidar_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& meas_package) {
  MatrixXd Zsig(3, n_sigma_);
  for (int i = 0; i < n_sigma_; ++i) {
    double px  = Xsig_pred_(0, i);
    double py  = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double rho = std::sqrt(px*px + py*py);
    // if both px and py are zero, atan2 is undefined, so we skip updating
    if (rho < 1e-6) return;
    double phi = std::atan2(py, px);
    double rhod = (px*std::cos(yaw)*v + py*std::sin(yaw)*v)/rho;
    Zsig(0, i) = rho;
    Zsig(1, i) = phi;
    Zsig(2, i) = rhod;
  }
  // Compute mean
  VectorXd z_pred(3);
  z_pred.fill(0.);
  for (int i = 0; i < n_sigma_; ++i) {
    z_pred += weights_(i)*Zsig.col(i);
  }
  // Compute covariance matrix
  MatrixXd S(3, 3);
  S.fill(0.);
  VectorXd diff(3); 
  for (int i = 0; i < n_sigma_; ++i) {
    diff = Zsig.col(i) - z_pred;
    diff(1) = std::fmod(diff(1), PI2);
    S += weights_(i)*diff*diff.transpose();
  }
  S(0, 0) += std_radr_*std_radr_;
  S(1, 1) += std_radphi_*std_radphi_;
  S(2, 2) += std_radrd_*std_radrd_;
  // Finish the update and calculate the normalized innovation squared
  double NIS = Update(meas_package, z_pred, S, Zsig);
  ++n_radar_;
  if (NIS > chi95_radar_) ++n_over95_radar_;
  std::cout << "Percentage of RADAR NIS observations over 95 chi^2 percentile: " << double(n_over95_radar_)*100./n_radar_<< std::endl;
}

double UKF::Update(const MeasurementPackage& meas_package, 
    const VectorXd& z_pred, const MatrixXd& S, const MatrixXd& Zsig) {
  int n_meas = z_pred.size();
  MatrixXd T(n_x_, n_meas);
  // Compute the cross-covariance matrix
  T.fill(0.);
  VectorXd xdiff(n_x_);
  VectorXd zdiff(n_meas);
  for (int i = 0; i < n_sigma_; ++i) {
    xdiff = Xsig_pred_.col(i) - x_;
    xdiff(3) = std::fmod(xdiff(3), PI2);
    zdiff = Zsig.col(i) - z_pred;
    // If the size of the measurement is 3, then we have a RADAR measurement
    // so we normalize the angle
    if (n_meas == 3) zdiff(1) = std::fmod(zdiff(1), PI2);
    T += weights_(i)*xdiff*zdiff.transpose();
  }
  MatrixXd Sinv = S.inverse();
  MatrixXd K = T*Sinv;
  auto& z = meas_package.raw_measurements_;
  zdiff = z - z_pred;
  if (n_meas == 3) zdiff(1) = std::fmod(zdiff(1), PI2);
  x_ += K*zdiff;
  P_ -= K*S*K.transpose();
  // Return NIS
  return zdiff.transpose()*Sinv*zdiff;
}
