#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
//  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  //use_radar_ = false;

  n_x_ = 5; // state dim
  n_aug_ = 7; // aug dim
  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.55555;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  lambda_ = 3. - n_aug_;
  int n_sigma = 1+2*n_aug_;
  weights_ = VectorXd(n_sigma);
  weights_.fill(0.5/(lambda_+n_aug_));
  weights_(0) = lambda_/(lambda_+n_aug_);
  Xsig_pred_ = MatrixXd(n_x_, n_sigma);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (is_initialized_ == false) {
    x_ = VectorXd(n_x_);
    P_ = MatrixXd::Identity(n_x_, n_x_);
    P_(2, 2) = 1000.;
    P_(3, 3) = 1000.;
    P_(4, 4) = 1000.;
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
  std::cout << "Before update" << std::endl;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    std::cout << "Update Radar" << std::endl;
    UpdateRadar(meas_package);
  } else {
    std::cout << "Update Lidar" << std::endl;
    UpdateLidar(meas_package);
  }
}

MatrixXd UKF::GenerateSigmaPoints(const VectorXd& x_aug, const MatrixXd& P_aug) {
  double coef = std::sqrt(n_aug_+lambda_);
  MatrixXd sqrtP = P_aug.llt().matrixL();
  MatrixXd sigma_pts(n_aug_, 1 + 2*n_aug_);
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
  int n_sigma = 1 + 2*n_aug_;
  VectorXd delta_x(n_x_);
  delta_x.fill(0.);
  VectorXd delta_noise(n_x_);
  delta_noise.fill(0.);
  for (int i = 0; i < n_sigma; ++i) {
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
  } 
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  PredictSigmaPoints(delta_t);
  int n_sigma = 1 + 2*n_aug_;
  x_.fill(0.);
  for (int i = 0; i < n_sigma; ++i) {
    x_ += weights_(i)*Xsig_pred_.col(i);
  }
  P_.fill(0.);
  for (int i = 0; i < n_sigma; ++i) {
    VectorXd diff = Xsig_pred_.col(i) - x_;
    diff(3) = std::fmod(diff(3), 2.*3.14159265359);
    P_ += weights_(i)*diff*diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_sigma = 1 + 2*n_aug_;
  MatrixXd Zsig = Xsig_pred_.topLeftCorner(2, n_sigma);
  VectorXd z = x_.head(2);
  MatrixXd S(2, 2);
  S.fill(0.);
  VectorXd diff(2); 
  for (int i = 0; i < n_sigma; ++i) {
    diff = Zsig.col(i) - z;
    S += weights_(i)*diff*diff.transpose();
  }
  S(0, 0) += std_laspx_*std_laspx_;
  S(1, 1) += std_laspy_*std_laspy_;
  double error = Update(meas_package, z, S, Zsig);
  ++n_lidar_;
  if (error > chi95_lidar_) ++n_over95_lidar_;
  std::cout << "LIDAR: " << double(n_over95_lidar_)/n_lidar_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_sigma = 1 + 2*n_aug_;
  MatrixXd Zsig(3, n_sigma);
  for (int i = 0; i < n_sigma; ++i) {
    double px  = Xsig_pred_(0, i);
    double py  = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double rho = std::sqrt(px*px + py*py);
    double phi = std::atan2(py, px);
    double rhod = (px*std::cos(yaw)*v + py*std::sin(yaw)*v)/rho;
    Zsig(0, i) = rho;
    Zsig(1, i) = phi;
    Zsig(2, i) = rhod;
  }
  VectorXd z(3);
  z.fill(0.);
  for (int i = 0; i < n_sigma; ++i) {
    z += weights_(i)*Zsig.col(i);
  }
  MatrixXd S(3, 3);
  S.fill(0.);
  VectorXd diff(3); 
  for (int i = 0; i < n_sigma; ++i) {
    diff = Zsig.col(i) - z;
    diff(1) = std::fmod(diff(1), 2.*3.14159265359);
    S += weights_(i)*diff*diff.transpose();
  }
  S(0, 0) += std_radr_*std_radr_;
  S(1, 1) += std_radphi_*std_radphi_;
  S(2, 2) += std_radrd_*std_radrd_;
  //std::cout << Xsig_pred_ << std::endl;
  //std::cout << z << std::endl;
  //std::cout << S << std::endl;
  double error = Update(meas_package, z, S, Zsig);
  std::cout << ": " << error << std::endl;
  ++n_radar_;
  if (error > chi95_radar_) ++n_over95_radar_;
  std::cout << "RADAR: " << double(n_over95_radar_)/n_radar_<< std::endl;
}

double UKF::Update(const MeasurementPackage& meas_package, 
    const VectorXd& z, const MatrixXd& S, const MatrixXd& Zsig) {
  int n_meas = z.size();
  int n_sigma = 1 + 2*n_aug_;
  MatrixXd T(n_x_, n_meas);
  T.fill(0.);
  const double PI2 = 2.*3.14159265359;
  VectorXd xdiff(n_x_);
  VectorXd zdiff(n_meas);
  for (int i = 0; i < n_sigma; ++i) {
    xdiff = Xsig_pred_.col(i) - x_;
    xdiff(3) = std::fmod(xdiff(3), PI2);
    zdiff = Zsig.col(i) - z;
    if (n_meas == 3) zdiff(1) = std::fmod(zdiff(1), PI2);
    T += weights_(i)*xdiff*zdiff.transpose();
  }
  MatrixXd Sinv = S.inverse();
  MatrixXd K = T*Sinv;
  zdiff = meas_package.raw_measurements_ - z;
  if (n_meas == 3) zdiff(1) = std::fmod(zdiff(1), PI2);
  x_ += K*zdiff;
  P_ -= K*S*K.transpose();
  return zdiff.transpose()*Sinv*zdiff;
}
