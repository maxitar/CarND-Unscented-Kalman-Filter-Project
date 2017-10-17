#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;
  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;
  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;
  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  ///* state covariance matrix
  MatrixXd P_;
  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;
  ///* time when the state is true, in us
  long long time_us_;
  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;
  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;
  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;
  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;
  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;
  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;
  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;
  ///* Weights of sigma points
  VectorXd weights_;
  ///* State dimension
  int n_x_;
  ///* Augmented state dimension
  int n_aug_;
  ///* Number of sigma points
  int n_sigma_;
  ///* Sigma point spreading parameter
  double lambda_;


  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const MeasurementPackage& meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const MeasurementPackage& meas_package);
private:
  const double PI2 = 2.*3.14159265359;
  /**
   * Generate the sigma points
   * @param x_aug The augmented mean state at k
   * @param P_aug The augmented covariance matrix at k
   */
  MatrixXd GenerateSigmaPoints(const VectorXd& x_aug, const MatrixXd& P_aug);
  /** 
   * Predict the sigma points
   * @param delta_t The time between the measurements at k+1 and k
   */
  void PredictSigmaPoints(double delta_t);
  /**
   * Combine common update operations for Lidar and Radar
   * @param meas_package The measurement at k+1
   * @param z The predicted measurement mean at k+1
   * @param S The predicted covariance at k+1
   * @param Zsig The individual predcited measurement points at k+1
   */
  double Update(const MeasurementPackage& meas_package, const VectorXd& z_pred, const MatrixXd& S, const MatrixXd& Zsig);
  ///* Number of radar measurements
  int n_radar_ = 0;
  ///* Number of radar measurements outside of 95-percentile
  int n_over95_radar_ = 0;
  ///* Number of lidar measurements
  int n_lidar_ = 0;
  ///* Number of lidar measurements outside of 95-percentile
  int n_over95_lidar_ = 0;
  ///* 95-percentile of chi^2 distribution with 3 dofs
  const double chi95_radar_ = 7.815;
  ///* 95-percentile of chi^2 distribution with 2 dofs
  const double chi95_lidar_ = 5.991;
};

#endif /* UKF_H */
