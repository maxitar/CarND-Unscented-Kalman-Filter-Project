#include <iostream>
#include <fstream>
#include <math.h>
#include "ukf.h"
#include "tools.h"

int file_input(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << "Need input text file" << std::endl;
    return 1;
  }
  // Create a Kalman Filter instance
  UKF ukf;

  // used to compute the RMSE later
  Tools tools;
  std::vector<VectorXd> estimations;
  std::vector<VectorXd> ground_truth;

  std::ifstream input_file(argv[1]);
  if (!input_file) { std::cout << "Couldn't open input file" << std::endl; return 1; }

  std::ofstream output_file("output.csv");
  if (!output_file) { std::cout << "Couldn't open output file" << std::endl; return 1; }
  output_file << "Type\tTime\tpx_est\tpx_gt\tpy_est\tpy_gt\tvx_est\tvx_gt\tvy_est\tvy_gt\tyaw_est\tyaw_gt\tyawd_est\tyawd_gt\tNIS\n";
  std::string sensor_measurment;
  MeasurementPackage meas_package;
  long long timestamp_0 = -1;
  while (std::getline(input_file, sensor_measurment)) { 
    std::istringstream iss(sensor_measurment);
    long long timestamp;

    // reads first element from the current line
    std::string sensor_type;
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0) {
      if (ukf.is_initialized_ == true && ukf.use_laser_ == false) {
        VectorXd RMSE = tools.CalculateRMSE(estimations, ground_truth);
        continue;
      }
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float px;
      float py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
    } else if (sensor_type.compare("R") == 0) {
      if (ukf.is_initialized_ == true && ukf.use_radar_ == false) {
        VectorXd RMSE = tools.CalculateRMSE(estimations, ground_truth);
        continue;
      }
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float theta;
      float ro_dot;
      iss >> ro;
      iss >> theta;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro,theta, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
    }
    if (timestamp_0 == -1) timestamp_0 = timestamp;
    double px_gt;
    double py_gt;
    double vx_gt;
    double vy_gt;
    double yaw_gt;
    double yawd_gt;
    iss >> px_gt;
    iss >> py_gt;
    iss >> vx_gt;
    iss >> vy_gt;
    iss >> yaw_gt;
    iss >> yawd_gt;

    //Call ProcessMeasurment(meas_package) for Kalman filter
    double NIS = ukf.ProcessMeasurement(meas_package);    	  

    double px = ukf.x_(0);
    double py = ukf.x_(1);
    double v   = ukf.x_(2);
    double yaw = ukf.x_(3);
    double yawd = ukf.x_(4);
    double vx = cos(yaw)*v;
    double vy = sin(yaw)*v;

    double time = (timestamp - timestamp_0)/1e6;
    output_file << sensor_type << '\t' << time << '\t' <<
                   px << '\t' << px_gt << '\t' <<
                   py << '\t' << py_gt << '\t' <<
                   vx << '\t' << vx_gt << '\t' <<
                   vy << '\t' << vy_gt << '\t' <<
                   yaw << '\t' << yaw_gt << '\t' <<
                   yawd << '\t' << yawd_gt << '\t' << NIS << '\n';
  }
  output_file.close();
  input_file.close();
}
