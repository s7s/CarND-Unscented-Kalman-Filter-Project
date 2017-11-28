#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;

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
  //state covariance matrix P
  P_ << 1, 0, 0, 0,0,
        0, 1, 0, 0,0,
        0, 0, 1, 0,0,
        0, 0, 0, 1,0,
        0, 0, 0, 0,1;


  x_ << 1, 1, 0, 0,0;

  n_x_=x_.size();

  n_aug_=x_.size()+2;

  lambda_=3-n_x_;


  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //create vector for weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0)=lambda_/(lambda_+n_aug_);
  weights_.tail(2*n_aug_)=VectorXd::Ones(2*n_aug_)*(0.5/(lambda_+n_aug_));


  long x_size = x_.size();
  I_ = MatrixXd::Identity(x_size, x_size);

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
float UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    time_us_=meas_package.timestamp_;
    cout << "UKF: " << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_(0)=meas_package.raw_measurements_(0) *cos(meas_package.raw_measurements_(1));
      x_(1)=-meas_package.raw_measurements_(0) *sin(meas_package.raw_measurements_(1));
      x_(3)=atan2(x_(1),x_(0));
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_(0)=meas_package.raw_measurements_(0);
      x_(1)=meas_package.raw_measurements_(1);
      x_(3)=atan2(x_(1),x_(0));
    }
    ////long x_size = ekf_.x_.size();
    ////ekf_.I_ = MatrixXd::Identity(x_size, x_size);

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return 0;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  float NIS;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    if(use_radar_){
      NIS=UpdateRadar(meas_package);
    }
  }
  else {
    // Laser updates
    if(use_laser_){
      NIS=UpdateLidar(meas_package);
    }
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
return NIS;

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


  //generate augmented sigma points

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug= VectorXd::Zero(n_aug_);
  x_aug.head(n_x_)=x_;
  P_aug=MatrixXd::Zero(n_aug_,n_aug_);
  P_aug.block(0,0,n_x_,n_x_)=P_;
  P_aug(n_x_,n_x_)=std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1)=std_yawdd_*std_yawdd_;
  MatrixXd A=P_aug.llt().matrixL();

  Xsig_aug.col(0)=x_aug;
  for(int i=0;i<n_aug_;i++){
    Xsig_aug.col(i+1)=x_aug+(sqrt(lambda_+n_aug_)*A.col(i));
    Xsig_aug.col(i+n_aug_+1)=x_aug-(sqrt(lambda_+n_aug_)*A.col(i));
  }


  //predict sigma points

  double delta_t_2 =delta_t*delta_t;
  for (int i=0;i<2*n_aug_+1;i++){
    if(Xsig_aug(4,i)==0){
      Xsig_pred_(0,i)=Xsig_aug(0,i)+Xsig_aug(2,i)*cos(Xsig_aug(3,i))*delta_t+0.5*delta_t_2*cos(Xsig_aug(3,i))*Xsig_aug(5,i);
      Xsig_pred_(1,i)=Xsig_aug(1,i)+Xsig_aug(2,i)*sin(Xsig_aug(3,i))*delta_t+0.5*delta_t_2*sin(Xsig_aug(3,i))*Xsig_aug(5,i);
    }
    else{
      Xsig_pred_(0,i)=Xsig_aug(0,i)+(Xsig_aug(2,i)/Xsig_aug(4,i))*(sin(Xsig_aug(3,i)+Xsig_aug(4,i)*delta_t)-sin(Xsig_aug(3,i)))+0.5*delta_t_2*cos(Xsig_aug(3,i))*Xsig_aug(5,i);
      Xsig_pred_(1,i)=Xsig_aug(1,i)+(Xsig_aug(2,i)/Xsig_aug(4,i))*(-cos(Xsig_aug(3,i)+Xsig_aug(4,i)*delta_t)+cos(Xsig_aug(3,i)))+0.5*delta_t_2*sin(Xsig_aug(3,i))*Xsig_aug(5,i);
    }
    Xsig_pred_(2,i)=Xsig_aug(2,i)+0+delta_t*Xsig_aug(5,i);
    Xsig_pred_(3,i)=Xsig_aug(3,i)+Xsig_aug(4,i)*delta_t+0.5*delta_t_2*Xsig_aug(6,i);
    Xsig_pred_(4,i)=Xsig_aug(4,i)+0+delta_t*Xsig_aug(6,i);
  }


  //predict mean and covarience


  x_=Xsig_pred_*weights_;

  VectorXd x_diff = VectorXd(n_x_);
  P_=MatrixXd::Zero(n_x_,n_x_);
  for (int i=0;i<2*n_aug_+1;i++){
    x_diff=Xsig_pred_.col(i)-x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    P_=P_+weights_(i)*(x_diff)*(x_diff).transpose();
  }


}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
float UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
 float NIS;
 MatrixXd H(2,5);
  H << 1, 0, 0, 0,0,
      0, 1, 0, 0,0;
  VectorXd z= meas_package.raw_measurements_;
  VectorXd y=z-H*x_;

  MatrixXd S=H*P_*H.transpose();
  S(0,0)+=std_laspx_*std_laspx_;
  S(1,1)+=std_laspy_*std_laspy_;
  MatrixXd K=P_*H.transpose()*S.inverse();
  x_=x_+K*y;

  P_ = (I_ - K * H) * P_;
  NIS=y.transpose()*S.inverse()*y;
  return NIS;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
float UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //create matrix for sigma points in measurement space
  float NIS;
  int n_z = 3;
  VectorXd z= meas_package.raw_measurements_;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  VectorXd z_diff = VectorXd(n_z);

  for(int i=0 ;i<2*n_aug_+1;i++){
    Zsig(0,i)=sqrt(Xsig_pred_(0,i)*Xsig_pred_(0,i)+Xsig_pred_(1,i)*Xsig_pred_(1,i));
    Zsig(1,i)=atan2(Xsig_pred_(1,i),Xsig_pred_(0,i));
    Zsig(2,i)=((Xsig_pred_(0,i)*cos(Xsig_pred_(3,i))+Xsig_pred_(1,i)*sin(Xsig_pred_(3,i)))*Xsig_pred_(2,i))/Zsig(0,i);

  }

  z_pred=Zsig*weights_;
  S=MatrixXd::Zero(n_z,n_z);
  for(int i=0;i<2*n_aug_+1;i++){
    z_diff=Zsig.col(i)-z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S=S+weights_(i)*(z_diff)*(z_diff).transpose();
  }
  S(0,0)+=std_radr_*std_radr_;
  S(1,1)+=std_radphi_*std_radphi_;
  S(2,2)+=std_radrd_*std_radrd_;



  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc=MatrixXd::Zero(n_x_, n_z);
  VectorXd x_diff = VectorXd(n_x_);

  for(int i=0;i<2*n_aug_+1;i++){
    z_diff=Zsig.col(i)-z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    x_diff=Xsig_pred_.col(i)-x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc=Tc+weights_(i)*(x_diff)*(z_diff).transpose();
  }
  MatrixXd K=Tc*S.inverse();
  VectorXd y=z-z_pred;
  while (y(1)> M_PI) y(1)-=2.*M_PI;
  while (y(1)<-M_PI) y(1)+=2.*M_PI;

  x_=x_+K*y;

  P_=P_-K*S*K.transpose();
  NIS=y.transpose()*S.inverse()*y;
  return NIS;
}
