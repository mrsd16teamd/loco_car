#ifndef IMU_ODOMPUB_H
#define IMU_ODOMPUB_H

#include "sensor_msgs/Imu.h"
#include "tf/transform_datatypes.h"
#include <tf/transform_broadcaster.h>
#include <nav_msgs/Odometry.h>
#include <cmath>

struct ImuData
{
  // Store (possibly filtered) state of IMU
  
  ImuData(): orientation({0,0,0}), angular_velocity({0,0,0}),
    linear_acceleration({0,0,0}) { }
  float orientation[3]; // rpy[radians], pre-processed onboard with Madgwick's
  float angular_velocity[3];     // [rad/s]
  float linear_acceleration[3];  // [m/s^2]
};

struct State2d
{
  State2d() : x(0), y(0), t(0), vx(0), vy(0), w(0) {}
  float x, y, t, vx, vy, w;
};

class ImuOdometryPublisher
{
public:
  ImuOdometryPublisher();
  ~ImuOdometryPublisher();

private:
  ImuData latest_imu;
  State2d current_state;

  float lp_alpha_;  // Parameter for low-pass-ness. 0.0=passthru, 0.9 is LP
  float lp_beta_;   // Derived from lp_alpha_
  float threshold_; // Threshold for imu acceleration to be non-zero (not noise)
  float inactive_timeout_; // Timer to zero out velocity.

  ros::Time last_active_time; // for imu
  ros::Time last_update_time; // for imu

  ros::Subscriber imu_sub;
  ros::Publisher odom_pub;
  tf::TransformBroadcaster odom_broadcaster;

  //All steps called from SensorCallback
  void SensorCallback(const sensor_msgs::Imu::ConstPtr& msg);
  void UpdateIMU(const sensor_msgs::Imu::ConstPtr& msg);
  void CalculateOdometry();
  void PublishOdometry();
};

#endif IMU_ODOMPUB_H
