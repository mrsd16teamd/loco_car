#include "ros/ros.h"
#include "sensor_msgs/Imu.h"
#include "tf/transform_datatypes.h"
#include <tf/transform_broadcaster.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Twist.h>
#include <cmath>
#include <string>

#define LOWPASS false

class ImuData
{
public:
  ImuData() : orientation{ 0, 0, 0 }, angular_velocity{ 0, 0, 0 }, linear_acceleration{ 0, 0, 0 }
  {
  }
  float orientation[3];  // in roll, pitch yaw [radians]
  float angular_velocity[3];     //[rad/s]
  float linear_acceleration[3];  //[m/s^2]
};

class DriftChecker
{
public:
  DriftChecker();
  ~DriftChecker();

  // void InitializeParams();
  //void CalculateOdometry();

private:
  float wheelbase = 0.255;
  float drift_threshold_ = 0.1;

  ImuData latest_imu;
  float wz_expected;

  // Filter parameters
  float lp_alpha_;
  float lp_beta_;
  float threshold_;
  float inactive_timeout_;
  ros::Time last_active_time;
  ros::Time last_update_time;

  ros::Subscriber imu_sub;
  ros::Subscriber cmd_sub;
  //ros::Publisher odom_pub;

  void SensorCallback(const sensor_msgs::Imu::ConstPtr& msg);
  void CommandCallback(const geometry_msgs::Twist::ConstPtr& msg);
};

DriftChecker::DriftChecker()
{
  ros::NodeHandle nh;

  wz_expected = 0;

  imu_sub = nh.subscribe("imu", 1, &DriftChecker::SensorCallback, this);
  cmd_sub = nh.subscribe("cmd_vel", 1, &DriftChecker::CommandCallback, this);

  last_update_time = ros::Time::now();
  last_active_time = ros::Time::now();

  inactive_timeout_ = 1.0;
  drift_threshold_ = 0.3;

  ROS_INFO("started drift checker");
}

DriftChecker::~DriftChecker()
{
}

void DriftChecker::CommandCallback(const geometry_msgs::Twist::ConstPtr& msg)
{
  //update latest CommandCallback
  wz_expected = 1/wheelbase * msg->angular.z * msg->linear.x;
  //ROS_INFO("wz_expected: %f", wz_expected);

  last_update_time = ros::Time::now();

  //TODO account for latency between command and actuation?
}

void DriftChecker::SensorCallback(const sensor_msgs::Imu::ConstPtr& msg)
{
  // update latest_imu, call publish

  // orientation is originally in geometry_msgs/Quaternion
  tf::Quaternion quat;
  tf::quaternionMsgToTF(msg->orientation, quat);

  double roll, pitch, yaw;
  tf::Matrix3x3(quat).getRPY(roll, pitch, yaw);

  if (!LOWPASS)
  {
    latest_imu.orientation[0] = roll;
    latest_imu.orientation[1] = pitch;
    latest_imu.orientation[2] = yaw;

    latest_imu.angular_velocity[0] = msg->angular_velocity.x;
    latest_imu.angular_velocity[1] = msg->angular_velocity.y;
    latest_imu.angular_velocity[2] = msg->angular_velocity.z;

    latest_imu.linear_acceleration[0] = msg->linear_acceleration.x;
    latest_imu.linear_acceleration[1] = msg->linear_acceleration.y;
    latest_imu.linear_acceleration[2] = msg->linear_acceleration.z;
  }
  else
  {
    // Low pass filter
    latest_imu.orientation[0] = lp_alpha_ * roll + lp_beta_ * latest_imu.orientation[0];
    latest_imu.orientation[1] = lp_alpha_ * pitch + lp_beta_ * latest_imu.orientation[1];
    latest_imu.orientation[2] = lp_alpha_ * yaw + lp_beta_ * latest_imu.orientation[2];

    latest_imu.angular_velocity[0] = lp_alpha_ * (msg->angular_velocity.x) + lp_beta_ * latest_imu.angular_velocity[0];
    latest_imu.angular_velocity[1] = lp_alpha_ * (msg->angular_velocity.y) + lp_beta_ * latest_imu.angular_velocity[1];
    latest_imu.angular_velocity[2] = lp_alpha_ * (msg->angular_velocity.z) + lp_beta_ * latest_imu.angular_velocity[2];

    latest_imu.linear_acceleration[0] =
        lp_alpha_ * (msg->linear_acceleration.x) + lp_beta_ * (latest_imu.linear_acceleration[0]);
    latest_imu.linear_acceleration[1] =
        lp_alpha_ * (msg->linear_acceleration.y) + lp_beta_ * (latest_imu.linear_acceleration[1]);
    latest_imu.linear_acceleration[2] =
        lp_alpha_ * (msg->linear_acceleration.z) + lp_beta_ * (latest_imu.linear_acceleration[2]);
  }

   ROS_INFO("wz_sen: %f, wz_exp: %f",latest_imu.angular_velocity[2], wz_expected);

   // check for drift
   ros::Time current_time = ros::Time::now();
   if ((current_time.toSec() - last_active_time.toSec()) < inactive_timeout_){
     //no expectations
   }
   else{
    // ROS_INFO("drifting?");
     if(abs(wz_expected - latest_imu.angular_velocity[2]) > drift_threshold_){
       ROS_INFO("DRIFTING!!");
     }
     else{
       ROS_INFO("looks like you're not drifting");
     }
   }
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "drift_checker");

  DriftChecker drift_checker;

  ros::spin();

  return 0;
}
