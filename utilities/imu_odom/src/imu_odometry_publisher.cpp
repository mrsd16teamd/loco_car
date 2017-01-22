/*
 * imu_odometry_publisher
 * Reads IMU messages and integrates them to provide odometry frame for AMCL.
 * Will drift over time, but just need to provide continuous reference frame.
 * Digital low-pass filter available for accelerometer readings.
 * Tested on Razor 9DOF IMU
 *
 * Author: Kazu Otani, kazuotani14@gmail.com
 */

#include "imu_odometry_publisher.h"
#include "ros/ros.h"

#define LOWPASS true

ImuOdometryPublisher::ImuOdometryPublisher()
{
  ros::NodeHandle nh;

  // Initialize parameters from yaml file
  if (nh.getParam("imu/lp_alpha", lp_alpha_) && nh.getParam("imu/threshold", threshold_) &&
      nh.getParam("odom/inactive_timeout", inactive_timeout_))
  {
    lp_beta_ = 1 - lp_alpha_;
    ROS_INFO("Loaded imu_odometry parameters:");
    ROS_INFO("alpha=%f, threshold=%f, inactive_t=%f", lp_alpha_, threshold_, inactive_timeout_);
  }
  else
  {
    ROS_INFO("Some or all parameters not found! Resorting to defaults.");
    lp_alpha_ = 0.9;
    threshold_ = 0.1;
    inactive_timeout_ = 0.5;
  }

  imu_sub = nh.subscribe("imu", 1, &ImuOdometryPublisher::SensorCallback, this);
  odom_pub = nh.advertise<nav_msgs::Odometry>("odom", 50);

  last_update_time = ros::Time::now();
  last_active_time = ros::Time::now();
}

ImuOdometryPublisher::~ImuOdometryPublisher()
{
}

void ImuOdometryPublisher::PublishOdometry()
{
  // Publish stuff

  // Since all odometry is 6DOF we need a quaternion created from yaw
  geometry_msgs::Quaternion odom_quat = tf::createQuaternionMsgFromYaw(current_state.t);

  ros::Time current_time = ros::Time::now();

  // TF transform
  geometry_msgs::TransformStamped odom_tf;
  odom_tf.header.stamp = current_time;
  odom_tf.header.frame_id = "odom";
  odom_tf.child_frame_id = "base_link";
  odom_tf.transform.translation.x = current_state.x;
  odom_tf.transform.translation.y = current_state.y;
  odom_tf.transform.translation.z = 0.0;
  odom_tf.transform.rotation = odom_quat;

  // send the transform
  odom_broadcaster.sendTransform(odom_tf);

  //  Odometry message
  nav_msgs::Odometry odom_msg;
  odom_msg.header.stamp = current_time;
  odom_msg.header.frame_id = "odom";

  // Set the position
  odom_msg.pose.pose.position.x = current_state.x;
  odom_msg.pose.pose.position.y = current_state.y;
  odom_msg.pose.pose.orientation = odom_quat;

  // Set the velocity
  odom_msg.child_frame_id = "base_link";
  odom_msg.twist.twist.linear.x = current_state.vx;
  odom_msg.twist.twist.linear.y = current_state.vy;
  odom_msg.twist.twist.angular.z = current_state.w;

  // Publish the message
  odom_pub.publish(odom_msg);
} // PublishOdometry

void ImuOdometryPublisher::CalculateOdometry()
{
  // Update current_state

  ros::Time current_time = ros::Time::now();
  float h = current_time.toSec() - last_update_time.toSec();

  // Update position, heading
  current_state.x += current_state.vx * h;
  current_state.y += current_state.vy * h;
  current_state.t += current_state.w * h;

  // Update velocities, compensating for tilt
  float ax = latest_imu.linear_acceleration[0] * cos(latest_imu.orientation[1]) +
             latest_imu.linear_acceleration[2] * sin(latest_imu.orientation[1]);
  float ay = latest_imu.linear_acceleration[1] * cos(-latest_imu.orientation[0]) +
             latest_imu.linear_acceleration[2] * sin(-latest_imu.orientation[0]);

  current_state.vx += ax * h;
  current_state.vy += ay * h;
  current_state.w = latest_imu.angular_velocity[2];
  // w_z is approximate; it doesn't take into account orientation

  // Check if robot is moving above threshold amount
  // So small static noise doesn't move odom frame when robot is stationary.
  if ((std::abs(ax) > threshold_) || (std::abs(ay) > threshold_) ||
      (std::abs(latest_imu.angular_velocity[2]) > threshold_))
  {
    last_active_time = ros::Time::now();
  }

  // Check if robot has been stationary -> zero velocity to reduce drift.
  if ((current_time.toSec() - last_active_time.toSec()) > inactive_timeout_)
  {
    current_state.vx = current_state.vy = current_state.w = 0;
  }

  last_update_time = current_time;
} // CalculateOdometry

void ImuOdometryPublisher::UpdateIMU(const sensor_msgs::Imu::ConstPtr& msg){
  // Update latest_imu data

  // Orientation is originally in geometry_msgs/Quaternion. convert
  tf::Quaternion quat;
  tf::quaternionMsgToTF(msg->orientation, quat);

  double roll, pitch, yaw;
  tf::Matrix3x3(quat).getRPY(roll, pitch, yaw);

  latest_imu.orientation[0] = roll;
  latest_imu.orientation[1] = pitch;
  latest_imu.orientation[2] = yaw;

  latest_imu.angular_velocity[0] = msg->angular_velocity.x;
  latest_imu.angular_velocity[1] = msg->angular_velocity.y;
  latest_imu.angular_velocity[2] = msg->angular_velocity.z;

  // Accelerometer is the most problematic; may need to be low-pass filtered
  if (!LOWPASS)
  {
    latest_imu.linear_acceleration[0] = msg->linear_acceleration.x;
    latest_imu.linear_acceleration[1] = msg->linear_acceleration.y;
    latest_imu.linear_acceleration[2] = msg->linear_acceleration.z;
  }
  else
  {
    // Low pass filter
    latest_imu.linear_acceleration[0] =
        lp_alpha_ * (msg->linear_acceleration.x) + lp_beta_ * (latest_imu.linear_acceleration[0]);
    latest_imu.linear_acceleration[1] =
        lp_alpha_ * (msg->linear_acceleration.y) + lp_beta_ * (latest_imu.linear_acceleration[1]);
    latest_imu.linear_acceleration[2] =
        lp_alpha_ * (msg->linear_acceleration.z) + lp_beta_ * (latest_imu.linear_acceleration[2]);
  }
} // UpdateIMU

void ImuOdometryPublisher::SensorCallback(const sensor_msgs::Imu::ConstPtr& msg)
{
  UpdateIMU(msg);
  CalculateOdometry();
  PublishOdometry();
} // SensorCallback

int main(int argc, char** argv)
{
  ros::init(argc, argv, "imu_odometry_publisher");

  ImuOdometryPublisher imu_odometry_publisher;

  ros::spin();

  return 0;
}
