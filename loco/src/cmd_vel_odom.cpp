
#include "ros/ros.h"
#include <math.h>
#include "tf/transform_datatypes.h"
#include <tf/transform_broadcaster.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Twist.h>

class cmdVelOdomPublisher
{
public:
  cmdVelOdomPublisher(): x(0), y(0), theta(0), vx(0), vy(0), w(0){
    ros::NodeHandle nh;

    cmd_sub = nh.subscribe("cmd_vel", 1, &cmdVelOdomPublisher::UpdateOdom, this);
    odom_pub = nh.advertise<nav_msgs::Odometry>("odom", 50);

    last_update_time = ros::Time::now();

    ROS_INFO("Started cmd_vel odom publisher");
  }

private:
  static const double wheelbase = 0.255;
  static const double cmd_timesteps = 0.1;

  double x, y, theta, vx, vy, w;

  ros::Subscriber cmd_sub;
  ros::Publisher odom_pub;
  tf::TransformBroadcaster odom_broadcaster;

  ros::Time last_update_time;

  void UpdateOdom(const geometry_msgs::Twist::ConstPtr& msg);
  void PublishOdom();
};

void cmdVelOdomPublisher::UpdateOdom(const geometry_msgs::Twist::ConstPtr& msg)
{
  //Read commanded velocity, and calculate expected angular velocity
  float wz_expected = 0.0;
  if (msg->linear.x > 0.1)
    wz_expected = 1/wheelbase *tan( msg->angular.z) * msg->linear.x;

  // Guess time until next odom update
  ros::Time current_time = ros::Time::now();
  double dt = current_time.toSec() - last_update_time.toSec(); //timestep
  if (dt>1.0) dt = cmd_timesteps;
  last_update_time = ros::Time::now();

  // Update odom frame
  x += msg->linear.x * cos(theta) * dt;
  y += msg->linear.x * sin(theta) * dt;
  theta += wz_expected * dt;

  vx = msg->linear.x;
  w  = wz_expected;

  // Publish
  PublishOdom();
}

void cmdVelOdomPublisher::PublishOdom()
{
  // Since all odometry is 6DOF we need a quaternion created from yaw
  geometry_msgs::Quaternion odom_quat = tf::createQuaternionMsgFromYaw(theta);

  ros::Time current_time = ros::Time::now();

  // TF transform
  geometry_msgs::TransformStamped odom_tf;
  odom_tf.header.stamp = current_time;
  odom_tf.header.frame_id = "odom";
  odom_tf.child_frame_id = "base_link";
  odom_tf.transform.translation.x = x;
  odom_tf.transform.translation.y = y;
  odom_tf.transform.translation.z = 0.0;
  odom_tf.transform.rotation = odom_quat;

  // send the transform
  odom_broadcaster.sendTransform(odom_tf);

  //  Odometry message
  nav_msgs::Odometry odom_msg;
  odom_msg.header.stamp = current_time;
  odom_msg.header.frame_id = "odom";

  // Set the position
  odom_msg.pose.pose.position.x = x;
  odom_msg.pose.pose.position.y = y;
  odom_msg.pose.pose.orientation = odom_quat;

  // Set the velocity. Velocity can just be zero.
  odom_msg.child_frame_id = "base_link";
  odom_msg.twist.twist.linear.x = vx;
  odom_msg.twist.twist.linear.y = vy;
  odom_msg.twist.twist.angular.z = w;

  // Publish the message
  odom_pub.publish(odom_msg);
} // PublishOdometry


int main(int argc, char** argv)
{
  ros::init(argc, argv, "cmd_vel_odom");

  cmdVelOdomPublisher cmd_vel_odom;

  ros::spin();

  return 0;
}
