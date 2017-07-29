//
// MIT License
//
// Copyright (c) 2017 MRSD Team D - LoCo
// The Robotics Institute, Carnegie Mellon University
// http://mrsdprojects.ri.cmu.edu/2016teamd/
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include "ros/ros.h"
#include "sensor_msgs/Imu.h"
#include "tf/transform_datatypes.h"
#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>
#include <nav_msgs/Odometry.h>
#include <cmath>

class ScanmatchAugment
{
public:
  ScanmatchAugment();
  ~ScanmatchAugment();

  // void InitializeParams();
  void CalculateOdometry();
  void PublishOdometry();

private:
  float wz_latest;
  tf::StampedTransform laserodom_tf;
  float alpha_ = 0.8;

  float x;
  float y;
  float theta;

  ros::Time last_imu_time;
  ros::Time last_laserodom_time;

  ros::Subscriber imu_sub;
  tf::TransformListener laserodom_listener;

  ros::Publisher odom_pub;
  tf::TransformBroadcaster odom_broadcaster;

  void SensorCallback(const sensor_msgs::Imu::ConstPtr& msg);
  void PublishOdometry(int &x, int &y, int &t);
};

ScanmatchAugment::ScanmatchAugment()
{
  ros::NodeHandle nh;

  //TODO add params here

  imu_sub = nh.subscribe("imu", 1, &ScanmatchAugment::SensorCallback, this);


  odom_pub = nh.advertise<nav_msgs::Odometry>("odom", 50);

  last_imu_time = ros::Time::now();  // TODO this isn't right.. where should it go?
  last_laserodom_time = ros::Time::now();
}

ScanmatchAugment::~ScanmatchAugment()
{
}

void ScanmatchAugment::PublishOdometry()
{
  // Publish stuff

  // since all odometry is 6DOF we need a quaternion created from yaw
  //TODO get t from method input
  geometry_msgs::Quaternion odom_quat = tf::createQuaternionMsgFromYaw(theta);

  ros::Time current_time = ros::Time::now();

  // Publish TF transform
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
  // ROS_INFO("sent base_link to odom transform");

  // Publish odometry message
  nav_msgs::Odometry odom_msg;
  odom_msg.header.stamp = current_time;
  odom_msg.header.frame_id = "odom";

  // set the position
  odom_msg.pose.pose.position.x = x;
  odom_msg.pose.pose.position.y = y;
  odom_msg.pose.pose.orientation = odom_quat;

  // set the velocity
  odom_msg.child_frame_id = "base_link";
  odom_msg.twist.twist.linear.x = 0.0;
  odom_msg.twist.twist.linear.y = 0.0;
  odom_msg.twist.twist.angular.z = wz_latest;
  //TODO check what happens if I just add angular velocity to odom message
  //      without complementary filter.
  //      Nothing should happen; odeometry model just samples around odom pose

  // publish the message
  odom_pub.publish(odom_msg);
}

void ScanmatchAugment::CalculateOdometry()
{
  // Update current_state, publish
  x = laserodom_tf.getOrigin().x();
  y = laserodom_tf.getOrigin().y();

  tf::Quaternion quat;
  quat = laserodom_tf.getRotation();
  double roll, pitch,yaw;
  tf::Matrix3x3(quat).getRPY(roll, pitch, yaw);
  theta = 0; //TODO how to get heading here?

  //TODO calculate time here. sample at fixed intervals?
  float h = 0.1;

  theta = alpha_*(theta + wz_latest*h) + (1-alpha_)*theta;

  PublishOdometry();
}

void ScanmatchAugment::SensorCallback(const sensor_msgs::Imu::ConstPtr& msg)
{
  wz_latest = msg->angular_velocity.z;

  try{
    //TODO check where remapping is occuring, change if necessary
    laserodom_listener.lookupTransform("/scanmatch_odom", "/base_link",
                             ros::Time(0), laserodom_tf);
  }
  catch (tf::TransformException ex){
    ROS_ERROR("%s",ex.what());
    ros::Duration(1.0).sleep();
  }

  CalculateOdometry();
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "scanmatch_augment");

  ScanmatchAugment imu_odometry_publisher;

  ros::spin();

  return 0;
}
