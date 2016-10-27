#include "ros/ros.h"
#include "sensor_msgs/Imu.h"

void imuCallback(const sensor_msgs::Imu::ConstPtr& msg){
  //ROS_INFO("Received IMU data.");

  //TODO put odometry stuff here

  //read: gyros, accel

  double ax = msg->linear_acceleration.x;
  double ay = msg->linear_acceleration.y;
  double w  = msg->angular_velocity.z;

  ROS_INFO("IMU data: ax=%f, ay=%f, w=%f",ax,ay,w);

}

int main(int argc, char** argv){
  ros::init(argc, argv, "imu_listner");
  ros::NodeHandle nh;

  ros::Subscriber imu_sub = nh.subscribe("fake_imu",1000,imuCallback);

  ros::spin();

  return 0;

}
