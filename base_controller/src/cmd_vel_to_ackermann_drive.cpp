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

#include <ros/ros.h>
#include <math.h>
#include <geometry_msgs/Twist.h>
#include <ackermann_msgs/AckermannDriveStamped.h>

class VelMsgConverter
{
public:
  VelMsgConverter()
  {
    // Initialize the publisher and subscriber
    vel_pub = n_.advertise<ackermann_msgs::AckermannDriveStamped>("/cmd_vel_ack", 10);
    vel_sub = n_.subscribe("/cmd_vel",10,&VelMsgConverter::velCallback,this);
  }

  private:
    ros::NodeHandle n_;
    ros::Publisher vel_pub;
    ros::Subscriber vel_sub;

    // Initialize parameters
    float wheelbase = 0.255;
    /*
    if (n_.getParam("chassis/wheelbase",wheelbase)) {
      n_.getParam("chassis/wheelbase",wheelbase);
      ROS_INFO("Loaded chassis parameters");
      ROS_INFO("wheelbase=%f", wheelbase);
    }
    else
      ROS_INFO("parameters not found!");
      */

/*
  float vel_to_steering_angle(float v, double omega) {
    if(omega==0 || v==0)
      return 0;
    float radius = v/omega;
    float steering_angle = atan(wheelbase/radius);
    if(v<0){
	steering_angle *= -1;
    }
    return steering_angle;
  }
*/

  float vel_to_steering_angle(float v, double omega) {
     return omega;
  }

  void velCallback(const geometry_msgs::Twist& twist_msg) {
    ackermann_msgs::AckermannDriveStamped ack_msg;

    float wheelbase = 0.255;
    float lin_vel = twist_msg.linear.x;
    ack_msg.header.stamp = ros::Time::now();
    ack_msg.header.frame_id = "base_link";

    ack_msg.drive.speed = lin_vel;
    ack_msg.drive.steering_angle = vel_to_steering_angle(lin_vel, twist_msg.angular.z);

    vel_pub.publish(ack_msg);

  }
};//End of class

int main(int argc, char** argv) {
  ros::init(argc, argv, "cmd_vel_to_ackermann_drive");
  VelMsgConverter ackermann_drive_pub;
  ros::spin();
  return 0;
}
