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
#include "tf/transform_datatypes.h"
#include "move_base_msgs/MoveBaseActionGoal.h"
#include <tf/transform_listener.h>
#include "std_msgs/Float64.h"
#include <cmath>
#include <string>
#include "geometry_msgs/PoseWithCovarianceStamped.h"
#include "geometry_msgs/Quaternion.h"


ros::Publisher pose_pub;

int main(int argc, char** argv)
{
  ros::init(argc, argv, "amcl_echo_node");

  ros::NodeHandle nh;

  pose_pub = nh.advertise<geometry_msgs::PoseWithCovarianceStamped>("amcl_pose_echo",1);

  tf::TransformListener listener;

  // TODO change this to wait until amcl_pose is publishing
  ROS_INFO("amcl_echo: Waiting for everything to start.");
  ros::Duration(5).sleep();
  ROS_INFO("Publishing amcl_pose_echo now.");

  ros::Rate r(40);

  while(ros::ok()){
    ros::spinOnce();

    double current_x, current_y;

    tf::StampedTransform transform;
    tf::Quaternion q;

    try{
       listener.lookupTransform("map", "base_link", ros::Time(0), transform);
       current_x = transform.getOrigin().x();
       current_y = transform.getOrigin().y();
       q = transform.getRotation();

       //ROS_INFO("Currently at: %f , %f ", current_x, current_y);
    }
    catch (tf::TransformException ex)
    {
      ROS_ERROR("amcl_echo _node: map to base_link transform not found");
      ros::Duration(1).sleep();
      continue;
    }

  	geometry_msgs::PoseWithCovarianceStamped pose_msg;
  	pose_msg.pose.pose.position.x = current_x;
    pose_msg.pose.pose.position.y = current_y;
    pose_msg.pose.pose.position.z = 0.0;

    double yaw = getYaw(q);
    geometry_msgs::Quaternion q_msg = tf::createQuaternionMsgFromYaw(yaw);
    pose_msg.pose.pose.orientation = q_msg;

  	pose_msg.header.stamp = ros::Time::now();
	pose_msg.header.frame_id = "map";
  	pose_pub.publish(pose_msg);

    r.sleep();
  }

  return 0;
}
