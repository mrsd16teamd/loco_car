/*
Executes most recent trajectory given by iLQR.
Sends feedback to action client (planner) in terms of "almost done", "done".
*/

#include "iLQR_executer.h"

iLQR_Executer::iLQR_Executer()
{
  ros::NodeHandle nh;
  cmd_pub = nh.advertise<geometry_msgs::Twist>("cmd_vel", 3);
  ROS_INFO("Started iLQR executer node.");
}

// provides action to execute plans
void iLQR::ExecuteTrajectory(int steps);

//     geometry_msgs::Twist msg;
//     msg.linear.x = vx;
//     msg.angular.z = steer;
//
//     ROS_INFO("Sending vx: %f, steer: %f", vx, steer);
//
//     cmd_pub.publish(msg);
//     ros::spinOnce();

int main(int argc, char** argv)
{
  ros::init(argc, argv, "iLQR_executer");

  iLQR_Executer executer;

  return 0;
}
