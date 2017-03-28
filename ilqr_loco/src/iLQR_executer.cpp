/*
Executes most recent trajectory given by iLQR.
Sends feedback to action client (planner) in terms of "almost done", "done".
*/

#include "ros/ros.h"
#include <actionlib/server/simple_action_server.h>
#include <ilqr_loco/TrajExecAction.h">

#include <math.h>
#include <geometry_msgs/Twist.h>
#include <fstream>
#include <string>
#include <sstream>


class iLQR_Executer
{
private:
  ros::NodeHandle nh;
  actionlib::SimpleActionServer<ilqr_loco::TrajExecAction> as;
  std::string traj_action;

  // create messages that are used to published feedback/result
  ilqr_loco::TrajExecFeedback feedback_;
  ilqr_loco::TrajExecResult result_;

  ros::Publisher cmd_pub;

  void ExecuteTrajectory(int steps);

public:
  iLQR_Executer();
  void executeCB(const ilqr_loco::TrajExecGoalConstPtr &goal);
};

iLQR_Executer::iLQR_Executer()
{
  as(nh, name, boost::bind(&TrajExecAction::executeCB, this, _1), false), traj_action(name)
  {
    as.start();
  }

  cmd_pub = nh.advertise<geometry_msgs::Twist>("cmd_vel", 3);
  //TODO we can also subscribe to state information here for local control around traj

  ROS_INFO("Started iLQR executer node.");
}



// provides action to execute plans
void iLQR_Executer::ExecuteTrajectory(int steps){

};

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
