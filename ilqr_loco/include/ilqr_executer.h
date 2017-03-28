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


public:
  iLQR_Executer();
  void execute_trajectory(const ilqr_loco::TrajExecGoalConstPtr &goal);
};
