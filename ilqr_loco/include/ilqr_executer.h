#include "ros/ros.h"
#include <actionlib/server/simple_action_server.h>
#include <ilqr_loco/TrajExecAction.h>

#include <math.h>
#include <geometry_msgs/Twist.h>
#include <fstream>
#include <string>
#include <sstream>

class iLQR_Executer
{
protected:
  ros::NodeHandle nh;
  actionlib::SimpleActionServer<ilqr_loco::TrajExecAction> as;
  std::string traj_action;

  // create messages that are used to published feedback/result
  ilqr_loco::TrajExecFeedback feedback;
  ilqr_loco::TrajExecResult result;

  ros::Publisher cmd_pub;

public:
  iLQR_Executer():
    as(nh, "traj_executer", boost::bind(&iLQR_Executer::execute_trajectory, this,
    _1), false), traj_action("traj_executer")
    {
      as.start();
      cmd_pub = nh.advertise<geometry_msgs::Twist>("cmd_vel", 3);
      ROS_INFO("Started iLQR executer node.");
    }

  void execute_trajectory(const ilqr_loco::TrajExecGoalConstPtr &goal);
};
