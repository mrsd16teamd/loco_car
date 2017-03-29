#ifndef _TRAJ_SERVER_H_
#define _TRAJ_SERVER_H_

#include "ros/ros.h"
#include <actionlib/server/simple_action_server.h>
#include <ilqr_loco/TrajExecAction.h>

#include <math.h>
#include <geometry_msgs/Twist.h>
#include <fstream>
#include <string>
#include <sstream>

//TODO is there a way for server to know who called it?

class TrajServer
{
public:
  TrajServer():
    as(nh, "traj_executer", boost::bind(&TrajServer::execute_trajectory, this,
    _1), false), traj_action("traj_executer"), old_msg_thres(0.5)
    {
      as.start();
      cmd_pub = nh.advertise<geometry_msgs::Twist>("cmd_vel", 3);
      ROS_INFO("Started iLQR executer node. Send me actions!");
    }

  void execute_trajectory(const ilqr_loco::TrajExecGoalConstPtr &goal);

private:
  ros::NodeHandle nh;
  actionlib::SimpleActionServer<ilqr_loco::TrajExecAction> as;
  std::string traj_action;

  // If a command was planned to be executed more than this many seconds, then
  // client will ignore it and move to processing next command.
  // Arbitrarily set to 0.5sec in constructor above.
  double old_msg_thres;

  // create messages that are used to published feedback/result
  ilqr_loco::TrajExecFeedback feedback;
  ilqr_loco::TrajExecResult result;

  ros::Publisher cmd_pub;
};

#endif
