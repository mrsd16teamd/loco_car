#ifndef _TRAJ_SERVER_H_
#define _TRAJ_SERVER_H_

#include "ros/ros.h"
#include <actionlib/server/simple_action_server.h>
#include <ilqr_loco/TrajExecAction.h>

#include <math.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/PoseStamped.h>
#include <tf/transform_datatypes.h>
#include <nav_msgs/Path.h>
#include <fstream>
#include <string>
#include <sstream>

#include "try_get_param.h"

//TODO is there a way for server to know who called it?

class TrajServer
{
public:
  TrajServer():
    as(nh, "traj_server", boost::bind(&TrajServer::execute_trajectory, this,
    _1), false), traj_action("traj_server"), cur_yaw_(0), cur_integral_(0),
    prev_error_(0), dt(0.02)
    {
	    ROS_INFO("Starting traj server.");

      as.start();
      cmd_pub = nh.advertise<geometry_msgs::Twist>("cmd_vel", 3);
      path_pub = nh.advertise<nav_msgs::Path>("path", 1);
      state_sub  = nh.subscribe("odometry/filtered", 1, &TrajServer::stateCb, this);

      LoadParams();
      ROS_INFO("Started iLQR executer node. Send me actions!");
    }

private:
  ros::NodeHandle nh;
  actionlib::SimpleActionServer<ilqr_loco::TrajExecAction> as;

  std::string traj_action;

  // If a command was planned to be executed more than this many seconds, then
  // client will ignore it and move to processing next command.
  double old_msg_thres;

  // create messages that are used to publish feedback/result
  ilqr_loco::TrajExecFeedback feedback_;
  ilqr_loco::TrajExecResult result_;

  //PID control
  double cur_yaw_;
  double cur_integral_;
  double prev_error_;
  float kp_, ki_, kd_, dt;

  ros::Publisher cmd_pub;
  ros::Publisher path_pub;
  ros::Subscriber state_sub;

  void LoadParams();
  void execute_trajectory(const ilqr_loco::TrajExecGoalConstPtr &goal);
  geometry_msgs::Twist pid_correct_yaw(geometry_msgs::Twist orig_twist, nav_msgs::Odometry state);
  void stateCb(const nav_msgs::Odometry &msg);
  void PublishPath(const ilqr_loco::TrajExecGoalConstPtr &goal);

};

#endif
