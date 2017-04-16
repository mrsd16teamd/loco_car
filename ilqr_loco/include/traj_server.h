#ifndef _TRAJ_SERVER_H_
#define _TRAJ_SERVER_H_

#include "ros/ros.h"
#include "Eigen/Core"
#include <std_msgs/Float64MultiArray.h>
#include <actionlib/server/simple_action_server.h>
#include <ilqr_loco/TrajExecAction.h>
#include <tf/transform_datatypes.h>

#include <math.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>
#include <fstream>
#include <string>
#include <sstream>


//TODO is there a way for server to know who called it?

class TrajServer
{
public:
  TrajServer():
    as(nh, "traj_server", boost::bind(&TrajServer::execute_trajectory, this,
    _1), false), traj_action("traj_server")
    {
      as.start();

      cmd_pub = nh.advertise<geometry_msgs::Twist>("cmd_vel", 3);
      path_pub = nh.advertise<nav_msgs::Path>("path", 1);
      state_sub_  = nh.subscribe("odometry/filtered", 1, &TrajServer::stateCb, this);
      LoadParams();

      cur_state.resize(8);
      L.resize(2,8);
      x.resize(8);
      dx.resize(8);

      ROS_INFO("Started iLQR executer node. Send me actions!");
    }

private:
  ros::NodeHandle nh;
  actionlib::SimpleActionServer<ilqr_loco::TrajExecAction> as;
  std::string traj_action;
  ros::Subscriber state_sub_;

  // If a command was planned to be executed more than this many seconds, then
  // client will ignore it and move to processing next command.
  double old_msg_thres;
  Eigen::VectorXd cur_state;
  Eigen::Vector2d l;
  Eigen::MatrixXd L;
  Eigen::VectorXd x;
  Eigen::VectorXd dx;
  Eigen::Vector2d u;
  Eigen::Vector2d last_u;

  // create messages that are used to published feedback/result
  ilqr_loco::TrajExecFeedback feedback;
  ilqr_loco::TrajExecResult result;

  ros::Publisher cmd_pub;
  ros::Publisher path_pub;
  ros::Subscriber state_sub;

  void LoadParams();
  void execute_trajectory(const ilqr_loco::TrajExecGoalConstPtr &goal);
  void stateCb(const nav_msgs::Odometry &msg);
  void FillVecFromOdom(const nav_msgs::Odometry &odom, Eigen::VectorXd &v);
  void FillVecFromTwist(const geometry_msgs::Twist &twist, Eigen::Vector2d &v);
  void FillTwistFromVec(geometry_msgs::Twist &twist, const Eigen::Vector2d &v);
};

#endif
