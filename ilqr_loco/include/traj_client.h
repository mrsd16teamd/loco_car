#ifndef _TRAJ_CLIENT_H_
#define _TRAJ_CLIENT_H_

#include "ros/ros.h"
#include <actionlib/client/simple_action_client.h>
#include <actionlib/client/terminal_state.h>
#include <ilqr_loco/TrajExecAction.h>

#include <math.h>
#include <geometry_msgs/Twist.h>
#include <std_msgs/Float32MultiArray.h>

#include <fstream>
#include <string>
#include <sstream>

class TrajClient
{
public:
  TrajClient(): ac("traj_executer", true)
  {
    state_sub  = nh.subscribe("odometry/filtered", 1, &TrajClient::stateCb, this);
    obs_sub = nh.subscribe("ccs", 1, &TrajClient::obsCb, this);

    ROS_INFO("Waiting for action server to start.");
    ac.waitForServer(); //will wait for infinite time
    ROS_INFO("Action server started.");
    Plan();
  }

private:
  ros::NodeHandle nh;
  ros::Publisher cmd_pub;
  ros::Subscriber state_sub;
  ros::Subscriber obs_sub;

  actionlib::SimpleActionClient<ilqr_loco::TrajExecAction> ac;

  nav_msgs::Odometry most_recent_state;
  std_msgs::Float32MultiArray obs_pos;

  void Plan();
  ilqr_loco::TrajExecGoal GenerateTrajectory();
  void SendTrajectory(ilqr_loco::TrajExecGoal &goal);

  void activeCb();
  void feedbackCb(const ilqr_loco::TrajExecFeedbackConstPtr& feedback);
  void doneCb(const actionlib::SimpleClientGoalState& state,
              const ilqr_loco::TrajExecResultConstPtr& result);

  void stateCb(const nav_msgs::Odometry &msg);
  void obsCb(const std_msgs::Float32MultiArray &msg);
};

#endif
