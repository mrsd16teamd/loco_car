#ifndef _RAMP_CLIENT_H_
#define _RAMP_CLIENT_H_

#include "ros/ros.h"
#include <actionlib/client/simple_action_client.h>
#include <actionlib/client/terminal_state.h>
#include <ilqr_loco/TrajExecAction.h>

// Do we need these?
#include <string>

class RampPlanner
{
public:
  void Plan();
  bool flag_;

  RampPlanner() : ac("traj_executer", true), flag_(false), T_(0) {
    state_sub  = nh.subscribe("odometry/filtered", 1, &RampPlanner::stateCb, this);

    ROS_INFO("Waiting for action server to start.");
    ac.waitForServer(); //will wait for infinite time
    ROS_INFO("Action server started.");
    Plan();
  }

protected:
  ros::NodeHandle nh_;
  ros::Subscriber state_sub_;
  actionlib::SimpleActionClient<ilqr_loco::TrajExecAction> ac_;

  nav_msgs::Odometry cur_state_;
  nav_msgs::Odometry prev_state_;

  geometry_msgs::Vector3 cur_rpy_;
  geometry_msgs::Vector3 prev_rpy_;

  // geometry_msgs::Point cur_pos_;
  // geometry_msgs::Point prev_pos_;

  const float kp_ = 0.1;
  const float ki_ = 0;
  const float kd_ = 0.05;
  const float target_vel_ = 3;
  const float accel_ = 3;
  const float timestep = 0.5;

  float cur_integral_;
  float prev_integral_;
  float prev_error_;
  int T_;

  void stateCb(const nav_msgs::Odometry &msg);
  void activeCb();
  void feedbackCb(const ilqr_loco::TrajExecFeedbackConstPtr& feedback);
  void doneCb(const actionlib::SimpleClientGoalState& state,
              const ilqr_loco::TrajExecResultConstPtr& result);

  ilqr_loco::TrajExecGoal GenerateTrajectory(nav_msgs::Odometry prev_state_,
                                             nav_msgs::Odometry cur_state_);
  void SendTrajectory(ilqr_loco::TrajExecGoal &goal);

};

#endif
