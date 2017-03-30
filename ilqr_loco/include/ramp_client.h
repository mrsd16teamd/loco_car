#ifndef _RAMP_CLIENT_H_
#define _RAMP_CLIENT_H_

#include <ros/ros.h>
#include <tf/transform_listener.h>
#include <tf/transform_datatypes.h>
#include <actionlib/client/simple_action_client.h>
#include <actionlib/client/terminal_state.h>
#include <ilqr_loco/TrajExecAction.h>


// Do we need these?
#include <string>
#include <iostream>

class RampPlanner
{
public:
  void Plan();
  bool flag_;
  bool init_flag_;

  RampPlanner() : ac_("traj_executer", true) {
    state_sub_  = nh_.subscribe("odometry/filtered", 1, &RampPlanner::stateCb, this);

    ROS_INFO("Waiting for action server to start.");
    ac_.waitForServer(); //will wait for infinite time
    ROS_INFO("Action server started.");

    flag_ = false;
    cur_integral_ = 0.0;
    prev_error_ = 0.0;
    T_ = 0;
    init_flag_ = false;
    start_time_ = ros::Time::now();
  }

protected:
  ros::NodeHandle nh_;
  ros::Subscriber state_sub_;
  actionlib::SimpleActionClient<ilqr_loco::TrajExecAction> ac_;

  nav_msgs::Odometry cur_state_;
  nav_msgs::Odometry prev_state_;

  static const float kp_ = 0.45;
  static const float ki_ = 0.05;
  static const float kd_ = 0.1;
  static const float target_vel_ = 3;
  static const float accel_ = 3;
  static const float timestep_ = 0.02;
  static const float timeout_ = 0.5;
  
  ros::Time start_time_;
  double cur_integral_;
  double prev_error_;
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
