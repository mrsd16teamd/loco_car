#ifndef _TRAJ_CLIENT_H_
#define _TRAJ_CLIENT_H_

#include <vector>
#include <math.h>

#include "ilqr_planner.h"
#include <ilqr_loco/TrajExecAction.h>

#include <ros/ros.h>

#include <tf/transform_listener.h>
#include <tf/transform_datatypes.h>
#include <actionlib/client/simple_action_client.h>
#include <actionlib/client/terminal_state.h>

#include <geometry_msgs/Twist.h>
#include <std_msgs/Float32MultiArray.h>
#include <geometry_msgs/Point.h>
#include <std_msgs/Char.h>

#define PI 3.1415926535

class TrajClient
{
public:
  bool switch_flag_;
  bool ramp_goal_flag_;

  TrajClient();

protected:
  ros::NodeHandle nh_;
  ros::Subscriber state_sub_;
  ros::Subscriber obs_sub_;
  ros::Subscriber mode_sub_;
  actionlib::SimpleActionClient<ilqr_loco::TrajExecAction> ac_;

  nav_msgs::Odometry start_state_;
  nav_msgs::Odometry cur_state_;
  nav_msgs::Odometry prev_state_;
  geometry_msgs::Point obs_pos_;
  std::vector<double> desired_state_; // Not used anywhere?

  //Constants for rampup planner
  static const float kp_ = 0.45;
  static const float ki_ = 0.05;
  static const float kd_ = 0.1;
  static const float target_vel_ = 3;
  static const float accel_ = 3;
  static const float timestep_ = 0.02;
  static const float timeout_ = 2.5;

  ros::Time start_time_;
  double cur_integral_;
  double prev_error_;
  double cur_vel_;
  int T_;
  bool state_estimate_received_;
  int mode_;

  void rampPlan();
  ilqr_loco::TrajExecGoal rampGenerateTrajectory(nav_msgs::Odometry prev_state_,
                                                 nav_msgs::Odometry cur_state_);

  void ilqgPlan();
  ilqr_loco::TrajExecGoal ilqgGenerateTrajectory(nav_msgs::Odometry cur_state);

  void SendTrajectory(ilqr_loco::TrajExecGoal &goal);

  void activeCb();
  void feedbackCb(const ilqr_loco::TrajExecFeedbackConstPtr& feedback);
  void doneCb(const actionlib::SimpleClientGoalState& state,
              const ilqr_loco::TrajExecResultConstPtr& result);

  void stateCb(const nav_msgs::Odometry &msg);
  void obsCb(const geometry_msgs::PointStamped &msg);
  void modeCb(const geometry_msgs::Point &msg);
  void RampAndiLQR();
};

#endif
