#ifndef _TRAJ_CLIENT_H_
#define _TRAJ_CLIENT_H_

#include <vector>
#include <math.h>

#include <ilqr_loco/TrajExecAction.h>

#include <ros/ros.h>
#include <tf/transform_listener.h>
#include <tf/transform_datatypes.h>
#include <actionlib/client/simple_action_client.h>
#include <actionlib/client/terminal_state.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/Point.h>

extern "C"{
  #include "iLQG.h"
  #include "iLQG_plan.h"
}

#define PI 3.1415926535

class TrajClient
{
public:
  bool ramp_goal_flag_;

  TrajClient();

protected:
  // ROS Handles
  ros::NodeHandle nh_;
  ros::Subscriber state_sub_;
  ros::Subscriber obs_sub_;
  ros::Subscriber mode_sub_;
  actionlib::SimpleActionClient<ilqr_loco::TrajExecAction> ac_;

  // ROS Parameters
  int T_horizon_;
  std::vector<double> init_control_seq_;
  std::vector<double> x_des_;
  tOptSet Opt;
  double dt;

  // iLQR Opt.p: Car Params
  double g_;
  double L_;
  double m_;
  double b_;
  double a_;
  double G_f_;
  double G_r_;
  double c_x_;
  double c_a_;
  double Iz_;
  double mu_;
  double mu_s_;
  std::vector<double> limThr_;
  std::vector<double> limSteer_;

  // iLQR Opt.p: Costs
  std::vector<double> cu_;        // Control cost
  std::vector<double> cdu_;       // Change of control cost
  std::vector<double> cf_;        // Final cost
  std::vector<double> pf_;        // Final cost smooth terms
  std::vector<double> cx_;        // Running cost of position
  std::vector<double> cdx_;       // Running cost of velocity
  std::vector<double> px_;        // Running cost smooth terms
  double cdrift_;                 // Drift cost
  double k_pos_;                  // Obstacle pos cost
  double k_vel_;                  // Obstacle vel cost
  double d_thres_;                // Obstacle threshold

  // Helper variables
  int T_;                         // Sequence ID number (starts from 0, in lifetime of client)
  int mode_;                      // Operation mode from keyboard teleop
  ros::Time start_time_;          // Operation start time
  bool state_estimate_received_;  // Initial estimate flag
  bool obs_received_;

  // State variables
  nav_msgs::Odometry start_state_;
  nav_msgs::Odometry cur_state_;
  nav_msgs::Odometry prev_state_;
  geometry_msgs::Point obs_pos_;

  // Ramp up
  double cur_integral_;
  double prev_error_;
  double cur_vel_;

  //Constants for rampup planner
  float kp_, ki_, kd_;
  float accel_;
  float target_vel_;
  float timeout_;
  double timestep_;

  //iLQR parameters
  float mpc_timeout_;
  float goal_threshold_;
  double ilqr_tolFun_;
  double ilqr_tolConstraint_;
  double ilqr_tolGrad_;
  int ilqr_max_iter_;
  int ilqr_regType_;
  int ilqr_debug_level_;

  void LoadParams();
	void LoadCarParams();
  void LoadCostParams();
	void SetOptParams(tOptSet *o);
  void LoadOpt();

  void rampPlan();
  ilqr_loco::TrajExecGoal rampGenerateTrajectory(nav_msgs::Odometry prev_state_,
                                                 nav_msgs::Odometry cur_state_);

  void ilqrPlan();
  ilqr_loco::TrajExecGoal ilqgGenerateTrajectory(nav_msgs::Odometry cur_state);
  void iLQR_gen_traj(nav_msgs::Odometry &x_cur, std::vector<double> &u_init, std::vector<double> &x_des,
                     geometry_msgs::Point &obstacle_pos, int T, tOptSet* o, ilqr_loco::TrajExecGoal &goal);
  void ilqrMPC();
  double DistToGoal();

  void RampAndiLQR();
	void SendZeroCommand();
  void SendTrajectory(ilqr_loco::TrajExecGoal &goal);

  void stateCb(const nav_msgs::Odometry &msg);
  void obsCb(const geometry_msgs::PointStamped &msg);
  void modeCb(const geometry_msgs::Point &msg);

  void FillGoalMsgHeader(ilqr_loco::TrajExecGoal &goal);
  void FillTwistMsg(geometry_msgs::Twist &twist, double lin_x, double ang_z);
  void FillOdomMsg(nav_msgs::Odometry &odom, double x, double y,
                   double yaw, double Ux, double Uy, double w);

   void activeCb();
   void feedbackCb(const ilqr_loco::TrajExecFeedbackConstPtr& feedback);
   void doneCb(const actionlib::SimpleClientGoalState& state,
               const ilqr_loco::TrajExecResultConstPtr& result);
};

#endif
