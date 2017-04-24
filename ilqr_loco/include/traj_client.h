#ifndef _TRAJ_CLIENT_H_
#define _TRAJ_CLIENT_H_

#include <vector>
#include <math.h>

#include <ilqr_loco/TrajExecAction.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Point.h>
#include <sensor_msgs/LaserScan.h>
#include "try_get_param.h"

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
  ros::NodeHandle nh;
  ros::Subscriber state_sub_;
  ros::Subscriber obs_sub_;
  ros::Subscriber mode_sub_;
  ros::Subscriber scan_sub_;
  ros::Publisher predicted_state_pub_;
  ros::Publisher obs_pos_pub_;
  actionlib::SimpleActionClient<ilqr_loco::TrajExecAction> ac_;
  tf::TransformListener tf_listener_;

  // ilqr parameters and saved data
  int T_horizon_;
  std::vector<double> init_control_seq_;
  std::vector<double> u_seq_saved_;
  std::vector<double> x_traj_saved_;
  std::vector<double> x_des_;
  tOptSet Opt;

  // iLQR Opt.p: Car Params
  double g_, L_, m_, b_, a_, G_f_, G_r_, c_x_, c_a_, Iz_, mu_, mu_s_;
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
  bool found_obstacle_;

  // State variables
  nav_msgs::Odometry start_state_;
  nav_msgs::Odometry cur_state_;
  nav_msgs::Odometry prev_state_;
  geometry_msgs::Point obs_pos_;

  // Ramp up
  double cur_integral_;
  double prev_error_;
  double cur_vel_;
  double ramp_start_y_;

  //Constants for rampup planner
  float kp_, ki_, kd_, kp_y_;
  float accel_;
  float target_vel_;
  float pre_ramp_vel_;
  float pre_ramp_time_;
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
  std::vector<double> replan_times_;
  double replan_rate_;
  int step_on_last_traj_;
  int use_extrapolate_;
  double extrapolate_dt_;
  int use_pid_;

  // Obstacle detector parameters
  float obs_dist_thres_; //[m]f
  float obs_percent_thres_;
  float scan_front_angle_;
  float scan_min_index_, scan_max_index_;

  void LoadParams();
  void LoadCarParams();
  void LoadCostParams();
  void SetOptParams(tOptSet *o);
  void LoadOpt();

  void rampPlan();
  ilqr_loco::TrajExecGoal rampGenerateTrajectory(nav_msgs::Odometry prev_state_,
                                                 nav_msgs::Odometry cur_state_);

  void Plan();
  void PlanFromCurrentStateILQR();
  void PlanFromExtrapolatedILQR();
  ilqr_loco::TrajExecGoal GenTrajILQR(nav_msgs::Odometry &x_cur, std::vector<double> &u_init,
          std::vector<double> &x_des, geometry_msgs::Point &obstacle_pos);
  void MpcILQR();
  void FixedRateReplanILQR();
  double DistToGoal();
  nav_msgs::Odometry ExtrapolateState(const nav_msgs::Odometry &state);

  void SendZeroCommand();
  void SendSwerveCommand();
  void SendTrajectory(ilqr_loco::TrajExecGoal &goal);
  void SendInitControlSeq();

  void stateCb(const nav_msgs::Odometry &msg);
  void modeCb(const geometry_msgs::Point &msg);
  void scanCb(const sensor_msgs::LaserScanConstPtr &msg);

  //  void activeCb();
  void feedbackCb(const ilqr_loco::TrajExecFeedbackConstPtr& feedback);
  //  void doneCb(const actionlib::SimpleClientGoalState& state,
              //  const ilqr_loco::TrajExecResultConstPtr& result);

  void FillGoalMsgHeader(ilqr_loco::TrajExecGoal &goal);
  void FillTwistMsg(geometry_msgs::Twist &twist, double lin_x, double ang_z);
  void FillOdomMsg(nav_msgs::Odometry &odom, double x, double y,
                   double yaw, double Ux, double Uy, double w);

  void InitDetector();
  bool TransformLaserToMap(geometry_msgs::PointStamped &pos_laser_frame, geometry_msgs::PointStamped &pos_map_frame);
  void InsertFakeObs();
  void ReactToObstacle();

  double clamp(double val, double min_val, double max_val)
  {
    return std::max(min_val, std::min(val, max_val));
  }
};

#endif
