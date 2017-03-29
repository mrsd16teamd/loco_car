#ifndef _ILQR_PLANNER_H_
#define _ILQR_PLANNER_H_

#include <vector>

#include <ilqr_loco/TrajExecAction.h>

// #include "iLQR_mpc.c" //TODO integrate generated c-code
//iLQR_mpc(double x_cur[10], double x_des[6], double obs[2], int T);

// Note that the inputs to this function can be whatever is convenient for client
ilqr_loco::TrajExecGoal iLQR_gen_traj(nav_msgs::Odometry x_cur, std::vector<double> x_des,
                                        std_msgs::Float32MultiArray obstacle_pos, int T)
{
  // TODO pre-process inputs as necessary, put into C-style arrays

  // TODO convert orientation  from quaternion
  double theta = ___;

  double x_current[10] = {x_cur.pose.pose.position.x, x_cur.pose.pose.position.y, theta,
                          x_cur.twist.twist.linear.x, x_cur.twist.twist.linear.y,
                          x_cur.twist.twist.angular.z,
                          x_cur.twist.twist.linear.x, 0, 0, 0}

  // TODO what is the output type?
  sometype outputs = iLQR_mpc(x_current, x_desired, obs, T);

  ilqr_loco::TrajExecGoal goal;
  // TODO put output from c-code into action message

  return goal;
}

#endif
