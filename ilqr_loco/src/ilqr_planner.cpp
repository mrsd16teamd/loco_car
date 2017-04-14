#ifndef _ILQR_PLANNER_H_
#define _ILQR_PLANNER_H_

#include <vector>
#include <ctime>

#include <ilqr_loco/TrajExecAction.h>
#include <tf/transform_datatypes.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Point.h>

#include "traj_client.h"

extern "C"{
  #include "iLQG.h"
  #include "iLQG_plan.h"
}

// Note that the inputs to this function can be whatever is convenient for client
void TrajClient::iLQR_gen_traj(nav_msgs::Odometry &x_cur, std::vector<double> &u_init, std::vector<double> &x_des,
                               geometry_msgs::Point &obstacle_pos, int T, tOptSet *o, ilqr_loco::TrajExecGoal &goal)
{
  //Pre-process inputs - put them in format that C-code wants
  double theta = tf::getYaw(x_cur.pose.pose.orientation);

  double x0[10] = {x_cur.pose.pose.position.x, x_cur.pose.pose.position.y, theta,
                   x_cur.twist.twist.linear.x, x_cur.twist.twist.linear.y,
                   x_cur.twist.twist.angular.z,
                   x_cur.twist.twist.linear.x, 0, 0, 0};

  double* xDes = &x_des[0]; //std::vector trick to convert vector to C-style array
  double* u0 = &u_init[0];
  double Obs[2] = {(double)obstacle_pos.x,(double)obstacle_pos.y};

  int N = T+1;
  int n = 10; //state size
  int m = 2;  //control size

  //Run iLQR trajectory generation
  struct trajectory Traj;
  Traj.x = (double *) malloc(n*N*sizeof(double));
  Traj.u = (double *) malloc(m*(N-1)*sizeof(double));

  // traj[0]: states, traj[1]: controls
  plan_trajectory(x0,u0,xDes,Obs,T,o,&Traj);

  // TODO find better way that doesn't copy twice
  std::vector<double> u_sol(Traj.u, Traj.u+N);
  u_init = u_sol;

  //TODO bring this back!
  //Put states and controls into format that action client wants.
  // goal.traj.states.reserve(N);
  // goal.traj.commands.reserve(N);

  for(int i=0; i<N; i++) {
   	nav_msgs::Odometry odom;
    FillOdomMsg(odom, Traj.x[i*n+0], Traj.x[i*n+1], Traj.x[i*n+2],
                      Traj.x[i*n+3], Traj.x[i*n+4], Traj.x[i*n+5]);
  	goal.traj.states.push_back(odom);
  }

  for(int i=0; i<N-1; i++) {
  	geometry_msgs::Twist twist;
    FillTwistMsg(twist, double(Traj.u[i*m+0]), double(Traj.u[i*m+1]));
  	goal.traj.commands.push_back(twist);
  }

  // append zero command to stop vehicle
  geometry_msgs::Twist twist;
  FillTwistMsg(twist, 0.0, 0.0);
  goal.traj.commands.push_back(twist);
}

#endif
