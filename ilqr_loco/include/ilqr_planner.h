#ifndef _ILQR_PLANNER_H_
#define _ILQR_PLANNER_H_

#include <vector>

#include <ilqr_loco/TrajExecAction.h>
#include <tf/transform_datatypes.h>
#include <nav_msgs/Odometry.h>
#include <std_msgs/Float32MultiArray.h>

extern "C"{
  #include "iLQG_mpc.h"
}

//iLQR_mpc(double x_cur[10], double x_des[6], double obs[2], int T);

// Note that the inputs to this function can be whatever is convenient for client
ilqr_loco::TrajExecGoal iLQR_gen_traj(nav_msgs::Odometry x_cur, std::vector<double> x_des,
                                        std_msgs::Float32MultiArray obstacle_pos, int T)
{
  //Pre-process inputs - put them in format that C-code wants
  double theta = tf::getYaw(x_cur.pose.pose.orientation);

  double x0[10] = {x_cur.pose.pose.position.x, x_cur.pose.pose.position.y, theta,
                          x_cur.twist.twist.linear.x, x_cur.twist.twist.linear.y,
                          x_cur.twist.twist.angular.z,
                          x_cur.twist.twist.linear.x, 0, 0, 0};

  double* xDes = &x_des[0];
  double Obs[2] = {(double)obstacle_pos.data[0],(double)obstacle_pos.data[1]};

  int N = T+1;
  int n = 10;
  int m = 2;

  //Run iLQR trajectory generation
  struct trajectory Traj;
  Traj.x = (double *) malloc(n*N*sizeof(double));
  Traj.u = (double *) malloc(m*(N-1)*sizeof(double));
  // traj[0]: states, traj[1]: controls
  plan_trajectory(x0,xDes,Obs,T,&Traj);

  //Post-process iLQR trajectory - put states and controls into format that action client wants.
  ilqr_loco::TrajExecGoal goal;	
  for(int i=0; i<N; i++) {
   	nav_msgs::Odometry odom;

    odom.pose.pose.position.x = Traj.x[i*n+0]; // x
    odom.pose.pose.position.y = Traj.x[i*n+1]; // y
    odom.pose.pose.position.z = 0.0;

    geometry_msgs::Quaternion odom_quat = tf::createQuaternionMsgFromYaw(Traj.x[i*n+2]); // phi
    odom.twist.twist.linear.x = Traj.x[i*n+3]; // Ux
 	odom.twist.twist.linear.x = Traj.x[i*n+4]; // Uy
 	odom.twist.twist.angular.z = Traj.x[i*n+5]; // r
	goal.traj.states.push_back(odom);
  }

  for(int i=0; i<N-1; i++) {
  	geometry_msgs::Twist twist;
  	twist.linear.x = Traj.u[i*m+0];
  	twist.angular.z = Traj.u[i*m+1];
  	goal.traj.commands.push_back(twist);
  }

  return goal;
}

#endif
