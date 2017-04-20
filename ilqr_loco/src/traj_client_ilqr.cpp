#include "traj_client.h"

ilqr_loco::TrajExecGoal TrajClient::GenTrajILQR(nav_msgs::Odometry &x_cur, std::vector<double> &u_init,
                                  std::vector<double> &x_des, geometry_msgs::Point &obstacle_pos)
{
  ROS_INFO("Generating iLQG trajectory.");
  ilqr_loco::TrajExecGoal goal;
  FillGoalMsgHeader(goal);
  goal.traj.timestep = timestep_;

  // ROS_INFO("Start state (before prediction): %f, %f, %f, %f, %f, %f",
  //         cur_state_.pose.pose.position.x, cur_state_.pose.pose.position.y, theta,
  //         cur_state_.twist.twist.linear.x, cur_state_.twist.twist.linear.y, cur_state_.twist.twist.angular.z);

  //Pre-process inputs - put them in format that C-code wants
  double theta = tf::getYaw(x_cur.pose.pose.orientation);

  // TODO figure out good way to initialize previous steering
  double x0[10] = {x_cur.pose.pose.position.x, x_cur.pose.pose.position.y, theta,
                   x_cur.twist.twist.linear.x, x_cur.twist.twist.linear.y,
                   x_cur.twist.twist.angular.z,
                   x_cur.twist.twist.linear.x, 0.1, 0, 0};

  double* xDes = &x_des[0]; //std::vector trick to convert vector to C-style array
  double* u0 = &u_init[0];
  double Obs[2] = {(double)obstacle_pos.x, (double)obstacle_pos.y};

  int N = T_horizon_+1;
  int n = 10; //state size
  int m = 2;  //control size

  //Run iLQR trajectory generation
  struct trajectory Traj;
  Traj.x = (double *) malloc(n*N*sizeof(double));
  Traj.u = (double *) malloc(m*(N-1)*sizeof(double));

  // traj[0]: states, traj[1]: controls
  plan_trajectory(x0, u0, xDes, Obs, T_horizon_, &Opt, &Traj);

  // TODO find better way that doesn't copy twice
  std::vector<double> u_sol(Traj.u, Traj.u+(2*T_horizon_));
  u_init = u_sol;

  // TODO save trajectory here
  // std::vector<double> x_sol(Traj.x, Traj.x+(n*N));
  // x_traj_saved_ = x_sol;

  //Put states and controls into format that action client wants.
  goal.traj.states.reserve(N);
  goal.traj.commands.reserve(N);

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

  return goal;
}

void TrajClient::PlanFromCurrentStateILQR()
{
  ilqr_loco::TrajExecGoal goal = GenTrajILQR(cur_state_, u_seq_saved_, x_des_, obs_pos_);

  if (mode_==7) // turn on pid heading corrections during server execution
  	goal.traj.mode = 1;
  SendTrajectory(goal);
}

void TrajClient::MpcILQR()
{
  start_time_ = ros::Time::now();

  ilqr_loco::TrajExecGoal goal;

  while( (DistToGoal() > goal_threshold_) && (ros::Time::now() - start_time_ > ros::Duration(mpc_timeout_)) )
  {
    ROS_INFO("Receding horizon iteration #%d", T_);

    ilqr_loco::TrajExecGoal goal = GenTrajILQR(cur_state_, u_seq_saved_, x_des_, obs_pos_);
    // TODO do some quick checks on trajectory?

    SendTrajectory(goal);

    ros::spinOnce(); // to pick up new state estimates
    T_++;
    ROS_INFO("DistToGoal: %f", DistToGoal());
  }
  SendZeroCommand();
}

void TrajClient::SparseReplanILQR()
{
  start_time_ = ros::Time::now();

  bool plan_next_ = true;
  while (ros::ok())
  {
    if(plan_next_)
    {
      ROS_INFO("Replan #%d", T_);
      PlanFromCurrentStateILQR();
      plan_next_ = false;
      if (++T_ >= replan_times_.size())
        ROS_INFO("Done planning.");
    }

    if(ros::Time::now() - start_time_ > ros::Duration(replan_times_[T_]) && T_ < replan_times_.size())
      plan_next_ = true;

    if(ros::Time::now() - start_time_ > ros::Duration(mpc_timeout_))
    {
      ROS_INFO("iLQR sparse replan timed out.");
      SendZeroCommand();
      break;
    }
    ros::spinOnce(); // to pick up new state estimates
  }
}

nav_msgs::Odometry TrajClient::ExtrapolateState(const nav_msgs::Odometry &state)
{
  nav_msgs::Odometry extrapolated = state;

  double dt = 0.1; //[s] TODO make this parameter, or function of T_horizon_ and max_iter_

  extrapolated.pose.pose.position.x += (dt*extrapolated.twist.twist.linear.x);
  extrapolated.pose.pose.position.y += (dt*extrapolated.twist.twist.linear.y);

  double theta = tf::getYaw(extrapolated.pose.pose.orientation);
  theta += dt*extrapolated.twist.twist.angular.z;
  extrapolated.pose.pose.orientation = tf::createQuaternionMsgFromYaw(theta);

  return extrapolated;
}

double TrajClient::DistToGoal()
{
  return sqrt( pow((x_des_[0]- cur_state_.pose.pose.position.x), 2) +
               pow((x_des_[1]- cur_state_.pose.pose.position.y), 2) );
}
