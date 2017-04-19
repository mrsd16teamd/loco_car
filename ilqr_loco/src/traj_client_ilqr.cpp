#include "traj_client.h"

double TrajClient::DistToGoal()
{
  ROS_INFO("x_des: %f, %f. x_cur: %f, %f", x_des_[0], x_des_[1], cur_state_.pose.pose.position.x, cur_state_.pose.pose.position.y);
  return sqrt( pow((x_des_[0]- cur_state_.pose.pose.position.x), 2) +
               pow((x_des_[1]- cur_state_.pose.pose.position.y), 2) );
}

// change this to take x_cur in vector form
// need a helper to translate odometry to vector, and revisit all the callers <--- good amount of work
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
  std::copy(Traj.u, Traj.u+m*(N-1), u_init.begin());
  ROS_INFO("T: %d, size of u_init: %d", T, int(u_init.size()));

  x_seq_saved_.resize(n*N,0);
  std::copy(Traj.x, Traj.x+n*N, x_seq_saved_.begin());
  ROS_INFO("T: %d, size of x_seq_saved: %d", T, int(x_seq_saved_.size()));

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

// merge this with ilqr_gen_traj
ilqr_loco::TrajExecGoal TrajClient::ilqgGenerateTrajectory(nav_msgs::Odometry cur_state)
{
  ROS_INFO("Generating iLQG trajectory.");
  ilqr_loco::TrajExecGoal goal;
  FillGoalMsgHeader(goal);
  goal.traj.timestep = timestep_;

  double theta = tf::getYaw(cur_state.pose.pose.orientation);
  ROS_INFO("Start state: %f, %f, %f, %f, %f, %f",
            cur_state.pose.pose.position.x, cur_state.pose.pose.position.y, theta,
            cur_state.twist.twist.linear.x, cur_state.twist.twist.linear.y, cur_state.twist.twist.angular.z);
  ROS_INFO("Obs pos: %f, %f", obs_pos_.x, obs_pos_.y);

  iLQR_gen_traj(cur_state, u_seq_saved_, x_des_, obs_pos_, T_horizon_, &Opt, goal);
  ++T_;

  return goal;
}

void TrajClient::ilqrPlan()
{
  if (mode_==2) {
    double theta = tf::getYaw(cur_state_.pose.pose.orientation);
    ROS_INFO("Start state (before prediction): %f, %f, %f, %f, %f, %f",
            cur_state_.pose.pose.position.x, cur_state_.pose.pose.position.y, theta,
            cur_state_.twist.twist.linear.x, cur_state_.twist.twist.linear.y, cur_state_.twist.twist.angular.z);

    cur_state_.pose.pose.position.x += (execution_delay_*cur_state_.twist.twist.linear.x);
    // this is extra work, maybe use if statements in iLQR_gen_traj would be better

    theta += execution_delay_*cur_state_.twist.twist.angular.z;
    cur_state_.pose.pose.orientation = tf::createQuaternionMsgFromYaw(theta);
  }

  ilqr_loco::TrajExecGoal goal = ilqgGenerateTrajectory(cur_state_);
  if (mode_==7) // turn on pid heading corrections during server execution
  	goal.traj.mode = 1;
  SendTrajectory(goal);
}

void TrajClient::ilqrMPC()
{
  start_time_ = ros::Time::now();

  ilqr_loco::TrajExecGoal goal;

  //HACK HERE
  FillGoalMsgHeader(goal);
goal.traj.timestep = timestep_;
  for(int i=0; i<u_seq_saved_.size()/2; i++) {
	geometry_msgs::Twist twist;
	FillTwistMsg(twist, u_seq_saved_[2*i], u_seq_saved_[2*i+1]);
	goal.traj.commands.push_back(twist);
  }
for(int i=0; i<u_seq_saved_.size()+1; i++) {
  goal.traj.states.push_back(cur_state_);
}
  SendTrajectory(goal);
  ROS_INFO("Sent first swerve.");
////
  int iter_count = 0;
  bool goal_achieved = false;

  while(!goal_achieved)
  {
    if(ros::Time::now() - start_time_ > ros::Duration(mpc_timeout_))
    {
      ROS_INFO("iLQR timed out.");
      SendZeroCommand();
      break;
    }

    ROS_INFO("Receding horizon iteration #%d", iter_count);
    ROS_INFO("u0[0]: %f", u_seq_saved_[0]);

    goal = ilqgGenerateTrajectory(cur_state_);

    // TODO do some quick checks on trajectory

    SendTrajectory(goal);

    ROS_INFO("DistToGoal: %f", DistToGoal());
    if (DistToGoal() < goal_threshold_) {
    	ROS_INFO("Reached goal point.");
	  	SendZeroCommand();
      break;
    }

    ros::spinOnce(); // to pick up new state estimates
    iter_count++;
  }
}

void TrajClient::ilqrSparseReplan()
{
  start_time_ = ros::Time::now();
  std::vector<double> replan_times;
  if (target_vel_ <= 1)
    replan_times = {0.0, 0.5, 1.0, 1.5}; // for static: {0.0, 0.75, 1.5};
  else if (target_vel_ <= 2)
    replan_times = {0.0, 0.75, 1.5}; // for static: {0.0, 0.75, 1.5};
  else if (target_vel_ <= 4)
    replan_times = {0.0, 0.75, 1.5}; // for static: {0.0, 0.75, 1.5};

  bool plan_next_ = true;

  ilqr_loco::TrajExecGoal goal;

  while (ros::ok())
  {
    if(plan_next_)
    {
      ROS_INFO("Replan #%d", T_);
      
      if(T_>0) {	
      	double theta = tf::getYaw(cur_state_.pose.pose.orientation);

      	ROS_INFO("Start state (before prediction): %f, %f, %f",
        cur_state_.pose.pose.position.x, cur_state_.pose.pose.position.y, theta);

      	int step_offset = int(execution_delay_/timestep_);
	  	int n = 10; //state size
      	std::vector<double> state_offset(3,0);
      	for(int i=0; i<3; i++) {
      		state_offset[i] = x_seq_saved_[step_offset*n+i] - x_seq_saved_[i];
      	}

      	ROS_INFO("Predicted offset: %f, %f, %f",
        state_offset[0], state_offset[1], state_offset[2]);
      	
      	cur_state_.pose.pose.position.x += state_offset[0];
      	cur_state_.pose.pose.position.y += state_offset[1];

      	
      	theta += state_offset[2];
    	cur_state_.pose.pose.orientation = tf::createQuaternionMsgFromYaw(theta);
      }

      ilqrPlan();
      plan_next_ = false;
      if (T_ >= replan_times.size())
      {
        ROS_INFO("Done planning.");
      }
    }
    if(ros::Time::now() - start_time_ > ros::Duration(replan_times[T_]) && T_ < replan_times.size())
    {
      plan_next_ = true;
      ROS_INFO("Next plan.");
    }
    if(ros::Time::now() - start_time_ > ros::Duration(mpc_timeout_))
    {
      ROS_INFO("iLQR timed out.");
      SendZeroCommand();
      break;
    }
    ros::spinOnce(); // to pick up new state estimates
  }

}
