#include "traj_client.h"

double TrajClient::DistToGoal()
{
  ROS_INFO("x_des: %f, %f. x_cur: %f, %f", x_des_[0], x_des_[1], cur_state_.pose.pose.position.x, cur_state_.pose.pose.position.y);
  return sqrt( pow((x_des_[0]- cur_state_.pose.pose.position.x), 2) +
               pow((x_des_[1]- cur_state_.pose.pose.position.y), 2) );
}

// Calls iLQG_plan.c to generate new trajectory
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
  ilqr_loco::TrajExecGoal goal = ilqgGenerateTrajectory(cur_state_);
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
  std::vector<double> replan_times = {0.0, 0.75, 1.5};
  bool plan_next_ = true;

  ilqr_loco::TrajExecGoal goal;

  int iter_count = 0;

  while (ros::ok())
  {
    if(plan_next_)
    {
      ROS_INFO("Replan #%d", iter_count);
      ilqrPlan();
      plan_next_ = false;
      if (++iter_count >= replan_times.size())
      {
        ROS_INFO("Done planning.");
      }
    }
    if(ros::Time::now() - start_time_ > ros::Duration(replan_times[iter_count]) && iter_count < replan_times.size())
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
