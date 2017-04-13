#include "traj_client.h"

double TrajClient::DistToGoal()
{
  return sqrt( pow((x_des_[0]- cur_state_.pose.pose.position.x), 2) +
               pow((x_des_[1]- cur_state_.pose.pose.position.x), 2) );
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

  iLQR_gen_traj(cur_state, init_control_seq_, x_des_, obs_pos_, T_horizon_, &Opt, goal);
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

  int iter_count = 0;
  bool goal_achieved = false;

  while(!goal_achieved)
  {
    if(ros::Time::now() - start_time_ < ros::Duration(mpc_timeout_))
    {
      ROS_INFO("iLQR timed out.");
      SendZeroCommand();
      return;
    }

    ROS_INFO("Receding horizon iteration #%d", iter_count);

    goal = ilqgGenerateTrajectory(cur_state_);

    // TODO do some quick checks on trajectory

    SendTrajectory(goal);

    ROS_INFO("DistToGoal: %f", DistToGoal());
    if (DistToGoal() < goal_threshold_) {
      ROS_INFO("Reached goal point.");
			SendZeroCommand();
      return;
    }

    ros::spinOnce(); // to pick up new state estimates
    iter_count++;
  }
}
