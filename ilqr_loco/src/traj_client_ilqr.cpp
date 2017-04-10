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

  iLQR_gen_traj(cur_state, init_control_seq_, x_des_, obs_pos_, T_horizon_, goal);
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
  ilqr_loco::TrajExecGoal goal;

  if(ros::Time::now() - start_time_ < ros::Duration(mpc_timeout_))
  {
    ROS_INFO("iLQR timed out.");
    return;;
  }

  bool goal_achieved = false;

  while(!goal_achieved)
  {
    ilqr_loco::TrajExecGoal goal = ilqgGenerateTrajectory(cur_state_);
    SendTrajectory(goal);

    if (DistToGoal() < goal_threshold_){
      ROS_INFO("Reached goal point.");
      ilqr_loco::TrajExecGoal stop;
      geometry_msgs::Twist twist;
      FillTwistMsg(twist, 0.0, 0.0);
      stop.traj.commands.push_back(twist);
      SendTrajectory(stop);
      return;
    }

    ros::spinOnce();
  }
}
