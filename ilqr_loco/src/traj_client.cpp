#include "traj_client.h"

void TrajClient::activeCb() {}

void TrajClient::feedbackCb(const ilqr_loco::TrajExecFeedbackConstPtr& feedback) { }

void TrajClient::doneCb(const actionlib::SimpleClientGoalState& state,
  const ilqr_loco::TrajExecResultConstPtr& result) {}

void TrajClient::stateCb(const nav_msgs::Odometry &msg)
{
  most_recent_state = msg;
}
void TrajClient::obsCb(const std_msgs::Float32MultiArray &msg)
{
  obs_pos = msg;
}

// Calls iLQR_mpc.c to generate new trajectory
ilqr_loco::TrajExecGoal TrajClient::GenerateTrajectory()
{
  ROS_INFO("Generating trajectory.");
  ilqr_loco::TrajExecGoal goal;

  //TODO use sensor feedback here
  //TODO replace this with actual trajectory planners; see below for possible implementation
  // goal = ramp_planner.GenerateTrajectory();
  // goal = iLQR_gen_traj(x_current, x_desired, obstacle_pos, T);
  //    see ilqr_planner.h

  ros::Time begin = ros::Time::now();
  
  goal.traj.timestep = 0.05;


  // For now, just hard-coded commands
  int traj_length = 20;
  for (int i=0; i<traj_length; i++)
  {
    geometry_msgs::Twist msg;
    msg.linear.x = 0.5;
    msg.angular.z = 0.5;
    goal.traj.commands.push_back(msg);
  }
}

void TrajClient::SendTrajectory(ilqr_loco::TrajExecGoal &goal)
{
  ac.sendGoal(goal,
              boost::bind(&TrajClient::doneCb, this, _1, _2),
              boost::bind(&TrajClient::activeCb, this),
              boost::bind(&TrajClient::feedbackCb, this, _1));
}


// This function gets called when server comes online!
void TrajClient::Plan()
{
  ilqr_loco::TrajExecGoal goal = GenerateTrajectory();
  SendTrajectory(goal);
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "test_client");
  TrajClient client;
  ros::spin();

  return 0;
}
