#include "traj_client.h"

void TrajClient::FillGoalMsgHeader(ilqr_loco::TrajExecGoal &goal)
{
  goal.traj.header.seq = T_;
  goal.traj.header.stamp = ros::Time::now();
  goal.traj.header.frame_id = "base_link";
}

void TrajClient::FillTwistMsg(geometry_msgs::Twist &twist, double lin_x, double ang_z)
{
  twist.linear.x = lin_x;
  twist.angular.z = ang_z;
}

void TrajClient::FillOdomMsg(nav_msgs::Odometry &odom, double x, double y,
                             double yaw, double Ux, double Uy, double w)
{
  odom.pose.pose.position.x = x;
  odom.pose.pose.position.y = y;
  odom.pose.pose.position.z = 0.0;

  geometry_msgs::Quaternion odom_quat = tf::createQuaternionMsgFromYaw(yaw);
  odom.pose.pose.orientation = odom_quat;

  odom.twist.twist.linear.x = Ux;
  odom.twist.twist.linear.y = Uy;
  odom.twist.twist.angular.z = w;
}

void TrajClient::SendZeroCommand()
{
  mode_ = 0;

  ilqr_loco::TrajExecGoal end_goal;

  FillGoalMsgHeader(end_goal);
  geometry_msgs::Twist control_msg;
  FillTwistMsg(control_msg, 0, 0);

  end_goal.traj.commands.push_back(control_msg);
  end_goal.traj.states.push_back(cur_state_);
  SendTrajectory(end_goal);
  ROS_INFO("Sent zero command.");
}

void TrajClient::SendTrajectory(ilqr_loco::TrajExecGoal &goal)
{
  // ROS_INFO("Sending trajectory.");
  // ac_.sendGoal(goal);

  ac_.sendGoal(goal,
               NULL,
               NULL,
              //  boost::bind(&TrajClient::doneCb, this, _1, _2),
              //  boost::bind(&TrajClient::activeCb, this),
               boost::bind(&TrajClient::feedbackCb, this, _1));
}
