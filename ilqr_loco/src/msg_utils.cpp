//
// MIT License
//
// Copyright (c) 2017 MRSD Team D - LoCo
// The Robotics Institute, Carnegie Mellon University
// http://mrsdprojects.ri.cmu.edu/2016teamd/
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include "traj_client.h"

void TrajClient::FillGoalMsgHeader(ilqr_loco::TrajExecGoal &goal)
{
  goal.traj.header.seq = T_;
  goal.traj.header.stamp = ros::Time::now();
  goal.traj.header.frame_id = "base_link";
  goal.traj.timestep = timestep_;
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

void TrajClient::SendInitControlSeq()
{
  ilqr_loco::TrajExecGoal goal;
  FillGoalMsgHeader(goal);
  geometry_msgs::Twist control_msg;

  for (int i=0; i<(init_control_seq_.size()/2); i++)
  {
    FillTwistMsg(control_msg, init_control_seq_[2*i], init_control_seq_[(2*i)+1]);
    goal.traj.commands.push_back(control_msg);
    goal.traj.states.push_back(cur_state_);
  }
  FillTwistMsg(control_msg, 0, 0);
  goal.traj.commands.push_back(control_msg);
  goal.traj.states.push_back(cur_state_);

  SendTrajectory(goal);
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
