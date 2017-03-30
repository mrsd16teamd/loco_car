#include "ramp_client.h"

void RampPlanner::stateCb(const nav_msgs::Odometry &msg) {
  prev_state_ = cur_state_;
  cur_state_ = msg;
  init_flag_ = true;
}

void RampPlanner::activeCb() {}

void RampPlanner::feedbackCb(const ilqr_loco::TrajExecFeedbackConstPtr& feedback) {}

void RampPlanner::doneCb(const actionlib::SimpleClientGoalState& state,
                        const ilqr_loco::TrajExecResultConstPtr& result) {}

ilqr_loco::TrajExecGoal RampPlanner::GenerateTrajectory(nav_msgs::Odometry prev_state,
                                                        nav_msgs::Odometry cur_state) {

  ilqr_loco::TrajExecGoal goal;
  goal.traj.header.seq = T_;
  goal.traj.header.stamp = ros::Time::now();
  goal.traj.header.frame_id = "/base_link";

  double dt = (cur_state.header.stamp).toSec() - (prev_state.header.stamp).toSec();

  // Get current orientation
  // tf::Quaternion quat(cur_state.pose.pose.orientation.x, cur_state.pose.pose.orientation.y,
  //                     cur_state.pose.pose.orientation.z, cur_state.pose.pose.orientation.w);

  double roll, pitch, yaw;
  // tf::Matrix3x3 m(quat);
  // m.getRPY(roll, pitch, yaw);
  // prev_rpy_ = cur_rpy_;
  // cur_rpy_.x = roll;
  // cur_rpy_.y = pitch;
  // cur_rpy_.z = yaw;

  yaw = tf::getYaw(cur_state.pose.pose.orientation);

  ROS_INFO("yaw = %f",yaw);

  // PID control for vehicle heading
  float error = 0 - yaw;
  prev_integral_ = cur_integral_;
  cur_integral_ += error*dt;
  float output = kp_*error + ki_*cur_integral_ + kd_*(error-prev_error_)/dt;
  ROS_INFO("P = %f, I = %f, D = %f",kp_*error, ki_*cur_integral_, kd_*(error-prev_error_)/dt);
  ROS_INFO("PID output = %f",output);
  prev_error_ = error;

  // Generate goal
  float v = cur_state.twist.twist.linear.x + accel_*dt;

  geometry_msgs::Twist control_msg;

  control_msg.linear.x = v<3.0 ? v : 3.0;
  control_msg.angular.z = output;
  goal.traj.commands.push_back(control_msg);
  goal.traj.states.push_back(cur_state);
  ++T_;

  flag_ = v>=3.0 ? true : false;  // Ramp completion flag

  return goal;
}

void RampPlanner::SendTrajectory(ilqr_loco::TrajExecGoal &goal) {
  ac_.sendGoal(goal,
              boost::bind(&RampPlanner::doneCb, this, _1, _2),
              boost::bind(&RampPlanner::activeCb, this),
              boost::bind(&RampPlanner::feedbackCb, this, _1));
}

void RampPlanner::Plan() {
  ROS_INFO("Waiting for initial odometry information...");
  while (!init_flag_)
    ros::spinOnce();

  ROS_INFO("Starting speed ramp!");
  ros::Time start_time = ros::Time::now();
  ros::Duration timeout(3.0);
  while(ros::Time::now() - start_time < timeout) {
    ilqr_loco::TrajExecGoal goal = RampPlanner::GenerateTrajectory(prev_state_, cur_state_);
    RampPlanner::SendTrajectory(goal);
    ros::spinOnce();
  }
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "ramp_client");
  RampPlanner ramp;
  ros::spin();

  return 0;
}
