#include "ramp_client.h"

void RampPlanner::stateCb(const nav_msgs::Odometry &msg) {
  prev_state_ = cur_state_;
  cur_state_ = msg;
  if (!end_flag_)
    RampPlanner::Plan();
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
  double error = 0 - yaw;
  cur_integral_ += error*dt;
  double output = kp_*error + (std::abs(cur_integral_)<1 ? ki_*cur_integral_ : 0) + (dt>0.01 ? kd_*(error-prev_error_)/dt : 0);
  ROS_INFO("P = %f | I = %f | D = %f",kp_*error, ki_*cur_integral_, kd_*(error-prev_error_)/dt);
  ROS_INFO("PID output = %f",output);
  prev_error_ = error;

  // Generate goal
  double v = cur_state.twist.twist.linear.x + accel_*dt +0.4;

  geometry_msgs::Twist control_msg;

  control_msg.linear.x = v<target_vel_ ? v : target_vel_;
  control_msg.angular.z = output;
  goal.traj.commands.push_back(control_msg);
  goal.traj.states.push_back(cur_state);
  ++T_;

  goal_flag_ = v>=target_vel_ ? true : false;  // Ramp completion flag

  return goal;
}

void RampPlanner::SendTrajectory(ilqr_loco::TrajExecGoal &goal) {
  ac_.sendGoal(goal,
              boost::bind(&RampPlanner::doneCb, this, _1, _2),
              boost::bind(&RampPlanner::activeCb, this),
              boost::bind(&RampPlanner::feedbackCb, this, _1));
}

void RampPlanner::Plan() {
 
  if(ros::Time::now() - start_time_ < ros::Duration(timeout_)) {
    ilqr_loco::TrajExecGoal goal = RampPlanner::GenerateTrajectory(prev_state_, cur_state_);
    RampPlanner::SendTrajectory(goal);
  }
  else {
    // Stop car after ramp timeout
    ROS_INFO("Timeout exceeded, stopping car");
    ilqr_loco::TrajExecGoal end_goal;
    geometry_msgs::Twist control_msg;
    control_msg.linear.x = 0.0;
    control_msg.angular.z = 0.0;
    end_goal.traj.commands.push_back(control_msg);
    end_goal.traj.commands.push_back(control_msg);
    end_goal.traj.states.push_back(cur_state_);
    end_goal.traj.states.push_back(cur_state_);
    RampPlanner::SendTrajectory(end_goal);
    end_flag_ = true;
  }
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "ramp_client");
  RampPlanner ramp;
  ros::spin();

  return 0;
}
