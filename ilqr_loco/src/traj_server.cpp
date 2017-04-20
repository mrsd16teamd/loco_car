/*
Executes most recent trajectory given by iLQR.
Sends feedback to action client (planner) in terms of "almost done", "done".
*/

#include "traj_server.h"

void TrajServer::LoadParams()
{
    TRYGETPARAM("old_msg_discard_thres", old_msg_thres)
    TRYGETPARAM("kp_heading", kp_)
    TRYGETPARAM("ki_heading", ki_)
    TRYGETPARAM("kd_heading", kd_)
}

void TrajServer::PublishPath(const ilqr_loco::TrajExecGoalConstPtr &goal)
{
  nav_msgs::Path path_msg;
  path_msg.header.stamp = goal->traj.header.stamp;
  path_msg.header.frame_id = "map";
  std::vector<geometry_msgs::PoseStamped> poses(goal->traj.states.size());

  for (int i=0; i < goal->traj.commands.size(); i++)
  {
    poses.at(i).pose.position.x = goal->traj.states[i].pose.pose.position.x;
    poses.at(i).pose.position.y = goal->traj.states[i].pose.pose.position.y;
    poses.at(i).pose.orientation = goal->traj.states[i].pose.pose.orientation;
  }
  path_msg.poses = poses;
  path_pub.publish(path_msg);
  ros::spinOnce();
}

void TrajServer::stateCb(const nav_msgs::Odometry &msg)
{
  cur_yaw_  = tf::getYaw(msg.pose.pose.orientation);
}

geometry_msgs::Twist TrajServer::pid_correct_yaw(geometry_msgs::Twist orig_twist, nav_msgs::Odometry state)
{
  double yaw_des = tf::getYaw(state.pose.pose.orientation);
  double orig_steer = orig_twist.angular.z;

  // Correct steering to compensate for yaw error
  double error = yaw_des - cur_yaw_;
  cur_integral_ += error*dt;
  double steer = orig_steer + kp_*error + (std::abs(cur_integral_)<1 ? ki_*cur_integral_ : 0) + (dt>0.01 ? kd_*(error-prev_error_)/dt : 0);
  prev_error_ = error;

  geometry_msgs::Twist new_twist;
  new_twist.linear.x = orig_twist.linear.x;
  new_twist.angular.z = steer;

  return new_twist;
}

// provides action to execute plans
void TrajServer::execute_trajectory(const ilqr_loco::TrajExecGoalConstPtr &goal){
  // TODO? check that states and commands are right length
  ROS_INFO("%s: Received trajectory.", traj_action.c_str());

  bool success = true;

  double timestep = goal->traj.timestep;
  double traj_start_time = (goal->traj.header.stamp).toSec();

  ros::Rate loop_rate(1.0/timestep);

  // ROS_INFO("Executing trajectory in mode %d", goal->traj.mode); // TODO print client name

  for (int i=0; i < goal->traj.commands.size(); i++)
  {
    // check that preempt has not been requested by the client
    if (as.isPreemptRequested() || !ros::ok())
    {
      ROS_INFO("%s: Preempted at %dth command.", traj_action.c_str(), i);
      as.setPreempted();
      success = false;
      break;
    }
    // check that commands in plan are not too old
    else if ((ros::Time::now().toSec() - (traj_start_time + (i*timestep))) > old_msg_thres)
    {
     ROS_INFO("%s: Ignoring old command.", traj_action.c_str());
     continue;
    }
    else
    {
      if (goal->traj.mode == 1) {
        geometry_msgs::Twist pid_twist = pid_correct_yaw(goal->traj.commands[i], goal->traj.states[i]);
        cmd_pub.publish(pid_twist);
      }
      else {
        cmd_pub.publish(goal->traj.commands[i]);
      }
      ros::spinOnce();

      int steps_left = goal->traj.commands.size() - i;
      if (steps_left>1)
        loop_rate.sleep();
    }
  }

  PublishPath(goal);   // For visualization

  if(success)
  {
    ROS_INFO("%s: Finished publishing trajectory", traj_action.c_str());
    result.done = 1;
    as.setSucceeded(result);
  }
}


int main(int argc, char** argv)
{
  ros::init(argc, argv, "traj_server");

  TrajServer executer;
  ros::spin();

  return 0;
}
