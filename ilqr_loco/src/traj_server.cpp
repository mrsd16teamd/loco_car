/*
Executes most recent trajectory given by iLQR.
Sends feedback to action client (planner) in terms of "almost done", "done".
*/

#include "traj_server.h"

void TrajServer::LoadParams()
{
  try
	{
    nh.getParam("old_msg_discard_thres", old_msg_thres);
  }
  catch(...)
	{
    ROS_ERROR("Please put all params into yaml file, and load it.");
  }
}

void TrajServer::ClampControls(Eigen::Vector2d &u)
{
  u(0) = std::min(throttle_lims(1), std::max(u(0), throttle_lims(0)));
  u(1) = std::min(steering_lims(1), std::max(u(1), steering_lims(0)));
}

void TrajServer::FillVecFromOdom(const nav_msgs::Odometry &odom, Eigen::VectorXd &v)
{
  v(0) = odom.pose.pose.position.x;
  v(1) = odom.pose.pose.position.y;
  v(2) = tf::getYaw(odom.pose.pose.orientation);
  v(3) = odom.twist.twist.linear.x;
  v(4) = odom.twist.twist.linear.y;
  v(5) = odom.twist.twist.angular.z;
  v(6) = last_u(0); // TODO figure out what to do here
  v(7) = last_u(1);
}

void TrajServer::FillVecFromTwist(const geometry_msgs::Twist &twist, Eigen::Vector2d &v)
{
  v(0) = twist.linear.x;
  v(1) = twist.angular.z;
}

void TrajServer::FillTwistFromVec(geometry_msgs::Twist &twist, const Eigen::Vector2d &v)
{
  twist.linear.x = v(0);
  twist.angular.z = v(1);
}

void TrajServer::stateCb(const nav_msgs::Odometry &odom)
{
  FillVecFromOdom(odom, cur_state);
}

// provides action to execute plans
void TrajServer::execute_trajectory(const ilqr_loco::TrajExecGoalConstPtr &goal){
  // TODO? check that states and commands are right length

  bool success = true;

  double timestep = goal->traj.timestep;
  double traj_start_time = (goal->traj.header.stamp).toSec();

  // ros::Rate loop_rate(1/timestep);
  ros::Rate loop_rate(50); //TODO make this work without hard-coding?

  ROS_INFO("Executing trajectory."); // TODO print client name

  for (int i=0; i < goal->traj.commands.size(); i++)
  {
    // check that preempt has not been requested by the client
    if (as.isPreemptRequested() || !ros::ok())
    {
      ROS_INFO("%s: Preempted", traj_action.c_str());
      as.setPreempted();
      success = false;
      break;
    }
    // check that commands in plan are not too old
    else if ((ros::Time::now().toSec() - (traj_start_time + (i*timestep))) > old_msg_thres)
    {
     ROS_INFO("Ignoring old command.");
     continue;
    }
    else
    {
      if (goal->traj.mode == 0) // old behavior - executing control sequence
      {
        cmd_pub.publish(goal->traj.commands[i]);
        ros::spinOnce();

        int steps_left = goal->traj.commands.size() - i;
        if (steps_left>1)
          loop_rate.sleep();
        // feedback.steps_left =  steps_left;
        // as.publishFeedback(feedback);
      }
      else if (goal->traj.mode == 1) // new behavior - execute control policy
      {
        // Get x, u, l, L into eigen matrices from messages
        FillVecFromOdom(goal->traj.states[i], x);
        if (i==0){
          x(6) = x(7) = 0;
        }
        else{
          x(6) = goal->traj.commands[i-1].linear.x;
          x(7) = goal->traj.commands[i-1].angular.z;
        }
        FillVecFromTwist(goal->traj.commands[i], u);

        l(0) = goal->traj.l.data[(2*i)];
        l(1) = goal->traj.l.data[(2*i)+1];
        for (int j=0; j<2; j++) {
          for (int k=0; k<8; k++) {
            L(j,k) = goal->traj.L.data[(i*16)+(j*2)+k]; // TODO change this index
          }
        }

        // Adjust control with control gains
        // Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
        // std::cout << "cur_state: " << cur_state.format(CommaInitFmt) << '\n';
        // std::cout << "l: " << l.format(CommaInitFmt) << '\n';
        // std::cout << "L: " << L.format(CommaInitFmt) << '\n';
        // std::cout << "u: " << u.format(CommaInitFmt) << '\n';
        u = u - l;
        // std::cout << "u-l: " << u.format(CommaInitFmt) << '\n';
        dx = cur_state - x;
        // std::cout << "x: " << x.format(CommaInitFmt) << '\n';
        // std::cout << "dx: " << u.format(CommaInitFmt) << '\n';
        u = u + L*dx;
        // std::cout << "u: " << u.format(CommaInitFmt) << "\n\n\n";

        // Publish command and save the commands as the previous cmds
        // TODO clamp or boxQP control inputs
        ClampControls(u);
        geometry_msgs::Twist control;
        FillTwistFromVec(control, u);
        last_u = u;
        cmd_pub.publish(control);
        ros::spinOnce();

        int steps_left = goal->traj.commands.size() - i;
        if (steps_left>1)
          loop_rate.sleep();
      }
    }
  }

  // For visualization
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

  if(success)
  {
    ROS_INFO("Finished publishing trajectory");
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
