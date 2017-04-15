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
      cmd_pub.publish(goal->traj.commands[i]);
      ros::spinOnce();

      int steps_left = goal->traj.commands.size() - i;
      if (steps_left>1)
        loop_rate.sleep();
      // feedback.steps_left =  steps_left;
      // as.publishFeedback(feedback);
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
