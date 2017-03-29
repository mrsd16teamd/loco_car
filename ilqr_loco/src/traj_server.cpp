/*
Executes most recent trajectory given by iLQR.
Sends feedback to action client (planner) in terms of "almost done", "done".
*/

#include "traj_server.h"

// provides action to execute plans
void TrajServer::execute_trajectory(const ilqr_loco::TrajExecGoalConstPtr &goal){
  // TODO? check that states and commands are right length

  bool success = true;
  ROS_INFO("Executing trajectory."); // TODO print client name

  double timestep = goal->traj.timestep;
  double traj_start_time = (goal->traj.header.stamp).toSec();

  for (int i=0; i < goal->traj.commands.size(); i++)
  {
    double now = ros::Time::now().toSec();
    double cmd_planned_time = traj_start_time + (i*timestep);

    // check that preempt has not been requested by the client
    if (as.isPreemptRequested() || !ros::ok())
    {
      ROS_INFO("%s: Preempted", traj_action.c_str());
      as.setPreempted();
      success = false;
      break;
    }
    // check that commands in plan are not too old
    else if ((now - cmd_planned_time) > old_msg_thres)
    {
      ROS_INFO("Ignoring old command.");
      continue;
    }
    else{
      ROS_INFO("Publishing command: %f, %f", goal->traj.commands[i].linear.x, goal->traj.commands[i].angular.z);
      cmd_pub.publish(goal->traj.commands[i]);
      ros::spinOnce();

      feedback.steps_left =  goal->traj.commands.size() - i;
      as.publishFeedback(feedback);
    }
  }

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
