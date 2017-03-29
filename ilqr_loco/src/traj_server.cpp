/*
Executes most recent trajectory given by iLQR.
Sends feedback to action client (planner) in terms of "almost done", "done".
*/

#include "traj_server.h"

// provides action to execute plans
void TrajServer::execute_trajectory(const ilqr_loco::TrajExecGoalConstPtr &goal){
  // TODO check that states and commands are right length

  bool success = true;
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
    else{
      ROS_INFO("Publishing command: %f, %f", goal->traj.commands[i].linear.x, goal->traj.commands[i].angular.z);
      cmd_pub.publish(goal->traj.commands[i]);
      ros::spinOnce();

      ROS_INFO("Some of states: %f, %f", goal->traj.states[i].pose.pose.position.x,
          goal->traj.states[i].twist.twist.linear.x);

      feedback.steps_left = 1; //TODO change this
      as.publishFeedback(feedback);
      ros::Duration(1.0).sleep();
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
