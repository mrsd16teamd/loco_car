/*
Executes most recent trajectory given by iLQR.
Sends feedback to action client (planner) in terms of "almost done", "done".
*/

#include "ilqr_executer.h"

// provides action to execute plans
void iLQR_Executer::execute_trajectory(const ilqr_loco::TrajExecGoalConstPtr &goal){
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
      // TODO instead of making new message, send from goal message
      geometry_msgs::Twist msg = goal->traj.commands[i];
      ROS_INFO("Publishing %f and %f", msg.linear.x, msg.angular.z);
      cmd_pub.publish(msg);
      ros::spinOnce();

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

};


int main(int argc, char** argv)
{
  ros::init(argc, argv, "ilqr_executer");

  iLQR_Executer executer;
  ros::spin();

  return 0;
}
