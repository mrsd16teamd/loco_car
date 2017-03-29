/*
Continuously calls iLQR to generate trajectories, based on changing states,
desired states, and obstacle positions.
*/

#include "ilqr_planner.h"

iLQR_Planner::iLQR_Planner()
{
  ros::NodeHandle nh;
  cmd_pub = nh.advertise<geometry_msgs::Twist>("cmd_vel", 3);
  // TODO subscriber for state
  // TODO subscriber for obstacle position

  actionlib::SimpleActionServer<actionlib_tutorials::FibonacciAction> as_;
}

void iLQR_Planner::GetNewSensorInfo()
{
  // get feedback on state from server
}

// Calls iLQR_mpc.c to generate new trajectory
bool iLQR_Planner::GenerateTrajectory(double x_cur[10], double x_des[6], double obs[2], int T)
{

}

void iLQR_Planner::SendTrajectory()
{

}

// Main function. Reads new state estimate and obstacle position, generates trajectory
// based on that, and sends trajectory to executer.
void iLQR_Planner::Plan()
{

}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "iLQR_planner");

  iLQR_Planner planner;

  return 0;
}
