/*
Continuously calls iLQR to generate trajectories, based on changing states,
desired states, and obstacle positions.
*/

#include "ros/ros.h"
#include <actionlib/server/simple_action_server.h>
#include <actionlib/client/terminal_state.h>
#include <ilqr_loco/TrajExecAction.h">

#include <math.h>
#include <geometry_msgs/Twist.h>

#include <fstream>
#include <string>
#include <sstream>
// #include "iLQR_mpc.c" //TODO integrate generated c-code


class iLQR_Planner
{
public:
  iLQR_Planner();

private:
  ros::Publisher cmd_pub;
  // TODO subscriber for state
  // TODO subscriber for obstacle position
  // TODO most recent state estimate
  // TODO most recent obstacle position

  void GetNewSensorInfo();
  void SendTrajectory();
  bool GenerateTrajectory(double x_cur[10], double x_des[6], double obs[2], int T);
  void Plan();
};


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
