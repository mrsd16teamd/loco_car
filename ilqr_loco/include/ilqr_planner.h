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

class iLQR_planner
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
