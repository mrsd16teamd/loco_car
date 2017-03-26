#include "ros/ros.h"

#include <math.h>
#include <geometry_msgs/Twist.h>

#include <fstream>
#include <string>
#include <sstream>


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
