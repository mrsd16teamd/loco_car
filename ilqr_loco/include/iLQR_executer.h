#include "ros/ros.h"

#include <math.h>
#include <geometry_msgs/Twist.h>

#include <fstream>
#include <string>
#include <sstream>


class iLQR_Executer
{
public:
  iLQR_Executer();

private:
  ros::Publisher cmd_pub;

  void ExecuteTrajectory(int steps);
};
