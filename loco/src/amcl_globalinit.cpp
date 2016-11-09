#include "ros/ros.h"
#include <cstdlib>
#include "std_srvs/Empty.h"
#include "geometry_msgs/PoseArray.h"

int main(int argc, char **argv)
{
  ros::init(argc, argv, "amcl_globalinit");

  ros::NodeHandle n;
  ros::topic::waitForMessage<geometry_msgs::PoseArray>("particlecloud",ros::Duration(5));
  ros::ServiceClient client = n.serviceClient<std_srvs::Empty>("global_localization");
  std_srvs::Empty srv;
  client.call(srv);
  return 0;
}
