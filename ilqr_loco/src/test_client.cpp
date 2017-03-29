#include "ros/ros.h"
#include <actionlib/server/simple_action_server.h>
#include <actionlib/client/terminal_state.h>
#include <ilqr_loco/TrajExecAction.h>

#include <math.h>
#include <geometry_msgs/Twist.h>

#include <fstream>
#include <string>
#include <sstream>
// #include "iLQR_mpc.c" //TODO integrate generated c-code

class TestClient
{
public:
  TestClient();

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

TestClient::TestClient()
{
  ros::NodeHandle nh;
  // TODO subscriber for state
  // TODO subscriber for obstacle position

  actionlib::SimpleActionServer<ilqr_loco::TrajExecAction> as;
}

void TestClient::GetNewSensorInfo()
{
  // get feedback on state from server
}

// Calls iLQR_mpc.c to generate new trajectory
bool TestClient::GenerateTrajectory(double x_cur[10], double x_des[6], double obs[2], int T)
{

}

void TestClient::SendTrajectory()
{

}

// Main function. Reads new state estimate and obstacle position, generates trajectory
// based on that, and sends trajectory to executer.
void TestClient::Plan()
{

}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "TestClient");

  TestClient client;

  return 0;
}
