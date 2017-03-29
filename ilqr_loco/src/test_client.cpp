#include "ros/ros.h"
#include <actionlib/client/simple_action_client.h>
#include <actionlib/client/terminal_state.h>
#include <ilqr_loco/TrajExecAction.h>

#include <math.h>
#include <geometry_msgs/Twist.h>

#include <fstream>
#include <string>
#include <sstream>
// #include "iLQR_mpc.c" //TODO integrate generated c-code

typedef actionlib::SimpleActionClient<ilqr_loco::TrajExecAction> Client;
//
// class TestClient
// {
// public:
//   TestClient();
//   void Plan();
//
// private:
//   ros::NodeHandle nh;
//   ros::Publisher cmd_pub;
//   actionlib::SimpleActionClient<ilqr_loco::TrajExecAction> ac;
//
//   void GetNewSensorInfo();
//   void SendTrajectory();
//   bool GenerateTrajectory(double x_cur[10], double x_des[6], double obs[2], int T);
// };
//
// TestClient::TestClient()
// {
//   ac = ac1;
//     // TODO subscriber for state
//   // TODO subscriber for obstacle position
//
// }
//
// void TestClient::GetNewSensorInfo()
// {
//   // get feedback on state from server
// }
//
// // Calls iLQR_mpc.c to generate new trajectory
// bool TestClient::GenerateTrajectory(double x_cur[10], double x_des[6], double obs[2], int T)
// {
//
// }
//
// void TestClient::SendTrajectory()
// {
//
// }
//
// // Main function. Reads new state estimate and obstacle position, generates trajectory
// // based on that, and sends trajectory to executer.
// void TestClient::Plan()
// {
//   ROS_INFO("Waiting for action server to start.");
//   // wait for the action server to start
//   ac.waitForServer(); //will wait for infinite time
// }

int main(int argc, char** argv)
{
  ros::init(argc, argv, "test_client");

  Client ac("traj_executer", true);
  ROS_INFO("Waiting for action server to start.");
  ac.waitForServer(); //will wait for infinite time

  ROS_INFO("Action server started, sending goal.");
  ilqr_loco::TrajExecGoal goal;
  int traj_length = 20;
  for (int i=0; i<traj_length; i++)
  {
    geometry_msgs::Twist msg;
    msg.linear.x = 0.5;
    msg.angular.z = 0.5;
    goal.traj.commands.push_back(msg);
  }

  ac.sendGoal(goal);

  ros::spin();

  // TestClient client;
  // client.Plan();

  return 0;
}
