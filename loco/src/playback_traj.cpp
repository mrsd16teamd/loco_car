//
// MIT License
//
// Copyright (c) 2017 MRSD Team D - LoCo
// The Robotics Institute, Carnegie Mellon University
// http://mrsdprojects.ri.cmu.edu/2016teamd/
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// 

/*
Plays back pre-recorded control input sequence.
Format is "u [vx] [steer] [dt]"
  vx : m/s
  steer: radians
  dt : sec
*/

//TODO
// take filepath as argument so we don't need to recompile to change file

#include "ros/ros.h"
#include <math.h>
#include <geometry_msgs/Twist.h>

#include <fstream>
#include <string>
#include <sstream>

class Player
{
public:
  Player();

private:
  ros::Publisher cmd_pub;

  void playControls(std::ifstream& controls);
};

Player::Player()
{
  ros::NodeHandle nh;

  cmd_pub = nh.advertise<geometry_msgs::Twist>("cmd_vel", 3);

  // Try to open command input file
  // TODO change this so we don't need to compile to change control file
  std::string filepath_str;
  if (nh.getParam("/playback_traj/file_path", filepath_str)){
    ROS_INFO("Playing %s", filepath_str.c_str());
  }
  else{
    ROS_INFO("Failed to read param.");
    return;
    // filepath_str = "/home/parallels/ros/src/loco/controlinputs/controls.txt";
  }
  const char* filepath = filepath_str.c_str();

  std::ifstream controls(filepath);
  if (!controls.is_open()){
    ROS_INFO("Failed to open control file.");
    return;
  }

  ROS_INFO("Started cmd_vel playback node. Sending commands in 3 seconds...");
  ros::Duration(3).sleep();

  // Call playback function
  playControls(controls);
}

void Player::playControls(std::ifstream& controls){
  std::string line;
  while (getline(controls,line) && ros::ok())
  {
    float vx, steer, dt;
    std::string first;
    std::stringstream data(line);

    data >> first >> vx >> steer >> dt;

    if(first=="end") {
      ROS_INFO("Finished playback.");
      return;
    }
    else if (first != "u"){
      ROS_INFO("Invalid input.");
      return;
    }

    geometry_msgs::Twist msg;
    msg.linear.x = vx;
    msg.angular.z = steer;

    ROS_INFO("Sending vx: %f, steer: %f", vx, steer);

    cmd_pub.publish(msg);
    ros::spinOnce();

    ros::Duration(dt).sleep();
  }
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "player_traj");

  Player player;

  ROS_INFO("Done.");

  return 0;
}
