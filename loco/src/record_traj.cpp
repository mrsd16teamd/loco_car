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
Subscribes to cmd_vel topic and records what it sees, along with timesteps
between commands. Records are saved to txt file for playback later
Format is "u [vx] [steer] [dt]"
  vx : m/s
  steer: radians
  dt : sec
*/

// TODO
// add 'end'

#include "ros/ros.h"
#include <math.h>
#include <geometry_msgs/Twist.h>

#include <iostream>
#include <fstream>
#include <string>

class Recorder
{
public:
  Recorder();
  // ~Recorder(); TODO custom sigint handler to write 'end'?

private:
  std::ofstream output_file;
  ros::Time last_cmd_time;


  ros::Subscriber cmd_sub;
  void CommandCallback(const geometry_msgs::Twist::ConstPtr& msg);
};

Recorder::Recorder()
{
  ros::NodeHandle nh;
  cmd_sub = nh.subscribe("cmd_vel", 1, &Recorder::CommandCallback, this);

  last_cmd_time = ros::Time::now();

  // Open file to write to
  output_file.open("/home/parallels/ros/src/loco/controlinputs/recorded.txt");
  if (!output_file.is_open()){
    ROS_INFO("Failed to open output file.");
    return;
  }

  ROS_INFO("Started cmd_vel recorder.");
}

void Recorder::CommandCallback(const geometry_msgs::Twist::ConstPtr& msg)
{
  if(output_file.is_open())
  {
    float vx, steer, dt;

    vx = msg->linear.x;
    steer = msg->angular.z;

    ros::Time current_time = ros::Time::now();
    dt = (current_time-last_cmd_time).toSec();
    if(dt>0.3)
      dt = 0.0;
    last_cmd_time = current_time;

    output_file << "u " << vx << " " << steer << " " << dt << '\n';
    output_file.flush();
    std::cout << "Heard " << vx << " " << steer << " dt:" << dt << '\n';
  }
  else
  {
    ROS_INFO("Unable to write to output file.");
  }
}


int main(int argc, char** argv)
{
  ros::init(argc, argv, "record_traj");

  Recorder recorder;

  ros::spin();

  return 0;
}
