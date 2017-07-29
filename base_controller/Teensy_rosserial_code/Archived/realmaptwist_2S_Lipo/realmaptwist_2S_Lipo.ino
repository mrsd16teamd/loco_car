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

#if (ARDUINO >= 100)
 #include <Arduino.h>
#else
 #include <WProgram.h>
#endif

#include <Servo.h>
#include <ros.h>
#include <std_msgs/UInt16.h>
#include <std_msgs/String.h>
#include <geometry_msgs/Twist.h>
//#include <ackermann_msgs/AckermannDriveStamped.h>

#define led_pin 13
#define kill_pin 20
#define disable_pin 21
#define esc_pin 22
#define servo_pin 23

ros::NodeHandle  nh;

Servo servo;
Servo esc;

double x, w;
long steer, throttle;
char buf[200];
unsigned long last_received;
const unsigned long timeout = 1000; //timeout in ms before resetting steering and throttle to 0

bool disabled = 0;
bool kill = 0;


double mapf(double x, double in_min, double in_max, double out_min, double out_max)
{
    return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

/*
void cmd_vel_cb(const ackermann_msgs::AckermannDriveStamped& cmd_msg){
  x = cmd_msg.drive.speed;
  w = cmd_msg.drive.steering_angle;
  last_received = millis();

}*/
void cmd_vel_cb(const geometry_msgs::Twist& cmd_msg){
  x = cmd_msg.linear.x;
  w = cmd_msg.angular.z;
  last_received = millis();

}
//ros::Subscriber<ackermann_msgs::AckermannDriveStamped> sub("cmd_vel_ack", cmd_vel_cb);
ros::Subscriber<geometry_msgs::Twist> sub("cmd_vel", cmd_vel_cb);

std_msgs::String out_msg;
ros::Publisher teensy("teensy", &out_msg);

void setup(){
  pinMode(led_pin, OUTPUT);
  pinMode(esc_pin, OUTPUT);
  pinMode(servo_pin, OUTPUT);
  pinMode(disable_pin, INPUT);
//  pinMode(kill_pin, INPUT);
  attachInterrupt(disable_pin, disable_ISR, CHANGE);
//  attachInterrupt(kill_pin, kill_ISR, CHANGE);
  nh.getHardware()->setBaud(115200);
  nh.initNode();
  nh.subscribe(sub);
  nh.advertise(teensy);

  servo.attach(servo_pin,1000,2000); //attach it to pin A9/23
  esc.attach(esc_pin,1000,2000); //attach it to pin A8/22

}

void loop(){

  unsigned long elapsed = millis() - last_received;

  if (elapsed > timeout) {
    x = 0;
    w = 0;
  }

  nh.spinOnce();
  String out;
  out +=  "Throttle: " + String(x) + ", " + String(throttle) + '\t' + "Steering: " + String(w) + ", " + String(steer) + '\t' + "Disabled: " + String(disabled) + "\t Elapsed: " + elapsed ;
  out.toCharArray(buf,200);
  out_msg.data = buf;
  teensy.publish( &out_msg );

  if (!disabled) {

    steer = mapf(w, 0.8, -0.8, 1160,1960); //maxes out at +/- 0.8
    servo.attach(servo_pin,1000,2000);

//    servo.writeMicroseconds(steer);
//    if (abs(steer-1500) <= 10) {
//          steer = 1500;
//    }

    if (x>0) {
      //       throttle = mapf(x, 0, 0.4, 1535, 1560); //hand tuned values. default to 1500, 2000 if problems
      if( x < 0.394) {
        throttle = 1550; //minimum to start moving
      }

      else if (x > 1.55 && (w > 0.60 || w < - 0.60)) {   //if turbo mode on teleop is activated, make sure steering is high and spin fast enough so it can drift
//      else if (x > 1.55) {
        throttle = 1666;
      }

      else if  (x > 1.5) {
        throttle = 1570;  //maximum that's safe to go
      }

      else {

        throttle = (float(x)*20.08)+1542.06; //emperically determined equation for forward motion
      }
    }

    else if (x < 0) {
//      throttle = mapf(x, -0.4, 0, 1410, 1430); //hand tuned values. default to 1500, 2000 if problems
      if( x > -0.394) {
        throttle = 1433; //minimum to start moving
      }
      else if (x < -1.55 && (w > 0.60 || w < - 0.60)) {  //if turbo mode on teleop is activated, make sure steering is high and spin fast enough so it can drift
//     else if (x < -1.55) {
        throttle = 1324;
    }
      else if  (x < -1.5) {
        throttle = 1424;  //maximum that's safe to go
      }

      else {


          throttle = (float(x)*10.08)+1437.032; //emperically determined equation for forward motion transformed for backward motion UNSTABLE
      }

    }

    else {
      throttle = 1500;
    }

    esc.writeMicroseconds(throttle);
    servo.writeMicroseconds(steer);

    if (elapsed > 500 && steer == 1560){
      servo.detach();
    }

    digitalWrite(led_pin, LOW);
  }

  else {  //when disabled

    throttle = 1500;
    steer = 1560;
    servo.writeMicroseconds(1560);
    servo.detach();
    esc.writeMicroseconds(1500);
    digitalWrite(led_pin, HIGH);

  }

  delay(10);
}

void disable_ISR() {

  disabled = digitalRead(disable_pin);
  w = 0;
  x = 0;
  throttle = 1500;
  steer = 1500;
  esc.writeMicroseconds(1500);
  servo.writeMicroseconds(1560);
  servo.detach();

}
