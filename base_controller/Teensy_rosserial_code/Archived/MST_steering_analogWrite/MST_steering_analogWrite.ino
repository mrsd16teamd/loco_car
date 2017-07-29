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
#include <std_msgs/Float64.h>
#include <std_msgs/String.h>
//#include <geometry_msgs/Twist.h>
//#include <ackermann_msgs/AckermannDriveStamped.h>

#define led_pin 13
#define on_pin 20
#define off_pin 19
#define disable_pin 21
#define esc_pin 22
#define servo_pin 23

ros::NodeHandle  nh;

double x;
double w = 0.22;
double actual_w, actual_deg;
const float pi = 3.14159;
const int freq = 100;
//long steer_zero = 567;
long steer_zero = 9076; // for 16 bit analogWrite
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


//void cmd_vel_cb(const geometry_msgs::Twist& cmd_msg){
//  x = cmd_msg.linear.x;
//  w = cmd_msg.angular.z;
//  last_received = millis();
//
//}

void cmd_vel_cb(const std_msgs::Float64& steering_angle){
//  x = cmd_msg.linear.x;
//  w = cmd_msg.angular.z;

  w = steering_angle.data;
//  servo.attach(servo_pin,885,1885);

  last_received = millis();

}

//ros::Subscriber<geometry_msgs::Twist> sub("cmd_vel", cmd_vel_cb);
ros::Subscriber<std_msgs::Float64> sub("/commands/servo/position", cmd_vel_cb);

std_msgs::String out_msg;
ros::Publisher teensy("teensy", &out_msg);

void setup(){
  pinMode(led_pin, OUTPUT);
  pinMode(esc_pin, OUTPUT);
  pinMode(servo_pin, OUTPUT);
  pinMode(disable_pin, INPUT);
  pinMode(on_pin, OUTPUT);
  pinMode(off_pin, OUTPUT);
  attachInterrupt(disable_pin, disable_ISR, CHANGE);
//  attachInterrupt(kill_pin, kill_ISR, CHANGE);
  nh.getHardware()->setBaud(115200);
  nh.initNode();
  nh.subscribe(sub);
  nh.advertise(teensy);

//  servo.attach(servo_pin,885,1885); //attach it to pin A9/23
//  esc.attach(esc_pin,1000,2000); //attach it to pin A8/22

  analogWriteFrequency(servo_pin, freq);
//  analogWriteResolution(12);
  analogWriteResolution(16);
  // just to show it's alive, LED pin pulse can be used to turn on the VESC as well
  digitalWrite(led_pin, HIGH);
  delay(100);
  digitalWrite(led_pin, LOW);
  delay(100);
  digitalWrite(led_pin, HIGH);
  delay(100);
  digitalWrite(led_pin, LOW);

  // send pulse to turn ON VESC
  digitalWrite(on_pin, HIGH);
  delay(10);
  digitalWrite(on_pin, LOW);
}

void loop(){

  unsigned long elapsed = millis() - last_received;

  nh.spinOnce();
  String out;
  actual_w = w-0.22;
  actual_deg = actual_w*180.0/pi;
  out +=  "Steering: " + String(actual_w) + ", " + String(actual_deg) + ", " + String(steer) + '\t' + "Disabled: " + String(disabled) + "\t Elapsed: " + elapsed ;
  out.toCharArray(buf,200);
  out_msg.data = buf;
  teensy.publish( &out_msg );

  if (!disabled) {


      steer = mapf(w, 0.9977, -0.5577, 7728 , 10408); //maxes out at +/- 0.77 rads = +/- 44.56 degs, with offset of 0.22 rad from center

    steer = long(steer);

    if (elapsed > timeout && w == 0.22){
      analogWrite(servo_pin, 0);
    }

    else{
      analogWrite(servo_pin, steer);
    }
    digitalWrite(led_pin, LOW);




  }

  else {  //when disabled

    steer = steer_zero;

    analogWrite(servo_pin, steer_zero);
    digitalWrite(led_pin, HIGH);
  }

  delay(10);
}

void disable_ISR() {

  disabled = digitalRead(disable_pin);

  if (disabled) {
    digitalWrite(off_pin, HIGH);
    delay(10);
    digitalWrite(off_pin, LOW);

  }

  else {
    digitalWrite(on_pin, HIGH);
    delay(10);
    digitalWrite(on_pin, LOW);

  }




//  throttle = 1500;
  steer = steer_zero;
//  esc.writeMicroseconds(throttle);
//  servo.writeMicroseconds(steer);
//  servo.detach();


}
