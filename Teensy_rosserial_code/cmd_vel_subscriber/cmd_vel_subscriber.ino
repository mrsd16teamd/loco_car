/*
 * rosserial Servo Control Example
 *
 * This sketch demonstrates the control of hobby R/C servos
 * using ROS and the arduiono
 * 
 * For the full tutorial write up, visit
 * www.ros.org/wiki/rosserial_arduino_demos
 *
 * For more information on the Arduino Servo Library
 * Checkout :
 * http://www.arduino.cc/en/Reference/Servo
 */

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

ros::NodeHandle  nh;

Servo servo;
Servo esc;

double x, y, z, w;
long st, th;
char buf[200];
bool reverse_tap = 0;

double mapf(double x, double in_min, double in_max, double out_min, double out_max)
{
    return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}


void cmd_vel_cb( const geometry_msgs::Twist& cmd_msg){
  x = cmd_msg.linear.x;
  y = cmd_msg.linear.y;
//  z = cmd_msg.linear.z;
  w = cmd_msg.angular.z;

//  st = mapf(w, -1.0, 1.0, 0, 65535); //maxes out at +/- 0.8
//  analogWrite(23, st);
  
  st = mapf(w, -1.0, 1.0, 1000,2000); //maxes out at +/- 0.8
  servo.writeMicroseconds(st); //set servo PWM pulse high microseconds, should be from 1000-2000  
 
//  th = mapf(x, -1.0, 1.0, 0, 65535); //minimum fwd speed at ~0.15-0.8max, min reverse spd at -0.20, max -0.8
//  analogWrite(22, th);

//  th = mapf(x, -1.0, 1.0, 1000,2000); //minimum fwd speed at ~0.15-0.8max, min reverse spd at -0.20, max -0.8
  if (x > 0) {
    th = mapf(x, 0, 1.0, 1580, 1836); //hand tuned values. default to 1500, 2000 if problems
    esc.writeMicroseconds(th);
    reverse_tap = 1;
  }
  else if ( x <0) {
    th = mapf(x, -1.0, 0, 962 , 1390); //hand tuned values. default to 1000, 1500 if problems
    if (reverse_tap) {
      esc.writeMicroseconds(1350);
      delay(100);
      esc.writeMicroseconds(1500);
      delay(80);
      reverse_tap = 0;
    }

    
    esc.writeMicroseconds(th);
  }

  else esc.writeMicroseconds(1500);
  
//  th = mapf(x, -1.0, 1.0, 1000,2000); //minimum fwd speed at ~0.15-0.8max, min reverse spd at -0.20, max -0.8
//  esc.writeMicroseconds(th);
 
  digitalWrite(13, HIGH-digitalRead(13));  //toggle led  
}


ros::Subscriber<geometry_msgs::Twist> sub("cmd_vel", cmd_vel_cb);

std_msgs::String out_msg;
ros::Publisher teensy("teensy", &out_msg);

void setup(){
  pinMode(13, OUTPUT);
  pinMode(23, OUTPUT);
  pinMode(22, OUTPUT);
  nh.initNode();
  nh.subscribe(sub);
  nh.advertise(teensy);
  
  servo.attach(23,1000,2000); //attach it to pin A9/23
  esc.attach(22,1000,2000); //attach it to pin A8/22
//  analogWriteResolution(16);
//  Serial.begin(9600);
}

void loop(){
  nh.spinOnce();
  String out;
  out +=  "Throttle: " + String(x) + ", " + String(th) + '\t' + "Steering: " + String(w) + ", " + String(st) ;
  out.toCharArray(buf,200);
  out_msg.data = buf;
  teensy.publish( &out_msg );
  

  
//  Serial.println(out);
  delay(20);
}
