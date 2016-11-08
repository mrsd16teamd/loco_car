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

#define led_pin 13
#define kill_pin 20
#define disable_pin 21
#define esc_pin 22
#define servo_pin 23

ros::NodeHandle  nh;

Servo servo;
Servo esc;

double x, y, z, w;
long st, th;
char buf[200];
bool reverse_tap = 0;

bool disabled = 0;
bool kill = 0;

double mapf(double x, double in_min, double in_max, double out_min, double out_max)
{
    return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}


void cmd_vel_cb( const geometry_msgs::Twist& cmd_msg){
  x = cmd_msg.linear.x;
  y = cmd_msg.linear.y;
//  z = cmd_msg.linear.z;
  w = cmd_msg.angular.z;  

 if (!disabled) {
  digitalWrite(led_pin, HIGH-digitalRead(led_pin));  //toggle led  
 }
}


ros::Subscriber<geometry_msgs::Twist> sub("cmd_vel", cmd_vel_cb);

std_msgs::String out_msg;
ros::Publisher teensy("teensy", &out_msg);

void setup(){
  pinMode(led_pin, OUTPUT);
  pinMode(esc_pin, OUTPUT);
  pinMode(servo_pin, OUTPUT);
  pinMode(disable_pin, INPUT);
  pinMode(kill_pin, INPUT);
  attachInterrupt(disable_pin, disable_ISR, CHANGE);
  attachInterrupt(kill_pin, kill_ISR, CHANGE);
  nh.initNode();
  nh.subscribe(sub);
  nh.advertise(teensy);
  
  servo.attach(servo_pin,1000,2000); //attach it to pin A9/23
  esc.attach(esc_pin,1000,2000); //attach it to pin A8/22
//  analogWriteResolution(16);
//  Serial.begin(9600);
}

void loop(){
 
  
  nh.spinOnce();
  String out;
  out +=  "Throttle: " + String(x) + ", " + String(th) + '\t' + "Steering: " + String(w) + ", " + String(st) + '\t' + "Disabled: " + String(disabled) ;
  out.toCharArray(buf,200);
  out_msg.data = buf;
  teensy.publish( &out_msg );
  
if (!disabled) {
  st = mapf(w, -1.0, 1.0, 1000,2000); //maxes out at +/- 0.8
  servo.writeMicroseconds(st); 

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
}

else {
    digitalWrite(led_pin, LOW);  //toggle led  
    th = 1500;
    st = 1500;
    servo.writeMicroseconds(1500);
    esc.writeMicroseconds(1500);
}
  
  
  delay(20);
}

void disable_ISR() {
  disabled = digitalReadFast(disable_pin);
}

void kill_ISR() {
  while (1) {
    digitalWrite(led_pin, HIGH-digitalRead(led_pin));  //toggle led  
    delay(1000);
  }
}
