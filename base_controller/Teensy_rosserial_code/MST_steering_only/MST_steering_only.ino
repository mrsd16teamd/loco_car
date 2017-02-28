/* 
 *  MRSD Team D '16 - Team LoCo - Evasive Maneuvers and Drifting for Autonomous Vehicles - Teensy Firmware
 *  This sets up the Teensy as a ROS subscriber for steering command inputs and also a ROS publisher on the /teensy topic for diagnostics
 *  The teensy then sends out PWM output to the servo to control the steering
*/

#if (ARDUINO >= 100)
#include <Arduino.h>
#else
#include <WProgram.h>
#endif

#include <Servo.h>
#include <ros.h>
//#include <std_msgs/UInt16.h>
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

Servo servo;


//steering angle input in rads
double steer; 
//double steer_zero = 0.22;
//double steer_min = 0.99;
//double steer_max = -0.55;
double steer_zero = 0;          //maxes out at +/- 0.77 rads = +/- 44.65 degs
double steer_min = 0.7777;
double steer_max = -0.7777;

//pwm value to write to the servo, zeros at 1385 in current steering arm config
long pwm_val;       
long pwm_zero = 1385;
long pwm_min = 1135;
long pwm_max = 1635;
//long pwm_min = 885;
//long pwm_max = 1885;


char buf[200];
unsigned long last_received;
const unsigned long timeout = 500; //timeout in ms before resetting steering to 0

bool disabled = 0;

double mapf(double x, double in_min, double in_max, double out_min, double out_max)
{
  return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}


void servo_cb(const std_msgs::Float64& steering_angle) {
  steer = steering_angle.data;
  servo.attach(servo_pin, pwm_min , pwm_max);
  last_received = millis();

}

ros::Subscriber<std_msgs::Float64> sub("/commands/servo/position", servo_cb);

std_msgs::String out_msg;
ros::Publisher teensy("teensy", &out_msg);

void setup() {
  pinMode(led_pin, OUTPUT);
  pinMode(servo_pin, OUTPUT);
  pinMode(on_pin, OUTPUT);
  pinMode(off_pin, OUTPUT);
  pinMode(disable_pin, INPUT);
  
  attachInterrupt(disable_pin, disable_ISR, CHANGE);
  
  nh.getHardware()->setBaud(115200);
  nh.initNode();
  nh.subscribe(sub);
  nh.advertise(teensy);

  servo.attach(servo_pin, pwm_min, pwm_max); //attach it to pin A9/23


  // just to show it's alive
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

void loop() {

  unsigned long elapsed = millis() - last_received;

  if (elapsed > timeout && pwm_val == pwm_zero) { //if no new commands are received, detach servo to prevent stressing it out
    servo.detach();
  }

  nh.spinOnce();
  String out;

  out += "Steering Angle: " + String(steer) +" rads | ";
  out += String(steer*RAD_TO_DEG, 2) + " degs \t";
  out += "PWM: " + String(pwm_val);
//  out += "Disabled: " + String(disabled)                    + '\t';
//  out += "Elapsed: " + elapsed;
  out.toCharArray(buf, 200);
  out_msg.data = buf;
  teensy.publish( &out_msg );

  if (!disabled) {

    pwm_val = mapf(steer, steer_min, steer_max, pwm_min, pwm_max); 
    servo.writeMicroseconds(pwm_val);

    digitalWrite(led_pin, LOW);
  }

  else {  //when disabled
    
    pwm_val = pwm_zero;
    servo.writeMicroseconds(pwm_val);
    servo.detach();

    digitalWrite(led_pin, HIGH);
  }

  delay(10);
}

void disable_ISR() {

  disabled = digitalRead(disable_pin);

  if (disabled) {
    digitalWrite(off_pin, HIGH);  //send pulse to MOSFET switch OFF pin to turn off VESC
    delay(10);
    digitalWrite(off_pin, LOW);

  }

  else {
    digitalWrite(on_pin, HIGH);  //send pulse to MOSFET switch ON pin to turn on VESC
    delay(10);
    digitalWrite(on_pin, LOW);

  }

  steer = steer_zero;
  pwm_val = pwm_zero;

}


