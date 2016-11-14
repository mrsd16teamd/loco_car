#include <ros/ros.h>
#include <math.h>
#include <geometry_msgs/Twist.h>
#include <ackermann_msgs/AckermannDriveStamped.h>

class VelMsgConverter
{
public:
  VelMsgConverter()
  {
    // Initialize the publisher and subscriber
    vel_pub = n_.advertise<ackermann_msgs::AckermannDriveStamped>("/cmd_vel_ack", 1);
    vel_sub = n_.subscribe("/cmd_vel",1,&VelMsgConverter::velCallback,this);
  }

  private:
    ros::NodeHandle n_;
    ros::Publisher vel_pub;
    ros::Subscriber vel_sub;

    // Initialize parameters
    float wheelbase = 0.255;
    /*
    if (n_.getParam("chassis/wheelbase",wheelbase)) {
      n_.getParam("chassis/wheelbase",wheelbase);
      ROS_INFO("Loaded chassis parameters");
      ROS_INFO("wheelbase=%f", wheelbase);
    }
    else
      ROS_INFO("parameters not found!");
      */


  float vel_to_steering_angle(float v, double omega) {
    if(omega==0 || v==0)
      return 0;
    float radius = v/omega;
    float steering_angle = atan(wheelbase/radius);
    return steering_angle;
  }

  void velCallback(const geometry_msgs::Twist& twist_msg) {
    ackermann_msgs::AckermannDriveStamped ack_msg;

    float wheelbase = 0.255;
    float lin_vel = twist_msg.linear.x;
    ack_msg.header.stamp = ros::Time::now();
    ack_msg.header.frame_id = "base_link";

    ack_msg.drive.speed = lin_vel;
    ack_msg.drive.steering_angle = vel_to_steering_angle(lin_vel, twist_msg.angular.z);

    vel_pub.publish(ack_msg);

  }
};//End of class

int main(int argc, char** argv) {
  ros::init(argc, argv, "cmd_vel_to_ackermann_drive");
  VelMsgConverter ackermann_drive_pub;
  ros::spin();
  return 0;
}
