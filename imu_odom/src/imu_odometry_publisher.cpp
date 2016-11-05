#include "ros/ros.h"
#include "sensor_msgs/Imu.h"
#include "tf/transform_datatypes.h"
#include <tf/transform_broadcaster.h>
#include <nav_msgs/Odometry.h>
#include <cmath>
#include <string>

#define LOWPASS false

class ImuData{
  public:
    ImuData(): orientation{0,0,0},angular_velocity{0,0,0},
      linear_acceleration{0,0,0} { }
    float orientation[3];          //in roll, pitch yaw [radians]
    float angular_velocity[3];     //[rad/s]
    float linear_acceleration[3];  //[m/s^2]
};

class State2d{
  public:
    State2d(): x(0), y(0), t(0), vx(0), vy(0), w(0) { }
    float x, y, t, vx, vy, w;
};

class ImuOdometryPublisher {
  public:
    ImuOdometryPublisher();
    ~ImuOdometryPublisher();

    //void InitializeParams();
    void CalculateOdometry();

  private:
    ImuData latest_imu;
    State2d current_state;

    float lp_alpha_;
    float lp_beta_;
    float threshold_;
    float inactive_timeout_;

    ros::Time last_active_time;
    ros::Time last_update_time;

    ros::Subscriber imu_sub;
    ros::Publisher odom_pub;

    void SensorCallback(const sensor_msgs::Imu::ConstPtr& msg);
    void PublishOdometry();
};

ImuOdometryPublisher::ImuOdometryPublisher() {
  ros::NodeHandle nh;

  //Initialize parameters
  if (nh.getParam("imu/lp_alpha", lp_alpha_ ) && nh.getParam("imu/threshold", threshold_)
    && nh.getParam("odom/inactive_timeout", inactive_timeout_) ){
      nh.getParam("imu/lp_alpha", lp_alpha_);
      nh.getParam("imu/threshold", threshold_);
      nh.getParam("odom/inactive_time", inactive_timeout_);
      lp_beta_ = 1-lp_alpha_;
      ROS_INFO("Loaded imu_odometry parameters");
      ROS_INFO("alpha=%f, threshold=%f, inactive_t=%f",lp_alpha_, threshold_, inactive_timeout_);
  }
  else{ ROS_INFO("parameters not found!"); }

  imu_sub = nh.subscribe("imu",1, &ImuOdometryPublisher::SensorCallback, this);
  odom_pub = nh.advertise<nav_msgs::Odometry>("odom",50);

  last_update_time = ros::Time::now(); //TODO this isn't right.. where should it go?
  last_active_time = ros::Time::now();
}

ImuOdometryPublisher::~ImuOdometryPublisher() { }

void ImuOdometryPublisher::PublishOdometry(){
  //Publish stuff
  ros::Time current_time = ros::Time::now();

  //since all odometry is 6DOF we need a quaternion created from yaw
  geometry_msgs::Quaternion odom_quat = tf::createQuaternionMsgFromYaw(current_state.t);

  //Publish TF transform
  tf::TransformBroadcaster odom_broadcaster;
  geometry_msgs::TransformStamped odom_tf;
  odom_tf.header.stamp = current_time;
  odom_tf.header.frame_id = "odom";
  odom_tf.child_frame_id = "base_link";
  odom_tf.transform.translation.x = current_state.x;
  odom_tf.transform.translation.y = current_state.y;
  odom_tf.transform.rotation = odom_quat;

  //send the transform
  odom_broadcaster.sendTransform(odom_tf);


  // Publish odometry message
  nav_msgs::Odometry odom_msg;
  odom_msg.header.stamp = current_time;
  odom_msg.header.frame_id = "odom";

  //set the position
  odom_msg.pose.pose.position.x = current_state.x;
  odom_msg.pose.pose.position.y = current_state.y;
  odom_msg.pose.pose.orientation = odom_quat;

  //set the velocity
  odom_msg.child_frame_id = "base_link";
  odom_msg.twist.twist.linear.x = current_state.vx;
  odom_msg.twist.twist.linear.y = current_state.vy;
  odom_msg.twist.twist.angular.z = current_state.w;

  //publish the message
  odom_pub.publish(odom_msg);
}

void ImuOdometryPublisher::CalculateOdometry(){
  //Update current_state, publish

  ros::Time current_time = ros::Time::now();
  float h = current_time.toSec() - last_update_time.toSec();

  //update position, heading
  current_state.x += current_state.vx * h;
  current_state.y += current_state.vy * h;
  current_state.t += current_state.w  * h;

  //update velocities
  float ax = latest_imu.linear_acceleration[0]*cos(latest_imu.orientation[1])
    + latest_imu.linear_acceleration[2]*sin(latest_imu.orientation[1]);
  float ay = latest_imu.linear_acceleration[1]*cos(-latest_imu.orientation[0])
    + latest_imu.linear_acceleration[2]*sin(-latest_imu.orientation[0]);

  current_state.vx += ax * h;
  current_state.vy += ay * h;
  current_state.w = latest_imu.angular_velocity[2];
  //w_z is approximate; it doesn't take into account orientation

  ROS_INFO("roll= %f pitch=%f yaw=%f",latest_imu.orientation[0], latest_imu.orientation[1], latest_imu.orientation[2]);
  ROS_INFO("x= %f y=%f t= %f",current_state.x, current_state.y, current_state.t);
  ROS_INFO("vx= %f vy=%f w= %f",current_state.vx, current_state.vy, current_state.w);
  ROS_INFO("ax= %f ay=%f",ax, ay);

  //check ifrobot is (close to) not moving
  if((std::abs(ax)>threshold_) || (std::abs(ay)>threshold_) ||
    (std::abs(latest_imu.angular_velocity[2])>threshold_))
  {
    last_active_time = ros::Time::now();
    //ROS_INFO("updated last_active_time");
  }

  if((current_time.toSec() - last_active_time.toSec()) > inactive_timeout_){
    current_state.vx = current_state.vy = current_state.w = 0;
  }

  last_update_time = current_time;

  PublishOdometry();

}

void ImuOdometryPublisher::SensorCallback(const sensor_msgs::Imu::ConstPtr& msg) {
  // update latest_imu, call publish

  // orientation is originally in geometry_msgs/Quaternion
  tf::Quaternion quat;
  tf::quaternionMsgToTF(msg->orientation,quat);

  double roll, pitch, yaw;
  tf::Matrix3x3(quat).getRPY(roll, pitch, yaw);

  if(!LOWPASS){
    latest_imu.orientation[0] = roll;
    latest_imu.orientation[1] = pitch;
    latest_imu.orientation[2] = yaw;

    latest_imu.angular_velocity[0] = msg->angular_velocity.x;
    latest_imu.angular_velocity[1] = msg->angular_velocity.y;
    latest_imu.angular_velocity[2] = msg->angular_velocity.z;

    latest_imu.linear_acceleration[0] = msg->linear_acceleration.x;
    latest_imu.linear_acceleration[1] = msg->linear_acceleration.y;
    latest_imu.linear_acceleration[2] = msg->linear_acceleration.z;
  }
  else{
    //Low pass filter
    latest_imu.orientation[0] = lp_alpha_*roll  + lp_beta_*latest_imu.orientation[0];
    latest_imu.orientation[1] = lp_alpha_*pitch + lp_beta_*latest_imu.orientation[1];
    latest_imu.orientation[2] = lp_alpha_*yaw   + lp_beta_*latest_imu.orientation[2];

    latest_imu.angular_velocity[0] = lp_alpha_*(msg->angular_velocity.x) +
      lp_beta_*latest_imu.angular_velocity[0];
    latest_imu.angular_velocity[1] = lp_alpha_*(msg->angular_velocity.y) +
      lp_beta_*latest_imu.angular_velocity[1];
    latest_imu.angular_velocity[2] = lp_alpha_*(msg->angular_velocity.z) +
      lp_beta_*latest_imu.angular_velocity[2];

    latest_imu.linear_acceleration[0] = lp_alpha_*(msg->linear_acceleration.x) +
      lp_beta_*(latest_imu.linear_acceleration[0]);
    latest_imu.linear_acceleration[1] = lp_alpha_*(msg->linear_acceleration.y) +
      lp_beta_*(latest_imu.linear_acceleration[1]);
    latest_imu.linear_acceleration[2] = lp_alpha_*(msg->linear_acceleration.z) +
      lp_beta_*(latest_imu.linear_acceleration[2]);
  }

  //ROS_INFO("imu: %f",latest_imu.linear_acceleration[0]);

  CalculateOdometry();
}

int main(int argc, char** argv) {
  ros::init(argc, argv, "imu_odometry_publisher");

  ImuOdometryPublisher imu_odometry_publisher;

  ros::spin();

  return 0;
}
