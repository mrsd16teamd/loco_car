#include "ros/ros.h"
#include "sensor_msgs/Imu.h"
#include <tf/transform_broadcaster.h>
#include <nav_msgs/Odometry.h>

struct imuData {
  double ax, ay, w;
} imu_data;

void imuCallback(const sensor_msgs::Imu::ConstPtr& msg){
  //ROS_INFO("Received IMU data.");

  //TODO add transform from imu_frame to base_link
  //define base_link as midpoint of rect defined by wheels?
  //needs to make static_transform_publisher for this and laser to base_link

  imu_data.ax = msg->linear_acceleration.x;
  imu_data.ay = msg->linear_acceleration.y;
  imu_data.w  = msg->angular_velocity.z;

}

int main(int argc, char** argv){
  ros::init(argc, argv, "odometry_publisher");

  ros::NodeHandle n;
  ros::Publisher odom_pub = n.advertise<nav_msgs::Odometry>("odom",50);
  tf::TransformBroadcaster odom_broadcaster;
  ros::Subscriber imu_sub = n.subscribe("imu",1000,imuCallback);

  double x = 0.0; double y = 0.0; double th = 0.0;
  double vx = 0.0; double vy = 0.0; double vth = 0.0;

  ros::Time current_time, last_time;
  current_time = ros::Time::now(); last_time = ros::Time::now();
  double no_move_time = current_time.toSec();

  ros::Rate r(50.0); //TODO make this faster later [Hz]
  while(n.ok()){

    ros::spinOnce();               // check for incoming messages

    // Assuming x is forward, y is right, z is up. acc[m/s^2], w[rad/s]
    //integrate new sensor readings

    current_time = ros::Time::now();

    //compute odometry in a typical way given the velocities of the robot
    double dt = (current_time - last_time).toSec();

    vx += imu_data.ax*dt;
    vy += imu_data.ay*dt;
    vth = imu_data.w;

    // Before integrating odometry, check if robot has been inactive for a while
    if( (imu_data.ax==0.0) && (imu_data.ay==0.0) && (imu_data.w==0.0) ) {
	no_move_time = no_move_time;
    }
    else{
	no_move_time = current_time.toSec();
    }
    if((current_time.toSec() - no_move_time)>0.3){
	vx = 0;
	vy = 0;
	vth = 0;
    }

    //ROS_INFO("current time = %f, no_move_time = %f", current_time.toSec(), no_move_time);

    double delta_x = (vx * cos(th) - vy * sin(th)) * dt;
    double delta_y = (vx * sin(th) + vy * cos(th)) * dt;
    double delta_th = vth * dt;

    x += delta_x;
    y += delta_y;
    th += delta_th;

    ROS_INFO("x=%.2f, y=%.2f, th=%.2f, vx=%.2f, vy=%.2f, vth=%.2f", x,y,th,vx,vy,vth);

    //since all odometry is 6DOF we'll need a quaternion created from yaw
    geometry_msgs::Quaternion odom_quat = tf::createQuaternionMsgFromYaw(th);

    //first, we'll publish the transform over tf
    geometry_msgs::TransformStamped odom_trans;
    odom_trans.header.stamp = current_time;
    odom_trans.header.frame_id = "odom";
    odom_trans.child_frame_id = "base_link";

    odom_trans.transform.translation.x = x;
    odom_trans.transform.translation.y = y;
    odom_trans.transform.translation.z = 0.0;
    odom_trans.transform.rotation = odom_quat;

    //send the transform
    odom_broadcaster.sendTransform(odom_trans);

    //next, we'll publish the odometry message over ROS
    nav_msgs::Odometry odom;
    odom.header.stamp = current_time;
    odom.header.frame_id = "odom";

    //set the position
    odom.pose.pose.position.x = x;
    odom.pose.pose.position.y = y;
    odom.pose.pose.position.z = 0.0;
    odom.pose.pose.orientation = odom_quat;

    //set the velocity
    odom.child_frame_id = "base_link";
    odom.twist.twist.linear.x = vx;
    odom.twist.twist.linear.y = vy;
    odom.twist.twist.angular.z = vth;

    //publish the message
    odom_pub.publish(odom);

    last_time = current_time;
    r.sleep();
  }
}

// class State {
// public:
//   double x,y,th,vx,vy,w;
//   State() : x(0), y(0), th(0), vx(0), vy(0), w(0) {}
// } state;
// right way to do it with callback is here: http://answers.ros.org/question/59725/publishing-to-a-topic-via-subscriber-callback-function/
