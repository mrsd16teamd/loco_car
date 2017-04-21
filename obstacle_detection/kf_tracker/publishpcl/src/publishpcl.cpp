#include <ros/ros.h>
#include <laser_assembler/AssembleScans.h>
#include "laser_geometry/laser_geometry.h"
#include <algorithm>
#include <ros/assert.h>
#include <tf2/LinearMath/Transform.h>
#include <tf/transform_listener.h>
#include <sensor_msgs/PointCloud.h>
#include <tf/message_filter.h>
#include <message_filters/subscriber.h>
#include <vector>
#include <math.h>
#include <ctime>

class Scan2Cloud {
public:
  Scan2Cloud();
  void scanCallback(const sensor_msgs::LaserScan::ConstPtr& scan);
private:
  ros::NodeHandle node_;
  laser_geometry::LaserProjection projector_;
  tf::TransformListener tfListener_;

  ros::Publisher point_cloud_publisher_;
  ros::Subscriber scan_sub_;

  sensor_msgs::PointCloud2 cloud_;
  double front_angle; // currently 90 degrees
  double new_angle_min;
  double new_angle_max;
  int min_index;
  int max_index;

  sensor_msgs::LaserScan scan_front;
};


Scan2Cloud::Scan2Cloud()
{
  scan_sub_ = node_.subscribe<sensor_msgs::LaserScan> ("scan", 1, &Scan2Cloud::scanCallback, this);
  point_cloud_publisher_ = node_.advertise<sensor_msgs::PointCloud2> ("scan_cloud", 1, false);

  try{
    node_.getParam("scan_clip_angle", front_angle);
  }
  catch(...){
    ROS_ERROR("Need param scan_clip_angle!");
    ros::shutdown();
  }

  new_angle_min = -front_angle/2;
  new_angle_max = front_angle/2;

  double old_angle = 2.35619449615;
  double increment = 0.00436332309619;

  min_index = floor((old_angle-new_angle_max) / increment);
  max_index = 1080 - min_index;


}



void Scan2Cloud::scanCallback(const sensor_msgs::LaserScan::ConstPtr& scan){

  clock_t start = clock();
  //pre-process scans, chopping to front sub-section
  scan_front = *scan;
  scan_front.angle_min = new_angle_min;
  scan_front.angle_max = new_angle_max;
  std::vector<float> new_data(&scan->ranges[min_index],&scan->ranges[max_index]);
  scan_front.ranges = new_data;

  projector_.projectLaser(scan_front, cloud_);
  point_cloud_publisher_.publish(cloud_);


    std::cout<<"publishing time"<<std::endl;
    std::cout<<seconds<<std::endl;
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "scan2cloud_publisher");

  Scan2Cloud scan2cloud_publisher;

  ros::spin();

  return 0;
}
