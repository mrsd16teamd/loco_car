#include <ros/ros.h>
#include <nodelet/nodelet.h>
#include <laser_assembler/AssembleScans.h>
#include "laser_geometry.h"
#include <algorithm>
#include <ros/assert.h>
#include <tf2/LinearMath/Transform.h>
#include <tf/transform_listener.h>
#include <sensor_msgs/PointCloud.h>
#include <tf/message_filter.h>
#include <message_filters/subscriber.h>
#include <vector>
#include <math.h>

#include <iostream>

#include <pluginlib/class_list_macros.h>

namespace publishpcl_nodelet // @(namespace)
{
  class Scan2Cloud : public nodelet::Nodelet // @(NodeletClass)
  {
  public:
    Scan2Cloud() = default;

  private:
    laser_geometry::LaserProjection projector_;
    tf::TransformListener tfListener_;

    ros::Publisher pcl_pub_;
    ros::Subscriber scan_sub_;
    ros::Timer timer_;

    double front_angle_; // currently 90 degrees
    double new_angle_min_;
    double new_angle_max_;
    int min_index_;
    int max_index_;

    sensor_msgs::PointCloud2 cloud_;
    sensor_msgs::LaserScan scan_front_;


    virtual void onInit(){
      ros::NodeHandle nh = getNodeHandle();
      ros::NodeHandle& private_nh = getPrivateNodeHandle();

      scan_sub_ = nh.subscribe<sensor_msgs::LaserScan> ("scan", 1, &Scan2Cloud::scanCallback, this);
      pcl_pub_ = private_nh.advertise<sensor_msgs::PointCloud2> ("scan_cloud", 1, false);

      if (nh.hasParam("scan_clip_angle")){
        nh.getParam("scan_clip_angle", front_angle_);
      }
      else{
        NODELET_INFO_STREAM("Need param scan_clip_angle!");
        ros::shutdown();
      }

      new_angle_min_ = -front_angle_/2;
      new_angle_max_ = front_angle_/2;

      double old_angle = 2.35619449615;
      double increment = 0.00436332309619;

      min_index_ = floor((old_angle-new_angle_max_) / increment);
      max_index_ = 1080 - min_index_;
    //   std::cout << "min_index_ = " << min_index_ << std::endl;
    //   std::cout << "max_index_ = " << max_index_ << std::endl;
    };

    void scanCallback(const sensor_msgs::LaserScan::ConstPtr scan) {
      scan_front_ = *scan;
      scan_front_.angle_min = new_angle_min_;
      scan_front_.angle_max = new_angle_max_;
      scan_front_.ranges.assign(scan->ranges.begin()+min_index_, scan->ranges.begin()+max_index_);
      projector_.projectLaser(scan_front_, cloud_);
      pcl_pub_.publish(cloud_);
    }
  };
} // namespace publishpcl_nodelet

PLUGINLIB_DECLARE_CLASS(publishpcl_nodelet, Scan2Cloud, publishpcl_nodelet::Scan2Cloud, nodelet::Nodelet);
