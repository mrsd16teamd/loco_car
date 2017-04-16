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

#include <pluginlib/class_list_macros.h>

namespace publishpcl_nodelet // @(namespace)
{
  class Scan2Cloud : public nodelet::Nodelet // @(NodeletClass)
  {
  public:
    Scan2Cloud();

  private:
    laser_geometry::LaserProjection projector_;
    tf::TransformListener tfListener_;

    ros::Publisher pcl_pub_;
    ros::Subscriber scan_sub_;
    ros::Timer timer_;

    double front_angle; // currently 90 degrees
    double new_angle_min;
    double new_angle_max;
    int min_index;
    int max_index;

    sensor_msgs::PointCloud2 cloud_;
    sensor_msgs::LaserScan scan_front_;


    virtual void onInit(){
      ros::NodeHandle nh = getNodeHandle();
      ros::NodeHandle& private_nh = getPrivateNodeHandle();

      // timer_ = nh.createTimer(ros::Duration(1.0), boost::bind(& Scan2Cloud::timerCb, this));
      scan_sub_ = nh.subscribe<sensor_msgs::LaserScan> ("scan", 1, &Scan2Cloud::scanCallback, this);
      pcl_pub_ = private_nh.advertise<sensor_msgs::PointCloud2> ("scan_cloud", 1, false);

      try{
        nh.getParam("scan_clip_angle", front_angle);
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
    };

    void timerCb(const ros::TimerEvent& event){
    // Using timers is the preferred 'ROS way' to manual threading
    NODELET_INFO_STREAM("The time is now " << event.current_real);
    }

    void scanCallback(const sensor_msgs::LaserScan::ConstPtr& scan) {
      scan_front_ = *scan;
      scan_front_.angle_min = new_angle_min;
      scan_front_.angle_max = new_angle_max;
      std::vector<float> new_data(&scan->ranges[min_index],&scan->ranges[max_index]);
      scan_front_.ranges = new_data;

      projector_.projectLaser(scan_front_, cloud_);
      pcl_pub_.publish(cloud_);
    }
  };
} // namespace publishpcl

PLUGINLIB_DECLARE_CLASS(publishpcl_nodelet, Scan2Cloud, publishpcl_nodelet::Scan2Cloud, nodelet::Nodelet);
