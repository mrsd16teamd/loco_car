#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <std_msgs/Float32MultiArray.h>
#include <std_msgs/Int32MultiArray.h>
#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/LaserScan.h>
#include <tf/transform_listener.h>

tf::TransformListener *tran;
float obstacle_thres = 1.25; //[m]
ros::Publisher cc_pos;
float obs_dist = 0;
bool found_obs = false;

void scan_cb(const sensor_msgs::LaserScanConstPtr &msg)
{
  if(found_obs) return;

  int n_scans_close_enough = 0;
  int size_scan = msg->ranges.size();

  for(int i=0; i<size_scan; i++)
  {
    if (msg->ranges[i]<obstacle_thres)
    {
      n_scans_close_enough++;
      obs_dist = msg->ranges[i];
    }
  }

  float percent_scans_close = float(n_scans_close_enough)/float(size_scan);

  geometry_msgs::PointStamped cluster_pos_localframe;
  geometry_msgs::PointStamped cluster_pos_mapframe;

  //TODO fill, transform, sendcluster_center_pos here
  if(percent_scans_close>0.2)
  {
    cluster_pos_localframe.point.x = obs_dist;
    cluster_pos_localframe.point.y = 0;
    cluster_pos_localframe.point.z = 0;
    cluster_pos_localframe.header.frame_id = "laser";
    cluster_pos_mapframe.header.frame_id = "map";
    try
    {
      tf::StampedTransform transform;
      tran->waitForTransform("map", "laser", ros::Time::now(), ros::Duration(0.1));
      tran->lookupTransform("map", "laser", ros::Time(0), transform);
      tran->transformPoint("map", cluster_pos_localframe, cluster_pos_mapframe);
      cc_pos.publish(cluster_pos_mapframe);
      found_obs = true;
      ROS_INFO("Found obstacle, published at %f %f.", cluster_pos_mapframe.point.x, cluster_pos_mapframe.point.y);
    }
    catch(...)
    {
      ROS_INFO("Transform not available.");
    }
  }
  else
  {
    cluster_pos_mapframe.point.x = 999;
    cluster_pos_mapframe.point.y = 999;
    cluster_pos_mapframe.point.z = 999;
    cluster_pos_mapframe.header.frame_id = "map";
    cc_pos.publish(cluster_pos_mapframe);
  }
}

int main(int argc, char **argv) {
  ros::init(argc, argv, "obstacle_detector");
  ros::NodeHandle nh;

  ros::Subscriber sub = nh.subscribe("scan_front", 1, scan_cb);

  tf::TransformListener lr(ros::Duration(10));
  tran = &lr;

  cc_pos = nh.advertise<geometry_msgs::PointStamped>("cluster_center", 100); // clusterCenter1

  ros::spin();

  return 0;
}
