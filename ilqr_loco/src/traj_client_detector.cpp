#include "traj_client.h"

// TODO put this all in a NaiveObstacleDetector class
// put that class within TrajClient - make sure it has access to TrajClient's members

void TrajClient::InitDetector()
{
  double new_angle_min = -scan_front_angle_/2;
  double new_angle_max = scan_front_angle_/2;
  double old_angle = 2.35619449615;
  double increment = 0.00436332309619;

  scan_min_index_ = floor((old_angle-new_angle_max) / increment);
  scan_max_index_ = 1080 - scan_min_index_;
}

bool TrajClient::TransformLaserToMap(geometry_msgs::PointStamped &pos_laser_frame,
										                    geometry_msgs::PointStamped &pos_map_frame)
{
  try
  {
    tf::StampedTransform transform;
    tf_listener_.waitForTransform("map", "laser", ros::Time::now(), ros::Duration(0.01));
    tf_listener_.lookupTransform("map", "laser", ros::Time(0), transform);
    tf_listener_.transformPoint("map", pos_laser_frame, pos_map_frame);
    return true;
  }
  catch (...)
  {
    ROS_INFO("Obstacle detector: Map to laser transform not available.");
    return false;
  }
}

void TrajClient::scanCb(const sensor_msgs::LaserScanConstPtr &msg)
{
  // ROS_INFO("start obstacle detector");

  if(found_obstacle_) return;

  // Look for obstacles in slice of scan
  int n_scans_close_enough = 0;
  double obs_dist;

  for(int i=scan_min_index_; i<scan_max_index_; i++)
  {
    if (msg->ranges[i]<obs_dist_thres_)
    {
      n_scans_close_enough++;
      obs_dist = msg->ranges[i];
    }
  }

  float percent_scans_close = float(n_scans_close_enough)/float(scan_max_index_ - scan_min_index_);

  // Fill data, transform, send cluster_center_pos here
  if(percent_scans_close > obs_percent_thres_)
  {
    ROS_INFO("Obstacle detector: Found obstacle.");
    geometry_msgs::PointStamped obs_pos_localframe;
    geometry_msgs::PointStamped obs_pos_mapframe;
    obs_pos_localframe.header.frame_id = "laser";
    obs_pos_mapframe.header.frame_id = "map";

    double yaw = tf::getYaw(cur_state_.pose.pose.orientation);
    obs_pos_localframe.point.x = obs_dist * cos(yaw);
    obs_pos_localframe.point.y = 0;

    if (TransformLaserToMap(obs_pos_localframe, obs_pos_mapframe))
    {
      obs_pos_.x = obs_pos_mapframe.point.x;
      obs_pos_.y = 0; // hack to deal with heading uncertainty
      obs_pos_mapframe.point.y = 0; // hack to deal with heading uncertainty
      found_obstacle_ = true;
      obs_pos_pub_.publish(obs_pos_mapframe);
  	  ROS_INFO("Transformed obs pos.");

      if (!reacted_to_obstacle_)
        ReactToObstacle();

        ros::spinOnce();
    }
  }
}

void TrajClient::InsertFakeObs()
{
  found_obstacle_ = true;

  geometry_msgs::PointStamped obs_pos_localframe;
  geometry_msgs::PointStamped obs_pos_mapframe;
  obs_pos_localframe.header.frame_id = "laser";
  obs_pos_mapframe.header.frame_id = "map";

  obs_pos_localframe.point.x = obs_dist_thres_;
  obs_pos_localframe.point.y = 0;

  TransformLaserToMap(obs_pos_localframe, obs_pos_mapframe);
  obs_pos_.x = obs_pos_mapframe.point.x;
  obs_pos_.y = 0; // hack to deal with heading uncertainty
  obs_pos_mapframe.point.y = 0; // hack to deal with heading uncertainty
  obs_pos_pub_.publish(obs_pos_mapframe);

  if (!reacted_to_obstacle_)
    ReactToObstacle();

  ros::spinOnce();
  // ROS_INFO("Obstacle_detector: Published fake obstacle at %f, %f", cluster_pos_mapframe.point.x, cluster_pos_mapframe.point.y);
}

void TrajClient::ResetObstacle()
{
  reacted_to_obstacle_ = false;
  found_obstacle_ = false;

  geometry_msgs::PointStamped reset_obs_pos;
  reset_obs_pos.header.frame_id = "map";
  reset_obs_pos.point.x = 999;
  reset_obs_pos.point.y = 0;
  obs_pos_pub_.publish(reset_obs_pos);
  ros::spinOnce();
}
