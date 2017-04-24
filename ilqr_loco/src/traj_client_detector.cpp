

bool TrajClient::transform_laser_to_map(geometry_msgs::PointStamped &pos_laser_frame,
										                    geometry_msgs::PointStamped &pos_map_frame)
{
  try
  {
    tf::StampedTransform transform;
    tf_listener.waitForTransform("map", "laser", ros::Time::now(), ros::Duration(0.01));
    tf_listener.lookupTransform("map", "laser", ros::Time(0), transform);
    tf_listener.transformPoint("map", pos_laser_frame, pos_map_frame);
    return true;
  }
  catch (...)
  {
    ROS_INFO("Naive_obstacle_detector: Map to laser transform not available.");
    return false;
  }
}

void TrajClient::scanCb(const sensor_msgs::LaserScanConstPtr &msg)
{
  if(found_obs) return;

  // Look for obstacles in slice of scan
  int n_scans_close_enough = 0;

  for(int i=min_index; i<max_index; i++)
  {
    if (msg->ranges[i]<obstacle_thres)
    {
      n_scans_close_enough++;
      obs_dist = msg->ranges[i];
    }
  }

  float percent_scans_close = float(n_scans_close_enough)/float(max_index - min_index);

  geometry_msgs::PointStamped cluster_pos_localframe;
  geometry_msgs::PointStamped cluster_pos_mapframe;
  cluster_pos_localframe.header.frame_id = "laser";
  cluster_pos_mapframe.header.frame_id = "map";

  // Fill data, transform, send cluster_center_pos here
  if(percent_scans_close > percent_thres)
  {
    // ROS_INFO("Naive_obstacle_detector: Found obstacle.");
    found_obs = true;

    cluster_pos_localframe.point.x = obs_dist;
    cluster_pos_localframe.point.y = 0;
    cluster_pos_localframe.point.z = 0;

    if (transform_laser_to_map(cluster_pos_localframe, cluster_pos_mapframe))
    {
      cluster_pos_mapframe.point.y = 0; // hack to deal with heading uncertainty
      cc_pos.publish(cluster_pos_mapframe);
    }
  }
  else
  {
    cluster_pos_mapframe.point.x = 999;
    cluster_pos_mapframe.point.y = 0;
    cluster_pos_mapframe.point.z = 0;
    cc_pos.publish(cluster_pos_mapframe);
  }
}

void TrajClient::insert_fake_obs()
{
  found_obs = true;

  geometry_msgs::PointStamped cluster_pos_localframe;
  geometry_msgs::PointStamped cluster_pos_mapframe;
  cluster_pos_localframe.header.frame_id = "laser";
  cluster_pos_mapframe.header.frame_id = "map";

  cluster_pos_localframe.point.x = obstacle_thres;
  cluster_pos_localframe.point.y = 0;
  cluster_pos_localframe.point.z = 0;

  transform_laser_to_map(cluster_pos_localframe, cluster_pos_mapframe);
  cluster_pos_mapframe.point.y = 0; // hack to deal with heading uncertainty
  cc_pos.publish(cluster_pos_mapframe);
  // ROS_INFO("Naive_obstacle_detector: Published fake obstacle at %f, %f", cluster_pos_mapframe.point.x, cluster_pos_mapframe.point.y);
}
