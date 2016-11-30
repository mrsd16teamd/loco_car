#include "ros/ros.h"
#include "tf/transform_datatypes.h"
#include "move_base_msgs/MoveBaseActionGoal.h"
#include <tf/transform_listener.h>
#include "std_msgs/Float64.h"
#include <cmath>
#include <string>

float last_goal_x = 0.0;
float last_goal_y = 0.0;
ros::Subscriber goal_sub;
ros::Publisher dist_pub;

void GoalCallback(const move_base_msgs::MoveBaseActionGoal::ConstPtr& msg){
  ROS_INFO("dist2goal_node: received goal");
  last_goal_x = msg->goal.target_pose.pose.position.x;
  last_goal_y = msg->goal.target_pose.pose.position.y;
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "dist2goal_node");

  ros::NodeHandle nh;

  goal_sub = nh.subscribe("move_base/goal", 1, GoalCallback);
  dist_pub = nh.advertise<std_msgs::Float64>("dist2goal_metric", 50);

  tf::TransformListener listener;

  ros::Rate r(1);

  while(ros::ok()){
    ros::spinOnce();

    float current_x, current_y;

    tf::StampedTransform transform;

    try{
       listener.lookupTransform("map", "base_link", ros::Time(0), transform);
       current_x = transform.getOrigin().x();
       current_y = transform.getOrigin().y();

       //ROS_INFO("Currently at: %f , %f ", current_x, current_y);
    }
    catch (tf::TransformException ex){
      ROS_ERROR("dist2goal_node: map to base_link transform not found on this cycle");
    }

    float distance = sqrt( pow((last_goal_x-current_x),2) + pow((last_goal_y-current_y),2) );

    std_msgs::Float64 dist_msg;
    dist_msg.data = distance;
    dist_pub.publish(dist_msg);

    r.sleep();
  }

  ros::spin();

  return 0;
}
