#include <iostream>
#include <string.h>
#include <fstream>
#include <algorithm>
#include <iterator>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/video/video.hpp>
#include "opencv2/video/tracking.hpp"
#include <ros/ros.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include "pcl_ros/point_cloud.h"
#include <geometry_msgs/Point.h>
#include <std_msgs/Float32MultiArray.h>
#include <std_msgs/Int32MultiArray.h>

#include <sensor_msgs/PointCloud2.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/geometry.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/common/centroid.h>
 #include <tf/transform_listener.h>
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>
#include <limits>
#include <utility>
#include <pcl/registration/correspondence_estimation.h>
#include "Kf.h"
#include <utility>
#include <ctime>
float obstacle_thres;
int mode;
static bool found_obs = false;
ros::Publisher objID_pub;
ros::Publisher cc_pos;
//ros::Publisher markerPub1;

tf::TransformListener* tran;

void cluster_extraction (const sensor_msgs::PointCloud2ConstPtr& input, std::vector<pcl::PointIndices>& cluster_indices);

int DEBUGMODE = 0;



void cloud_cb (const sensor_msgs::PointCloud2ConstPtr& input)

{ 
 if(mode)
{
    if(found_obs) return;
}
    //initialize the clustercenter
    std_msgs::Float32MultiArray cluster_center;
    float xcoordinate(9.0f);
    float ycoordinate(0.0f);
    float zcoordinate(0.0f);
    static bool obstaclepresent(0);
    pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud (new pcl::PointCloud<pcl::PointXYZ>);
    //initialize the clustercenter
    std_msgs::Float32MultiArray clustermapframe;
    // initialize the clusterindices to be passed to cluster
    std::vector<pcl::PointIndices> cluster_indices;
    // do clustering 


    pcl::fromROSMsg (*input, *input_cloud);
    clock_t start = clock();
    cluster_extraction (input, cluster_indices);



/*
    if ( DEBUGMODE ==1 )
    {    std:cout<<cluster_indices.size()<<std::endl;

	}*/
    std::vector<pcl::PointIndices>::const_iterator it;
    std::vector<int>::const_iterator pit;
    // Vector of cluster pointclouds
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr > cluster_vec;

    // Cluster centroids
    std::vector<pcl::PointXYZ> clusterCentroids;

    cluster_center.data.push_back(0);
    cluster_center.data.push_back(0);
    cluster_center.data.push_back(0); 

    for(it = cluster_indices.begin(); it != cluster_indices.end(); ++it)
    {
        float x=0.0; float y=0.0;
        int numPts=0;
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_cluster (new pcl::PointCloud<pcl::PointXYZ>);
        for(pit = it->indices.begin(); pit != it->indices.end(); pit++)
        {
            cloud_cluster->points.push_back(input_cloud->points[*pit]);
            x+=input_cloud->points[*pit].x;
            y+=input_cloud->points[*pit].y;
            numPts++;
	    

                  
        }

        pcl::PointXYZ centroid;
        centroid.x=x/numPts;
        centroid.y=y/numPts;
        centroid.z=0.0;
 
        if (centroid.x < obstacle_thres && centroid.x > -1*obstacle_thres && centroid.y > -1*obstacle_thres && centroid.y < obstacle_thres && centroid.x != 0 && centroid.y != 0 )
        {
             
            xcoordinate = centroid.x;
            ycoordinate = centroid.y;
	    zcoordinate = centroid.z;
	    obstaclepresent = 1;

	    
        }


	
    };

    Eigen::Vector4f obstaclepoint;
    geometry_msgs::PointStamped laserframe;
    geometry_msgs::PointStamped mapframe;
    laserframe.header.frame_id = "laser";  

    mapframe.header.frame_id  = "map";
   
 

    tf::StampedTransform transform;

   // setting up laserframe coordinates from cluster centers
   laserframe.point.x = xcoordinate;
   laserframe.point.y = ycoordinate;
   laserframe.point.z = zcoordinate;

    // converting laserframe to mapframe 
    try
    {
	    tran->waitForTransform("map","laser",ros::Time::now(), ros::Duration(0.01));
            tran->lookupTransform("map","laser", ros::Time(0), transform);
            tran->transformPoint("map",laserframe, mapframe);
    }
    catch (tf::TransformException& ex) 
    {
    }




   if ( DEBUGMODE == 1)
    {

    
        std::cout<<"worldcoordiante"<<std::endl;
        std::cout<<mapframe.point.x<<std::endl;
        std::cout<<mapframe.point.y<<std::endl;
        std::cout<<mapframe.point.z<<std::endl;

        std::cout<<"laserframe"<<std::endl;
        std::cout<<laserframe.point.x<<std::endl;
        std::cout<<laserframe.point.y<<std::endl;

    }

    //setting up the marker 
/*
    visualization_msgs::MarkerArray clusterMarkers1;


        visualization_msgs::Marker m1;

        m1.id=0;
        m1.type=visualization_msgs::Marker::CUBE;
        m1.header.frame_id="/map";
        m1.scale.x=0.3;         m1.scale.y=0.3;         m1.scale.z=0.3;
        m1.action=visualization_msgs::Marker::ADD;
        m1.color.a=1.0;
        m1.color.r=1;
        m1.color.g=0;
        m1.color.b=0;
*/
     
// checking if the cluster was found near the give radious of robot 

/*

   if (  abs(laserframe.point.x)  != 9.0  &&  abs(laserframe.point.y)  != 0.0)
   {





       clustermapframe.data.push_back(mapframe.point.x);
       clustermapframe.data.push_back(mapframe.point.y);
       clustermapframe.data.push_back(mapframe.point.z);
    }


   else 
    {
       mapframe.point.x = xcoordinate;
       mapframe.point.y  = ycoordinate;
       mapframe.point.z  = zcoordinate;
    }


       m1.pose.position.x=mapframe.point.x;
       m1.pose.position.y=mapframe.point.y;
       m1.pose.position.z=mapframe.point.z;

*/

   if( obstaclepresent == 0 )
   {    
       mapframe.point.x = 999 ;
       mapframe.point.y  = 999 ;
       mapframe.point.z  = 0;
      // cc_pos.publish(mapframe);

   } 

  if( obstaclepresent == 1 )
   {    
        found_obs = true;
//	clusterMarkers1.markers.push_back(m1);  
	cc_pos.publish(mapframe);
//	markerPub1.publish(clusterMarkers1);


   } 

      
  obstaclepresent = 0;
   

    clock_t end = clock();
    float seconds = (float)(end - start) / CLOCKS_PER_SEC;

 //   std::cout<<"clustering time"<<std::endl;
 //   std::cout<<seconds<<std::endl;
    



};


void mode_cb(const geometry_msgs::Point &msg)
{
  if(msg.x == 8){
    ROS_INFO("Looking for obstacle again.");
    found_obs = false;
  }
}


int main(int argc, char** argv)
{




 
    ros::init (argc,argv,"naive_detector");
    ros::NodeHandle nh;
    std::cout<<"About to setup callback tracker\n";
    tf::TransformListener lr(ros::Duration(10));
    tran=&lr;

  try{
    nh.getParam("naive_obstacle_dist_thres", obstacle_thres);
   nh.getParam("setforcontinuousdetection", mode);
  }
  catch(...){
    ROS_ERROR("Need param obstacle_thres");
    ros::shutdown();
  }


    ros::Subscriber sub = nh.subscribe ("scan_cloud", 1, cloud_cb);
  ros::Subscriber mode_sub = nh.subscribe("client_command", 1, mode_cb);
    cc_pos=nh.advertise<geometry_msgs::PointStamped>("cluster_center",100);//clusterCenter1

   // markerPub1= nh.advertise<visualization_msgs::MarkerArray> ("viz1",1);


    ros::spin();

}

