#include "kf_tracker.h"

using namespace std;
using namespace cv;

// tf::TransformListener listener;
// ros::Publisher objID_pub;

tf::TransformListener *tran;

// KF init
int stateDim = 4; // [x,y,v_x,v_y]//,w,h]
int measDim = 2;  // [z_x,z_y,z_w,z_h]
int ctrlDim = 0;

int n_clusters = 1;
std::vector<cv::KalmanFilter> KF_vec;
cv::KalmanFilter KF_start(stateDim, measDim, ctrlDim, CV_32F);

ros::Publisher cc_pos;
// ros::Publisher markerPub;
ros::Publisher markerPub1;

std::vector<geometry_msgs::Point> prevClusterCenters;

cv::Mat state(stateDim, 1, CV_32F);
cv::Mat_<float> measurement(2, 1);

std::vector<int> objID; // Output of the data association using KF

bool firstFrame = true;

void KFT(const std_msgs::Float32MultiArray ccs) {
  // First predict, to update the internal statePre variable
  std::vector<cv::Mat> pred;
  for (int i = 0; i < n_clusters; i++) {
    pred.push_back(KF_vec[i].predict());
  }

  // Get measurements
  // Extract the position of the clusters forom the multiArray. To check if the
  // data
  // coming in, check the .z (every third) coordinate and that will be 0.0
  std::vector<geometry_msgs::Point> clusterCenters; // clusterCenters

  int i = 0;
  for (std::vector<float>::const_iterator it = ccs.data.begin();
       it != ccs.data.end(); it += 3) {
    geometry_msgs::Point pt;
    pt.x = *it;
    pt.y = *(it + 1);
    pt.z = *(it + 2);

    clusterCenters.push_back(pt);
  }

  std::vector<geometry_msgs::Point> KFpredictions;
  i = 0;
  for (std::vector<cv::Mat>::iterator it = pred.begin(); it != pred.end();
       it++) {
    geometry_msgs::Point pt;
    pt.x = (*it).at<float>(0);
    pt.y = (*it).at<float>(1);
    pt.z = (*it).at<float>(2);

    KFpredictions.push_back(pt);
  }

  // Find the cluster that is more probable to be belonging to a given KF.
  objID.clear(); // Clear the objID vector
  objID.resize(n_clusters); // Allocate default elements so that [i] doesnt
                            // segfault. Should be done better

  // Copy clusterCentres for modifying it and preventing multiple assignments of
  // the same ID
  std::vector<geometry_msgs::Point> copyOfClusterCenters(clusterCenters);

  std::vector<std::vector<float> > distMat;

  for (int filterN = 0; filterN < n_clusters; filterN++) {
    std::vector<float> distVec;
    for (int n = 0; n < n_clusters; n++) {
      distVec.push_back(
          euclidean_distance(KFpredictions[filterN], copyOfClusterCenters[n]));
    }
    distMat.push_back(distVec);
  }

  for (int clusterCount = 0; clusterCount < n_clusters; clusterCount++) {
    // 1. Find min(distMax)==> (i,j);
    std::pair<int, int> minIndex(findIndexOfMin(distMat));

    // 2. objID[i]=clusterCenters[j]; counter++
    objID[minIndex.first] = minIndex.second;

    // 3. distMat[i,:]=10000; distMat[:,j]=10000
    distMat[minIndex.first] = std::vector<float>(
        n_clusters, 10000.0); // Set the row to a high number.
    for (int row = 0; row < distMat.size();
         row++) // set the column to a high number
    {
      distMat[row][minIndex.second] = 10000.0;
    }
  }

  visualization_msgs::MarkerArray clusterMarkers;
  for (int i = 0; i < n_clusters; i++) {
    visualization_msgs::Marker m;

    m.id = i;
    m.type = visualization_msgs::Marker::CUBE;
    m.header.frame_id = "/laser";
    m.scale.x = 0.3;
    m.scale.y = 0.3;
    m.scale.z = 0.3;
    m.action = visualization_msgs::Marker::ADD;
    m.color.a = 1.0;
    m.color.r = i % 2 ? 1 : 0;
    m.color.g = i % 3 ? 1 : 0;
    m.color.b = i % 4 ? 1 : 0;

    // geometry_msgs::Point clusterC(clusterCenters.at(objID[i]));
    geometry_msgs::Point clusterC(KFpredictions[i]);
    m.pose.position.x = clusterC.x;
    m.pose.position.y = clusterC.y;
    m.pose.position.z = clusterC.z;

    clusterMarkers.markers.push_back(m);
  }

  prevClusterCenters = clusterCenters;

  // markerPub.publish(clusterMarkers);

  std_msgs::Int32MultiArray obj_id;
  for (std::vector<int>::iterator it = objID.begin(); it != objID.end(); it++)
    obj_id.data.push_back(*it);

  // Publish the object IDs
  // objID_pub.publish(obj_id);

  // convert clusterCenters from geometry_msgs::Point to floats
  std::vector<std::vector<float> > cc;
  for (int i = 0; i < n_clusters; i++) {
    vector<float> pt;
    pt.push_back(clusterCenters[objID[i]].x);
    pt.push_back(clusterCenters[objID[i]].y);
    pt.push_back(clusterCenters[objID[i]].z);

    cc.push_back(pt);
  }

  // The update phase
  std::vector<cv::Mat> measmat_vec(n_clusters);
  for (int i=0; i<n_clusters; i++)
  {
    cv::Mat measmat = cv::Mat(2, 1, CV_32F, &cc[i]);
    measmat_vec[i] = measmat;
  }
}

// If this is the first frame, initialize kalman filters for the clustered objects
void init_kf(const sensor_msgs::PointCloud2ConstPtr &input)
{
  // Initialize n Kalman Filters; Assuming n max objects in the dataset.
  // Could be made generic by creating a Kalman Filter only when a new object is
  // detected

  for (int i = 0; i < n_clusters; i++) {
    KF_vec[i].transitionMatrix = (Mat_<float>(4, 4) << dx, 0, 1, 0, 0, dy, 0, 1,
                                  0, 0, dvx, 0, 0, 0, 0, dvy);
    cv::setIdentity(KF_vec[i].measurementMatrix);
    setIdentity(KF_vec[i].processNoiseCov, Scalar::all(sigmaP));

    // Meas noise cov matrix R
    cv::setIdentity(KF_vec[i].measurementNoiseCov, cv::Scalar(sigmaQ)); // 1e-1
  }

  // Process the point cloud
  pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud(
      new pcl::PointCloud<pcl::PointXYZ>);
  pcl::PointCloud<pcl::PointXYZ>::Ptr clustered_cloud(
      new pcl::PointCloud<pcl::PointXYZ>);

  /* Creating the KdTree from input point cloud*/
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(
      new pcl::search::KdTree<pcl::PointXYZ>);

  pcl::fromROSMsg(*input, *input_cloud);

  tree->setInputCloud(input_cloud);

  std::vector<pcl::PointIndices> cluster_indices;
  pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
  ec.setClusterTolerance(0.04);
  ec.setMinClusterSize(50);
  ec.setMaxClusterSize(200);
  ec.setSearchMethod(tree);
  ec.setInputCloud(input_cloud);

  /* Extract the clusters out of pc and save indices in cluster_indices.*/
  ec.extract(cluster_indices);

  std::vector<pcl::PointIndices>::const_iterator it;
  std::vector<int>::const_iterator pit;
  // Vector of cluster pointclouds
  std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> cluster_vec;
  // Cluster centroids
  std::vector<pcl::PointXYZ> clusterCentroids;

  for (it = cluster_indices.begin(); it != cluster_indices.end(); ++it) {

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_cluster(
        new pcl::PointCloud<pcl::PointXYZ>);
    float x = 0.0;
    float y = 0.0;
    int numPts = 0;
    for (pit = it->indices.begin(); pit != it->indices.end(); pit++) {

      cloud_cluster->points.push_back(input_cloud->points[*pit]);
      x += input_cloud->points[*pit].x;
      y += input_cloud->points[*pit].y;
      numPts++;
    }

    pcl::PointXYZ centroid;
    centroid.x = x / numPts;
    centroid.y = y / numPts;
    centroid.z = 0.0;

    cluster_vec.push_back(cloud_cluster);

    // Get the centroid of the cluster
    clusterCentroids.push_back(centroid);
  }

  // Ensure at least 6 clusters exist to publish (later clusters may be empty)
  while (cluster_vec.size() < n_clusters) {
    pcl::PointCloud<pcl::PointXYZ>::Ptr empty_cluster(
        new pcl::PointCloud<pcl::PointXYZ>);
    empty_cluster->points.push_back(pcl::PointXYZ(0, 0, 0));
    cluster_vec.push_back(empty_cluster);
  }

  while (clusterCentroids.size() < n_clusters) {
    pcl::PointXYZ centroid;
    centroid.x = 0.0;
    centroid.y = 0.0;
    centroid.z = 0.0;

    clusterCentroids.push_back(centroid);
  }

  // Set initial state
  for (int i = 0; i < n_clusters; i++) {
    KF_vec[i].statePre.at<float>(0) = clusterCentroids.at(i).x;
    KF_vec[i].statePre.at<float>(1) = clusterCentroids.at(i).y;
    KF_vec[i].statePre.at<float>(2) = 0; // initial v_x
    KF_vec[i].statePre.at<float>(3) = 0; // initial v_y
  }

  firstFrame = false;

  for (int i = 0; i < n_clusters; i++) {
    geometry_msgs::Point pt;
    pt.x = clusterCentroids.at(i).x;
    pt.y = clusterCentroids.at(i).y;
    prevClusterCenters.push_back(pt);
  }
}

void cloud_cb(const sensor_msgs::PointCloud2ConstPtr &input) {
  if (firstFrame) {
    init_kf(input);
  }
  else
  {
    pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud(
        new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr clustered_cloud(
        new pcl::PointCloud<pcl::PointXYZ>);

    /* Creating the KdTree from input point cloud*/
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(
        new pcl::search::KdTree<pcl::PointXYZ>);

    pcl::fromROSMsg(*input, *input_cloud);

    tree->setInputCloud(input_cloud);

    /* Here we are creating a vector of PointIndices, which contains the actual
    * index
    * information in a vector<int>. The indices of each detected cluster are
    * saved here.
    * Cluster_indices is a vector containing one instance of PointIndices for
    * each detected
    * cluster. Cluster_indices[0] contain all indices of the first cluster in
    * input point cloud.
    */
    std::vector<pcl::PointIndices> cluster_indices;
    pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
    ec.setClusterTolerance(0.04);
    ec.setMinClusterSize(50);
    ec.setMaxClusterSize(200);
    ec.setSearchMethod(tree);
    ec.setInputCloud(input_cloud);

    /* Extract the clusters out of pc and save indices in cluster_indices.*/
    ec.extract(cluster_indices);

    /* To separate each cluster out of the vector<PointIndices> we have to
    * iterate through cluster_indices, create a new PointCloud for each
    * entry and write all points of the current cluster in the PointCloud.
    */

    std::vector<pcl::PointIndices>::const_iterator it;
    std::vector<int>::const_iterator pit;
    // Vector of cluster pointclouds
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> cluster_vec;

    // Cluster centroids
    std::vector<pcl::PointXYZ> clusterCentroids;

    for (it = cluster_indices.begin(); it != cluster_indices.end(); ++it) {
      float x = 0.0;
      float y = 0.0;
      int numPts = 0;
      pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_cluster(
          new pcl::PointCloud<pcl::PointXYZ>);
      for (pit = it->indices.begin(); pit != it->indices.end(); pit++) {
        cloud_cluster->points.push_back(input_cloud->points[*pit]);

        x += input_cloud->points[*pit].x;
        y += input_cloud->points[*pit].y;
        numPts++;
      }

      pcl::PointXYZ centroid;
      centroid.x = x / numPts;
      centroid.y = y / numPts;
      centroid.z = 0.0;

      if (centroid.x < 4 && centroid.x > -4 && centroid.y > -0.75 &&
          centroid.y < 0.75 && centroid.x != 0 && centroid.y != 0) {
        cluster_vec.push_back(cloud_cluster);

        // Get the centroid of the cluster
        clusterCentroids.push_back(centroid);
      }
    }

    // Ensure at least n clusters exist to publish (later clusters may be empty)
    while (cluster_vec.size() < n_clusters) {
      pcl::PointCloud<pcl::PointXYZ>::Ptr empty_cluster(
          new pcl::PointCloud<pcl::PointXYZ>);
      empty_cluster->points.push_back(pcl::PointXYZ(0, 0, 0));
      cluster_vec.push_back(empty_cluster);
    }

    while (clusterCentroids.size() < n_clusters) {
      pcl::PointXYZ centroid;
      centroid.x = 0.0;
      centroid.y = 0.0;
      centroid.z = 0.0;

      clusterCentroids.push_back(centroid);
    }

    std_msgs::Float32MultiArray cc;

    Eigen::Vector4f obstaclepoint;
    geometry_msgs::PointStamped bill;
    geometry_msgs::PointStamped m;
    bill.header.frame_id = "laser";
    for (int i = 0; i < n_clusters; i++) {
      {
        cc.data.push_back(clusterCentroids.at(i).x);
        cc.data.push_back(clusterCentroids.at(i).y);
        cc.data.push_back(clusterCentroids.at(i).z);

        if (clusterCentroids.at(i).x < 4 && clusterCentroids.at(i).y < 0.7 &&
            clusterCentroids.at(i).y > -0.7 &&
            clusterCentroids.at(i).x > -4 && clusterCentroids.at(i).x != 0 &&
            clusterCentroids.at(i).y != 0) {
          obstaclepoint[0] = clusterCentroids.at(i).x;
          obstaclepoint[1] = clusterCentroids.at(i).y;
          obstaclepoint[2] = clusterCentroids.at(i).z;
        }
      }
    }

    tf::StampedTransform transform;

    bill.point.x = obstaclepoint[0];
    bill.point.y = obstaclepoint[1];
    bill.point.z = obstaclepoint[2];

    try {
      tran->waitForTransform("/map", "/laser", ros::Time::now(),
                             ros::Duration(0.01));
      tran->lookupTransform("/map", "/laser", ros::Time(0), transform);
      tran->transformPoint("map", bill, m);
    } catch (tf::TransformException &ex) {
    }

    visualization_msgs::MarkerArray clusterMarkers1;
    visualization_msgs::Marker m1;

    m1.id = 0;
    m1.type = visualization_msgs::Marker::CUBE;
    m1.header.frame_id = "/map";
    m1.scale.x = 0.3;
    m1.scale.y = 0.3;
    m1.scale.z = 0.3;
    m1.action = visualization_msgs::Marker::ADD;
    m1.color.a = 1.0;
    m1.color.r = 1;
    m1.color.g = 0;
    m1.color.b = 0;

    std_msgs::Float32MultiArray cctemp;
    if (abs(obstaclepoint[0]) > 0.00001 && abs(obstaclepoint[1]) > 0.000001) {
      m1.pose.position.x = m.point.x;
      m1.pose.position.y = m.point.y;
      m1.pose.position.z = m.point.z;

      clusterMarkers1.markers.push_back(m1);

      cctemp.data.push_back(m.point.x);
      cctemp.data.push_back(m.point.y);
      cctemp.data.push_back(m.point.z);
    }

    // Publish cluster mid-points.
    cc_pos.publish(cctemp);
    markerPub1.publish(clusterMarkers1);

    KFT(cc);
  } // else
} // cloud_cb

int main(int argc, char **argv) {
  ros::init(argc, argv, "KFTracker");
  ros::NodeHandle nh;

  for(int i=0; i<n_clusters; i++)
  {
    KF_vec.push_back(KF_start);
  }

  ros::Subscriber sub = nh.subscribe("scan_cloud", 1, cloud_cb);

  tf::TransformListener lr(ros::Duration(10));
  tran = &lr;

  cc_pos = nh.advertise<std_msgs::Float32MultiArray>("cluster_center", 100); // clusterCenter1
  markerPub1 = nh.advertise<visualization_msgs::MarkerArray>("viz1", 1);

  ros::spin();
}
