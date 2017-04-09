

void publish_cloud(ros::Publisher& pub, pcl::PointCloud<pcl::PointXYZ>::Ptr cluster, const sensor_msgs::PointCloud2ConstPtr& input){
    sensor_msgs::PointCloud2::Ptr clustermsg (new sensor_msgs::PointCloud2);
    pcl::toROSMsg (*cluster , *clustermsg);
    clustermsg->header.frame_id = input->header.frame_id;
    clustermsg->header.stamp = ros::Time::now();
    pub.publish (*clustermsg);

}


double euclidean_distance(geometry_msgs::Point& p1, geometry_msgs::Point& p2)
{
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}


std::pair<int,int> findIndexOfMin(std::vector<std::vector<float> > distMat)
{
   // cout<<"findIndexOfMin cALLED\n";
    std::pair<int,int>minIndex;
    float minEl=std::numeric_limits<float>::max();
  //  cout<<"minEl="<<minEl<<"\n";
    for (int i=0; i<distMat.size();i++)
        for(int j=0;j<distMat.at(0).size();j++)
        {
            if( distMat[i][j]<minEl)
            {
                minEl=distMat[i][j];
                minIndex=std::make_pair(i,j);

            }

        }
  //  cout<<"minIndex="<<minIndex.first<<","<<minIndex.second<<"\n";
    return minIndex;
}
