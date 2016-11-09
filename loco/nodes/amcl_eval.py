#!/usr/bin/env python
import rospy
import geometry_msgs.msg
import roslaunch
from math import pow
from numpy import mean

err_list = []

def launch_stuff():
    rospy.set_param('use_sim_time',True)
    package = 'rosbag'
    executable = 'play'
    bag_args = '/home/parallels/ros_bags/fvelab_amcl1.bag --topics /scan /map /map_metadata /scanmatch_odom /tf /clock -r 5'
    node = roslaunch.core.Node(package, executable)
    node.args=bag_args
    rosbag_launch = roslaunch.scriptapi.ROSLaunch()
    rosbag_launch.start()
    process = rosbag_launch.launch(node)

    package = 'loco'
    executable = 'amcl_score.py'
    node = roslaunch.core.Node(package, executable)
    node.args='output="screen"'
    eval_launch = roslaunch.scriptapi.ROSLaunch()
    eval_launch.start()
    process = eval_launch.launch(node)

def processPose(data):
    cov_xx = data.pose.covariance[0]
    cov_yy = data.pose.covariance[7]
    cov_tt = data.pose.covariance[35]
    cov_magn = (pow(cov_xx,2)+pow(cov_yy,2)+pow(cov_tt,2))**(1./3.)
    err_list.append(cov_tt)
    rospy.loginfo("avg: %s", mean(err_list))

def listener():
    # In ROS, nodes are uniquely named. If two nodes with the same
    # node are launched, the previous one is kicked off. The
    # anonymous=True flag means that rospy will choose a unique
    # name for our 'listener' node so that multiple listeners can
    # run simultaneously.
    rospy.init_node('launcher', anonymous=True)
    # rospy.Subscriber("amcl_pose", geometry_msgs.msg.PoseWithCovarianceStamped, processPose)
    # rospy.loginfo("launched stuff")
    rospy.spin()

if __name__ == '__main__':
    # listener()
    launch_stuff()
    listener()
    # rospy.spin()
