#!/usr/bin/env python
import rospy
import geometry_msgs.msg
import roslaunch
from math import pow
from numpy import mean
import std_msgs.msg
import sys

err_list = []

def processPose(data):
    cov_xx = data.pose.covariance[0]
    cov_yy = data.pose.covariance[7]
    cov_tt = data.pose.covariance[35]
    cov_magn = (pow(cov_xx,2)+pow(cov_yy,2)+pow(cov_tt,2))**(1./3.)

    err_list.append(cov_magn)
    # rospy.loginfo("%s", mean(err_list))
    sys.stdout.write('\r')
    sys.stdout.write(str(mean(err_list)))
    sys.stdout.flush()

def listener():
    # In ROS, nodes are uniquely named. If two nodes with the same
    # node are launched, the previous one is kicked off. The
    # anonymous=True flag means that rospy will choose a unique
    # name for our 'listener' node so that multiple listeners can
    # run simultaneously.
    rospy.init_node('listener', anonymous=True)
    rospy.Subscriber("amcl_pose", geometry_msgs.msg.PoseWithCovarianceStamped, processPose)
    # rospy.loginfo("initialized listener")

    rospy.spin()


if __name__ == '__main__':
    listener()
