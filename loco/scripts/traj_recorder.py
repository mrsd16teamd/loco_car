#!/usr/bin/env python

import rospy
from geometry_msgs.msg import Twist
import os

record_file = "record.txt"    # This should be a txt file
filepath = os.path.abspath(record_file)
try:
    file = open(filepath, "w")
except IOError:
    rospy.logerr("Could not open '%s' to write.", command_file)

def callback(data):
    vx = data.linear.x
    steer = data.angular.z

    rospy.loginfo("I heard vx: %f, steer: %f",vx,steer)
    command = "u " + str(vx) + " " + str(steer)
    file.write(command + '\n')

def recorder():
    rospy.init_node('traj_recorder')
    rospy.Subscriber("cmd_vel", Twist, callback)

    rospy.loginfo("Started control input recorder.")

    rospy.spin()

if __name__ == '__main__':
    try:
        recorder()
    except rospy.ROSInterruptException:
        pass
