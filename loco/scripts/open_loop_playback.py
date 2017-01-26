#!/usr/bin/env python

import rospy
from geometry_msgs.msg import Twist
import os

command_file = "controls.txt"
filepath = os.path.abspath(command_file)

def talker():
    pub = rospy.Publisher('cmd_vel', Twist, queue_size=1)
    rospy.init_node('open_loop_playback', anonymous=True)
    try:
        file = open(filepath, "r")
    except IOError:
        rospy.logerr("'%s' not found.", command_file)
        return 0
    line1 = file.readline().rstrip().split()
    assert(line1[0]=="dt")
    dt = float(line1[1])/1000 # in [ms], change to [s]
    rospy.loginfo("Starting open loop command play back at " + str(1/dt) + "Hz.")
    rospy.sleep(2.0)

    rate = rospy.Rate(1/dt) # [Hz]. Change this according to dt

    while not rospy.is_shutdown():
        # Read file line by line
        line = file.readline().rstrip().split()
        if (line[0]=="end"):
            rospy.loginfo("Finished executing control inputs.")
            break
        assert(line[0]=="u")
        vx = float(line[1])
        steer = float(line[2])
        rospy.loginfo("vx: " + str(vx) + ", steer: " + str(steer))

        # Fill twist message
        msg = Twist()
        msg.linear.x = vx
        msg.angular.z = steer

        # rospy.loginfo(msg)
        pub.publish(msg)
        rate.sleep()

if __name__ == '__main__':
    try:
        talker()
    except rospy.ROSInterruptException:
        pass
