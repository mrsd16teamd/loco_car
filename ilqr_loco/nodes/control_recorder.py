#!/usr/bin/env python
import rospy
from geometry_msgs.msg import Twist



def callback(data):
    rospy.loginfo(rospy.get_caller_id() + "I heard %s", data.data)

def listener():
    rospy.init_node('control_recorder', anonymous=True)
    rospy.Subscriber("cmd_vel_out", Twist, callback)
    rospy.spin()

if __name__ == '__main__':
    listener()
