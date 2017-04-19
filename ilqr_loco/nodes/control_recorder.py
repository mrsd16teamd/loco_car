#!/usr/bin/env python
# -*- coding: utf-8 -*-

import rospy
from geometry_msgs.msg import Twist


u_seq = []
filename = '0mps_T20_20Hz.yaml'

def callback(msg):
	print("I heard: " + str(msg.linear.x) + ", " + str(msg.angular.z))
	u_seq.append(msg.linear.x)
	u_seq.append(msg.angular.z)

def listener():
    rospy.init_node('control_recorder', anonymous=True)
    rospy.Subscriber("cmd_vel_out", Twist, callback)
    rospy.spin()

def write_to_file():
	n_commands = len(u_seq)/2
	with open(filename, "w") as f:
		target.write("T_horizon: " + str(n_commands) + "\n")
		target.write("init_control_seq: [")
		for i in range(n_commands):
			target.write(str(u_seq[2*i] ) + ", " + str(u_seq[2*i+1]) )
			if i != (n_commands-1):
				target.write(", ")
		target.write("]")
	print("Done writing to file.")

if __name__ == '__main__':
	print("Writing to: " + filename)
	listener()
	raw_input("Hit enter when you're done.")
	write_to_file()
