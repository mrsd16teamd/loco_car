#!/usr/bin/env python

import rospy
from geometry_msgs.msg import Twist
import sys, select, termios, tty


u_seq = []
filename = 'most_recent_recording.yaml'
filepath = '/home/ubuntu/catkin_ws/src/ilqr_loco/control_seq/' + filename


def callback(msg):
	if (msg.linear.x < 0.05 and msg.angular.z < 0.05):
		return
	else:
		print("I heard: " + str(msg.linear.x) + ", " + str(msg.angular.z))
	u_seq.append(msg.linear.x)
	u_seq.append(msg.angular.z)

def listener():
	print("Writing to: " + filepath)
	print("Go!")
	print("Hit ctrl+c when you're done.")

	rospy.init_node('control_recorder', anonymous=True)
	rospy.Subscriber("cmd_vel_out", Twist, callback)

	rospy.spin()

def write_to_file():
	n_commands = len(u_seq)/2
	with open(filepath, "w+") as f:
		f.write("T_horizon: " + str(n_commands) + "\n")
		f.write("init_control_seq: [\n")
		for i in range(n_commands):
			f.write(str(u_seq[2*i] ) + ", " + str(u_seq[2*i+1]))
			if i != (n_commands-1):
				f.write(", ")
			f.write("\n")
		f.write("]")
	print("Done writing to file")

if __name__ == '__main__':
	rospy.on_shutdown(write_to_file)
	listener()
