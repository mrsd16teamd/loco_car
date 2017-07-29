#
# MIT License
#
# Copyright (c) 2017 MRSD Team D - LoCo
# The Robotics Institute, Carnegie Mellon University
# http://mrsdprojects.ri.cmu.edu/2016teamd/
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

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
