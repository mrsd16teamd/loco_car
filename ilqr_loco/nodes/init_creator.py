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

import sys, select, termios, tty
import numpy as np

u_seq = []
filename = 'generated.yaml'
filepath = '/home/ubuntu/catkin_ws/src/ilqr_loco/control_seq/' + filename

T = 50
vel = 1.0

def fill_sequence():
	for i in range(T):
		throttle = vel + 0.1*np.random.normal()
		steer = 0.2*np.random.normal()
		u_seq.append(throttle)
		u_seq.append(steer)

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
	fill_sequence()
	write_to_file()
