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
