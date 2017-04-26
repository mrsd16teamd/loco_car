#!/usr/bin/env python
import roslib; roslib.load_manifest('keyboard_command')
import rospy

from geometry_msgs.msg import Point
from std_msgs.msg import Char

import sys, select, termios, tty
from collections import OrderedDict

msg = """
Reading from the keyboard  and Publishing to Point!
---------------------------
CTRL-C to quit
"""

# Change the two dictionaries before to add/modify instructions
commandBindings = {
		'a': 1,
		'b': 2,
		'c': 3,
		'e': 5,
		'f': 6,
		'g': 7,
        'h': 13,
        'i': 14,
		'j': 15,
		'r': 8,
		'z': 10,
		'o': 12
	     }

commands = {
	'a': 'ramp (-> brake)',
	'b': 'iLQR open-loop from static',
	'c': 'iLQR mpc from static',
	'e': 'ramp -> iLQR open-loop',
	'f': 'ramp -> iLQR mpc',
	'g': 'ramp -> playback',
    'h': 'ramp -> playback -> iLQR open-loop',
    'i': 'ramp -> playback -> iLQR mpc',
	'j': 'ramp without brake',
	'z': 'execute initial control sequence',
	'r': 'reset obs',
	'o': 'insert fake obstacle 1m in front of robot'
}
instructions = OrderedDict(sorted(commands.items(), key=lambda t: t[0]))

def getKey():
	tty.setraw(sys.stdin.fileno())
	select.select([sys.stdin], [], [], 0)
	key = sys.stdin.read(1)
	termios.tcsetattr(sys.stdin, termios.TCSADRAIN, settings)
	return key

if __name__=="__main__":
	settings = termios.tcgetattr(sys.stdin)
	pub = rospy.Publisher('client_command', Point, queue_size = 1)
	rospy.init_node('keyboard_command')
	print("----------")
	print("Keyboard command running! Give me commands. \nMake sure these commands are synced with traj_client.")
	for key in instructions.keys():
		print(key + ': ' + instructions[key])
	print("----------")

	try:
		while(1):
			key = getKey()

			if key in commandBindings.keys():
				command = commandBindings[key]
				print "Key: ", key, " - ", instructions[key], " (", commandBindings[key] , ")"
			else:
				command = 0
				print "Key: ", key, " NO COMMAND"

				if (key == '\x03'):
					break

			cmd = Point()
			cmd.x = command;
			pub.publish(cmd)

	except:
		print 'Error.'

	finally:
		cmd = Point()
		cmd.x = 0
		pub.publish(cmd)

    	termios.tcsetattr(sys.stdin, termios.TCSADRAIN, settings)
