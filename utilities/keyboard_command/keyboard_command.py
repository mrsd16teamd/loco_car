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
		'd': 4,
		'e': 5,
		'f': 6,
		'g': 7,
		'r': 8,
		'k': 9
	     }

commands = {
	'a': 'ramp',
	'b': 'iLQR open-loop from static',
	'c': 'ramp and iLQR open-loop',
	'd': 'ramp and iLQR rh',
	'e': 'iLQR rh from static',
	'f': 'iLQR sparse replanning from static',
	'g': 'iLQR open-loop with PID corrections from static'
	'r': 'reset obs',
	'k': 'kill client'
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
				print "Key: ", key, " - ", instructions[key]
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
