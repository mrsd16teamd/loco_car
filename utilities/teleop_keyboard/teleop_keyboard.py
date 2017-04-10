#!/usr/bin/env python
import roslib; roslib.load_manifest('teleop_keyboard')
import rospy

from geometry_msgs.msg import Point
from std_msgs.msg import Char

import sys, select, termios, tty

msg = """
Reading from the keyboard  and Publishing to Point!
---------------------------
CTRL-C to quit
"""

commandBindings = {
		'a': 1,
		'b': 2,
		'c': 3,
		'd': 4,
		'r': 8,
		'k': 9
	     }

instructions = {
	'a': 'ramp',
	'b': 'iLQR static',
	'c': 'ramp and iLQR open-loop',
	'd': 'ramp and iLQR mpc',
	'r': 'reset obs',
	'k': 'kill client'
}
#DONT CHANGE THESE! TOO MUCH WORK


def getKey():
	tty.setraw(sys.stdin.fileno())
	select.select([sys.stdin], [], [], 0)
	key = sys.stdin.read(1)
	termios.tcsetattr(sys.stdin, termios.TCSADRAIN, settings)
	return key

if __name__=="__main__":
	settings = termios.tcgetattr(sys.stdin)
	pub = rospy.Publisher('client_command', Point, queue_size = 1)
	rospy.init_node('teleop_keyboard')
	print("----------")
	print("Teleop keyboard running! Give me commands. \nMake sure these commands are synced with traj_client.\n" +
		"a: ramp\nb: iLQR static\nc: ramp and iLQR open-loop\nd: ramp and iLQR mpc\nr: reset obs\nk: kill client")
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
