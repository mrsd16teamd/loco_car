#!/usr/bin/env python
import roslib; roslib.load_manifest('teleop_keyboard')
import rospy

from geometry_msgs.msg import Point
# from std_msgs import Char

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
		'k': 8,
		'r': 9
	     }



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
	print("Teleop keyboard running! Give me commands.\n" +
		"a: ramp\nb: iLQR static\nc: ramp and iLQR\nr: reset obs\nk: kill client")
	print("----------")

	try:
		while(1):
			key = getKey()

			if key in commandBindings.keys():
				command = commandBindings[key]
				print "Key: ", key, " Command: ", command
			else:
				command = 0
				print "Key: ", key, " Command: ", command

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
