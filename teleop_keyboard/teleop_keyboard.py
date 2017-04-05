#!/usr/bin/env python
import roslib; roslib.load_manifest('teleop_keyboard')
import rospy

from geometry_msgs.msg import Point
# from std_msgs import Char

import sys, select, termios, tty

msg = """
Reading from the keyboard  and Publishing to Point!
---------------------------
1 - mode 1
2 - mode 2
3 - mode 3

anything else : stop

CTRL-C to quit
"""

commandBindings = {
		'a': 1,
		'b': 2,
		'c': 3,
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

	try:
		while(1):
			key = getKey()
			print "------"
			print "Key pressed:", key

			if key in commandBindings.keys():
				command = commandBindings[key]
				print "Command: ", command
			else:
				command = 0
				print "Command: ", command

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
