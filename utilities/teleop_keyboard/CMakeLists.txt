cmake_minimum_required(VERSION 2.8.3)
project(teleop_keyboard)

find_package(catkin REQUIRED std_msgs)

catkin_package()

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
catkin_install_python(PROGRAMS
   teleop_keyboard.py
   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
