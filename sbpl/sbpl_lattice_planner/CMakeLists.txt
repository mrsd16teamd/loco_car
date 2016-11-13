cmake_minimum_required(VERSION 2.8.3)
project(sbpl_lattice_planner)
set(PROJECT_NAME sbpl_lattice_planner)

##############################################################################
# Find dependencies
##############################################################################

find_package(catkin REQUIRED COMPONENTS
  roscpp costmap_2d nav_core pluginlib geometry_msgs nav_msgs message_generation sbpl)

find_package(PkgConfig REQUIRED)

#find_package(sbpl REQUIRED)

include_directories(include ${catkin_INCLUDE_DIRS})

##############################################################################
# Generate messages
##############################################################################

add_message_files(FILES SBPLLatticePlannerStats.msg)
generate_messages(DEPENDENCIES geometry_msgs)

##############################################################################
# Define package
##############################################################################

#MESSAGE(SBPL_LIBRARIES ${SBPL_LIBRARIES})
# Since sbpl is found through pkg-config, we need to manually add its
# libraries and include directories into catkin_package
catkin_package(
#  INCLUDE_DIRS include ${sbpl_lattice_planner_INCLUDE_DIRS} ${sbpl_INCLUDE_DIRS}
  INCLUDE_DIRS include
  LIBRARIES ${sbpl_lattice_planner_LIBRARIES}
  CATKIN_DEPENDS ${THIS_PACKAGE_ROS_DEPS} message_runtime geometry_msgs
#  DEPENDS
)

##############################################################################
# Build
##############################################################################

add_library(${PROJECT_NAME} SHARED src/sbpl_lattice_planner.cpp)
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})
add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}_generated_messages_cpp ${PROJECT_NAME}_gencpp)

##############################################################################
# Install
##############################################################################

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
