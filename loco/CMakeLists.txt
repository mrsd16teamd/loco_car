cmake_minimum_required(VERSION 2.8.3)
project(loco)
find_package(catkin REQUIRED
  roscpp
  std_msgs
  tf
  move_base_msgs
  geometry_msgs
  sensor_msgs
  nav_msgs
  rospy
  message_generation
  )

# For actions
#find_package(catkin REQUIRED genmsg actionlib_msgs actionlib)
#add_action_files(DIRECTORY action FILES
  # DoDishes.action TODO add our own here
#)

include_directories(
  ${catkin_INCLUDE_DIRS}
)

#add_message_files(
#   FILES
#   Trajectory.msg
#)

# generate_messages(
#    DEPENDENCIES
#   std_msgs
#   nav_msgs
#   geometry_msgs
#   actionlib_msgs
#)

catkin_package(
  CATKIN_DEPENDS message_runtime
)

# Manually add executables
add_executable(amcl_globalinit src/amcl_globalinit.cpp)
target_link_libraries(amcl_globalinit ${catkin_LIBRARIES})

add_executable(dist2goal src/dist2goal.cpp)
target_link_libraries(dist2goal ${catkin_LIBRARIES})

add_executable(cmd_vel_odom src/cmd_vel_odom.cpp)
target_link_libraries(cmd_vel_odom ${catkin_LIBRARIES})

add_executable(playback_traj src/playback_traj)
target_link_libraries(playback_traj ${catkin_LIBRARIES})

add_executable(record_traj src/record_traj)
target_link_libraries(record_traj ${catkin_LIBRARIES})

add_executable(amcl_echo src/amcl_echo)
target_link_libraries(amcl_echo ${catkin_LIBRARIES})
