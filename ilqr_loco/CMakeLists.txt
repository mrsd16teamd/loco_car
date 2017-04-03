cmake_minimum_required(VERSION 2.8.3)
project(ilqr_loco)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  genmsg
  actionlib_msgs
  actionlib
  nav_msgs
  geometry_msgs
  loco_msgs
)
find_package(Boost REQUIRED COMPONENTS system)

add_action_files(
   DIRECTORY action
   FILES TrajExec.action
)

generate_messages(
   DEPENDENCIES
   std_msgs
   nav_msgs
   geometry_msgs
   actionlib_msgs
   loco_msgs
)

catkin_package(
  CATKIN_DEPENDS roscpp std_msgs message_runtime actionlib_msgs loco_msgs
  INCLUDE_DIRS include
  LIBRARIES ilqr_loco
  # DEPENDS system_lib
)

include_directories(
  include ${catkin_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}
)

file(GLOB_RECURSE C_SOURCES "src/*.c")
add_definitions(-DPRNT=printf -DDEBUG_BACKPASS=1 -DDEBUG_FORWARDPASS=1 -DFULL_DDP=0 -O3)

add_executable(traj_client src/traj_client.cpp ${C_SOURCES})
target_link_libraries(traj_client ${catkin_LIBRARIES})
add_dependencies(traj_client ${ilqr_loco_EXPORTED_TARGETS})

add_executable(ramp_client src/ramp_client.cpp)
target_link_libraries(ramp_client ${catkin_LIBRARIES})
add_dependencies(ramp_client ${ilqr_loco_EXPORTED_TARGETS})

add_executable(traj_server src/traj_server.cpp)
target_link_libraries(traj_server ${catkin_LIBRARIES})
add_dependencies(traj_server ${ilqr_loco_EXPORTED_TARGETS})

add_executable(test_client src/test_client.cpp)
target_link_libraries(test_client ${catkin_LIBRARIES})
add_dependencies(test_client ${ilqr_loco_EXPORTED_TARGETS})
