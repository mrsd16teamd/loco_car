cmake_minimum_required(VERSION 2.8.3)
project(ilqr_loco)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

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
  INCLUDE_DIRS include include/iLQG
  LIBRARIES ilqr_loco
)

include_directories(
  include include/iLQG ${catkin_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}
)

file(GLOB_RECURSE C_SOURCES "src/iLQG/*.c")
add_definitions(-DPRNT=printf -DDEBUG_ILQG=1 -DDEBUG_BACKPASS=0 -DDEBUG_FORWARDPASS=0 -DFULL_DDP=0 -DMULTI_THREADED=0 -O3)

add_executable(traj_client src/traj_client.cpp
													 src/ilqr_planner.cpp
													 src/traj_client_ramp.cpp
													 src/traj_client_ilqr.cpp
													 src/traj_client_params.cpp
                           src/msg_utils.cpp
												   ${C_SOURCES})
target_link_libraries(traj_client ${catkin_LIBRARIES})
add_dependencies(traj_client ${ilqr_loco_EXPORTED_TARGETS})

add_executable(traj_server src/traj_server.cpp)
target_link_libraries(traj_server ${catkin_LIBRARIES})
add_dependencies(traj_server ${ilqr_loco_EXPORTED_TARGETS})
