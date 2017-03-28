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

#add_executable(ilqr_planner src/ilqr_planner.cpp)
#target_link_libraries(ilqr_planner ${catkin_LIBRARIES})
#add_dependencies(ilqr_planner ${ilqr_loco_EXPORTED_TARGETS} )

add_executable(ilqr_executer src/ilqr_executer.cpp)
target_link_libraries(ilqr_executer ${catkin_LIBRARIES})
add_dependencies(ilqr_executer ${ilqr_loco_EXPORTED_TARGETS} )
