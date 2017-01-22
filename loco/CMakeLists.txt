cmake_minimum_required(VERSION 2.8.3)
project(loco)
find_package(catkin REQUIRED roscpp std_msgs tf move_base_msgs geometry_msgs sensor_msgs nav_msgs)


include_directories(
  ${catkin_INCLUDE_DIRS}
)
catkin_package()




add_executable(amcl_globalinit src/amcl_globalinit.cpp)
target_link_libraries(amcl_globalinit ${catkin_LIBRARIES})

add_executable(dist2goal src/dist2goal.cpp)
target_link_libraries(dist2goal ${catkin_LIBRARIES})

add_executable(cmd_vel_odom src/cmd_vel_odom.cpp)
target_link_libraries(cmd_vel_odom ${catkin_LIBRARIES})
