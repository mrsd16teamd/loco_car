cmake_minimum_required(VERSION 2.8.3)
project(loco)
find_package(catkin REQUIRED roscpp std_msgs)
catkin_package()

add_executable(amcl_globalinit src/amcl_globalinit.cpp)
target_link_libraries(amcl_globalinit ${catkin_LIBRARIES})
