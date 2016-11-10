cmake_minimum_required(VERSION 2.8.3)
project(base_controller)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_BUILD_TYPE Release) 

find_package(catkin REQUIRED COMPONENTS
  ackermann_msgs
  roscpp
  geometry_msgs
)

# Node: catkin_package macro is required
catkin_package(
    CATKIN_DEPENDS
    )


###########
## Build ##
###########

include_directories(include)
include_directories(
  ${catkin_INCLUDE_DIRS}
)

add_executable(cmd_vel_to_ackermann_drive src/cmd_vel_to_ackermann_drive.cpp)
target_link_libraries(cmd_vel_to_ackermann_drive ${catkin_LIBRARIES} )


