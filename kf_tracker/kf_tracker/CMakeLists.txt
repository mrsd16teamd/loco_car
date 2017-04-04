cmake_minimum_required(VERSION 2.8.3)
project(kf_tracker)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  pcl_ros
  roscpp
  sensor_msgs
)
find_package( OpenCV REQUIRED )

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES kf_tracker
#  CATKIN_DEPENDS pcl_ros roscpp sensor_msgs
#  DEPENDS system_lib
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
# include_directories(include)
include_directories(
  ${catkin_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  include
)


## Declare a cpp library
# add_library(kf_tracker
#   src/${PROJECT_NAME}/kf_tracker.cpp
# )

## Declare a cpp executable
# add_executable(kf_tracker_node src/kf_tracker_node.cpp)
add_executable( tracker src/kf_tracker.cpp )
target_link_libraries ( tracker ${OpenCV_LIBRARIES} ${catkin_LIBRARIES})

## Add cmake target dependencies of the executable/library
## as an example, message headers may need to be generated before nodes
# add_dependencies(kf_tracker_node kf_tracker_generate_messages_cpp)

## Specify libraries to link a library or executable target against
# target_link_libraries(kf_tracker_node
#   ${catkin_LIBRARIES}
# )
