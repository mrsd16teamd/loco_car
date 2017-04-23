cmake_minimum_required(VERSION 2.8.3)
project(publishpcl_nodelet)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(catkin REQUIRED COMPONENTS
  nodelet
  roscpp
  rospy
  std_msgs
  sensor_msgs
  tf
  laser_geometry
  laser_assembler
  tf2
)

## Setup include directories
include_directories(
  include
  include/Eigen
  ${catkin_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  ${Boost_LIBRARIES}
)

catkin_package(
  LIBRARIES nodelet_pcl
  CATKIN_DEPENDS
  nodelet
  roscpp
  rospy
  std_msgs
  sensor_msgs
  tf
  laser_geometry
  laser_assembler
  tf2
)

## Create the nodelet tutorial library
add_library(nodelet_pcl src/scan2cloud.cpp)
target_link_libraries(nodelet_pcl ${catkin_LIBRARIES} ${tf_LIBRARIES})
if(catkin_EXPORTED_LIBRARIES)
  add_dependencies(nodelet_pcl ${catkin_EXPORTED_LIBRARIES})
endif()

## Mark the nodelet library for installations
install(TARGETS nodelet_pcl
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})

## Mark other files for installation (e.g. launch and bag files, etc.)
install(FILES nodelet_pcl.xml scan2cloud.launch #plus_default.yaml
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})
