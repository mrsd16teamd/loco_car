cmake_minimum_required(VERSION 2.8.3)
project(actionlib_tutorials)

find_package(catkin REQUIRED COMPONENTS roscpp actionlib actionlib_msgs)
find_package(Boost REQUIRED COMPONENTS system)

add_action_files(
  DIRECTORY action
  FILES Fibonacci.action
)

generate_messages(
  DEPENDENCIES actionlib_msgs std_msgs
)

catkin_package(
  CATKIN_DEPENDS actionlib_msgs
)

include_directories(include ${catkin_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

add_executable(fibonacci_server src/fibonacci_server.cpp)
add_executable(fibonacci_client src/fibonacci_client.cpp)
add_executable(fibonacci_client2 src/fibonacci_client2.cpp)

target_link_libraries(
  fibonacci_server
  ${catkin_LIBRARIES}
)

target_link_libraries(
  fibonacci_client
  ${catkin_LIBRARIES}
)

target_link_libraries(
  fibonacci_client2
  ${catkin_LIBRARIES}
)


add_dependencies(
  fibonacci_server
  ${actionlib_tutorials_EXPORTED_TARGETS}
)

add_dependencies(
  fibonacci_client
  ${actionlib_tutorials_EXPORTED_TARGETS}
)

add_dependencies(
  fibonacci_client2
  ${actionlib_tutorials_EXPORTED_TARGETS}
)
