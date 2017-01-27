# generated from genmsg/cmake/pkg-genmsg.cmake.em

message(STATUS "ackermann_msgs: 2 messages, 0 services")

set(MSG_I_FLAGS "-Iackermann_msgs:/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg;-Istd_msgs:/opt/ros/indigo/share/std_msgs/cmake/../msg")

# Find all generators
find_package(gencpp REQUIRED)
find_package(genlisp REQUIRED)
find_package(genpy REQUIRED)

add_custom_target(ackermann_msgs_generate_messages ALL)

# verify that message/service dependencies have not changed since configure



get_filename_component(_filename "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDrive.msg" NAME_WE)
add_custom_target(_ackermann_msgs_generate_messages_check_deps_${_filename}
  COMMAND ${CATKIN_ENV} ${PYTHON_EXECUTABLE} ${GENMSG_CHECK_DEPS_SCRIPT} "ackermann_msgs" "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDrive.msg" ""
)

get_filename_component(_filename "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDriveStamped.msg" NAME_WE)
add_custom_target(_ackermann_msgs_generate_messages_check_deps_${_filename}
  COMMAND ${CATKIN_ENV} ${PYTHON_EXECUTABLE} ${GENMSG_CHECK_DEPS_SCRIPT} "ackermann_msgs" "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDriveStamped.msg" "ackermann_msgs/AckermannDrive:std_msgs/Header"
)

#
#  langs = gencpp;genlisp;genpy
#

### Section generating for lang: gencpp
### Generating Messages
_generate_msg_cpp(ackermann_msgs
  "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDrive.msg"
  "${MSG_I_FLAGS}"
  ""
  ${CATKIN_DEVEL_PREFIX}/${gencpp_INSTALL_DIR}/ackermann_msgs
)
_generate_msg_cpp(ackermann_msgs
  "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDriveStamped.msg"
  "${MSG_I_FLAGS}"
  "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDrive.msg;/opt/ros/indigo/share/std_msgs/cmake/../msg/Header.msg"
  ${CATKIN_DEVEL_PREFIX}/${gencpp_INSTALL_DIR}/ackermann_msgs
)

### Generating Services

### Generating Module File
_generate_module_cpp(ackermann_msgs
  ${CATKIN_DEVEL_PREFIX}/${gencpp_INSTALL_DIR}/ackermann_msgs
  "${ALL_GEN_OUTPUT_FILES_cpp}"
)

add_custom_target(ackermann_msgs_generate_messages_cpp
  DEPENDS ${ALL_GEN_OUTPUT_FILES_cpp}
)
add_dependencies(ackermann_msgs_generate_messages ackermann_msgs_generate_messages_cpp)

# add dependencies to all check dependencies targets
get_filename_component(_filename "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDrive.msg" NAME_WE)
add_dependencies(ackermann_msgs_generate_messages_cpp _ackermann_msgs_generate_messages_check_deps_${_filename})
get_filename_component(_filename "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDriveStamped.msg" NAME_WE)
add_dependencies(ackermann_msgs_generate_messages_cpp _ackermann_msgs_generate_messages_check_deps_${_filename})

# target for backward compatibility
add_custom_target(ackermann_msgs_gencpp)
add_dependencies(ackermann_msgs_gencpp ackermann_msgs_generate_messages_cpp)

# register target for catkin_package(EXPORTED_TARGETS)
list(APPEND ${PROJECT_NAME}_EXPORTED_TARGETS ackermann_msgs_generate_messages_cpp)

### Section generating for lang: genlisp
### Generating Messages
_generate_msg_lisp(ackermann_msgs
  "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDrive.msg"
  "${MSG_I_FLAGS}"
  ""
  ${CATKIN_DEVEL_PREFIX}/${genlisp_INSTALL_DIR}/ackermann_msgs
)
_generate_msg_lisp(ackermann_msgs
  "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDriveStamped.msg"
  "${MSG_I_FLAGS}"
  "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDrive.msg;/opt/ros/indigo/share/std_msgs/cmake/../msg/Header.msg"
  ${CATKIN_DEVEL_PREFIX}/${genlisp_INSTALL_DIR}/ackermann_msgs
)

### Generating Services

### Generating Module File
_generate_module_lisp(ackermann_msgs
  ${CATKIN_DEVEL_PREFIX}/${genlisp_INSTALL_DIR}/ackermann_msgs
  "${ALL_GEN_OUTPUT_FILES_lisp}"
)

add_custom_target(ackermann_msgs_generate_messages_lisp
  DEPENDS ${ALL_GEN_OUTPUT_FILES_lisp}
)
add_dependencies(ackermann_msgs_generate_messages ackermann_msgs_generate_messages_lisp)

# add dependencies to all check dependencies targets
get_filename_component(_filename "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDrive.msg" NAME_WE)
add_dependencies(ackermann_msgs_generate_messages_lisp _ackermann_msgs_generate_messages_check_deps_${_filename})
get_filename_component(_filename "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDriveStamped.msg" NAME_WE)
add_dependencies(ackermann_msgs_generate_messages_lisp _ackermann_msgs_generate_messages_check_deps_${_filename})

# target for backward compatibility
add_custom_target(ackermann_msgs_genlisp)
add_dependencies(ackermann_msgs_genlisp ackermann_msgs_generate_messages_lisp)

# register target for catkin_package(EXPORTED_TARGETS)
list(APPEND ${PROJECT_NAME}_EXPORTED_TARGETS ackermann_msgs_generate_messages_lisp)

### Section generating for lang: genpy
### Generating Messages
_generate_msg_py(ackermann_msgs
  "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDrive.msg"
  "${MSG_I_FLAGS}"
  ""
  ${CATKIN_DEVEL_PREFIX}/${genpy_INSTALL_DIR}/ackermann_msgs
)
_generate_msg_py(ackermann_msgs
  "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDriveStamped.msg"
  "${MSG_I_FLAGS}"
  "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDrive.msg;/opt/ros/indigo/share/std_msgs/cmake/../msg/Header.msg"
  ${CATKIN_DEVEL_PREFIX}/${genpy_INSTALL_DIR}/ackermann_msgs
)

### Generating Services

### Generating Module File
_generate_module_py(ackermann_msgs
  ${CATKIN_DEVEL_PREFIX}/${genpy_INSTALL_DIR}/ackermann_msgs
  "${ALL_GEN_OUTPUT_FILES_py}"
)

add_custom_target(ackermann_msgs_generate_messages_py
  DEPENDS ${ALL_GEN_OUTPUT_FILES_py}
)
add_dependencies(ackermann_msgs_generate_messages ackermann_msgs_generate_messages_py)

# add dependencies to all check dependencies targets
get_filename_component(_filename "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDrive.msg" NAME_WE)
add_dependencies(ackermann_msgs_generate_messages_py _ackermann_msgs_generate_messages_check_deps_${_filename})
get_filename_component(_filename "/home/ubuntu/catkin_ws/src/utilities/ackermann_msgs/msg/AckermannDriveStamped.msg" NAME_WE)
add_dependencies(ackermann_msgs_generate_messages_py _ackermann_msgs_generate_messages_check_deps_${_filename})

# target for backward compatibility
add_custom_target(ackermann_msgs_genpy)
add_dependencies(ackermann_msgs_genpy ackermann_msgs_generate_messages_py)

# register target for catkin_package(EXPORTED_TARGETS)
list(APPEND ${PROJECT_NAME}_EXPORTED_TARGETS ackermann_msgs_generate_messages_py)



if(gencpp_INSTALL_DIR AND EXISTS ${CATKIN_DEVEL_PREFIX}/${gencpp_INSTALL_DIR}/ackermann_msgs)
  # install generated code
  install(
    DIRECTORY ${CATKIN_DEVEL_PREFIX}/${gencpp_INSTALL_DIR}/ackermann_msgs
    DESTINATION ${gencpp_INSTALL_DIR}
  )
endif()
add_dependencies(ackermann_msgs_generate_messages_cpp std_msgs_generate_messages_cpp)

if(genlisp_INSTALL_DIR AND EXISTS ${CATKIN_DEVEL_PREFIX}/${genlisp_INSTALL_DIR}/ackermann_msgs)
  # install generated code
  install(
    DIRECTORY ${CATKIN_DEVEL_PREFIX}/${genlisp_INSTALL_DIR}/ackermann_msgs
    DESTINATION ${genlisp_INSTALL_DIR}
  )
endif()
add_dependencies(ackermann_msgs_generate_messages_lisp std_msgs_generate_messages_lisp)

if(genpy_INSTALL_DIR AND EXISTS ${CATKIN_DEVEL_PREFIX}/${genpy_INSTALL_DIR}/ackermann_msgs)
  install(CODE "execute_process(COMMAND \"/usr/bin/python\" -m compileall \"${CATKIN_DEVEL_PREFIX}/${genpy_INSTALL_DIR}/ackermann_msgs\")")
  # install generated code
  install(
    DIRECTORY ${CATKIN_DEVEL_PREFIX}/${genpy_INSTALL_DIR}/ackermann_msgs
    DESTINATION ${genpy_INSTALL_DIR}
  )
endif()
add_dependencies(ackermann_msgs_generate_messages_py std_msgs_generate_messages_py)
