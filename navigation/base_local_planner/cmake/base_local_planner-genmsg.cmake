# generated from genmsg/cmake/pkg-genmsg.cmake.em

message(STATUS "base_local_planner: 1 messages, 0 services")

set(MSG_I_FLAGS "-Ibase_local_planner:/home/ubuntu/catkin_ws/src/navigation/base_local_planner/msg;-Istd_msgs:/opt/ros/indigo/share/std_msgs/cmake/../msg")

# Find all generators
find_package(gencpp REQUIRED)
find_package(genlisp REQUIRED)
find_package(genpy REQUIRED)

add_custom_target(base_local_planner_generate_messages ALL)

# verify that message/service dependencies have not changed since configure



get_filename_component(_filename "/home/ubuntu/catkin_ws/src/navigation/base_local_planner/msg/Position2DInt.msg" NAME_WE)
add_custom_target(_base_local_planner_generate_messages_check_deps_${_filename}
  COMMAND ${CATKIN_ENV} ${PYTHON_EXECUTABLE} ${GENMSG_CHECK_DEPS_SCRIPT} "base_local_planner" "/home/ubuntu/catkin_ws/src/navigation/base_local_planner/msg/Position2DInt.msg" ""
)

#
#  langs = gencpp;genlisp;genpy
#

### Section generating for lang: gencpp
### Generating Messages
_generate_msg_cpp(base_local_planner
  "/home/ubuntu/catkin_ws/src/navigation/base_local_planner/msg/Position2DInt.msg"
  "${MSG_I_FLAGS}"
  ""
  ${CATKIN_DEVEL_PREFIX}/${gencpp_INSTALL_DIR}/base_local_planner
)

### Generating Services

### Generating Module File
_generate_module_cpp(base_local_planner
  ${CATKIN_DEVEL_PREFIX}/${gencpp_INSTALL_DIR}/base_local_planner
  "${ALL_GEN_OUTPUT_FILES_cpp}"
)

add_custom_target(base_local_planner_generate_messages_cpp
  DEPENDS ${ALL_GEN_OUTPUT_FILES_cpp}
)
add_dependencies(base_local_planner_generate_messages base_local_planner_generate_messages_cpp)

# add dependencies to all check dependencies targets
get_filename_component(_filename "/home/ubuntu/catkin_ws/src/navigation/base_local_planner/msg/Position2DInt.msg" NAME_WE)
add_dependencies(base_local_planner_generate_messages_cpp _base_local_planner_generate_messages_check_deps_${_filename})

# target for backward compatibility
add_custom_target(base_local_planner_gencpp)
add_dependencies(base_local_planner_gencpp base_local_planner_generate_messages_cpp)

# register target for catkin_package(EXPORTED_TARGETS)
list(APPEND ${PROJECT_NAME}_EXPORTED_TARGETS base_local_planner_generate_messages_cpp)

### Section generating for lang: genlisp
### Generating Messages
_generate_msg_lisp(base_local_planner
  "/home/ubuntu/catkin_ws/src/navigation/base_local_planner/msg/Position2DInt.msg"
  "${MSG_I_FLAGS}"
  ""
  ${CATKIN_DEVEL_PREFIX}/${genlisp_INSTALL_DIR}/base_local_planner
)

### Generating Services

### Generating Module File
_generate_module_lisp(base_local_planner
  ${CATKIN_DEVEL_PREFIX}/${genlisp_INSTALL_DIR}/base_local_planner
  "${ALL_GEN_OUTPUT_FILES_lisp}"
)

add_custom_target(base_local_planner_generate_messages_lisp
  DEPENDS ${ALL_GEN_OUTPUT_FILES_lisp}
)
add_dependencies(base_local_planner_generate_messages base_local_planner_generate_messages_lisp)

# add dependencies to all check dependencies targets
get_filename_component(_filename "/home/ubuntu/catkin_ws/src/navigation/base_local_planner/msg/Position2DInt.msg" NAME_WE)
add_dependencies(base_local_planner_generate_messages_lisp _base_local_planner_generate_messages_check_deps_${_filename})

# target for backward compatibility
add_custom_target(base_local_planner_genlisp)
add_dependencies(base_local_planner_genlisp base_local_planner_generate_messages_lisp)

# register target for catkin_package(EXPORTED_TARGETS)
list(APPEND ${PROJECT_NAME}_EXPORTED_TARGETS base_local_planner_generate_messages_lisp)

### Section generating for lang: genpy
### Generating Messages
_generate_msg_py(base_local_planner
  "/home/ubuntu/catkin_ws/src/navigation/base_local_planner/msg/Position2DInt.msg"
  "${MSG_I_FLAGS}"
  ""
  ${CATKIN_DEVEL_PREFIX}/${genpy_INSTALL_DIR}/base_local_planner
)

### Generating Services

### Generating Module File
_generate_module_py(base_local_planner
  ${CATKIN_DEVEL_PREFIX}/${genpy_INSTALL_DIR}/base_local_planner
  "${ALL_GEN_OUTPUT_FILES_py}"
)

add_custom_target(base_local_planner_generate_messages_py
  DEPENDS ${ALL_GEN_OUTPUT_FILES_py}
)
add_dependencies(base_local_planner_generate_messages base_local_planner_generate_messages_py)

# add dependencies to all check dependencies targets
get_filename_component(_filename "/home/ubuntu/catkin_ws/src/navigation/base_local_planner/msg/Position2DInt.msg" NAME_WE)
add_dependencies(base_local_planner_generate_messages_py _base_local_planner_generate_messages_check_deps_${_filename})

# target for backward compatibility
add_custom_target(base_local_planner_genpy)
add_dependencies(base_local_planner_genpy base_local_planner_generate_messages_py)

# register target for catkin_package(EXPORTED_TARGETS)
list(APPEND ${PROJECT_NAME}_EXPORTED_TARGETS base_local_planner_generate_messages_py)



if(gencpp_INSTALL_DIR AND EXISTS ${CATKIN_DEVEL_PREFIX}/${gencpp_INSTALL_DIR}/base_local_planner)
  # install generated code
  install(
    DIRECTORY ${CATKIN_DEVEL_PREFIX}/${gencpp_INSTALL_DIR}/base_local_planner
    DESTINATION ${gencpp_INSTALL_DIR}
  )
endif()
add_dependencies(base_local_planner_generate_messages_cpp std_msgs_generate_messages_cpp)

if(genlisp_INSTALL_DIR AND EXISTS ${CATKIN_DEVEL_PREFIX}/${genlisp_INSTALL_DIR}/base_local_planner)
  # install generated code
  install(
    DIRECTORY ${CATKIN_DEVEL_PREFIX}/${genlisp_INSTALL_DIR}/base_local_planner
    DESTINATION ${genlisp_INSTALL_DIR}
  )
endif()
add_dependencies(base_local_planner_generate_messages_lisp std_msgs_generate_messages_lisp)

if(genpy_INSTALL_DIR AND EXISTS ${CATKIN_DEVEL_PREFIX}/${genpy_INSTALL_DIR}/base_local_planner)
  install(CODE "execute_process(COMMAND \"/usr/bin/python\" -m compileall \"${CATKIN_DEVEL_PREFIX}/${genpy_INSTALL_DIR}/base_local_planner\")")
  # install generated code
  install(
    DIRECTORY ${CATKIN_DEVEL_PREFIX}/${genpy_INSTALL_DIR}/base_local_planner
    DESTINATION ${genpy_INSTALL_DIR}
  )
endif()
add_dependencies(base_local_planner_generate_messages_py std_msgs_generate_messages_py)
