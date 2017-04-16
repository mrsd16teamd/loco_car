#ifndef _TRYGETPARAM_H_
#define _TRYGETPARAM_H_

// Assumptions
// * all nodes that use this have node handle named nh
// * name is a string
// * variable is a member variable name
#define TRYGETPARAM(name, variable) if(!nh.getParam(name, variable)){ROS_ERROR("MISSING PARAMETER: %s.", name); ros::shutdown();}

#endif
