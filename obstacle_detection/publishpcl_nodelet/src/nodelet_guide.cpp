// NOdelet everything


#include "ros/ros.h"
#include "nodelet/nodelet.h"

namespace @(namespace)
{

  class @(NodeletClass) : public nodelet::Nodelet
  {
  public:
  @(NodeletClass)();

  private:
  virtual void onInit(){
  nh = getNodeHandle();
  private_nh = getPrivateNodeHandle();
  timer_ = nh.createTimer(ros::Duration(1.0), boost::bind(& @(NodeletClass)::timerCb, this, _1));
  sub_ = nh.subscribe("incoming_chatter", 10, boost::bind(& @(NodeletClass)::messageCb, this, _1));
  pub_ = private_nh.advertise<std_msgs::String>("outgoing_chatter", 10);
  };

  void timerCb(const ros::TimerEvent& event){
  // Using timers is the preferred 'ROS way' to manual threading
  NODELET_INFO_STREAM("The time is now " << event.current_real);
  }

  // must use a ConstPtr callback to use zero-copy transport
  void messageCb(const std_msgs::StringConstPtr message){

  // can republish the old message no problem, since we're not modifying it
  pub_.publish(message);

  std_msgs::String new_message;
  new_message.data = message.data + " fizz buzz";
  pub_.publish(new_message);

  // we can't modify any messages after they've been published, unless we want our subscribers to get VERY confused
  // new_message.data = "can't do this!";
   }

  ros::Subscriber sub_;
  ros::Publisher pub_;
  ros::Timer timer_;
  };

} // namespace @(namespace)

PLUGINLIB_DECLARE_CLASS(@(package), @(NodeletClass), @(namespace)::@(NamespaceClass), nodelet::Nodelet);






//////////////////////////////////////////////////////////////////////

// NODE CODE

////////////////////////////////////////////////////////////////////////////

#include "ros/ros.h"
#include "nodelet/loader.h"

int main(int argc, char **argv){
  ros::init(argc, argv, "@(node)");
  nodelet::Loader nodelet;
  nodelet::M_string remap(ros::names::getRemappings());
  nodelet::V_string nargv;
  std::string nodelet_name = ros::this_node::getName();
  nodelet.load(nodelet_name, "@(package)/@(nodelet)", remap, nargv);
  ros::spin();
  return 0;
  }
