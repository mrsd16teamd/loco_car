#include "traj_client.h"

#define ILQRDEBUG 0

TrajClient::TrajClient(): ac_("traj_server", true), mode_(0), T_(0),
                          cur_integral_(0), prev_error_(0), cur_vel_(0.5)
{
  state_sub_  = nh.subscribe("odometry/filtered", 1, &TrajClient::stateCb, this);
  obs_sub_ = nh.subscribe("cluster_center", 1, &TrajClient::obsCb, this);
  mode_sub_ = nh.subscribe("client_command", 1, &TrajClient::modeCb, this);

  state_estimate_received_ = false;
  obs_received_ = false;
  ramp_goal_flag_ = false;

  LoadParams();

  ROS_INFO("Waiting for action server to start.");
  ac_.waitForServer(); //will wait for infinite time
  ROS_INFO("Action client started. Send me commands from keyboard_command!");
}

void TrajClient::stateCb(const nav_msgs::Odometry &msg)
{
  prev_state_ = cur_state_;
  cur_state_ = msg;

  state_estimate_received_ = true;

  if (T_ == 0) {
    cur_integral_ = 0;
    prev_error_ = 0;
    start_state_ = cur_state_;
	  ramp_start_y_ = start_state_.pose.pose.position.y;
  }

  if (mode_==1 || (mode_==3 && !obs_received_) || (mode_==4 && !obs_received_) )
  {
    rampPlan();
  }
}

void TrajClient::obsCb(const geometry_msgs::PointStamped &msg)
{
  if (msg.point.x < 100) {
    ROS_INFO("Received obstacle message.");
    obs_pos_.x = msg.point.x;
    obs_pos_.y = msg.point.y;
    obs_received_ = true;

    if (mode_==1) {
      SendZeroCommand();
      mode_ = 0;
    }
    else if (mode_==2 || mode_==3 || mode_==7){
      ilqrPlan();
      mode_ = 0;
    }
    else if (mode_==4 || mode_==5){
      ilqrMPC();
    }
    else if (mode_==6){
      ilqrSparseReplan();
    }
  }
}

void TrajClient::modeCb(const geometry_msgs::Point &msg)
{
  int command = msg.x;

  #if ILQRDEBUG
  state_estimate_received_ = true;
  #endif

  if (!state_estimate_received_){
    ROS_INFO("Haven't received state info yet.");
    return;
  }

  switch (command)
  {
    //DONT CHANGE THESE! TOO MUCH WORK
    case 1: {
      ROS_INFO("Mode 1: ramp. If I see an obstacle, I'll brake!");
      T_ = 0;
      mode_ = 1;
      start_time_ = ros::Time::now();
      break;
      // wait for stateCb to ramp
    }
    case 2: {

      #if ILQRDEBUG
      obs_pos_.x = 2.599635;
      obs_pos_.y = 0.365210;

      cur_state_.pose.pose.position.x = 1.826;
      cur_state_.pose.pose.position.y = 0.340;
      double theta = 0.0032;
      cur_state_.pose.pose.orientation = tf::createQuaternionMsgFromYaw(theta);
      cur_state_.twist.twist.linear.x = 0.062;
      cur_state_.twist.twist.linear.y = -0.009;
      cur_state_.twist.twist.angular.z = 0.00023;
      ilqrPlan();
      #endif

      ROS_INFO("Mode 2: iLQR from static initial conditions.");

      T_ = 0;
      mode_ = 2;
	  u_seq_saved_ = init_control_seq_;
      obs_received_ = false;
      break;
		  // wait for obsCb to plan
    }
    case 3: {
      ROS_INFO("Mode 3: ramp -> iLQR open loop.");
      T_ = 0;
      mode_ = 3;
      obs_received_ = false;
      start_time_ = ros::Time::now();
      break;
      //wait for stateCb to ramp
    }
    case 4: {
      ROS_INFO("Mode 4: ramp -> receding horizon iLQR.");
      T_ = 0;
      mode_ = 4;
      obs_received_ = false;
      start_time_ = ros::Time::now();
      break;
      //wait for stateCb to ramp
    }
    case 5: {
      ROS_INFO("Mode 5: Receding horizon iLQR from static initial conditions.");
      T_ = 0;
      mode_ = 5;
      obs_received_ = false;
      break;
      //wait for obsCb to plan
    }
    case 6: {
      ROS_INFO("Mode 6: iLQR with sparse replanning from static.");

      #if ILQRDEBUG
      obs_pos_.x = 2.599635;
      obs_pos_.y = 0.365210;

      cur_state_.pose.pose.position.x = 1.826;
      cur_state_.pose.pose.position.y = 0.340;
      double theta = 0.0032;
      cur_state_.pose.pose.orientation = tf::createQuaternionMsgFromYaw(theta);
      cur_state_.twist.twist.linear.x = 0.062;
      cur_state_.twist.twist.linear.y = -0.009;
      cur_state_.twist.twist.angular.z = 0.00023;
      ilqrSparseReplan();
      #endif

      T_ = 0;
      mode_ = 6;
      obs_received_ = false;
      break;
    }
    case 7: {
      ROS_INFO("Mode 7: iLQR w/ pid corrections from static initial conditions.");
      mode_ = 7;
	    u_seq_saved_ = init_control_seq_;
      obs_received_ = false;
      break;
    }
    case 8: {
      ROS_INFO("Resetting obs_received_ to false.");
      obs_received_ = false;
      break;
    }
    case 9: {
      ROS_INFO("Killing node.");
      mode_ = 0;
      ros::shutdown();
      break;
    }
    default: {
      ROS_INFO("Please enter valid command.");
    }
  }
}

void TrajClient::SendTrajectory(ilqr_loco::TrajExecGoal &goal)
{
  ROS_INFO("Sending trajectory.");
  ac_.sendGoal(goal);
              //  ,boost::bind(&TrajClient::doneCb, this, _1, _2),
              //  boost::bind(&TrajClient::activeCb, this),
              //  boost::bind(&TrajClient::feedbackCb, this, _1));
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "traj_client");
  TrajClient client;
  ros::spin();

  return 0;
}
