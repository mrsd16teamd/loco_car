#include "traj_client.h"

TrajClient::TrajClient(): ac_("traj_executer", true)
{
  state_sub_  = nh_.subscribe("odometry/filtered", 1, &TrajClient::stateCb, this);
  obs_sub_ = nh_.subscribe("cluster_center", 1, &TrajClient::obsCb, this);
  mode_sub_ = nh_.subscribe("client_command", 1, &TrajClient::modeCb, this);

  ROS_INFO("Waiting for action server to start.");
  ac_.waitForServer(); //will wait for infinite time
  ROS_INFO("Action server started. Send me commands from teleop_keyboard!");

  state_estimate_received_ = false;
  obs_received_ = false;
  ramp_goal_flag_ = false;
  mode_ = 0;
  T_ = 0;

  cur_integral_ = 0.0;
  prev_error_ = 0.0;
  cur_vel_ = 0.5;
}

void TrajClient::LoadParams()
{
  try
  {
    // Get parameters from ROS Param server
    nh_.getParam("timestep", timestep_);

    nh_.getParam("kp_ramp", kp_);
    nh_.getParam("ki_ramp", ki_);
    nh_.getParam("kd_ramp", kd_);
    nh_.getParam("accel_ramp", accel_);
    nh_.getParam("target_vel_ramp", target_vel_);
    nh_.getParam("timeout_ramp", timeout_);

    nh_.getParam("T_horizon", T_horizon_);
    nh_.getParam("init_control_seq", init_control_seq_);
    nh_.getParam("X_des", x_des_);
  }
  catch(...)
  {
    ROS_ERROR("Please put all params into yaml file, and load it.");
  }
}

void TrajClient::stateCb(const nav_msgs::Odometry &msg)
{
  prev_state_ = cur_state_;
  cur_state_ = msg;

  //turn velocities into body frame
  double theta = tf::getYaw(cur_state_.pose.pose.orientation);
  double old_vx = cur_state_.twist.twist.linear.x;
  double old_vy = cur_state_.twist.twist.linear.y;
  cur_state_.twist.twist.linear.x = cos(theta)*old_vx + sin(theta)*old_vy;
  cur_state_.twist.twist.linear.y = cos(theta+PI/2)*old_vx + sin(theta+PI/2)*old_vy;

  state_estimate_received_ = true;
  if (mode_==1 || (mode_==3 && !obs_received_) ){
    rampPlan();
  }
}

void TrajClient::obsCb(const geometry_msgs::PointStamped &msg)
{
  if (msg.point.x < 100) //Note: This is just cuz detector pubs 999 when sees nothing
  {
    obs_pos_.x = msg.point.x;
    obs_pos_.y = msg.point.y;
    ROS_INFO("Received obstacle message: x = %f, y = %f", obs_pos_.x, obs_pos_.y);
    obs_received_ = true;
    if (mode_==3){
      ilqrPlan();
    }
  }
}

void TrajClient::modeCb(const geometry_msgs::Point &msg)
{
  int command = msg.x;
  ROS_INFO("Received command %d", command);
  if (!state_estimate_received_){
    ROS_INFO("Haven't received state info yet.");
    return;
  }
  switch (command)
  {
    //DONT CHANGE THESE! TOO MUCH WORK
    case 1: { //ramp
      start_time_ = ros::Time::now();
      mode_ = 1;
      break;
    }
    case 2: { //iLQR static
      if(!obs_received_){
        ROS_INFO("Haven't received obstacle info yet.");
        break;
      }
      mode_=2;
      ilqrPlan();
      break;
    }
    case 3: { //ramp and iLQR open loop
      obs_received_ = false;
      start_time_ = ros::Time::now();
      mode_=3;
      break;
      //wait for next stateCb
    }
    case 4: { //ramp and iLQR closed loop
      //TODO fill this
    }
    case 8: { //reset obs
      obs_received_ = false;
      break;
    }
    case 9: { //kill client
      ROS_INFO("Killing node.");
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
  ros::init(argc, argv, "test_client");
  TrajClient client;
  ros::spin();

  return 0;
}
