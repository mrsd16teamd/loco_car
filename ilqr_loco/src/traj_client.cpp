#include "traj_client.h"

// Constructor
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
    nh_.getParam("T_horizon", T_horizon_);
    nh_.getParam("init_control_seq", init_control_seq_);
    nh_.getParam("X_des", x_des_);

    nh_.getParam("timestep", timestep_);

    nh_.getParam("kp_ramp", kp_);
    nh_.getParam("ki_ramp", ki_);
    nh_.getParam("kd_ramp", kd_);
    nh_.getParam("accel_ramp", accel_);
    nh_.getParam("target_vel_ramp", target_vel_);
    nh_.getParam("timeout_ramp", timeout_);
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

//////////////////////// GENERATE RAMP START //////////////////////////////////

ilqr_loco::TrajExecGoal TrajClient::rampGenerateTrajectory(nav_msgs::Odometry prev_state,
                                                           nav_msgs::Odometry cur_state) {

  ilqr_loco::TrajExecGoal goal;
  FillGoalMsgHeader(goal);

  double dt = (cur_state.header.stamp).toSec() - (prev_state.header.stamp).toSec();
  double yaw = tf::getYaw(cur_state.pose.pose.orientation);

  // PID control for vehicle heading
  double error = 0 - yaw;
  cur_integral_ += error*dt;
  double output = kp_*error + (std::abs(cur_integral_)<1 ? ki_*cur_integral_ : 0) + (dt>0.01 ? kd_*(error-prev_error_)/dt : 0);
  prev_error_ = error;

  // Generate goal
  cur_vel_ += accel_*dt;
  double v = cur_state.twist.twist.linear.x + accel_*dt + 0.75;
  v = cur_vel_<target_vel_ ? cur_vel_ : target_vel_;

  geometry_msgs::Twist control_msg;
  FillTwistMsg(control_msg, v, output);
  goal.traj.commands.push_back(control_msg);

  nav_msgs::Odometry state_msg;
  double expected_x = start_state_.pose.pose.position.x + 0.5*accel_*start_time_.toSec()*start_time_.toSec();
  double expected_y = start_state_.pose.pose.position.y;

  FillOdomMsg(state_msg, expected_x, expected_y, 0, v, 0, 0); //0s: yaw, vy, w
  goal.traj.states.push_back(state_msg);

  ++T_;
  ramp_goal_flag_ = v>=target_vel_ ? true : false;  // Ramp completion flag

  return goal;
}


void TrajClient::rampPlan() {

  if(ros::Time::now() - start_time_ < ros::Duration(timeout_) && (mode_==1 || (mode_==3 && !obs_received_))) {
    ilqr_loco::TrajExecGoal goal = rampGenerateTrajectory(prev_state_, cur_state_);
    SendTrajectory(goal);
  }
  else {
    // Stop car after ramp timeout
    ROS_INFO("Timeout exceeded, stopping car");
    ilqr_loco::TrajExecGoal end_goal;

    geometry_msgs::Twist control_msg;
    FillTwistMsg(control_msg, 0, 0);

    end_goal.traj.commands.push_back(control_msg);
    end_goal.traj.states.push_back(cur_state_);
    SendTrajectory(end_goal);

    obs_received_ = true;
    mode_ = 0;
  }
}

//////////////////////// GENERATE RAMP END //////////////////////////////////

//////////////////////// GENERATE iLQG START //////////////////////////////////

// Calls iLQG_mpc.c to generate new trajectory
ilqr_loco::TrajExecGoal TrajClient::ilqgGenerateTrajectory(nav_msgs::Odometry cur_state)
{
  ROS_INFO("Generating iLQG trajectory.");
  ilqr_loco::TrajExecGoal goal;
  FillGoalMsgHeader(goal);
  goal.traj.timestep = timestep_;

  double theta = tf::getYaw(cur_state.pose.pose.orientation);

  ROS_INFO("Start state: %f, %f, %f, %f, %f, %f",
            cur_state.pose.pose.position.x, cur_state.pose.pose.position.y, theta,
            cur_state.twist.twist.linear.x, cur_state.twist.twist.linear.y, cur_state.twist.twist.angular.z);
  ROS_INFO("Obs pos: %f, %f", obs_pos_.x, obs_pos_.y);

  iLQR_gen_traj(cur_state, init_control_seq_, x_des_, obs_pos_, T_horizon_, goal);
  ++T_;

  return goal;
}

void TrajClient::ilqrPlan()
{
  ilqr_loco::TrajExecGoal goal = ilqgGenerateTrajectory(cur_state_);
  SendTrajectory(goal);
}

//////////////////////// Message utils //////////////////////////////////

void TrajClient::FillGoalMsgHeader(ilqr_loco::TrajExecGoal &goal)
{
  goal.traj.header.seq = T_;
  goal.traj.header.stamp = ros::Time::now();
  goal.traj.header.frame_id = "/base_link";
}

void TrajClient::FillTwistMsg(geometry_msgs::Twist &twist, double lin_x, double ang_z)
{
  twist.linear.x = lin_x;
  twist.angular.z = ang_z;
}

void TrajClient::FillOdomMsg(nav_msgs::Odometry &odom, double x, double y,
                             double yaw, double Ux, double Uy, double w)
{
  odom.pose.pose.position.x = x;
  odom.pose.pose.position.y = y;
  odom.pose.pose.position.z = 0.0;

  geometry_msgs::Quaternion odom_quat = tf::createQuaternionMsgFromYaw(yaw);
  odom.pose.pose.orientation = odom_quat;

  odom.twist.twist.linear.x = Ux;
  odom.twist.twist.linear.y = Uy;
  odom.twist.twist.angular.z = w;
}

//////////////////////// GENERATE iLQG END //////////////////////////////////


int main(int argc, char** argv)
{
  ros::init(argc, argv, "test_client");
  TrajClient client;
  ros::spin();

  return 0;
}
