#include "traj_client.h"

TrajClient::TrajClient(): ac_("traj_executer", true), mode_(0), T_(0),
                          cur_integral_(0), prev_error_(0), cur_vel_(0.5)
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

	LoadParams();
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
    nh_.getParam("timeout_ilqr_mpc", mpc_timeout_);
    nh_.getParam("stop_goal_threshold", goal_threshold_);

    LoadOpt();
  }
  catch(...)
  {
    ROS_ERROR("Please put all params into yaml file, and load it.");
  }
}

void TrajClient::LoadCarParams()
{
  nh_.getParam("Opt_car_param/g", g_);
  nh_.getParam("Opt_car_param/L", L_);
  nh_.getParam("Opt_car_param/m", m_);
  nh_.getParam("Opt_car_param/b", b_);
  nh_.getParam("Opt_car_param/c_x", c_x_);
  nh_.getParam("Opt_car_param/c_a", c_a_);
  nh_.getParam("Opt_car_param/Iz", Iz_);
  nh_.getParam("Opt_car_param/mu", mu_);
  nh_.getParam("Opt_car_param/mu_s", mu_s_);

  a_ = L_ - b_;
  G_f_ = m_*g_*b_/L_;
  G_r_ = m_*g_*a_/L_;
}

void TrajClient::LoadCostParams()
{
  nh_.getParam("Opt_cost/cu", cu_);
  nh_.getParam("Opt_cost/cdu", cdu_);
  nh_.getParam("Opt_cost/cf", cf_);
  nh_.getParam("Opt_cost/pf", pf_);
  nh_.getParam("Opt_cost/cx", cx_);
  nh_.getParam("Opt_cost/cdx", cdx_);
  nh_.getParam("Opt_cost/px", px_);
  nh_.getParam("Opt_cost/cdrift", cdrift_);
  nh_.getParam("Opt_cost/k_pos", k_pos_);
  nh_.getParam("Opt_cost/k_vel", k_vel_);
  nh_.getParam("Opt_cost/d_thres", d_thres_);
}

void TrajClient::LoadOpt()
{
  LoadCarParams();
  LoadCostParams();

  Opt = INIT_OPTSET;
  standard_parameters(&Opt);
  Opt.p= (double **) malloc(n_params*sizeof(double *));
  Opt.p[0] = assignPtrVal(&G_f_,1);
  Opt.p[1] = assignPtrVal(&G_r_,1);;
  Opt.p[2] = assignPtrVal(&Iz_,1);;
  // [3] Obs
  Opt.p[4] = assignPtrVal(&a_,1);
  Opt.p[5] = assignPtrVal(&b_,1);
  Opt.p[6] = assignPtrVal(&c_a_,1);
  Opt.p[7] = assignPtrVal(&c_x_,1);
  Opt.p[8] = assignPtrVal(&cdrift_,1);
  Opt.p[9] = assignPtrVal(&cdu_[0],2);
  Opt.p[10] = assignPtrVal(&cdx_[0],3);
  Opt.p[11] = assignPtrVal(&cf_[0],6);
  Opt.p[12] = assignPtrVal(&cu_[0],2);
  Opt.p[13] = assignPtrVal(&cx_[0],3);
  Opt.p[14] = assignPtrVal(&d_thres_,1);
  Opt.p[15] = assignPtrVal(&timestep_,1);
  Opt.p[16] = assignPtrVal(&k_pos_,1);
  Opt.p[17] = assignPtrVal(&k_vel_,1);
  Opt.p[18] = assignPtrVal(&limSteer_[0],2);
  Opt.p[19] = assignPtrVal(&limThr_[0],2);
  Opt.p[20] = assignPtrVal(&m_,1);
  Opt.p[21] = assignPtrVal(&mu_,1);
  Opt.p[22] = assignPtrVal(&mu_s_,1);
  Opt.p[23] = assignPtrVal(&pf_[0],6);
  Opt.p[24] = assignPtrVal(&px_[0],3);
  // [25] xDes

  char *err_msg, *fname;
  fname = "max_iter";
  double max_iter = 100;

  err_msg = setOptParam(&Opt, fname, &max_iter, 1);
  if(err_msg) {
      printf("Dimagree error, Error setting optimization parameter '%s': %s.\n", fname, err_msg);
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

  if (mode_==1 || (mode_==3 && !obs_received_) || (mode_==4 && !obs_received_) )
  {
    rampPlan();
  }
}

void TrajClient::obsCb(const geometry_msgs::PointStamped &msg)
{
  if (msg.point.x < 100) //Note: This is just cuz detector pubs 999 when it sees nothing
  {
    obs_pos_.x = msg.point.x;
    obs_pos_.y = msg.point.y;
    ROS_INFO("Received obstacle message: x = %f, y = %f", obs_pos_.x, obs_pos_.y);
    obs_received_ = true;

    if (mode_==1){
      SendZeroCommand();
    }
    else if (mode_==2 || mode_==3){
      ilqrPlan();
      mode_ = 0;
    }
    else if (mode_==4){
      ilqrMPC();
    }
  }
}

void TrajClient::modeCb(const geometry_msgs::Point &msg)
{
  int command = msg.x;

  if (!state_estimate_received_){
    ROS_INFO("Haven't received state info yet.");
    return;
  }

  switch (command)
  {
    //DONT CHANGE THESE! TOO MUCH WORK
    case 1: {
      ROS_INFO("Mode 1: ramp. If I see an obstacle, I'll brake!");
      mode_ = 1;
      start_time_ = ros::Time::now();
      // wait for stateCb to ramp
      break;
    }
    case 2: {
      ROS_INFO("Mode 2: iLQR from static initial conditions.");
      mode_ = 2;
		  // wait for obsCb to plan
      break;
    }
    case 3: {
      ROS_INFO("Mode 3: ramp -> iLQR open loop.");
      mode_ = 3;
      obs_received_ = false;
      start_time_ = ros::Time::now();
      break;
      //wait for stateCb to ramp
    }
    case 4: {
      ROS_INFO("Mode 4: ramp -> iLQR closed loop.");
      mode_ = 4;
      obs_received_ = false;
      //wait for stateCb to ramp
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
