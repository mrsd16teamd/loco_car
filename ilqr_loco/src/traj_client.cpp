#include "traj_client.h"

#define ILQRDEBUG 0
#define DUMMYOBS obs_pos_.x = 999; obs_pos_.y = 0.365210;
#define DUMMYOBSSTATE {obs_pos_.x = 2.599635; obs_pos_.y = 0.365210; cur_state_.pose.pose.position.x = 1.826; cur_state_.pose.pose.position.y = 0.340; double theta = 0.0032; cur_state_.pose.pose.orientation = tf::createQuaternionMsgFromYaw(theta); cur_state_.twist.twist.linear.x = 0.062; cur_state_.twist.twist.linear.y = -0.009; cur_state_.twist.twist.angular.z = 0.00023;}


TrajClient::TrajClient(): ac_("traj_server", true), mode_(0), T_(0),
                          cur_integral_(0), prev_error_(0), step_on_last_traj_(0),
                          tf_listener_(nh, ros::Duration(5), true)
{
  state_sub_  = nh.subscribe("odometry/filtered", 1, &TrajClient::stateCb, this);
  mode_sub_ = nh.subscribe("client_command", 1, &TrajClient::modeCb, this);
  scan_sub_ = nh.subscribe("scan", 1, &TrajClient::scanCb, this);

  predicted_state_pub_ = nh.advertise<nav_msgs::Odometry>("odometry/predicted", 1);
  obs_pos_pub_ = nh.advertise<geometry_msgs::PointStamped>("cluster_center", 1);

  state_estimate_received_ = false;
  found_obstacle_ = false;
  ramp_goal_flag_ = false;
  reacted_to_obstacle_ = false;

  LoadParams();
  InitDetector();
  FillInitControlSeq();

  ROS_INFO("Waiting for action server to start.");
  ac_.waitForServer(); //will wait for infinite time
  ROS_INFO("Action client started. Send me commands from keyboard_command!");
}

void TrajClient::stateCb(const nav_msgs::Odometry &msg)
{
  state_estimate_received_ = true;

  prev_state_ = cur_state_;
  cur_state_ = msg;

  if (T_ == 0)
  {
    cur_integral_ = 0;
    prev_error_ = 0;
    start_state_ = cur_state_;
    ramp_start_y_ = start_state_.pose.pose.position.y;
    u_seq_saved_ = init_control_seq_;
    start_time_ = ros::Time::now();
  }

  if (mode_==1 || ((mode_==5 || mode_==6 || mode_==7 || mode_==13 || mode_==14) && !found_obstacle_))
    rampPlan();
}

void TrajClient::ReactToObstacle()
{
    ROS_INFO("Reacting to obstacle!");
    if (mode_ == 1)
      SendZeroCommand(); //brake
    else if (mode_==2 || mode_==5)
      Plan();
    else if (mode_==3 || mode_==6)
      MpcILQR();
    else if (mode_==7 || mode_==13 || mode_==14 || mode_==15) {
      SendInitControlSeq();
      if (mode_==13)
        Plan();
      else if (mode_==14)
        MpcILQR();
    }

    reacted_to_obstacle_ = true;
    mode_ = 0;
}

void TrajClient::modeCb(const geometry_msgs::Point &msg)
{
  int command = msg.x;

  #if ILQRDEBUG
  state_estimate_received_ = true;
  #endif

  if (!state_estimate_received_)
  {
    ROS_INFO("Haven't received state info yet.");
    return;
  }

  T_ = 0;

  switch (command)
  {
    case 1: ROS_INFO("Mode 1: ramp. If I see an obstacle, I'll brake!");
			mode_ = 1;
            break;
            // wait for stateCb to ramp
    case 2: ROS_INFO("Mode 2:iLQR open-loop from static.");
            mode_ = 2;
            #if ILQRDEBUG
                DUMMYOBS
                u_seq_saved_ = init_control_seq_;
                PlanFromCurrentStateILQR();
                mode_ = 0;
            #endif
            break;
            // wait for obsCb to plan
    case 3: ROS_INFO("Mode 3: iLQR mpc from static");
            mode_ = 3;
            #if ILQRDEBUG
                DUMMYOBS
                MpcILQR();
                mode_ = 0;
            #endif
            break;
            //wait for stateCb to ramp
    case 5: ROS_INFO("Mode 5: ramp -> iLQR open-loop.");
            mode_ = 5;
            break;
            //wait for obsCb to plan
    case 6: ROS_INFO("Mode 6: ramp -> iLQR mpc");
            mode_ = 6;
            break;
	  case 7: ROS_INFO("Mode 7: ramp -> playback");
            mode_ = 7;
            break;
    case 8: ROS_INFO("Resetting obstacle.");
			      ResetObstacle();
            break;
    case 10: ROS_INFO("Play back initial control sequence");
			       mode_ = 10;
             SendInitControlSeq();
             break;
    case 12: ROS_INFO("Placing fake obstacle %f meters ahead.", obs_dist_thres_);
             InsertFakeObs();
             break;
    case 13: ROS_INFO("ramp -> playback -> iLQR open-loop.");
			 mode_ = 13;
             break;
    case 14: ROS_INFO("ramp -> playback -> iLQR mpc");
    	       mode_ = 14;
             break;

    default: ROS_INFO("Please enter valid command.");
  }
}

void TrajClient::feedbackCb(const ilqr_loco::TrajExecFeedbackConstPtr& feedback)
{
  // Keeps track of progress of TrajAction server along most recently sent trajectory
  step_on_last_traj_ = feedback->step;
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "traj_client");
  TrajClient client;
  ros::spin();

  return 0;
}
