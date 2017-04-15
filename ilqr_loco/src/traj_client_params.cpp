#include "traj_client.h"

void TrajClient::LoadParams()
{
  try
  {
    ROS_INFO("Loading parameters.");

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

    nh_.getParam("ilqr_tolFun", ilqr_tolFun_);
    nh_.getParam("ilqr_tolConstraint", ilqr_tolConstraint_);
    nh_.getParam("ilqr_tolGrad", ilqr_tolGrad_);
    nh_.getParam("ilqr_max_iter", ilqr_max_iter_);
    nh_.getParam("ilqr_regType", ilqr_regType_);
    nh_.getParam("ilqr_debug_level", ilqr_debug_level_);

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
  nh_.getParam("Opt_car_param/limSteer", limSteer_);
  nh_.getParam("Opt_car_param/limThr", limThr_);

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

// changed to pass by reference to apply and keep edit
void TrajClient::SetOptParams(tOptSet *o)
{
    o->tolFun= ilqr_tolFun_;
    o->tolConstraint= ilqr_tolConstraint_;
    o->tolGrad= ilqr_tolGrad_;
    o->max_iter= ilqr_max_iter_;
    o->regType= ilqr_regType_;
    o->debug_level= ilqr_debug_level_;

    // double default_alpha[]= {1.0, 0.3727594, 0.1389495, 0.0517947, 0.0193070, 0.0071969, 0.0026827, 0.0010000};
    // o->alpha= default_alpha;
    o->n_alpha= 8;
    o->lambdaInit= 1;
    o->dlambdaInit= 1;
    o->lambdaFactor= 1.6;
    o->lambdaMax= 0.0000000001;
    o->lambdaMin= 0.000001;
    o->zMin= 0.0;
    o->w_pen_init_l= 1.0;
    o->w_pen_init_f= 1.0;
    o->w_pen_max_l= INF;
    o->w_pen_max_f= INF;
    o->w_pen_fact1= 4.0; // 4...10 Bertsekas p. 123
    o->w_pen_fact2= 1.0;
}

void TrajClient::LoadOpt()
{
  LoadCarParams();
  LoadCostParams();

  Opt = INIT_OPTSET;

  SetOptParams(&Opt);

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

  char *err_msg;
  double max_iter = 100;

  err_msg = setOptParam(&Opt, "max_iter", &max_iter, 1);
  if(err_msg) {
      printf("Dimagree error, Error setting optimization parameter '%s': %s.\n", "max_iter", err_msg);
  }
}
