#ifndef _ILQR_PLANNER_H_
#define _ILQR_PLANNER_H_

#include <ilqr_loco/TrajExecAction.h>

// #include "iLQR_mpc.c" //TODO integrate generated c-code
//iLQR_mpc(double x_cur[10], double x_des[6], double obs[2], int T);

// Note that the inputs to this function can be whatever is convenient for client
ilqr_loco::TrajExecGoal iLQR_gen_traj(x_current, x_desired, obstacle_pos, T)
{
  // TODO pre-process inputs as necessary, put into C-style arrays
  double x_cur[10] = { }
  double x_des[6]  = { }
  double obs[2]    = { }

  // TODO what is the output type?
  sometype outputs = iLQR_mpc(x_cur, x_des, obs, T);

  ilqr_loco::TrajExecGoal goal;
  // TODO put output from c-code into action message

  return goal;
}

#endif 
