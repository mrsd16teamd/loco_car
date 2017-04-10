#ifndef ILQG_MPC_H
#define ILQG_MPC_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "iLQG.h"
#include "printMat.h"
#include "matMult.h"

struct trajectory {
    double* x;  // states
    double* u;  // controls
};

double* assignPtrVal(double* values, int numVal);
void init_params(tOptSet *o, double* xDes, double* Obs);
void plan_trajectory(double* x0, double* u0, double* xDes, double* Obs, int T,  tOptSet *o, struct trajectory* Traj);

#endif
