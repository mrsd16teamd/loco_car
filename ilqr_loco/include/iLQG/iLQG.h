#ifndef ILQG_H
#define ILQG_H

#ifndef FULL_DDP
    #define FULL_DDP 1
#endif

#include "iLQG_problem.h"

#define INIT_OPTSET {0, 0, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, NULL, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0, NULL, NULL, NULL, {0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, NULL}

#ifndef PRNT
#define PRNT printf
#endif

#ifndef MULTI_THREADED
#define MULTI_THREADED 0
#endif
#if MULTI_THREADED
#include <pthread.h>
#endif

#ifndef NUMBER_OF_THREADS
#define NUMBER_OF_THREADS 1
#endif

#define DO_PREFIX1(VAL)  1 ## VAL
#define PREFIX1(VAL)     DO_PREFIX1(VAL)


typedef struct paramDesc {
  char *name;
  int size;
  int is_var;
} tParamDesc;

typedef struct optSet {
    int n_hor;
    int debug_level;
    double *x0, new_cost, cost, dcost, lambda, g_norm, expected;
    double **p;
    const double *alpha;
    int n_alpha;
    double lambdaMax;
    double lambdaMin;
    double lambdaInit;
    double dlambdaInit;
    double lambdaFactor;
    int max_iter;
    double tolGrad;
    double tolFun;
    double tolConstraint;
    double zMin;
    int regType;
    int iterations;
    int *log_linesearch;
    double *log_z;
    double *log_cost;
    double dV[2];

    double w_pen_l;
    double w_pen_f;
    double w_pen_max_l;
    double w_pen_max_f;
    double w_pen_init_l;
    double w_pen_init_f;
    double w_pen_fact1;
    double w_pen_fact2;

    traj_t *nominal;
    traj_t *candidates[NUMBER_OF_THREADS];

    traj_t trajectories[NUMBER_OF_THREADS+1];

    multipliers_t multipliers;
} tOptSet;

void printParams(double **p, int k);
void standard_parameters(tOptSet *o);
int iLQG(tOptSet *o);
char *setOptParam(tOptSet *o, const char *name, const double *value, const int n);
int forward_pass(traj_t *c, tOptSet *o, double alpha, double *csum, int cost_only);
void makeCandidateNominal(tOptSet *o, int idx);
int calc_derivs(tOptSet *o);
int init_opt(tOptSet *o);
int update_multipliers(tOptSet *o, int init);
int get_g_size();
int calcG(double g[], trajEl_t *t, int k, double *p[]);

static double max(double a, double b) {
    return (a > b)? a: b;
}

static double min(double a, double b) {
    return (a < b)? a: b;
}
extern int n_params;
extern int n_vars;
extern tParamDesc *paramdesc[];

#if MULTI_THREADED
extern pthread_mutex_t step_mutex;
extern pthread_cond_t  next_step_condition;
extern int step_calc_done;
#endif


#endif /* ILQG_H */
