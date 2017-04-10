// MATLAB Mex function wrapper for iLQG algorithm
// Copyright (c) 2016 Jens Geisler
#include "iLQG_plan.h"

double* assignPtrVal(double* values, int numVal) {
    double* temp = (double *) malloc(numVal*sizeof(double));
    memcpy(temp, values, numVal*sizeof(double));

    return temp;
}

void init_params(tOptSet *o, double* xDes, double* Obs)
{
    o->p[3] = assignPtrVal(Obs,2);
    o->p[25] = assignPtrVal(xDes,6);
}

void plan_trajectory(double* x0, double* u0, double* xDes, double* Obs, int T, tOptSet *o, struct trajectory* Traj)
{
    // dims
    int N, n, m, m_, n_, si, i, k;
    // inputs
    tOptSet o= INIT_OPTSET;
    int dims[3];

    //outputs
    double *x_nom, *u_nom; //, *l, *L;

    // aux
    char *err_msg, *fname;
    clock_t begin, end;

    n= sizeof(x0)/sizeof(x0[0]);  // length of state vector
    m= 2; // number of inputs
    N= T+1;

    // inputs
    o->x0= x0; //double *
    u_nom= u0;  // double **
    o->n_hor= N-1;

    // Set optimization parameters
    // standard_parameters(o);
    // fname = "max_iter";
    // double max_iter = 100;

    // err_msg = setOptParam(o, fname, &max_iter, 1);
    // if(err_msg) {
    //     printf("Dimagree error, Error setting optimization parameter '%s': %s.\n", fname, err_msg);
    // }

    // Set model and problem parameters
    init_params(o, xDes, Obs);

		// TODO store o in traj_client so we can skip above

    // outputs
    double success[1];
    double new_cost[1];

    // aux
    for(i= 0; i<NUMBER_OF_THREADS+1; i++)
        o->trajectories[i].t= (trajEl_t *) malloc(sizeof(trajEl_t)*(N-1));
    o->multipliers.t= (multipliersEl_t *) malloc(sizeof(multipliersEl_t)*N);

    printf("Set const vars\n");
    if(!init_opt(o)) {
        success[0]= 0;
        new_cost[0]= o->cost;
    } else {
        printf("Initializing trajectory\n");
        for(k= 0; k<N-1; k++)
            for(i= 0; i<N_U; i++)
                o->nominal->t[k].u[i]= u_nom[MAT_IDX(i, k, N_U)];
        if(!forward_pass(o->candidates[0], o, 0.0, &o->cost, 0)) {
            printf("forward_pass failed\n");
            success[0]= 0;
            new_cost[0]= o->cost;
        } else {
            makeCandidateNominal(o, 0);

            printf("Starting iLQG\n");
            begin = clock();
            success[0]= iLQG(o);
            end = clock();
            printf("Time for iLQG: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
            for(k= 0; k<N-1; k++)
                for(i= 0; i<N_X; i++)
                    Traj->x[MAT_IDX(i, k, N_X)]= o->nominal->t[k].x[i];
            for(i= 0; i<N_X; i++)
                Traj->x[MAT_IDX(i, N-1, N_X)]= o->nominal->f.x[i];

            for(k= 0; k<N-1; k++)
                for(i= 0; i<N_U; i++)
                    Traj->u[MAT_IDX(i, k, N_U)]= o->nominal->t[k].u[i];

            new_cost[0]= o->cost;
        }
    }
    printf("runs to here");

    for(i= 0; i<NUMBER_OF_THREADS+1; i++)
        free(o->trajectories[i].t);
}
