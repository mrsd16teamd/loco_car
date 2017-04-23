/* File generated form template iLQG_problem.tem on 2017-04-23 17:10:08-04:00. Do not edit! */

#ifndef ILQG_PROBLEM_H
#define ILQG_PROBLEM_H

#include <math.h>
// #include "mex.h"
#ifndef  HAVE_OCTAVE
// #include "matrix.h"
#endif

// #define isNANorINF(v) (mxIsNaN(v) || mxIsInf(v))
// #define INF mxGetInf()
#define isNANorINF(v) (isnan(v) || isinf(v))
#define INF INFINITY

#define N_X 9 
#define N_U 2 

#define sizeofQxx 45 
#define sizeofQuu 3 
#define sizeofQxu 18 

typedef struct {
    double x[N_X];
    double u[N_U];
    double lower[N_U];
    double upper[N_U];
    double lower_sign[N_U];
    double upper_sign[N_U];
    double lower_hx[N_X*N_U];
    double upper_hx[N_X*N_U];
    
    double l[N_U];
    double L[N_U*N_X]; 
    double c;
    double cx[N_X];
    double cxx[sizeofQxx];
    double cu[N_U];
    double cuu[sizeofQuu];
    double cxu[sizeofQxu];
    double fx[N_X*N_X];
    double fu[N_X*N_U];
#if FULL_DDP
    double fxx[N_X*sizeofQxx];
    double fuu[N_X*sizeofQuu];
    double fxu[N_X*sizeofQxu];
#endif

} trajEl_t;

typedef struct {
    double x[N_X];

    double c;
    double cx[N_X];
    double cxx[sizeofQxx];


} trajFin_t;

typedef struct {
    trajEl_t* t;
    trajFin_t f;
} traj_t;

typedef struct {

} multipliersEl_t;

typedef struct {

} multipliersFin_t;

typedef struct {
    multipliersEl_t* t;
    multipliersFin_t f;
} multipliers_t;

#endif // ILQG_PROBLEM_H
