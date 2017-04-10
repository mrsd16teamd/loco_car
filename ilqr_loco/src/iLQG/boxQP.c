// Box constrained quadratic optimizer c implementation of http://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization by Yuval Tassa
// Copyright (c) 2016 Jens Geisler
//
// BIBTeX:
// @INPROCEEDINGS{
// author={Tassa, Y. and Mansard, N. and Todorov, E.},
// booktitle={Robotics and Automation (ICRA), 2014 IEEE International Conference on},
// title={Control-Limited Differential Dynamic Programming},
// year={2014}, month={May}, doi={10.1109/ICRA.2014.6907001}}

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "iLQG.h"
#include "boxQP.h"
#include "cholesky.h"
#include "matMult.h"
#include "printMat.h"

#ifndef MOD_CHOL
#define MOD_CHOL 0
#endif

#ifndef DEBUG_BOXQP
#define DEBUG_BOXQP 0
#else
    #if PREFIX1(DEBUG_BOXQP)==1
    #define DEBUG_BOXQP 1
    #endif
#endif

#define TRACE(x) do { if (DEBUG_BOXQP) PRNT x; } while (0)
#define printVec_(x) do { if (DEBUG_BOXQP) printVec x; } while (0)
#define printMat_(x) do { if (DEBUG_BOXQP) printMat x; } while (0)
   


int boxQP(double *H, const double *g, const double *lower, const double *upper, double *x, double *Hfree, double *L, double *grad, double *grad_clamped, double *search, int *is_clamped, int *n_free_, double *invHfree, const int n) {
    int i, j, iter, i_free, j_free;
    double oldvalue= 0.0;
    int res= 0;
    double gnorm= 0.0;
    int nfactor= 0;
    double trace= 0.0;
    double value, d1, mod_chol_regu, sdotg;
    int all_clamped, clamps_changed, n_free, was_clamped;
    double step, vc;
    int nstep;
    double *xc;
    
    const int maxIter           = 100;    // maximum number of iterations
    const double minGrad        = 1e-8;   // minimum norm of non-fixed gradient
    const double minRelImprove  = 1e-8;   // minimum relative improvement
    const double stepDec        = 0.6;    // factor for decreasing stepsize
    const double minStep        = 1e-22;  // minimal stepsize for linesearch
    const double Armijo         = 0.1;    // Armijo parameter (fraction of linear improvement required)
    
    memset(Hfree, 0, sizeof(double)*(n*(n+1))/2);

    for(i= 0; i<n; i++) {
        // clamp to limits
        if(x[i]>upper[i]) x[i]= upper[i];
        if(x[i]<lower[i]) x[i]= lower[i];

        is_clamped[i]= 0; // initial values
    }

#if MOD_CHOL
    if(mod_chol(H, n, L, (int *)grad_clamped /* reuse of memory */, search /* reuse of memory */)>0.0)
        perm_tri_square(L, H, (int *)grad_clamped, n);
#endif
    
    // initial objective value
    value= 0.0;
    for(i= 0; i<n; i++) {
        d1= 0.0;
        for(j= 0; j<n; j++)
            d1+= H[SYMTRI_MAT_IDX(i, j)]*x[j]; // TODO: in case of MOD_CHOL save some cycles and use L instead of H
        
        value+= x[i]*(g[i] + 0.5*d1);
    }
    
    for(iter= 0; iter<maxIter; iter++) {
        if(iter>0 && (oldvalue-value) < minRelImprove*fabs(oldvalue))
            return 4;

        oldvalue= value;
        printVec_((&oldvalue, 1, "oldvalue"));
        
        all_clamped= 1;
        clamps_changed= 0;
        n_free= 0;
        gnorm= 0.0;
        for(i= 0; i<n; i++) {
            // get gradient
            // grad     = g + H*x;
            d1= 0.0;
            for(j= 0; j<n; j++)
                d1+= H[SYMTRI_MAT_IDX(i, j)]*x[j];

            grad[i]= g[i] + d1;

            was_clamped= is_clamped[i];
            if(x[i]<=lower[i] && grad[i]>0)
                is_clamped[i]= 1;
            else if(x[i]>=upper[i] && grad[i]<0)
                is_clamped[i]= 2;
            else {
                is_clamped[i]= 0;
                all_clamped= 0;
                gnorm+= grad[i]*grad[i];
                n_free++;
            }
            if((!was_clamped) != (!is_clamped[i]))
                clamps_changed= 1;
        }
        n_free_[0]= n_free;
        
//         for(i= 0; i<n; i++) search[i]= !is_clamped[i];
//         printVec_((search, n, "free"));
        
        TRACE(("gnorm= %g\n", sqrt(gnorm)));
        
        if(all_clamped)
            return 6;
        
        // [Hfree, indef]  = chol(H(free,free));
        if(iter==0 || clamps_changed) {
            for(j= 0, j_free= 0; j<n; j++) { // cols
                if(!is_clamped[j]) {
                    for(i= 0, i_free= 0; i<=j; i++) { // rows
                        if(!is_clamped[i]) {
                            Hfree[UTRI_MAT_IDX(i_free, j_free)]= H[UTRI_MAT_IDX(i, j)];
                            i_free++;
                        }
                    }
                    j_free++;
                }
            }
            if(!cholesky_tri(Hfree, n_free, L)) {
                return -1;
            }
            cholesky_tri_inv(L, invHfree, n_free, search /* reuse of memory */);
            nfactor++;
        }
        
        // gnorm= sqrt(gnorm);
        if(gnorm<minGrad*minGrad)
            return 5;
        

        // get search direction
        // grad_clamped   = g  + H*(x.*clamped);
        for(i= 0, i_free= 0; i<n; i++) {
            if(!is_clamped[i]) {
                d1= 0.0;
                for(j= 0;  j<n; j++) {
                    if(is_clamped[j]) {
                        d1+= H[SYMTRI_MAT_IDX(i, j)] * x[j];
                    }
                }
                grad_clamped[i_free]= g[i] + d1;
                i_free++;
            }
        }
        // search(free)   = -Hfree\(Hfree'\grad_clamped(free)) - x(free);
        for(i= 0, i_free= 0; i<n; i++) {
            if(!is_clamped[i]) {
                search[i]= -x[i];
                for(j_free= 0;  j_free<n_free; j_free++)
                    search[i]-= invHfree[SYMTRI_MAT_IDX(i_free, j_free)] * grad_clamped[j_free];

                i_free++;
            } else
                search[i]= 0.0;
        }
//         cholesky_solve_tri(Hfree, grad_clamped, search, n_free);
//         for(i= n-1, i_free= n_free-1; i>=0; i--) {
//             if(!is_clamped[i]) {
//                 search[i]= -search[i_free] - x[i];
//                 i_free--;
//             } else
//                 search[i]= 0.0;
//         }
        
        // check for descent direction
        // sdotg          = sum(search.*grad);
        sdotg= 0.0;
        for(i= 0; i<n; i++)
            sdotg+= search[i]*grad[i];
        
        if(sdotg>=0.0) {
            printTri(H, n, "H");
            return -2;
        }
        
        // armijo linesearch
        step= 1.0;
        nstep= 0;
        xc= grad; // reuse memory
        while(1) {
            // xc    = clamp(x+step*search);
            for(i= 0; i<n; i++) {
                xc[i]= x[i] + step*search[i];
                if(xc[i]>upper[i]) xc[i]= upper[i];
                if(xc[i]<lower[i]) xc[i]= lower[i];
            }
            // vc    = xc'*g + 0.5*xc'*H*xc;
            vc= 0.0;
            for(i= 0; i<n; i++) {
                d1= 0.0;
                for(j= 0; j<n; j++)
                    d1+= H[SYMTRI_MAT_IDX(i, j)]*xc[j];

                vc+= xc[i]*(g[i] + 0.5*d1);
            }
            
            if(((vc-oldvalue)/(step*sdotg)) >= Armijo)
                break;
            
            step= step*stepDec;
            if(step<minStep)
                return 2;
            
            nstep++;
        }
        
        // accept candidate
        for(i= 0; i<n; i++)
            x[i]= xc[i];
        value= vc;
        
        TRACE(("iter %-3d  value % -9.5g |g| %-9.3g  reduction %-9.3g  steps %-2d  n_clamped %d nfactor %d\n", iter+1, vc, sqrt(gnorm), oldvalue-vc, nstep, n-n_free, nfactor));
    }
    
    return 1;
}
