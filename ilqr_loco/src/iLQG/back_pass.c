// iLQG back pass c implementation of http://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization by Yuval Tassa
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

//#include "mex.h"
#ifndef  HAVE_OCTAVE
//#include "matrix.h"
#endif

#include "back_pass.h"
#include "iLQG.h"
#include "matMult.h"
#include "boxQP.h"
#include "printMat.h"

#ifndef DEBUG_BACKPASS
#define DEBUG_BACKPASS 0
#else
    #if PREFIX1(DEBUG_BACKPASS)==1
    #define DEBUG_BACKPASS 1
    #endif
#endif

#define TRACE(x) do { if (DEBUG_BACKPASS) PRNT x; } while (0)
#define printVec_(x) do { if (DEBUG_BACKPASS) printVec x; } while (0)
#define printMat_(x) do { if (DEBUG_BACKPASS) printMat x; } while (0)
   
int back_pass(tOptSet *o) {
    double d1, g_norm_i, g_norm_max, g_norm_sum;
    int k, i_, j_, k_, l_, i_free, j_free, k_free, m_free, qpRes;
    int N= o->n_hor;
    double invHfree[sizeofQuu], L[sizeofQuu];
    double R[sizeofQuu], grad[N_U], grad_clamped[N_U], search[N_U];
    int is_clamped[N_U];
    double Vx[N_X], Vxx[sizeofQxx], Qx[N_X], Qu[N_U], Qxx[sizeofQxx];
    double Qxu_reg[sizeofQxu], Qxu[sizeofQxu];
    double QuuF[sizeofQuu], Quu[sizeofQuu];
    double dummy[N_X*N_X];
    trajEl_t *t= o->nominal->t + N - 1;
    trajFin_t *f= &o->nominal->f;
    
    g_norm_sum= 0.0;

    o->dV[0]= 0.0;
    o->dV[1]= 0.0;
    
#if MULTI_THREADED  
    pthread_mutex_lock(&step_mutex);
        while(step_calc_done>N)
            pthread_cond_wait(&next_step_condition, &step_mutex);
    pthread_mutex_unlock(&step_mutex);
    if(step_calc_done<0)
        return 2;
#endif

    memcpy(Vx, f->cx, sizeof(double)*N_X);
    memcpy(Vxx, f->cxx, sizeof(double)*sizeofQxx);

    for(k= N-1; k>=0; k--, t--) {
#if MULTI_THREADED  
        pthread_mutex_lock(&step_mutex);
            while(step_calc_done>k)
                pthread_cond_wait(&next_step_condition, &step_mutex);
        pthread_mutex_unlock(&step_mutex);
        if(step_calc_done<0)
            return 2;
#endif
//         TRACE(("k: %d\n", k));
//         TRACE(("Qu=\n"));
        // Qu  = cu(:,i)      + fu(:,:,i)'*Vx(:,i+1);
        memcpy(Qu, t->cu, sizeof(double)*N_U);
        addMulVec(Qu, Vx, t->fu, N_X, N_U);

//         TRACE(("Qx=\n"));
        // Qx  = cx(:,i)      + fx(:,:,i)'*Vx(:,i+1);
        memcpy(Qx, t->cx, sizeof(double)*N_X);
        addMulVec(Qx, Vx, t->fx, N_X, N_X);

//         TRACE(("Qxu=\n"));
        // Qux = cxu(:,:,i)'  + fu(:,:,i)'*Vxx(:,:,i+1)*fx(:,:,i);
        memcpy(Qxu, t->cxu, sizeof(double)*sizeofQxu);
        addMul2Tri(Qxu, Vxx, t->fx, N_X, N_X, t->fu, N_X, N_U, dummy);
        // fxuVx = vectens(Vx(:,i+1),fxu(:,:,:,i));
        // Qux   = Qux + fxuVx;
#if FULL_DDP
        for(j_= 0; j_<N_X*N_U; j_++) { // x, u
            d1= 0.0;
            for(i_= 0, k_= 0; i_<N_X; i_++, k_+= N_X*N_U) // f
                d1+= Vx[i_]*t->fxu[j_+k_];
            Qxu[j_]+= d1;
        }
#endif
//         TRACE(("Quu=\n"));
        // Quu = cuu(:,:,i)   + fu(:,:,i)'*Vxx(:,:,i+1)*fu(:,:,i);
        memcpy(Quu, t->cuu, sizeof(double)*sizeofQuu);
        addSquareTri(Quu, Vxx, t->fu, N_X, N_U, dummy);
        // fuuVx = vectens(Vx(:,i+1),fuu(:,:,:,i));
        // Quu   = Quu + fuuVx;
#if FULL_DDP
        for(j_= 0; j_<sizeofQuu; j_++) { // u, u
            d1= 0.0;
            for(i_= 0, k_= 0; i_<N_X; i_++, k_+= sizeofQuu) // f
                d1+= Vx[i_]*t->fuu[j_+k_];
            Quu[j_]+= d1;
        }                
#endif
        
//         TRACE(("Qxx=\n"));
        // Qxx = cxx(:,:,i)   + fx(:,:,i)'*Vxx(:,:,i+1)*fx(:,:,i);
        memcpy(Qxx, t->cxx, sizeof(double)*sizeofQxx);
        addSquareTri(Qxx, Vxx, t->fx, N_X, N_X, dummy);

        // Qxx = Qxx + vectens(Vx(:,i+1),fxx(:,:,:,i));
#if FULL_DDP
        for(j_= 0; j_<sizeofQxx; j_++) {// x, x
            d1= 0.0;
            for(i_= 0, k_= 0; i_<N_X; i_++, k_+= sizeofQxx) // f
                d1+= Vx[i_]*t->fxx[j_+k_];
            Qxx[j_]+= d1;
        }
#endif
        
//         TRACE(("regularization\n"));
        memcpy(QuuF, Quu, sizeof(double)*sizeofQuu);
        memcpy(Qxu_reg, Qxu, sizeof(double)*sizeofQxu);
        if(o->regType==2) {
//             TRACE(("type 2\n"));
            for(j_= 0; j_<N_U; j_++)
                for(i_= 0; i_<=j_; i_++) {
                    d1= 0.0;
                    for(k_= 0; k_<N_U; k_++)
                        d1+= t->fu[SYMTRI_MAT_IDX(k_, i_)]*t->fu[SYMTRI_MAT_IDX(k_, j_)];

                    QuuF[UTRI_MAT_IDX(i_, j_)]+= d1*o->lambda;
                }

            for(i_= 0; i_<N_X; i_++)
                for(j_= 0; j_<N_U; j_++) {
                    d1= 0.0;
                    for(k_= 0; k_<N_X; k_++)
                        d1+= t->fx[MAT_IDX(k_, i_, N_X)]*t->fu[MAT_IDX(k_, j_, N_U)];

                    Qxu_reg[MAT_IDX(i_, j_, N_X)]+= d1*o->lambda;
                }
        }
        if(o->regType==1) {
//             TRACE(("type 1\n"));
            for(i_= 0; i_<N_U; i_++) QuuF[UTRI_MAT_IDX(i_, i_)]+= o->lambda;          
        }

        // solve Quadratic Program
//         TRACE(("boxQP\n"));
        if(k==o->n_hor-1)
            memset(t->l, 0, sizeof(double)*N_U);
        else
            memcpy(t->l, (t+1)->l, sizeof(double)*N_U);

        if((qpRes= boxQP(QuuF, Qu, t->lower, t->upper, t->l, R, L, grad, grad_clamped, search, is_clamped, &m_free, invHfree, N_U))<1) {
            TRACE(("@k= %d: qpRes= %d \n", k, qpRes));
            return 1;
        }
        
        // Lfree        = -R\(R'\Qux_reg(free,:));
        // K_i(free,:)   = Lfree;
        memset(t->L, 0, sizeof(double)*N_U*N_X);
        for(i_= 0, i_free= 0; i_<N_U; i_++) {
            if(!is_clamped[i_]) {
                for(j_= 0, j_free= 0; j_<N_U; j_++) {
                    if(!is_clamped[j_]) {
                        for(l_= 0; l_<N_X; l_++)
                            t->L[MAT_IDX(i_, l_, N_U)]-= invHfree[SYMTRI_MAT_IDX(i_free, j_free)] * Qxu_reg[MAT_IDX(l_, j_, N_X)];
                        
                        j_free++;
                    } else {
                        d1= 0.0;
                        for(k_= 0, k_free= 0; k_<N_U; k_++) {
                            if(!is_clamped[k_]) {
                                d1-= invHfree[SYMTRI_MAT_IDX(i_free, k_free)] * QuuF[SYMTRI_MAT_IDX(k_, j_)];
                                k_free++;
                            }
                        }
                        for(l_= 0; l_<N_X; l_++)
                            t->L[MAT_IDX(i_, l_, N_U)]-= d1*((is_clamped[j_]==1)? t->lower_sign[j_]*t->lower_hx[MAT_IDX(l_, j_, N_X)]: t->upper_sign[j_]*t->upper_hx[MAT_IDX(l_, j_, N_X)]);
                    }
                }
                i_free++;
            } else {
                for(l_= 0; l_<N_X; l_++)
                    t->L[MAT_IDX(i_, l_, N_U)]-= (is_clamped[i_]==1)? t->lower_sign[i_]*t->lower_hx[MAT_IDX(l_, i_, N_X)]: t->upper_sign[i_]*t->upper_hx[MAT_IDX(l_, i_, N_X)];
            }
        }
        
        
        // dV          = dV + [k_i'*Qu  .5*k_i'*Quu*k_i];
        for(i_= 0; i_<N_U; i_++) {
            o->dV[0]+= Qu[i_]*t->l[i_];
        }
        for(i_= 0; i_<N_U; i_++) {
            d1= 0.0;
            for(j_= 0; j_<N_U; j_++) {
                d1+= t->l[j_]*Quu[SYMTRI_MAT_IDX(j_, i_)];
            }
            o->dV[1]+= 0.5*t->l[i_]*d1;
        }
//         TRACE(("dV= %g, %g, %g, %g, %g\n", o->dV[0], o->dV[1], t->l[0], Qu[0], Qu[1]));

//         TRACE(("Vx=\n"));
        // Vx(:,i)     = Qx  + K_i'*Quu*k_i + K_i'*Qu  + Qux'*k_i;
        memcpy(Vx, Qx, sizeof(double)*N_X);
        addMul2Tri(Vx, Quu, t->L, N_U, N_X, t->l, N_U, 1, dummy);
        for(i_= 0; i_<N_X; i_++) // row Vx
            for(j_= 0; j_<N_U; j_++) // col o->L' & row Qu
                Vx[i_]+= t->L[MAT_IDX(j_, i_, N_U)]*Qu[j_];

        for(i_= 0; i_<N_X; i_++) // row Vx
            for(j_= 0; j_<N_U; j_++) // col Qxu & row o->l
                Vx[i_]+= Qxu[MAT_IDX(i_, j_, N_X)]*t->l[j_];
        
    
//         TRACE(("Vxx=\n"));
        // Vxx(:,:,i)  = Qxx + K_i'*Quu*K_i + K_i'*Qux + Qux'*K_i;
        memcpy(Vxx, Qxx, sizeof(double)*sizeofQxx);
        addSquareTri(Vxx, Quu, t->L, N_U, N_X, dummy);
        for(i_= 0; i_<N_X; i_++)
            for(j_= 0; j_<N_X; j_++)
                for(k_= 0; k_<N_U; k_++) {
                    d1= t->L[MAT_IDX(k_, i_, N_U)]*Qxu[MAT_IDX(j_, k_, N_X)];
                    if(i_==j_) d1*= 2.0;
                    
                    Vxx[SYMTRI_MAT_IDX(i_, j_)]+= d1;
                }
        

//         TRACE(("g_norm=\n"));
        // g_norm= mean(max(abs(l) ./ (abs(u)+1),[],1));
        g_norm_max= 0.0;
        for(i_= 0; i_<N_U; i_++) {
            g_norm_i= fabs(t->l[i_]) / (fabs(t->u[i_])+1.0);
            if(g_norm_i>g_norm_max) g_norm_max= g_norm_i;
        }
        g_norm_sum+= g_norm_max;
    }
    
    o->g_norm= g_norm_sum/((double)(o->n_hor-1));
    
    return 0;
}
