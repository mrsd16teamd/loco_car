#include <math.h>
#include "cholesky.h"
#include "printMat.h"
#include "matMult.h"

int cholesky_tri(const double *A, int n, double *L) {
    int i, j, k;
    double s;
    
    for(i = 0; i < n; i++)
        for(j = 0; j < (i+1); j++) {
            s= 0;
            for(k = 0; k < j; k++)
                s += L[UTRI_MAT_IDX(k, i)] * L[UTRI_MAT_IDX(k, j)];

            s= A[UTRI_MAT_IDX(j, i)] - s;
            if(i == j) {
                if(s<=0.0)
                    return 0;
                L[UTRI_MAT_IDX(j, i)]= sqrt(s);
            } else
//                L[i * n + j]= s / L[j * n + j];
                L[UTRI_MAT_IDX(j, i)]= 1.0 / L[UTRI_MAT_IDX(j, j)] * s;
        }
 
    return 1;
}

void cholesky_solve_tri(const double *L_, const double *b, double *x, int n) {
    int i, k;
    double sum;

    // Solve L*y= b;
    for(k= 0; k<n; k++) {
        sum= b[k];
        for(i= 0; i<k; i++)
            sum-= x[i]*L_[UTRI_MAT_IDX(i, k)];
        
        x[k]= sum/L_[UTRI_MAT_IDX(k, k)];
    }

    // Solve L'*X= Y;
    for(k= n-1; k >= 0; k--) {
       sum= x[k];
       for(i= k+1; i<n; i++)
             sum-= x[i]*L_[UTRI_MAT_IDX(k, i)];
       x[k]= sum/L_[UTRI_MAT_IDX(k, k)];
    }
}

void cholesky_tri_inv(const double *L_, double *invA, const int n, double *x) {
    int i, k, l;

    for(l= 0; l<n; l++) {
        x[l]= 1.0;
        for(k= l+1; k<n; k++) x[k]= 0.0;

        // Solve L*y= b;
        for(k= l; k<n; k++) {
            for(i= l; i<k; i++)
                x[k]-= x[i]*L_[UTRI_MAT_IDX(i, k)];
            x[k]/= L_[UTRI_MAT_IDX(k, k)];
        }

        // Solve L'*X= Y;
        for(k= n-1; k >= l; k--) {
            for(i= k+1; i<n; i++)
                 x[k]-= x[i]*L_[UTRI_MAT_IDX(k, i)];
            x[k]/= L_[UTRI_MAT_IDX(k, k)];
            
            invA[UTRI_MAT_IDX(l, k)]= x[k];
        }
    }
}

double fmax(double a, double b) {
    if(a>b)
        return a;
    else
        return b;
}

void switch_row_and_colMAT(double *A, int n, int i, int j) {
    double tmp;
    int k;

    for(k= 0; k<j; k++) {
        tmp= A[UTRI_MAT_IDX(k, j)];
        A[UTRI_MAT_IDX(k, j)]= A[UTRI_MAT_IDX(k, i)];
        A[UTRI_MAT_IDX(k, i)]= tmp;
    }
    
    for(k= j+1; k<i; k++) {
        tmp= A[UTRI_MAT_IDX(j, k)];
        A[UTRI_MAT_IDX(j, k)]= A[UTRI_MAT_IDX(k, i)];
        A[UTRI_MAT_IDX(k, i)]= tmp;
    } 

    for(k= i+1; k<n; k++) {
        tmp= A[UTRI_MAT_IDX(j, k)];
        A[UTRI_MAT_IDX(j, k)]= A[UTRI_MAT_IDX(i, k)];
        A[UTRI_MAT_IDX(i, k)]= tmp;
    }

    tmp= A[UTRI_MAT_IDX(j, j)];
    A[UTRI_MAT_IDX(j, j)]= A[UTRI_MAT_IDX(i, i)];
    A[UTRI_MAT_IDX(i, i)]= tmp;
}

void switch_row_and_colVEC(int *b, int i, int j) {
    int id;
    id= b[i];
    b[i]= b[j];
    b[j]= id;
}

void jthIteration(double *L_, int n, int j) {
    int i, k;

    L_[UTRI_MAT_IDX(j, j)]= sqrt(L_[UTRI_MAT_IDX(j, j)]);
    for(i=j+1; i<n; i++) {
        L_[UTRI_MAT_IDX(j, i)]/= L_[UTRI_MAT_IDX(j, j)];
        for(k=j+1; k<=i; k++) {
            L_[UTRI_MAT_IDX(k, i)]-= L_[UTRI_MAT_IDX(j, i)]*L_[UTRI_MAT_IDX(j, k)];
        }
    }
}

double mod_chol(double *A, int n, double *E, int *P, double *g) {
    int phaseone= 1;
    int i, j, k;
    double tau= pow(2.22044604925031e-16, 1./3.);
    double taubar= pow(2.22044604925031e-16, 2.0/3.0);
    double mu= 0.1;
    double gamma= 0.0;
    double tmp= 0.;
    double delta= 0.;
    double deltaprev= 0.;
    double normj;
    double lambda_hi;
    double lambda_lo;

    if(n==1) {
        delta= (taubar * fabs(A[0])) - A[0];
        if(delta>0.0) E[0]= delta; else E[0]= 0.;
        if(A[0]==0.0) E[0]= taubar;
        A[0]= sqrt(A[0]+E[0]);
        P[0]= 0;
        return E[0];
    }

    for(i= 0; i<n; i++) {
        P[i]= i;
        E[i]= 0.0;
    }
    for(i= 0; i<n; i++) {
        tmp= fabs(A[UTRI_MAT_IDX(i, i)]);
        if(tmp>gamma) gamma= tmp;
        if(A[UTRI_MAT_IDX(i, i)]<0.0) phaseone= 0;
    }

    /* Phase one, A potentially positive-definite */
    j= 0;
    while((j<n) && phaseone) {
        double tmp_max= A[UTRI_MAT_IDX(j, j)];
        double tmp_min= A[UTRI_MAT_IDX(j, j)];
        int id= j;

        // find extrema
        for(i= j+1; i<n; i++) {
            if(tmp_max<A[UTRI_MAT_IDX(i, i)]) {
                tmp_max= A[UTRI_MAT_IDX(i, i)];
                id= i;
            }
            if(tmp_min>A[UTRI_MAT_IDX(i, i)])
                tmp_min= A[UTRI_MAT_IDX(i, i)];
        }
//      if(tmp_max<taubar*gamma || tmp_min<tau*tmp_max) {
        if(tmp_max<taubar*gamma || tmp_min<-mu*tmp_max) {
            phaseone= 0;
            break; //go to phasetwo
        } else {
            /* Pivot on maximum diagonal of remaining submatrix */
            if(id!=j) {
                //switch rows and cols of id and j of A
                switch_row_and_colMAT(A, n, id, j);
                switch_row_and_colVEC(P, id, j);
            }

            tmp_min= 0.0;
            for(i=j+1; i<n; i++) {
                tmp= A[UTRI_MAT_IDX(i, i)] - A[UTRI_MAT_IDX(j, i)]*A[UTRI_MAT_IDX(j, i)]/A[UTRI_MAT_IDX(j, j)];
                if(tmp_min>tmp) tmp_min= tmp;
            }
//          if(tmp_min<tau*gamma) {
            if(tmp_min<-mu*gamma) {
                phaseone= 0;
                break;
            } else {// perform jth iteration of factorization
                jthIteration(A, n, j);
                j++;
            }
        }//end of if
    }//end of while

    /* Phase two, A not positive-definite */
    if(!phaseone && (j==(n-1))) {

        delta= -A[UTRI_MAT_IDX(n-1, n-1)] + fmax(tau*A[UTRI_MAT_IDX(n-1, n-1)]/(tau-1.), taubar*gamma);
        A[UTRI_MAT_IDX(n-1, n-1)]+= delta;
        A[UTRI_MAT_IDX(n-1, n-1)]= sqrt(A[UTRI_MAT_IDX(n-1, n-1)]);
        E[n-1]= delta;
        deltaprev= delta;
    }

    if(!phaseone && (j<(n-1))) {
        // k= number of iterations performed in phase one
        k= j - 1;

        /* Caculate lower Gerschgorin bound */
        for(i=k+1; i<n; i++) {
            g[i]= A[UTRI_MAT_IDX(i, i)];
            for(j=k+1; j<=i-1; j++) {
                g[i]-= fabs(A[UTRI_MAT_IDX(j, i)]);
            }
            for(j=i+1; j<n; j++) {
                g[i]-= fabs(A[UTRI_MAT_IDX(i, j)]);
            }
        }
        /* Modified Cholesky Decomposition */
        for(j=k+1; j<n-2; j++) {
            // Pivot on maximum lower Gerschgorin bound estimate
            int id= j;
            tmp= g[id];
            for(i=j+1; i<n; i++) {
                if(tmp<g[i]) {
                    tmp= g[i];
                    id= i;
                }
            }
            if(id!=j) {
                switch_row_and_colMAT(A, n, id, j);
                switch_row_and_colVEC(P, id, j);
                tmp= g[id];
                g[id]= g[j];
                g[j]= tmp;
            }
            //Calculate E[j, j] and add to diagonal
            normj= 0.;
            for(i=j+1; i<n; i++) {
                normj+= fabs(A[UTRI_MAT_IDX(j, i)]);
            }
            delta= fmax(0.0, fmax(fmax(normj, taubar*gamma)-A[UTRI_MAT_IDX(j, j)], deltaprev));
            if(delta>0) {
                A[UTRI_MAT_IDX(j, j)]+= delta;
                deltaprev= delta;
                E[j]= delta;
            }

            //Update Gerschgorin bound estimates
            if(A[UTRI_MAT_IDX(j, j)]!=normj) {
                tmp= 1.0 - normj/A[UTRI_MAT_IDX(j, j)];
                for(i=j+1; i<n; i++) {
                    g[i]+= fabs(A[UTRI_MAT_IDX(j, i)])*tmp;
                }
            }
            // Perform jth iteration of factorization
            jthIteration(A, n, j);
        }
        //Final 2x2 submatrix
        tmp= sqrt((A[UTRI_MAT_IDX(n-2, n-2)]-A[UTRI_MAT_IDX(n-1, n-1)])*(A[UTRI_MAT_IDX(n-2, n-2)]-A[UTRI_MAT_IDX(n-1, n-1)]) +4.0*A[UTRI_MAT_IDX(n-1, n-2)]*A[UTRI_MAT_IDX(n-1, n-2)]);
        lambda_hi= ((A[UTRI_MAT_IDX(n-2, n-2)]+A[UTRI_MAT_IDX(n-1, n-1)]) + tmp)*0.5;
        lambda_lo= ((A[UTRI_MAT_IDX(n-2, n-2)]+A[UTRI_MAT_IDX(n-1, n-1)]) - tmp)*0.5;
        delta= fmax(fmax(0.0, -lambda_lo + fmax(tau*(lambda_hi-lambda_lo)/(1.0-tau), taubar*gamma)), deltaprev);
        if(delta>0) {
            A[UTRI_MAT_IDX(n-2, n-2)]+= delta;
            A[UTRI_MAT_IDX(n-1, n-1)]+= delta;
            deltaprev= delta;
            E[n-2]= delta;
            E[n-1]= delta;
        }
        A[UTRI_MAT_IDX(n-2, n-2)]= sqrt(A[UTRI_MAT_IDX(n-2, n-2)]);
        A[UTRI_MAT_IDX(n-1, n-2)]/= A[UTRI_MAT_IDX(n-2, n-2)];
        A[UTRI_MAT_IDX(n-1, n-1)]= sqrt(A[UTRI_MAT_IDX(n-1, n-1)]-A[UTRI_MAT_IDX(n-1, n-2)]*A[UTRI_MAT_IDX(n-1, n-2)]);
    }
    return deltaprev;
}

void mod_chol_solve(const double *L_, const int *P, const double *b, double *x, int n, double *y) {
    int i, k;

    for(k= 0; k<n; k++) x[k]= b[P[k]];

    // Solve L*y= b;
    for(k= 0; k<n; k++) {
        for(i= 0; i<k; i++)
            x[k]-= x[i]*L_[UTRI_MAT_IDX(i, k)];
        x[k]/= L_[UTRI_MAT_IDX(k, k)];
    }

    // Solve L'*X= Y;
    for(k= n-1; k >= 0; k--) {
       for(i= k+1; i<n; i++)
             x[k]-= x[i]*L_[UTRI_MAT_IDX(k, i)];
       x[k]/= L_[UTRI_MAT_IDX(k, k)];
    }

    //Solve P'*y= x
    for(k= 0; k<n; k++) y[P[k]]= x[k];
    for(k= 0; k<n; k++) x[k]= y[k];
}

void mod_chol_inv(const double *L_, const int *P, double *invA, int n, double *x) {
    int i, k, l, cp, rp, tp;

    for(l= 0; l<n; l++) {
        x[l]= 1.0;
        for(k= l+1; k<n; k++) x[k]= 0.0;

        // Solve L*y= b; die Einträge in x sind nach P sortiert
        for(k= l; k<n; k++) {
            for(i= l; i<k; i++)
                x[k]-= x[i]*L_[UTRI_MAT_IDX(i, k)];
            x[k]/= L_[UTRI_MAT_IDX(k, k)];
        }

        // Solve L'*X= Y;
        for(k= n-1; k >= l; k--) {
            for(i= k+1; i<n; i++)
                 x[k]-= x[i]*L_[UTRI_MAT_IDX(k, i)];
            x[k]/= L_[UTRI_MAT_IDX(k, k)];

            rp= P[k]; cp= P[l]; if(rp>cp) { tp= rp; rp= cp; cp= tp; }
            invA[UTRI_MAT_IDX(rp, cp)]= x[k];
        }
    }
}

void perm_tri_square(const double *L, double *H, const int *P, const int n) {
    int r, c, i, r_, c_, i_;
    
    for(c= 0; c<n; c++) {
        for(r= 0; r<=c; r++) {
            c_= P[c];
            r_= P[r];
            if(r_>c_) {
                i_= r_;
                r_= c_;
                c_= i_;
            }
            H[UTRI_MAT_IDX(r_, c_)]= 0.0;
            for(i= 0; i<=r; i++)
                H[UTRI_MAT_IDX(r_, c_)]+= L[UTRI_MAT_IDX(i, c)]*L[UTRI_MAT_IDX(i, r)];
        }
    }
}
