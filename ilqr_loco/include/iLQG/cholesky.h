#ifndef CHOLESKY_H
#define CHOLESKY_H

int cholesky_tri(const double *A, int n, double *L);
void cholesky_tri_solve(const double *L_, const double *b, double *x, int n);
void cholesky_tri_inv(const double *L_, double *invA, const int n, double *x);

double mod_chol(double *A, int n, double *E, int *P, double *g);
void mod_chol_solve(const double *L_, const int *P, const double *b, double *x, int n, double *y);
void mod_chol_inv(const double *L_, const int *P, double *invA, int n, double *x);
void perm_tri_square(const double *L, double *H, const int *P, const int n);

#endif // CHOLESKY_H
