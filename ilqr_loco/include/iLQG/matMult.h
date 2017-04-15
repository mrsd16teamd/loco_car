#ifndef MATMUL_H
#define MATMUL_H

#define MAT_IDX(r, c, nr) ((r) + (c)*(nr))
#define MAT_IDX3(r, c, b, nr, nc) ((r) + (c)*(nr) + (b)*(nr)*(nc))
#define MAT_IDX4(r, c, i3, i4, nr, nc, n3) ((r) + (c)*(nr) + (i3)*(nr)*(nc) + (i4)*(nr)*(nc)*(n3))

#define UTRI_MAT_IDX(r, c) (((c)*((c)+1))/2 + (r)) // row first order
#define SYMTRI_MAT_IDX(r, c) ((r>c)? UTRI_MAT_IDX(c, r): UTRI_MAT_IDX(r, c))

void addMulVec(double base[], const double a[], const double b[], const int n_r, const int n_c);
void addSquareTri(double base[], const double b[], const double a[], const int n_r, const int n_c, double ba[]);
void addMul2Mat(double base[], const double b[], const double a[], const int n_ra, const int n_ca, const double c[], const int n_rc, const int n_cc, double bc[]);
void addMul2Tri(double base[], const double b[], const double a[], const int n_ra, const int n_ca, const double c[], const int n_rc, const int n_cc, double bc[]);

#endif // MATMUL_H
