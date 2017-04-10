#ifndef BOXQP_H
#define BOXQP_H

int boxQP(double *H, const double *g, const double *lower, const double *upper, double *x, double *Hfree, double *L, double *grad, double *grad_clamped, double *search, int *is_clamped, int *n_free_, double *invHfree, const int n);

#endif /* BOXQP_H */
