#ifndef PRINTMAT_H
#define PRINTMAT_H

#include <string.h>
#include <stdio.h>

void printVec(const double *A, const int n, const char *nm);
void printTri(const double *A, const int n, const char *nm);
void printMat(const double *A, const int n, const int m, const char *nm);

#endif // PRINTMAT_H
