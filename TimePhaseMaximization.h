
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double swap, tempr;
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void Sum_Extreme(double *a, double *b, double *Sn, int n, double *delt, double *pshift, double Tobs, int imin, int imax, double t0);
