//
//  LISAPathfinder.h
//  
//
//  Created by Tyson Littenberg on 12/1/15.
//
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void MomentOfInertia(double **I);

void draw_face(double *x, gsl_rng *seed);

//void map_face_to_3D(double x, double y, int face, double *xp, double *yp, double *zp);

int which_face(double x, double y);

void rotate_face_to_3D(double *x0, double *xp);

void rotate_sky_to_3D(double *omega, int face, double *costheta, double *phi);

void rotation(double *r, double *k, double theta, double *kp);

int check_impact(double costheta, double phi, int face);

