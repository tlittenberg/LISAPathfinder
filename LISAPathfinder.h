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

void MomentOfInertia(double ***I);

void draw_face(double *x, gsl_rng *seed);
void draw_r(double *r, gsl_rng *seed);

//void map_face_to_3D(double x, double y, int face, double *xp, double *yp, double *zp);

int which_face(double x, double y);
int which_face_r(double *r);

void rotate_face_to_3D(double *x0, double *xp);

void rotate_sky_to_3D(double *omega, int face, double *costheta, double *phi);

void rotation(double *r, double *k, double theta, double *kp);

int check_impact_cube(double costheta, double phi, int face);
int check_impact(double costheta, double phi, int face);



void draw_octagon(struct Spacecraft *lpf, double *r, gsl_rng *seed);
void draw_side(struct Spacecraft *lpf, double *r, int face, gsl_rng *seed);
int check_side(struct Spacecraft *lpf, double *r);
int which_side(struct Spacecraft *lpf, double *r);

void write_octagon();


void get_normal(double *n, int face);
void get_edge(struct Spacecraft *lpf, double *x0, double *xf, int face);

void face2map(struct Spacecraft *lpf, double *r, double *x);
void map2face(struct Spacecraft *lpf, double *r, double *x);


