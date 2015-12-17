//
//  LISAPathfinder.c
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

#include "LISAPathfinder.h"


void MomentOfInertia(double **I)
{
  double b = 1.0;

  I[0][0] = b*b/6.;
  I[0][1] = 0.0;
  I[0][2] = 0.0;

  I[1][0] = 0.0;
  I[1][1] = b*b/6.;
  I[1][2] = 0.0;

  I[2][0] = 0.0;
  I[2][1] = 0.0;
  I[2][2] = b*b/6.;

}

void draw_face(double *x, gsl_rng *seed)
{
  double b=1.0;
  double xmin = 0;
  double xmax = 4.*b;
  double ymin = 0;
  double ymax = 3.*b;

  //pick where on the surface
  x[0] = xmin + (xmax-xmin)*gsl_rng_uniform(seed);
  x[1] = ymin + (ymax-ymin)*gsl_rng_uniform(seed);
}

int which_face(double x, double y)
{
  double b=1.0;
  double xmin = 0;
  double xmax = 4.*b;
  double ymin = 0;
  double ymax = 3.*b;

  int i;

  int face = -1;

  for(i=0; i<6; i++)
  {
    switch(i)
    {
        //sides
      case 0:

        xmin=0.0;
        xmax=b;
        ymin=b;
        ymax=2*b;

        break;
      case 1:

        xmin=b;
        xmax=2*b;
        ymin=b;
        ymax=2*b;

        break;
      case 2:

        xmin=2*b;
        xmax=3*b;
        ymin=b;
        ymax=2*b;

        break;
      case 3:

        xmin=3*b;
        xmax=4*b;
        ymin=b;
        ymax=2*b;

        break;

        //bottom
      case 4:

        xmin=2*b;
        xmax=3*b;
        ymin=0;
        ymax=b;

        break;

        //top
      case 5:

        xmin=2*b;
        xmax=3*b;
        ymin=2*b;
        ymax=3*b;

        break;

      default:
        continue;
    }

    if(x>xmin && x<=xmax && y>ymin && y<=ymax) face=i;
  }

  return face;
}


void rotate_face_to_3D(double *x0, double *xp)
{
  int i;
  double b=1.0;

  double theta;

  double x = x0[0];
  double y = x0[1];
  double z = 0;

  double *k  = malloc(3*sizeof(double));
  double *kp = malloc(3*sizeof(double));
  double *r  = malloc(3*sizeof(double));

  int face = which_face(x0[0],x0[1]);

  //get point into cartesian coordinates
  k[0] = x;
  k[1] = y;
  k[2] = z;

  //shift down so face one is at the origin
  k[1] -= b;


  //rotate about x 90 degrees transforming y->z
  r[0] = 1;
  r[1] = 0;
  r[2] = 0;
  theta = M_PI/2.0;

  rotation(r,k,theta,kp);
  for(i=0; i<3; i++) k[i]=kp[i];

  if(face>0)
  {
    //translate face 1 to origin
    k[0] -= b;

    //rotate about x=1 z 90 degrees
    r[0] = 0;
    r[1] = 0;
    r[2] = 1;
    theta = M_PI/2.0;

    rotation(r,k,theta,kp);

    //shift back to origin
    kp[0] += b;

    for(i=0; i<3; i++) k[i]=kp[i];
  }

  if(face>1)
  {
    //translate face 2 to origin
    k[0] -= b;
    k[1] -= b;

    //rotate about x=1 z 90 degrees
    r[0] = 0;
    r[1] = 0;
    r[2] = 1;
    theta = M_PI/2.0;

    rotation(r,k,theta,kp);

    //shift back to origin
    kp[0] += b;
    kp[1] += b;

    for(i=0; i<3; i++) k[i]=kp[i];
  }

  if(face==3)
  {
    //translate face 3 to origin
    k[1] -= b;

    //rotate about x=1 z 90 degrees
    r[0] = 0;
    r[1] = 0;
    r[2] = 1;
    theta = M_PI/2.0;

    rotation(r,k,theta,kp);

    //shift back to origin
    kp[1] += b;

    for(i=0; i<3; i++) k[i]=kp[i];
  }

  //fold up bottom
  if(face==4)
  {
    //translate face 3 to origin
    k[1] -= b;

    //rotate up by x -90 degrees
    r[0] = 1;
    r[1] = 0;
    r[2] = 0;
    theta = -M_PI/2.0;

    rotation(r,k,theta,kp);

    //shift back to origin
    kp[1] += b;

    for(i=0; i<3; i++) k[i]=kp[i];
  }

  //fold down top
  if(face==5)
  {
    //translate face 3 to origin
    k[1] -= b;
    k[2] -= b;

    //rotate down by x 90 degrees
    r[0] = 1;
    r[1] = 0;
    r[2] = 0;
    theta = M_PI/2.0;

    rotation(r,k,theta,kp);

    //shift back to origin
    kp[1] += b;
    kp[2] += b;

    for(i=0; i<3; i++) k[i]=kp[i];
  }

  xp[0] = kp[0];
  xp[1] = kp[1];
  xp[2] = kp[2];

  free(k);
  free(kp);
  free(r);

}

void rotate_sky_to_3D(double *omega, int face, double *costheta, double *phi)
{
  int i;

  double theta;

  double x = omega[0];
  double y = omega[1];
  double z = omega[2];

  double *k  = malloc(3*sizeof(double));
  double *kp = malloc(3*sizeof(double));
  double *r  = malloc(3*sizeof(double));

  //get point into cartesian coordinates
  k[0] = x;
  k[1] = y;
  k[2] = z;


  //rotate about x 90 degrees transforming y->z
  r[0] = 1;
  r[1] = 0;
  r[2] = 0;
  theta = M_PI/2.0;

  rotation(r,k,theta,kp);
  for(i=0; i<3; i++) k[i]=kp[i];

  if(face>0)
  {

    //rotate about x=1 z 90 degrees
    r[0] = 0;
    r[1] = 0;
    r[2] = 1;
    theta = M_PI/2.0;

    rotation(r,k,theta,kp);

    for(i=0; i<3; i++) k[i]=kp[i];
  }

  if(face>1)
  {

    //rotate about x=1 z 90 degrees
    r[0] = 0;
    r[1] = 0;
    r[2] = 1;
    theta = M_PI/2.0;

    rotation(r,k,theta,kp);

    for(i=0; i<3; i++) k[i]=kp[i];
  }

  if(face==3)
  {

    //rotate about x=1 z 90 degrees
    r[0] = 0;
    r[1] = 0;
    r[2] = 1;
    theta = M_PI/2.0;

    rotation(r,k,theta,kp);

    for(i=0; i<3; i++) k[i]=kp[i];
  }

  //fold up bottom
  if(face==4)
  {

    //rotate up by x 90 degrees
    r[0] = 1;
    r[1] = 0;
    r[2] = 0;
    theta = -M_PI/2.0;

    rotation(r,k,theta,kp);

    for(i=0; i<3; i++) k[i]=kp[i];

  }

  //fold up bottom
  if(face==5)
  {

    //rotate up by x 90 degrees
    r[0] = 1;
    r[1] = 0;
    r[2] = 0;
    theta = M_PI/2.0;

    rotation(r,k,theta,kp);

    for(i=0; i<3; i++) k[i]=kp[i];
  }

  double phitemp = atan2(kp[1],kp[0]);
  while(phitemp<0.0)phitemp+=2.0*M_PI;

  *costheta = kp[2];
  *phi = phitemp;

  free(k);
  free(kp);
  free(r);


}


void rotation(double *r, double *k, double theta, double *kp)
{
  int i;

  /*
   rotation axis
   */
  double normalize;
  double n[3];

  normalize=0.0;
  for(i=0; i<3; i++)
  {
    n[i] = r[i];
    normalize += n[i]*n[i];
  }
  normalize = 1./sqrt(normalize);
  for(i=0; i<3; i++) n[i] *= normalize;

  /*
   rotation angle
   */
  double cosomega = cos(theta);
  double sinomega = sin(theta);
  double c1momega = 1.0 - cosomega;


  /*
   rotate k' = Rk
   */
  kp[0] = (c1momega*n[0]*n[0] + cosomega)     *k[0] + (c1momega*n[0]*n[1] - sinomega*n[2])*k[1] + (c1momega*n[0]*n[2] + sinomega*n[1])*k[2];
  kp[1] = (c1momega*n[0]*n[1] + sinomega*n[2])*k[0] + (c1momega*n[1]*n[1] + cosomega)     *k[1] + (c1momega*n[1]*n[2] - sinomega*n[0])*k[2];
  kp[2] = (c1momega*n[0]*n[2] - sinomega*n[1])*k[0] + (c1momega*n[1]*n[2] + sinomega*n[0])*k[1] + (c1momega*n[2]*n[2] + cosomega)     *k[2];
  
  
}


int check_impact(double costheta, double phi, int face)
{
  double n[3];
  double k[3];

  if(face==-1)return 1;

  //get norm to the current face
  switch(face)
  {
    case 0:

      n[0] =  0.;
      n[1] = -1.;
      n[2] =  0.;

      break;

    case 1:

      n[0] =  1.;
      n[1] =  0.;
      n[2] =  0.;
      
      break;

    case 2:

      n[0] =  0.;
      n[1] =  1.;
      n[2] =  0.;

      break;

    case 3:

      n[0] = -1.;
      n[1] =  0.;
      n[2] =  0.;

      break;


    case 4:

      n[0] =  0.;
      n[1] =  0.;
      n[2] = -1.;
      break;


    case 5:

      n[0] =  0.;
      n[1] =  0.;
      n[2] =  1.;
      break;


    default:
      return 1;
      break;
  }

  //compute line-of-sight vector k
  //get angle w.r.t. normal
  k[0] = sin(acos(costheta))*cos(phi);
  k[1] = sin(acos(costheta))*sin(phi);
  k[2] = costheta;

  //compute dot product between n and k
  double nk = n[0]*k[0] + n[1]*k[1] + n[2]*k[2];

  if(nk > 0) return 0;
  else       return 1;
  
}












