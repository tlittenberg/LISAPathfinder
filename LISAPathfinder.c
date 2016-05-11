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

#include "LPF.h"
#include "BayesLine.h"
#include "Subroutines.h"
#include "LISAPathfinder.h"


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


void draw_r(double *r, gsl_rng *seed)
{
  //draw direction from COM
  double costheta = -1.0 + 2.0*gsl_rng_uniform(seed);
  double phi      =  2.0*M_PI*gsl_rng_uniform(seed);
  
  
  //get unit vector
  double *rhat=malloc(3*sizeof(double));
  rhat[0] = sin(acos(costheta))*cos(phi);
  rhat[1] = sin(acos(costheta))*sin(phi);
  rhat[2] = costheta;
  
  //get face
  int face = which_face_r(rhat);
  
  //now with rhat and face calculate |r|
  double n[3]    = {0,0,0};
  double nhat[3] = {0,0,0};
  
  switch(face)
  {
    case 0:
      n[0] =  1.;
      nhat[0] =  0.5;
      break;
    case 1:
      n[1] =  1.;
      nhat[1] = 0.5;
      break;
    case 2:
      n[0] = -1.;
      nhat[0] = -0.5;
      break;
    case 3:
      n[1] = -1.;
      nhat[1] = -0.5;
      break;
    case 4:
      n[2] =  1.;
      nhat[2] = 0.5;
      break;
    case 5:
      n[2] = -1.;
      nhat[2] = -0.5;
      break;
    default:
      break;
  }
  
  int i;
  double ndotx = 0.5;
  double ndotr = 0.0;
  
  for(i=0; i<3; i++)
  {
    ndotr += n[i]*rhat[i];
    ndotx += n[i]*nhat[i];
  }
  
  double rmag = ndotx/ndotr;
  
  for(i=0; i<3; i++) r[i] = rmag*rhat[i];
  
  free(rhat);
  
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
  double *n=malloc(3*sizeof(double));
  double k[3];
  
  if(face==-1)
  {
    free(n);
    return 1;
  }

  //get norm to the current face
  get_normal(n,face);
  
  //compute line-of-sight vector k
  //get angle w.r.t. normal
  k[0] = sin(acos(costheta))*cos(phi);
  k[1] = sin(acos(costheta))*sin(phi);
  k[2] = costheta;
  
  //compute dot product between n and k
  double nk = n[0]*k[0] + n[1]*k[1] + n[2]*k[2];
  
  free(n);
  
  if(nk > 0) return 0;
  else       return 1;
  
}

int check_side(struct Spacecraft *lpf, double *r)
{
  //spacecraft dimensions
  //double A = lpf->x[2][1] - lpf->x[2][0];
  double D = lpf->x[0][0] - lpf->x[4][0];
  double H = lpf->x[2][1] - lpf->x[6][1];
  

  double *x0 = malloc(3*sizeof(double));
  double *xf = malloc(3*sizeof(double));
  int face = which_face_r(r);

  if(face<8)get_edge(lpf, x0, xf, face);
  else
  {
    x0[0] =  lpf->x[1][0];//0.921;
    x0[1] =  lpf->x[1][1];//
    
    xf[0] =  lpf->x[2][0];//0.521;
    xf[1] =  lpf->x[2][1];//0.903;
  }

  double m = -fabs((xf[1] - x0[1])/(xf[0] - x0[0]));

  int flag = 0;
  
  switch(face)
  {
    case 0:
      if(r[1]>xf[1] || r[1]<x0[1] || r[2]<0.0 || r[2]>xf[2]) flag = 1;
      //printf("x0={%g,%g,%g}, xf={%g,%g,%g}, r={%g,%g,%g}, flag=%i\n",x0[0],x0[1],x0[2],xf[0],xf[1],xf[2],r[0],r[1],r[2],flag);
      break;
    case 1:
      //if(fabs(r[0])>D*0.5 || fabs(r[1])>H*0.5 || fabs(r[1])-m*fabs(r[0])>fabs(x0[1]) - m*fabs(x0[0])) flag = 1;
      //printf("x0={%g,%g,%g}, xf={%g,%g,%g}, r={%g,%g,%g}, flag=%i\n",x0[0],x0[1],x0[2],xf[0],xf[1],xf[2],r[0],r[1],r[2],flag);
      break;
    case 2:
      if(r[0]<xf[0] || r[0]>x0[0] || r[2]<0.0 || r[2]>xf[2]) flag = 1;
      break;
    case 3:
      break;
    case 4:
      if(r[1]<xf[1] || r[1]>x0[1] || r[2]<0.0 || r[2]>xf[2]) flag = 1;
      //if(flag==1) printf("failed face 4: %g %g -> %g,%g %g,%g\n",r[1],r[2],xf[1],x0[1],0.0,xf[2]);
      break;
    case 5:
      break;
    case 6:
      if(r[0]>xf[0] || r[0]<x0[0] || r[2]<0.0 || r[2]>xf[2]) flag = 1;
      break;
    case 7:
      break;
    case 8:
      if(fabs(r[0])>D*0.5 || fabs(r[1])>H*0.5 || fabs(r[1])-fabs(x0[1]) > m*(fabs(r[0])-fabs(x0[0])) ) flag = 1;
      break;
    case 9:
      if(fabs(r[0])>D*0.5 || fabs(r[1])>H*0.5 || fabs(r[1])-fabs(x0[1]) > m*(fabs(r[0])-fabs(x0[0])) ) flag = 1;
      break;
    default:
      break;
  }
  
  
  free(x0);
  free(xf);
  
  return flag;
}

int which_side(struct Spacecraft *lpf, double *r)
{
  double A = lpf->H;//lpf->x[2][1] - lpf->x[2][0];
  //double D = lpf->x[0][0] - lpf->x[4][0];
  //double H = lpf->x[2][1] - lpf->x[6][1];

  int i,face;
  double *n=malloc(3*sizeof(double));
  double nk;
  double nk_max = -1.0;
  
  double rhat[3];
  double rnorm=0;
  for(i=0; i<3; i++)
  {
    rnorm += r[i]*r[i];
  }
  rnorm=sqrt(rnorm);
  for(i=0; i<3; i++) rhat[i] = r[i]/sqrt(rnorm);
  
  if(r[2]==A) face = 9;
  else if(r[2]==0) face = 8;
  else
  {
    
    //get norm to the current face
    for(i=0; i<8; i++)
    {
      switch(i)
      get_normal(n,i);
      
      //compute dot product between n and k
      nk = n[0]*rhat[0] + n[1]*rhat[1] + n[2]*rhat[2];
      if(nk > nk_max)
      {
        nk_max = nk;
        face = i;
      }
    }
  }
  free(n);
  return face;
}



void draw_octagon(struct Spacecraft *lpf, double *r, gsl_rng *seed)
{

  double A = lpf->H;
  double D = lpf->W;//lpf->x[0][0] - lpf->x[4][0];
  double H = lpf->D;//lpf->x[2][1] - lpf->x[6][1];

//  lpf->W = lpf->x[0][0] - lpf->x[4][0];
//  lpf->D = lpf->x[2][1] - lpf->x[6][1];

  //spacecraft dimensions
  //double a = 1.0;
  //double D = a*(1.0 + sqrt(2.0));

  double x0[2],xf[2];
  x0[0] =  lpf->x[1][0];//0.921;
  x0[1] =  lpf->x[1][1];//

  xf[0] =  lpf->x[2][0];//0.521;
  xf[1] =  lpf->x[2][1];//0.903;

  //line connecting
  // y - y0 = m*(x-x0)
  // y = mx - mx0 + y0
  // y - mx < y0 - mx0
  double m = (xf[1] - x0[1])/(xf[0] - x0[0]);



  double q[3];
  
  int pass = 0;
  
  double x,y,z;
  
  while(!pass)
  {
    x = gsl_rng_uniform(seed)*D/2.0;
    y = gsl_rng_uniform(seed)*H/2.0;
    z = A;


    if( y-m*x < x0[1] - m*x0[0] ) pass = 1;
  }
  
  
  int i;
  for(i=0; i<2; i++)
  {
    q[i] = 1.0;
    if(gsl_rng_uniform(seed)<0.5)q[i] *= -1.0;
  }
  
  r[0]=q[0]*x;
  r[1]=q[1]*y;
  r[2]=z;

}

void draw_side(struct Spacecraft *lpf, double *r, int face, gsl_rng *seed)
{
  
  double *x0,*xf;
  x0 = malloc(3*sizeof(double));
  xf = malloc(3*sizeof(double));
  
  //int face = (int)floor(gsl_rng_uniform(seed)*8.0);
  get_edge(lpf, x0, xf, face);
  
  if(fabs(xf[0]-x0[0]) > fabs(xf[1]-x0[1]))
  {
    //find x
    r[0] = x0[0] + (xf[0] - x0[0])*gsl_rng_uniform(seed);
    //solve for y
    r[1] = x0[1] + ((xf[1] - x0[1])/(xf[0] - x0[0]))*(r[0]-x0[0]);
  }
  else
  {
    //find y
    r[1] = x0[1] + (xf[1] - x0[1])*gsl_rng_uniform(seed);
    //solve for x
    r[0] = x0[0] + ((xf[0] - x0[0])/(xf[1] - x0[1]))*(r[1]-x0[1]);
  }

//
//
//  //find x
//  r[0] = x0[0] + (xf[0] - x0[0])*gsl_rng_uniform(seed);
//
//  if(face==0 || face==2 || face==4 || face== 6)
//    r[1] = x0[1] + (xf[1] - x0[1])*gsl_rng_uniform(seed);
//
//  else
//    r[1] = ((xf[1] - x0[1])/(xf[0] - x0[0]))*(r[0] - x0[0]) + x0[1];

  //z
  r[2] = x0[2] + (xf[2] - x0[2])*gsl_rng_uniform(seed);

  free(x0);
  free(xf);
}

int which_face_r(double *r)
{
  int i,face;
  double *n=malloc(3*sizeof(double));
  double nk;
  double nk_max = -1.0;

  if(r[2]>=SC_H)    face = 9;
  else if (r[2]<=0) face = 8;
  else
  {
    
    //get norm to the current face
    for(i=0; i<8; i++)
    {
      get_normal(n,i);
      
      //compute dot product between n and k
      nk = n[0]*r[0] + n[1]*r[1] + n[2]*r[2];
      if(nk > nk_max)
      {
        nk_max = nk;
        face = i;
      }
      
    }
  }
  free(n);
  
  return face;
}


void write_octagon()
{
  double a = 1.0;
  double D = a*(1+sqrt(2.));
  a/=2.0;
  D/=2.0;
  FILE *out;
  //draw bottom
  out = fopen("bottom_temp.dat","w");
  fprintf(out,"%lg %lg %lg\n",D,a,0.0);
  fprintf(out,"%lg %lg %lg\n",a,D,0.0);
  fprintf(out,"%lg %lg %lg\n",-a,D,0.0);
  fprintf(out,"%lg %lg %lg\n",-D,a,0.0);
  fprintf(out,"%lg %lg %lg\n",-D,-a,0.0);
  fprintf(out,"%lg %lg %lg\n",-a,-D,0.0);
  fprintf(out,"%lg %lg %lg\n",a,-D,0.0);
  fprintf(out,"%lg %lg %lg\n",D,-a,0.0);
  fprintf(out,"%lg %lg %lg\n",D,a,0.0);
  fclose(out);
  
  //draw top
  out = fopen("top_temp.dat","w");
  fprintf(out,"%lg %lg %lg\n",D,a,2*a);
  fprintf(out,"%lg %lg %lg\n",a,D,2*a);
  fprintf(out,"%lg %lg %lg\n",-a,D,2*a);
  fprintf(out,"%lg %lg %lg\n",-D,a,2*a);
  fprintf(out,"%lg %lg %lg\n",-D,-a,2*a);
  fprintf(out,"%lg %lg %lg\n",-a,-D,2*a);
  fprintf(out,"%lg %lg %lg\n",a,-D,2*a);
  fprintf(out,"%lg %lg %lg\n",D,-a,2*a);
  fprintf(out,"%lg %lg %lg\n",D,a,2*a);
  fclose(out);
  
  //draw sides
  out = fopen("sides_temp.dat","w");
  //all the way around the bottom
  fprintf(out,"%lg %lg %lg\n",D,a,0.0);
  fprintf(out,"%lg %lg %lg\n",a,D,0.0);
  fprintf(out,"%lg %lg %lg\n",-a,D,0.0);
  fprintf(out,"%lg %lg %lg\n",-D,a,0.0);
  fprintf(out,"%lg %lg %lg\n",-D,-a,0.0);
  fprintf(out,"%lg %lg %lg\n",-a,-D,0.0);
  fprintf(out,"%lg %lg %lg\n",a,-D,0.0);
  fprintf(out,"%lg %lg %lg\n",D,-a,0.0);
  fprintf(out,"%lg %lg %lg\n",D,a,0.0);
  //step up
  fprintf(out,"%lg %lg %lg\n",D,a,2*a);
  //all the way around the top
  fprintf(out,"%lg %lg %lg\n",a,D,2*a);
  fprintf(out,"%lg %lg %lg\n",-a,D,2*a);
  fprintf(out,"%lg %lg %lg\n",-D,a,2*a);
  fprintf(out,"%lg %lg %lg\n",-D,-a,2*a);
  fprintf(out,"%lg %lg %lg\n",-a,-D,2*a);
  fprintf(out,"%lg %lg %lg\n",a,-D,2*a);
  fprintf(out,"%lg %lg %lg\n",D,-a,2*a);
  fprintf(out,"%lg %lg %lg\n",D,a,2*a);
  fprintf(out,"%lg %lg %lg\n",a,D,2*a);
  //step down
  fprintf(out,"%lg %lg %lg\n",a,D,0.0);
  //work your way around the sides
  fprintf(out,"%lg %lg %lg\n",-a,D,0.0);
  fprintf(out,"%lg %lg %lg\n",-a,D,2*a);
  fprintf(out,"%lg %lg %lg\n",-D,a,2*a);
  fprintf(out,"%lg %lg %lg\n",-D,a,0.0);
  fprintf(out,"%lg %lg %lg\n",-D,-a,0.0);
  fprintf(out,"%lg %lg %lg\n",-D,-a,2*a);
  fprintf(out,"%lg %lg %lg\n",-a,-D,2*a);
  fprintf(out,"%lg %lg %lg\n",-a,-D,0.0);
  fprintf(out,"%lg %lg %lg\n",a,-D,0.0);
  fprintf(out,"%lg %lg %lg\n",a,-D,2*a);
  fprintf(out,"%lg %lg %lg\n",D,-a,2*a);
  fprintf(out,"%lg %lg %lg\n",D,-a,0.0);
  fprintf(out,"%lg %lg %lg\n",D,a,0.0);
  
}



void get_normal(double *n, int face)
{
  /************************************************/
  /*                                              */
  /*    xyz componants for spacecraft normals     */
  /*                                              */
  /************************************************/
  
  double e1 = 0.865769;
  double e2 = 0.500444;
  
  switch(face)
  {
    case 0:
      
      n[0] =  1.;
      n[1] =  0.;
      n[2] =  0.;
      
      break;
      
    case 1:
      
      n[0] =  e1;
      n[1] =  e2;
      n[2] =  0.;
      
      break;
      
    case 2:
      
      n[0] =  0.;
      n[1] = 1.;
      n[2] =  0.;
      
      break;
      
    case 3:
      
      n[0] = -e1;
      n[1] = e2;
      n[2] =  0.;
      
      break;
      
      
    case 4:
      
      n[0] = -1.;
      n[1] =  0.;
      n[2] =  0.;
      break;
      
      
    case 5:
      
      n[0] = -e1;
      n[1] = -e2;
      n[2] =  0.;
      break;
      
    case 6:
      
      n[0] =  0.;
      n[1] =  -1.;
      n[2] =  0.;
      break;
      
    case 7:
      
      n[0] =  e1;
      n[1] =  -e2;
      n[2] =  0.;
      break;
      
    case 8:
      
      n[0] =  0.;
      n[1] =  0.;
      n[2] =  -1.;
      break;
      
    case 9:
      
      n[0] =  0.;
      n[1] =  0.;
      n[2] =  1.;
      break;
      
    default:
      break;
  }
  
}


void get_edge(struct Spacecraft *lpf, double *x0, double *xf, int face)
{
  /************************************************/
  /*                                              */
  /*    xyz coordinates for spacecraft corners    */
  /*                                              */
  /************************************************/
//  double d = 0.831;
//  double w = 1.842;
//  double h = 1.806;

  int i;
  if(face>-1 && face < 8)
  {
    for(i=0; i<2; i++)
    {
      x0[i] = lpf->x[face][i];
      xf[i] = lpf->x[face+1][i];
    }
    x0[2] = 0.0;
    xf[2] = lpf->H;
  }
  else
  {
    for(i=0; i<2; i++)
    {
      x0[i] = lpf->x[1][i];
      xf[i] = lpf->x[2][i];
    }
    x0[2] = 0.0;
    xf[2] = lpf->H;
  }
  /*
  switch(face)
  {
    case 0:
      x0[0] =  0.921;
      x0[1] = -0.211;
      x0[2] =  0;
      
      xf[0] =  0.921;
      xf[1] =  0.211;
      xf[2] =  0.831;
      break;
    case 1:
      x0[0] =  0.921;
      x0[1] =  0.211;
      x0[2] =  0;
      
      xf[0] =  0.521;
      xf[1] =  0.903;
      xf[2] =  0.831;
      break;
    case 2:
      x0[0] =  0.521;
      x0[1] =  0.903;
      x0[2] =  0;
      
      xf[0] = -0.521;
      xf[1] =  0.903;
      xf[2] =  0.831;
      break;
    case 3:
      x0[0] = -0.521;
      x0[1] =  0.903;
      x0[2] =  0;
      
      xf[0] = -0.921;
      xf[1] =  0.211;
      xf[2] =  0.831;
      break;
    case 4:
      x0[0] = -0.921;
      x0[1] =  0.211;
      x0[2] =  0;
      
      xf[0] = -0.921;
      xf[1] = -0.211;
      xf[2] =  0.831;
      break;
    case 5:
      x0[0] = -0.921;
      x0[1] = -0.211;
      x0[2] =  0;
      
      xf[0] = -0.521;
      xf[1] = -0.903;
      xf[2] =  0.831;
      break;
    case 6:
      x0[0] = -0.521;
      x0[1] = -0.903;
      x0[2] =  0;
      
      xf[0] =  0.521;
      xf[1] = -0.903;
      xf[2] =  0.831;
      break;
    case 7:
      x0[0] =  0.521;
      x0[1] = -0.903;
      x0[2] =  0;
      
      xf[0] =  0.906;
      xf[1] = -0.217;
      xf[2] =  0.831;
      break;
    default:

      //face 1 is the default (top right quadrant of decks
      x0[0] =  0.921;
      x0[1] =  0.211;
      x0[2] =  0;

      xf[0] =  0.521;
      xf[1] =  0.903;
      xf[2] =  0.831;

      break;
  }
   */
  
}

void face2map(struct Spacecraft *lpf, double *r, double *x)
{
  //which face of spacecraft was hit?
  int face = which_face_r(r);
  
  //get normal vector for that face
  double *n = malloc(3*sizeof(double));
  get_normal(n,face);
  

  //get height on face
  double z = r[2];
  
  int i;

  double A = lpf->H;//lpf->x[2][1] - lpf->x[2][0];
  double D = lpf->x[0][0] - lpf->x[4][0];
  //double H = lpf->x[2][1] - lpf->x[6][1];


  
  //    if(face==6)
  //    {
  //        if(face==6||face==7)printf("%i r=(%g,%g,%g), x=(%g,%g), n=(%g,%g)\n",face,r[0],r[1],r[2],x[0],x[1],n[0],n[1]);
  //    }
  
  if(face==8) //bottom
  {
    for(i=0; i<2; i++) x[i] = r[i];
  }
  else if(face==9) //top
  {
    
    x[0] = A + D - r[0];
    x[1] = r[1];
    
  }
  else if(face<8)
  {
    for(i=0; i<2; i++)x[i] = r[i] + z*n[i];
  }
  
  
  //    if(face==8 || face==9)
  //    {
  //        printf("face=%i (%g,%g,%g):-->",face,r[0],r[1],r[2]);
  //        map2face(r, x);
  //        face = which_face_r(r);
  //        printf("%i (%g,%g,%g)\n",face,r[0],r[1],r[2]);
  //    }
  
  
}

void map2face(struct Spacecraft *lpf, double *r, double *x)
{
  int i;
  
  //spacecraft dimensions
  //double a = 1.0;
  //double D = a*(1.0 + sqrt(2.0));

  double x0[2],xf[2];
  x0[0] =  lpf->x[1][0];//0.921;
  x0[1] =  lpf->x[1][1];//
  
  xf[0] =  lpf->x[2][0];//0.521;
  xf[1] =  lpf->x[2][1];//0.903;

  //line connecting
  // y - y0 = m*(x-x0)
  // y = mx - mx0 + y0
  // y - mx < y0 - mx0
  double m = (xf[1] - x0[1])/(xf[0] - x0[0]);


  double A = lpf->H;//lpf->x[2][1] - lpf->x[2][0];
  double D = lpf->x[0][0] - lpf->x[4][0];
  double H = lpf->x[2][1] - lpf->x[6][1];

  //top face is way to the right, so that's easy
  if(x[0] > D/2.0 + A)
  {
    r[0] = A + D - x[0];
    r[1] = x[1];
    r[2] = A;
  }
  //bottom face is also a tad bit easier
  else if( fabs(x[0])<D*0.5 && fabs(x[1])<H*0.5 && fabs(x[1])<=m*(fabs(x[0])-x0[0]) + x0[1] )
  {
    r[0] = x[0];
    r[1] = x[1];
    r[2] = 0;
  }
  //the sides are really tricky
  else
  {
    //find intersection of x and edge
    double mx = x[1]/x[0]; //slope of map vector
    
    //find edge of octagon intersected by x
    double *x0 = malloc(3*sizeof(double));
    double *xf = malloc(3*sizeof(double));
    int face = which_face_r(r);
    get_edge(lpf, x0, xf, face);
    
    //find slope of edge
    double mp;
    if(face==0 || face==4) mp = 1.e32; //
    else                   mp = (xf[1] - x0[1])/(xf[0] - x0[0]);
    
    //find intersection point
    double *xi = malloc(2*sizeof(double));
    xi[0] = (xf[1] - mp*xf[0])/(mx - mp);
    xi[1] = mx*xi[0];
    
    //get norm to surface
    double *n = malloc(2*sizeof(double));
    get_normal(n, face);
    
    //find |z|
    double zmag = 0.0;
    for(i=0; i<2; i++) zmag += (x[i] - xi[i])*n[i];
    
    //find rx,ry
    for(i=0; i<2; i++) r[i] = x[i] - zmag*n[i];
    
    //add rz
    r[2] = zmag;
    
    free(x0);
    free(xf);
    free(xi);
    free(n);
  }
}



