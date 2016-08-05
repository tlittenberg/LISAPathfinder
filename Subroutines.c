/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


#include "BayesLine.h"
#include "Subroutines.h"
#include "LISAPathfinder.h"
#include "TimePhaseMaximization.h"
#include "LPF.h"

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define eps 1.2e-7
#define RNMX (1.0-eps)
double swap, tempr;
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

/* ********************************************************************************** */
/*                                                                                    */
/*                                  Data structures                                   */
/*                                                                                    */
/* ********************************************************************************** */


/* ********************************************************************************** */
/*                                                                                    */
/*                                    MCMC tools                                      */
/*                                                                                    */
/* ********************************************************************************** */

//static void system_pause()
//{
//  printf("Press Any Key to Continue\n");
//  getchar();
//}

void bayesline_mcmc(struct Data *data, struct Model **model, struct BayesLineParams ***bayesline, int *index, double beta, int ic)
{
  int ifo,i;

  int N  = data->N;
  int NI = data->DOF;


  //pointers to structures
  struct Model *model_x = model[index[ic]];
  struct BayesLineParams **bl_x = bayesline[index[ic]];

  double *r = malloc(2*N*sizeof(double));

  //update PSD for each interferometer
  for(ifo=0; ifo<NI; ifo++)
  {
    //copy over current multi-template & current residual
    for(i=0; i<2*N; i++)
    {
      r[i] = data->d[ifo][i]-model_x->s[ifo][i];
    }

    //re-run Markovian, full spectrum, full model part of BayesLine
    BayesLineRJMCMC(bl_x[ifo], r, model_x->Snf[ifo], model_x->invSnf[ifo], model_x->SnS[ifo], N*2, 10, beta, 1);

    for(i=0; i<N; i++)
    {
      model_x->Snf[ifo][i]*=2.0;
      model_x->SnS[ifo][i] = model_x->Snf[ifo][i];
      model_x->invSnf[ifo][i] = 1./model_x->Snf[ifo][i];
    }
  }

  free(r);
}

void ptmcmc(struct Model **model, double *temp, int *index, gsl_rng *r, int NC, int mc)
{
  /*
   NC is passed to PTMCMC redundendtly because the main code
   can reduce the number of chains based on the temperature
   of the hottest chain
   */
  int a, b, ic;
  int olda, oldb;
  
  double heat1, heat2;
  double logL1, logL2;
  double dlogL;
  double H;
  double alpha;
  double beta;
  
  double A[NC],S[NC];
  
  //b = (int)(ran2(seed)*((double)(chain->NC-1)));
  for(b=NC-1; b>0; b--)
  {
    a = b - 1;
    
    olda = index[a];
    oldb = index[b];
    
    heat1 = temp[a];
    heat2 = temp[b];
    
    logL1 = model[olda]->logL;
    logL2 = model[oldb]->logL;
    
    //Hot chains jump more rarely
    if(gsl_rng_uniform(r)<1.0)///heat1)
    {
      dlogL = logL2 - logL1;
      H  = (heat2 - heat1)/(heat2*heat1);
      
      alpha = exp(dlogL*H);
      beta  = gsl_rng_uniform(r);
      
      if(alpha >= beta)
      {
        index[a] = oldb;
        index[b] = olda;
        A[a]=1;
      }
      else A[a]=0;
    }
  }
  
  double nu=30;
  double t0=100;
  
  for(ic=1; ic<NC-1; ic++)
  {
    S[ic] = log(temp[ic] - temp[ic-1]);
  }
  
  ic=0;
  for(ic=1; ic<NC-1; ic++)
  {
    S[ic] += (A[ic-1] - A[ic])*(t0/((double)mc+t0))/nu;
    
    temp[ic] = temp[ic-1] + exp(S[ic]);
    
  }//end loop over ic
  //  for(ic=0; ic<NC; ic++)printf("%.2g ",temp[ic]);
  //  printf("\n");
  
}


void draw_impact_point(struct Data *data, struct Spacecraft *lpf, struct Source *source, gsl_rng *seed)
{
//  double a = 1.0;
//  double A = 2.*(1.0 + sqrt(2.))*a*a;
  //TODO: Areas need to be computed from LPF.h
  double area[10];
  area[0] = 0.350682;
  area[1] = 0.66421;
  area[2] = 0.865902;
  area[3] = 0.66421;
  area[4] = 0.350682;
  area[5] = 0.66421;
  area[6] = 0.865902;
  area[7] = 0.66421;
  area[8] = 1.99823;
  area[9] = 1.99823;

  //pick a side, any side
  source->face = (int)floor(10*gsl_rng_uniform(seed));
  if(source->face == 9)
  {
    draw_octagon(lpf, source->r, seed);
    source->r[2] = source->r[2];
  }
  else if(source->face == 8)
  {
    draw_octagon(lpf, source->r, seed);
    source->r[2] = 0.0;
  }
  else if(gsl_rng_uniform(seed) < area[source->face]/area[9])
  {
    draw_side(lpf, source->r, source->face, seed);
  }
  else source->face=-1;

  //printf("trial side=%i ",source->face);
  //map side to x-y plane
  if(source->face > -1) face2map(lpf, source->r, source->map);
  
  
  //pick sky location
  source->costheta = -1.0 + 2.0*gsl_rng_uniform(seed);
  source->phi      = gsl_rng_uniform(seed)*2.0*M_PI;
//  source->costheta = -0.0763315;//-1.0 + 2.0*gsl_rng_uniform(seed);
//  source->phi      = 3.43195;//gsl_rng_uniform(seed)*2.0*M_PI;

  
  //flag if face is not exposed to sky location
  if(check_impact(source->costheta, source->phi, source->face)) source->face = -1;

  //momentum and impact time
  //printf("fudge\n");
//  source->P  = gsl_ran_exponential(seed,20);
//  source->t0 = gsl_rng_uniform(seed)*data->T;

}

void proposal(struct Flags *flags, struct Data *data, struct Spacecraft *lpf, struct Model *model, struct Model *trial, gsl_rng *r, int *reject, int nmax, int *drew_prior)
{
  int n;

    //model->logQ = 0.0;
    //trial->logQ = 0.0;

  //Always update spacecraft parameters
  detector_proposal(data, model, trial, r);
  
    //MCMC proposal
    if(gsl_rng_uniform(r)<0.99 || !flags->rj)
    {
        for(n=0; n<model->N; n++){
            //printf("(source=%i) ",n);
            //printf("in-proposal: n=%i\n",n);
            if(flags->use_spacecraft==0)
                impact_proposal(data, lpf, model->source[n], trial->source[n], r);
            else
                impact_proposal_sc(data, lpf, model->source[n], trial->source[n], r, &drew_prior[n]);
        }
    }
    //RJ proposal
    else dimension_proposal(flags, data, lpf, model, trial, r, nmax, reject);
  
}

/*  John's version of this with some fixes especially to make draw consistent with a meanaingful prior.
void draw_impact_point(struct Data *data, struct Spacecraft *lpf, struct Source *source, gsl_rng *seed)
{
//  double a = 1.0;
//  double A = 2.*(1.0 + sqrt(2.))*a*a;
  double area[10];
  area[0] = 0.350682;
  area[1] = 0.66421;
  area[2] = 0.865902;
  area[3] = 0.66421;
  area[4] = 0.350682;
  area[5] = 0.66421;
  area[6] = 0.865902;
  area[7] = 0.66421;
  area[8] = 2.77331;
  area[9] = 2.77331;
  
  source->face=-1;
  while(source->face<0){
    //pick a side, any side
    source->face = (int)floor(10*gsl_rng_uniform(seed));
    if(source->face == 9)
      {
	draw_octagon(lpf, source->r, seed);
	source->r[2] = source->r[2];
      }
    else if(source->face == 8)
      {
	draw_octagon(lpf, source->r, seed);
	source->r[2] = 0.0;
      }
    else if(gsl_rng_uniform(seed) < area[source->face]/area[9])
      {
	draw_side(lpf, source->r, source->face, seed);
      }
    else source->face=-1;
    
    //printf("trial side=%i ",source->face);
    //map side to x-y plane
    if(source->face > -1) face2map(lpf, source->r, source->map);
    
    //pick sky location
    //source->costheta = -1.0 + 2.0*gsl_rng_uniform(seed);
    //source->phi      = gsl_rng_uniform(seed)*2.0*M_PI;
    double cos2thetaF= gsl_rng_uniform(seed);
    double costhetaF=sqrt(cos2thetaF);
    double phiF=gsl_rng_uniform(seed)*2.0*M_PI;
    face_sky_to_body_sky(source->face,costhetaF,phiF,&source->costheta,&source->phi);
    while(source->phi > 2.0*M_PI) source->phi -= 2.0*M_PI;
    while(source->phi < 0)        source->phi += 2.0*M_PI;
    
    //flag if face is not exposed to sky location
    static int count=0;
    if(check_impact(source->costheta, source->phi, source->face)){//
      if(source->face>-1)printf("This should (almost?) never fail. (count=%i),costhetaF=%g\n",count,costhetaF);
      source->face = -1;
    }
    count++;
    //momentum and impact time
    //printf("fudge\n");
    source->P  = gsl_ran_exponential(seed,20);
    source->t0 = gsl_rng_uniform(seed)*data->T;
  }
}
*/

void draw_impact_point_sc(struct Data *data, struct Spacecraft *lpf, struct Source *source, gsl_rng *seed)
{
  enum{ nface=10};
  int i;
  static double cum_area[nface];
  static double total_area;
  static int have_cum_area=0;
  //compute and store the cumulative area
  if(have_cum_area==0){
    have_cum_area=1;
    double sum=0;
    for(i=0;i<nface;i++){
      double area=lpf->faces[i]->area;
      sum+=area;
      printf("area[%i]=%g, cum=%g\n",i,area,sum);
      cum_area[i]=sum;
    }
    total_area=sum;
  }

  //draw a face distrib by area
  double x = gsl_rng_uniform(seed);
  source->face=-1;
  for(i=0;i<nface;i++)
    if(cum_area[i]>=x*total_area){
      source->face=i;
      break;
    }
  int iface=-1;

  //draw an impact point
  double rface[2];
  while(iface!=source->face){//since faces aren't nec. rectangles, keep drawing until it sticks
    iface=source->face;
    double x=gsl_rng_uniform(seed)*lpf->faces[iface]->rmax[0];
    double y=gsl_rng_uniform(seed)*lpf->faces[iface]->rmax[1];
    //printf(" drew (%g,%g) on face %i\n",x,y,iface);
    rface[0]=x;
    rface[1]=y;  
    adjust_face(lpf,&iface,rface,seed);
  }
  face2body(lpf,iface,rface,source->r);
  
  //pick sky location
  double cos2thetaF= gsl_rng_uniform(seed);
  double costhetaF=sqrt(cos2thetaF);
  double phiF=gsl_rng_uniform(seed)*2.0*M_PI;
  face_sky_to_body_sky(lpf,source->face,costhetaF,phiF,&source->costheta,&source->phi);
  while(source->phi > 2.0*M_PI) source->phi -= 2.0*M_PI;
  while(source->phi < 0)        source->phi += 2.0*M_PI;
  
  //flag if face is not exposed to sky location
  static int count=0;
  if(incidence(lpf,source->face,source->costheta, source->phi)<0){
    //printf("face=%i, cth=%g, phi=%g, inc=%g\n",source->face,source->costheta,source->phi,incidence(lpf,source->costheta, source->phi, source->face));
    if(source->face>-1)printf("This should (almost?) never fail. (count=%i),costhetaF=%g\n",count,costhetaF);
    source->face = -1;
  }
  count++;
  //printf("drew impact point on face=%i\n",source->face);
}

void draw_impactor(struct Data *data, struct Source *source, gsl_rng *seed)
{
  //momentum and impact time
  source->P  = gsl_ran_exponential(seed,10000);
  source->t0 = gsl_rng_uniform(seed)*data->T;
}



void detector_proposal(struct Data *data, struct Model *model, struct Model *trial, gsl_rng *r)
{
  int i;
  //trial->mass = model->mass + gsl_ran_ugaussian(r)*1.0;     // kg
  for(i=0; i<3; i++)
  {
    trial->Ais[i]  = model->Ais[i];//  + gsl_ran_ugaussian(r)*2.0e-11; // m*Hz^-1/2
    trial->Ath[i]  = model->Ath[i];//  + gsl_ran_ugaussian(r)*1.0e-10; // N*Hz^-1/2
    trial->Ars[i]  = model->Ars[i];//  + gsl_ran_ugaussian(r)*1.0e-9; // rad*Hz^-1/2
  }
}

void impact_proposal(struct Data *data, struct Spacecraft *lpf, struct Source *model, struct Source *trial, gsl_rng *r)
{
  //uniform
  if(gsl_rng_uniform(r)<0.5)
  {
    draw_impact_point(data,lpf,trial,r);

    //10% of time also draw time & amplitude from prior
    if(gsl_rng_uniform(r)<0.0)
    {
      draw_impactor(data, trial, r);
    }

  }
  //gaussian
  else
  {
    trial->P  = model->P  + gsl_ran_ugaussian(r)*10.0;
    trial->t0 = model->t0 + gsl_ran_ugaussian(r)*0.25;
    
    trial->phi      = model->phi + gsl_ran_ugaussian(r)*0.01;
    trial->costheta = model->costheta + gsl_ran_ugaussian(r)*0.01;
    
    //map to 2D
    trial->r[0] = model->r[0];
    trial->r[1] = model->r[1];
    trial->r[2] = model->r[2];

    model->face = which_face_r(trial->r);

    trial->face = model->face;

    //printf("face=%i, r0={%g,%g,%g} ",trial->face, trial->r[0],trial->r[1], trial->r[2]);

    face2map(lpf, model->r, model->map);

//    if(model->face==4)printf("x={%g,%g} ",model->r[0],model->r[1]);

    //perturb map coordinates
    if(gsl_ran_ugaussian(r)<0.5)
    {
      trial->map[0] = model->map[0] + gsl_ran_ugaussian(r)*0.01;
      trial->map[1] = model->map[1] + gsl_ran_ugaussian(r)*0.01;
      //printf("gaussian small\n");
    }
    else
    {
      trial->map[0] = model->map[0] + gsl_ran_ugaussian(r)*0.05;
      trial->map[1] = model->map[1] + gsl_ran_ugaussian(r)*0.05;
      //printf("gaussian big\n");
    }
    
    //map back to 3D
    map2face(lpf, trial->r, trial->map);
    
    //TODO: gaussian proposal is not robust to changing faces
    trial->face = which_face_r(trial->r);

//    if(trial->face==4)printf("rf={%g,%g,%g} face=%i\n", trial->r[0],trial->r[1], trial->r[2],trial->face);

    int flag = 0;
    if(check_impact(trial->costheta, trial->phi, trial->face))
    {flag = 1;
//    if(model->face == 4){
//      printf("jump failed condition %i\n",flag);
//      printf("   trial face=%i, r0={%g,%g,%g} sky={%g,%g}\n",trial->face, trial->r[0],trial->r[1], trial->r[2],trial->costheta, trial->phi);
//    }
    }
    else if(trial->face!=model->face)
    {
      flag = 2;
//      if(model->face == 1){
//      printf("%i->%i\n",model->face,trial->face);
//      printf("jump failed condition %i\n",flag);
//      printf("   model face=%i, r0={%g,%g,%g}\n",model->face, model->r[0],model->r[1], model->r[2]);
//      printf("   trial face=%i, r0={%g,%g,%g}\n",trial->face, trial->r[0],trial->r[1], trial->r[2]);
//      }
    }
    else if(check_side(lpf, trial->r))
    {
      flag = 3;
      //if(trial->face==1)printf("jump failed condition %i\n",flag);
    }
    if(flag>0)
    {
      //if(trial->face==1)printf("jump failed condition %i\n",flag);
      trial->face = -1;
    }
    //get 3D position
    //    rotate_face_to_3D(trial->map, trial->r);
    
    while(trial->phi > 2.0*M_PI) trial->phi -= 2.0*M_PI;
    while(trial->phi < 0)        trial->phi += 2.0*M_PI;
    
    while(trial->costheta >  1.0) trial->costheta =  2.0-trial->costheta;
    while(trial->costheta < -1.0) trial->costheta = -2.0-trial->costheta;
    
  }
  
  //    if(trial->face > -1)
  //    {
  //  printf("face=%i (%g,%g,%g):-->",trial->face,trial->r[0],trial->r[1],trial->r[2]);
  //  face2map(trial->r, trial->map);
  //
  //        trial->face = which_face_r(trial->r);
  //
  //    map2face(trial->r, trial->map);
  //  printf("%i (%g,%g,%g)\n",trial->face,trial->r[0],trial->r[1],trial->r[2]);
  //    }
  
}
  
void impact_proposal_sc(struct Data *data, struct Spacecraft *lpf, struct Source *model, struct Source *trial, gsl_rng *r, int *drew_prior)
{
  
  //uniform
  if(gsl_rng_uniform(r)<0.5)
  {
    draw_impact_point_sc(data,lpf,trial,r);

    //10% of time also draw time & amplitude from prior
    if(gsl_rng_uniform(r)<0.0) //TODO: Fix within-model prior-draw proposal
    {
      draw_impactor(data, trial, r);
    }

    //printf("draw\n");
    *drew_prior=1;
  }
  //gaussian
  else
  {
    double *rface=malloc(2*sizeof(double));

    //printf("gauss\n");
    trial->P  = model->P  + gsl_ran_ugaussian(r)*10.0;
    trial->t0 = model->t0 + gsl_ran_ugaussian(r)*0.25;
    
    trial->phi      = model->phi + gsl_ran_ugaussian(r)*0.01;
    trial->costheta = model->costheta + gsl_ran_ugaussian(r)*0.01;
    
    trial->r[0] = model->r[0];
    trial->r[1] = model->r[1];
    trial->r[2] = model->r[2];

    trial->face = model->face;

    //printf("face=%i, r0={%g,%g,%g} ",trial->face, trial->r[0],trial->r[1], trial->r[2]);
    //debug check
    //body2face(lpf,model->face,trial->r,rface);      

    int iface=model->face;
    //printf("\n\n\nGaussian step:\n");
    body2face(lpf,iface,trial->r,rface);      
    //debug check
    //face2body(lpf,iface,rface,trial->r);
    
    double scale=0.05;
    if(gsl_ran_ugaussian(r)<0.5)scale=0.01;
    double stepx=gsl_ran_ugaussian(r)*scale;
    double stepy=gsl_ran_ugaussian(r)*scale;
    double deltascale=sqrt(stepx*stepx+stepy*stepy);
    rface[0] += stepx;
    rface[1] += stepy;

    adjust_face(lpf,&iface,rface,r);
    face2body(lpf,iface,rface,trial->r);
    trial->face=iface;
    
    //debug check
    if(0 && iface>=0){
      //sanity check on results for testing:
      double dx=model->r[0] - trial->r[0];
      double dy=model->r[1] - trial->r[1];
      double dz=model->r[2] - trial->r[2];
      double delta=sqrt(dx*dx+dy*dy+dz*dz)/deltascale;
      if(delta>1.0001||delta<0.707){	  
	printf("--Problem: Step too %s model->r=(%g,%g,%g)  trial->r=(%g,%g,%g)  deltascale=%g  delta=%g\n",(delta>1?"large":"small"),
	       model->r[0],model->r[1],model->r[2],trial->r[0],trial->r[1],trial->r[2],deltascale,delta);
	body2face(lpf,model->face,model->r,rface);
	printf("  Model face[%i],rface=(%g,%g)\n",model->face,rface[0],rface[1]);
	body2face(lpf,trial->face,trial->r,rface);
	printf("  Trial face[%i],rface=(%g,%g)\n",trial->face,rface[0],rface[1]);
      }
    }
    
    while(trial->phi > 2.0*M_PI) trial->phi -= 2.0*M_PI;
    while(trial->phi < 0)        trial->phi += 2.0*M_PI;
    
    while(trial->costheta >  1.0) trial->costheta =  2.0-trial->costheta;
    while(trial->costheta < -1.0) trial->costheta = -2.0-trial->costheta;

    *drew_prior=0;
    //printf("returning from proposal trial->face=%i\n",trial->face);
  }

}

void dimension_proposal(struct Flags *flags, struct Data *data, struct Spacecraft *lpf, struct Model *model, struct Model *trial, gsl_rng *r, int Nmax, int *test)
{
  int n,kill;
  
  //birth
  if(gsl_rng_uniform(r)<0.5)
  {
    if(model->N==Nmax) *test = 1;
    else
    {
      if(flags->use_spacecraft==0)draw_impact_point(data,lpf,trial->source[trial->N],r);
      else draw_impact_point_sc(data,lpf,trial->source[trial->N],r);
      draw_impactor(data, trial->source[trial->N], r);
      //trial->logQ += log(gsl_ran_exponential_pdf(trial->source[trial->N]->P,1000));
      trial->N = model->N+1;
    }
  }
  //death
  else
  {
    if(model->N==0) *test = 1;
    else
    {
      //choose which to kill
      kill = (int)floor(gsl_rng_uniform(r)*(double)model->N);
      //model->logQ += log(gsl_ran_exponential_pdf(model->source[kill]->P,1000));
      for(n=kill; n<model->N-1; n++) copy_source(model->source[n+1], trial->source[n], data->DOF);
      trial->N = model->N-1;
    }
  }
}

void logprior(struct Data *data, struct Model *model, struct Model *injection)
{
  //model->logP = log_mass_prior(injection->mass,model->mass);
  
  int n;
  for(n=0; n<model->N; n++)
  {
    
    if(model->source[n]->t0 < 0.0 || model->source[n]->t0 > data->T) model->logP = -1.0e60;
    
    if(model->source[n]->P < 0.0) model->logP = -1.0e60;
  }
}

void logprior_sc(struct Data *data, struct Spacecraft *lpf, struct Model *model, struct Model *injection, int *drew_prior)
{
    //model->logP = log_mass_prior(injection->mass,model->mass);
    ///JGB:Seems that the prior was not initialized.  How did this work?
    model->logP=0;

    int n;
    for(n=0; n<model->N; n++)
    {

        if(model->source[n]->t0 < 0.0 || model->source[n]->t0 > data->T) model->logP = -1.0e60;

        if(model->source[n]->P < 0.0) model->logP = -1.0e60;
        //else model->logP += log(gsl_ran_exponential_pdf(model->source[n]->P,1000));

        ///JGB:Prior on the face/impact-dir, but not if the trial was drawn from prior/
        if(drew_prior[n]==0){
            int iface=model->source[n]->face;
            if(iface>=0){
                model->logP+=log(lpf->faces[iface]->area/lpf->lpf_area);
                double inc_cos=incidence(lpf,iface,model->source[n]->costheta,model->source[n]->phi);
                if(inc_cos>0)
                    model->logP+=log(inc_cos);
                else
                    model->logP=-1e60;
            } else model->logP=-1e60;
            
        }
    }
}

//just for debugging
void check_incidence(struct Spacecraft *lpf,struct Model *model){
  int i;
  for(i=0;i<model->N;i++){
    int iface=model->source[i]->face;
    double cth_inc=-1;
    if(iface>=0){
      double cth=model->source[i]->costheta;
      double phi=model->source[i]->phi;
      cth_inc=incidence(lpf,iface,cth,phi);
    }
    if(cth_inc<0)printf(" invalid incidence (cth_inc=%g) for source[%i] on face %i\n",cth_inc,i,iface);
  }
}  


double log_mass_prior(double m0, double m)
{
  return -0.5*(m0-m)*(m0-m) - 0.5*log(2.0*M_PI);
}

/* ********************************************************************************** */
/*                                                                                    */
/*                            Waveform basis functions                                */
/*                                                                                    */
/* ********************************************************************************** */

void LPFImpulseResponse(double **h, struct Data *data, struct Spacecraft *lpf, struct Source *source)
{
  int i,j,n;
  double *P;
  double *r = source->r;
  double *e = malloc(3*sizeof(double));
  double *d = malloc(3*sizeof(double));
  double *t = malloc(3*sizeof(double));

  double sinphi   = sin(source->phi);
  double cosphi   = cos(source->phi);
  double costheta = source->costheta;
  double sintheta = sin(acos(source->costheta));
  
  double *h_norm = malloc(data->N*2*sizeof(double));


  //calculate unit vector to source
  e[0] = -sintheta*cosphi;
  e[1] = -sintheta*sinphi;
  e[2] = -costheta;

  //calculate normalized waveform
  SineGaussianFourier(h_norm, source->t0, 1.0, data->N, 0, data->T);

  P = malloc(data->DOF*sizeof(double));
  

  for(i=0; i<3; i++) P[i] = e[i]*source->P*PC/lpf->M;

  if(data->DOF>3)
  {
    //get vector from proof mass to impact location
    for(i=0; i<3; i++) d[i] = r[i] - lpf->RB[i];

    //get torque direction about proof mass t = (r-R) x e
    crossproduct(d,e,t);


    //get angular acceleration omega = (I^-1)t
    for(i=0; i<3; i++)
    {
      P[i+3] = 0.0;
      for(j=0; j<3; j++) P[i+3] += t[j]*lpf->invI[0][i][j]*source->P*PC;
    }

    //add linear momentum at test mass
    for(i=0; i<3; i++) d[i] = lpf->RTM[0][i] - lpf->RB[i];

    //t = (rTB - rB) x h_theta
    crossproduct(d,P+3,t);
    for(i=0; i<3; i++) P[i] += t[i];

  }


  for(i=0; i<data->DOF; i++)
  {
    //SineGaussianFourier(h[i], source->t0, P[i], data->N, 0, data->T);
    for(n=0; n<data->N*2; n++)
    {
      h[i][n] = P[i]*h_norm[n];
    }
  }
  
  free(h_norm);
  free(P);
  free(e);
  free(d);
  free(t);
}

void SineGaussianFourier(double *hs, double t0, double P, int N, int flag, double Tobs)
{
  double f0, Q, sf, sx, Amp;
  double phi, f;
  double tau;
  double re,im;
  double invTobs = 1.0/Tobs;
  
  int i,istart,istop,imin,imax,even,odd;
  
  double TPI = 2.0*M_PI;
  double fmin = 0.0;
  
  
  //t0  = sigpar[0];
  f0  = (double)N/2.0/Tobs;//sigpar[1];
  Q   = 0.0;//sigpar[2];
    Amp = 1.0;//P*PC;//sigpar[3];
  phi = 0.0;//sigpar[4];
  
  tau = Q/(TPI*f0);
  
  imin = 0;
  imax = N;
  
  istart = 0;
  istop  = N;
  
  
  /* Use recursion relationship for cos(Phase) & sin(Phase) */
  
  //initial values of exp(iPhase)
  double phase = TPI*fmin*t0;
  double cosPhase_m  = cos(phase-phi);
  double sinPhase_m  = sin(phase-phi);
  double cosPhase_p  = cos(phase+phi);
  double sinPhase_p  = sin(phase+phi);
  
  //incremental values of exp(iPhase)
  double dim = -sin(TPI*t0*invTobs);
  double dre = sin(0.5*(TPI*t0*invTobs));
  dre = -2.0*dre*dre;
  
  double amplitude = Amp*sqrt(invTobs);
  double pi2tau2   = M_PI*M_PI*tau*tau;
  double Q2        = Q*Q/f0;
  
  for(i = istart; i < imin; i++)
  {
    even = 2*i;
    odd = even + 1;
    
    if(flag==0) hs[even] = hs[odd] = 0.0;
    
  }
  
  int iend = imax;
  for(i = imin; i < iend; i++)
  {
    even = 2*i;
    odd = even + 1;
    f = (double)(i)*invTobs;
    
    if(flag==0) hs[even] = hs[odd] = 0.0;
    
    sf = amplitude*expf(-pi2tau2*(f-f0)*(f-f0));
    sx = expf(-Q2*f);
    re =  sf*(cosPhase_m+sx*cosPhase_p);//cos(phi-TPI*f*(t0-Tobs));
    im = -sf*(sinPhase_m+sx*sinPhase_p);//sin(phi-TPI*f*(t0-Tobs));
    
    switch(flag)
    {
      case 1: // Add new wavelet to hs[]
        hs[even] += re;
        hs[odd]  += im;
        break;
      case -1: // Remove new wavelet from hs[]
        hs[even] -= re;
        hs[odd]  -= im;
        break;
      case 0:  // Replace hs[] with new wavelet
        hs[even] = re;
        hs[odd]  = im;
        break;
      default:
        fprintf(stderr,"Unsupported SineGaussian flag\n");
        abort();
        break;
    }
    
    /* Now update re and im for the next iteration. */
    recursive_phase_evolution(dre, dim, &cosPhase_m, &sinPhase_m);
    recursive_phase_evolution(dre, dim, &cosPhase_p, &sinPhase_p);
  }
  
  for(i = iend; i < istop; i++)
  {
    even = 2*i;
    odd = even + 1;
    
    if(flag==0) hs[even] = hs[odd] = 0.0;
    
  }
}

void recursive_phase_evolution(double dre, double dim, double *cosPhase, double *sinPhase)
{
  /* Update re and im for the next iteration. */
  double cosphi = *cosPhase;
  double sinphi = *sinPhase;
  
  double newRe = cosphi + cosphi*dre - sinphi*dim;
  double newIm = sinphi + cosphi*dim + sinphi*dre;
  
  *cosPhase = newRe;
  *sinPhase = newIm;
  
}

/* ********************************************************************************** */
/*                                                                                    */
/*                                Likelihood Functions                                */
/*                                                                                    */
/* ********************************************************************************** */

void max_loglikelihood(struct Data *data, struct Spacecraft *lpf, struct Model *model)
{
  int i,k,n;
  
  double **Snf = malloc(3*sizeof(double *));
  for(i=0; i<3; i++) Snf[i] = malloc(data->N*sizeof(double));
  
  // calculate noise model
  for(i=0; i<data->N; i++)data->f[i] = (double)i*data->df;
  Sn(data, lpf, model, Snf);
  
  double **r = malloc(3*sizeof(double *));
  for(k=0; k<3; k++) r[k] = malloc(data->N*2*sizeof(double));
  
  // initialize residual
  for(i=0; i<2*data->N; i++)
  {
    for(k=0; k<3; k++)
    {
      r[k][i] = data->d[k][i];
    }
  }
  
  
  // maximized parameters
  double dt=0.0;
  //double dphi=0.0;
  
  
  for(n=0; n<model->N; n++)
  {
    // calculate signal model
    //SineGaussianFourier(model->s, model->source[n]->t0, model->source[n]->P, data->N, 0, data->T);
    LPFImpulseResponse(model->s, data, lpf, model->source[n]);
    
    
    // find maximum time shift
    //Sum_Extreme(r, model->s, Snf, 2*data->N, &dt, &dphi, data->T, model->source[n]->t0);
    
    model->source[n]->t0 += dt;
    
    // re-calculate signal model
    //SineGaussianFourier(model->s, model->source[n]->t0, model->source[n]->P, data->N, 0, data->T);
    LPFImpulseResponse(model->s, data, lpf, model->source[n]);
    
    // form up residual for next parameter set
    for(i=0; i<2*data->N; i++)
    {
      for(k=0; k<data->DOF; k++) r[k][i] -= model->s[k][i];
    }
    
  }
  
  for(i=0; i<data->DOF; i++)free(Snf[i]);
  free(Snf);
  for(k=0; k<data->DOF; k++)free(r[k]);
  free(r);
}

double loglikelihood(struct Data *data, struct Spacecraft *lpf, struct Model *model, struct Flags *flags)
{
  //return 0;
  int i,k,n,re,im;
  double logL = 0.0;
  
  // calculate noise model
  double **Snf = malloc(data->DOF*sizeof(double *));
  for(i=0; i<data->DOF; i++) Snf[i] = malloc(data->N*sizeof(double));
  for(i=0; i<data->N; i++)data->f[i] = (double)i*data->df;
  //Sn(data, lpf, model, Snf);
  for(k=0; k<data->DOF; k++) for(i=0; i<data->N; i++) Snf[k][i] = model->Snf[k][i];


  double **r = malloc(data->DOF*sizeof(double *));
  for(k=0; k<data->DOF; k++) r[k] = malloc(data->N*2*sizeof(double));
  
  
  model->logL = 0.0;
  
  // initialize residual
  for(i=0; i<2*data->N; i++)
  {
    for(k=0; k<data->DOF; k++)
    {
      r[k][i] = data->d[k][i];
    }
  }
  
  
  // calculate residual of signal model
  for(n=0; n<model->N; n++)
  {
    //SineGaussianFourier(data->s, model->source[n]->t0, model->source[n]->P, data->N, 0, data->T);
    
    //bail if we are not on the surface of LPF
    if(model->source[n]->face==-1)
    {
      for(k=0; k<data->DOF; k++)
      {
        free(r[k]);
        free(Snf[k]);
      }
      free(r);
      free(Snf);
      
      return -1.0e60;
    }
    
    
    LPFImpulseResponse(data->s, data, lpf, model->source[n]);
    
    for(i=0; i<data->N; i++)
    {
      re = 2*i;
      im = re+1;
      
      for(k=0; k<data->DOF; k++)
      {
        r[k][re] -= data->s[k][re];
        r[k][im] -= data->s[k][im];
      }
    }
    
  }
  
  if(!flags->prior)
  {
    for(k=0; k<data->DOF; k++)
    //for(k=0; k<3; k++)
    {
      logL += -0.5*fourier_nwip(data->imin, data->imax, r[k], r[k], Snf[k]) + loglike_normalization(data->imin, data->imax, Snf[k]);
    }
  }
  
  for(k=0; k<data->DOF; k++)
  {
    free(r[k]);
    free(Snf[k]);
  }
  free(r);
  free(Snf);
  
  
  return logL;
  
}

double loglike_normalization(int imin, int imax, double *Sn)
{
  int i;
  double norm;
  
  norm = 0.0;
  
  for(i=imin; i<imax; i++)
  {
    norm -= log(Sn[i])*0.5;
  }
  
  return(norm);
}

/* ********************************************************************************** */
/*                                                                                    */
/*                             Instrument noise routines                              */
/*                                                                                    */
/* ********************************************************************************** */

double InertialSensorNoise(double f, double A, double M)
{
  //A is in units of m * Hz^-1/2
  //M is in units of kg
  
  double twopif = 2.0 * M_PI * f;
  double noise = A * M * twopif * twopif;
  return noise*noise;
}

double ThrusterNoise(double f, double A)
{
  //A is in units of Newtons * Hz^-1/2
  return A*A;
}

double AngularSensingNoise(double f, double A, double I)
{
  //A is in units of m * Hz^-1/2
  //M is in units of kg
  
  double twopif = 2.0 * M_PI * f;
  double noise = A * I * twopif * twopif;
  return noise*noise;
}


void Sn(struct Data *data, struct Spacecraft *lpf, struct Model *model, double **Snf)
{
  int i,n;
  for(i=0; i<data->DOF; i++)
  {
    for(n=0; n<data->N; n++)
    {
      //linear d.o.f.
      if(i<3)Snf[i][n] = InertialSensorNoise(data->f[n], model->Ais[i],lpf->M) + ThrusterNoise(data->f[n],model->Ath[i]);
      
      //angular d.o.f.
      else   Snf[i][n] = AngularSensingNoise(data->f[n], model->Ars[i-3],lpf->I[0][0][0]) + ThrusterNoise(data->f[n],model->Ath[i-3]*0.5);
    }
  }
}


void setup_psd_histogram(struct Data *data, struct Spacecraft *lpf, struct Model *model, struct PSDposterior *psd)
{
  int i,j;
  
  /* Create arrays to compute PSD 2D histogram */
  double **Snf = malloc(3*sizeof(double *));
  for(i=0; i<3; i++) Snf[i] = malloc(data->N*sizeof(double));
  Sn(data, lpf, model, Snf);
  double logSn;
  double Snfmin =  1.0e60;
  double Snfmax = -1.0e60;
  
  psd->Nx = 300;
  psd->Ny = 300;
  for(i=0; i<data->N; i++)
  {
    
    logSn = log10(sqrt(Snf[0][i]));
    
    if(logSn<Snfmin) Snfmin = logSn;
    if(logSn>Snfmax) Snfmax = logSn;
    
  }
  
  psd->ymin = Snfmin - 0.1*(Snfmax-Snfmin)/2.0;
  psd->ymax = Snfmax + 0.1*(Snfmax-Snfmin)/2.0;
  psd->dy = (psd->ymax-psd->ymin)/(double)psd->Ny;
  
  psd->xmin = log10(data->fmin);
  psd->xmax = log10(data->fmax);
  psd->dx = (psd->xmax-psd->xmin)/(double)psd->Nx;
  
  
  psd->histogram = malloc(psd->Nx*sizeof(double *));
  for(i=0; i<psd->Nx; i++)
  {
    psd->histogram[i] = malloc(psd->Ny*sizeof(double));
    for(j=0; j<psd->Ny; j++) psd->histogram[i][j] = 0.0;
  }
  
}

void populate_psd_histogram(struct Data *data, struct Spacecraft *lpf, struct Model *model, int MCMCSTEPS, struct PSDposterior *psd)
{
  int i,iy;
  double f;
  double logSn;
  
  for(i=0; i<psd->Nx; i++)
  {
    f = pow(10,psd->xmin+i*psd->dx);
    
    logSn = log10(sqrt(InertialSensorNoise(f, model->Ais[0],lpf->M) + ThrusterNoise(f,model->Ath[0])));
    
    iy = (int)floor((logSn-psd->ymin)/psd->dy);
    
    psd->histogram[i][iy] += 1.0/((double)MCMCSTEPS/2.0);
  }
}

void print_time_domain_waveforms(char filename[], double *h, int N, double *Snf, double eta, double Tobs, int imin, int imax, double tmin, double tmax)
{
  /*
   TODO: invFFT comes out time-reversed.
   LAL uses exp(-i2pift) for their Fourier transforms
   Why doesn't drealft(ty-1,N,1) work for inverse FFT?
   In the meantime, hacked invFFT by changing sign of imaginary part
   */
  int i;
  double x,t;
  FILE *waveout = fopen(filename,"w");

  double *ht = malloc(N*sizeof(double));
  for(i=0; i<N; i++)ht[i]=h[i];

  ht[0] = 0.0;
  ht[1] = 0.0;
  for(i=1; i< N/2; i++)
  {
    if(i>imin && i<imax)
    {
      x = sqrt(Snf[i]*eta);
      ht[2*i] /= x;
      ht[2*i+1] /= x;
      //printf("x=%g, ht[%i/%i] = %g + i%g\n",x,i,4096,ht[2*i],ht[2*i+1]);
    }
    else
    {
      ht[2*i]=ht[2*i+1]=0.0;
    }
  }

  drealft(ht-1,N,-1);
  double norm = 0.5*sqrt((double)N);

  double t0 = 0.0;//Tobs-2.0;

  for(i=0; i<N; i++)
  {
    t = (double)(i)/(double)(N)*Tobs;
    if(t>=tmin && t<tmax)fprintf(waveout,"%e %e\n", t, ht[i]/norm);
  }

  free(ht);

  fclose(waveout);
}


/* ********************************************************************************** */
/*                                                                                    */
/*                                    Math tools                                      */
/*                                                                                    */
/* ********************************************************************************** */

double snr(struct Data *data, struct Spacecraft *lpf, struct Model *model)
{
  int i,n,k;
  
  double **Snf = malloc(data->DOF*sizeof(double *));
  for(i=0; i<data->DOF; i++) Snf[i] = malloc(data->N*sizeof(double));
  
  // calculate noise model
  Sn(data, lpf, model, Snf);
  
  
  // calculate signal model
  double SNR = 0.0;
  
  for(n=0; n<model->N; n++)
  {
    //SineGaussianFourier(model->s, model->source[n]->t0, model->source[n]->P, data->N, 0, data->T);
    LPFImpulseResponse(model->s, data, lpf, model->source[n]);
    
    for(k=0; k<data->DOF; k++) SNR += fourier_nwip(data->imin, data->imax, model->s[k], model->s[k], Snf[k]);
  }
  
  for(i=0; i<data->DOF; i++)free(Snf[i]);
  free(Snf);
  
  return sqrt(SNR);
  
}

double fourier_nwip(int imin, int imax, double *a, double *b, double *Sn)
{
  int i, j, k;
  double arg, product;
  double ReA, ReB, ImA, ImB;
  
  arg = 0.0;
  
  for(i=imin; i<imax; i++)
  {
    j = i * 2;
    k = j + 1;
    ReA = a[j]; ImA = a[k];
    ReB = b[j]; ImB = b[k];
    product = ReA*ReB + ImA*ImB;
    arg  += product/Sn[i];
  }
  return(2.0*(arg));
  
}

void crossproduct(double *b, double *c, double *a)
{
  a[0] = b[1]*c[2] - b[2]*c[1];
  a[1] = b[2]*c[0] - b[0]*c[2];
  a[2] = b[0]*c[1] - b[1]*c[0];
}

void matrix_multiply(double **A, double **B, double **C, int N)
{
  int i,j,k;

  for(i=0; i<N; i++)
  {
    for(j=0; j<N; j++)
    {
      C[i][j] = 0.0;
      for(k=0; k<N; k++)
      {
        C[i][j] += A[i][k]*B[k][j];
      }
    }
  }
}

void matrix_invert(double **A, double **invA, int N)
{
  int i,j;

  // Find eigenvectors and eigenvalues
  gsl_matrix *GSLmatrx = gsl_matrix_alloc(N,N);
  gsl_matrix *GSLinvrs = gsl_matrix_alloc(N,N);

  for(i=0; i<N; i++)for(j=0; j<N; j++) gsl_matrix_set(GSLmatrx,i,j,A[i][j]);

  // sort and put them into evec
  gsl_permutation * permutation = gsl_permutation_alloc(N);
  gsl_linalg_LU_decomp(GSLmatrx, permutation, &i);
  gsl_linalg_LU_invert(GSLmatrx, permutation, GSLinvrs);

  //unpack arrays from gsl inversion
  for(i=0; i<N; i++)
  {
    for(j=0; j<N; j++) invA[i][j] = gsl_matrix_get(GSLinvrs,i,j);
  }

//  //test inversion
//  double **I = malloc(N*sizeof(double *));
//  for(i=0; i<N; i++) I[i] = malloc(N*sizeof(double));
//
//  matrix_multiply(A,invA,I,N);
//  for(i=0; i<N; i++)
//  {
//    for(j=0; j<N; j++)
//    {
//      printf("%.0f ",I[i][j]);
//    }
//    printf("\n");
//  }
//
//  for(i=0; i<N; i++) free(I[i]);
//  free(I);

}

void dfour1(double data[], unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi, swap;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m > 1 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

void drealft(double data[], unsigned long n, int isign)
{
  void dfour1(double data[], unsigned long nn, int isign);
  unsigned long i,i1,i2,i3,i4,np3;
  double c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;

  theta=3.141592653589793/(double) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    dfour1(data,n>>1,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for (i=2;i<=(n>>2);i++) {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    dfour1(data,n>>1,-1);
  }
}

/* ********************************************************************************** */
/*                                                                                    */
/*                        Data handling and injection routines                        */
/*                                                                                    */
/* ********************************************************************************** */

void simulate_injection(struct Data *data, struct Spacecraft *lpf, struct Model *injection)
{
  int n,i,k;
  int re,im;
  
  for(n=0; n<injection->N; n++)
  {
    
    LPFImpulseResponse(data->s, data, lpf, injection->source[n]);
    
    for(i=0; i<data->N; i++)
    {
      re = 2*i;
      im = re+1;
      
      for(k=0; k<data->DOF; k++)
      {
        data->d[k][re] += data->s[k][re];
        data->d[k][im] += data->s[k][im];
      }
    }
  }
  
  
  FILE *dataFile = fopen("data.dat","w");
  
  for(i=0; i<data->N; i++)
  {
    re = 2*i;
    im = re+1;
    
    data->f[i] = (double)i*data->df;
    
    fprintf(dataFile,"%lg ",data->f[i]);
    for(k=0; k<data->DOF; k++)fprintf(dataFile,"%lg ",data->d[k][re]*data->d[k][re] + data->d[k][im]*data->d[k][im]);
    fprintf(dataFile,"\n");
  }
  
  fclose(dataFile);
}

void simulate_data(struct Data *data)
{
  int i,k;
  int re,im;
  
  for(i=0; i<data->N; i++)
  {
    re = 2*i;
    im = re+1;
    
    data->f[i] = (double)i*data->df;
    
    for(k=0; k<data->DOF; k++)
    {
      data->d[k][re] = data->n[k][re];
      data->d[k][im] = data->n[k][im];
    }
  }
  
}

void simulate_noise(struct Data *data, struct Spacecraft *lpf, struct Model *injection, gsl_rng *r)
{
  int i,k;
  int re,im;
  
  FILE *dataFile = fopen("psd.dat","w");
  
  double **Snf = malloc(data->DOF*sizeof(double *));
  for(i=0; i<data->DOF; i++) Snf[i] = malloc(data->N*sizeof(double));
  
  // calculate noise model
  for(i=0; i<data->N; i++)data->f[i] = (double)i*data->df;
  Sn(data, lpf, injection, Snf);
  
  
  for(i=0; i<data->N; i++)
  {
    re = 2*i;
    im = re+1;
    
    fprintf(dataFile,"%lg  ",data->f[i]);
    
    for(k=0; k<data->DOF; k++)
    {
      data->n[k][re] = 0.5*gsl_ran_ugaussian(r)*sqrt(Snf[k][i]);
      data->n[k][im] = 0.5*gsl_ran_ugaussian(r)*sqrt(Snf[k][i]);
      fprintf(dataFile,"%lg %lg %lg %lg ",data->n[k][re]*data->n[k][re] + data->n[k][im]*data->n[k][im],InertialSensorNoise(data->f[i], injection->Ais[k],lpf->M),ThrusterNoise(data->f[i],injection->Ath[k]),Snf[k][i]);
      
    }
    fprintf(dataFile,"\n");
  }
  
  fclose(dataFile);
  
  for(i=0; i<data->DOF; i++)free(Snf[i]);
  free(Snf);
  
}

/* ********************************************************************************** */
/*                                                                                    */
/*                           Memory (de)allocation routines                           */
/*                                                                                    */
/* ********************************************************************************** */

void copy_source(struct Source *source, struct Source *copy, int DOF)
{
  copy->P        = source->P;
  copy->t0       = source->t0;
  copy->costheta = source->costheta;
  copy->phi      = source->phi;
  copy->face     = source->face;
  
  int i;
  for(i=0; i<2; i++) copy->map[i] = source->map[i];
  
  for(i=0; i<3; i++)
  {
    copy->r[i]     = source->r[i];
    copy->n[i]     = source->n[i];
    copy->k[i]     = source->k[i];
    copy->omega[i] = source->omega[i];
  }
}

void copy_model(struct Model *model, struct Model *copy, int Ndata, int DOF)
{

  //copy->mass = model->mass;
  copy->logL = model->logL;
  copy->logP = model->logP;
  
  copy->N = model->N;
  
  int n;
  for(n=0; n<model->N; n++) copy_source(model->source[n],copy->source[n], DOF);
  
  int i,k;

  for(k=0; k<3; k++)
  {
    copy->Ais[k]  = model->Ais[k];
    copy->Ath[k]  = model->Ath[k];
    copy->Ars[k]  = model->Ars[k];
  }

//  //moment of inertia tensors
//  for(i=0; i<2; i++) // test masses
//  {
//    for(j=0; j<3; j++) // x
//    {
//      for(k=0; k<3; k++) //y
//      {
//        copy->I[i][j][k]    = model->I[i][j][k];
//        copy->invI[i][j][k] = model->invI[i][j][k];
//      }
//    }
//  }

  
  for(k=0; k<DOF; k++)
  {
    for(i=0; i<Ndata*2; i++)
    {
      copy->s[k][i] = model->s[k][i];
    }
  }


  for(k=0; k<DOF; k++)
  {
    for(i=0; i<Ndata; i++)
    {
      copy->Snf[k][i]    = model->Snf[k][i];
      copy->SnS[k][i]    = model->SnS[k][i];
      copy->invSnf[k][i] = model->invSnf[k][i];
    }
  }

}

void initialize_source(struct Source *source)
{
  source->map   = malloc(2*sizeof(double));
  source->r     = malloc(3*sizeof(double));
  source->n     = malloc(3*sizeof(double));
  source->k     = malloc(3*sizeof(double));
  source->omega = malloc(3*sizeof(double));
}


void initialize_model(struct Model *model, int Ndata, int D, int DOF)
{
  int i;
  model->Snf    = malloc(DOF*sizeof(double *));
  model->invSnf = malloc(DOF*sizeof(double *));
  model->SnS    = malloc(DOF*sizeof(double *));
  for(i=0; i<DOF; i++)
  {
    model->Snf[i]    = malloc(Ndata*sizeof(double));
    model->invSnf[i] = malloc(Ndata*sizeof(double));
    model->SnS[i]    = malloc(Ndata*sizeof(double));
  }
  model->s    = malloc(DOF*sizeof(double *));
  for(i=0; i<DOF; i++) model->s[i] = malloc(Ndata*2*sizeof(double));
  
  model->Ais  = malloc(3*sizeof(double));
  model->Ath  = malloc(3*sizeof(double));
  model->Ars  = malloc(3*sizeof(double));
  
  model->N = D;
  model->source = malloc(model->N*sizeof(struct Source*));
  for(i=0; i<model->N; i++)
  {
    model->source[i] = malloc(sizeof(struct Source));
    initialize_source(model->source[i]);
  }

//  model->I    = malloc(2*sizeof(double **));
//  model->invI = malloc(2*sizeof(double **));
//  for(j=0; j<2; j++) model->I[j]    = malloc(3*sizeof(double *));
//  for(j=0; j<2; j++) model->invI[j] = malloc(3*sizeof(double *));
//  for(j=0; j<2; j++) for(i=0; i<3; i++) model->I[j][i]    = malloc(N*sizeof(double));
//  for(j=0; j<2; j++) for(i=0; i<3; i++) model->invI[j][i] = malloc(N*sizeof(double));
}

void initialize_bayesline(struct BayesLineParams **bayesline, struct Data *data, double **psd)
{
  int i,ifo;

  int N    = (int)(data->T*(data->fmax-data->fmin))+1;
  int imin = (int)(data->T*data->fmin);

  for(ifo=0; ifo<data->DOF; ifo++)
  {
    bayesline[ifo] = malloc(sizeof(struct BayesLineParams));
    bayesline[ifo]->priors = malloc(sizeof(BayesLinePriors));
    //bayesline[ifo]->priors->invsigma = malloc((int)(Tobs*(fmax-fmin))*sizeof(double));
    bayesline[ifo]->priors->upper = malloc(N*sizeof(double));
    bayesline[ifo]->priors->lower = malloc(N*sizeof(double));
    bayesline[ifo]->priors->mean  = malloc(N*sizeof(double));
    bayesline[ifo]->priors->sigma = malloc(N*sizeof(double));

    //Set BayesLine priors based on the data channel being used
    bayesline[ifo]->priors->SAmin = 1.0e-18;
    bayesline[ifo]->priors->SAmax = 1.0;
    bayesline[ifo]->priors->LQmin = 1.0e2;
    bayesline[ifo]->priors->LQmax = 1.0e8;
    bayesline[ifo]->priors->LAmin = 1.0e-18;
    bayesline[ifo]->priors->LAmax = 1.0;

    // set default flags
    //bayesline[ifo]->constantLogLFlag = data->constantLogLFlag;

    //Use initial estimate of PSD to set priors
    for(i=0; i<N; i++)
    {
      //bayesline[ifo]->priors->invsigma[i] = 1./(psd[ifo][i+(int)(Tobs*fmin)]*1.0);
      bayesline[ifo]->priors->sigma[i] = psd[ifo][i+imin];
      bayesline[ifo]->priors->mean[i]  = psd[ifo][i+imin];
      bayesline[ifo]->priors->lower[i] = psd[ifo][i+imin]/10.;
      bayesline[ifo]->priors->upper[i] = psd[ifo][i+imin]*10.;

    }

    //Allocates arrays and sets constants for for BayesLine
    BayesLineSetup(bayesline[ifo], data->d[ifo], data->fmin, data->fmax, data->dt, data->T);
  }
}

void free_source(struct Source *source)
{
  free(source->map);
  free(source->r);
  free(source->n);
  free(source->k);
  free(source->omega);
  free(source);
}

void free_model(struct Model *model, int N, int D, int DOF)
{
  int i;
  
  for(i=0; i<N; i++) free_source(model->source[i]);
  free(model->Ath);
  free(model->Ais);
  for(i=0; i<DOF; i++) free(model->s[i]);
  free(model->s);

//  for(j=0; j<2; j++) for(i=0; i<3; i++) free(model->I[j][i]);
//  for(j=0; j<2; j++) for(i=0; i<3; i++) free(model->invI[j][i]);
//  for(j=0; j<2; j++) free(model->I[j]);
//  for(j=0; j<2; j++) free(model->invI[j]);
//  free(model->I);
//  free(model->invI);


  free(model);
}

