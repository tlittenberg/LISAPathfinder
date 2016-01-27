/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Subroutines.h"
#include "LISAPathfinder.h"
#include "TimePhaseMaximization.h"

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

void draw_impact_point_cube(struct Data *data, struct Source *source, gsl_rng *seed)
{
  //pick where on the surface
  draw_face(source->map, seed);

  //pick sky location
  source->costheta = -1.0 + 2.0*gsl_rng_uniform(seed);
  source->phi      = gsl_rng_uniform(seed)*2.0*M_PI;

  //store incident face
  source->face = which_face(source->map[0],source->map[1]);

  /*
  //get angle w.r.t. normal
  source->omega[0] = sin(acos(source->cosalpha))*cos(source->beta);
  source->omega[1] = sin(acos(source->cosalpha))*sin(source->beta);
  source->omega[2] = source->cosalpha;

   //get true source position
   //rotate_sky_to_3D(source->omega, source->face, &source->costheta, &source->phi);
   */

  //get 3D position
  rotate_face_to_3D(source->map, source->r);

  //flag if face is not exposed to sky location
  if(check_impact(source->costheta, source->phi, source->face)) source->face = -1;


  //momentum and impact time
  source->P  = gsl_ran_exponential(seed,20);
  source->t0 = gsl_rng_uniform(seed)*data->T;
  
}

void draw_impact_point(struct Data *data, struct Source *source, gsl_rng *seed)
{
   double a = 1.0;
   double A = 2.*(1.0 + sqrt(2.))*a*a;

   //pick a side, any side
   source->face = (int)floor(10*gsl_rng_uniform(seed));
   if(source->face == 9)
   {
      draw_octagon(source->r, seed);
      source->r[2] = -fabs(source->r[2]);
   }
   else if(source->face == 8)
   {
      draw_octagon(source->r, seed);
      source->r[2] = fabs(source->r[2]);
   }
   else if(gsl_rng_uniform(seed) < a*a/A)
   {
      draw_side(source->r, source->face, seed);
   }
   else source->face=-1;


   //pick sky location
   source->costheta = -1.0 + 2.0*gsl_rng_uniform(seed);
   source->phi      = gsl_rng_uniform(seed)*2.0*M_PI;


   //flag if face is not exposed to sky location
   if(check_impact(source->costheta, source->phi, source->face)) source->face = -1;


   //momentum and impact time
   source->P  = gsl_ran_exponential(seed,20);
   source->t0 = gsl_rng_uniform(seed)*data->T;
   
}

void proposal(struct Data *data, struct Model *model, struct Model *trial, gsl_rng *r, int *reject)
{
  int n;


  // Always update spacecraft parameters
  detector_proposal(data, model, trial, r);

  //MCMC proposal
  if(gsl_rng_uniform(r)<0.9)
  {
    for(n=0; n<model->N; n++) impact_proposal(data, model->source[n], trial->source[n], r);
  }
  //RJ proposal
  else dimension_proposal(data, model, trial, r, 10, reject);

}

void detector_proposal(struct Data *data, struct Model *model, struct Model *trial, gsl_rng *r)
{
  int i;
  trial->mass = model->mass + gsl_ran_ugaussian(r)*1.0;     // kg
  for(i=0; i<3; i++)
  {
    trial->Ais[i]  = model->Ais[i]  + gsl_ran_ugaussian(r)*2.0e-11; // m*Hz^-1/2
    trial->Ath[i]  = model->Ath[i]  + gsl_ran_ugaussian(r)*1.0e-10; // N*Hz^-1/2
      trial->Ars[i]  = model->Ars[i]  + gsl_ran_ugaussian(r)*1.0e-9; // rad*Hz^-1/2
  }
}

void impact_proposal(struct Data *data, struct Source *model, struct Source *trial, gsl_rng *r)
{

  //uniform
  if(gsl_rng_uniform(r)<0.1)
  {
    draw_impact_point(data,trial,r);

  }
  //gaussian
  else
  {
     draw_impact_point(data,trial,r);

    trial->P  = model->P  + gsl_ran_ugaussian(r)*1.0;
    trial->t0 = model->t0 + gsl_ran_ugaussian(r)*0.25;

    trial->phi      = model->phi + gsl_ran_ugaussian(r)*0.1;
    trial->costheta = model->costheta + gsl_ran_ugaussian(r)*0.1;

//    trial->map[0]   = model->map[0] + gsl_ran_ugaussian(r)*0.1;
//    trial->map[1]   = model->map[1] + gsl_ran_ugaussian(r)*0.1;
//
//    trial->face = which_face(trial->map[0],trial->map[1]);
    if(check_impact(trial->costheta, trial->phi, trial->face)) trial->face = -1;

    //get 3D position
//    rotate_face_to_3D(trial->map, trial->r);

    while(trial->phi > 2.0*M_PI) trial->phi -= 2.0*M_PI;
    while(trial->phi < 0)        trial->phi += 2.0*M_PI;

    while(trial->costheta >  1.0) trial->costheta =  2.0-trial->costheta;
    while(trial->costheta < -1.0) trial->costheta = -2.0-trial->costheta;

  }
}

void dimension_proposal(struct Data *data, struct Model *model, struct Model *trial, gsl_rng *r, int Nmax, int *test)
{
  int n,kill;

  //birth
  if(gsl_rng_uniform(r)<0.5)
  {
    if(model->N==Nmax) *test = 1;
    else
    {
      draw_impact_point(data,trial->source[trial->N],r);
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
      for(n=kill; n<model->N-1; n++) copy_source(model->source[n+1], trial->source[n], data->DOF);
      trial->N = model->N-1;
    }
  }
}

void logprior(struct Data *data, struct Model *model, struct Model *injection)
{
  model->logP = log_mass_prior(injection->mass,model->mass);

  int n;
  for(n=0; n<model->N; n++)
  {

    if(model->source[n]->t0 < 0.0 || model->source[n]->t0 > data->T) model->logP = -1.0e60;

    if(model->source[n]->P < 0.0) model->logP = -1.0e60;
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

void LPFImpulseResponse(double **h, struct Data *data, struct Source *source)
{
  int i,n;
  double *P;
  double *r = source->r;

  double sinphi   = sin(source->phi);
  double cosphi   = cos(source->phi);
  double costheta = source->costheta;
  double sintheta = sin(acos(source->costheta));

  double *h_norm = malloc(data->N*2*sizeof(double));

  //calculate normalized waveform
  SineGaussianFourier(h_norm, source->t0, 1.0, data->N, 0, data->T);


  P = malloc(data->DOF*sizeof(double));

  P[0] = -source->P*sintheta*cosphi;
  P[1] = -source->P*sintheta*sinphi;
  P[2] = -source->P*costheta;

  if(data->DOF>3)
  {
    double *L = malloc(3*sizeof(double));
    crossproduct(r,P,L);
    P[3] = L[0];
    P[4] = L[1];
    P[5] = L[2];

    free(L);
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
  Amp = P*PC;//sigpar[3];
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
    re = sf*(cosPhase_m+sx*cosPhase_p);//cos(phi-TPI*f*(t0-Tobs));
    im = sf*(sinPhase_m+sx*sinPhase_p);//sin(phi-TPI*f*(t0-Tobs));

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

void max_loglikelihood(struct Data *data, struct Model *model)
{
  int i,k,n;

  double **Snf = malloc(3*sizeof(double *));
  for(i=0; i<3; i++) Snf[i] = malloc(data->N*sizeof(double));

  // calculate noise model
  for(i=0; i<data->N; i++)data->f[i] = (double)i*data->df;
  Sn(data, model, Snf);

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
    LPFImpulseResponse(model->s, data, model->source[n]);


    // find maximum time shift
    //Sum_Extreme(r, model->s, Snf, 2*data->N, &dt, &dphi, data->T, model->source[n]->t0);

    model->source[n]->t0 += dt;

    // re-calculate signal model
    //SineGaussianFourier(model->s, model->source[n]->t0, model->source[n]->P, data->N, 0, data->T);
    LPFImpulseResponse(model->s, data, model->source[n]);

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

double loglikelihood(struct Data *data, struct Model *model)
{
  int i,k,n,re,im;
  double logL = 0.0;

  // calculate noise model
  double **Snf = malloc(data->DOF*sizeof(double *));
  for(i=0; i<data->DOF; i++) Snf[i] = malloc(data->N*sizeof(double));
  for(i=0; i<data->N; i++)data->f[i] = (double)i*data->df;
  Sn(data, model, Snf);

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


    LPFImpulseResponse(data->s, data, model->source[n]);

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

  for(k=0; k<data->DOF; k++)
  {
    re = data->N;
    im = re+1;

    //printf("signal[%i]=%g\n",k,data->s[k][re]*data->s[k][re]+data->s[k][im]*data->s[k][im]);
    //printf("data[%i]=%g\n",k,data->d[k][re]*data->d[k][re]+data->d[k][im]*data->d[k][im]);
    //printf("residual[%i]=%g\n",k,r[k][re]*r[k][re]+r[k][im]*r[k][im]);
    //printf("psd[%i]=%g\n",k,Snf[k][data->N/2]);

    logL += -0.5*fourier_nwip(data->imin, data->imax, r[k], r[k], Snf[k]) + loglike_normalization(data->imin, data->imax, Snf[k]);
    //printf("data->DOF[%i] = %g\n",k,logL);
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


void Sn(struct Data *data, struct Model *model, double **Snf)
{
  int i,n;
  for(i=0; i<data->DOF; i++)
  {
    for(n=0; n<data->N; n++)
    {
      //linear d.o.f.
      if(i<3)Snf[i][n] = InertialSensorNoise(data->f[n], model->Ais[i],model->mass) + ThrusterNoise(data->f[n],model->Ath[i]);

      //angular d.o.f.
      else   Snf[i][n] = AngularSensingNoise(data->f[n], model->Ars[i-3],model->I) + ThrusterNoise(data->f[n],model->Ath[i-3]*0.5);
    }
  }
}


void setup_psd_histogram(struct Data *data, struct Model *model, struct PSDposterior *psd)
{
  int i,j;

  /* Create arrays to compute PSD 2D histogram */
  double **Snf = malloc(3*sizeof(double *));
  for(i=0; i<3; i++) Snf[i] = malloc(data->N*sizeof(double));
  Sn(data, model, Snf);
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

void populate_psd_histogram(struct Data *data, struct Model *model, int MCMCSTEPS, struct PSDposterior *psd)
{
  int i,iy;
  double f;
  double logSn;

  for(i=0; i<psd->Nx; i++)
  {
    f = pow(10,psd->xmin+i*psd->dx);

    logSn = log10(sqrt(InertialSensorNoise(f, model->Ais[0],model->mass) + ThrusterNoise(f,model->Ath[0])));

    iy = (int)floor((logSn-psd->ymin)/psd->dy);

    psd->histogram[i][iy] += 1.0/((double)MCMCSTEPS/2.0);
  }
}

/* ********************************************************************************** */
/*                                                                                    */
/*                                    Math tools                                      */
/*                                                                                    */
/* ********************************************************************************** */

double snr(struct Data *data, struct Model *model)
{
  int i,n,k;

  double **Snf = malloc(data->DOF*sizeof(double *));
  for(i=0; i<data->DOF; i++) Snf[i] = malloc(data->N*sizeof(double));

  // calculate noise model
  Sn(data, model, Snf);


  // calculate signal model
  double SNR = 0.0;

  for(n=0; n<model->N; n++)
  {
    //SineGaussianFourier(model->s, model->source[n]->t0, model->source[n]->P, data->N, 0, data->T);
    LPFImpulseResponse(model->s, data, model->source[n]);

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

/* ********************************************************************************** */
/*                                                                                    */
/*                        Data handling and injection routines                        */
/*                                                                                    */
/* ********************************************************************************** */

void simulate_injection(struct Data *data, struct Model *injection)
{
  int n,i,k;
  int re,im;

  for(n=0; n<injection->N; n++)
  {

    LPFImpulseResponse(data->s, data, injection->source[n]);

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

void simulate_noise(struct Data *data, struct Model *injection, gsl_rng *r)
{
  int i,k;
  int re,im;

  FILE *dataFile = fopen("psd.dat","w");

  double **Snf = malloc(data->DOF*sizeof(double *));
  for(i=0; i<data->DOF; i++) Snf[i] = malloc(data->N*sizeof(double));

  // calculate noise model
  for(i=0; i<data->N; i++)data->f[i] = (double)i*data->df;
  Sn(data, injection, Snf);


  for(i=0; i<data->N; i++)
  {
    re = 2*i;
    im = re+1;

    fprintf(dataFile,"%lg  ",data->f[i]);

    for(k=0; k<data->DOF; k++)
    {
      data->n[k][re] = 0.5*gsl_ran_ugaussian(r)*sqrt(Snf[k][i]);
      data->n[k][im] = 0.5*gsl_ran_ugaussian(r)*sqrt(Snf[k][i]);
      fprintf(dataFile,"%lg %lg %lg %lg ",data->n[k][re]*data->n[k][re] + data->n[k][im]*data->n[k][im],InertialSensorNoise(data->f[i], injection->Ais[k],injection->mass),ThrusterNoise(data->f[i],injection->Ath[k]),Snf[k][i]);

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

  for(i=0; i<DOF; i++)
  {
    copy->r[i]     = source->r[i];
    copy->n[i]     = source->n[i];
    copy->k[i]     = source->k[i];
    copy->omega[i] = source->omega[i];
  }
}

void copy_model(struct Model *model, struct Model *copy, int N, int DOF)
{
  copy->I    = model->I;
  copy->mass = model->mass;
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


  for(k=0; k<DOF; k++)
  {
    for(i=0; i<N*2; i++)
    {
      copy->s[k][i] = model->s[k][i];
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

void initialize_model(struct Model *model, int N, int D, int DOF)
{
  int i;
  model->s    = malloc(DOF*sizeof(double *));
  for(i=0; i<DOF; i++) model->s[i] = malloc(N*2*sizeof(double));

  model->Ais  = malloc(3*sizeof(double));
  model->Ath  = malloc(3*sizeof(double));
  model->Ars  = malloc(3*sizeof(double));

  model->N = D;
  model->source = malloc(model->N*sizeof(struct Source*));
  for(i=0; i<N; i++)
  {
    model->source[i] = malloc(sizeof(struct Source));
    initialize_source(model->source[i]);
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

  free(model);
}

