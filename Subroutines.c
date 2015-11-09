/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Subroutines.h"
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

void proposal(struct Data *data, struct Model *model, struct Model *trial, gsl_rng *r)
{
  int n;

  trial->mass = model->mass + gsl_ran_ugaussian(r)*1.0;     // kg
  trial->Ais  = model->Ais  + gsl_ran_ugaussian(r)*2.0e-11; // m*Hz^-1/2
  trial->Ath  = model->Ath  + gsl_ran_ugaussian(r)*1.0e-10; // N*Hz^-1/2

  for(n=0; n<trial->N; n++)
  {
    //uniform
    if(gsl_rng_uniform(r)<0.1)
    {
      trial->source[n]->P    = gsl_rng_uniform(r)*100;
      trial->source[n]->t0   = gsl_rng_uniform(r)*data->T;

    }
    //gaussian
    else
    {
      trial->source[n]->P    = model->source[n]->P  + gsl_ran_ugaussian(r)*1.0;
      trial->source[n]->t0   = model->source[n]->t0 + gsl_ran_ugaussian(r)*1.0;
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

  int iend = imax+1;
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
  int i;

  double *Snf = malloc(data->N*sizeof(double));

  // calculate noise model
  for(i=0; i<data->N; i++)
  {
    Snf[i] = Sn(data->f[i], model);
  }

  // calculate signal model
  SineGaussianFourier(model->s, model->source[0]->t0, model->source[0]->P, data->N, 0, data->T);

  // experiment with maximizaiton
  double dt=0.0;
  double dphi=0.0;
  fprintf(stdout,"2*data->N = %i\n",2*data->N);fflush(stdout);
  Sum_Extreme(data->d, model->s, Snf, 2*data->N, &dt, &dphi, data->T, 0, 2*data->N, model->source[0]->t0);
  model->source[0]->t0 += dt;


  free(Snf);
  
}

void loglikelihood(struct Data *data, struct Model *model)
{
  int i,n,re,im;

  double *Snf = malloc(data->N*sizeof(double));
  double *r   = malloc(data->N*2*sizeof(double));

  // initialize residual
  for(i=0; i<2*data->N; i++)
  {
    r[i] = data->d[i];
  }

  // calculate noise model
  for(i=0; i<data->N; i++)
  {
    Snf[i] = Sn(data->f[i], model);
  }

  // calculate residual of signal model
  for(n=0; n<model->N; n++)
  {
    SineGaussianFourier(data->s, model->source[n]->t0, model->source[n]->P, data->N, 0, data->T);

    for(i=0; i<data->N; i++)
    {
      re = 2*i;
      im = re+1;

      r[i] -= data->s[re];
      r[i] -= data->s[im];

    }
  }

  model->logL =  -0.5*fourier_nwip(data->imin, data->imax, r, r, Snf) + loglike_normalization(data->imin, data->imax, Snf);

  free(r);
  free(Snf);
  
}

double loglike_normalization(int imin, int imax, double *Sn)
{
  int i;
  double norm;

  norm = 0.0;

  for(i=imin; i<imax; i++)
  {
    norm -= log(Sn[i])/2.0;
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

double Sn(double f, struct Model *model)
{
  return InertialSensorNoise(f, model->Ais,model->mass) + ThrusterNoise(f,model->Ath);
}

/* ********************************************************************************** */
/*                                                                                    */
/*                                    Math tools                                      */
/*                                                                                    */
/* ********************************************************************************** */

double snr(struct Data *data, struct Model *model)
{
  int i,n;

  double *Snf = malloc(data->N*sizeof(double));

  // calculate noise model
  for(i=0; i<data->N; i++)
  {
    Snf[i] = Sn(data->f[i], model);
  }

  // calculate signal model
  double SNR = 0.0;

  for(n=0; n<model->N; n++)
  {
    SineGaussianFourier(model->s, model->source[n]->t0, model->source[n]->P, data->N, 0, data->T);
    SNR += fourier_nwip(data->imin, data->imax, model->s, model->s, Snf);
  }
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

/* ********************************************************************************** */
/*                                                                                    */
/*                        Data handling and injection routines                        */
/*                                                                                    */
/* ********************************************************************************** */

void simulate_injection(struct Data *data, struct Model *injection)
{
  int n,i;
  int re,im;

  for(n=0; n<injection->N; n++)
  {
    SineGaussianFourier(data->s, injection->source[n]->t0, injection->source[n]->P, data->N, 0, data->T);

    for(i=0; i<data->N; i++)
    {
      re = 2*i;
      im = re+1;

      data->d[re] += data->s[re];
      data->d[im] += data->s[im];

    }
  }


  FILE *dataFile = fopen("data.dat","w");

  for(i=0; i<data->N; i++)
  {
    re = 2*i;
    im = re+1;

    data->f[i] = (double)i*data->df;

    fprintf(dataFile,"%lg %lg\n",data->f[i],data->d[re]*data->d[re] + data->d[im]*data->d[im]);
  }

  fclose(dataFile);
}

void simulate_data(struct Data *data)
{
  int i;
  int re,im;

  FILE *dataFile = fopen("data.dat","w");

  for(i=0; i<data->N; i++)
  {
    re = 2*i;
    im = re+1;

    data->f[i] = (double)i*data->df;

    data->d[re] = data->n[re];
    data->d[im] = data->n[im];

    fprintf(dataFile,"%lg %lg\n",data->f[i],data->d[re]*data->d[re] + data->d[im]*data->d[im]);

  }

  fclose(dataFile);
}

void simulate_noise(struct Data *data, struct Model *injection, gsl_rng *r)
{
  int i;
  int re,im;
  double Snf;

  FILE *dataFile = fopen("noise.dat","w");

  for(i=0; i<data->N; i++)
  {
    re = 2*i;
    im = re+1;

    data->f[i] = (double)i*data->df;

    Snf = Sn(data->f[i], injection);

    data->n[re] = 0.5*gsl_ran_ugaussian(r)*sqrt(Snf);
    data->n[im] = 0.5*gsl_ran_ugaussian(r)*sqrt(Snf);

    fprintf(dataFile,"%lg %lg %lg %lg %lg\n",data->f[i],data->n[re]*data->n[re] + data->n[im]*data->n[im],InertialSensorNoise(data->f[i], injection->Ais,injection->mass),ThrusterNoise(data->f[i],injection->Ath),Snf);

  }
  
  fclose(dataFile);
}

/* ********************************************************************************** */
/*                                                                                    */
/*                           Memory (de)allocation routines                           */
/*                                                                                    */
/* ********************************************************************************** */

void copy_model(struct Model *model, struct Model *copy, int N)
{
  copy->mass = model->mass;
  copy->Ais  = model->Ais;
  copy->Ath  = model->Ath;
  copy->logL = model->logL;
  copy->logP = model->logP;

  int n;
  for(n=0; n<model->N; n++)
  {
    copy->source[n]->P    = model->source[n]->P;
    copy->source[n]->t0   = model->source[n]->t0;
  }
  int i;
  for(i=0; i<N*2; i++)
  {
    copy->s[i] = model->s[i];
  }
}
