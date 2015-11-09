/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Subroutines.h"
#include "TimePhaseMaximization.h"

/* ============================  MAIN PROGRAM  ============================ */

int main()
{
  /* declare variables */
  int ic,n,mc;
  int accept;
  int MCMCSTEPS;

  double H;
  double alpha;

  struct Data  *data      = malloc(sizeof(struct Data));
  struct Model *trial     = malloc(sizeof(struct Model));
  struct Model *injection = malloc(sizeof(struct Model));
  injection->N = trial->N = 1;
  injection->source = malloc(injection->N*sizeof(struct Source*));
  trial->source     = malloc(injection->N*sizeof(struct Source*));
  for(n=0; n<injection->N; n++)
  {
    injection->source[n] = malloc(sizeof(struct Source));
    trial->source[n]     = malloc(sizeof(struct Source));
  }

  struct Chain *chain =  malloc(sizeof(struct Chain));
  chain->N = 10;
  chain->temp  = malloc(sizeof(double)*chain->N);
  chain->model = malloc(sizeof(struct Model*)*chain->N);
  for(ic=0; ic<chain->N; ic++)
  {
    chain->temp[ic] = (double)ic;
    chain->model[ic] = malloc(sizeof(struct Model));
    chain->model[ic]->N  = 2;
    chain->model[ic]->source = malloc(chain->model[ic]->N*sizeof(struct Source*));

    for(n=0; n<chain->model[ic]->N; n++)
    {
      chain->model[ic]->source[n] = malloc(sizeof(struct Source));
    }
  }

  /* set up GSL random number generator */
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc (T);
  gsl_rng_env_setup();

  /* Initialize data structure */
  data->T  = 2048.;
  data->dt = 0.5;
  data->df = 1.0/data->T;
  data->N  = (int)(data->T/data->dt)/2;

  data->fmin = 1.0e-4; //Hz
  data->fmax = (double)data->N/data->T;  //Hz

  data->imin = (int)floor(data->fmin*data->T);
  data->imax = (int)floor(data->fmax*data->T);

  data->d = malloc(data->N*2*sizeof(double));
  data->s = malloc(data->N*2*sizeof(double));
  data->n = malloc(data->N*2*sizeof(double));
  data->f = malloc(data->N*sizeof(double));

  /* Simulate noise data */
  injection->mass = 422.0; //kg
  injection->Ais  = 2.0e-9; // m*Hz^-1/2
  injection->Ath  = 1.0e-8; // N*Hz^-1/2

  /* Simulate source data */
  injection->source[0]->t0   = data->T/4.0;
  injection->source[0]->P    = 8.0;

//  injection->source[1]->t0   = data->T/8.0;
//  injection->source[1]->P    = 8.0;

  injection->s = malloc(data->N*2*sizeof(double));
  trial->s     = malloc(data->N*2*sizeof(double));
  for(ic=0; ic<chain->N; ic++)
    chain->model[ic]->s = malloc(data->N*2*sizeof(double));

  simulate_noise(data, injection, r);
  simulate_data(data);
  simulate_injection(data,injection);
  loglikelihood(data, injection);

  printf("Injected parameters:   {%lg,%lg}\n",injection->source[0]->t0 ,injection->source[0]->P );
  printf("SNR of injection = %g\n",snr(data,injection));

  /* Initialize model */
  trial->N = 1;
  for(ic=0; ic<chain->N; ic++)
  {
    chain->model[ic]->N = 1;
    proposal(data,injection,chain->model[ic],r);

    logprior(data, chain->model[ic], injection);
    loglikelihood(data, chain->model[ic]);
  }

  FILE *outfile;


  /* set up distribution */

  /* set up MCMC run */
  accept    = 0;
  MCMCSTEPS = 1000000;
  outfile = fopen("chain.dat","w");
  fprintf(outfile,"#dlogL mass dAx dAf t0 P\n");

  /* Here is the MCMC loop */
  for(mc=0;mc<MCMCSTEPS;mc++)
  {
    //copy x to y
    copy_model(chain->model[0], trial, data->N);

    //choose new parameters for y
    proposal(data, chain->model[0], trial, r);

    //compute maximized likelihood
    //max_loglikelihood(data, trial);

    //compute new likelihood
    loglikelihood(data, trial);

    //compute new prior
    logprior(data, trial, injection);

    //compute Hastings ratio
    H     = trial->logL - chain->model[0]->logL + trial->logP - chain->model[0]->logP;
    alpha = log(gsl_rng_uniform(r));

    //adopt new position w/ probability H
    if(H>alpha)
    {
      copy_model(trial, chain->model[0], data->N);
      accept++;
    }

    //output to file
    fprintf(outfile,"%lg %lg %lg %lg ",chain->model[0]->logL-injection->logL,chain->model[0]->mass,(injection->Ais-chain->model[0]->Ais)/injection->Ais,(injection->Ath-chain->model[0]->Ath)/injection->Ath);
    for(n=0; n<chain->model[0]->N; n++)fprintf(outfile,"%lg %lg ", chain->model[0]->source[n]->t0,chain->model[0]->source[n]->P);
    fprintf(outfile,"\n");
  }
  printf("acceptance rate = %g\n",(double)accept/(double)MCMCSTEPS);

	return 0;
}




























