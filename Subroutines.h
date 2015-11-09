/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/* ********************************************************************************** */
/*                                                                                    */
/*                                  Data structures                                   */
/*                                                                                    */
/* ********************************************************************************** */

#define PC 3.6e-8

struct Source
{
  double t0;
  double P;
};

struct Model
{
  double mass;
  double Ais;
  double Ath;
  double logL;
  double logP;

  int N;

  struct Source **source;

  double *s;
};

struct Chain
{
  int N;
  struct Model **model;
  double *temp;
  int *index;

};

struct Data
{
  int N;
  int imin;
  int imax;
  double T;
  double dt;
  double df;
  double fmax;
  double fmin;
  double *d;
  double *n;
  double *s;
  double *f;
};

/* ********************************************************************************** */
/*                                                                                    */
/*                                    MCMC tools                                      */
/*                                                                                    */
/* ********************************************************************************** */

void proposal(struct Data *data, struct Model *model, struct Model *trial, gsl_rng *r);

void logprior(struct Data *data, struct Model *model, struct Model *injection);

double log_mass_prior(double m0, double m);

/* ********************************************************************************** */
/*                                                                                    */
/*                            Waveform basis functions                                */
/*                                                                                    */
/* ********************************************************************************** */

void SineGaussianFourier(double *hs, double t0, double P, int N, int flag, double Tobs);

void recursive_phase_evolution(double dre, double dim, double *cosPhase, double *sinPhase);

/* ********************************************************************************** */
/*                                                                                    */
/*                                Likelihood Functions                                */
/*                                                                                    */
/* ********************************************************************************** */

void max_loglikelihood(struct Data *data, struct Model *model);

void loglikelihood(struct Data *data, struct Model *model);

double loglike_normalization(int imin, int imax, double *Sn);

/* ********************************************************************************** */
/*                                                                                    */
/*                             Instrument noise routines                              */
/*                                                                                    */
/* ********************************************************************************** */

double InertialSensorNoise(double f, double A, double M);

double ThrusterNoise(double f, double A);

double Sn(double f, struct Model *model);

/* ********************************************************************************** */
/*                                                                                    */
/*                                    Math tools                                      */
/*                                                                                    */
/* ********************************************************************************** */

double snr(struct Data *data, struct Model *model);

double fourier_nwip(int imin, int imax, double *a, double *b, double *Sn);

/* ********************************************************************************** */
/*                                                                                    */
/*                        Data handling and injection routines                        */
/*                                                                                    */
/* ********************************************************************************** */

void simulate_data(struct Data *data);

void simulate_injection(struct Data *data, struct Model *injection);

void simulate_noise(struct Data *data, struct Model *injection, gsl_rng *r);

/* ********************************************************************************** */
/*                                                                                    */
/*                           Memory (de)allocation routines                           */
/*                                                                                    */
/* ********************************************************************************** */

void copy_model(struct Model *model, struct Model *copy, int N);
