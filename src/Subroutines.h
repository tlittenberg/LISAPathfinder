/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Spacecraft.h"

/* ********************************************************************************** */
/*                                                                                    */
/*                                  Data structures                                   */
/*                                                                                    */
/* ********************************************************************************** */

#define PC 3.6e-8
//#define DOF 3

struct Source
{
  //time of impact
  double t0;

  //total momentum
  double P;

  //2D location of impact on the surface
  double *map;

  //3D location of impact
  double *r;

  //sky location
  double costheta;
  double phi;
  double *omega; //line of site

  //ID for which face
  int face;

  //norm of face
  double *n;

  //sky location w.r.t. norm
  double cosalpha;
  double beta;
  double *k;

};

struct Model
{
  double *Ais; //inertial sensing noise
  double *Ath; //thruster noise
  double *Ars; //rotational sensing noise
  double **Snf;
  double **invSnf;
  double **SnS;
  double logL; //likelihod
  double logP; //prior
  double logQ; //proposal
  double snr;  //snr

  /* Impact parameters */

  int N;       //number of impacts

  struct Source **source;

  /* Instrument response */
  double **s;
};

struct Data
{
  int N;
  int NFFT;
   int DOF;
  int grs;
  int imin;
  int imax;
  long seed;
  long nseed;
  long iseed;
  double T;
  double dt;
  double df;
  double fmax;
  double fmin;
  double tmin;
  double tmax;
  double **d;
  double **n;
  double **s;
  double *f;

  double **t_density;
  double *t_density_max;

    char path[128];
    char gps[16];
    char duration[8];
};

struct PSDposterior
{
  int Nx;
  int Ny;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double dx;
  double dy;
  double **histogram;
};

struct Flags
{
  int verbose;
  int prior;
  int rj;
    int use_spacecraft;
    int simdata;
    int realdata;

};


/* ********************************************************************************** */
/*                                                                                    */
/*                                    MCMC tools                                      */
/*                                                                                    */
/* ********************************************************************************** */

void bayesline_mcmc(struct Data *data, struct Model **model, struct BayesLineParams ***bayesline, int *index, double beta, int ic);

void ptmcmc(struct Model **model, double *temp, int *index, gsl_rng *r, int NC, int mc);

void proposal(struct Flags *flags, struct Data *data, struct Spacecraft *lpf, struct Model *model, struct Model *trial, gsl_rng *r, int *reject, int nmax, int *drew_prior);

void dimension_proposal(struct Flags *flags, struct Data *data, struct Spacecraft *lpf, struct Model *model, struct Model *trial, gsl_rng *r, int Nmax, int *test);

void detector_proposal(struct Data *data, struct Model *model, struct Model *trial, gsl_rng *r);

void impact_proposal(struct Data *data, struct Spacecraft *lpf, struct Source *model, struct Source *trial, gsl_rng *r);
void impact_proposal_sc(struct Data *data, struct Spacecraft *lpf, struct Source *model, struct Source *trial, gsl_rng *r, int *drew_prior);

void draw_impact_point(struct Data *data, struct Spacecraft *lpf, struct Source *source, gsl_rng *seed);
void draw_impact_point_sc(struct Data *data, struct Spacecraft *lpf, struct Source *source, gsl_rng *seed);
void draw_impactor(struct Data *data, struct Source *source, gsl_rng *seed);

void logprior(struct Data *data, struct Model *model, struct Model *injection);
void logprior_sc(struct Data *data,struct Spacecraft *lpf,  struct Model *model, struct Model *injection, int *drew_prior);

double log_mass_prior(double m0, double m);

void find_impacts(double *h, int N, double *Snf, double eta, double Tobs, int imin, int imax, double tmin, double tmax, double *t_density);

/* ********************************************************************************** */
/*                                                                                    */
/*                            Waveform basis functions                                */
/*                                                                                    */
/* ********************************************************************************** */

void LPFImpulseResponse(double **h, struct Data *data, struct Spacecraft *lpf, struct Source *source);

void SineGaussianFourier(double *hs, double t0, double P, int N, int flag, double Tobs);

void recursive_phase_evolution(double dre, double dim, double *cosPhase, double *sinPhase);

/* ********************************************************************************** */
/*                                                                                    */
/*                                Likelihood Functions                                */
/*                                                                                    */
/* ********************************************************************************** */

void max_loglikelihood(struct Data *data, struct Spacecraft *lpf, struct Model *model);

double loglikelihood(struct Data *data, struct Spacecraft *lpf, struct Model *model, struct Flags *flags);

double loglike_normalization(int imin, int imax, double *Sn);

/* ********************************************************************************** */
/*                                                                                    */
/*                             Instrument noise routines                              */
/*                                                                                    */
/* ********************************************************************************** */

double InertialSensorNoise(double f, double A, double M);

double AngularSensingNoise(double f, double A, double I);

double ThrusterNoise(double f, double A);

void Sn(struct Data *data, struct Spacecraft *lpf, struct Model *model, double **Snf);

void setup_psd_histogram(struct Data *data, struct Spacecraft *lpf, struct Model *model, struct PSDposterior *psd);

void populate_psd_histogram(struct Data *data, struct Spacecraft *lpf, struct Model *model, int MCMCSTEPS, struct PSDposterior *psd);

void print_power_spectra(char filename[], double *d, double *h, int N, double *Snf, double Tobs, int imin, int imax);

void print_time_domain_waveforms(char filename[], double *h, int N, double *Snf, double eta, double Tobs, int imin, int imax, double tmin, double tmax);

/* ********************************************************************************** */
/*                                                                                    */
/*                                    Math tools                                      */
/*                                                                                    */
/* ********************************************************************************** */

double snr(struct Data *data, struct Spacecraft *lpf, struct Model *model);

double fourier_nwip(int imin, int imax, double *a, double *b, double *Sn);

void crossproduct(double *b, double *c, double *a);

void matrix_multiply(double **A, double **B, double **C, int N);

void matrix_invert(double **A, double **invA, int N);

void check_incidence(struct Spacecraft *lpf,struct Model * model);

void dfour1(double data[], unsigned long nn, int isign);

void drealft(double data[], unsigned long n, int isign);

/* ********************************************************************************** */
/*                                                                                    */
/*                        Data handling and injection routines                        */
/*                                                                                    */
/* ********************************************************************************** */

void simulate_data(struct Data *data);

void simulate_injection(struct Data *data, struct Spacecraft *lpf, struct Model *injection);

void simulate_noise(struct Data *data, struct Spacecraft *lpf, struct Model *injection, gsl_rng *r);

/* ********************************************************************************** */
/*                                                                                    */
/*                           Memory (de)allocation routines                           */
/*                                                                                    */
/* ********************************************************************************** */

void copy_source(struct Source *source, struct Source *copy, int DOF);

void copy_model(struct Model *model, struct Model *copy, int N, int DOF);

void initialize_source(struct Source *source);

void initialize_model(struct Model *model, int N, int D, int DOF);

void initialize_bayesline(struct BayesLineParams **bayesline, struct Data *data, double **psd);

void free_source(struct Source *source);

