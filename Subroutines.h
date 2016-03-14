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
  double logL; //likelihod
  double logP; //prior

  /* Impact parameters */

  int N;       //number of impacts

  struct Source **source;

  /* Instrument response */
  double **s;
};

struct Spacecraft
{
  /* Spacecraft parameters */
  double **R;        //position of proof masses
  double **x;        //position of spacecraft corners
  double ***I;       //moment of inertia (ultimately a tensor)
  double ***invI;    //inverse moment of inertia (ultimately a tensor)
  double M;          //mass
  double W;          //spacecraft width  (x-direction, m)
  double D;          //spacecraft depth  (y-direction, m)
  double H;          //spacecraft height (z-direciton, m)

};

struct Data
{
  int N;
   int DOF;
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
  double **d;
  double **n;
  double **s;
  double *f;
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
};


/* ********************************************************************************** */
/*                                                                                    */
/*                                    MCMC tools                                      */
/*                                                                                    */
/* ********************************************************************************** */

void ptmcmc(struct Model **model, double *temp, int *index, gsl_rng *r, int NC, int mc);

void proposal(struct Flags *flags, struct Data *data, struct Spacecraft *lpf, struct Model *model, struct Model *trial, gsl_rng *r, int *reject);

void dimension_proposal(struct Data *data, struct Spacecraft *lpf, struct Model *model, struct Model *trial, gsl_rng *r, int Nmax, int *test);

void detector_proposal(struct Data *data, struct Model *model, struct Model *trial, gsl_rng *r);

void impact_proposal(struct Data *data, struct Spacecraft *lpf, struct Source *model, struct Source *trial, gsl_rng *r);

void draw_impact_point(struct Data *data, struct Spacecraft *lpf, struct Source *source, gsl_rng *seed);

void logprior(struct Data *data, struct Model *model, struct Model *injection);

double log_mass_prior(double m0, double m);

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

void free_source(struct Source *source);

void initialize_spacecraft(struct Spacecraft *spacecraft);

