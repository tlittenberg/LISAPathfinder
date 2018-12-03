/* Compile
 gcc power_law_mcmc.c -lm -lgsl -o power_law_mcmc
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>


#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
#define MC 10000000
#define DS 1000

static void printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}

struct Data
{
  int size;
  double *mean;
  double *variance;
  double *sigma;
  FILE *file;
};

struct Model
{
  int size;
  double alpha;
  double beta;
  double xmin;
  double xknee;
  double logL;
};

void loglikelihood(struct Data *data, struct Model *model, gsl_rng *seed);
double powerlaw_pdf(double alpha, double x, double xmin);
double broken_powerlaw_pdf(double alpha, double beta, double x, double xmin, double xknee);
void copy_model(struct Model *origin, struct Model *copy);

int main (int argc, char **argv)
{
  /*
   Parse command line
   */
  if(argc!=3)
  {
    printf("Usage: \n ./power_law_mcmc [data file] [D]\n");
    exit(0);
  }
  
  /*
   Read data
   */
  
  //set up data structure
  struct Data *data = malloc(sizeof(struct Data));
  data->file = fopen(argv[1],"r");
  
  //how big is the file?
  double junk;
  data->size = 0;
  while(!feof(data->file))
  {
    fscanf(data->file,"%lg%lg",&junk,&junk);
    data->size++;
  }
  data->size--;
  rewind(data->file);
  
  //store data
  data->mean     = malloc(data->size*sizeof(double));
  data->variance = malloc(data->size*sizeof(double));
  data->sigma    = malloc(data->size*sizeof(double));
  for(int n=0; n<data->size; n++)
  {
    fscanf(data->file,"%lg%lg",&data->mean[n],&data->variance[n]);
    data->sigma[n] = sqrt(data->variance[n]);
  }
  
  
  /*
   Set up MCMC
   */

  //RNG
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *seed = gsl_rng_alloc(T);
  gsl_rng_env_setup();
  gsl_rng_set (seed, time(0));

  //power law model parameters
  struct Model *model_x = malloc(sizeof(struct Model));
  struct Model *model_y = malloc(sizeof(struct Model));

  //initialize model likelihood
  model_x->size  = atoi(argv[2]);
  model_x->alpha = 1.5;
  model_x->beta  = 2.0;
  model_x->xmin  = 1.0;
  model_x->xknee = 10.0;
  loglikelihood(data, model_x, seed);
  copy_model(model_x,model_y);
  

  /*
   The MCMC Loop
   */
  
  //Chain productes
  FILE *chain = fopen("chain.dat","w");
  double logLmax=-INFINITY;
  double **fit = malloc(1000*sizeof(double *));
  for(int i=0; i<1000; i++) fit[i] = malloc((MC/DS)*sizeof(double *));
  

  //Here goes...
  for(int mc=0; mc<MC; mc++)
  {
    
    //proposal
    model_y->alpha = model_x->alpha + gsl_ran_ugaussian(seed)*0.2;
    model_y->beta  = model_y->alpha;
    model_y->xknee = model_x->xknee;
    if(model_x->size>1)
      model_y->beta  = model_x->beta  + gsl_ran_ugaussian(seed)*1.0;
    if(model_x->size>2)
      model_y->xknee = model_x->xknee + gsl_ran_ugaussian(seed)*5.0;

    //prior
    double logP = 0.0;
    if(model_y->beta<1.0) logP=-INFINITY;
    if(model_y->xknee<5 || model_y->xknee>100) logP=-INFINITY;
      
    //likelihood
    loglikelihood(data, model_y, seed);

    //Hasting's ratio
    double logH  = model_y->logL - model_x->logL + logP;

//    printf("logH=%g, logP=%g\n",logH,logP);
//    printf("   y: logL=%g,alpha=%g,beta=%g,knee=%g\n",model_y->logL,model_y->alpha,model_y->beta, model_y->xknee);
//    printf("   x: logL=%g,alpha=%g,beta=%g,knee=%g\n",model_x->logL,model_x->alpha,model_x->beta, model_x->xknee);

    if(logH > log(gsl_rng_uniform(seed)) && model_y->alpha)
    {
      copy_model(model_y,model_x);
      if(model_x->logL>logLmax) logLmax = model_x->logL;
    }
    
    //output chain file
    if(mc%DS==0)
    {
      fprintf(chain,"%lg %lg %lg %lg\n",model_x->logL,model_x->alpha,model_x->beta, model_x->xknee);
      for(int i=0; i<1000; i++)
        fit[i][mc/DS] = broken_powerlaw_pdf(model_x->alpha, model_x->beta,(double)i, model_x->xmin,model_x->xknee);
    }
    
    printProgress((double)mc/(double)(MC));
  }
  printProgress(1.0);
  printf("\n");

  fclose(chain);
  
  /*
   Bayesian Information Criterion
   */
  FILE *bicFile = fopen("bic.dat","w");
  fprintf(stdout,"BIC for model %i: %g\n",model_x->size,log(data->size)*(double)model_x->size - 2.*logLmax);
  fprintf(bicFile,"%i %g\n",model_x->size,log(data->size)*(double)model_x->size - 2.*logLmax);
  fclose(bicFile);
  
  /*
   Get credible intervals of model
   */
  FILE *fitFile = fopen("fit.dat","w");
  for(int i=0; i<1000; i++)
  {
    printProgress((double)i/(double)(1000));

    gsl_sort(fit[i],1,MC/DS);
    
    fprintf(fitFile,"%lg ",(double)i);
    fprintf(fitFile,"%lg ",gsl_stats_quantile_from_sorted_data (fit[i], 1, MC/DS, 0.05));
    fprintf(fitFile,"%lg ",gsl_stats_quantile_from_sorted_data (fit[i], 1, MC/DS, 0.25));
    fprintf(fitFile,"%lg ",gsl_stats_quantile_from_sorted_data (fit[i], 1, MC/DS, 0.50));
    fprintf(fitFile,"%lg ",gsl_stats_quantile_from_sorted_data (fit[i], 1, MC/DS, 0.75));
    fprintf(fitFile,"%lg ",gsl_stats_quantile_from_sorted_data (fit[i], 1, MC/DS, 0.95));
    fprintf(fitFile,"\n");
  }
  fclose(fitFile);
  printProgress(1.0);
  printf("\n");


  
  return 0;
}


void loglikelihood(struct Data *data, struct Model *model, gsl_rng *seed)
{
  double p;
  model->logL = 0.0;
  
  if(model->beta<1.0)
    model->logL = -INFINITY;
  else{
    
    for(int n=0; n<data->size; n++)
    {
      //fudge convolution over momentum errors
      p=0.0;
      
      for(int i=0; i<100; i++)
      {
        double x = data->mean[n] + gsl_ran_ugaussian(seed)*data->sigma[n];
        if(x>model->xmin)
        {
          //p = powerlaw_pdf(model->alpha, data->mean[n], model->xmin)*exp(-model->beta*data->mean[n]);
          //p = powerlaw_pdf(model->alpha, data->mean[n], model->xmin);
          p += broken_powerlaw_pdf(model->alpha, model->beta,x, model->xmin,model->xknee)/10;
          //model->logL += -0.5*p*p/data->variance[n];
          
        }
      }
      if(p!=0.)model->logL += log(p);
    }
  }
  
}

double powerlaw_pdf(double alpha, double x, double xmin)
{
  return ((alpha-1)/xmin)*pow((x/xmin),-alpha);
}

double broken_powerlaw_pdf(double alpha, double beta, double x, double xmin, double xknee)
{
  double knee = xknee;//pow(10,xknee);
  double norm = ((1./(1.-alpha))*(pow(knee,1.-alpha) - pow(xmin,1.-alpha))) - ((1./(1.-beta))*(pow(knee,1.-alpha)) );
  double A = 1./norm;
  double B = A*pow(knee,beta-alpha);
  
  if(x<=knee)
    return A*pow(x,-alpha);
  else
    return B*pow(x,-beta);
}

void copy_model(struct Model *origin, struct Model *copy)
{
  copy->alpha = origin->alpha;
  copy->beta  = origin->beta;
  copy->xmin  = origin->xmin;
  copy->xknee = origin->xknee;
  copy->logL  = origin->logL;
}
