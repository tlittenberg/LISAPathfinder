/*
 Compile:
 gcc -O3 -lm -lgsl -o pvalue_quantiles pvalue_quantiles.c
 
 Run:
 ./pvalue_quantiles
 
 Plot:
 gnuplot pvalue_quantiles.gpi
 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>

double draw_from_pdf(gsl_rng *seed, double alpha)
{
  /*
   Use the cdf to draw random numbers
   from the desired pdf:
   
   cdf(x) = int_{-infinity}^x pdf(x') dx'
   
   U[0,1] = cdf(x)
   
   Solve for x
   */
  
  //p(x) = U[0,1]
  //return gsl_rng_uniform(seed);
  
  //p(x) = -x^(-2) for x > 1
  //double alpha = 2.1;
  return pow( 1. - (gsl_rng_uniform(seed)), 1./(1.-alpha) );
}

int main( int argc, char* argv[] )
{
  int N = 42;       //Number of samples
  int M = 10000;  //Number of Monte Carlo realizations
  int D = 1000;     //Number of chain samples

  //set up & initialize random number generator
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *seed = gsl_rng_alloc(T);
  gsl_rng_env_setup();
  gsl_rng_set (seed, 150914);
  
  double *pvalues  = malloc(M*D*sizeof(double));
  double **samples = malloc(M*D*sizeof(double *));
  for(int m=0; m<M*D; m++) samples[m] = malloc(N*sizeof(double));
  
  FILE *alphaFile = fopen("alpha.dat","r");
  double *alpha = malloc(D*sizeof(double));
  for(int n=0; n<D; n++) fscanf(alphaFile,"%lg",&alpha[n]);
  fclose(alphaFile);
  
  // Monte Carlo M realizations of N samples from the cdf
  for(int d=0; d<D; d++)
  {
    if(d%(D/100)==0)
    {
      printf ("\r Monte Carlo simulations %.0f finished",100.*(double)(d+1)/(double)D);
      fflush(stdout);
    }
    for(int m=0; m<M; m++)
    {
      //let me know you're working
      //      if(m%(M/100)==0)printf ("\r Monte Carlo simulations %.0f finished",100.*(double)(m+1)/(double)M);
      //      fflush(stdout);
      
      //get random realization of desired distribution
      for(int n=0; n<N; n++)
      {
        samples[d*M+m][n] = draw_from_pdf(seed,alpha[d]);
      }
      //sort samples to make cdf of each Monte Carlo realization
      gsl_sort(samples[d*M+m],1,N);
    }
  }
  printf("\n");


  FILE *outfile = fopen("pvalue_quantiles.dat","w");
  
  //step through every p-value n/N)
  for(int n=0; n<N; n++)
  {
    //let me know you're working
    printf ("\r Computing quantiles %.0f finished",100.*(double)(n+1)/(double)N);
    fflush(stdout);

    //get all Monte carlo results at a partiuclar p-value
    for(int m=0; m<D*M; m++) pvalues[m] = samples[m][n];
    
    //sort all the values at particular p-value to get quantiles
    gsl_sort(pvalues,1,D*M);
    
    // "3 sigma" quantiles
    double L3 = gsl_stats_quantile_from_sorted_data (pvalues, 1, D*M, .0015);
    double U3 = gsl_stats_quantile_from_sorted_data (pvalues, 1, D*M, .9985);

    // "2 sigma" quantiles
    double L2 = gsl_stats_quantile_from_sorted_data (pvalues, 1, D*M, .0250);
    double U2 = gsl_stats_quantile_from_sorted_data (pvalues, 1, D*M, .9750);
    
    // "1 sigma" quantiles
    double L1 = gsl_stats_quantile_from_sorted_data (pvalues, 1, D*M, .1590);
    double U1 = gsl_stats_quantile_from_sorted_data (pvalues, 1, D*M, .8410);

    //print quantiles to file
    fprintf(outfile,"%.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",(double)(n+1)/N, L3, U3, L2, U2, L1, U1);
    
  }
  printf("\n");
  
  //clean up after yourself
  fclose(outfile);
  for(int m=0; m<D*M; m++) free(samples[m]);
  free(samples);
  free(pvalues);
  
  return 0;
}
