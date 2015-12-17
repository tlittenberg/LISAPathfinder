//
//  test.c
//  
//
//  Created by Tyson Littenberg on 12/2/15.
//
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Subroutines.h"
#include "LISAPathfinder.h"
#include "TimePhaseMaximization.h"

int main()
{

  int i;

  /* set up GSL random number generator */
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc (T);
  gsl_rng_env_setup();

  struct Source *source = malloc(sizeof(struct Source));
  initialize_source(source);

  /* Initialize data structure */
  struct Data  *data = malloc(sizeof(struct Data));

  data->T  = 1024.;
  data->dt = 0.5;
  data->df = 1.0/data->T;
  data->N  = (int)(data->T/data->dt)/2;

  data->fmin = 1.0e-4; //Hz
  data->fmax = (double)data->N/data->T;  //Hz

  data->imin = (int)floor(data->fmin*data->T);
  data->imax = (int)floor(data->fmax*data->T);

  data->d = malloc(3*sizeof(double *));
  data->s = malloc(3*sizeof(double *));
  data->n = malloc(3*sizeof(double *));

  for(i=0; i<3; i++)
  {
    data->d[i] = malloc(data->N*2*sizeof(double));
    data->s[i] = malloc(data->N*2*sizeof(double));
    data->n[i] = malloc(data->N*2*sizeof(double));
  }

  data->f = malloc(data->N*sizeof(double));

  for(i=0; i<1000000; i++)
  {
    draw_impact_point(data,source,r);


    if(source->face==4 || source->face==5)
    //if(source->face>=0)
    {
      printf("%lg %lg ",source->map[0],source->map[1]);
      printf("%lg %lg %lg ",source->r[0],source->r[1],source->r[2]);
      printf("%lg %lg ",source->costheta,source->phi);
      printf("%i\n",source->face);

    }
  }
  
  return 0;
}