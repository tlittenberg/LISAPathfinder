/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <getopt.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "LPF.h"
#include "Subroutines.h"
#include "LISAPathfinder.h"
#include "TimePhaseMaximization.h"

/* ============================  MAIN PROGRAM  ============================ */

void parse(int argc, char **argv, struct Data *data, struct Flags *flags);


int main(int argc, char **argv)
{

  /* declare variables */
  int i,ic,n,mc;
  int accept;
  int MCMCSTEPS;
  int BURNIN;
  int NC;

  double H;
  double alpha;


   /* Structure for run flags */
   struct Flags *flags = malloc(sizeof(struct Flags));

  /* Initialize spacecraft structure */
  struct Spacecraft *lpf = malloc(sizeof(struct Spacecraft));

  initialize_spacecraft(lpf);

  // Spacecraft mass
  lpf->M = EOM_SC_M;

  // Housing 1 Geometry
  lpf->R[0][0] = EOM_H1SC_X;
  lpf->R[0][1] = EOM_H1SC_Y;
  lpf->R[0][2] = EOM_H1SC_Z;
  // Housing 2 Geometry
  lpf->R[1][0] = EOM_H2SC_X;
  lpf->R[1][1] = EOM_H2SC_Y;
  lpf->R[1][2] = EOM_H2SC_Z;

  // Spacecraft Geometry
  lpf->x[0][0] = SC_BOT_CORNER_6_X;
  lpf->x[0][1] = SC_BOT_CORNER_6_Y;
  lpf->x[1][0] = SC_BOT_CORNER_5_X;
  lpf->x[1][1] = SC_BOT_CORNER_5_Y;
  lpf->x[2][0] = SC_BOT_CORNER_4_X;
  lpf->x[2][1] = SC_BOT_CORNER_4_Y;
  lpf->x[3][0] = SC_BOT_CORNER_3_X;
  lpf->x[3][1] = SC_BOT_CORNER_3_Y;
  lpf->x[4][0] = SC_BOT_CORNER_2_X;
  lpf->x[4][1] = SC_BOT_CORNER_2_Y;
  lpf->x[5][0] = SC_BOT_CORNER_1_X;
  lpf->x[5][1] = SC_BOT_CORNER_1_Y;
  lpf->x[6][0] = SC_BOT_CORNER_8_X;
  lpf->x[6][1] = SC_BOT_CORNER_8_Y;
  lpf->x[7][0] = SC_BOT_CORNER_7_X;
  lpf->x[7][1] = SC_BOT_CORNER_7_Y;
  lpf->x[8][0] = SC_BOT_CORNER_6_X;
  lpf->x[8][1] = SC_BOT_CORNER_6_Y;

  // Spacecraft dimensions
  lpf->H = SC_H;
  lpf->W = lpf->x[0][0] - lpf->x[4][0];
  lpf->D = lpf->x[2][1] - lpf->x[6][1];

  // Moment of inertia about H1
  lpf->I[0][0][0] = EOM_SC_IH1_XX;
  lpf->I[0][0][1] = EOM_SC_IH1_XY;
  lpf->I[0][0][2] = EOM_SC_IH1_XZ;

  lpf->I[0][1][0] = EOM_SC_IH1_YX;
  lpf->I[0][1][1] = EOM_SC_IH1_YY;
  lpf->I[0][1][2] = EOM_SC_IH1_YZ;

  lpf->I[0][2][0] = EOM_SC_IH1_ZX;
  lpf->I[0][2][1] = EOM_SC_IH1_ZY;
  lpf->I[0][2][2] = EOM_SC_IH1_ZZ;

  // Moment of inertia about H2
  lpf->I[1][0][0] = EOM_SC_IH2_XX;
  lpf->I[1][0][1] = EOM_SC_IH2_XY;
  lpf->I[1][0][2] = EOM_SC_IH2_XZ;

  lpf->I[1][1][0] = EOM_SC_IH2_YX;
  lpf->I[1][1][1] = EOM_SC_IH2_YY;
  lpf->I[1][1][2] = EOM_SC_IH2_YZ;

  lpf->I[1][2][0] = EOM_SC_IH2_ZX;
  lpf->I[1][2][1] = EOM_SC_IH2_ZY;
  lpf->I[1][2][2] = EOM_SC_IH2_ZZ;

  matrix_invert(lpf->I[0], lpf->invI[0], 3);
  matrix_invert(lpf->I[1], lpf->invI[1], 3);


  /* Initialize data structure */
  struct Data  *data = malloc(sizeof(struct Data));

  data->T  = 1024.;
  data->dt = 0.5;
  data->df = 1.0/data->T;
  data->N  = (int)(data->T/data->dt)/2;

   parse(argc, argv, data, flags);

  data->fmin = 1.0e-4; //Hz
  data->fmax = (double)data->N/data->T;  //Hz

  data->imin = (int)floor(data->fmin*data->T);
  data->imax = (int)floor(data->fmax*data->T);

  data->d = malloc(data->DOF*sizeof(double *));
  data->s = malloc(data->DOF*sizeof(double *));
  data->n = malloc(data->DOF*sizeof(double *));

  for(i=0; i<data->DOF; i++)
  {
    data->d[i] = malloc(data->N*2*sizeof(double));
    data->s[i] = malloc(data->N*2*sizeof(double));
    data->n[i] = malloc(data->N*2*sizeof(double));
  }

  data->f = malloc(data->N*sizeof(double));

   /* set up GSL random number generator */
   const gsl_rng_type *T = gsl_rng_default;
   gsl_rng *nr = gsl_rng_alloc (T);
   gsl_rng *r  = gsl_rng_alloc (T);
   gsl_rng_env_setup();
   gsl_rng_set (r, data->seed);
   gsl_rng_set (nr, data->nseed);

  /* Simulate noise data */
  struct Model *injection = malloc(sizeof(struct Model));
  initialize_model(injection,data->N,6, data->DOF);


  for(i=0; i<3; i++)
  {
    injection->Ais[i] = 2.0e-9; // m*Hz^-1/2
    injection->Ath[i] = 1.0e-8; // N*Hz^-1/2
    injection->Ars[i] = 2.0e-7; // rad*Hz^-1/2
  }

  /* Simulate source data */
   struct Source *source;

    injection->N = 1;
  for(n=0; n<injection->N; n++)
  {
    source = injection->source[n];
    source->face = -1;
    while(source->face ==-1) draw_impact_point(data, lpf, source, r);
      source->P = 20;
    printf("hit on face %i\n",source->face);
  }

  simulate_noise(data, lpf, injection, nr);
  simulate_data(data);
  simulate_injection(data,lpf,injection);
  injection->logL = loglikelihood(data, lpf, injection, flags);

  printf("Injected parameters:   \n");
  FILE *injfile = fopen("injection.dat","w");
  for(n=0; n<injection->N; n++)printf("     {%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg}\n"              ,injection->source[n]->t0, injection->source[n]->P, injection->source[n]->costheta, injection->source[n]->phi, injection->source[n]->map[0], injection->source[n]->map[1], injection->source[n]->r[0], injection->source[n]->r[1], injection->source[n]->r[2]);
  for(n=0; n<injection->N; n++)fprintf(injfile,"%lg %lg %lg %lg %lg %lg %lg %lg %lg\n",injection->source[n]->t0, injection->source[n]->P, injection->source[n]->costheta, injection->source[n]->phi, injection->source[n]->map[0], injection->source[n]->map[1], injection->source[n]->r[0], injection->source[n]->r[1], injection->source[n]->r[2]);

  fclose(injfile);
  printf("SNR of injection = %g\n",snr(data, lpf, injection));

  /* Initialize parallel chains */
  NC = 15;
  int *index = malloc(NC*sizeof(double));
  double *temp = malloc(NC*sizeof(double));
  double dT = 1.5;
  temp[0] = 1.0;
  index[0]=0;
  for(ic=1; ic<NC; ic++)
  {
    temp[ic]=temp[ic-1]*dT;
    index[ic]=ic;
  }
  temp[NC-1]=1.e6;

  /* Initialize model */
  struct Model **model = malloc(NC*sizeof(struct Model*));

  struct Model *trial  = malloc(sizeof(struct Model));
  initialize_model(trial, data->N, 10, data->DOF);


  for(ic=0; ic<NC; ic++)
  {
    model[ic] = malloc(sizeof(struct Model));
    initialize_model(model[ic],data->N,10, data->DOF);

    copy_model(injection, model[ic], data->N, data->DOF);


    detector_proposal(data,injection,model[ic],r);

    for(n=0; n<model[ic]->N; n++)
    {
      model[ic]->source[n]->P  = gsl_rng_uniform(r)*100;
      model[ic]->source[n]->t0 = gsl_rng_uniform(r)*data->T;

    }

    logprior(data, model[ic], injection);
    model[ic]->logL = loglikelihood(data, lpf, model[ic], flags);
  }

  /* set up distribution */
  //struct PSDposterior *psd = NULL;
  //setup_psd_histogram(data, injection, psd);

  /* set up MCMC run */
  accept    = 0;
  MCMCSTEPS = 1000000;
  BURNIN    = MCMCSTEPS/100;//1000;

  char filename[128];

  FILE *noisechain;
  sprintf(filename,"noisechain.dat");
  noisechain = fopen(filename,"w");
  //fprintf(noisechain,"#dlogL mass dAi[x] dAth[x] dAi[y] dAth[y] dAi[z] dAth[z]\n");

  FILE *impactchain;
  sprintf(filename,"impactchain.dat");
  impactchain = fopen(filename,"w");
  //fprintf(impactchain,"#dlogL N t0[0] P[0] costheta[0] phi[0] ... \n");

  FILE *logLchain;
  sprintf(filename,"logLchain.dat");
  logLchain = fopen(filename,"w");
  //fprintf(logLchain,"#dlogL[0] dlogL[1] ... T[0] T[1]...\n");

  int reject;

  /* Here is the MCMC loop */
  for(mc=0;mc<MCMCSTEPS;mc++)
  {

    for(ic=0; ic<NC; ic++)
    {
      for(n=0; n<10; n++)
      {
        reject=0;

        //copy x to y
        copy_model(model[index[ic]], trial, data->N, data->DOF);

        //choose new parameters for y
        proposal(flags, data, lpf, model[index[ic]], trial, r, &reject);
//        if(ic==0 && trial->source[0]->face==1)printf("draw={%g,%g,%g}\n", trial->source[0]->r[0],trial->source[0]->r[1],trial->source[0]->r[2]);



        //compute maximized likelihood
        //if(mc<BURNIN) max_loglikelihood(data, trial);
        if(reject) continue;
        else
        {
          //compute new likelihood
          trial->logL = loglikelihood(data, lpf, trial, flags);

          //compute new prior
          logprior(data, trial, injection);

          //compute Hastings ratio
          H     = (trial->logL - model[index[ic]]->logL)/temp[index[ic]] + trial->logP - model[index[ic]]->logP;
          alpha = log(gsl_rng_uniform(r));

//          if(ic==0 && trial->source[0]->face==1)printf("H=%g, logLy=%g, logLx=%g\n",H,trial->logL,model[index[ic]]->logL);

          //adopt new position w/ probability H
          if(H>alpha)
          {
            copy_model(trial, model[index[ic]], data->N, data->DOF);
            accept++;
          }
          source=trial->source[0];
//          if(trial->logL>-1e60)fprintf(stdout,"%lg %lg %lg %lg %lg %lg %i %lg %lg %lg\n", trial->logL-injection->logL ,source->P,source->map[0], source->map[1], source->costheta,source->phi,source->face, source->r[0], source->r[1], source->r[2]);

        }//Metropolis-Hastings
      }//Loop over inter-model updates
    }//Loop over chains


    ptmcmc(model, temp, index, r, NC, mc);

    //cute PSD histogram
    //if(mc>MCMCSTEPS/2) populate_psd_histogram(data, model[index[0]], MCMCSTEPS, psd);

    //print chain files

    //impact parameters
    ic = index[0];
    fprintf(noisechain,"%lg ",model[ic]->logL-injection->logL);
    //fprintf(noisechain,"%lg ",model[ic]->mass);
    for(i=0; i<3; i++)
    {
      fprintf(noisechain,"%lg %lg %lg ",(injection->Ais[i]-model[ic]->Ais[i])/injection->Ais[i],(injection->Ath[i]-model[ic]->Ath[i])/injection->Ath[i],(injection->Ars[i]-model[ic]->Ars[i])/injection->Ars[i]);
    }
    fprintf(noisechain,"\n");
    
    for(n=0; n<model[ic]->N; n++)
    {
      source = model[ic]->source[n];
      face2map(lpf, source->r,source->map);
      which_face_r(source->r);
      fprintf(impactchain,"%lg ",model[ic]->logL-injection->logL);
      fprintf(impactchain,"%i ",model[ic]->N);
      fprintf(impactchain,"%lg %lg %lg %lg %lg %lg %i %lg %lg %lg\n", source->t0,source->P,source->map[0], source->map[1], source->costheta,source->phi,source->face, source->r[0], source->r[1], source->r[2]);
    }fflush(impactchain);
    
    if(flags->verbose)
    {
      for(n=0; n<model[ic]->N; n++)
      {
        source = model[ic]->source[n];
        fprintf(stdout,"%lg ",model[ic]->logL-injection->logL);
        fprintf(stdout,"%i ",model[ic]->N);
        fprintf(stdout,"%lg %lg %lg %lg %lg %lg %i %lg %lg %lg\n", source->t0,source->P,source->map[0], source->map[1], source->costheta,source->phi,source->face, source->r[0], source->r[1], source->r[2]);
      }
    }

    //parallel tempering
    for(ic=0; ic<NC; ic++) fprintf(logLchain,"%lg ",model[index[ic]]->logL-injection->logL);
    for(ic=0; ic<NC; ic++) fprintf(logLchain,"%lg ",temp[ic]);
    fprintf(logLchain,"\n");





  }//MCMC
  printf("acceptance rate = %g\n",(double)accept/(double)MCMCSTEPS);


//  FILE *PSDfile = fopen("psdhistogram.dat","w");
//  for(i=0; i<psd->Nx; i++)
//  {
//    for(j=0; j<psd->Ny; j++)
//    {
//      fprintf(PSDfile,"%lg %lg %lg\n",psd->xmin + i*psd->dx, psd->ymin + j*psd->dy, psd->histogram[i][j]);
//    }
//  }
//  fclose(PSDfile);

  return 0;
}



static void print_usage() {
   printf("\n");
   printf("Usage: \n");
   printf("REQUIRED:\n");
   printf("  -d | --dof     : degrees of freedom (3 or 6) \n");
   printf("  -n | --nseed   : seed for noise simulation \n");
   printf("  -s | --seed    : seed for MCMC \n");
   printf("OPTIONAL:\n");
   printf("  -f | --fixd    : fixed dimension (no RJ)  \n");
   printf("  -h | --help    : usage information        \n");
   printf("  -p | --prior   : sample prior             \n");
   printf("  -v | --verbose : enable verbose output    \n");
   printf("EXAMPLE:\n");
   printf("./mcmc --dof 6 --seed 1234 --nseed 1234\n");
   printf("\n");
   exit(EXIT_FAILURE);
}



void parse(int argc, char **argv, struct Data *data, struct Flags *flags)
{
   data->DOF = 6;
   data->seed  = 1234;
   data->nseed = 1234;
   flags->verbose = 0;
   flags->prior = 0;
   flags->rj = 1;

   if(argc==1) print_usage();

   //Specifying the expected options
   static struct option long_options[] = {
      {"dof",     required_argument, 0,  'd' },
      {"fixd",    no_argument,       0,  'f' },
      {"help",    no_argument,       0,  'h' },
      {"nseed",   required_argument, 0,  'n' },
      {"seed",    required_argument, 0,  's' },
      {"verbose", no_argument,       0,  'v' },
      {"prior",   no_argument,       0,  'p' },
      {0,         0,                 0,   0  }
   };

   int opt=0;
   int long_index =0;

   //Loop through argv string and pluck out arguments
   while ((opt = getopt_long_only(argc, argv,"apl:b:",
                             long_options, &long_index )) != -1) {
      switch (opt) {
         case 'd' : data->DOF = atoi(optarg);
            break;
        case 'f' : flags->rj = 0;
          break;
         case 'h' :
            print_usage();
            exit(EXIT_FAILURE);
            break;
         case 'p' : flags->prior = 1;
              break;
        case 'n' : data->nseed = atoi(optarg);
          break;
         case 's' : data->seed = atoi(optarg);
            break;
         case 'v' : flags->verbose = 1;
            break;
         default: print_usage();
            exit(EXIT_FAILURE);
      }
   }

   //Report on set parameters
   fprintf(stdout,"***** RUN SETTINGS *****\n");
   fprintf(stdout,"Number of data channles ... %i\n",data->DOF);
   fprintf(stdout,"Random number seed ........ %li\n",data->seed);
   fprintf(stdout,"******* RUN FLAGS ******\n");
   if(flags->verbose)fprintf(stdout,"Verbose flag .............. ENABLED \n");
   else              fprintf(stdout,"Verbose flag .............. DISABLED\n");
   fprintf(stdout,"\n");

}
















