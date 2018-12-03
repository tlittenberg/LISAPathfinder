/*
 * gcc histograms.c -lm -o histograms
 *   -makes marginalized pdfs from MCMC chains
 *	-ignores first two columns as iteration and log-likelihood
 */

/***************************  REQUIRED LIBRARIES  ***************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>



#define NR_END 1
#define FREE_ARG char*

/*************  PROTOTYPE DECLARATIONS FOR INTERNAL FUNCTIONS  **************/
double gaussian(double x, double mean, double sigma);

/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char* argv[])
{
	/* -------------  DECLARATION OF VARIABLES  ------------- */
	int i, j, k, n;
	int count, bin;
	int Nparam;
	int NB;

	int col=0;
	double cmin=0;
	double cmax=0;

	double norm;
	double max;
	double bw;
	double x;
  int oversample=10;

	char filename[100];

	FILE* infile;
	FILE* outfile;


	if(argc<4) { printf("./histograms chainfile Nparam burn-in [downsample] [Nbin]\n"); return 0;}

	Nparam	= atoi(argv[2]);
	NB		= atoi(argv[3]);

  int DS=1;
  if(argc>=5) DS=atoi(argv[4]);

	int Nbin=101;
	if(argc==6)  Nbin = atoi(argv[5]);

  /*
	int Nhist[Nparam];
	double Parameters[Nparam];
	double Psigma[Nparam];
	double Pwidth[Nparam];
	double Pmean[Nparam];
	double Pmode[Nparam];
	double Pmax[Nparam];
	double Pmin[Nparam];
	double Pmle[Nparam];
	double KDEsigma[Nparam];
  double KDEwidth[Nparam];
   */
	int *Nhist=(int *)malloc(sizeof(int)*Nparam);//;
	double *Parameters=(double *)malloc(sizeof(double)*Nparam);//[Nparam];
	double *Psigma=(double *)malloc(sizeof(double)*Nparam);//;
	double *Pwidth=(double *)malloc(sizeof(double)*Nparam);//;
	double *Pmean=(double *)malloc(sizeof(double)*Nparam);//;
	double *Pmode=(double *)malloc(sizeof(double)*Nparam);//;
	double *Pmax=(double *)malloc(sizeof(double)*Nparam);//;
	double *Pmin=(double *)malloc(sizeof(double)*Nparam);//;
	double *Pmle=(double *)malloc(sizeof(double)*Nparam);//;
	double *KDEsigma=(double *)malloc(sizeof(double)*Nparam);//;
  double *KDEwidth=(double *)malloc(sizeof(double)*Nparam);//;
	double *KDEmin=(double *)malloc(sizeof(double)*Nparam);//;
  double *KDEmax=(double *)malloc(sizeof(double)*Nparam);//;

	/* initialize arrays */
	for(i=0; i<Nparam; i++)
	{
		Pmean[i] = 0.0;
		Pmin[i] = 1.0e60;
		Pmax[i] =-1.0e60;
	}


	/*  open chain file and find max and min for each parameter  */
	sprintf(filename,"%s",argv[1]);
	infile = fopen(filename,"r");
	char burnrows[10000];
	n = count = 0;
	while(!feof(infile))
	{
    if(n<NB)fgets(burnrows,10000,infile);
    else
    {
			for(i=0; i<Nparam; i++)
			{
        fscanf(infile,"%lg",  &Parameters[i]);
        if(n%DS==0)
        {
          Pmean[i] += Parameters[i];
          if(Parameters[i] < Pmin[i])  Pmin[i] = Parameters[i];
          if(Parameters[i] > Pmax[i])  Pmax[i] = Parameters[i];
        }
			}
      if(n%DS==0)count++;
    }
    n++;
	}
  int NSAMP = n-1-NB;
	rewind(infile);

  double **cdf=(double **)malloc(sizeof(double*)*Nparam);
  for(i=0; i<Nparam; i++) cdf[i] = (double *)malloc(sizeof(double)*NSAMP);

	for(i=0; i<Nparam; i++)
	{
		Nhist[i] = Nbin;
		Pwidth[i] = (Pmax[i]-Pmin[i])/(double)(Nhist[i] - 1);
    KDEwidth[i]=Pwidth[i]/(double)oversample;
	}

	int Nmax=0;
	for(i=0; i<Nparam; i++)
	{
		if(Nhist[i] > Nmax) Nmax = Nhist[i]+1;
	}

	//double histograms[Nparam][Nmax];
  double **histograms = (double **)malloc(sizeof(double *)*Nparam);
  for(i=0; i<Nparam; i++) histograms[i]=(double *)malloc(sizeof(double)*Nmax);
	for(i=0; i<Nparam; i++) for(j=0; j<Nmax; j++) histograms[i][j] = 0.0;

	for(i=0; i<Nparam; i++)
	{
		Pmean[i] /= (double)(count);
		Psigma[i] = 0.0;
	}

	/* calculate varience */
	n = 0;
	while(!feof(infile))
	{
    if(n<NB)fgets(burnrows,10000,infile);
    else
    {
      for(i=0; i<Nparam; i++)
      {
        fscanf(infile,"%lg",  &Parameters[i]);
        if(n%DS==0)Psigma[i] += (Parameters[i]-Pmean[i])*(Parameters[i]-Pmean[i]);
//        printf("%i %lg %lg\n",n-NB,cdf[i][n-NB],Parameters[i]);
        cdf[i][n-NB] = Parameters[i];
      }
    }
    n++;
  }

  rewind(infile);


  /* sort chains for cdfs */
  for(i=0; i<Nparam; i++) gsl_sort(cdf[i],1,NSAMP);

  
	/*  assign chain samples into bins for histograms  */
	norm = 1.0/((double)(count));

	n = 0;
	while(!feof(infile))
	{
    if(n<NB)fgets(burnrows,10000,infile);
    else
    {
      for(i=0; i<Nparam; i++)
      {
        fscanf(infile,"%lg",  &Parameters[i]);
        if(n%DS==0)
        {
          bin = (int)( (Parameters[i] - Pmin[i])/Pwidth[i] );
          if((bin >= 0 && bin < Nhist[i])) histograms[i][bin] += norm/Pwidth[i];
        }
      }
		}
    n++;
	}
	rewind(infile);

	for(i=0; i<Nparam; i++)
	{
		max = -1.0e60;
		for(j=0; j<Nhist[i]; j++)
		{
			if(histograms[i][j] > max)
			{
				max = histograms[i][j];
				bin = j;
			}
		}
		Pmode[i] = Pmin[i] + (double)bin*Pwidth[i];
	}

	/*  calculate std. deviation */
	printf("   Parameter  |  min  |  max  |  mean  |  mode | sigma\n");
	for(i=0; i<Nparam; i++)
	{
		Psigma[i] /= (double)(count-1);
		Psigma[i] = sqrt(Psigma[i]);
    KDEsigma[i] = Psigma[i]*pow((double)(count),-1./5.)*1.06*2;
		printf("   %i | %g | %g | %g | %g | %g\n", i+1, Pmin[i], Pmax[i],Pmean[i], Pmode[i], Psigma[i]);
	}

	/*  write histogram file  */
	sprintf(filename,"pdf_%s",argv[1]);
	outfile = fopen(filename,"w");
	fprintf(outfile,"# ");
	for(i=0; i<Nparam; i++)fprintf(outfile,"%.12g ",Pmode[i]);
	fprintf(outfile,"\n");
	for(i=0; i<Nmax; i++)
	{
		for(j=0;j<Nparam; j++)
		{
			x = Pmin[j] + (double)i*Pwidth[j];
			fprintf(outfile,"%.12e %e ", x, histograms[j][i]);
		}
		fprintf(outfile,"\n");
	}
	fclose(outfile);

	printf("  -Created file %s\n",filename);



  /*  write cdf file  */
  sprintf(filename,"cdf_%s",argv[1]);
  outfile = fopen(filename,"w");
  for(i=0; i<NSAMP; i++)
  {
    fprintf(outfile,"%.12g ",(double)i/(double)NSAMP);
    for(j=0;j<Nparam; j++)
    {
      fprintf(outfile,"%.12g ", cdf[j][i]);
    }
    fprintf(outfile,"\n");
  }
  fclose(outfile);

  printf("  -Created file %s\n",filename);



  //double KDEhistogram[Nparam][Nbin*oversample];
  double **KDEhistogram = (double **)malloc(Nparam*sizeof(double *));
  for(i=0; i<Nparam; i++) KDEhistogram[i] = (double *)malloc(Nbin*oversample*sizeof(double));


  for(i=0;i<Nparam;i++)
  {
    KDEmin[i] = Pmin[i]-3.0*KDEsigma[i];
    KDEmax[i] = Pmax[i]+3.0*KDEsigma[i];
    KDEwidth[i] = (KDEmax[i]-KDEmin[i])/(double)(Nbin*oversample);
    for(j=0;j<Nbin*oversample;j++)KDEhistogram[i][j]=0.0;
  }
  int jmin,jmax;
  n=0;
  while(!feof(infile))
	{
    if(n<NB)fgets(burnrows,10000,infile);
    else
    {

      for(i=0; i<Nparam; i++)
      {
        fscanf(infile,"%lg",  &Parameters[i]);
        if(n%DS==0)
        {
          jmin = (int)((Parameters[i] - KDEmin[i] - 3.0*KDEsigma[i])/(KDEwidth[i]));
          jmax = (int)((Parameters[i] - KDEmin[i] + 3.0*KDEsigma[i])/(KDEwidth[i]));
          for(j=jmin;j<jmax;j++)
          {
            KDEhistogram[i][j] += sqrt(2.)*gaussian(KDEmin[i]+(double)j*KDEwidth[i],Parameters[i],KDEsigma[i])/((double)count);
          }
        }
      }
		}
    n++;
	}
	fclose(infile);

  int min90flag;
  int max90flag;

  double min90;
  double max90;

  double total;
  double subtotal;

  sprintf(filename,"stats_%s",argv[1]);
	outfile = fopen(filename,"w");

  for(j=0;j<Nparam; j++)
  {
    min90flag=max90flag=0;
    total = 0.0;
    subtotal = 0.0;
    for(i=0; i<Nbin*oversample; i++) total+=KDEhistogram[j][i];
    for(i=0; i<Nbin*oversample; i++)
    {
      subtotal += KDEhistogram[j][i];
      if(subtotal/total>0.05 && min90flag==0)
      {
        min90 = KDEmin[j] + (double)i*KDEwidth[j];
        min90flag=1;
      }
      if(subtotal/total>0.95 && max90flag==0)
      {
        max90 = KDEmin[j] + (double)i*KDEwidth[j];
        max90flag=1;
      }
    }

    fprintf(outfile,"%i %g %g %g %g %g %g %g\n", j+1, Pmin[j], Pmax[j],Pmean[j], Pmode[j], Psigma[j],min90,max90);
    printf("Parameters %i 90%% credible interval:  [%.12g %.12g]\n",j,min90,max90);
	}
	fclose(outfile);
  
  
	/*  write kde file  */
	sprintf(filename,"kde_%s",argv[1]);
	outfile = fopen(filename,"w");
	fprintf(outfile,"# ");
	for(i=0; i<Nparam; i++)fprintf(outfile,"%.12g ",Pmode[i]);
	fprintf(outfile,"\n");


	for(i=0; i<Nbin*oversample; i++)
	{
		for(j=0;j<Nparam; j++)
		{
			x = KDEmin[j] + (double)(i)*KDEwidth[j];
			fprintf(outfile,"%.12e %e ", x, KDEhistogram[j][i]);
		}
		fprintf(outfile,"\n");
	}
	fclose(outfile);

	printf("  -Created file %s\n",filename);

	return 0;
}
double gaussian(double x, double mean, double sigma)
{
  return exp(-(x-mean)*(x-mean)/sigma/sigma)/sqrt(2.*M_PI*sigma*sigma);
}

