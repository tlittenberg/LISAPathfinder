
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "TimePhaseMaximization.h"

static void dfour1(double data[], unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi, swap;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m > 1 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

static void drealft(double data[], unsigned long n, int isign)
{
  void dfour1(double data[], unsigned long nn, int isign);
  unsigned long i,i1,i2,i3,i4,np3;
  double c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;

  theta=3.141592653589793/(double) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    dfour1(data,n>>1,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for (i=2;i<=(n>>2);i++) {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    dfour1(data,n>>1,-1);
  }
}

static void max_array_element(double *max, int *index, double *array, int n)
{
  int i;

  int j;
  double temp = -1.0e60;

  for(i = 0; i < n; i++)
  {
    if(array[i] > temp)
    {
      temp = array[i];
      j = i;
    }
  }

  *index = j;
  *max = temp;

}

static void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, double *psd, int n)
{
  int nb2, i, l, k;

  nb2 = n / 2;

  corr[0]  = 0.0;
  corr[1]  = 0.0;
  corrf[0] = 0.0;
  corrf[1] = 0.0;

  for (i=1; i < nb2; i++)
  {
    l=2*i;
    k=l+1;

    corr[l]	 =  ( data1[l]*data2[l] + data1[k]*data2[k]) / psd[i];
    corr[k]	 = -( data1[k]*data2[l] - data1[l]*data2[k]) / psd[i];
    corrf[l] =  ( data1[l]*data2[k] - data1[k]*data2[l]) / psd[i];
    corrf[k] = -( data1[k]*data2[k] + data1[l]*data2[l]) / psd[i];
  }

  drealft(corr-1, n, -1);
  drealft(corrf-1, n, -1);

}

void Sum_Extreme(double *a, double *b, double *Sn, int n, double *delt, double *pshift, double Tobs, double t0)
{
  double max=0.0;
  int i;
  int index = 0;
  double *AC, *AF;
  double *corr;
  double *corrcos, *corrsin;

  double newtime = t0;

  // this version simultaneously maximizes across the Network

  AC = malloc(sizeof(double)*n);
  AF = malloc(sizeof(double)*n);

  corr    = malloc(sizeof(double)*n);
  corrcos = malloc(sizeof(double)*n);
  corrsin = malloc(sizeof(double)*n);

  for(i = 0; i < n; i++)
  {
    AC[i] = 0.0;
    AF[i] = 0.0;
    corr[i] = 0.0;
    corrcos[i] = 0.0;
    corrsin[i] = 0.0;
  }

  phase_blind_time_shift(AC, AF, a, b, Sn, n);

  for(i = 0; i < n; i++)
  {
    corrcos[i] += AC[i];
    corrsin[i] += AF[i];
  }

  for(i = 0; i < n; i++) corr[i] = sqrt(corrcos[i]*corrcos[i]+corrsin[i]*corrsin[i]);


  *delt=0.0;
  newtime = -1.0;//t0+*delt;

  while(newtime<0.0 || newtime>Tobs)
  {
    max_array_element(&max, &index, corr, n);

    if(index < (n/2))
      *delt = ( (double)index/(double)n )*Tobs;
    else if(index >= n/2)
      *delt = ( (double)(index - n)/(double)n )*Tobs;

    max = 2.0*(max);
    *pshift = atan2(corrsin[index],corrcos[index]);

    newtime = *delt + t0;
    corr[index]=-1.0;
  }

  free(AC);
  free(AF);
  free(corr);
  free(corrcos);
  free(corrsin);
  
}