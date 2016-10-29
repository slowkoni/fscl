#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <math.h>

#include <kmacros.h>

#include "fscl.h"

static double ascprob_subsample(int k, int d, int min_obs, int n) {
  int i;
  double no_asc;

  no_asc = 0.;
  for(i=0;i<min_obs;i++) 
    no_asc += exp(lchoose(k,d-i) + lchoose(n-k,i)) + 
      exp(lchoose(n-k,d-i) + lchoose(k,i));
  no_asc /= exp(lchoose(n,d));

  return 1.0 - no_asc;
  /*  return 1.0 - ((exp(lchoose(k,d)) + exp(lchoose(n-k, d))) /
      exp(lchoose(n,d))); */
}

double *ascbias_adjust_background(double *bsf, int n, int asc_depth, 
				  int min_obs) {
  double *asc, *adj, asc_sum, adj_sum;
  int i;

  asc = alloca(sizeof(double)*(n+1));
  asc[0] = asc[n] = 0.;
  asc_sum = 0.;
  for(i=1;i<n;i++) {
    asc[i] = ascprob_subsample(i, asc_depth, min_obs, n);
    asc_sum += asc[i];
  }
  for(i=1;i<n;i++)
    asc[i] /= asc_sum;

  MA(adj, sizeof(double)*(n+1));
  adj[0] = adj[n] = 0.;
  adj_sum = 0.;
  for(i=1;i<n;i++) {
    adj[i] = bsf[i]/asc[i];
    adj_sum += adj[i];
  }

  for(i=1;i<n;i++)
    adj[i] /= adj_sum;

#if 0
  /* debugging block to assess if ascertainmnet bias correction to the 
     background frequency spectrum estimated from the data is producing
     something close to the 1/x expectation of the standard constant size
     neutral model. 

     NOTE: Adjustment appears to be producing correct spectra except for the
           upper tail where the relative error is as much as 50% -- not sure
	   if this will affect much as this adjusted background is used to
	   model the frequency spectrum of snps prior to the sweep, and high
	   derived frequency SNPs are most likely going to be fixed. If the
	   presence of fixed derived sites are not used, these will have
	   little impact. */
  {
    double *tmp[3];
    int sample_depths[3];
    scan_t scan_obj;

    sample_depths[0] = sample_depths[1] = sample_depths[2] = n;

    scan_obj.n_depths = 3;
    scan_obj.sample_depths = sample_depths;

    tmp[0] = bsf;
    tmp[1] = adj;

    tmp[2] = alloca(sizeof(double)*(n+1));
    tmp[2][0] = tmp[2][n] = 0;
    adj_sum = 0.;
    for(i=1;i<n;i++) {
      tmp[2][i] = 1/(double) i;
      adj_sum += tmp[2][i];
    }
    for(i=1;i<n;i++)
      tmp[2][i] /= adj_sum;

    output_background_fs("ascbias-background-adj.bs", &scan_obj, tmp);
  }
#endif

  return adj;

}

void ascbias_adjust_expect(double *fsp, int n, int min_obs, int d) {
  int i;
  double asc_sum;

  asc_sum = 0;
  for(i=0;i<=n;i++)
    asc_sum += fsp[i]*ascprob_subsample(i, d, min_obs, n);
    
  //  fsp[0] = fsp[n] = 0.;
  for(i=0;i<=n;i++)
    fsp[i] = fsp[i]*ascprob_subsample(i, d, min_obs, n)/asc_sum;

}

#if 0

static void ascbias_fsp_adjustment(double *fsp, int n_snps, int n, int d) {
  double *asc, asc_sum, inv;
  int i;

  asc = alloca(sizeof(double)*(n+1));
  asc_sum = 0;
  for(i=1;i<n;i++) {
    asc[i] = ascprob_subsample(i, d, n);
    asc_sum += fsp[i]*n_snps/asc[i];
  }
  inv = 1.0 - fsp[0] + fsp[n];

  for(i=1;i<n;i++)
    fsp[i] = fsp[i]*n_snps/asc[i]*(inv/asc_sum);

}


/* from pjh & spline routines */
    if (asc_depth > 0) {
      double inv;
      asc_sum = 0.;
      for(f=1;f<sample_size;f++) {
	asc[f] = ascprob_subsample(f, asc_depth, sample_size);
	asc_sum += p[f]*asc[f];
      }
      inv = 1.0 - p[0] - p[sample_size];
      
      for(f=1;f<sample_size;f++) {
       	p[f] = (p[f]*asc[f])*(inv/asc_sum);
      }
      
    }
#endif
