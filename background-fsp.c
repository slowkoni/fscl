
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <kmacros.h>

#include <math.h>
#include <float.h>
#include <unistd.h>

#include "fscl.h"

static double log_fact(int n) {
  static int n_max = 0;
  static double *lookup_table = NULL;
  int i;

  if (n < 0) return -DBL_MAX;

  if (n == 0 || n == 1) return 0.;

  if (n > n_max) {
    RA(lookup_table, sizeof(double)*(n+1));
    if (n_max == 0) {
      lookup_table[0] = 0;
      n_max = 1;
    }
    for(i=n_max;i<=n;i++)
      lookup_table[i] = lookup_table[i-1] + log(i);
    n_max = n;
  }

  return lookup_table[n];
}

double lchoose(int n, int k) {

  if (n == 0 && k == 0) return 0.;
  if (k > n || n == 0) return -DBL_MAX;
  return log_fact(n) - log_fact(k) - log_fact(n-k);
}


static double **neutral_spectra(scan_t *scan_obj) {
  int n_invariant, n_fixed, n_segregating, i, k, m;
  double segregating_sum, **fsp;

  n_invariant = 0;
  n_fixed = 0;
  for(i=0;i<scan_obj->n_snps;i++) {
    if (scan_obj->snps[i].obs_freq == 0) n_invariant++;
    if (scan_obj->snps[i].obs_freq == 
	scan_obj->sample_depths[scan_obj->snps[i].depth_p])
      n_fixed++;
  }

  MA(fsp, sizeof(double *)*scan_obj->n_depths);
  for(i=0;i<scan_obj->n_depths;i++) {
    m = scan_obj->sample_depths[i];
    MA(fsp[i], sizeof(double)*(m+1));
    fsp[i][0] = n_invariant;
    fsp[i][m] = n_fixed;

    segregating_sum = 0.;
    for(k=1;k<m;k++)
      segregating_sum += 1/(double) k;

    n_segregating = scan_obj->n_snps - n_fixed - n_invariant;
    for(k=1;k<m;k++)
      fsp[i][k] = (1./(double) k)/segregating_sum * n_segregating;
    for(k=0;k<=m;k++)
      fsp[i][k] /= (double) scan_obj->n_snps;
  }

  return fsp;
}

static void hypergeometric_upsample(double *fsp, double weight, int N, int k, 
				    int n) __attribute__((unused));
static void hypergeometric_upsample(double *fsp, double weight,
				    int N, int k, int n) {
  int m;

  if (weight == 0.) return;
  /* N: total balls in urn
     k,n: number of white balls in a sample of n < N drawn
     n-k: number of black balls drawn in a sample of n
     m: number of white balls in the urn. m >= k && N - m >= n - k */
  m = k;
  while(N - m >= n - k) {
    fsp[m] += exp(lchoose(m, k) + lchoose(N-m,n-k) - lchoose(N,n))*weight;
    m++;
  }

}

static void hypergeometric_downsample_fsp(double *d_fsp, double *fsp, 
					  int n, int N, int include_invariant) {
  int m, k;

  if (include_invariant) {
    for(m=0;m<=N;m++) {
      for(k=0;k<=m && k<=n;k++)
	d_fsp[k] += exp(lchoose(m,k) + lchoose(N-m, n-k) -lchoose(N, n))*fsp[m];
    }
  } else {
    for(m=1;m<=N;m++) {
      for(k=1;k<=m && k<n;k++)
	d_fsp[k] += exp(lchoose(m,k) + lchoose(N-m, n-k) -lchoose(N, n))*fsp[m];
    }
  }

}

static void binomial_sampling(double *fsp, double weight, int m, int k, 
			      int n) __attribute__((unused));
static void binomial_sampling(double *fsp, double weight, int m, int k, 
			      int n) {
  int i;
  double p, q, b;

  if (weight == 0.) return;

  if (k == 0) {
    fsp[0] += weight;
    return;
  }
  if (k == n) {
    fsp[m] += weight;
    return;
  }

  p = log(k/(double) n);
  q = log(1. - k/(double) n);

  i = (int) (m*(k/(double) n));
  while(i>=0) {
    b = exp(lchoose(m, i) + i * p + (m-i) * q);
    if (b == 0.) break;
    fsp[i] += b*weight;
    i--;
  }
  i = (int) (m*(k/(double) n)) + 1;
  while(i<=m) {
    b = exp(lchoose(m, i) + i * p + (m-i) * q);
    if (b == 0.) break;
    fsp[i] += b*weight;
    i++;
  }
}

static double **load_spectra(scan_t *scan_obj, char *bs_fname) {
  int i, j, depth;
  char *p, *q;
  char inputline[131072];
  double **fsp;
  FILE *f;

  MA(fsp, sizeof(double *)*scan_obj->n_depths);
  for(i=0;i<scan_obj->n_depths;i++)
    fsp[i] = NULL;

  f = fopen(bs_fname, "r");
  if (f == NULL) 
    logmsg(MSG_FATAL,"Can't background frequency spectrum file \"%s\" (%s)\n",
	   bs_fname, strerror(errno));

  while(fgets(inputline, 8192, f) != NULL) {
    CHOMP(inputline);
    if (inputline[0] == 0) continue;
    if (inputline[0] == '#') continue;

    p = inputline;
    q = strsep(&p, " \t");
    depth = atoi(q);

    while(i<scan_obj->n_depths && scan_obj->sample_depths[i] != depth) i++;

    if (i == scan_obj->n_depths) {
      logmsg(MSG_STATUS,"Frequency spectrum for %d classes not required by "
	     "snp data.\n", depth);
      continue;
    }

    j = 0;
    while((q = strsep(&p, " \t"))) {
      if (j % 128 == 0) RA(fsp[i], sizeof(double)*(j + 128));
      fsp[i][j++] = strtod(q, NULL);
    }
    if (j != depth) 
      logmsg(MSG_FATAL, "Error: Frequency spectrum on line %d indicates "
	     "%d classes but %d were found.\n", depth, j);

  }

  for(i=0;i<scan_obj->n_depths;i++) {
    if (fsp[i] == NULL) {
      logmsg(MSG_FATAL, "Error: data requires background frequency spectrum "
	     "for sample depth %d, not found %s\n", scan_obj->sample_depths[i],
	     bs_fname);
    }
  }

  return fsp;
}

double **background_fsp(scan_t *scan_obj, int force_neutral_spectrum,
			char *background_fsfname, int include_invariant) {
  int m, k, i, depth, max_depth;
  int __attribute__((unused))max_depth_p;
  double fsp_sum, **fsp, *tmp_fsp, wa, wd;

  if (force_neutral_spectrum) return neutral_spectra(scan_obj);
  if (background_fsfname) return load_spectra(scan_obj, background_fsfname);

  fprintf(stderr,"Estimating background site frequency spectrum....   \n");
  MA(fsp, sizeof(double *)*scan_obj->n_depths);
  max_depth_p = -1;
  max_depth = -1000;
  for(m=0;m<scan_obj->n_depths;m++) {
    MA(fsp[m], sizeof(double)*(scan_obj->sample_depths[m]+1));
    for(k=0;k<=scan_obj->sample_depths[m];k++) fsp[m][k] = 0.;

    if (scan_obj->sample_depths[m] > max_depth) {
      max_depth = scan_obj->sample_depths[m];
      max_depth_p = m;
    }
  }
  log_fact(max_depth+1);
  fprintf(stderr,"%d distinct sample depths observed. Maximum sample depth is %d haplotypes.\n", scan_obj->n_depths, max_depth);
  fprintf(stderr,"log(%d!) = %1.1f\n", max_depth, log_fact(max_depth)); 
  
  MA(tmp_fsp, sizeof(double)*(max_depth+1));
  for(k=0;k<=max_depth;k++) tmp_fsp[k] = 0.;
  for(i=0;i<scan_obj->n_snps;i++) {
    //    if (scan_obj->snps[i].folded) continue;
    depth = scan_obj->sample_depths[scan_obj->snps[i].depth_p];
    if (scan_obj->snps[i].folded) {
      if (scan_obj->snps[i].obs_freq == 0) {
	wa = 1;
	wd = 0;
      } else if (scan_obj->snps[i].obs_freq == depth) {
	wa = 0;
	wd = 1;
      } else {
	wa = 1./(scan_obj->snps[i].obs_freq);
	wd = 1./(depth - scan_obj->snps[i].obs_freq);
      }
    } else {
      wd = 1.;
      wa = 0.;
    }

    if (depth == max_depth) {
      tmp_fsp[scan_obj->snps[i].obs_freq] += wa/(wa + wd);
      tmp_fsp[depth - scan_obj->snps[i].obs_freq] += wd/(wa + wd);
    } else {
      //      hypergeometric_upsample(tmp_fsp, wd/(wa + wd), max_depth,
      //		      scan_obj->snps[i].obs_freq, depth);
      //hypergeometric_upsample(tmp_fsp, wa/(wa + wd), max_depth,
      //			      depth - scan_obj->snps[i].obs_freq, depth);
    }
  }

  fsp_sum = 0.;
  for(k=0;k<=max_depth;k++) fsp_sum += tmp_fsp[k];
  for(k=0;k<=max_depth;k++) tmp_fsp[k] /= fsp_sum;
  fprintf(stderr,"Total SNPs observed at max depth %d is %1.1f (%1.1f%%)\n", max_depth, fsp_sum, fsp_sum/(double) scan_obj->n_snps * 100.);  

  for(m=0;m<scan_obj->n_depths;m++) {
    if (isatty(2)) fprintf(stderr,"\r");
    fprintf(stderr,"Estimating frequency spectrum for sample depth %6d  (%1.1f%%)", scan_obj->sample_depths[m], (m+1)/(double) scan_obj->n_depths * 100.);
    if (!isatty(2)) fprintf(stderr,"\n");
    depth = scan_obj->sample_depths[m];
    //    if (depth < max_depth) {
      hypergeometric_downsample_fsp(fsp[m], tmp_fsp, depth,
				    max_depth, include_invariant);
    
      fsp_sum = 0.;
      for(k=0;k<=depth;k++) fsp_sum += fsp[m][k];
      for(k=0;k<=depth;k++) fsp[m][k] /= fsp_sum;
      //   }
  }
  if (isatty(2)) fprintf(stderr,"\n");
  fprintf(stderr,"\nDone estimating background frequency spectra.\n");

#if 0
  for(m=0;m<scan_obj->n_depths;m++) {
    depth = scan_obj->sample_depths[m];
    MA(fsp[m], sizeof(double)*(depth+1));

    for(k=0;k<=depth;k++) fsp[m][k] = 0.;

    for(i=0;i<scan_obj->n_snps;i++) {
      if (scan_obj->snps[i].folded) {
	if (scan_obj->snps[i].obs_freq == 0) {
	  wa = 1;
	  wd = 0;
	} else if (scan_obj->snps[i].obs_freq == depth) {
	  wa = 1;
	  wd = 0;
	} else {
	  wa = 1./(scan_obj->snps[i].obs_freq);
	  wd = 1./(depth - scan_obj->snps[i].obs_freq);
	}
      } else {
	wd = 1.;
	wa = 0.;
      }

      if (scan_obj->snps[i].depth_p == m) {
	fsp[m][scan_obj->snps[i].obs_freq] += wd/(wa+wd);
	fsp[m][depth - scan_obj->snps[i].obs_freq] += wa/(wa+wd);
      } else {
	binomial_sampling(fsp[m], wd/(wa+wd), depth, 
			  scan_obj->snps[i].obs_freq, 
			  scan_obj->sample_depths[scan_obj->snps[i].depth_p]);
	if (wa > 0.)
	  binomial_sampling(fsp[m], wa/(wa+wd), depth, 
			  scan_obj->sample_depths[scan_obj->snps[i].depth_p] -
			  scan_obj->snps[i].obs_freq, 
			  scan_obj->sample_depths[scan_obj->snps[i].depth_p]);
      }
    }

    fsp_sum = 0.;
    for(k=0;k<=scan_obj->sample_depths[m];k++) fsp_sum += fsp[m][k];
    for(k=0;k<=scan_obj->sample_depths[m];k++) fsp[m][k] /= fsp_sum;
  }    
#endif

  return fsp;
}

void output_background_fs(char *fname, scan_t *scan_obj, double **fsp) {
  FILE *f;
  int i, j;

  f = fopen(fname, "w");
  if (f == NULL) {
    fprintf(stderr,"Can't open background frequency spectrum file \"%s\" for "
	    "output. (%s)\n", fname, strerror(errno));
    exit(-1);
  }

  for(i=0;i<scan_obj->n_depths;i++) {
    fprintf(f,"%d",scan_obj->sample_depths[i]);
    for(j=0;j<=scan_obj->sample_depths[i];j++) fprintf(f,"\t%1.6f",fsp[i][j]);

    fprintf(f,"\n");
  }

  fclose(f);
}


