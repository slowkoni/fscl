
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <math.h>
#include <float.h>
#include <omp.h>

#include "kmacros.h"
#include "fscl.h"

static omp_lock_t thread_lock;
extern int spline_pts;
double log_ad_step;

double log_fact(int n) {
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

double spline_interpolate(spline_t *spf, double x) {
  int i;
  double sp_y;

  i = (x - LOG_AD_MIN)/log_ad_step;
  if (i >= spf->n) i = spf->n-1;
  if (i < 0) i = 0;
 
  sp_y = x*(spf->coef[i][0]*x*x + spf->coef[i][1]*x + spf->coef[i][2]) + 
    spf->coef[i][3];

  return sp_y;
}


static void solve_linear_system(double *b, double **m, double *v, int n) {
  int i, j ,k;
  double f;

  for(i=0;i<n;i++) {
    if (fabs(m[i][i]) < 1e-20) {
      int max;
      j = i+1;
      max = i;
      while(j<n) {
	if (fabs(m[j][i])>0 && 
	    (max == i || fabs(fabs(m[j][i])-1) < fabs(fabs(m[max][i])-1))) 
	     max = j;
	j++;
      }

      if (max == i)
	logmsg(MSG_FATAL,"Ill conditioned matrix while trying to estimate "
	       "spline functions to approximate sweep model likelihoods."); 
       
      
      for(k=0;k<n;k++) {
	m[i][k] += m[max][k];
      }
      v[i] += v[max];
    }

    f = m[i][i];
    for(k=i;k<i+8 && k<n;k++)
      m[i][k] /= f;

    v[i] /= f;

    for(j=i+1;j<i+8 && j<n;j++) {
      f = m[j][i];
      for(k=i;k<i+8 && k<n;k++)
	m[j][k] = m[j][k] - m[i][k]*f;
      v[j] = v[j] - v[i]*f;
    }
  }

  i = n - 1;
  while(i>=0) {
    if (fabs(m[i][i]) < 1e-10) {
      logmsg(MSG_WARN,"Warning: setting a spline coefficient %d to zero",i);
      b[i--] = 0;
      continue;
    }
    b[i] = v[i];
    for(k=i+1;k<i+8 && k<n;k++) {
      b[i] -= m[i][k]*b[k];
    }
    i--;
  }

}

static spline_t *estimate_spline(double *x, double *y, int n) {
  double **m, *v, *b;
  int i, j, k;
  spline_t *spline_func;

  MA(m, sizeof(double *)*(4*(n+1)));
  MA(v, sizeof(double)*(4*(n+1)));
  MA(b, sizeof(double)*(4*(n+1)));

  for(i=0;i<4*(n+1);i++) {
    MA(m[i], sizeof(double)*(4*(n+1)));
    for(j=0;j<4*(n+1);j++)
      m[i][j] = 0.;
    v[i] = 0.;
    b[i] = 0.;
  }

  
  m[0][0] = 6*x[0];
  m[0][1] = 2;

  for(i=1,j=0,k=0;k<n-1;i+=4,j+=4,k++) {
    m[i  ][j  ] = x[k]*x[k]*x[k];
    m[i  ][j+1] = x[k]*x[k];
    m[i  ][j+2] = x[k];
    m[i  ][j+3] = 1.;
    v[i] = y[k];

    m[i+1][j  ] = x[k+1]*x[k+1]*x[k+1];
    m[i+1][j+1] = x[k+1]*x[k+1];
    m[i+1][j+2] = x[k+1];
    m[i+1][j+3] = 1;
    v[i+1] = y[k+1]; 

    m[i+2][j  ] = 3*x[k+1]*x[k+1];
    m[i+2][j+1] = 2*x[k+1];
    m[i+2][j+2] = 1;
    m[i+2][j+3] = 0;
    m[i+2][j+4] = -3*x[k+1]*x[k+1];
    m[i+2][j+5] = -2*x[k+1];
    m[i+2][j+6] = -1;
    m[i+2][j+7] = 0;
    v[i+2] = 0;
    
    m[i+3][j  ] = 6*x[k+1];
    m[i+3][j+1] = 2;
    m[i+3][j+2] = 0;
    m[i+3][j+3] = 0;
    m[i+3][j+4] = -6*x[k+1];
    m[i+3][j+5] = -2;
    m[i+3][j+6] = 0;
    m[i+3][j+7] = 0;
    v[i+3] = 0;

  }

  m[i  ][j  ] = x[k]*x[k]*x[k];
  m[i  ][j+1] = x[k]*x[k];
  m[i  ][j+2] = x[k];
  m[i  ][j+3] = 1;
  v[i] = y[k];

  m[i+1][j  ] = x[n]*x[n]*x[n];
  m[i+1][j+1] = x[n]*x[n];
  m[i+1][j+2] = x[n];
  m[i+1][j+3] = 1;
  v[i+1] = y[n];

  m[i+2][j  ] = 6*x[n];
  m[i+2][j+1] = 2;
  m[i+2][j+2] = 0;
  m[i+2][j+3] = 0;
  v[i+2] = 0;

  solve_linear_system(b, m, v, 4*n);

  MA(spline_func, sizeof(spline_t));
  MA(spline_func->knot_points, sizeof(double)*(n+1));
  MA(spline_func->coef, sizeof(double *)*n);
  MA(spline_func->coef[0], sizeof(double)*4*n);

  for(i=1;i<n;i++)
    spline_func->coef[i] = spline_func->coef[i-1] + 4;


  spline_func->n = n;
  for(i=0;i<n;i++) {
    spline_func->knot_points[i] = x[i];
    for(k=0;k<4;k++) {
      spline_func->coef[i][k] = b[i*4 + k];
    }
  }
  spline_func->knot_points[n] = x[n];
  
  for(i=0;i<4*n;i++)
    free(m[i]);
  free(m);
  free(v);
  free(b);
  return spline_func;
}

static void sm_free(sm_ptable_t *sm_p) __attribute__((unused));
static void sm_free(sm_ptable_t *sm_p) {
  int i;
  for(i=1;i<sm_p->sample_size;i++) {
    free(sm_p->spline_func[i]->knot_points);
    free(sm_p->spline_func[i]->coef[0]);
    free(sm_p->spline_func[i]->coef);
    free(sm_p->spline_func[i]);
  }
  free(sm_p->spline_func);
  free(sm_p);
}


static double p_kescape(int k, int n, double ad) {

  if (k == 0) return exp(-n*ad);
  return exp(lchoose(n, k) + k * log(1.0 - exp(-ad)) - (n-k)*ad);
}

#define EXACT_KNOTS (100)
static double **pjh_use_splines(double *fsp, int n) {
  int n_knots, i, j, h, s;
  double spline_interval;
  double *xjh, *yjh, **pjh;
  spline_t *spf;

  n_knots = 500;//N_SPLINE_KNOTS;
  MA(xjh, sizeof(double)*(n_knots));
  MA(yjh, sizeof(double)*(n_knots));

  MA(pjh, sizeof(double *)*(n + 1));
  MA(pjh[0], sizeof(double)*((n+1)*(n+1)));
  for(j=0;j<=n;j++) {
    pjh[j] = pjh[0] + j*(n+1);

    cr_logmsg(MSG_DEBUG2,"Building pjh array %1.1f%%", j/(double) n * 100.0);
    if (n - j > n_knots) {
      spline_interval = (n - j - EXACT_KNOTS*2) /(double) (n_knots - EXACT_KNOTS*2);
      for(s=0;s<EXACT_KNOTS;s++)
	xjh[s] = j+s;
      while(s<n_knots-EXACT_KNOTS) {
	xjh[s] = (j+EXACT_KNOTS) + rint((s-EXACT_KNOTS)*spline_interval);
	s++;
      }
      while(s<n_knots-1) {
	xjh[s] = n - (n_knots - s);
	s++;
      }
      xjh[s] = n;

      for(s=0;s<n_knots;s++) {
	yjh[s] = 0.;

	for(i=j;i<=n;i++) {
	  yjh[s] += fsp[i]*exp(lchoose(i,j) +
			       lchoose(n - i, xjh[s] - j) -
			       lchoose(n, xjh[s]));
	}

      }
      spf = estimate_spline(xjh, yjh, n_knots-1);

      for(h=0;h<j;h++)
	pjh[j][h] = 0.;
      for(;h<=n;h++) {
	pjh[j][h] = spline_interpolate(spf, (double) h);
	if (isnan(pjh[j][h])) abort();
	if (pjh[j][h] < 0.) pjh[j][h] = 0.;
      }
      free(spf->knot_points);
      free(spf->coef[0]);
      free(spf->coef);
      free(spf);
    } else {
      //      for(h=0;h<j;h++)
      //	pjh[j][h] = 0.;
      for(h=0;h<=n;h++) {
	pjh[j][h] = 0.;
	for(i=j;i<=n;i++) {
	  pjh[j][h] += fsp[i]*exp(lchoose(i,j) +
			     lchoose(n - i, h - j) -
			     lchoose(n, h));
	}
      }
    }
  }
  logmsg(MSG_DEBUG2,"");
  free(xjh);
  free(yjh);

  return pjh;
}

#ifdef CUDA
void cuda_calculate_pjh(double *, double *, double *, int);
double *cuda_pjh_init_tables(int);
static double *cu_lft = NULL;
#endif

sm_ptable_t compute_sweep_model_fsp(double *fsp, int sample_size,
				    int asc_depth, int asc_min_freq,
				    int ascbias_background_only,
				    int include_invariant) {
  int b, i, j, k, q, f;
  double **pjh, **pbk, *x, **y, **fy, *p, p_sum;
  sm_ptable_t sm_ptable;
  static int warned = 0;

  log_ad_step = (LOG_AD_MAX - LOG_AD_MIN)/(spline_pts + 1.);

  if (0 && sample_size > 800) {
    omp_set_lock(&thread_lock);
    if (warned == 0) {
      logmsg(MSG_WARN,"sample depths are > 800, using splines to approximiate pjh "
	     "pre-computation.");
      warned = 1;
    }
    omp_unset_lock(&thread_lock);    
    pjh = pjh_use_splines(fsp, sample_size);
  } else {
    struct timeval stopwatch;
    gettimeofday(&stopwatch, NULL);
    fprintf(stderr,"Computing pjh[][] for sample size %d\n", sample_size);
    MA(pjh, sizeof(double *)*(sample_size + 1));
    MA(pjh[0], sizeof(double)*((sample_size+1)*(sample_size+1)));
    for(j=0;j<=sample_size;j++)
       pjh[j] = pjh[0] + j*(sample_size+1);
#ifndef CUDA
    for(j=0;j<=sample_size;j++) {
      int h;
      
      for(h=0;h<=sample_size;h++) {
	pjh[j][h] = 0.;

	for(i=j;i<=sample_size;i++) {
	  pjh[j][h] += fsp[i]*exp(lchoose(i,j) + 
				  lchoose(sample_size - i, h - j) -
				  lchoose(sample_size, h));
	}
      }
    }
#else
    cuda_calculate_pjh(pjh[0], fsp, cu_lft, sample_size);
#endif
    
    fprintf(stderr,"Done computing pjh[][] for sample size %d (%1.1f sec)\n",
	    sample_size, elapsed_time(&stopwatch));
  }

  MA(pbk, sizeof(double *)*(sample_size+1));
  MA(pbk[0], sizeof(double)*((sample_size+1)*(sample_size+1)));
  for(b=0;b<=sample_size;b++) {
    pbk[b] = pbk[0] + b*(sample_size+1);
    for(k=0;k<sample_size;k++) {
      pbk[b][k] = 0.;

      /* Suppose <k> lineages escape the sweep and <n - k> lineages do not. 
	 (That is, <n - k> lineages coalesce at the point in time the swept
	 allele arose, but <k> lineages coalesce at times prior to the sweep 
	 event). Assuming there is no coalescence
	 among the escaped lineages prior to (more recently than) to the MRCA 
	 of the swept lineages, then there are <k+1> lineages at the time just 
	 prior to the sweep. Therefore, if <j> of the <k+1> lineages at this 
	 point are mutant, the probability that the MRCA of the swept 
	 lineages is of the mutant type is <j>/(<k+1>).

	 pjh[j][h] is the probability of <j> mutant lineages in a sample of 
	 size <h> < <sample size>, given the background frequency spectrum 
	 observed for a sample of size <sample size>

	 if there are <b> observed mutant lineages in the present day sample
	 of size <sample size>, and <k> lineages escaped the sweep, then in 
	 the ancestoral sample of size <k+1> at the time the swept mutation 
	 arose, there were either b - (n - k) + 1 mutant lineages if the MRCA 
	 of the swept lineages is mutant (eg, the MRCA of the swept lineage 
	 plus whatever else is not amongst the <n-k> lineages derived from 
	 the swept lineage), or just <b> mutant lineages but in sample of 
	 size <k+1>. <pbk[b][k]> below holds the probability of observing <b>
	 mutant alleles in the present day sample, given <k> lineages escaped
	 and the frequency spectrum of the ancestoral sample given by 
	 <pjh[][]>

	 To get the expectation for present day frequecies, given a sweep of
	 strength <whatever> at a genetic distance of <whatever> from the snp
	 of interest, we must integrate over the distribution of the number of
	 escaped lineages, given the strength and genetic distance
	 parameters. This is captured by the composite parameter 
	 alpha*distance, expressed here as <ad>.
      */ 

      q = b - (sample_size - k) + 1;
      if (q > 0) {
	pbk[b][k] += pjh[q][k+1]*(q/(double) (k+1));
      }
      if (b < k+1) {
	pbk[b][k] += pjh[b][k+1]*((k+1-b)/(double) (k+1));
      }
    }
  }

  free(pjh[0]);
  free(pjh);


  sm_ptable.sample_size = sample_size;
  MA(sm_ptable.spline_func, sizeof(spline_t *)*(sample_size+1));
  MA(sm_ptable.fspline_func, sizeof(spline_t *)*(sample_size/2 + 1));

  MA(x, sizeof(double)*(spline_pts+1));
  MA(y, sizeof(double *)*(sample_size+1));
  MA(fy, sizeof(double *)*(sample_size/2 + 1));
  for(f=0;f<=sample_size;f++)
    MA(y[f], sizeof(double)*(spline_pts+1));
  for(f=0;f<=sample_size-f;f++)
    MA(fy[f], sizeof(double)*(spline_pts+1));

  MA(p, sizeof(double)*(sample_size+1));

  for(i=0; i<=spline_pts; i++) {
    double log_ad, ad;

    log_ad = LOG_AD_MIN + i*log_ad_step;
    ad = exp(log_ad);

    p_sum = 0.;
    for(f=0;f<=sample_size;f++) {
      p[f] = p_kescape(sample_size, sample_size, ad)*fsp[f];
      for(k=0;k<sample_size;k++) {
	p[f] += p_kescape(k, sample_size, ad)*pbk[f][k];
      }
      p_sum += p[f];
    }
    if (!include_invariant) {
      p_sum -= p[0] + p[sample_size];
      p[0] = p[sample_size] = 0.;
    }
    for(f=0;f<=sample_size;f++) p[f] /= p_sum;

    if (asc_depth > 0 && ascbias_background_only == 0)
      ascbias_adjust_expect(p, sample_size, asc_min_freq, asc_depth);

    for(f=0;f<=sample_size;f++) {
      if (p[f] == 0.) y[f][i] = log(DBL_MIN);
      else y[f][i] = log(p[f]);
    }

    for(f=0;f<sample_size - f;f++) {
      if (p[f] + p[sample_size - f] == 0.) fy[f][i] = log(DBL_MIN);
      else fy[f][i] = log(p[f] + p[sample_size - f]);
    }
    if (f == sample_size - f) {
      if (p[f] == 0.) fy[f][i] = log(DBL_MIN);
      else fy[f][i] = log(p[f]);
    }

    x[i] = log_ad;
    logmsg(MSG_DEBUG1,"%5d %5d %5.2f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f", sample_size, i, x[i], p[0], p[1], p[2], p[sample_size/2], p[sample_size-2], p[sample_size-1],p[sample_size]);
  }
  logmsg(MSG_DEBUG1,"bgrnd %d %5.3f %5.3f %5.3f %5.3f", sample_size, fsp[0], fsp[1], fsp[sample_size-1], fsp[sample_size]);

  sm_ptable.pbk = pbk;
  sm_ptable.fsp = fsp;

  for(f=0;f<=sample_size;f++)
    sm_ptable.spline_func[f] = estimate_spline(x, y[f], spline_pts);
  for(f=0;f<=sample_size - f;f++)
    sm_ptable.fspline_func[f] = estimate_spline(x, fy[f], spline_pts);

  for(f=0;f<=sample_size;f++)
    free(y[f]);
  for(f=0;f<sample_size - f;f++)
    free(fy[f]);

  free(x);
  free(y);
  free(fy);
  free(p);

  return sm_ptable;
}

sm_ptable_t *compute_sweep_model_tables(scan_t *scan_obj, double **fsp,
					int asc_depth, int asc_min_freq,
					int ascbias_background_only,
					int include_invariant) {
  int i;
  sm_ptable_t *sm_p;
  double *asc;
  int n_complete, max_depth;
  
  omp_init_lock(&thread_lock);
  n_complete = 0;
  MA(sm_p, sizeof(sm_ptable_t)*scan_obj->n_depths);
  max_depth = scan_obj->sample_depths[0];
  for(i=0; i < scan_obj->n_depths; i++)
    if (max_depth < scan_obj->sample_depths[i]) max_depth = scan_obj->sample_depths[i];
  
#ifdef CUDA
  cu_lft = cuda_pjh_init_tables(max_depth+1);
#endif
  
#pragma omp parallel for schedule(dynamic, 1)
  for(i=0;i<scan_obj->n_depths;i++) {
    if (asc_depth > 0) {
      asc = ascbias_adjust_background(fsp[i], scan_obj->sample_depths[i],
				      asc_depth, asc_min_freq);
    } else {
      asc = fsp[i];
    }
    sm_p[i] = compute_sweep_model_fsp(asc, scan_obj->sample_depths[i],
				      asc_depth, asc_min_freq, 
				      ascbias_background_only,
				      include_invariant);
    omp_set_lock(&thread_lock);
    n_complete++;
    cr_logmsg(MSG_STATUS,"Computing sweep models for all sample depths - %1.1f%% ",
	      n_complete/(double) scan_obj->n_depths * 100.);
    omp_unset_lock(&thread_lock);
  }
  logmsg(MSG_STATUS,"");
  omp_destroy_lock(&thread_lock);
  return sm_p;
}



