#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <math.h>
#include <float.h>

#include <kmacros.h>
#include "fscl.h"

static double *log_table = NULL;

void init_log_table(void) {
  int i;

  MA(log_table, sizeof(double)*0x10000);
  for(i=1;i<=0xFFFF;i++)
    log_table[i] = log(i);

  /* zero corresponds to the case where the sweep target is directly on top of
     a SNP. For this case, we assume the SNP (or the target) is at the
     adjacent base, so we use this hack to make the log of a "zero distance"
     the same as the quantity for 1bp apart */
  log_table[0] = 0.; 
}

/* the log() function in the standard math library is rather slow, but we need
   to take the log of integer values hundreds of times per each likelihood
   evaluation. This is a simple optimization for integer valued arguments with
   a speed-space-accuracy trade off. A table of pre-computed values is used
   for small integer arguments, and these values plus a constant is used for
   higher values by first dividing (left shift here) by an appropriate
   value. The approximation is due to the low-order bits being ignored for
   larger valued arguments. The difference between this function and log() is
   less than <0.01%

   Note that this function is only used for taking the log of the (integer) 
   distance between two SNPs. */
static double logt(int d) {

  if (d > 0xFFFFFF) return 11.783502069519070 + log_table[d >> 16];
  if (d > 0xFFFF) return 5.545177444479562 + log_table[d >> 8 ];
  return log_table[d];
}

static double p_kescape(int k, int n, double ad) __attribute__((unused));
static double p_kescape(int k, int n, double ad) {

  if (k == 0) return exp(-n*ad);
  return exp(lchoose(n, k) + k * log(1.0 - exp(-ad)) - (n-k)*ad);
}

#if 0
extern double log_ad_step;
inline static double spline_interpolate(spline_t *spf, double x) {
  int i, j, m;
  double sp_y, xsq;

#if 0
  i = 0;
  j = spf->n;
  while(j-i>1) {
    m = (i+j)/2;
    if (x < spf->knot_points[m]) {
      j = m;
    } else {
      i = m;
    }
  }
#endif
  i = (x - LOG_AD_MIN)/log_ad_step;
  if (i >= spf->n) i = spf->n-1;
  if (i < 0) i = 0;

  xsq = x*x;
  sp_y = spf->coef[i][0]*x*xsq + spf->coef[i][1]*xsq + spf->coef[i][2]*x +
    spf->coef[i][3];

  return sp_y;
}
#endif

static double snp_likelihood(snp_t *snp, double log_ad, sm_ptable_t *sm_p) {
  int f;
  sm_ptable_t *sm_fsp;
  double logl;
  
  //if (log_ad > LOG_AD_MAX) log_ad = LOG_AD_MAX;
  //if (log_ad < LOG_AD_MIN) log_ad = LOG_AD_MIN;

  sm_fsp = sm_p + snp->depth_p;
  f = snp->obs_freq;

  if (snp->folded) {
    logl = spline_interpolate(sm_fsp->fspline_func[f], log_ad);
  } else {
    logl = spline_interpolate(sm_fsp->spline_func[f], log_ad);
  }

  return logl - snp->null_logl;
}

static void sm_likelihood(scan_pt_t *result, snp_t *snps, sm_ptable_t *sm_p) {
  int i;
  double log_ad;

  /*  fprintf(stdout,"sm_likelihood search at %d, lalpha %1.2f:\t%d\t%d\t%d\n",
	  result->sweep_pos, result->lalpha, result->window_start, 
	  result->nearest_snp, result->window_end); */

  result->sm_logl = result->null_logl;
  log_ad = logt(abs(result->sweep_pos - snps[result->nearest_snp].pos)) +
    result->lalpha;
  if (log_ad > LOG_AD_MAX) {
    //result->sm_logl = result->null_logl;
    return;
  }
  result->sm_logl += snp_likelihood(snps + result->nearest_snp, log_ad, sm_p);

  i = result->nearest_snp - 1;
  while(i>=result->window_start) {
    log_ad = logt(result->sweep_pos - snps[i].pos) + result->lalpha;
    if (log_ad > LOG_AD_MAX) break;
    result->sm_logl += snp_likelihood(snps+i, log_ad, sm_p);
    i--;
  }
#if 0
  while(i>=result->window_start) {
    result->sm_logl += snps[i].null_logl;
    i--;
  }
#endif

  i = result->nearest_snp+1;
  while(i<=result->window_end) {
    log_ad = logt(snps[i].pos - result->sweep_pos) + result->lalpha;
    if (log_ad > LOG_AD_MAX) break;
    result->sm_logl += snp_likelihood(snps+i, log_ad, sm_p);
    i++;
  }
#if 0
  while(i<=result->window_end) {
    result->sm_logl += snps[i].null_logl;
    i++;
  }
#endif

}

#if 0

static void rsearch_maxalpha(scan_pt_t *result, snp_t *snps,
			     sm_ptable_t *sm_p, scan_pt_t l_result,
			     scan_pt_t r_result, int recursion_limit,
			     int depth) {
  scan_pt_t m1_result, m2_result;
  double s1, s2, s3;

  if ((recursion_limit && depth > recursion_limit) || 
      r_result.lalpha - l_result.lalpha < 1e-2) {
    //    if (exp(r_result.lalpha) - exp(l_result.lalpha) < 1e-3) {
    //  fprintf(stderr,"Exiting recursion by limited range of lalpha\n");
    //}
    if (l_result.sm_logl > r_result.sm_logl) {
      *result = l_result;
    } else {
      *result = r_result;
    }
    return;
  }
 
  //  fprintf(stdout,"rsearch called: %1.2f\t%1.2f\n", l_result.lalpha, r_result.lalpha);
  m1_result = m2_result = l_result;
  m1_result.lalpha = l_result.lalpha + (r_result.lalpha - l_result.lalpha)/3.0;
  m2_result.lalpha = l_result.lalpha + (r_result.lalpha - l_result.lalpha)*2.0/3.0;

  sm_likelihood(&m1_result, snps, sm_p);
  sm_likelihood(&m2_result, snps, sm_p);

  if (fabs(l_result.sm_logl  - m1_result.sm_logl)  < 1e-3 && 
      fabs(m1_result.sm_logl - m2_result.sm_logl)  < 1e-3 &&
      fabs(m2_result.sm_logl - r_result.sm_logl)   < 1e-3) {
    //fprintf(stderr,"Exiting recursion by matching likelihoods\n");
    *result = l_result;
    result->lalpha = (r_result.lalpha + l_result.lalpha)/2.0;
    return;
  }

  s1 = l_result.sm_logl + m1_result.sm_logl;
  s2 = m1_result.sm_logl + m2_result.sm_logl;
  s3 = m2_result.sm_logl + r_result.sm_logl;

  if (s1 > s2 && s1 > s3) {
    rsearch_maxalpha(result, snps, sm_p, l_result, m2_result,
		     recursion_limit, depth+1);
  } else if (s3 > s2) {
    rsearch_maxalpha(result, snps, sm_p, m1_result, r_result,
		     recursion_limit, depth+1);
  } else {
    rsearch_maxalpha(result, snps, sm_p, m1_result, m2_result,
		     recursion_limit, depth+1);
  }
}

void search_maxalpha(scan_pt_t *result, snp_t *snps, sm_ptable_t *sm_p) {
  scan_pt_t l_result, r_result, tmp_result;

  result->sm_logl = -DBL_MAX;
  l_result = r_result = tmp_result = *result;

  r_result.lalpha = LOG_AD_MIN;
  sm_likelihood(&r_result, snps, sm_p);
  while(r_result.lalpha < LOG_AD_MAX) {

    l_result = r_result;
    if (r_result.lalpha == LOG_AD_MIN) {
      r_result.lalpha = -15.0;
    } else if (r_result.lalpha < -10.0) {
      r_result.lalpha += 1.5;
    } else if (r_result.lalpha < -6.0) {
      r_result.lalpha += 0.5;
    } else if (r_result.lalpha < -2.0) {
      r_result.lalpha += 1.5;
    } else {
      r_result.lalpha += 2.0;
    }
    sm_likelihood(&r_result, snps, sm_p);

    //    fprintf(stdout,"rsearch top loop\n");
    rsearch_maxalpha(&tmp_result, snps, sm_p, l_result, r_result, 1,1);
    if (tmp_result.sm_logl > result->sm_logl) *result = tmp_result;
  }

  if (result->lalpha < -16.0) {
    l_result.lalpha = LOG_AD_MIN;
    r_result.lalpha = -16.0;
  } else if (result->lalpha < -10.0) {
    l_result.lalpha = result->lalpha - 1.0;
    r_result.lalpha = result->lalpha + 1.0;
  } else if (result->lalpha < -4.0) {
    l_result.lalpha = result->lalpha - 0.5;
    r_result.lalpha = result->lalpha + 0.5;
  } else if (result->lalpha < 0.0) {
    l_result.lalpha = result->lalpha - 0.1;
    r_result.lalpha = result->lalpha + 0.1;
  } else {
    l_result.lalpha = result->lalpha - 1.0;
    r_result.lalpha = result->lalpha + 1.0;
  }

  //  fprintf(stdout,"rsearch bottom loop - left\n");
  sm_likelihood(&l_result, snps, sm_p);
  //fprintf(stdout,"rsearch bottom loop - right\n");
  sm_likelihood(&r_result, snps, sm_p);
  rsearch_maxalpha(&tmp_result, snps, sm_p, l_result, r_result, 3, 0);

  //  if (tmp_result.sm_logl < tmp_result.null_logl) 
  //  tmp_result.sm_logl = tmp_result.null_logl;

  tmp_result.clr = tmp_result.sm_logl - tmp_result.null_logl;

  *result = tmp_result;
}

#else

void search_maxalpha(scan_pt_t *result, snp_t *snps, sm_ptable_t *sm_p) {
  double lalpha, step_size, left_endpoint, right_endpoint;
  scan_pt_t tmp_result, max_result;

  max_result = tmp_result = *result;
  max_result.sm_logl = -DBL_MAX;

  step_size = (LOG_AD_MAX - LOG_AD_MIN)/10.0;

  for(lalpha = LOG_AD_MIN; lalpha <= LOG_AD_MAX; lalpha += step_size) {
    tmp_result.lalpha = lalpha;
    sm_likelihood(&tmp_result, snps, sm_p);
    if (tmp_result.sm_logl > max_result.sm_logl) max_result = tmp_result;
  }

  left_endpoint = max_result.lalpha - step_size;
  if (left_endpoint < LOG_AD_MIN) left_endpoint = LOG_AD_MIN;
  right_endpoint = max_result.lalpha + step_size;
  if (right_endpoint > LOG_AD_MAX) right_endpoint = LOG_AD_MAX; 

  step_size = (right_endpoint - left_endpoint)/15.;

  for(lalpha = left_endpoint + step_size; lalpha < right_endpoint; 
      lalpha += step_size) {
    tmp_result.lalpha = lalpha;
    sm_likelihood(&tmp_result, snps, sm_p);
    if (tmp_result.sm_logl > max_result.sm_logl) max_result = tmp_result;
  }

  max_result.clr = max_result.sm_logl - max_result.null_logl;
  *result = max_result;
}
#endif

