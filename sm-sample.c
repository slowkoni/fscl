#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <math.h>
#include <float.h>

#include <time.h>
#include <gsl/gsl_rng.h>

#include <kmacros.h>
#include <cmdline-utils.h>

#include "fscl.h"

sm_ptable_t compute_sweep_model_fsp(double *fsp, int sample_size,
				    int asc_depth);

int spline_pts = 201;

static char *output_basename;
static int output_complete, bp_length, n_sweeps, asc_depth, sample_size;
static double alpha, snp_density, mb_length;

static option_t options[] = {
  { 'o', "output-basename", (void *) &output_basename, OPT_STR, 1, 1,
    "basename for output files" },
  { 0, "output-complete", (void *) &output_complete, OPT_FLAG, 0, 0,
    "output the full sample of snps as well as the ascertained and "
    "snp-density matched random sample" },
  { 'a', "alpha", (void *) &alpha, OPT_DBL, 1, 1, 
    "strength of selection at selected sites" },
  { 's', "snp-density", (void *) &snp_density, OPT_DBL, 1, 1,
    "number of SNPs per kb in master sample" },
  { 'd', "asc-depth", (void *) &asc_depth, OPT_INT, 1, 1,
    "depth of ascertainment sample" },
  { 'N', "sample-size", (void *) &sample_size, OPT_INT, 1, 1,
    "number of chromosomes in generated SNP frequency sample" },
  { 'l', "segment-length", (void *) &mb_length, OPT_DBL, 1, 1,
    "length of segment in megabases" },
  { 'n', "n-sweeps", (void *) &n_sweeps, OPT_INT, 1, 1,
    "number of sweeps positions to simulate" },
  { 0, NULL, NULL, 0, 0, 0, NULL }
};

static gsl_rng *rng;

typedef struct {
  int pos;
  double alpha;
} spos_t;

typedef struct {
  int pos;
  int f;
} sm_snp_t;

static void init_options(void) {

  rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, time(NULL));

  snp_density = 1.;
  mb_length = 10.;
  alpha = 1e-6;
  output_complete = 0;
  n_sweeps = 1;
  asc_depth = 2;
  sample_size = 128;
  output_basename = NULL;
}

static void verify_options(void) {
  if (output_basename == NULL) {
    fprintf(stderr,"Specify basename for output files with -o option.\n");
    exit(-1);
  }
  if (snp_density <= 0.) {
    fprintf(stderr,"Specify positive SNP/kb density with -s option.\n");
    exit(-1);
  }
  if (mb_length <= 0.) {
    fprintf(stderr,"Specify positive segment length in Mb with -l option.\n");
    exit(-1);
  }
  if (alpha <= 0.0) {
    fprintf(stderr,"Specify positive selection strength with -a option.\n");
    exit(-1);
  }
  if (asc_depth < 2) {
    fprintf(stderr,"Specify ascertainment sample size > 1 with -d option.\n");
    exit(-1);
  }
  if (sample_size < asc_depth) {
    fprintf(stderr,"Specify a sample size larger than the ascertainment sample"
	    " size with -N option.\n");
    exit(-1);
  }

  bp_length = mb_length*1e6;
}

static double *standard_neutral_spectrum(int n) {
  int i;
  double *fsp, fsp_sum;

  MA(fsp, sizeof(double)*(n+1));
  fsp[0] = fsp[n] = 0.;

  fsp_sum = 0.;
  for(i=1;i<n;i++) {
    fsp[i] = 1./i;
    fsp_sum += fsp[i];
  }

  for(i=1;i<n;i++) fsp[i] /= fsp_sum;

  return fsp;
}

/* the spos_t business is meant to allow alpha to be drawn from a distribution
   rather than the fixed single alpha being used at the moment */
static spos_t *place_sweeps(int n_sweeps, int bp_length, double alpha) {
  int i;
  double sweep_spacing;
  spos_t *spos;

  MA(spos, sizeof(spos_t)*n_sweeps);
  sweep_spacing = bp_length/(double) n_sweeps;

  for(i=0;i<n_sweeps;i++) {
    spos[i].pos = (int) (sweep_spacing*(i+0.5));
    spos[i].alpha = alpha;
  }

  return spos;
}

static spos_t *search_nearest_sweep(int pos, spos_t *spos, int n) {
  int i, j, m;

  i = 0;
  j = n;
  while((j-i)>1) {
    m = (i+j)/2;
    if (pos < spos[m].pos) {
      j = m;
    } else if (pos > spos[m].pos) {
      i = m;
    } else {
      i = j = m;
    }
  }
  if (j<n && pos - spos[i].pos > spos[j].pos - pos) return spos + j;
  return spos + i;
}
    
static int snp_compare(sm_snp_t *a, sm_snp_t *b) {
  return a->pos - b->pos;
}

static sm_snp_t *sample_snps(int n_snps, int bp_length, spos_t *spos,
			  int n_sweeps, sm_ptable_t *sm_p) {
  int i, f;
  double u, log_ad, *fsp, fsp_sum;
  spos_t *nearest_sweep;
  sm_snp_t *snps;

  MA(snps, sizeof(sm_snp_t)*n_snps);
  MA(fsp, sizeof(double)*(sample_size+1));
  for(i=0;i<n_snps;i++) {
    if (i % 100 == 0) 
      fprintf(stderr,"\rsnps sampled: %8d",i);
    snps[i].pos = gsl_rng_uniform(rng)*bp_length;
    if (n_sweeps > 0) {
      nearest_sweep = search_nearest_sweep(snps[i].pos, spos, n_sweeps);
      log_ad = log(nearest_sweep->alpha*abs(nearest_sweep->pos - snps[i].pos));
    } else {
      log_ad = LOG_AD_MAX;
    }

    if (log_ad < LOG_AD_MAX) {
      fsp_sum = 0.;
      for(f=1;f<sample_size;f++) {
	fsp[f] = exp(spline_interpolate(sm_p->spline_func[f],
					log_ad));
	fsp_sum += fsp[f];
      }
    } else {
      fsp_sum = 0.;
      for(f=1;f<sample_size;f++) {
	fsp[f] = sm_p->fsp[f];
	fsp_sum += fsp[f];
      }
    }
      

    f = 0;
    u = gsl_rng_uniform(rng)*fsp_sum;
    while(f<sample_size-1 && u>0.) {
      //      fprintf(stderr,"%d\t%d\t%g\n",i,f,u);
      f++;
      u -= fsp[f];
    }
    snps[i].f = f;
  }
  fprintf(stderr,"\rsnps sampled: %8d - done\n",i);
  qsort(snps, n_snps, sizeof(sm_snp_t), (void *) snp_compare);
  return snps;
}

static double ascprob_subsample(int k, int d, int n) {

  return 1.0 - ((exp(lchoose(k,d)) + exp(lchoose(n-k, d))) /
		exp(lchoose(n,d)));
}

static sm_snp_t *subsample_ascertainment(int *r_nasc, sm_snp_t *snps, 
					 int n_snps, int asc_depth) {
  int i, n_asc;
  double p;
  sm_snp_t *asc_snps;

  n_asc = 0;
  asc_snps = NULL;
  for(i=0;i<n_snps;i++) {
    p = ascprob_subsample(snps[i].f, asc_depth, sample_size);
    if (gsl_rng_uniform(rng) < p) {
      if (n_asc % SNP_ALLOC_STEP == 0)
	RA(asc_snps, sizeof(sm_snp_t)*(n_asc + SNP_ALLOC_STEP));
      asc_snps[n_asc] = snps[i];
      n_asc++;
    }
  }

  *r_nasc = n_asc;
  return asc_snps;
}

static void permute(int *p, int n) {
  int i, j, t;

  for(i=0;i<n;i++) {
    j = gsl_rng_uniform(rng)*n;
    t = p[i];
    p[i] = p[j];
    p[j] = t;
  }
}

static sm_snp_t *random_ascertainment(sm_snp_t *snps, int n_snps, int n_asc) {
  int i, *random_order;
  sm_snp_t *rnd_snps;

  MA(random_order, sizeof(int)*n_snps);
  for(i=0;i<n_snps;i++)
    random_order[i] = i;
  permute(random_order, n_snps);

  MA(rnd_snps, sizeof(sm_snp_t)*n_asc);
  for(i=0;i<n_asc;i++)
    rnd_snps[i] = snps[random_order[i]];

  qsort(rnd_snps, n_asc, sizeof(sm_snp_t), (void *) snp_compare);
  
  free(random_order);
  return rnd_snps;
}

static void output_snps(char *fname, sm_snp_t *snps, int n_snps) {
  int i;
  FILE *f;

  f = fopen(fname, "w");
  if (f == NULL) {
    fprintf(stderr,"Can't open output file \"%s\" (%s)\n", fname,
	    strerror(errno));
    return;
  }

  for(i=0;i<n_snps;i++)
    fprintf(f,"%d\t%d\t%d\t0\n",snps[i].pos, snps[i].f, sample_size);

  fclose(f);
}

int main(int argc, char *argv[]) {
  int n_snps, n_asc;
  double *fsp;
  char *complete_fname, *asc_fname, *rnd_fname;
  sm_snp_t *snps, *asc_snps, *rnd_snps;
  sm_ptable_t sm_p;
  spos_t *sweep_positions;

  init_options();
  cmdline_getoptions(options, argc, argv);
  verify_options();

  fsp = standard_neutral_spectrum(sample_size);
  sm_p = compute_sweep_model_fsp(fsp, sample_size, 0);

  sweep_positions = place_sweeps(n_sweeps, bp_length, alpha);

  n_snps = snp_density*(bp_length/1000.0);
  snps = sample_snps(n_snps, bp_length, sweep_positions, n_sweeps, &sm_p);

  asc_snps = subsample_ascertainment(&n_asc, snps, n_snps, asc_depth);
  rnd_snps = random_ascertainment(snps, n_snps, n_asc);

  (void) asprintf(&complete_fname,"%s-complete.sf", output_basename);
  (void) asprintf(&asc_fname,"%s-asc.sf", output_basename);
  (void) asprintf(&rnd_fname,"%s-rnd.sf", output_basename);

  if (output_complete) 
    output_snps(complete_fname, snps, n_snps);

  output_snps(asc_fname, asc_snps, n_asc);
  output_snps(rnd_fname, rnd_snps, n_asc);

  return 0;
}
