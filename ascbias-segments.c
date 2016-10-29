#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <kmacros.h>
#include <cmdline-utils.h>
#include <msparser.h>

static int mb_length, sample_size, asc_depth, n_snps, n_repl, double_hit;
static double rho_Mb, rho_std;
static char *ms_opts;

static option_t options[] = {
  { 'r', "rho-Mb", &rho_Mb, OPT_DBL, 1, 1, 
    "total segment rho" },
  { 0, "rho-stdv", &rho_std, OPT_DBL, 1, 1,
    "standard deviation of rho/Mb" },
  { 'l', "mb-length", &mb_length, OPT_INT, 1, 1,
    "length of chromosome in Mb" },
  { 'n', "sample-size", &sample_size, OPT_INT, 1, 1,
    "number of chromosomes in sample" },
  { 'd', "asc-depth", &asc_depth, OPT_INT, 1, 1,
    "depth of ascertainment sample" },
  { 's', "n-snps", &n_snps, OPT_INT, 1, 1,
    "total number of SNPs on chromosome" },
  { 'N', "n-repl", &n_repl, OPT_INT, 1, 1,
    "number of replicate samples to draw" },
  { 0, "double-hit", &double_hit, OPT_FLAG, 1, 1,
    "apply double-hit rule for ascertainment" },
  { 0, "ms-opts", &ms_opts, OPT_STR, 1, 1,
    "demographic model options for ms" },
  { 0, NULL, NULL, 0, 0, 0, NULL }
};

static gsl_rng *rng;

typedef struct {
  int freq;
  int pos;
  char *alleles;
} snp_t;

static void init_options(void) {
  rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, time(NULL));

  rho_Mb = 100.;
  rho_std = 0.;
  mb_length = 25;
  sample_size = 128;
  asc_depth = 0;
  n_snps = 4000;
  n_repl = 100;
  double_hit = 0;
  ms_opts = "";
}

static void verify_options(void) {

  if (rho_Mb < 0. || mb_length < 1 || sample_size < asc_depth || 
      sample_size < 2 || (asc_depth < 2 && asc_depth != 0) || n_snps < 2 ||
      n_repl < 1) {
    fprintf(stderr,"invalid settings for command line options used.\n");
    exit(-1);
  }
}

static int extract_snps(snp_t *snps, msblock_t *msb, int s_pos, 
			int segment_length, int sample_size, int asc_depth,
			int double_hit) {
  int i, j, k, d, ascertained;
  char *alleles;

  alleles = alloca(sizeof(char)*sample_size);
  k = 0;
  for(i=0;i<msb->n_poly;i++) {
    d = 0;
    ascertained = 0;
    for(j=0;j<asc_depth;j++) {
      if (msb->haplotypes[j][i] == '1') d++;
    }
    if (asc_depth == 0 || (double_hit == 0 && d > 0 && d < asc_depth) ||
	(double_hit == 1 && d > 1 && d < asc_depth - 1)) {
      ascertained = 1;
    }
    d = 0;
    for(j=0;j<msb->n_haplotypes;j++) {
      if (msb->haplotypes[j][i] == '1') d++;
      alleles[j] = msb->haplotypes[j][i];
    }
    if (ascertained) {
      snps[k].freq = d;
      snps[k].pos = s_pos + msb->positions[i]*segment_length;
      strncpy(snps[k].alleles, alleles, sample_size);
      k++;
    }
  }

  return k;
}

static void permute_snps(snp_t *snps, int n) {
  int i, j;
  snp_t t;

  for(i=0;i<n;i++) {
    j = (rand()/(RAND_MAX + 1.0))*n;
    t = snps[i];
    snps[i] = snps[j];
    snps[j] = t;
  }
}

static int snp_compare(snp_t *a, snp_t *b) {
  return a->pos - b->pos;
}

static void ms_output(snp_t *snps, int n_snps, int sample_size, 
		      double bp_length) {
  int i, j;

  fprintf(stdout,"\n//\n");
  fprintf(stdout,"segsites: %d\n",n_snps);
  fprintf(stdout,"positions:");
  for(i=0;i<n_snps;i++) {
    fprintf(stdout," %1.8e",snps[i].pos/bp_length);
  }
  fprintf(stdout,"\n");

  for(j=0;j<sample_size;j++) {
    for(i=0;i<n_snps;i++)
      putc(snps[i].alleles[j], stdout);
    putc('\n', stdout);
  }

}

static double elapsed_seconds(struct timeval then) {
  struct timeval now;

  gettimeofday(&now, NULL);
  return (now.tv_sec - then.tv_sec) + (now.tv_usec - then.tv_usec)/1e6;
}
  

static void draw_samples(int n_repl, int sample_size, int mb_length,
			 int n_snps, int asc_depth, int double_hit,
			 double rho_Mb, double rho_std, char *ms_opts) {
  int i, j, b, k, status;
  int ms_segments, ms_snps;
  double ms_length, mean_rho, stdv_rho, segment_rho, asc_factor;
  char *alleles, ms_cmd[1024];
  msblock_t *msb;
  snp_t *snps;
  struct timeval stopwatch;
  FILE *msf;

  ms_segments = mb_length/1;
  ms_length = mb_length / ms_segments;
  mean_rho = rho_Mb*ms_length;
  stdv_rho = rho_std*ms_length;

  fprintf(stdout,"ms %d %d -s %d -r %1.1f %d %s\n", sample_size, n_repl,
	  n_snps, rho_Mb*mb_length, (int) (mb_length*1e6), ms_opts);
  fprintf(stdout,"%d %d %d\n",rand(),rand(),rand());

  if (asc_depth > 0) {
    asc_factor = 1.0 + 12.0/asc_depth;
  } else {
    asc_factor = 1.0;
  }
  ms_snps = 0;
  snps = NULL;
  alleles = NULL;

  gettimeofday(&stopwatch, NULL);
  i = 0;
  while(i<n_repl) {
    if ((int) (n_snps / ms_segments * asc_factor + 1) != ms_snps) {
      ms_snps = n_snps / ms_segments * asc_factor + 1;
      RA(snps, sizeof(snp_t)*(ms_snps*(ms_segments+1)));
      RA(alleles, sizeof(char)*(ms_snps*(ms_segments+1)*sample_size));
      snps[0].alleles = alleles;
      for(j=1;j<ms_snps*ms_segments;j++)
	snps[j].alleles = snps[j-1].alleles + sample_size;
    }

    k = 0;
    for(b=0;b<ms_segments;b++) {
      //      segment_rho = mean_rho + gsl_ran_gaussian(rng, stdv_rho);
      //      if (segment_rho < 0.) segment_rho = 0.;
      segment_rho = (mean_rho - stdv_rho) + (2*stdv_rho)*b/(double)ms_segments;

      sprintf(ms_cmd,"ms %d 1 -s %d -r %1.1f %d %s", sample_size, 
	      ms_snps, segment_rho, (int) (ms_length*1e6), ms_opts);
      msf = msparser_execute(ms_cmd);

      msb = msparser_block();
      if (msb == NULL) {
	fprintf(stderr,"ms execution failed: %s\n", ms_cmd);
	exit(-1);
      }
      
      fprintf(stderr,"\rrepl %3d, segment %3d, snp %7d", i, b, k);
      k += extract_snps(snps + k, msb, b*ms_length*1e6, ms_length*1e6,
			  sample_size, asc_depth, double_hit);
      msparser_block_free(msb);

      fclose(msf);
      wait(&status);
    }
    fprintf(stderr,"\t%1.1f sec/repl\n", elapsed_seconds(stopwatch)/(i+1));


    if (k < n_snps) {
      asc_factor *= 2.0;
      continue;
    }
    if (k > n_snps*2.0) asc_factor = asc_factor * 0.67;
    permute_snps(snps, k);

    if (n_snps < k) k = n_snps;
    qsort(snps, k, sizeof(snp_t), (void *) snp_compare);

    ms_output(snps, k, sample_size, mb_length*1e6);
    i++;

  }

  free(alleles);
  free(snps);
}


int main(int argc, char *argv[]) {

  init_options();
  cmdline_getoptions(options, argc, argv);
  verify_options();

  draw_samples(n_repl, sample_size, mb_length, n_snps, asc_depth, double_hit,
	       rho_Mb, rho_std, ms_opts);

  return 0;
}
