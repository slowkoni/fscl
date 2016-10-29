#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef LINUX
#include <signal.h>
#include <execinfo.h>
#endif

#include <kmacros.h>
#include <cmdline-utils.h>

#include "fscl.h"

static char *snp_fname, *ms_fname, *bs_fname, *output_bs_fname;
static int force_neutral_spectrum, ms_segment_length, ms_folded;
static int ms_sample_first, ms_sample_size;
static int asc_depth, asc_min_freq, ascbias_background_only;
static int include_invariant, maximum_only;
static int small_grid_sp, large_grid_sp, dont_scan;
static int n_threads;
static double permute_nbp, alpha_factor, scan_width_mb;
static int verbosity_level;
char *prepend_label, *output_fname;
int spline_pts, n_permute;

static int eval_range, bp_resl;

static option_t options[] = {
  { 'f', "snpfile", &snp_fname, OPTION_STRINGTYPE, 1, 1, 
    "File name of file with SNP frequency data" },
  { 'd', "asc-depth", &asc_depth, OPTION_INTTYPE, 1, 1,
    "Depth of SNP ascertainment sample" },
  { 0, "asc-minimum-freq", &asc_min_freq, OPT_INT, 1, 1,
    "minimum number of observations of both alleles for SNP ascertainment"},
  { 'p', "n-permute", &n_permute, OPT_INT, 1, 1,
    "number of snp block permutations for p-value computations" },
  { 0, "permute-nbp", &permute_nbp, OPT_DBL, 1, 1,
    "probability for switching to a new snp block for permutations" },
  { 0, "n-threads", &n_threads, OPT_INT, 1, 1,
    "number of simultaneous threads for snp block permutations" },
  { 'a', "alpha-factor", &alpha_factor, OPT_DBL, 1, 1,
    "multiply 1/alpha by this factor to determine single sweep window size" },
  { 'g', "fine-grid-spacing", &small_grid_sp, OPTION_INTTYPE, 1, 1,
    "Spacing of candidate sweep points along the chromosome (in bp)" },
  { 'G', "coarse-grid-spacing", &large_grid_sp, OPTION_INTTYPE, 1, 1,
    "Size of coarse grid in which CLR maxima will be selected" },
  { 'w', "sweep-width", &scan_width_mb, OPT_DBL, 1, 1,
    "maximum width of sweep effect in scanning, in Mb" },

  { 'm', "msfile", &ms_fname, OPTION_STRINGTYPE, 1, 1,
    "Name of an ms output file" },
  { 0, "ms-segment-length", &ms_segment_length, OPTION_INTTYPE, 1, 1,
    "Length in bp of simulated ms segments (use with -m option only)" },
  { 0, "ms-folded", &ms_folded, OPTION_FLAGTYPE, 0, 0,
    "For ms input, treat all sites as folded" },
  { 0, "max-only", &maximum_only, OPTION_FLAGTYPE, 0, 0,
    "for ms input, output only the maximum CLR for each input block" },
  { 0, "ms-sample-first", &ms_sample_first, OPT_INT, 1, 1,
    "index of first chromosome in ms sample to analyze" },
  { 0, "ms-sample-size", &ms_sample_size, OPT_INT, 1, 1,
    "number of consecutive chromosomes in ms output to take as the sample" },

  { 0, "force-neutral-spectrum", &force_neutral_spectrum, OPTION_FLAGTYPE, 
    0, 0,"Do not estimate background spectrum from the data. Use sum(1/i)/i" },
  { 'b', "background-spectrum", &bs_fname, OPTION_STRINGTYPE,
    1, 1, "Load the background frequency spectrum from a file" },
  { 0, "output-bs", &output_bs_fname, OPT_STR, 1, 1,
    "write estimated background site-frequency spectra to file" },

  { 0, "include-invariant", &include_invariant, OPTION_FLAGTYPE, 0, 0,
    "Include invariant sites in analysis (default is to ignore them)"},
  { 0, "splines", &spline_pts, OPTION_INTTYPE, 1, 1,
    "Number of knot points to approximate sweep model w/ respect to alpha" },

  { 0, "prepend-label", &prepend_label, OPT_STR, 1, 1,
    "optional token to prepend to each line of the sweep scan output" },
  { 'v', "verbosity", &verbosity_level, OPT_INT, 1, 1,
    "verbosity level 0-4, default 3, debug 4 and above" },
  { 'o', "output-file", &output_fname, OPT_STR, 1, 1,
    "output file for scan results" },

  { 0, "no-scan", &dont_scan, OPT_FLAG, 0, 0, 
    "do not scan chromosome, compute background frequency spectrum only" },

  { 0, "ascbias-background-only", &ascbias_background_only, OPT_FLAG, 0, 0,
    "correct for ascertainment bias only in estimating the background site "
    "frequency spectrum" },

  { 0, NULL, NULL, 0, 0, 0, NULL }
};

gsl_rng *rng;

#ifdef LINUX
static void show_backtrace(int signal, struct sigcontext ctx) __attribute__((unused));
static void show_backtrace(int signal, struct sigcontext ctx) {
  void *bt_symbols[16];
  int n_symbols, i;
  char **strings;

  fprintf(stderr,"Segmentation fault -- stack backtrace follows:\n");

  n_symbols = backtrace(bt_symbols, 16);
  //  bt_symbols[1] = (void *) ctx.rip;

  strings = backtrace_symbols(bt_symbols, n_symbols);
  for(i=1;i<n_symbols;i++) {
    fprintf(stderr,"frame %d: %s\n", i, strings[i]);
  }

  exit(11);
}  
#endif

static void init_options() {

#ifdef LINUX
  //  signal(SIGSEGV, (void *) show_backtrace);
#endif

  //  srand(getpid() + time(NULL));
#warning using hard-coded fixed random seed
  srand(0xFD821A6);
  rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, 0xFD821A6);

  verbosity_level = MSG_STATUS;

  dont_scan = 0;
  snp_fname = NULL;
  asc_depth = 0;
  ascbias_background_only = 0;
  asc_min_freq = 1;
  n_permute = 0;
  permute_nbp = 0.1;
  n_threads = 1;

  ms_fname = NULL;
  ms_segment_length = 0;
  ms_folded = 0;
  maximum_only = 0;
  ms_sample_first = 0;
  ms_sample_size = 0;

  small_grid_sp = 1000;
  large_grid_sp = 100000;
  scan_width_mb = 1.0;

  force_neutral_spectrum = 0;
  bs_fname = NULL;
  output_fname = NULL;
  output_bs_fname = NULL;

  spline_pts = 200;
  include_invariant = 0;

  alpha_factor = 1.0;

  prepend_label = NULL;

  bp_resl = 128;
  eval_range = 8192;

  init_log_table();
}

static void validate_options() {
  int stop;

  if (verbosity_level < 0) verbosity_level = 0;
  configure_logmsg(verbosity_level);

  if (ascbias_background_only) logmsg(MSG_STATUS,"ascertainment bias correction to background site frequency spectrum only\n");

  stop = 0;
  if (spline_pts < 200) {
    logmsg(MSG_ERROR,"Error: must use at least 200 spline functions to "
	   "approximate sweep model\nlikelihood function.\n");
    stop = 1;
  }
  if (spline_pts > 500) {
    logmsg(MSG_WARN,"Warning: estimating %d spline functions will "
	   "significantly inflate execution\ntime with little gain "
	   "in accuracy...\n", spline_pts);
    if (spline_pts > 1000) {
      logmsg(MSG_WARN,"Warning: estimating %d spline functions will require "
	     "%3.1f Mb\nduring estimation...\n", spline_pts, 
	     (128.0*spline_pts*spline_pts)/(1000000.0));
    }
  }
  if (snp_fname == NULL && ms_fname == NULL) {
    logmsg(MSG_ERROR,"Error: input snp frequency file or ms file not "
	   "specified. Use -f option or -m option.\n");
    stop = 1;
  }
  if (snp_fname && ms_fname) {
    logmsg(MSG_ERROR,"Specify either a snp frequency file or an ms file, not "
	   "both.\n");
    stop = 1;
  }
  if (ms_segment_length && ms_fname==NULL) {
    logmsg(MSG_WARN,"Warning: --ms-segment-length option ignored if -m "
	   "option is not used.\n");
    ms_segment_length = 0;
  }

  if (asc_depth == 1 || asc_depth < 0) {
    logmsg(MSG_ERROR,"Error: if specified, ascertainment sample depth must be "
	   "at least 2.\n");
    stop = 1;
  }

  if (asc_depth >= 2 && asc_min_freq > 2*asc_depth) {
    logmsg(MSG_ERROR,"Error: SNP ascertainment is impossible with asc. sample"
	   " depth and asc. minimum allele frequency setting\n");
    stop=1;
  }
  if (asc_depth >= 2 && asc_min_freq == 0) asc_min_freq = 1;

  if (small_grid_sp < 1 && output_bs_fname == NULL) {
    logmsg(MSG_ERROR,"Error: specify sweep position grid spacing with -g "
	   "option (in bp).\n");
    stop = 1;
  }

  if (output_bs_fname == NULL && large_grid_sp % small_grid_sp != 0) {
    logmsg(MSG_ERROR,"Error: fine grid spacing must evenly divide coarse grid "
	   "spacing.\n");
    stop = 1;
  }

  if (stop) {
    if (verbosity_level <= MSG_FATAL)
      logmsg(MSG_FATAL, "Fatal errors have occurred, use -v 1 or greater to "
	     "see them.\n");
    exit(-1);
  }

}

static void scan_free(scan_t *scan_obj) {
  int i;

  if (scan_obj->scan_pts) free(scan_obj->scan_pts);
  if (scan_obj->sample_depths) free(scan_obj->sample_depths);
  if (scan_obj->snps) free(scan_obj->snps);
  for(i=0;i<scan_obj->n_chromosomes;i++)
    free(scan_obj->chr_limits[i].name);
  free(scan_obj->chr_limits);
  free(scan_obj);
}

int main(int argc, char *argv[]) {
  double **fsp;
  scan_t *scan_obj;
  sm_ptable_t *sm_p;

  init_options();
  cmdline_getoptions(options, argc, argv);
  validate_options();

  if (ms_fname) {
    scan_obj = ms_background(ms_fname, ms_segment_length, ms_folded,
			     ms_sample_first, ms_sample_size);
    fsp = background_fsp(scan_obj, force_neutral_spectrum, bs_fname,
			 include_invariant);

    if (output_bs_fname) output_background_fs(output_bs_fname, scan_obj, fsp);

    if (dont_scan == 0) {
      sm_p = compute_sweep_model_tables(scan_obj, fsp, asc_depth, asc_min_freq,
					ascbias_background_only,
					include_invariant);
      scan_free(scan_obj);

      ms_openfile(ms_fname);
      while((scan_obj = ms_next_block(ms_segment_length, ms_folded,
				      ms_sample_first, ms_sample_size))) {
	if (scan_obj->n_snps > 0) {
	  compute_snp_null_model(scan_obj, fsp);

	  scan_chromosome(scan_obj, sm_p, eval_range, bp_resl, large_grid_sp, 
			  n_threads);
	  if (n_permute > 0)
	    scan_permute(scan_obj, sm_p, n_permute, permute_nbp, alpha_factor, 
			 n_threads, eval_range, bp_resl, large_grid_sp);

	  scan_output(output_fname, scan_obj, maximum_only, n_permute, 
		      prepend_label);
	}
	scan_free(scan_obj);
      }
    } else {
      scan_free(scan_obj);
    }
  } else {
    scan_obj = load_snp_input(snp_fname, include_invariant);

    fsp = background_fsp(scan_obj, force_neutral_spectrum, bs_fname,
			 include_invariant);
    if (output_bs_fname) output_background_fs(output_bs_fname, scan_obj, fsp);
    
    if (dont_scan == 0) {
      sm_p = compute_sweep_model_tables(scan_obj, fsp, asc_depth, asc_min_freq,
					ascbias_background_only,
					include_invariant);
      compute_snp_null_model(scan_obj, fsp);

      scan_chromosome(scan_obj, sm_p, eval_range, bp_resl, large_grid_sp, 
		      n_threads);
      if (n_permute > 0) 
	scan_permute(scan_obj, sm_p, n_permute, permute_nbp, alpha_factor, 
		     n_threads, eval_range, bp_resl, large_grid_sp);

      scan_output(output_fname, scan_obj, maximum_only, n_permute, 
		  prepend_label);
    }
    scan_free(scan_obj);
  }

  return 0;
}

