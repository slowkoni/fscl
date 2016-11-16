#ifndef FSCL_H
#define FSCL_H

#include <stdint.h>

#if 1
typedef struct {
  int chr;
  int pos;
  double null_logl;
  int obs_freq;
  int depth_p;
  int folded; 
} snp_t;

#else
typedef struct {
  int pos;
  float null_logl;
  uint16_t obs_freq;
  uint8_t depth_p;
  uint8_t chr:7;
  uint8_t folded:1;
} snp_t;
#endif
typedef struct {
  int chr;
  char *name;
  int start_index;
  int n_snps;
  int start_pos;
  int bp_length;
} chr_limits_t;

typedef struct {
  int chr;
  int nearest_snp;
  int sweep_pos;
  int n_snps;
  int window_start;
  int window_end;
  double lalpha;
  double null_logl;
  double sm_logl;
  double clr;
  int permute_n;
  int permute_p;
  int permute_finished;
  int scan_running;
  float *permute_clr;
} scan_pt_t;

typedef struct {
  int n_snps;
  snp_t *snps;
  int n_depths;
  int *sample_depths;
  int n_scan_pts;
  scan_pt_t *scan_pts;
  chr_limits_t *chr_limits;
  int n_chromosomes;
} scan_t;

typedef struct {
  int n;
  double *knot_points;
  double **coef;
} spline_t;

typedef struct {
  spline_t **spline_func;
  spline_t **fspline_func;
  int sample_size;
  double **pbk;
  double *fsp;
} sm_ptable_t;

/* controls the range of the alpha*distance spline interpolation */
#define LOG_AD_MIN (-20.0)
#define LOG_AD_MAX (4.0)

#define N_SPLINE_KNOTS (200)

#define SNP_ALLOC_STEP (512*1024)

/* snp-input.c */
scan_t *load_snp_input(char *snp_fname, int include_invariant,
		       int minimum_obs_depth);

/* background-fsp.c */
double **background_fsp(scan_t *scan_obj, int force_neutral_spectrum,
			char *background_fsfname, int include_invariant);
void output_background_fs(char *fname, scan_t *scan_obj, double **fsp);
double lchoose(int n, int k);

/* sm-spline.c */
sm_ptable_t *compute_sweep_model_tables(scan_t *scan_obj, double **fsp,
					int asc_depth, int asc_min_freq,
					int ascbias_background_only,
					int include_invariant);
double spline_interpolate(spline_t *spf, double x);

/* sm-search.h */
void init_log_table(void);
void search_maxalpha(scan_pt_t *scan_pt, snp_t *snps, sm_ptable_t *sm_p);

/* scan-chromosome.c */
void compute_snp_null_model(scan_t *scan_obj, double **fsp);
void scan_chromosome(scan_t *scan_obj, sm_ptable_t *sm_p, int eval_range,
		     int bp_resl, int large_grid_sp, int n_threads);
void scan_permute(scan_t *scan_obj, sm_ptable_t *sm_p,
		  int n_permute, double permute_nbp, double alpha_factor, 
		  int n_threads, int eval_range, int bp_resl, 
		  int large_grid_sp, double scan_width_mb);
void scan_output(char *output_fname, scan_t *scan_obj, int maximum_only, 
		 int n_permute, char *prepend_label);

/* ms-input.c */
void ms_openfile(char *ms_fname);
scan_t *ms_background(char *ms_fname, int ms_segment_length, int ms_folded,
		      int ms_sample_first, int ms_sample_last);
scan_t *ms_next_block(int ms_segment_length, int ms_folded,
		      int ms_sample_first, int ms_sample_last);

/* asc-bias.c */
double *ascbias_adjust_background(double *bsf, int n, int asc_depth, 
				  int min_obs);
void ascbias_adjust_expect(double *fsp, int n, int min_obs, int d);

/* logmsg.c */
#define MAX_LOGMSG_LENGTH (4096)
enum { MSG_FATAL = 0, MSG_ERROR, MSG_WARN, MSG_STATUS, MSG_DEBUG1, MSG_DEBUG2 };

void configure_logmsg(int level);
void logmsg(int priority, volatile char *s, ...);
void cr_logmsg(int priority, volatile char *s, ...);

#endif
