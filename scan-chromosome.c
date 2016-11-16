#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>
#include <signal.h>
#include <math.h>
#include <float.h>
#include <pthread.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>

#include <kmacros.h>
#include "fscl.h"

extern gsl_rng *rng;

void compute_snp_null_model(scan_t *scan_obj, double **fsp) {
  snp_t *snps;
  int i, depth;

  snps = scan_obj->snps;
  for(i=0;i<scan_obj->n_snps;i++) {
    depth = scan_obj->sample_depths[snps[i].depth_p];
    if (snps[i].folded && snps[i].obs_freq != depth - snps[i].obs_freq) {
      snps[i].null_logl = log(fsp[snps[i].depth_p][snps[i].obs_freq] + 
			      fsp[snps[i].depth_p][depth - snps[i].obs_freq]);
    } else {
      snps[i].null_logl = log(fsp[snps[i].depth_p][snps[i].obs_freq]);
    }
  }
}

static int search_snppos(snp_t *snps, int n_snps, double sweep_pos) {
  int i, j, m;

  i = 0;
  j = n_snps;
  while(j-i>1) {
    m = (i+j)/2;
    if (snps[m].pos < sweep_pos) {
      i = m;
    } else {
      j = m;
    }
  }
  if (j == n_snps) return n_snps - 1;
  if ((sweep_pos - snps[i].pos) < (snps[j].pos - sweep_pos)) 
    return i;
  return j;
}

static void init_scan_result(scan_pt_t *scan_pt, int chr, snp_t *snps, 
			     chr_limits_t *limits, int eval_range, int pos) {
  int i, chm_start, chm_stop;

  scan_pt->chr = chr;
  scan_pt->nearest_snp = limits->start_index + 
    search_snppos(snps + limits->start_index, limits->n_snps, pos);

  i = scan_pt->nearest_snp;
  while(i < limits->n_snps && snps[i].pos == pos) {
    i++;
    pos++;
  }
  scan_pt->sweep_pos = pos;
  
  chm_start = limits->start_index;
  chm_stop  = limits->start_index + limits->n_snps - 1;

  if (scan_pt->nearest_snp - eval_range < chm_start) {
    scan_pt->window_start = chm_start;
    scan_pt->window_end   = chm_start + eval_range*2;
    if (scan_pt->window_end   >  chm_stop) scan_pt->window_end    = chm_stop;
  } 
  else if (scan_pt->nearest_snp + eval_range > chm_stop) {
    scan_pt->window_end   = chm_stop;
    scan_pt->window_start = chm_stop - eval_range * 2;
    if (scan_pt->window_start <  chm_start) scan_pt->window_start = chm_start;
  } 
  else {
    scan_pt->window_start = scan_pt->nearest_snp - eval_range;
    scan_pt->window_end   = scan_pt->nearest_snp + eval_range;
  }

  scan_pt->n_snps = scan_pt->window_end - scan_pt->window_start + 1;
  scan_pt->null_logl = 0.;
  for(i=scan_pt->window_start;i<=scan_pt->window_end;i++)
    scan_pt->null_logl += snps[i].null_logl;

  scan_pt->sm_logl = -DBL_MAX;
  scan_pt->lalpha = LOG_AD_MAX;
  scan_pt->permute_n = 0;
  scan_pt->permute_p = 0;
  scan_pt->permute_finished = 0;
}

static scan_pt_t
search_maxpos_recursive(scan_pt_t start_pt, scan_pt_t end_pt,
		        snp_t *snps, chr_limits_t *limits,
			int eval_range, int bp_resl, sm_ptable_t *sm_p) {
  scan_pt_t mid_pt;

  if (end_pt.sweep_pos - start_pt.sweep_pos <= bp_resl) 
    return start_pt.clr > end_pt.clr ? start_pt : end_pt;

  init_scan_result(&mid_pt, start_pt.chr, snps, limits, eval_range,
		   (start_pt.sweep_pos + end_pt.sweep_pos)/2);
  search_maxalpha(&mid_pt, snps, sm_p);

  if ((start_pt.clr + mid_pt.clr) >= (end_pt.clr + mid_pt.clr)) {
    return search_maxpos_recursive(start_pt, mid_pt, snps, limits, eval_range,
				   bp_resl, sm_p);
  } else {
    return search_maxpos_recursive(mid_pt, end_pt, snps, limits, eval_range,
				   bp_resl, sm_p);
  }
}


static scan_pt_t search_maxpos(int chr, int start_pos, int end_pos,
			       snp_t *snps, chr_limits_t *limits,
			       int eval_range, int bp_resl, sm_ptable_t *sm_p) {
  scan_pt_t start_pt, end_pt;

  init_scan_result(&start_pt, chr, snps, limits, eval_range, start_pos);
  search_maxalpha (&start_pt, snps, sm_p);

  init_scan_result(&end_pt,   chr, snps, limits, eval_range, end_pos);
  search_maxalpha (&end_pt,   snps, sm_p);

  return search_maxpos_recursive(start_pt, end_pt, snps, limits, eval_range, 
				 bp_resl, sm_p);
}


/*static double elapsed_time(struct timeval then) {
  struct timeval now;

  gettimeofday(&now, NULL);
  return (now.tv_sec - then.tv_sec) + (now.tv_usec - then.tv_usec)/1e6;
  }*/

typedef struct {
  scan_t *scan_obj;
  sm_ptable_t *sm_ptable;

  int large_grid_sp;
  int eval_range;
  int bp_resl;

  int scan_chm;
  int scan_pos;
  pthread_mutex_t scan_lock;
} scan_args_t;

static void *scan_thread(scan_args_t *args) {
  int last_chm, chm, start_pos, end_pos;
  scan_pt_t max_pt;

  scan_t *scan_obj;
  snp_t *snps;
  chr_limits_t *limits;

  scan_obj = args->scan_obj;
  snps = scan_obj->snps;

  last_chm = 0;
  limits = scan_obj->chr_limits;

  pthread_mutex_lock(&args->scan_lock);
  for(;;) {
    if (args->scan_chm == scan_obj->n_chromosomes) break;

    if (args->scan_pos >= scan_obj->chr_limits[args->scan_chm].bp_length) {
      args->scan_chm++;

      if (args->scan_chm == scan_obj->n_chromosomes) break;
      args->scan_pos = scan_obj->chr_limits[args->scan_chm].start_pos;
    }

    start_pos = args->scan_pos;
    chm = args->scan_chm;

    cr_logmsg(MSG_STATUS, "Scanning chromosome %s - %5.2f Mb", 
	      scan_obj->chr_limits[chm].name, args->scan_pos/1e6);

    args->scan_pos += args->large_grid_sp;
    pthread_mutex_unlock(&args->scan_lock);

    if (chm != last_chm) {
      limits = scan_obj->chr_limits + chm;
      last_chm = chm;
    }

    end_pos = start_pos + args->large_grid_sp;
    if (end_pos > scan_obj->chr_limits[chm].bp_length) 
      end_pos = scan_obj->chr_limits[chm].bp_length;

    max_pt = search_maxpos(chm, start_pos, end_pos, snps, limits,
			   args->eval_range, args->bp_resl, args->sm_ptable);

    pthread_mutex_lock(&args->scan_lock);

    max_pt.permute_clr = scan_obj->scan_pts[scan_obj->n_scan_pts].permute_clr;
    scan_obj->scan_pts[scan_obj->n_scan_pts++] = max_pt;
  }
  pthread_mutex_unlock(&args->scan_lock);

  return NULL;
}

static int scan_pt_pos_compare(scan_pt_t *a, scan_pt_t *b) {

  if (a->chr < b->chr) return -1;
  if (a->chr > b->chr) return  1;
  if (a->sweep_pos < b->sweep_pos) return -1;
  if (a->sweep_pos > b->sweep_pos) return  1;
  return 0;
}

#define CLR_NULL_DIST_SAVE (10000)
void scan_chromosome(scan_t *scan_obj, sm_ptable_t *sm_ptable, int eval_range,
		     int bp_resl, int large_grid_sp, int n_threads) {
  int i;
  pthread_t *threads;
  scan_args_t scan_args;

  scan_obj->n_scan_pts = 0;
  for(i=0;i<scan_obj->n_chromosomes;i++)
    scan_obj->n_scan_pts += scan_obj->chr_limits[i].bp_length/large_grid_sp + 1;
  MA(scan_obj->scan_pts, sizeof(scan_pt_t)*scan_obj->n_scan_pts);

  for(i=0;i<scan_obj->n_scan_pts;i++) {
    MA(scan_obj->scan_pts[i].permute_clr, sizeof(float)*CLR_NULL_DIST_SAVE);
  }
  
  scan_obj->n_scan_pts = 0;

  scan_args.scan_obj = scan_obj;
  scan_args.sm_ptable = sm_ptable;

  scan_args.large_grid_sp = large_grid_sp;
  scan_args.eval_range = eval_range;
  scan_args.bp_resl = bp_resl;

  scan_args.scan_chm = 0;
  scan_args.scan_pos = scan_obj->chr_limits[scan_args.scan_chm].start_pos;
  pthread_mutex_init(&scan_args.scan_lock, NULL);

  MA(threads, sizeof(pthread_t)*n_threads);
  for(i=0;i<n_threads;i++)
    pthread_create(threads + i, NULL, (void *) scan_thread, &scan_args);

  for(i=0;i<n_threads;i++) pthread_join(threads[i], NULL);
  pthread_mutex_destroy(&scan_args.scan_lock);
  free(threads);

  qsort(scan_obj->scan_pts, scan_obj->n_scan_pts, sizeof(scan_pt_t),
	(void *) scan_pt_pos_compare);

  logmsg(MSG_STATUS,"\nInitial scan finished.\n");

#if 0
  for(i=0;i<scan_obj->n_chromosomes;i++) {
    window_start_pos = 0;
    while(window_start_pos < scan_obj->chr_limits[i].bp_length) {
      cr_logmsg(MSG_STATUS, "Scanning chromosome %s - %5.2f Mb", 
		scan_obj->chr_limits[i].name, window_start_pos/1e6);
      pos = window_start_pos;
      init_scan_result(&max_pt, i, scan_obj->snps, 
		       scan_obj->chr_limits[i].start_index,
		       scan_obj->chr_limits[i].n_snps, pos);
      search_maxalpha(&max_pt, scan_obj->snps, scan_obj->n_snps, sm_ptable);
      pos += small_grid_sp;
      while(pos < scan_obj->chr_limits[i].bp_length && 
	    pos - window_start_pos < large_grid_sp) {
	init_scan_result(&scan_pt, i, scan_obj->snps, 
			 scan_obj->chr_limits[i].start_index,
			 scan_obj->chr_limits[i].n_snps, pos);
	search_maxalpha(&scan_pt, scan_obj->snps, scan_obj->n_snps, sm_ptable);
	if (scan_pt.clr > max_pt.clr) {
	  max_pt = scan_pt;
	}
	pos += small_grid_sp;
      }
      scan_obj->scan_pts[scan_obj->n_scan_pts++] = max_pt;
      window_start_pos = pos;
    }
  }
#endif	
#if 0
  scan_obj->n_scan_pts = 0;
  for(i=0;i<scan_obj->n_chromosomes;i++)
    scan_obj->n_scan_pts += scan_obj->chr_limits[i].bp_length/grid_sp + 1;

  MA(scan_obj->scan_pts, sizeof(scan_pt_t)*scan_obj->n_scan_pts);
  k = 0;
  for(i=0;i<scan_obj->n_chromosomes;i++) {
    j = 0;
    while(j<=scan_obj->chr_limits[i].bp_length) {
      init_scan_result(scan_obj->scan_pts + k, i, scan_obj->snps,  
		       scan_obj->chr_limits[i].start_index,
		       scan_obj->chr_limits[i].n_snps, j);
      k++;
      j += grid_sp;
    }
  }

  gettimeofday(&stop_watch, NULL);

  for(i=0;i<scan_obj->n_scan_pts;i++) {
    if (i % 100 == 0)
      cr_logmsg(MSG_STATUS,"Scanning chromosomes - %5.1f%%", 
		i/(double) scan_obj->n_scan_pts * 100.);
    search_maxalpha(scan_obj->scan_pts + i, scan_obj->snps, scan_obj->n_snps,
		    sm_ptable);
  }

  cr_logmsg(MSG_STATUS,"Scanning finished."); 
  logmsg(MSG_STATUS,"");
#endif

}

//#define LN_DEBUG
#ifdef LN_DEBUG
static FILE *rfile;
#endif

static void snp_block_permute(snp_t *p_snps, snp_t *snps, int n_snps, 
			      double permute_nbp, double scan_width_mb) {
  int i, j, k, __attribute__((unused))*already_permuted;
  snp_t tmp;

  memcpy(p_snps, snps, sizeof(snp_t)*n_snps);
  //  for(i=0;i<n_snps;i++) p_snps[i] = snps[i];

  //  CA(already_permuted, sizeof(int)*n_snps);

  i = 0;
  while(i<n_snps) {
    //    while(i<n_snps && already_permuted[i]) i++;
    //if (i == n_snps) break;
    j = rand()/(RAND_MAX + 1.0) * n_snps;
    k = j + (int) (-1.0/permute_nbp * log(rand()/(RAND_MAX + 1.0)));
    //k = 1 + gsl_ran_lognormal(rng, -0.2027, 1.3386);
#ifdef LN_DEBUG
    pthread_mutex_lock(&nd_update_lock);
    fprintf(rfile, "%d\n", k);
    pthread_mutex_unlock(&nd_update_lock);
#endif
    while(k < n_snps && snps[k].chr == snps[j].chr && snps[k].pos - snps[j].pos < scan_width_mb*1e6) k++;
    if (i + (k - j) >= n_snps) k = n_snps;
    if (k > n_snps) {
      j = n_snps - k;
      k = n_snps;
    }

    while(j<k && i < n_snps && j < n_snps) {
      tmp.obs_freq = p_snps[i].obs_freq;
      tmp.depth_p = p_snps[i].depth_p;
      tmp.folded = p_snps[i].folded;
      tmp.null_logl = p_snps[i].null_logl;

      p_snps[i].obs_freq = p_snps[j].obs_freq;
      p_snps[i].depth_p = p_snps[j].depth_p;
      p_snps[i].folded = p_snps[j].folded;
      p_snps[i].null_logl = p_snps[j].null_logl;

      p_snps[j].obs_freq = tmp.obs_freq;
      p_snps[j].depth_p = tmp.depth_p;
      p_snps[j].folded = tmp.folded;
      p_snps[j].null_logl = tmp.null_logl;

      //already_permuted[j] = 1;
      i++;
      j++;
    } 
  }

  //  for(i=0;i<n_snps;i++) p_snps[i].pos = snps[i].pos;
  //free(already_permuted);
}

typedef struct {
  scan_t *scan_obj;
  sm_ptable_t *sm_p;

  int n_permute;
  double permute_nbp;

  snp_t *p_snps;

  int eval_range;
  int bp_resl;
  int large_grid_sp;
  double scan_width_mb;
} thread_args_t;

static int global_permute, global_scan_pt, n_remaining, n_global_scan_pts;
static int *global_scan_pts;
static int n_remaining;
static pthread_mutex_t scan_pt_lock;
static pthread_barrier_t permute_barrier, run_barrier;

static void *scan_permute_thread(thread_args_t *args) {
  int i, j, k, n_snps, my_scan_pt, large_grid_sp;
  int __attribute__((unused)) small_grid_sp;
  int __attribute__((unused)) my_permute_index;
  scan_pt_t max_pt;

  /* actual arguments to this function (packed into args) because we
     can only pass a single (void *) to a thread starting function */
  scan_t *scan_obj;
  sm_ptable_t *sm_p;
  int n_permute;
  double permute_nbp, scan_width_mb;
  snp_t *p_snps;
  
  /* unpack full set of arguments from the structure used to pass through
     pthread_create() */
  scan_obj = args->scan_obj;
  sm_p = args->sm_p;

  n_permute = args->n_permute;
  permute_nbp = args->permute_nbp;
  scan_width_mb = args->scan_width_mb;

  p_snps = args->p_snps;
  large_grid_sp = args->large_grid_sp;
  n_snps = scan_obj->n_snps;
  
  my_permute_index = 0;
  usleep((int) (rand()/(RAND_MAX+1.0)*10000));
  for(;;) {
    if (pthread_barrier_wait(&permute_barrier)==PTHREAD_BARRIER_SERIAL_THREAD) {
      snp_block_permute(p_snps, scan_obj->snps, n_snps, permute_nbp, scan_width_mb);
      global_permute++;

      k = i = 0;
      while(i<n_global_scan_pts) {
	j = i;
	while(j<n_global_scan_pts &&
	      scan_obj->scan_pts[global_scan_pts[j]].permute_finished == 1) j++;
	if (j<n_global_scan_pts)
	  global_scan_pts[k++] = global_scan_pts[j];
	i = j + 1;
      }

      n_global_scan_pts = k;
      global_scan_pt = 0;

      cr_logmsg(MSG_STATUS,"Scanning snp block permutations... %7d (%d "
		"scan pts remaining)        ", global_permute, n_global_scan_pts);
    }
    pthread_barrier_wait(&permute_barrier);

    pthread_mutex_lock(&scan_pt_lock);
    if (n_global_scan_pts == 0) break;
    if (global_permute > n_permute) break;
    pthread_mutex_unlock(&scan_pt_lock);

    for(;;) {
      int chr, start_pos, __attribute__((unused))current_pos;

      pthread_mutex_lock(&scan_pt_lock);
      if (global_scan_pt == n_global_scan_pts) {
	pthread_mutex_unlock(&scan_pt_lock);
	break;
      }
      my_scan_pt = global_scan_pts[global_scan_pt++];
      pthread_mutex_unlock(&scan_pt_lock);

      chr = scan_obj->scan_pts[my_scan_pt].chr;
      start_pos = scan_obj->scan_pts[my_scan_pt].sweep_pos - 
	(scan_obj->scan_pts[my_scan_pt].sweep_pos % large_grid_sp);
      
      max_pt = search_maxpos(chr, start_pos, start_pos + large_grid_sp, p_snps,
			     scan_obj->chr_limits + chr, args->eval_range,
			     args->bp_resl, sm_p);

      if (max_pt.clr >= scan_obj->scan_pts[my_scan_pt].clr) {
	scan_obj->scan_pts[my_scan_pt].permute_p++;
	if (scan_obj->scan_pts[my_scan_pt].permute_p >= 20 &&
	    scan_obj->scan_pts[my_scan_pt].permute_p /(double)
	    scan_obj->scan_pts[my_scan_pt].permute_n >= rand()/(RAND_MAX + 1.0)) 
	  scan_obj->scan_pts[my_scan_pt].permute_finished = 1;
      }
      
      if (scan_obj->scan_pts[my_scan_pt].permute_n < CLR_NULL_DIST_SAVE)
	scan_obj->scan_pts[my_scan_pt].permute_clr[scan_obj->scan_pts[my_scan_pt].permute_n] = max_pt.clr;
      scan_obj->scan_pts[my_scan_pt].permute_n++;
      if (max_pt.clr < 0 || max_pt.clr > 1000000 || isnan(max_pt.clr)) {
	fprintf(stderr,"%d\t%d\t%g\t%1.3e\n",chr,start_pos,max_pt.clr, 
		exp(max_pt.lalpha));
      }

#if 0
#if 0
      chr = scan_obj->scan_pts[my_scan_pt].chr;
      start_pos = scan_obj->scan_pts[my_scan_pt].sweep_pos - 
	(scan_obj->scan_pts[my_scan_pt].sweep_pos % large_grid_sp);
      current_pos = start_pos;
      while(current_pos < start_pos + large_grid_sp) {
	init_scan_result(&pscan_pt, chr, p_snps, n_snps,
			 scan_obj->chr_limits[chr].start_index,
			 scan_obj->chr_limits[chr].n_snps, current_pos);
	search_maxalpha(&pscan_pt, p_snps, n_snps, sm_p);
	if (pscan_pt.clr >= scan_obj->scan_pts[my_scan_pt].clr) {
	  scan_obj->scan_pts[my_scan_pt].permute_p++;
	  if (scan_obj->scan_pts[my_scan_pt].permute_p >= 20)
	    scan_obj->scan_pts[my_scan_pt].permute_finished = 1;
	  break;
	}
	current_pos += small_grid_sp;
      }
      scan_obj->scan_pts[my_scan_pt].permute_n++;

#else
      pscan_pt = scan_obj->scan_pts[my_scan_pt];
      pscan_pt.null_logl = 0.;
      for(j=pscan_pt.window_start;j<=pscan_pt.window_end;j++)
	pscan_pt.null_logl += p_snps[j].null_logl;
      
      search_maxalpha(&pscan_pt, p_snps, n_snps, sm_p);

      scan_obj->scan_pts[my_scan_pt].permute_n++;
      if (scan_obj->scan_pts[my_scan_pt].clr <= pscan_pt.clr) 
	scan_obj->scan_pts[my_scan_pt].permute_p++;

      if (scan_obj->scan_pts[my_scan_pt].permute_p >= 20) 
	scan_obj->scan_pts[my_scan_pt].permute_finished = 1;
#endif
#endif
    }
  }

  pthread_mutex_unlock(&scan_pt_lock);
  return NULL;
}

static struct timeval stop_watch;
static scan_t *g_scan_obj;
extern char *output_fname;
extern int n_permute;
extern char *prepend_label;

static void output_clr_null_distribution(char *output_fname, 
					 scan_t *scan_obj);

static void interrupt_signal_handler(int signal) {

  if (elapsed_time(&stop_watch) < 10) {
    fprintf(stderr,"\nanother interrupt signal received, aborting "
	    "permutation\n");
    exit(-1);
  } else {
    scan_output(output_fname, g_scan_obj, 0, n_permute, prepend_label);
    output_clr_null_distribution(output_fname, g_scan_obj);
    gettimeofday(&stop_watch, NULL);
  }

}

static int scan_pt_clr_compare(scan_pt_t *a, scan_pt_t *b) __attribute__((unused));
static int scan_pt_clr_compare(scan_pt_t *a, scan_pt_t *b) {

  if (a->permute_finished == 0 && b->permute_finished == 1) return -1;
  if (a->permute_finished == 1 && b->permute_finished == 0) return 1;

  if (b->clr < a->clr) return -1;
  if (b->clr > a->clr) return  1;
  return 0;
}

void scan_permute(scan_t *scan_obj, sm_ptable_t *sm_p, 
		  int n_permute, double permute_nbp, 
		  double alpha_factor, 
		  int n_threads, int eval_range, int bp_resl, 
		  int large_grid_sp, double scan_width_mb) {
  int i;
  thread_args_t args;
  pthread_t *threads;
  struct sigaction int_sigaction;

#ifdef LN_DEBUG
  rfile = fopen("block-size.tsv", "w");
#endif

  /* Allow the user to get a current output at anytime during the permutation
     by pressing CTRL-C on the terminal. A second CTRL-C within 10 seconds will
     terminate the program */
  gettimeofday(&stop_watch, NULL);
  int_sigaction.sa_handler = interrupt_signal_handler;
  sigemptyset(&int_sigaction.sa_mask);
  int_sigaction.sa_flags = 0;
  g_scan_obj = scan_obj;

  sigaction(SIGINT, &int_sigaction, NULL);

  /* Since pthread_create() only allows a single (void *) to be passed to the
     thread starting function, we must pack the arguments into a struct and
     pass a pointer to the struct */
  args.scan_obj = scan_obj;
  args.sm_p = sm_p;

  args.permute_nbp = permute_nbp;
  args.n_permute = n_permute;

  args.large_grid_sp = large_grid_sp;
  args.eval_range = eval_range;
  args.bp_resl = bp_resl;
  args.scan_width_mb = scan_width_mb;
  
  MA(args.p_snps, sizeof(snp_t)*scan_obj->n_snps);

  MA(global_scan_pts, sizeof(int)*scan_obj->n_scan_pts);
  for(i=0;i<scan_obj->n_scan_pts;i++)
    global_scan_pts[i] = i;
  n_global_scan_pts = scan_obj->n_scan_pts;

  MA(threads, sizeof(pthread_t)*n_threads);
  pthread_mutex_init(&scan_pt_lock, NULL);
  pthread_barrier_init(&permute_barrier, NULL, n_threads);
  pthread_barrier_init(&run_barrier, NULL, n_threads);

  n_remaining = scan_obj->n_scan_pts;

  global_permute = -1;
  for(i=0;i<n_threads;i++) {
    pthread_create(threads + i, NULL, (void *) scan_permute_thread, 
		   &args);
  }

  for(i=0;i<n_threads;i++) pthread_join(threads[i], NULL);

  cr_logmsg(MSG_STATUS,"Scanning snp block permutations... finished.\n");
  pthread_mutex_destroy(&scan_pt_lock);
  pthread_barrier_destroy(&permute_barrier);
  pthread_barrier_destroy(&run_barrier);
  free(threads);

#ifdef LN_DEBUG
  fclose(rfile);
#endif
}

static int __attribute__((unused))float_compare(float *a, float *b) {
  if (*a < *b) return -1;
  if (*a > *b) return  1;
  return 0;
}

static int __attribute__((unused))double_compare(double *a, double *b) {
  if (*a < *b) return -1;
  if (*a > *b) return  1;
  return 0;
} 

void scan_output(char *output_fname, scan_t *scan_obj, int maximum_only, 
		 int n_permute, char *prepend_label) {
  int i;
  double max_clr;
  scan_pt_t *s_pt;
  char pos_str[40];
  FILE *out_f;

  if (output_fname) {
    out_f = fopen(output_fname, "w");
    if (out_f == NULL) {
      fprintf(stderr,"Can't open output file \"%s\" (%s)", output_fname,
	      strerror(errno));
      return;
    }
  } else {
    out_f = stdout;
  }

  s_pt = scan_obj->scan_pts;
  max_clr = scan_obj->scan_pts[0].clr;
  for(i=1;i<scan_obj->n_scan_pts;i++) {
    if (scan_obj->scan_pts[i].clr > max_clr) {
      max_clr = scan_obj->scan_pts[i].clr;
      s_pt = scan_obj->scan_pts + i;
    }
  }

  if (s_pt->sweep_pos > 1000000) {
    sprintf(pos_str, "chromosome %s %1.2f Mb", 
	    scan_obj->chr_limits[s_pt->chr].name, s_pt->sweep_pos/1e6);
  } else if (s_pt->sweep_pos > 2000) {
    sprintf(pos_str, "chromosome %s %1.2f Kb", 
	    scan_obj->chr_limits[s_pt->chr].name, s_pt->sweep_pos/1e3);
  } else {
    sprintf(pos_str, "chromosome %s %d bp", 
	    scan_obj->chr_limits[s_pt->chr].name, s_pt->sweep_pos);
  }
  logmsg(MSG_STATUS,"\rOutput complete -- maximum CLR of %g at %s "
	 "(alpha = %g)\n", max_clr, pos_str, exp(s_pt->lalpha));

  if (maximum_only) {
    if (prepend_label) fprintf(out_f,"%s\t",prepend_label);

    fprintf(out_f,"%s\t%d\t%1.2f\t%1.3e\t%d\t%d\t%d\n", 
	    scan_obj->chr_limits[s_pt->chr].name, s_pt->sweep_pos, max_clr,
	    exp(s_pt->lalpha), s_pt->n_snps, 
	    scan_obj->snps[s_pt->window_start].pos,
	    scan_obj->snps[s_pt->window_end].pos);
    return;
  }

  if (n_permute > 0) {
    double pvalue;
    
    for(i=0;i<scan_obj->n_scan_pts;i++) {
      s_pt = scan_obj->scan_pts + i;

      pvalue = (s_pt->permute_p + 0.5)/(double) (s_pt->permute_n + 0.5);

      if (prepend_label) fprintf(out_f,"%s\t",prepend_label);
      fprintf(out_f,"%s\t%d\t%1.2f\t%1.3e\t%d\t%d\t%1.3f\n", 
	      scan_obj->chr_limits[s_pt->chr].name, s_pt->sweep_pos, 
	      s_pt->clr, exp(s_pt->lalpha), s_pt->permute_p, s_pt->permute_n,
	      -log10(pvalue));
    }
  } else {
    for(i=0;i<scan_obj->n_scan_pts;i++) {
      s_pt = scan_obj->scan_pts + i;

      if (prepend_label) fprintf(out_f,"%s\t",prepend_label);
      fprintf(out_f,"%s\t%d\t%1.2f\t%1.3e\t%d\t%d\t%d\n", 
	      scan_obj->chr_limits[s_pt->chr].name, s_pt->sweep_pos, 
	      s_pt->clr, exp(s_pt->lalpha), s_pt->n_snps, 
	      scan_obj->snps[s_pt->window_start].pos, 
	      scan_obj->snps[s_pt->window_end].pos);
    }
  }
  
  if (output_fname) fclose(out_f);
}


static void output_clr_null_distribution(char *output_fname, 
					 scan_t *scan_obj) {
  int i, j;
  char *output_fname2;
  FILE *f;

  output_fname2 = alloca(strlen(output_fname)+20);
  sprintf(output_fname2, "%s-nulldist", output_fname);

  f = fopen(output_fname2, "w");
  if (f == NULL) {
    fprintf(stderr,"Can't open output file for CLR null distribution (%s)\n",
	    strerror(errno));
    return;
  }

  fprintf(f,"chr\tpos\tCLR\talpha\tp\tn");
  for(j=0;j<CLR_NULL_DIST_SAVE;j++)
    fprintf(f,"\t%1.4f",j/(double) CLR_NULL_DIST_SAVE);
  fprintf(f,"\n");

  for(i=0;i<scan_obj->n_scan_pts;i++) {
    scan_pt_t *s_pt;
    int n_pts;
    s_pt = scan_obj->scan_pts + i;

    if (CLR_NULL_DIST_SAVE < s_pt->permute_n)
      n_pts = CLR_NULL_DIST_SAVE;
    else
      n_pts = s_pt->permute_n;

    qsort(s_pt->permute_clr, n_pts, sizeof(float), (void *) float_compare);

    fprintf(f,"%s\t%d\t%1.3f\t%1.3e\t%d\t%d", 
	    scan_obj->chr_limits[s_pt->chr].name, s_pt->sweep_pos, 
	    s_pt->clr, exp(s_pt->lalpha), s_pt->permute_p, s_pt->permute_n);
    for(j=0;j<n_pts;j++)
      fprintf(f,"\t%1.2f",(double) s_pt->permute_clr[j]);
    fprintf(f,"\n");

  }

  fclose(f);
}
