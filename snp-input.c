#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <kmacros.h>
#include "fscl.h"

static int pos_compare(snp_t *a, snp_t *b) {
  if (a->chr < b->chr) return -1;
  if (a->chr > b->chr) return  1;
  if (a->pos < b->pos) return -1;
  if (a->pos > b->pos) return  1;
  return 0;
}

scan_t *load_snp_input(char *snp_fname, int include_invariant) {
  FILE *f;
  int i, j, obs_freq, sample_size, folded, line_no, pos, max_pos;
  int n_chromosomes, chr_index;
  char inputline[8192], chr_name[8192], **chr_names;
  scan_t *scan_obj;
  chr_limits_t *chr_limits;

  f = fopen(snp_fname, "r");
  if (f == NULL) {
    fprintf(stderr,"Can't open snp file \"%s\" (%s)\n",snp_fname,
	    strerror(errno));
    exit(-1);
  }

  MA(scan_obj, sizeof(scan_t));

  scan_obj->n_snps = 0;
  scan_obj->snps = NULL;
  scan_obj->n_depths = 0;
  scan_obj->sample_depths = NULL;
  scan_obj->n_scan_pts = 0;
  scan_obj->scan_pts = NULL;

  line_no = 0;
  max_pos = 0;
  n_chromosomes = 0;
  chr_index = -1;
  chr_names = NULL;
  while(fgets(inputline, 8192, f)) {
    line_no++;
    CHOMP(inputline);
    if (inputline[0] == 0) continue;
    if (inputline[0] == '#') continue;

    if (sscanf(inputline,"%s %d %d %d %d", chr_name, &pos, &obs_freq, 
	       &sample_size, &folded) != 5) {
      if (strcmp(inputline,"chromosome")!=0) {
	fprintf(stderr,"Can't parse SNP input at line %d: \"%s\"\n", line_no,
		inputline);
      }
      continue;
    }

    if (sample_size < 5) continue;
    if (!include_invariant && (obs_freq < 1 || obs_freq > sample_size-1))
      continue;

    if (chr_index == -1 || strcmp(chr_name, chr_names[chr_index])!=0) {
      chr_index=0;
      while(chr_index<n_chromosomes && 
	    strcmp(chr_names[chr_index], chr_name)!=0) chr_index++;
      if (chr_index == n_chromosomes) {
	if (n_chromosomes % 32 == 0)
	  RA(chr_names, sizeof(char *)*(n_chromosomes + 32));
	chr_names[n_chromosomes] = strdup(chr_name);
	chr_index = n_chromosomes;
	n_chromosomes++;
      }
    }

    if (pos > max_pos) max_pos = pos;


    if (scan_obj->n_snps % SNP_ALLOC_STEP == 0) {
      RA(scan_obj->snps, sizeof(snp_t)*(scan_obj->n_snps + SNP_ALLOC_STEP));
    }

    if (folded && obs_freq > (sample_size - obs_freq))
      obs_freq = sample_size - obs_freq;

    scan_obj->snps[scan_obj->n_snps].chr = chr_index;
    scan_obj->snps[scan_obj->n_snps].pos = pos;
    scan_obj->snps[scan_obj->n_snps].obs_freq = obs_freq;
    scan_obj->snps[scan_obj->n_snps].folded = folded;

    j = 0;
    while(j<scan_obj->n_depths && scan_obj->sample_depths[j]!=sample_size) j++;
    if (j == scan_obj->n_depths) {
      if (scan_obj->n_depths % 32 == 0) 
	RA(scan_obj->sample_depths, sizeof(int)*(scan_obj->n_depths + 32));
      scan_obj->sample_depths[scan_obj->n_depths++] = sample_size;
    }
    scan_obj->snps[scan_obj->n_snps].depth_p = j;
    scan_obj->n_snps++;
  }
  fclose(f);

  if (scan_obj->n_snps == 0) {
    fprintf(stderr,"No usable snps found in file \"%s\"\n",snp_fname);
    exit(1);
  }
  
  qsort(scan_obj->snps, scan_obj->n_snps, sizeof(snp_t), (void *) pos_compare);

  MA(chr_limits, sizeof(chr_limits_t)*n_chromosomes);
  i = 0;
  while(i<scan_obj->n_snps) {
    j = i;
    while(j<scan_obj->n_snps && scan_obj->snps[j].chr == scan_obj->snps[i].chr)
      j++;
    chr_limits[scan_obj->snps[i].chr].chr = scan_obj->snps[i].chr;
    chr_limits[scan_obj->snps[i].chr].start_index = i;
    chr_limits[scan_obj->snps[i].chr].n_snps = j - i;
    chr_limits[scan_obj->snps[i].chr].bp_length = scan_obj->snps[j-1].pos;
    chr_limits[scan_obj->snps[i].chr].name = chr_names[scan_obj->snps[i].chr];
    i = j;
  }    
  
  scan_obj->n_chromosomes = n_chromosomes;
  scan_obj->chr_limits = chr_limits;

  free(chr_names);
  return scan_obj;
}
