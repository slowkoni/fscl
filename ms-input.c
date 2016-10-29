#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <kmacros.h>
#include <msparser.h>

#include "fscl.h"

static FILE *msf = NULL;

void ms_openfile(char *ms_fname) {

  if (msf != NULL) fclose(msf);

  msf = fopen(ms_fname, "r");
  if (msf == NULL)
    logmsg(MSG_FATAL, "Can't open ms input file \"%s\" (%s)", ms_fname,
	   strerror(errno));

  msparser_setfile(msf);
}

scan_t *ms_background(char *ms_fname, int ms_segment_length, int ms_folded,
		      int ms_sample_first, int ms_sample_size) {
  int i, j, m, d, n_blocks, n_snps, block_sample_size;
  msblock_t *msb;
  scan_t *scan_obj;

  n_snps = 0;
  MA(scan_obj, sizeof(scan_t));
  scan_obj->snps = NULL;
  scan_obj->n_snps = 0;
  scan_obj->n_depths = 0;
  scan_obj->sample_depths = NULL;
  scan_obj->scan_pts = NULL;
  scan_obj->n_scan_pts = 0;

  ms_openfile(ms_fname);
  n_blocks = 0;
  while((msb = msparser_block())) {
    m = 0;
    if (ms_sample_size == 0) 
      block_sample_size = msb->n_haplotypes - ms_sample_first;
    else
      block_sample_size = ms_sample_size;

    while(m<scan_obj->n_depths && 
	  scan_obj->sample_depths[m] != block_sample_size) 
      m++;
    if (m == scan_obj->n_depths) {
      if (scan_obj->n_depths % 8 == 0)
	RA(scan_obj->sample_depths, sizeof(int)*(scan_obj->n_depths + 8));
      scan_obj->sample_depths[scan_obj->n_depths++] = block_sample_size;
    }

    

    for(i=0;i<msb->n_poly;i++) {
      if (n_snps % SNP_ALLOC_STEP == 0)
	RA(scan_obj->snps, sizeof(snp_t)*(n_snps + SNP_ALLOC_STEP));
      
      scan_obj->snps[n_snps].pos = msb->positions[i]*ms_segment_length +
	(n_blocks * ms_segment_length);

      scan_obj->snps[n_snps].obs_freq = 0;
      scan_obj->snps[n_snps].depth_p = m;
      d = 0;
      for(j=ms_sample_first;j<ms_sample_first + block_sample_size;j++)
	if (msb->haplotypes[j][i] == '1') d++;
      if (d == 0 || d == block_sample_size) continue;
	
      if (ms_folded) {
	if (d > scan_obj->sample_depths[m] - d)
	  scan_obj->snps[n_snps].obs_freq = scan_obj->sample_depths[m] - d;
	scan_obj->snps[n_snps].folded = 1;
      } else {
	scan_obj->snps[n_snps].obs_freq = d;
	scan_obj->snps[n_snps].folded = 0;
      }
      n_snps++;
    }
    
    msparser_block_free(msb);
    n_blocks++;
  }

  scan_obj->n_snps = n_snps;
  return scan_obj;
}

scan_t *ms_next_block(int ms_segment_length, int ms_folded, 
		      int ms_sample_first, int ms_sample_size) {
  int i, j, m, d, n_snps;
  msblock_t *msb;
  scan_t *scan_obj;

  n_snps = 0;
  MA(scan_obj, sizeof(scan_t));
  scan_obj->snps = NULL;
  scan_obj->n_snps = 0;
  scan_obj->n_depths = 0;
  scan_obj->sample_depths = NULL;
  scan_obj->scan_pts = NULL;
  scan_obj->n_scan_pts = 0;

  msb = msparser_block();
  if (msb == NULL) return NULL;

  if (ms_sample_size == 0) ms_sample_size = msb->n_haplotypes-ms_sample_first;

  m = 0;
  while(m<scan_obj->n_depths && 
	scan_obj->sample_depths[m] != msb->n_haplotypes) m++;
  if (m == scan_obj->n_depths) {
    if (scan_obj->n_depths % 8 == 0)
      RA(scan_obj->sample_depths, sizeof(int)*(scan_obj->n_depths + 8));
    scan_obj->sample_depths[m] = ms_sample_size;
  }


  for(i=0;i<msb->n_poly;i++) {
    if (n_snps % SNP_ALLOC_STEP == 0)
      RA(scan_obj->snps, sizeof(snp_t)*(n_snps + SNP_ALLOC_STEP));
    
    scan_obj->snps[n_snps].pos = msb->positions[i]*ms_segment_length;
    
    scan_obj->snps[n_snps].obs_freq = 0;
    scan_obj->snps[n_snps].depth_p = m;
    d = 0;
    for(j=ms_sample_first;j<ms_sample_first + ms_sample_size;j++)
      if (msb->haplotypes[j][i] == '1') d++;
    if (d == 0 || d == ms_sample_size) continue;

    
    if (ms_folded) {
      if (d > scan_obj->sample_depths[m] - d)
	scan_obj->snps[n_snps].obs_freq = scan_obj->sample_depths[m] - d;
      scan_obj->snps[n_snps].folded = 1;
    } else {
      scan_obj->snps[n_snps].obs_freq = d;
      scan_obj->snps[n_snps].folded = 0;
    }
    n_snps++;
  }    
  msparser_block_free(msb);

  scan_obj->n_snps = n_snps;
  return scan_obj;
}
