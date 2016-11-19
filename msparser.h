#ifndef MSPARSER_H
#define MSPARSER_H
#include "kmacros.h"

typedef struct tree {
  double length; /* Length of branch from parent to this node */
  int sample_index;
  struct tree *left;
  struct tree *right;
} gtree_t;

typedef struct {
  gtree_t *gtree;
  gtree_t *gspace;
  int segment_size;
} segment_t;

typedef struct {
  int n_segments;
  segment_t *segments;
  gtree_t **gspace;
  double prob;
  int n_poly;
  double *positions;
  int n_haplotypes;
  char **haplotypes;
} msblock_t;

typedef struct {
  int n;
  double Tw;
  double Tpi;
  double Th;
  double Dt;
  double Dfl;
  double H;
} sfs_summary_t;


#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
typedef void* yyscan_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif
yyscan_t msparser_setfile(FILE *);
msblock_t *msparser_block();
void msparser_block_free(msblock_t *);
FILE *msparser_execute(char *ms_cmd);
int *msblock_sfs(msblock_t *msb, int s_index, int n);
sfs_summary_t *sfs_summaries(int *sfs, int n);
double *msblock_fsbranch_lengths(msblock_t *msb, int s_index, int n);
#ifdef __cplusplus
}
#endif
#endif
