%{

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <math.h>
#include <sys/types.h>
#include <unistd.h>


#include "msparser.h"

  int msget_lineno(yyscan_t);
  char *msget_text(yyscan_t);

#define mslineno msget_lineno(ms_scanner)
#define mstext msget_text(ms_scanner)
  //extern int mslineno;
  //extern char *mstext;
int mslex(yyscan_t yyscanner, void *);
static void mserror(msblock_t **msb_unused, yyscan_t ms_scanner, int stemp_ptr,
		    int htemp_ptr, int dtemp_ptr, int stemp_size, 
		    int htemp_size, int dtemp_size, double *dtemp, segment_t *stemp, char **htemp, char *error_msg);
void mslex_init(yyscan_t yyscanner);

#define STEMP_ALLOC_STEP (2)
#define HTEMP_ALLOC_STEP (2)
#define DTEMP_ALLOC_STEP (2)
%}

%define api.pure
%lex-param {yyscan_t ms_scanner}
%parse-param {msblock_t **return_msb}
%parse-param {yyscan_t ms_scanner}

%parse-param {int stemp_ptr}
%parse-param {int htemp_ptr}
%parse-param {int dtemp_ptr}
%parse-param {int stemp_size}
%parse-param {int htemp_size}
%parse-param {int dtemp_size}

%parse-param {double *dtemp}
%parse-param {segment_t *stemp}
%parse-param {char **htemp}


%union {
  char      *strtype;
  int        itype;
  double     dtype;
  double    *dptr;
  segment_t *sptr;
  gtree_t   *gptr;
  char      *hptr;
  msblock_t *rptr;
}

%token MS_START
%token OPEN_BRACKET
%token CLOSE_BRACKET
%token OPEN_PAREN
%token COLON
%token CLOSE_PAREN
%token COMMA
%token SEMICOLON
%token <strtype>    INTEGER
%token <dtype>  DECIMAL
%token <itype>  SEGSITES_TOKEN
%token PROB
%token POSITIONS_TOKEN
%token <strtype> HAPLOTYPE
%token BLANK_LINE
%token SCAN_ERROR
%token UNMATCHED_LINE;

%token <strtype> MS_CMDLINE

%type <rptr> msblock
%type <gptr> t
%type <itype> segsites
%type <itype> integer
%type <dtype> decimal
%type <strtype> haplotype
%type <gptr> tline
%type <dtype> tlength
%type <itype> norecomb_segsize
%type <dtype> maybeprob

%debug

%%

s: msblock
   ;

msblock:   { *return_msb = NULL; YYACCEPT; }
          | ms_start trees segsites maybeprob positions haplotypes maybeblankline {
	   msblock_t *msblock;
           int i;

	      MA(msblock, sizeof(msblock_t));
	      msblock->n_segments = stemp_ptr;
	      if (msblock->n_segments > 0) {
		MA(msblock->segments, sizeof(segment_t)*stemp_ptr);
		for(i=0;i<stemp_ptr;i++) {
		  msblock->segments[i].gtree = stemp[i].gtree;
		  msblock->segments[i].segment_size = stemp[i].segment_size;
		}
	      } else {
	         msblock->segments = NULL;
              }

	      msblock->prob = $4;

	      if (dtemp_ptr != $3) {
	         fprintf(stderr,"Different number of positions than "
		                "segregating sites.\n");
              }

	      msblock->n_poly = dtemp_ptr;
	      if (msblock->n_poly > 0) {
	        MA(msblock->positions, sizeof(double)*dtemp_ptr);
	        for(i=0;i<dtemp_ptr;i++) {
	          msblock->positions[i] = dtemp[i];
	        }
	      } else {
	        msblock->positions = NULL;
              }

	      msblock->n_haplotypes = htemp_ptr;
	      if (msblock->n_haplotypes > 0) {
	         MA(msblock->haplotypes, sizeof(char *)*htemp_ptr);
	         for(i=0;i<htemp_ptr;i++) {
	            msblock->haplotypes[i] = htemp[i];
	         }
	      } else {
	         msblock->haplotypes = NULL;
              }
	      *return_msb = msblock;
	      YYACCEPT;
	   }
	   ;

ms_start: MS_START { htemp_ptr = 0; dtemp_ptr = 0; stemp_ptr = 0; }
	  ;

maybeprob:                { $$ = 0.; }
           | PROB decimal { $$ = $2; }
           ;

maybeblankline: 
		| BLANK_LINE
		;

trees:   {}
         | trees norecomb_segsize tline {
	   if (stemp_ptr == stemp_size) {
              stemp_size += STEMP_ALLOC_STEP;
	      RA(stemp, sizeof(segment_t)*stemp_size);
           }
	   stemp[stemp_ptr].gtree = $3;
	   stemp[stemp_ptr].segment_size = $2;
	   stemp_ptr++;
         }
	 ;

norecomb_segsize:   { $$ = 1; /* if run without recombination, no segment
                                 size token will appear -- the whole locus
                                 is thus one unit length */ }
                  | OPEN_BRACKET integer CLOSE_BRACKET { $$ = $2; } 
                  ;

tline:   t SEMICOLON { $$ = $1; }

t:   integer COLON tlength {
        gtree_t *g;

	MA(g, sizeof(gtree_t));
        g->sample_index = $1;
        g->length = $3;
        g->left = NULL;
        g->right = NULL;
	$$ = g;
     } 
   | OPEN_PAREN t COMMA t CLOSE_PAREN COLON tlength {
        gtree_t *g;

	MA(g, sizeof(gtree_t));
        g->sample_index = 0;
        g->length = $7;
        g->left = $2;
        g->right = $4;
        $$ = g;
     }
   | OPEN_PAREN t COMMA t CLOSE_PAREN {
        gtree_t *g;

	MA(g, sizeof(gtree_t));
        g->sample_index = 0;
        g->length = 0.;
        g->left = $2;
        g->right = $4;
        $$ = g;
   }
   ;

tlength:   integer { $$ = (double) $1; }
	 | decimal { $$ = $1; }
	 ;

segsites:  { $$ = 0; }
	   | SEGSITES_TOKEN integer { $$ = $2; }
	   ;

positions: 
	   | POSITIONS_TOKEN decimal_list;
	   ;

decimal_list: {  }
	      | decimal_list decimal { 
		  if (dtemp_ptr == dtemp_size) {
		    dtemp_size += DTEMP_ALLOC_STEP;
		    RA(dtemp, sizeof(double)*dtemp_size);
		  }
		  dtemp[dtemp_ptr++] = $2;
	      }
	      ;

decimal: DECIMAL { $$ = strtod(mstext, NULL); }

haplotypes:   {  }
	      |	haplotypes haplotype {
	          if (htemp_ptr == htemp_size) {
		    htemp_size += HTEMP_ALLOC_STEP;
		    RA(htemp, sizeof(char *)*htemp_size);
		  }
		  htemp[htemp_ptr++] = $2;
	      }
	      ;

haplotype: HAPLOTYPE { $$ = strdup(mstext); };

integer: INTEGER { $$ = atoi(mstext); }

%%

static void mserror(msblock_t **msb_unused, yyscan_t ms_scanner, int stemp_ptr,
		    int htemp_ptr, int dtemp_ptr, int stemp_size, 
		    int htemp_size, int dtemp_size, double *dtemp, segment_t *stemp, char **htemp, char *error_msg) {
       fprintf(stderr,"Parsing error (%s) at line %d of ms output\n",
       error_msg, mslineno);
       fprintf(stderr,"Current mstext: \"%s\"\n",mstext);
}

msblock_t *msparser_block(yyscan_t ms_scanner) {
  msblock_t *msb;
  int stemp_ptr, htemp_ptr, dtemp_ptr;
  int stemp_size = 0;
  int htemp_size = 0;
  int dtemp_size = 0;

  double *dtemp = NULL;
  segment_t *stemp = NULL;
  char **htemp = NULL;

//  msdebug = 1;
  msb = NULL;
  stemp_ptr = htemp_ptr = dtemp_ptr = 0;
  msparse(&msb, ms_scanner, stemp_ptr, htemp_ptr, dtemp_ptr, stemp_size, htemp_size, dtemp_size, dtemp, stemp, htemp);
  return msb;
}

static void gtree_free(gtree_t *g) {
  if (g == NULL) return;
  gtree_free(g->left);
  gtree_free(g->right);
  free(g);
}

FILE *msparser_execute(char *ms_cmd) {
  FILE *f;
  int pfd[2], n_args, k;
  char *p, **cmdargs;
  
  n_args = 0;
  p = ms_cmd;
  while(*p!=0) {
    n_args++;
    while(*p!=0 && *p!=' ') p++;
    while(*p!=0 && *p==' ') p++;
  }
  
  ms_cmd = strdup(ms_cmd);
  MA(cmdargs, sizeof(char *)*(n_args+1));
  k = 0;
  p = ms_cmd;
  while(*p!=0) {
    cmdargs[k++] = p;
    while(*p!=0 && *p!=' ') p++;
    if (*p != 0) {
      while(*p!=0 && *p==' ') { 
	*p = 0; 
	p++; 
      }
    }
  }
  cmdargs[k] = NULL;

  if (pipe(pfd)!=0) {
    fprintf(stderr,"Can't open pipe to ms (%s)\n",strerror(errno));
    exit(1);
  }

  if (fork()==0) {
    close(0);
    close(1);
    close(2);
    close(pfd[0]);
    dup2(pfd[1], 1);
    dup2(pfd[1], 2);
    execvp("ms", cmdargs);
    fprintf(stderr,"Can't execute ms (%s)\n",strerror(errno));
    exit(1);
  }

  close(pfd[1]);
  f = fdopen(pfd[0], "r");
  if (f == NULL) {
    fprintf(stderr,"Can't fdopen() pipe from ms (%s)\n", strerror(errno));
    exit(1);
  }

  free(cmdargs);
  free(ms_cmd);

  return f;
}

void msparser_block_free(msblock_t *msb) {
  int i;

  for(i=0;i<msb->n_haplotypes;i++) {
    free(msb->haplotypes[i]);
  }
  free(msb->haplotypes);
  free(msb->positions);
  for(i=0;i<msb->n_segments;i++)
    gtree_free(msb->segments[i].gtree);
  if (msb->n_segments) free(msb->segments);
  free(msb);
}

static int bfs_descend(gtree_t *node, int segment_size, int s_index, int n, 
		       double *bfs) {
  int s;

  if (node->left == NULL) {
    if (node->sample_index >= s_index && node->sample_index < s_index + n) {
      bfs[1] += node->length*segment_size;
    } else {
      bfs[0] += node->length*segment_size;
    }
    return 1;
  }

  s = bfs_descend(node->left, segment_size, s_index, n, bfs) +
    bfs_descend(node->right, segment_size, s_index, n, bfs);

  bfs[s] += node->length*segment_size;
  return s;
}

double *msblock_fsbranch_lengths(msblock_t *msb, int s_index, int n) {
  int i;
  double *bfs;

  MA(bfs, sizeof(double)*(n+1));
  for(i=0;i<=n;i++) bfs[i] = 0.;

  for(i=0;i<msb->n_segments;i++)
    bfs_descend(msb->segments[i].gtree, msb->segments[i].segment_size,
		s_index, n, bfs);
  
  return bfs;
}

int *msblock_sfs(msblock_t *msb, int s_index, int n) {
  int i, j, q, *sfs;

  MA(sfs, sizeof(int)*(n+1));
  for(i=0;i<=n;i++) 
    sfs[i] = 0;
  
  for(j=0;j<msb->n_poly;j++) {
    q = 0;
    for(i=s_index;i<=s_index+n && i<msb->n_haplotypes;i++)
      if (msb->haplotypes[i][j] == '1') q++;
    sfs[q]++;
  }

  return sfs;
}

sfs_summary_t *sfs_summaries(int *sfs, int n) {

  int i, s;
  double a1, a2, b1, b2,c1, c2, e1, e2, vd, ud, v;
  sfs_summary_t *ss;

  MA(ss, sizeof(sfs_summary_t));

  s = 0;
  a1 = a2 = 0.;
  for(i=1;i<n;i++) {
    s += sfs[i];

    /* Constants for Tajima's D and Fu and Li's D */
    a1 += 1.0/i;
    a2 += 1.0/(i*i);
  }

  ss->Tw = ss->Tpi = ss->Th = 0.;
  ss->Dt = ss->Dfl = ss->H = 0.;
  ss->n = s;

  if (s == 0) return ss;


  /* For variance of Tajima's D */
  b1 = (n+1)/((double) 3*(n-1));
  b2 = (2*(n*n + n + 3))/(double) (9*n*(n-1));
  c1 = b1 - 1.0/a1;
  c2 = b2 - (n+2.)/(a1*n) + a2/(a1*a1);
  e1 = c1/a1;
  e2 = c2/(a1*a1 + a2);

  /* For variance of Fu and Li's D */
  vd = 1. + (a1/(a2 + a1*a1))*((2*n*a1-4.*(n-1)-(n+1)*(n-2))/(double) 
			       ((n-1)*(n-2)));
  ud = a1 - 1. - vd;
  v =  ud*s + vd*s*s;

  ss->Tw = s/a1;
  for(i=1;i<n;i++) {
    ss->Tpi += i*(n-i)*sfs[i];
    ss->Th += i*i*sfs[i];
  }
  ss->Tpi /= (n*(n-1)/2.0);
  ss->Th = ss->Th*(2.0/(double) (n*(n-1)));

  ss->Dt = (ss->Tpi - ss->Tw)/sqrt(e1*s + e2*s*(s-1));
  ss->H = (ss->Tpi - ss->Th);
  ss->Dfl = (s - sfs[1]*a1)/sqrt(v);
  
  return ss;

}
 
 
