#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <float.h>

#include "kmacros.h"
#include "cuda-macros.h"

__device__ double cuda_lchoose(double *cu_lft, int n, int k) {

  if (k < 0) return -DBL_MAX;
  if (n == 0 && k == 0) return 0;
  if (k > n || n == 0) return -DBL_MAX;

  return cu_lft[n] - (cu_lft[k] + cu_lft[n-k]);
  //return lgamma((double) n+1) - (lgamma((double) k+1) + lgamma((double) (n-k+1)));
}

double lchoose(double *lft, int n, int k) {
  if (k < 0) return -DBL_MAX;
  if (n == 0 && k == 0) return 0;
  if (k > n || n == 0) return -DBL_MAX;

  return lft[n] - (lft[k] + lft[n-k]);
}

__global__ void cuda_pjh(double *pjh, double *fsp, double *cu_lft, int sample_size) {
  int i, j, h;
  j = blockDim.x * blockIdx.x + threadIdx.x;
  h = blockDim.y * blockIdx.y + threadIdx.y;

  if (j > sample_size) return;
  if (h > sample_size) return;

  for(i=j; i <= sample_size; i++) {
    pjh[j*(sample_size + 1) + h] += fsp[i]*exp(cuda_lchoose(cu_lft, i,j) + 
					       cuda_lchoose(cu_lft, sample_size - i, h - j) -
					       cuda_lchoose(cu_lft, sample_size, h));
  }
}

__global__ void cuda_pbk(double *pbk, double *pjh, int sample_size) {
  int b, k, q;

  b = blockDim.x * blockIdx.x + threadIdx.x;
  k = blockDim.y * blockIdx.y + threadIdx.y;
  
  if (b > sample_size) return;
  if (k >= sample_size) return;

  q = b - (sample_size - k) + 1;
  if (q > 0) {
    pbk[b*(sample_size+1) + k] += pjh[q*(sample_size+1) + k+1]*(q/(k+1.));
  }
  if (b < k+1) {
    pbk[b*(sample_size+1) + k] += pjh[b*(sample_size+1) + k+1]*((k+1.-b)/(k+1.));
  }
}
   
void cuda_calculate_pjh(double *pjh, double *fsp, double *cu_lft, int sample_size) {
  int n = sample_size;
  double *cu_pjh, *cu_pbk, *pbk;
  double *cu_fsp;
  

  CUDA_MA(cu_fsp, sizeof(double)*(n+1));
  CUDA_MEMCPY_TO(cu_fsp, fsp, sizeof(double)*(n+1));

  memset(pjh, 0x0, sizeof(double)*(n+1)*(n+1));
  CUDA_MA(cu_pjh, sizeof(double)*(n+1)*(n+1));
  CUDA_MEMCPY_TO(cu_pjh, pjh, sizeof(double)*(n+1)*(n+1));

  MA(pbk, (n+1)*(n+1), double);
  memset(pbk, 0x0, sizeof(double)*(n+1)*(n+1));
  CUDA_MA(cu_pbk, sizeof(double)*(n+1)*(n+1));
  CUDA_MEMCPY_TO(cu_pbk, pjh, sizeof(double)*(n+1)*(n+1));

  dim3 dimBlock(8,8);
  dim3 dimGrid((n+1)/dimBlock.x + 1, (n+1)/dimBlock.y + 1);
  cuda_pjh<<<dimGrid, dimBlock>>>(cu_pjh, cu_fsp, cu_lft, sample_size);

#if 1
  cuda_pbk<<<dimGrid, dimBlock>>>(cu_pbk, cu_pjh, sample_size);
#endif

  CUDA_MEMCPY_FROM(pjh, cu_pjh, sizeof(double)*(n+1)*(n+1));
  CUDA_MEMCPY_FROM(pbk, cu_pbk, sizeof(double)*(n+1)*(n+1));
  
  CUDA_FREE(cu_pbk);
  CUDA_FREE(cu_pjh);
  CUDA_FREE(cu_fsp);
}

static double *mk_log_factorial_table(double **r_lft, int n) {
  double *t, *cu_lft;
  int i;
  
  MA(t, n+1, double);

  t[0] = 0;
  t[1] = 0;
  for(i=2;i <= n; i++)
    t[i] = lgamma(i+1.);
  //    t[i] = t[i-1] + log(i);

  CUDA_MA(cu_lft, sizeof(double)*(n+1));
  CUDA_MEMCPY_TO(cu_lft, t, sizeof(double)*(n+1));
  if (r_lft != NULL) *r_lft = t;
  else free(t);
  return cu_lft;
}

double *cuda_pjh_init_tables(int n) {
  return mk_log_factorial_table(NULL, n);
}

#ifdef UNIT_TEST
static double *neutral_fsp(int n) {
  int i;
  double *fsp, s;
  
  MA(fsp, (n+1), double);
  
  fsp[0] = 0.;
  fsp[n] = 0.;
  s = 0.;
  for(i=1;i<n;i++) {
    fsp[i] = 1./i;
    s += fsp[i];
  }
  for(i=1;i<n;i++) {
    fsp[i] /= s;
  }

  return fsp;
}

void cpu_calculate_pjh(double **pjh, double *fsp, double *lft, int sample_size) {
  int i, j, h;
  
  for(j=0;j<=sample_size;j++) {      
    for(h=0;h<=sample_size;h++) {
      pjh[j][h] = 0.;
      
      for(i=j;i<=sample_size;i++) {
	pjh[j][h] += fsp[i]*exp(lchoose(lft, i,j) + 
				lchoose(lft, sample_size - i, h - j) -
				lchoose(lft, sample_size, h));
      }
    }
  }

}

int main(int argc, char *argv[]) {
  double **pjh_gpu, **pjh_cpu;
  double *cu_lft, *lft, *fsp;
  double cpu_ms, gpu_ms;
  int sample_size, i, j, h;
  struct timeval stopwatch_gpu, stopwatch_cpu;
  
  if (argc < 2) {
    fprintf(stderr,"\nSpecify sample size as first command line argument.\n");
    exit(-1);
  }
  sample_size = atoi(argv[1]);

  fsp = neutral_fsp(sample_size+1);
  cu_lft = mk_log_factorial_table(&lft, sample_size+1);

  MA(pjh_gpu, sample_size+1, double *);
  CA(pjh_gpu[0], (sample_size+1)*(sample_size+1), double);
  MA(pjh_cpu, sample_size+1, double *);
  CA(pjh_cpu[0], (sample_size+1)*(sample_size+1), double);
  for(i=1; i<=sample_size; i++) {
    pjh_gpu[i] = pjh_gpu[i-1] + (sample_size+1);
    pjh_cpu[i] = pjh_cpu[i-1] + (sample_size+1);
  }
  
  fprintf(stderr,"Computing pjh[][] for sample size %d -   ", sample_size);

  fprintf(stderr,"GPU: ");
  gettimeofday(&stopwatch_gpu, NULL);
  cuda_calculate_pjh(pjh_gpu[0], fsp, cu_lft, sample_size);
  gpu_ms = elapsed_time_ms(&stopwatch_gpu);
  fprintf(stderr," %1.1f ms", gpu_ms);

  fprintf(stderr,"\tCPU: ");
  gettimeofday(&stopwatch_cpu, NULL);
  cpu_calculate_pjh(pjh_cpu, fsp, lft, sample_size);
  cpu_ms = elapsed_time_ms(&stopwatch_cpu);
  fprintf(stderr," %1.1f ms\n", cpu_ms);

  double pjh_error = 0.;
  for(j=0; j<=sample_size; j++) {
    for(h=0; h<=sample_size; h++) {
      pjh_error += fabs(pjh_cpu[j][h] - pjh_gpu[j][h])/(pjh_cpu[j][h] + DBL_EPSILON);
      //      fprintf(stdout,"%d\t%d\t%1.5f\t%1.5f\t%1.5f\t%1.5f\n", j, h, fsp[h], pjh_cpu[j][h], pjh_gpu[j][h],
      //	      (pjh_cpu[j][h] - pjh_gpu[j][h])/(pjh_cpu[j][h]+DBL_EPSILON));
    }
  }
  pjh_error /= (sample_size+1)*(sample_size+1);
  fprintf(stderr,"Average %% deviance of GPU calculation from CPU calculation: %g\n",
	  pjh_error*100.);
  fprintf(stderr,"GPU is %1.1fX faster than CPU\n", cpu_ms/gpu_ms);
  
  free(pjh_gpu[0]);
  free(pjh_gpu);

  free(pjh_cpu[0]);
  free(pjh_cpu);

  CUDA_FREE(cu_lft);
  free(lft);
  free(fsp);
  return 0;
}
#endif
