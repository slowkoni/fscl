#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <float.h>

#include "kmacros.h"
#include "cuda-macros.h"

extern "C" {

__device__ double cuda_lchoose(double *cu_lft, int n, int k) {

  if (k < 0) return -DBL_MAX;
  if (n == 0 && k == 0) return 0;
  if (k > n || n == 0) return -DBL_MAX;

  return cu_lft[n] - (cu_lft[k] + cu_lft[n-k]);
  //return lgamma((double) n+1) - (lgamma((double) k+1) + lgamma((double) (n-k+1)));
}

  static double __attribute__((unused))lchoose(double *lft, int n, int k) {
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

void cuda_calculate_pbk(double *pjh, double *pbk, double **r_cu_pbk, double **r_cu_fsp,
			double *fsp, double *cu_lft, int sample_size) {
  int n = sample_size;
  double *cu_pjh, *cu_pbk;
  double *cu_fsp;
  

  CUDA_MA(cu_fsp, sizeof(double)*(n+1));
  CUDA_MEMCPY_TO(cu_fsp, fsp, sizeof(double)*(n+1));

  memset(pjh, 0x0, sizeof(double)*(n+1)*(n+1));
  CUDA_MA(cu_pjh, sizeof(double)*(n+1)*(n+1));
  CUDA_MEMCPY_TO(cu_pjh, pjh, sizeof(double)*(n+1)*(n+1));

  memset(pbk, 0x0, sizeof(double)*(n+1)*(n+1));
  CUDA_MA(cu_pbk, sizeof(double)*(n+1)*(n+1));
  CUDA_MEMCPY_TO(cu_pbk, pbk, sizeof(double)*(n+1)*(n+1));

  dim3 dimBlock(8,8);
  dim3 dimGrid((n+1)/dimBlock.x + 1, (n+1)/dimBlock.y + 1);
  cuda_pjh<<<dimGrid, dimBlock>>>(cu_pjh, cu_fsp, cu_lft, sample_size);
  CUDA_MEMCPY_FROM(pjh, cu_pjh, sizeof(double)*(n+1)*(n+1));

  {
    dim3 dimBlock(8,8);
    dim3 dimGrid((n+1)/dimBlock.x + 1, (n+1)/dimBlock.y + 1);
    cuda_pbk<<<dimGrid, dimBlock>>>(cu_pbk, cu_pjh, sample_size);
  }
  CUDA_MEMCPY_FROM(pbk, cu_pbk, sizeof(double)*(n+1)*(n+1));

  CUDA_FREE(cu_pjh);

  if (r_cu_pbk != NULL)
    *r_cu_pbk = cu_pbk;
  else
    CUDA_FREE(cu_pbk);
  if (r_cu_fsp != NULL)
    *r_cu_fsp = cu_fsp;
  else
    CUDA_FREE(cu_fsp);
  
}

__device__ double p_kescape_gpu(double *cu_lft, int k, int n, double ad) {
  if (k == 0) return exp(-n*ad);
  return exp(cuda_lchoose(cu_lft,n,k) + k*log(1.0 - exp(-ad)) - (n-k)*ad);
}

static double __attribute__((unused))p_kescape_cpu(double *lft, int k, int n, double ad) {
  if (k == 0) return exp(-n*ad);
  return exp(lchoose(lft,n,k) + k*log(1.0 - exp(-ad)) - (n-k)*ad);
}

__global__ void cuda_spf(double *cu_spf, double *cu_pbk, double *cu_fsp, double *cu_x,
			 double *cu_lft, int spline_pts, int sample_size) {
  int s, f, k;
  double ad;
  
  s = blockDim.x * blockIdx.x + threadIdx.x;
  f = blockDim.y * blockIdx.y + threadIdx.y;

  if (s > spline_pts || f > sample_size) return;
  
  ad = exp(cu_x[s]);
  
  cu_spf[s*(sample_size+1) + f] = p_kescape_gpu(cu_lft, sample_size, sample_size, ad)*cu_fsp[f];
  for(k=0; k < sample_size; k++)
    cu_spf[s*(sample_size+1) + f] +=
      p_kescape_gpu(cu_lft, k, sample_size, ad)*cu_pbk[f*(sample_size+1) + k];
  
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

void cpu_calculate_pbk(double **pjh, double **pbk, double *fsp, double *lft, int sample_size) {
  int i, j, h, b, k, q;
  
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

  for(b=0;b<=sample_size;b++) {
    for(k=0;k<sample_size;k++) {
      q = b - (sample_size - k) + 1;
      if (q > 0) {
	pbk[b][k] += pjh[q][k+1]*(q/(double) (k+1));
      }
      if (b < k+1) {
	pbk[b][k] += pjh[b][k+1]*((k+1-b)/(double) (k+1));
      }
    }
  }

}

void cuda_calculate_spline_pts(double *x, double *spf, double *cu_pbk,
			       double *cu_fsp, double *cu_lft,
			       int spline_pts, int sample_size) {
  double *cu_x, *cu_spf;
  
  memset(spf, 0x0, sizeof(double)*(spline_pts+1)*(sample_size+1));

  CUDA_MA(cu_x,   sizeof(double)*(spline_pts+1));
  CUDA_MA(cu_spf, sizeof(double)*(spline_pts+1)*(sample_size+1));

  CUDA_MEMCPY_TO(cu_x, x, sizeof(double)*(spline_pts+1));
  CUDA_MEMCPY_TO(cu_spf, spf, sizeof(double)*(spline_pts+1)*(sample_size+1));

  dim3 dimBlock(8,8);
  dim3 dimGrid((spline_pts + 1)/dimBlock.x + 1, (sample_size+1)/dimBlock.y + 1);
  cuda_spf<<<dimGrid, dimBlock>>>(cu_spf, cu_pbk, cu_fsp, cu_x, cu_lft, spline_pts,
				  sample_size);

  CUDA_MEMCPY_FROM(spf, cu_spf,  sizeof(double)*(spline_pts+1)*(sample_size+1));

  CUDA_FREE(cu_spf);
  CUDA_FREE(cu_x);
  CUDA_FREE(cu_pbk);
  CUDA_FREE(cu_fsp);

}
			       
void cpu_calculate_spline_pts(double *x, double **spf, double **pbk, double *fsp, double *lft,
			      int spline_pts, int sample_size) {
  int i, f, k;
  
  for(i=0; i<=spline_pts; i++) {
    double ad;

    ad = exp(x[i]);
    for(f=0;f<=sample_size;f++) {
      spf[i][f] = p_kescape_cpu(lft, sample_size, sample_size, ad)*fsp[f];
      //      fprintf(stderr,"%d\t%d\t%g\t%g\n",i,f,x[i],spf[i][f]);
      for(k=0;k<sample_size;k++) {
	double tmp;
	tmp = p_kescape_cpu(lft, k, sample_size, ad)*pbk[f][k];
	//	fprintf(stderr,"IFK: %d\t%d\t%g\t%g\t%g\n",i,f,x[i],pbk[f][k],tmp);
	spf[i][f] += tmp;
	//fprintf(stderr,"IFK: %d\t%d\t%g\t%g\t%g\n",i,f,x[i],pbk[f][k],spf[i][f]);
      }
      //fprintf(stderr,"%d\t%d\t%g\t%g\n",i,f,x[i],spf[i][f]);
    }
  }
}

#ifdef __cplusplus
}
#endif
#ifdef UNIT_TEST
#define LOG_AD_MIN (-20.0)
#define LOG_AD_MAX (4.0)
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

int main(int argc, char *argv[]) {
  double **pjh_gpu, **pjh_cpu, **pbk_gpu, **pbk_cpu, **spf_gpu, **spf_cpu;
  double *cu_lft, *cu_pbk, *cu_fsp, *lft, *fsp, *x;
  double cpu_ms, gpu_ms;
  int sample_size, i, j, h, b, k, spline_pts;
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

  MA(pbk_gpu, sample_size+1, double *);
  CA(pbk_gpu[0], (sample_size+1)*(sample_size+1), double);
  MA(pbk_cpu, sample_size+1, double *);
  CA(pbk_cpu[0], (sample_size+1)*(sample_size+1), double);
  for(i=1; i<=sample_size; i++) {
    pjh_gpu[i] = pjh_gpu[i-1] + (sample_size+1);
    pbk_gpu[i] = pbk_gpu[i-1] + (sample_size+1);
    pjh_cpu[i] = pjh_cpu[i-1] + (sample_size+1);
    pbk_cpu[i] = pbk_cpu[i-1] + (sample_size+1);
  }
  
  fprintf(stderr,"Computing pbk[][] for sample size %d -   ", sample_size);

  fprintf(stderr,"GPU: ");
  gettimeofday(&stopwatch_gpu, NULL);
  cuda_calculate_pbk(pjh_gpu[0], pbk_gpu[0], &cu_pbk, &cu_fsp, fsp, cu_lft, sample_size);
  gpu_ms = elapsed_time_ms(&stopwatch_gpu);
  fprintf(stderr," %1.1f ms", gpu_ms);

  fprintf(stderr,"\tCPU: ");
  gettimeofday(&stopwatch_cpu, NULL);
  cpu_calculate_pbk(pjh_cpu, pbk_cpu, fsp, lft, sample_size);
  cpu_ms = elapsed_time_ms(&stopwatch_cpu);
  fprintf(stderr," %1.1f ms\n", cpu_ms);

  double pjh_error = 0.;
  for(j=0; j<=sample_size; j++) {
    for(h=0; h<=sample_size; h++) {
      pjh_error += fabs(pjh_cpu[j][h] - pjh_gpu[j][h])/(pjh_cpu[j][h] + DBL_EPSILON);
      /*      fprintf(stdout,"%d\t%d\t%1.5f\t%1.5f\t%1.5f\t%1.5f\n", j, h, fsp[h], pjh_cpu[j][h], pjh_gpu[j][h],
      	      (pjh_cpu[j][h] - pjh_gpu[j][h])/(pjh_cpu[j][h]+DBL_EPSILON));*/
    }
  }
  pjh_error /= (sample_size+1)*(sample_size+1);

  double pbk_error = 0.;
  for(b=0; b <= sample_size; b++) {
    for(k=0; k <= sample_size; k++) {
      pbk_error += fabs(pbk_cpu[b][k] - pbk_gpu[b][k])/(pbk_cpu[b][k] + DBL_EPSILON);
      /*      fprintf(stdout,"%d\t%d\t%1.5f\t%1.5f\t%1.5f\t%1.5f%%\n", b, k, fsp[k],
	      pbk_cpu[b][k], pbk_gpu[b][k],
	      (pbk_cpu[b][k] - pbk_gpu[b][k])/(pbk_cpu[b][k]+DBL_EPSILON)*100.);*/
    }
  }
  pbk_error /= (sample_size+1)*(sample_size+1);
  fprintf(stderr,"pjh[][] average %% deviance of GPU calculation from CPU calculation: %g\n",
	  pjh_error);
  fprintf(stderr,"pbk[][] average %% deviance of GPU calculation from CPU calculation: %g\n",
	  pbk_error);
  fprintf(stderr,"GPU is %1.1fX faster than CPU\n", cpu_ms/gpu_ms);

  spline_pts = 200;
  MA(x, spline_pts + 1, double);
  MA(spf_cpu, (spline_pts + 1), double *);
  CA(spf_cpu[0], (spline_pts+1)*(sample_size+1), double);
  MA(spf_gpu, (spline_pts+1), double *);
  CA(spf_gpu[0], (spline_pts+1)*(sample_size+1), double);
  for(i=1; i <= spline_pts; i++) {
    spf_cpu[i] = spf_cpu[i-1] + (sample_size+1);
    spf_gpu[i] = spf_gpu[i-1] + (sample_size+1);
  }
  for(i=0; i <= spline_pts; i++)
    x[i] = LOG_AD_MIN + i*(LOG_AD_MAX - LOG_AD_MIN)/(double) (spline_pts+1);

  fprintf(stderr,"Computing spline pts for sample size %d -   ", sample_size);

  fprintf(stderr,"GPU: ");
  gettimeofday(&stopwatch_gpu, NULL);
  cuda_calculate_spline_pts(x, spf_gpu[0], cu_pbk, cu_fsp, cu_lft, spline_pts,sample_size); 
  gpu_ms = elapsed_time_ms(&stopwatch_gpu);
  fprintf(stderr," %1.1f ms", gpu_ms);

  fprintf(stderr,"\tCPU: ");
  gettimeofday(&stopwatch_cpu, NULL);
  cpu_calculate_spline_pts(x, spf_cpu, pbk_cpu, fsp, lft, spline_pts, sample_size);
  cpu_ms = elapsed_time_ms(&stopwatch_cpu);
  fprintf(stderr," %1.1f ms\n", cpu_ms);

  double spf_error = 0.;
  for(int f=0; f< sample_size; f++) {
    for(i=0; i <= spline_pts; i++) {
      //      fprintf(stderr,"%d\t%d\t%1.3f\t%g\t%g\n", f, i, x[i], spf_cpu[i][f], spf_gpu[i][f]);
      spf_error += fabs(spf_cpu[i][f] - spf_gpu[i][f])/(spf_cpu[i][f] + DBL_EPSILON);
    }
  }
  spf_error /= (spline_pts+1)*(sample_size+1);
  fprintf(stderr,"spf[][] average deviance of GPU calculation from CPU calculation: %g\n",
	  spf_error);
  fprintf(stderr,"GPU is %1.1fX faster than CPU\n", cpu_ms/gpu_ms);
  
  free(pjh_gpu[0]);
  free(pjh_gpu);

  free(pjh_cpu[0]);
  free(pjh_cpu);

  free(pbk_gpu[0]);
  free(pbk_gpu);

  free(pbk_cpu[0]);
  free(pbk_cpu);

  CUDA_FREE(cu_lft);
  free(lft);
  free(fsp);
  return 0;
}
#endif

