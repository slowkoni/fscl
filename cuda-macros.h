#ifndef CUDA_MACROS_H
#define CUDA_MACROS_H

static void cuda_exit(int exit_code) {
  cudaError_t cerr;

  if (exit_code) {
    fprintf(stderr,"Resetting CUDA device on abnormal termination (%d)\n", 
	    exit_code);
  }
  cerr = cudaDeviceReset();
  if (cerr != cudaSuccess) {
    fprintf(stderr, "Failed to reset CUDA device (%s)\n", 
	    cudaGetErrorString(cerr));
  }
  exit(exit_code);
}

#undef MA
#define MA(p,s,t)       {			\
    (p) = (t*) malloc(sizeof(t)*(s));		\
    if ((p) == NULL) {						   \
      fprintf(stderr,"Failed allocating memory at %s:%d (%1.1f Mb)\n",	\
	  __FILE__,__LINE__, (s)/1e6);				    \
      abort();							    \
    }								    \
}

#undef CA
#define CA(p,s,t)	 {		\
    (p) = (t*) malloc(sizeof(t)*s);	\
if ((p) == NULL) {  \
  fprintf(stderr,"Failed allocating memory at %s:%d (%1.1f Mb)\n", \
	  __FILE__,__LINE__, (s)/1e6);				    \
  abort();							    \
  }\
    memset((p), 0x0, (s)); \
}

#undef RA
#define RA(p,s,t)           {			\
    (p) = (t *) realloc((p),(s)*sizeof(t));	\
if ((p) == NULL) {  \
  fprintf(stderr,"Failed allocating memory at %s:%d (%1.1f Mb)\n", \
	  __FILE__,__LINE__, (s)/1e6);				    \
  abort();							    \
 }                                \
}

#define CUDA_MA(p, s) { cudaError_t cerr; cerr = cudaMalloc(&(p), (s)); if (cerr != cudaSuccess) { fprintf(stderr,"Can't allocate %d KB CUDA device memory at line %d (%s)\n", (int) ((s)/1024), __LINE__, cudaGetErrorString(cerr)); cuda_exit(-1); } }

#define CUDA_MA_HOST_MAP(cp, p, s) { cudaError_t cerr; cerr = cudaHostAlloc(&(p), (s), cudaHostAllocMapped); if (cerr != cudaSuccess) { fprintf(stderr,"Can't allocate %d KB host pinned memory at line %d (%s)\n", (int) ((s)/1024), __LINE__, cudaGetErrorString(cerr)); cuda_exit(-1); } cerr = cudaHostGetDevicePointer(&(cp), (p), 0); if (cerr != cudaSuccess) { fprintf(stderr,"Can't retrieve device pointer for host allocated pinned memory at line %d (%s)\n", __LINE__, cudaGetErrorString(cerr)); cuda_exit(-1); } }


#define CUDA_MA_HOST(p, s, flags) { cudaError_t cerr; cerr = cudaHostAlloc(&(p), (s), (flags)); if (cerr != cudaSuccess) { fprintf(stderr,"Can't allocate %d KB host (cudaHostAlloc()) memory at line %d (%s)\n", (int) ((s)/1024), __LINE__, cudaGetErrorString(cerr)); cuda_exit(-1); } }

#define CUDA_MEMCPY_TO(cp, p, s) { cudaError_t cerr; cerr = cudaMemcpy((cp), (p), (s), cudaMemcpyHostToDevice); if (cerr != cudaSuccess) { fprintf(stderr,"Can't copy memory to CUDA device at line %d (%s)\n", __LINE__, cudaGetErrorString(cerr)); cuda_exit(-1); } }

#define CUDA_MEMCPY_TOA(cp, p, s, stream) { cudaError_t cerr; cerr = cudaMemcpyAsync((cp), (p), (s), cudaMemcpyHostToDevice, stream); if (cerr != cudaSuccess) { fprintf(stderr,"Can't copy memory to CUDA device at line %d (%s)\n", __LINE__, cudaGetErrorString(cerr)); cuda_exit(-1); } }

#define CUDA_MEMCPY_FROM(p, cp, s) { cudaError_t cerr; cerr = cudaMemcpy((p), (cp), (s), cudaMemcpyDeviceToHost); if (cerr != cudaSuccess) { fprintf(stderr,"Can't copy memory from CUDA device at line %d (%s)\n", __LINE__, cudaGetErrorString(cerr)); cuda_exit(-1); } }

#define CUDA_MEMCPY_FROMA(p, cp, s, stream) { cudaError_t cerr; cerr = cudaMemcpyAsync((p), (cp), (s), cudaMemcpyDeviceToHost, stream); if (cerr != cudaSuccess) { fprintf(stderr,"Can't copy memory from CUDA device at line %d (%s)\n", __LINE__, cudaGetErrorString(cerr)); cuda_exit(-1); } }

#define CUDA_FREE(cp) { cudaError_t cerr; cerr = cudaFree((cp)); if (cerr != cudaSuccess) { fprintf(stderr,"Can't free CUDA device memory at line %d (%s)\n", __LINE__, cudaGetErrorString(cerr)); cuda_exit(-1); } }

#define CUDA_FREE_HOST(cp) { cudaError_t cerr; cerr = cudaFreeHost((cp)); if (cerr != cudaSuccess) { fprintf(stderr,"Can't free CUDA device memory at line %d (%s)\n", __LINE__, cudaGetErrorString(cerr)); cuda_exit(-1); } }

#define CUDA_MMAP(cp, p, s) { cudaError_t cerr; cerr = cudaHostRegister((p), (s), cudaHostRegisterMapped); if (cerr != cudaSuccess) { fprintf(stderr, "Can't register host memory for mapping to CUDA decvice at line %d (%s)\n", __LINE__, cudaGetErrorString(cerr)); cuda_exit(-1); } cudaHostGetDevicePointer(&(cp), (p), 0); if (cp == NULL) { fprintf(stderr, "NULL pointer returned from cudaHostGetDevicePointer() at line %d\n", __LINE__); cuda_exit(-1); } }



#endif
