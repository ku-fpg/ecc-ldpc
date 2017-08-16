#include "common.h"


// From
// http://on-demand.gputechconf.com/gtc/2013/presentations/S3174-Kepler-Shuffle-Tips-Tricks.pdf
__device__ __inline__ double shfl(double x, int lane) {
  // Split the double number into 2 32b registers.
  int lo, hi;
  asm volatile( "mov.b32 {%0,%1}, %2;" : "=r"(lo), "=r"(hi) : "d"(x));

  // Shuffle the two 32b registers.
  lo = __shfl(lo, lane);
  hi = __shfl(hi, lane);

  // Recreate the 64b number.
  asm volatile( "mov.b64 %0, {%1,%2};" : "=d"(x) : "r"(lo), "r"(hi));

  return x;
}

extern "C" __global__ void selfProduct(float_ty* mLet, float_ty* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  extern __shared__ double smem[];
  // int i = blockIdx.z;
  // int j = blockIdx.y*blockDim.y + threadIdx.y;
  int i = blockIdx.y*blockDim.y + threadIdx.y;
  int j = blockIdx.x*blockDim.x + threadIdx.x;

  // if (threadIdx.y == 0) {
  //   newMLet[(j*colCount) + i] = 1;
  // }

  double r = 1;

  double v = mLet[(j*colCount) + i];

  if (offsets[((j/sz)*colCount) + i] > -1) {
    smem[(threadIdx.x*colCount) + i] = v;
  } else {
    smem[(threadIdx.x*colCount) + i] = 1;
  }

  __syncthreads();

  if (offsets[((j/sz)*colCount) + i] > -1) {
    for (int k = 0; k < i; ++k) {
      r *= smem[(threadIdx.x*colCount) + k];//__shfl(orig, k);
    }

    for (int k = i+1; k < colCount; ++k) {
      r *= smem[(threadIdx.x*colCount) + k];//__shfl(orig, k);
    }

    newMLet[(j*colCount) + i] = r;
  } else {
    newMLet[(j*colCount) + i] = 1;
  }
}

