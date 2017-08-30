// #define NO_PARITY
#include "common.h"

// extern "C" __global__ void parityRowResults(bool* rowResults, float_ty* lam, int rowCount, int colCount, int sz, int* offsets) {
//   int startI = threadIdx.x;
//   int i      = startI;
//   int j      = blockIdx.y;

//   int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);

//   int count = __syncthreads_count(lamIx != -1 && hard(lam[lamIx]));

//   if (i == 0) {
//     rowResults[j] = count % 2 == 1;
//   }
// }

// // lam vector coordinates //
// extern "C" __global__ void checkParity(int* pop, bool* rowResults) {
//   int startJ = threadIdx.y;
//   int j      = startJ;

//   *pop = __syncthreads_or(rowResults[j]);
// }

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

extern "C" __global__ void selfProduct(float_ty* lam, float_ty* mLet, float_ty* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  extern __shared__ float_ty smem[];

  int i = blockIdx.y*blockDim.y + threadIdx.y;
  int j = blockIdx.x*blockDim.x + threadIdx.x;

  int globalIdx = threadIdx.x*colCount;

  double r = 1;

  double v = mLet[(j*colCount) + i];

  int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);

  if (lamIx > -1) {
    smem[globalIdx + i] = tanh(- ((lam[lamIx] - v)/2));
  } else {
    smem[globalIdx + i] = 1;
  }

  __syncthreads();

  if (lamIx > -1) {

    /* for (int k = 0; k < colCount; ++k) { */
    /*   if (k != i) { */
    /*     r *= smem[globalIdx + k]; */
    /*   } */
    /* } */

    for (int k = 0; k < i; ++k) {
      r *= smem[globalIdx + k];
    }

    for (int k = i+1; k < colCount; ++k) {
      r *= smem[globalIdx + k];
    }

    newMLet[(j*colCount) + i] = -2*atanh_(r);
  }
}

