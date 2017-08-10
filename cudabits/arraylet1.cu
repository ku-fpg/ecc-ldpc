#include "common.h"

// Arraylet matrix coordinates //
extern "C" __global__ void selfProduct(float_ty* mLet, float_ty* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.z;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int kStart = threadIdx.x*(colCount/blockDim.x);
  /* int k = blockIdx.x*blockDim.x + threadIdx.x; */

  if (kStart == 0) {
    newMLet[(j*colCount) + i] = 1;
  }
  __syncthreads();

  float_ty prod = 1;
  if (offsets[(j/sz)*colCount + i] > -1) {
    for (int k = kStart; k < (threadIdx.x+1)*(colCount/blockDim.x); ++k) {
      if (k != i && offsets[(j/sz)*colCount + k] > -1) {
        prod *= mLet[(j*colCount) + k];
      }
    }

    atomicMul(&newMLet[(j*colCount) + i], prod);
  }
}

extern "C" __global__ void selfProductRows(float_ty* mLet, float_ty* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  extern __shared__ float_ty smem[];

  int i = threadIdx.x;
  int j = blockIdx.x;

  smem[i] = mLet[(j*colCount) + i];
  __syncthreads();

  if (offsets[(j/sz)*colCount + i] > -1) {
    for (int s = 1; s < blockDim.x; s *= 2) {
      if (i % (2*s) == 0 && offsets[(j/sz)*colCount + i + sz] > -1) {
        smem[i] *= smem[i + s];
      }
      __syncthreads();
    }
    newMLet[(j*colCount) + i] = smem[0];
  }
}


