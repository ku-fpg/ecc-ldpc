// #define NO_PARITY
#include "common.h"

// Arraylet matrix coordinates //
extern "C" __global__ void updateLamT(float_ty* newLam, float_ty* mLetT, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x*blockDim.x + threadIdx.x;

  if (i < colCount*sz) {
    int blockI = i / sz;
    int modI   = i % sz;

    int c    = blockI;

    for (int j = 0; j < rowCount/sz; ++j) {
      int off    = sz - offsets[(j*colCount) + blockI];
      int localR = (i + off) % sz;

      int r    = (j * sz) + localR;

      if (off > -1) {
        /* atomicAdd(&newLam[i], newMLet[r * colCount + c]); */
        newLam[i] += mLetT[c * rowCount + r];
      }
    }
  }

  /* if (lamIx > -1) { */
  /*   atomicAdd(&newLam[lamIx], newMLet[(j*colCount) + i]); */
  /* } */
}

extern "C" __global__ void selfProduct(float_ty* lam, float_ty* mLet, float_ty* newMLet, float_ty* mLetT, int rowCount, int colCount, int sz, int* offsets) {
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

    for (int k = 0; k < i; ++k) {
      r *= smem[globalIdx + k];
    }

    for (int k = i+1; k < colCount; ++k) {
      r *= smem[globalIdx + k];
    }

    double updatedVal = -2*atanh_(r);
    newMLet[(j*colCount) + i] = updatedVal;
    mLetT[(i*rowCount) + j]   = updatedVal;
  }
}

