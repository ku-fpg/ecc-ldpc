#include "common.h"

// Arraylet matrix coordinates //
extern "C" __global__ void selfProduct(float_ty* mLet, float_ty* newMLet, int* nonzeroRows, int* nonzeroCols, int nonzeroCount, int rowCount, int colCount, int sz, int* offsets) {
  extern __shared__ float_ty prod[];
  /* int i = blockIdx.z; */
  /* int j = blockIdx.y*blockDim.y + threadIdx.y; */
  int j = nonzeroRows[blockIdx.x];
  int i = nonzeroCols[blockIdx.x];
  int kStart = threadIdx.x*(colCount/blockDim.x);
  int k = kStart;

  int prodIx = (threadIdx.y*colCount) + k;

  /* if (offsets[(j/sz)*colCount + i] == -1) { */
  /*   printf("== -1\t"); */
  /* } */

    if (k != i && offsets[(j/sz)*colCount + k] > -1) {
      prod[prodIx] = mLet[(j*colCount) + k];
    } else {
      prod[prodIx] = 1;
    }
    __syncthreads();

    int nearestPowerOfTwo = 32;
    double p;
    for (int s = nearestPowerOfTwo; s > 0; s >>= 1) {

      /* if (k + s < colCount) { */
      /*   p = prod[prodIx + s]; */
      /* } */
      /* __syncthreads(); */

      if (k + s < colCount) {
        prod[prodIx] *= prod[prodIx + s];
      }
      __syncthreads();
    }

    if (threadIdx.x == 0) {
      newMLet[(j*colCount) + i] = prod[threadIdx.y*colCount];
    }
  /* } */
}

extern "C" __global__ void selfProductRows(float_ty* mLet, float_ty* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  extern __shared__ float_ty smem[];

  int i = threadIdx.x;
  int j = blockIdx.x;
  /* int k = (threadIdx.x + 1)%blockDim.x; */

  smem[i] = mLet[(j*colCount) + i];
  __syncthreads();

  if (offsets[(j/sz)*colCount + i] > -1) {
    for (int s = 1; s < blockDim.x; s *= 2) {
      if (i % (2*s) == 0 && offsets[(j/sz)*colCount + i + sz] > -1) {
        smem[i] *= smem[i + s];
      }
      __syncthreads();
    }
    /* for (int s = blockDim.x/2; s > 0; s >>= 1) { */
    /*   if (i < s) { */
    /*     /1* smem[i] = smem[(i+1)%blockDim.x]*smem[i+s]; *1/ */
    /*     smem[i] *= smem[i+s]; */
    /*   } */
    /*   __syncthreads(); */
    /* } */
    newMLet[(j*colCount) + i] = smem[0]; // /mLet[(j*colCount) + i];
  }

  /* if (threadIdx.x == 0) { */
  /* atomicMul(&newMLet[(j*colCount) + i], smem[0]/mLet[(j*colCount) + i]); */
  /* } */
}

