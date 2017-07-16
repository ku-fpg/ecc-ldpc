// From http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#axzz4meEZrFDA
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}

__device__ void sumRow(double* mLet, double* initVals, double* rowSums, int rowCount, int colCount, int sz, int* runningRowPop, int* runningColPop, int* offsets) {
  int r = blockDim.x * blockIdx.x + threadIdx.x;
  int c = blockDim.y * blockIdx.y + threadIdx.y;

  atomicAdd(&rowSums[c], mLet[(c*rowCount) + r]);

  /* int iPrime = blockDim.x * blockIdx.x + threadIdx.x; */
  /* int iPrime = (i * sz) + */ 

  /* if (iPrime < rowCount) { */
  /*   rowSums[iPrime] = initVals[iPrime]; */
  /*   for (int j = 0; j < colCount; ++j) { */
  /*     /1* if (offsets[j] == 0) break; *1/ // TODO: See if this is more or less efficient */

  /*     /1* printf("(%d, %d)  ", j, iPrime); *1/ */
  /*     rowSums[iPrime] += mLet[(j*sz) + iPrime]; */
  /*   } */
  /*   /1* printf("%f ", rowSums[iPrime]); *1/ */
  /* } */
}

__device__ void sumCol(double* mLet, double* colSums, int rowCount, int colCount, int sz, int* offsets) {
}

__device__ double signum(double x) {
  return x < 0 ? -1 : 1;
}

// Arraylet matrix coordinates //
extern "C" __global__ void tanhTransform(double* mLet, double* newMLet, double* lam, int rowCount, int colCount, int sz, int* offsets) {
  int r = blockDim.x * blockIdx.x + threadIdx.x;
  int c = blockDim.y * blockIdx.y + threadIdx.y;

  int stride = rowCount;

  double prod = 1;

  if (offsets[c] != 0) {
    for (int j = 0; j < colCount; ++j) {
      double v = mLet[j*stride + r];
      if (j != c) {
        prod *= tanh(- ((lam[j] - v)/2));
      }
    }
  }

  newMLet[c*stride + r] = -2 * atanh(prod);
}

// Arraylet matrix coordinates //
extern "C" __global__ void updateLam(double* newLam, double* newMLet, int rowCount, int colCount, int sz) {
  int r = blockDim.x * blockIdx.x + threadIdx.x;
  int c = blockDim.y * blockIdx.y + threadIdx.y;

  atomicAdd(&newLam[r], newMLet[c*rowCount + r]);
}

__device__ bool hard(double v) {
  return v > 0 ? 1 : 0;
}

// lam vector coordinates //
extern "C" __global__ void checkParity(bool* result, double* newMLet, double* lamPrime, int rowCount, int colCount) {
  extern __shared__ int xored[];
  /* int r = blockDim.x * blockIdx.x + threadIdx.x; */
  int c = blockDim.x * blockIdx.x + threadIdx.x;

  bool hardLam = hard(lamPrime[c]);

  xored[c] = 0;

  for (int i = 0; i < colCount; ++i) {
    xored[c] = xored[c] != (hardLam != hard(newMLet[c*rowCount + i]));
  }

  /* int currAnsXored = __syncthreads_count(lamPrime[iPrime]>0) % 2; */
  result[0] = __syncthreads_or(xored[c] ? 1 : 0);
}



/* extern "C" __global__ void ldpc(double* mLet, int sz, int* runningRowPop, int* runningColPop, int* offsets, int rowCount, int colCount, int maxIterations, double* orig_lam, int orig_lam_len, int* result) { */
/*   extern __shared__ double lamPrime[]; */
/*   /1* __shared__ double mLetPrime[1408*32]; *1/ */

/*   result[blockIdx.x]   = orig_lam[blockIdx.x] > 0; */
/*   lamPrime[blockIdx.x] = orig_lam[blockIdx.x]; */
/*   __syncthreads(); */

/*   int currAns = 1; */
/*   int currAnsXored; */

/*   int i; */
/*   int iPrime = blockDim.x * blockIdx.x + threadIdx.x; */
/*   for (i = 0; i < maxIterations && !currAns; ++i) { */
/*     tanhTransform(mLet, lamPrime, rowCount, colCount, sz, runningRowPop, runningColPop, offsets); */
/*     __syncthreads(); */

/*     sumRow(mLet, orig_lam, lamPrime, rowCount, colCount, sz, runningRowPop, runningColPop, offsets); */
/*     __syncthreads(); */

/*     currAnsXored = __syncthreads_count(lamPrime[iPrime]>0) % 2; */
/*     currAns      = __syncthreads_or(currAnsXored); */
/*     __syncthreads(); // TODO: See if this is necessary */
/*   } */

/*   __syncthreads(); */

/*   if (currAns) { */
/*     result[blockIdx.x] = orig_lam[blockIdx.x] > 0; */
/*   } else { */
/*     result[blockIdx.x] = lamPrime[blockIdx.x] > 0; */
/*   } */
/* } */

