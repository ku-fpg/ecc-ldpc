extern "C" __global__ void final_result() {
   //result[threadIdx.x] = orig_lam[threadIdx.x] > 0;
}


#if 0
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
#endif

__device__ int currPop(int ix, int* runningPopulation) {
  return runningPopulation[ix+1] - runningPopulation[ix];
}

__device__ int skipAmt(int ix, int sz, int* runningPopulation) {
  return sz*runningPopulation[ix];
}

__device__ int indexR(int r, int sz, int* runningRowPop, int* runningColPop, int* offsets) {
  int blockR = r / sz;
  return skipAmt(blockR, sz, runningRowPop) + ((r + offsets[r]) % sz);
}

extern "C" __device__ void sumRow(double* mLet, double* initVals, double* rowSums, int rowCount, int colCount, int sz, int* runningRowPop, int* runningColPop, int* offsets) {
  int iPrime = blockDim.x * blockIdx.x + threadIdx.x;
  /* int iPrime = (i * sz) + */ 

  if (iPrime < rowCount) {
    rowSums[iPrime] = initVals[iPrime];
    for (int j = 0; j < colCount; ++j) {
      /* printf("(%d, %d)  ", j, iPrime); */
      rowSums[iPrime] += mLet[(j*sz) + iPrime];
    }
    /* printf("%f ", rowSums[iPrime]); */
  }
}

__global__ void sumCol(double* mLet, double* colSums, int rowCount, int colCount, int sz, int* runningRowPop, int* runningColPop, int* offsets) {
}

extern "C" __global__ void ldpc(double* mLet, int sz, int* runningRowPop, int* runningColPop, int* offsets, int rowCount, int colCount, int maxIterations, double* orig_lam, int orig_lam_len, int* result) {
  extern __shared__ double lamPrime[];

  result[blockIdx.x]   = orig_lam[blockIdx.x] > 0;
  lamPrime[blockIdx.x] = orig_lam[blockIdx.x];
  __syncthreads();

  int currAns = 1;
  int currAnsXored;

  int i;
  int iPrime = blockDim.x * blockIdx.x + threadIdx.x;
  for (i = 0; i < maxIterations && currAns; ++i) {
    sumRow(mLet, orig_lam, lamPrime, rowCount, colCount, sz, runningRowPop, runningColPop, offsets);


    currAnsXored = __syncthreads_count(lamPrime[iPrime]>0) % 2;
    currAns      = __syncthreads_or(currAnsXored);
    /* __syncthreads(); // TODO: See if this is necessary */
  }

  __syncthreads();

  if (!currAns) {
    result[blockIdx.x] = orig_lam[blockIdx.x] > 0;
  } else {
    result[blockIdx.x] = lamPrime[blockIdx.x] > 0;
  }
}

