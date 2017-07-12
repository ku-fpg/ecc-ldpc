extern "C" __global__ void final_result() {
   //result[threadIdx.x] = orig_lam[threadIdx.x] > 0;
}


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

__device__ void sumRow(double* mLet, double* rowSums, int rowCount, int colCount, int sz, int* runningRowPop, int* runningColPop, int* offsets) {
  if (blockIdx.x >= rowCount) return;

  double rowSum = 0;

  int i = blockIdx.x;
  for (int j = 0; j < currPop(i/sz, runningRowPop); ++j) {
    int r = i%sz;
    int c = j + (i/sz)*colCount;
    rowSum += mLet[(c*sz) + r];
  }

  atomicAdd(&rowSums[blockIdx.x], rowSum);
}

__device__ void sumCol(double* mLet, double* colSums, int rowCount, int colCount, int sz, int* runningRowPop, int* runningColPop, int* offsets) {

  double colSum = 0;

  int i = blockIdx.x;
  for (int j = 0; j < colCount; ++j) {
    int r = (j + offsets[j]) % sz;
    int c = (j / sz)*colCount + i;

    colSum += mLet[(c*sz) + r];
  }

  atomicAdd(&colSums[blockIdx.x], colSum);
}

extern "C" __global__ void ldpc(double* mLet, int sz, int* runningRowPop, int* runningColPop, int* offsets, int rowCount, int colCount, int maxIterations, double* orig_lam, int orig_lam_len, int* result) {
  extern __shared__ double lamPrime[];

  result[blockIdx.x] = orig_lam[blockIdx.x] > 0;
  __syncthreads();

  int currAns = 1;
  int currAnsXored;

  int i;
  for (i = 0; i < maxIterations && currAns; ++i) {

    /* for (int j = 0; j < rowCount; ++j) { */
    /*   mLet[(((j*sz) + j) * colCount) + blockIdx.x]; */
    /* } */
    /* __syncthreads(); */

    lamPrime[blockIdx.x] = orig_lam[blockIdx.x];
    sumCol(mLet, lamPrime, rowCount, colCount, sz, runningRowPop, runningColPop, offsets);
    /* for (int j = 0; j < rowCount; ++j) { */
    /*   /1* atomicAdd(&lamPrime[blockIdx.x] *1/ */
    /*   /1*          ,mLet[indexR(j, sz, runningRowPop, runningColPop, offsets)] *1/ */
    /*   /1*          ); *1/ */
    /* } */

    currAnsXored = __syncthreads_count(result[blockIdx.x]>0) % 2;
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

