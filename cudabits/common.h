typedef double float_ty;

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

__device__ double atomicMul(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val *
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
__device__ float atomicMul(float* address, float val)
{
    int* address_as_i = (int*)address;
    int old = *address_as_i, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_i, assumed,
                        __float_as_int(val *
                               __int_as_float(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __int_as_float(old);
}
__device__ double atomicAssign(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, val);

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}


__device__ float_ty signum(float_ty x) {
  if (x < 0) {
    return -1;
  } else if (x > 0) {
    return 1;
  } else {
    return 0;
  }
}

__device__ float_ty atanh_(float_ty x) {
  if (x == 1 || x == -1) {
    return signum(x) * 18.714973875118524;
  } else {
    return atanh(x);
  }
}

__device__ int lamIndex(int i, int j, int sz, int rowCount, int colCount, int* offsets) {
  int shift = i*sz;
  int off   = offsets[(j/sz)*colCount + i];
  if (off > -1) {
    return shift + ((off + j) % sz);
  } else {
    return -1;
  }
}

extern "C" __global__ void tanhTransform(float_ty* mLet, float_ty* lam, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;

  int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);
  if (lamIx > -1) {
    float_ty v = mLet[(j*colCount) + i];
    mLet[(j*colCount) + i] = tanh(- ((lam[lamIx] - v)/2));
  }
}

extern "C" __global__ void setToOne(float_ty* mLet, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x;
  int j = threadIdx.y;

  int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);
  /* if (lamIx > -1) { */
    mLet[(j*colCount) + i] = 1;
  /* } */
}

extern "C" __global__ void makeNonzeroMat(bool* nonzero, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x;
  int j = threadIdx.y;

  nonzero[(j*colCount) + i] = (offsets[(j/sz)*colCount + i] > -1);
}

extern "C" __global__ void insertOnes(float_ty* mLet, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x;
  int j = threadIdx.y;

  if (offsets[((j/sz)*colCount) + i] == -1) {
    mLet[(j*colCount) + i] = 1;
  }
}

extern "C" __global__ void atanhTransform(float_ty* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  /* int i = blockIdx.x; */
  /* int j = threadIdx.y; */
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;

  if (offsets[(j/sz)*colCount + i] > -1) {
    newMLet[(j*colCount) + i] = -2*atanh_(newMLet[(j*colCount)+i]);
  } else {
    newMLet[(j*colCount) + i] = 0;
  }
}

// Arraylet matrix coordinates //
extern "C" __global__ void updateLam(float_ty* newLam, float_ty* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);

  if (lamIx > -1) {
    atomicAdd(&newLam[lamIx], newMLet[(j*colCount) + i]);
  }
}

__device__ bool hard(float_ty v) {
  return v > 0;
}

#ifndef NO_PARITY

/* extern "C" __global__ void parityRowResults(int* rowResults, float_ty* lam, int rowCount, int colCount, int sz, int* offsets) { */
/*   extern __shared__ int paritySmem[]; */

/*   int i = blockIdx.x*blockDim.x + threadIdx.x; */
/*   int j = blockIdx.y*blockDim.y + threadIdx.y; */

/*   int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets); */

/*   // rowResults[j] = __syncthreads_count(lamIx != -1 && hard(lam[lamIx])); */

/*   paritySmem[(threadIdx.y*colCount) + i] = lamIx > -1 && hard(lam[lamIx]); */
/*   __syncthreads(); */

/*   int nearestPowerOfTwo = 32; */
/*   for (int s = nearestPowerOfTwo; s > 0; s >>= 1) { */
/*     // if (s < colCount) { */
/*       paritySmem[(threadIdx.y*colCount) + i] ^= paritySmem[(threadIdx.y*colCount) + i + s]; */
/*     // } */
/*     __syncthreads(); */
/*   } */

/*   if (i == 0) { */
/*     rowResults[j] = paritySmem[threadIdx.y*colCount]; */
/*   } */
/* } */

extern "C" __global__ void parityRowResults(int* done, float_ty* lam, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int jStart = blockIdx.y*blockDim.y + threadIdx.y;
  int j = jStart;

  if (i == 0 && j == 0) {
    *done = 0;
  }
  __syncthreads();

  if (!*done) {
    int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);

    int count = __syncthreads_count(lamIx > -1 && hard(lam[lamIx]));

    if (threadIdx.x == 0) {
      if (count % 2 == 1 == 1) {
        *done = 1;
      }
    }
  }
}


// lam vector coordinates //
extern "C" __global__ void checkParity(int* pop, int* rowResults) {
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  if (j == 0) {
    /* atomicAnd(pop, 0); */
    *pop = 0;
  }
  __syncthreads();

  int blockOr = __syncthreads_or(rowResults[j] % 2) != 0 ? 1 : 0;
  atomicOr(pop, blockOr);
}
#endif

