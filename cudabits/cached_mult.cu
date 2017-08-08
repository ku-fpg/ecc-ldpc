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
  int j = threadIdx.y;

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

// Arraylet matrix coordinates //
extern "C" __global__ void selfProduct(float_ty* mLet, float_ty* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  extern __shared__ float_ty prod[];
  int i = blockIdx.z;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int kStart = threadIdx.x; // *(colCount/blockDim.x);
  int k = kStart;

  if (offsets[(j/sz)*colCount + i] > -1) {
    if (k != i && offsets[(j/sz)*colCount + k] > -1) {
      prod[k] = mLet[(j*colCount) + k];
    } else {
      prod[k] = 1;
    }

    __syncthreads();

    int nearestPowerOfTwo = 32;
    if (k < blockDim.x/2) {
      for (int s = nearestPowerOfTwo; s > 0; s >>= 1) {
        if (k + s < blockDim.x) {
          prod[k] *= prod[k + s];
        }
        __syncthreads();
      }
    }

    if (k == 0) {
      newMLet[(j*colCount) + i] = prod[0];
    }
  }
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

extern "C" __global__ void parityRowResults(bool* rowResults, float_ty* lam, int rowCount, int colCount, int sz, int* offsets) {
  int startI = threadIdx.x;
  int i      = startI;
  int j      = blockIdx.y;

  int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);

  int count = __syncthreads_count(lamIx != -1 && hard(lam[lamIx]));

  if (i == 0) {
    rowResults[j] = count % 2 == 1;
  }
}

// lam vector coordinates //
extern "C" __global__ void checkParity(int* pop, bool* rowResults) {
  int startJ = threadIdx.y;
  int j      = startJ;

  *pop = __syncthreads_or(rowResults[j]);
}

