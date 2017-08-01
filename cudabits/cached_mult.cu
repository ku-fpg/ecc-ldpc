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


__device__ double signum(double x) {
  if (x < 0) {
    return -1;
  } else if (x > 0) {
    return 1;
  } else {
    return 0;
  }
}

__device__ double atanh_(double x) {
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

extern "C" __global__ void tanhTransform(double* mLet, double* lam, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x;
  int j = threadIdx.y;

  int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);
  if (lamIx > -1) {
    double v = mLet[(j*colCount) + i];
    mLet[(j*colCount) + i] = tanh(- ((lam[lamIx] - v)/2));
  }
}

extern "C" __global__ void setToOne(double* mLet, int rowCount, int colCount, int sz, int* offsets) {
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

extern "C" __global__ void insertOnes(double* mLet, bool* nonzero, int rowCount, int colCount, int sz) {
  int i = blockIdx.x;
  int j = threadIdx.y;

  if (!nonzero[(j*colCount) + i]) {
    mLet[(j*colCount) + i] = 1;
  }
}

// Arraylet matrix coordinates //
extern "C" __global__ void selfProduct(double* mLet, double* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.z;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int kStart = threadIdx.x*(colCount/blockDim.x);


  double prod = 1;
  if (offsets[(j/sz)*colCount + i] > -1) {
    for (int k = kStart; k < (threadIdx.x+1)*(colCount/blockDim.x); ++k) {

      if (k != i && offsets[(j/sz)*colCount + k] > -1) {
        prod *= mLet[(j*colCount) + k];
      }
    }

    atomicMul(&newMLet[(j*colCount) + i], prod);
  }
}

extern "C" __global__ void atanhTransform(double* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x;
  int j = threadIdx.y;

  if (offsets[(j/sz)*colCount + i] > -1) {
    newMLet[(j*colCount) + i] = -2*atanh_(newMLet[(j*colCount)+i]);
  } else {
    newMLet[(j*colCount) + i] = 0;
  }
}

// Arraylet matrix coordinates //
extern "C" __global__ void updateLam(double* newLam, double* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x;
  int j = threadIdx.y;
  int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);

  if (lamIx > -1) {
    atomicAdd(&newLam[lamIx], newMLet[(j*colCount) + i]);
  }
}

__device__ bool hard(double v) {
  return v > 0;
}

extern "C" __global__ void parityRowResults(bool* rowResults, double* lam, int rowCount, int colCount, int sz, int* offsets) {
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

