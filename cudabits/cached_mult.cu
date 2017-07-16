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

  if (prod == 1 || prod == -1) {
    newMLet[c*stride + r] = signum(prod) * 18.714973875118524;
  } else {
    newMLet[c*stride + r] = -2 * atanh(prod);
  }
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

