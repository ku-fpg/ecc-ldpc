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

__device__ double atanh_(double x) {
  if (x == 1 || x == -1) {
    return signum(x) * 18.714973875118524;
  } else {
    return atanh(x);
  }
}

__device__ int lamIndex(int i, int j, int sz, int rowCount, int colCount, int* offsets) {
  int shift = (i/sz)*sz;
  int off   = offsets[(j/sz)*colCount + i];
  if (off > -1) {
    return shift + ((off + j) % sz);
  } else {
    return -1;
  }
}

// Arraylet matrix coordinates //
extern "C" __global__ void tanhTransform(double* mLet, double* newMLet, double* lam, int rowCount, int colCount, int sz, int* offsets) {
  for (int i = 0; i < colCount; ++i) {
    for (int j = 0; j < rowCount; ++j) {
      if (offsets[(j/sz)*colCount + i] > -1) {
        double prod = 1;

        for (int k = 0; k < colCount; ++k) {
          double v = mLet[(j*colCount) + k];

          int lamIx = lamIndex(k, j, sz, rowCount, colCount, offsets);
          if (k != i) {
            if (lamIx > -1) {
              prod *= tanh(- ((lam[lamIx] - v)/2));
            }
          }
        }

        newMLet[(j*colCount) + i] = -2*atanh_(prod);
      }
    }
  }
}

// Arraylet matrix coordinates //
extern "C" __global__ void updateLam(double* newLam, double* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  /* printf("(rowCount,colCount)=(%d,%d)\t",rowCount, colCount*sz); */
    for (int i = 0; i < colCount; ++i) {
  for (int j = 0; j < rowCount; ++j) {
      newLam[lamIndex(i, j, sz, rowCount, colCount, offsets)] += newMLet[(j*colCount) + i];
    }
  }
}

__device__ bool hard(double v) {
  return v > 0;
}

// lam vector coordinates //
extern "C" __global__ void checkParity(bool* result, double* mLet, double* lam, int rowCount, int colCount, int sz, int* offsets) {
  *result = false;

    for (int j = 0; j < rowCount; ++j) {
    bool rowResult = false;
  for (int i = 0; i < colCount; ++i) {
      int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);

      if (lamIx > -1) {
        rowResult = (rowResult != hard(lam[lamIx]));
      }
    }
    if (rowResult) {
      *result = true;
      break;
    }
  }
}

