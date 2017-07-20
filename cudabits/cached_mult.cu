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
  int shift = i*sz;
  int off   = offsets[(j/sz)*colCount + i];
  if (off > -1) {
    return shift + ((off + j) % sz);
  } else {
    return -1;
  }
}

__device__ void absMinMax(double* minR, double* maxR, double x, double y) {
  if (abs(x) < abs(y)) {
    *minR = x;
    *maxR = y;
  } else {
    *minR = y;
    *maxR = x;
  }
}

__device__ void lit(double* worstValR, double* multResultR, double x) {
  if (x >= 1) {
    *worstValR   = 1;
    *multResultR = x;
  } else {
    *worstValR   = x;
    *multResultR = 1;
  }
}

__device__ void smult(double* worstValR, double* multResultR, double a, double b, double c, double d) {
  double minOfMins, maxOfMins;
  absMinMax(&minOfMins, &maxOfMins, a, c);

  *worstValR   = minOfMins;
  *multResultR = b * maxOfMins * d;
}

__device__ double sdiv(double a, double b, double c) {
  if (a == c) {
    return b;
  } else if (c == 0) {
#if 1
    printf("c==0\t");
#endif
    return a * b;
  } else {
    return a * (b/c);
  }
}

extern "C" __global__ void tanhMatrix(double* tanhMat, double* mLet, double* lam, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x;
  int j = blockIdx.y;

  if (offsets[(j/sz)*colCount + i] > -1) {
    int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);
    double v = mLet[(j*colCount) + i];
    if (lamIx > -1) {
      tanhMat[(j*colCount) + i] = tanh(- ((lam[lamIx] - v)/2));
    } else {
      tanhMat[(j*colCount) + i] = v;
    }
  } else {
    tanhMat[(j*colCount) + i] = 0;
  }
}

extern "C" __global__ void tanhMulted(double* worstVals, double* multResults, double* tanhMat, double* lam, int rowCount, int colCount, int sz, int* offsets) {
  int j = blockIdx.y;
  bool vZero = false;

  /* if (offsets[(j/sz)*colCount + i] > -1) { */
    double worstVal = 1;
    double multResult = 0;

    for (int k = 0; k < colCount; ++k) {
      double v = tanhMat[(j*colCount) + k];

      int lamIx = lamIndex(k, j, sz, rowCount, colCount, offsets);
      /* if (lamIx > -1) { */
      if (offsets[(j/sz)*colCount + k] > -1) {
        double tWorst, tMultResult;
        /* lit(&tWorst, &tMultResult, tanh(- ((lam[lamIx] - v)/2))); */
        lit(&tWorst, &tMultResult, v);
        /* printf("tWorst=%f\t",tWorst); */

        if (v == 0) {
          vZero = true;
        }
        smult(&worstVal, &multResult, worstVal, multResult, tWorst, tMultResult);
      }
    }
    worstVals[j]   = worstVal;
    multResults[j] = multResult;
    if (!vZero) {
      printf("!vZero ");
    }
    /* printf("worstVals[%d]=%f\t", j, worstVal); */
    /* printf("multResults[%d]=%f\t", j, multResult); */

  /* } */
}

// Arraylet matrix coordinates //
extern "C" __global__ void tanhTransform(double* mLet, double* newMLet, double* lam, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x;
  int j = blockIdx.y;
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

    /* double v = tanhMat[(j*colCount) + i]; */

    /* /1* if (worstVals[j] <= tanhMultResults[j]) { *1/ */
    /*   newMLet[(j*colCount) + i] = -2*atanh_(sdiv(worstVals[j], tanhMultResults[j], v)); */
    /* /1* } *1/ */
  }
}

// Arraylet matrix coordinates //
extern "C" __global__ void updateLam(double* newLam, double* newMLet, int rowCount, int colCount, int sz, int* offsets) {
  int i = blockIdx.x;
  int j = blockIdx.y;
  /* for (int i = 0; i < colCount; ++i) { */
  /*   for (int j = 0; j < rowCount; ++j) { */
      int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);
      if (lamIx > -1) {
        /* newLam[lamIx] += newMLet[(j*colCount) + i]; */
        atomicAdd(&newLam[lamIx], newMLet[(j*colCount) + i]);
      }
    /* } */
  /* } */
}

__device__ bool hard(double v) {
  return v > 0;
}

// lam vector coordinates //
extern "C" __global__ void checkParity(int* pop, double* mLet, double* lam, int rowCount, int colCount, int sz, int* offsets) {
  int startJ = threadIdx.y;
  int j      = startJ;

  bool rowResult = false;

  for (int i = 0; i < colCount; ++i) {
    int lamIx = lamIndex(i, j, sz, rowCount, colCount, offsets);

    if (lamIx > -1) {
      rowResult = (rowResult != hard(lam[lamIx]));
    }
  }

  atomicAdd(pop, (rowResult ? 1 : 0));
}

