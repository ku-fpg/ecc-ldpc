extern "C" __global__ void ldpc(double* mLet, double* offsets, int rowCount, int colCount, int maxIterations, double* orig_lam, int* result) {
  for (int i = 0; i < colCount; ++i) {
    result[i] = 0;
  }
}

