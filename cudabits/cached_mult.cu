extern "C" __global__ void ldpc(double* mLet, double* offsets, int rowCount, int colCount, int maxIterations, double* orig_lam, int orig_lam_len, int* result) {
  /* for (int i = 0; i < orig_lam_len; ++i) { */
  /*   result[i] = 0; */
  /* } */

  int iters;
  for (iters = 0; iters < maxIterations; ++iters) {

  }

  if (iters > 0) {
    for (int i = 0; i < orig_lam_len; ++i) {
      result[i] = orig_lam[i] > 0;
    }
  }
}

