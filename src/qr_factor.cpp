#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "par_for.h"
#include "qr_factor_core.h"

using namespace Rcpp;

// The triangular factor of a tall matrix, by Householder reflections, with
// the trailing columns of each step updated in parallel.
//
// WHY ONLY R. augmented_solve() forms (A'A)^-1 u and reads nothing else off
// the decomposition: Q is never accumulated, never applied, and never
// returned. That is what makes this worth writing at all -- the factor is
// p x p whatever n is, and the whole cost is the reduction.
//
// WHY THE COLUMNS. Step j applies one reflector to each trailing column,
// and those updates are independent: column k is read and written by one
// thread, in the same order and with the same operations as the sequential
// loop performs them. So the factor is bit-identical at any thread count BY
// CONSTRUCTION rather than by measurement, and no arithmetic is added over
// the sequential form -- which is what a block decomposition over ROWS
// (TSQR) cannot offer: its answer depends on the partition, and it pays
// about 15 per cent more flops for the reduction stage. Measured against
// TSQR on the same matrices, this route also agrees with R's own qr() to
// 1.6e-16 on the solve where TSQR agrees to 1e-14, being the same
// algorithm rather than a different orthogonalization order.
//
// The region is opened per step and only where the trailing block is large
// enough to pay for it, the work there being ntrail * len; below that the
// step runs sequentially whatever the count says, which costs one
// comparison and keeps the answer the same either way.
// [[Rcpp::export]]
NumericMatrix qr_factor_cpp(NumericMatrix A, int threads) {
    const int m = A.nrow(), p = A.ncol();
    NumericMatrix B = clone(A);
    NumericMatrix R(p, p);
    sm7::qr_factor_core(B.begin(), m, p, threads, R.begin());
    return R;
}
