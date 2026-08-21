#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "par_for.h"

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
    double* Bp = B.begin();
    std::vector<double> v(m > 0 ? m : 1);
    const int kmax = (m < p) ? m : p;

    // measured on this machine: the per-step region costs ~44 us to open,
    // and a trailing update runs at roughly one element per nanosecond, so
    // a step is worth splitting from about a quarter of a million elements
    const double kMinWork = 2.6e5;

    for (int j = 0; j < kmax; ++j) {
        const int len = m - j;
        double* col = Bp + (std::size_t) j * m + j;
        double nrm = 0.0;
        for (int i = 0; i < len; ++i) nrm += col[i] * col[i];
        nrm = std::sqrt(nrm);
        if (!R_FINITE(nrm) || nrm == 0.0) continue;
        const double alpha = (col[0] >= 0.0) ? -nrm : nrm;
        for (int i = 0; i < len; ++i) v[i] = col[i];
        v[0] -= alpha;
        double vtv = 0.0;
        for (int i = 0; i < len; ++i) vtv += v[i] * v[i];
        if (vtv == 0.0) continue;
        col[0] = alpha;
        for (int i = 1; i < len; ++i) col[i] = 0.0;

        const int ntrail = p - j - 1;
        if (ntrail <= 0) continue;
        const double* vp = v.data();
        const double inv = 2.0 / vtv;
        auto body = [&](std::size_t t) {
            const int k = j + 1 + (int) t;
            double* ck = Bp + (std::size_t) k * m + j;
            double s = 0.0;
            for (int i = 0; i < len; ++i) s += vp[i] * ck[i];
            s *= inv;
            for (int i = 0; i < len; ++i) ck[i] -= s * vp[i];
        };
        const double work = (double) ntrail * (double) len;
        sm7::par_for((std::size_t) ntrail, (work >= kMinWork) ? threads : 1,
                     body);
    }

    NumericMatrix R(p, p);
    for (int j = 0; j < p; ++j) {
        for (int i = 0; i <= j && i < kmax; ++i) {
            R(i, j) = Bp[(std::size_t) j * m + i];
        }
    }
    return R;
}
