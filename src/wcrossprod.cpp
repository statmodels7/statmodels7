#include <Rcpp.h>
#include "par_for.h"
using namespace Rcpp;

// The weighted cross product of the assembly, C = A' diag(w) B, decomposed
// over the ELEMENTS of C: element (j, k) is one dot product, accumulated in
// full by one thread in ascending row order with the weight multiplied onto
// A's entry first -- the same operations in the same order as
// crossprod(A * w, B) against the reference BLAS, whose dgemm accumulates a
// scalar in ascending order too. That is what licenses the bit-identity
// guarantee of n_threads(): no reduction is split across threads
// (piano_parallel.txt, section 0, measured difference exactly zero).
//
// The R wrapper (wcrossprod) holds the gate: only base dense matrices with
// a full-length weight reach this, and only above the measured work
// threshold; a sparse design keeps the Matrix route, where a threaded dense
// kernel would do nothing.

// X'v, one dot product per output element, ascending rows: the same
// operations in the same order as crossprod(X, v) against the reference
// BLAS. The caller precomputes v (e.g. w * r) exactly as the expression it
// replaces did, so the elementwise rounding is untouched.
// [[Rcpp::export]]
NumericVector xtv_cpp(NumericMatrix X, NumericVector v, int threads) {
    const int n = X.nrow(), p = X.ncol();
    NumericVector out(p);
    const double* Xp = X.begin();
    const double* vp = v.begin();
    double* op = out.begin();
    sm7::par_for(static_cast<std::size_t>(p), threads, [&](std::size_t j) {
        const double* x = Xp + j * n;
        double acc = 0.0;
        for (int i = 0; i < n; ++i) acc += x[i] * vp[i];
        op[j] = acc;
    });
    return out;
}

// w'(X * X), the column curvatures of a working model: element j is
// sum_i w_i * (x_ij * x_ij), with the square rounded first -- the same
// arithmetic as crossprod(w, X^2), without materializing the n x p square.
// [[Rcpp::export]]
NumericVector wxsq_cpp(NumericMatrix X, NumericVector w, int threads) {
    const int n = X.nrow(), p = X.ncol();
    NumericVector out(p);
    const double* Xp = X.begin();
    const double* wp = w.begin();
    double* op = out.begin();
    sm7::par_for(static_cast<std::size_t>(p), threads, [&](std::size_t j) {
        const double* x = Xp + j * n;
        double acc = 0.0;
        for (int i = 0; i < n; ++i) acc += wp[i] * (x[i] * x[i]);
        op[j] = acc;
    });
    return out;
}

// A'A, one dot product per element of the UPPER triangle, ascending rows,
// the lower triangle written as its mirror: half the work, which is what
// the reference dsyrk behind crossprod(A) also does -- a full-p^2 version
// measured no better than dsyrk at 8 workers. Each (j, k) pair is owned by
// exactly one task (the one holding its column k), so the mirrored write
// is disjoint too, and the mirror copies a value that equals the (k, j)
// accumulation exactly, the products and their order being the same.
// [[Rcpp::export]]
NumericMatrix xtx_cpp(NumericMatrix A, int threads) {
    const int n = A.nrow(), p = A.ncol();
    NumericMatrix C(p, p);
    const double* Ap = A.begin();
    double* Cp = C.begin();
    sm7::par_for(static_cast<std::size_t>(p), threads, [&](std::size_t k) {
        const double* b = Ap + k * n;
        for (std::size_t j = 0; j <= k; ++j) {
            const double* a = Ap + j * n;
            double acc = 0.0;
            for (int i = 0; i < n; ++i) acc += a[i] * b[i];
            Cp[k * p + j] = acc;
            if (j != k) Cp[j * p + k] = acc;
        }
    });
    return C;
}

// [[Rcpp::export]]
NumericMatrix wcrossprod_cpp(NumericMatrix A, NumericVector w,
                             NumericMatrix B, int threads) {
    const int n = A.nrow(), pa = A.ncol(), pb = B.ncol();
    NumericMatrix C(pa, pb);
    const double* Ap = A.begin();
    const double* Bp = B.begin();
    const double* wp = w.begin();
    double* Cp = C.begin();

    sm7::par_for(static_cast<std::size_t>(pa) * pb, threads,
                 [&](std::size_t t) {
        const std::size_t k = t / pa, j = t - k * pa;
        const double* a = Ap + j * n;
        const double* b = Bp + k * n;
        double acc = 0.0;
        for (int i = 0; i < n; ++i) acc += (a[i] * wp[i]) * b[i];
        Cp[k * pa + j] = acc;
    });
    return C;
}
