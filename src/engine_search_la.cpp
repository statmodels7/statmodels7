// joint_decrement() of the compiled outer search (engine_search.cpp), kept
// apart from RcppArmadillo: see engine_search_la.h.
// R evaluates these expressions one operation at a time; clang (arm64
// macOS) would otherwise fuse a*b + c into one rounding
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
#include <cfloat>
#include <vector>
#include "engine_search_la.h"
#define USE_FC_LEN_T
#define R_NO_REMAP
#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>

#ifndef FCONE
#define FCONE
#endif

// joint_decrement(g, A): eigen() with only.values (dsyevr) and solve()
// (dgesv with R's rcond test) as R calls them
double joint_decrement(const std::vector<double>& g, const std::vector<double>& A) {
    const int n = (int) g.size();
    if (n == 0) return 0.0;
    for (double x : g) if (!R_FINITE(x)) return NA_REAL;
    for (double x : A) if (!R_FINITE(x)) return NA_REAL;
    {
        // La_rs(x, only.values = TRUE)
        std::vector<double> a(A), w(n);
        char jobv = 'N', range = 'A', uplo = 'L';
        double vl = 0.0, vu = 0.0, abstol = 0.0;
        int il = 0, iu = 0, m = 0, info = 0, lwork = -1, liwork = -1, itmp = 0;
        double tmp = 0.0;
        std::vector<int> isuppz(2 * (size_t) n);
        double* z = nullptr;
        int ldz = 1;
        F77_CALL(dsyevr)(&jobv, &range, &uplo, &n, a.data(), &n, &vl, &vu, &il, &iu,
                         &abstol, &m, w.data(), z, &ldz, isuppz.data(), &tmp, &lwork,
                         &itmp, &liwork, &info FCONE FCONE FCONE);
        lwork = (int) tmp;
        liwork = itmp;
        std::vector<double> work(lwork);
        std::vector<int> iwork(liwork);
        F77_CALL(dsyevr)(&jobv, &range, &uplo, &n, a.data(), &n, &vl, &vu, &il, &iu,
                         &abstol, &m, w.data(), z, &ldz, isuppz.data(), work.data(),
                         &lwork, iwork.data(), &liwork, &info FCONE FCONE FCONE);
        if (info != 0) return NA_REAL;
        for (int i = 0; i < n; ++i) if (!(w[i] > 0.0)) return NA_REAL;
    }
    std::vector<double> v(g);
    {
        // La_solve(a, b, tol = .Machine$double.eps)
        std::vector<double> a(A);
        std::vector<int> ipiv(n);
        int nrhs = 1, info = 0;
        char one = '1';
        const double anorm = F77_CALL(dlange)(&one, &n, &n, A.data(), &n, nullptr FCONE);
        F77_CALL(dgesv)(&n, &nrhs, a.data(), &n, ipiv.data(), v.data(), &n, &info);
        if (info != 0) return NA_REAL;
        double rcond = 0.0;
        std::vector<double> work(4 * (size_t) n);
        std::vector<int> iwork(n);
        F77_CALL(dgecon)(&one, &n, a.data(), &n, &anorm, &rcond, work.data(), iwork.data(),
                         &info FCONE);
        if (rcond < DBL_EPSILON) return NA_REAL;
    }
    for (double x : v) if (!R_FINITE(x)) return NA_REAL;
    long double s = 0.0L;
    for (int i = 0; i < n; ++i) s += (long double) (g[i] * v[i]);
    return 0.5 * (double) s;
}
