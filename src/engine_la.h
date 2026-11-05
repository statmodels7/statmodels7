// Dense linear algebra of the compiled inner fit (piano_motore.md, phase 2),
// written as the R functions that iwls_fit() calls, through the same BLAS,
// LAPACK and LINPACK routines with the same arguments, so that the two routes
// produce the same numbers on the same matrix: qr() is dqrdc2 (qr.default),
// chol() is dpotrf on the upper triangle with the lower one zeroed (La_chol),
// eigen(symmetric = TRUE) is dsyevr on the lower triangle with abstol 0 and
// its output reversed (La_rs and eigen()), svd(nu = 0) is dgesdd with jobz
// "S" (La.svd), backsolve() and forwardsolve() are dtrsm on a copy of the
// right-hand side (do_backsolve). The sums R takes with sum(), colSums() and
// mean() accumulate in long double, as R's own do on x86-64 (D13); products
// that R sends to the BLAS (%*%, crossprod) accumulate in double in
// ascending order, as the reference BLAS does.
#ifndef SM7_ENGINE_LA_H
#define SM7_ENGINE_LA_H

#include <cfloat>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <vector>
#include <cmath>
#include <algorithm>

#ifndef FCONE
#define FCONE
#endif

namespace sm7e {

typedef long double ldouble;

// a dense column-major matrix
struct Dense {
    int nr = 0, nc = 0;
    std::vector<double> v;
    Dense() {}
    Dense(int r, int c) : nr(r), nc(c), v((std::size_t) r * c, 0.0) {}
    double& operator()(int i, int j) { return v[(std::size_t) j * nr + i]; }
    double operator()(int i, int j) const { return v[(std::size_t) j * nr + i]; }
    double* col(int j) { return v.data() + (std::size_t) j * nr; }
    const double* col(int j) const { return v.data() + (std::size_t) j * nr; }
};

// sum() of a double vector: long double accumulation, then R's conversion
inline double r_sum(const double* x, std::size_t n) {
    ldouble s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += x[i];
    if (s > DBL_MAX) return R_PosInf;
    if (s < -DBL_MAX) return R_NegInf;
    return (double) s;
}

// colSums(): long double per column (do_colsum)
inline double r_colsum(const double* x, std::size_t n) {
    ldouble s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += x[i];
    return (double) s;
}

// mean() of a double vector: R's two passes in long double (real_mean)
inline double r_mean(const double* x, std::size_t n) {
    ldouble s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += x[i];
    s /= n;
    if (R_FINITE((double) s)) {
        ldouble t = 0.0;
        for (std::size_t i = 0; i < n; ++i) t += (x[i] - s);
        s += t / n;
    }
    return (double) s;
}

// stats::median() of a finite vector
inline double r_median(std::vector<double> x) {
    const std::size_t n = x.size();
    if (n == 0) return NA_REAL;
    for (double v : x) if (ISNAN(v)) return NA_REAL;
    std::sort(x.begin(), x.end());
    const std::size_t half = (n + 1) / 2;
    if (n % 2 == 1) return x[half - 1];
    double two[2] = {x[half - 1], x[half]};
    return r_mean(two, 2);
}

// crossprod(X, v) for a dense X: one ascending dot product per column
inline void r_crossprod_vec(const double* X, int n, int m, const double* v,
                            double* out) {
    for (int j = 0; j < m; ++j) {
        const double* xj = X + (std::size_t) j * n;
        double t = 0.0;
        for (int i = 0; i < n; ++i) t = t + xj[i] * v[i];
        out[j] = t;
    }
}

// X %*% b for a dense X (dgemv "N"): y_i accumulated over the columns in
// ascending order
inline void r_matvec(const double* X, int n, int m, const double* b,
                     double* y) {
    for (int i = 0; i < n; ++i) y[i] = 0.0;
    for (int j = 0; j < m; ++j) {
        const double t = b[j];
        const double* xj = X + (std::size_t) j * n;
        for (int i = 0; i < n; ++i) y[i] = y[i] + t * xj[i];
    }
}

// do_backsolve(): solves with the k x k triangle of r, upper or lower, on a
// copy of x. Returns false where R signals "singular matrix in 'backsolve'".
inline bool r_trisolve(const Dense& r, bool upper, std::vector<double>& x) {
    const int k = r.nc;
    if (k == 0) return true;
    for (int i = 0; i < k; ++i) if (r(i, i) == 0.0) return false;
    const int one_i = 1;
    const double one = 1.0;
    int nrr = r.nr;
    int kk = k;
    F77_CALL(dtrsm)("L", upper ? "U" : "L", "N", "N", &kk, &one_i, &one,
                    r.v.data(), &nrr, x.data(), &kk FCONE FCONE FCONE FCONE);
    return true;
}

inline Dense transpose(const Dense& a) {
    Dense t(a.nc, a.nr);
    for (int j = 0; j < a.nc; ++j)
        for (int i = 0; i < a.nr; ++i) t(j, i) = a(i, j);
    return t;
}

// chol(A): the upper factor; false where R signals that a leading minor is
// not positive
inline bool r_chol(const Dense& A, Dense& U) {
    const int n = A.nr;
    U = A;
    for (int j = 0; j < n; ++j)
        for (int i = j + 1; i < n; ++i) U(i, j) = 0.0;
    int info = 0, m = n;
    F77_CALL(dpotrf)("U", &m, U.v.data(), &m, &info FCONE);
    return info == 0;
}

// backsolve(ch, forwardsolve(t(ch), u))
inline bool r_chol_solve(const Dense& U, std::vector<double>& u) {
    Dense L = transpose(U);
    if (!r_trisolve(L, false, u)) return false;
    return r_trisolve(U, true, u);
}

// eigen(A, symmetric = TRUE): values decreasing, vectors in the same order
inline bool r_eigen_sym(const Dense& A, bool only_values,
                        std::vector<double>& values, Dense& vectors) {
    int n = A.nr;
    values.assign(n, 0.0);
    if (n == 0) return true;
    std::vector<double> rx(A.v);
    std::vector<double> rv(n);
    std::vector<double> rz(only_values ? 1 : (std::size_t) n * n);
    std::vector<int> isuppz(2 * (std::size_t) n);
    const char* jobv = only_values ? "N" : "V";
    double vl = 0.0, vu = 0.0, abstol = 0.0, tmp = 0.0;
    int il = 0, iu = 0, m = 0, info = 0, lwork = -1, liwork = -1, itmp = 0;
    F77_CALL(dsyevr)(jobv, "A", "L", &n, rx.data(), &n, &vl, &vu, &il, &iu,
                     &abstol, &m, rv.data(), rz.data(), &n, isuppz.data(),
                     &tmp, &lwork, &itmp, &liwork, &info FCONE FCONE FCONE);
    if (info != 0) return false;
    lwork = (int) tmp;
    liwork = itmp;
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    F77_CALL(dsyevr)(jobv, "A", "L", &n, rx.data(), &n, &vl, &vu, &il, &iu,
                     &abstol, &m, rv.data(), rz.data(), &n, isuppz.data(),
                     work.data(), &lwork, iwork.data(), &liwork, &info
                     FCONE FCONE FCONE);
    if (info != 0) return false;
    for (int i = 0; i < n; ++i) values[i] = rv[n - 1 - i];
    if (!only_values) {
        vectors = Dense(n, n);
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                vectors(i, j) = rz[(std::size_t) (n - 1 - j) * n + i];
    }
    return true;
}

// svd(A, nu = 0): singular values d and v = t(vt)
inline bool r_svd(const Dense& A, std::vector<double>& d, Dense& v) {
    int n = A.nr, p = A.nc;
    const int np = std::min(n, p);
    std::vector<double> x(A.v);
    d.assign(np, 0.0);
    std::vector<double> u((std::size_t) n * np), vt((std::size_t) np * p);
    std::vector<int> iwork(8 * (std::size_t) np);
    int ldu = n, ldvt = np, lwork = -1, info = 0;
    double tmp = 0.0;
    F77_CALL(dgesdd)("S", &n, &p, x.data(), &n, d.data(), u.data(), &ldu,
                     vt.data(), &ldvt, &tmp, &lwork, iwork.data(), &info FCONE);
    if (info != 0) return false;
    lwork = (int) tmp;
    std::vector<double> work(lwork);
    F77_CALL(dgesdd)("S", &n, &p, x.data(), &n, d.data(), u.data(), &ldu,
                     vt.data(), &ldvt, work.data(), &lwork, iwork.data(), &info
                     FCONE);
    if (info != 0) return false;
    v = Dense(p, np);
    for (int j = 0; j < np; ++j)
        for (int i = 0; i < p; ++i) v(i, j) = vt[(std::size_t) i * np + j];
    return true;
}

// qr(A, tol): dqrdc2 as qr.default calls it; on return qr holds the
// compact factor, rank, qraux and the 1-based pivot
struct QRResult {
    Dense qr;
    int rank = 0;
    std::vector<double> qraux;
    std::vector<int> pivot;
};

inline bool r_qr(const Dense& A, double tol, QRResult& out) {
    int n = A.nr, p = A.nc;
    if (1.0 * n * p > 2147483647.0) return false;
    out.qr = A;
    out.qraux.assign(p, 0.0);
    out.pivot.resize(p);
    for (int j = 0; j < p; ++j) out.pivot[j] = j + 1;
    std::vector<double> work(2 * (std::size_t) std::max(p, 1));
    int rank = 0, ldx = n;
    F77_CALL(dqrdc2)(out.qr.v.data(), &ldx, &n, &p, &tol, &rank,
                     out.qraux.data(), out.pivot.data(), work.data());
    out.rank = rank;
    return true;
}

}  // namespace sm7e

#endif
