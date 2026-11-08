// The four covariance charts of the compiled engine at orders 0, 1 and 2
// (SPEC_penalty2 section 3, with the corrections of REVIEW_penalty2 and of
// REVIEW_penalty 3.1), written to reproduce R's numbers bit for bit:
//
//   code 1  log_cholesky(p)                       both orientations
//   code 2  compound_symmetry(p)                  both orientations
//   code 3  ar1(p)                                both orientations
//   code 4  dr_prod(p) with correlation_matrix(p) both orientations
//
// One function per derivative order: chart_order0(), chart_order1() and
// chart_order2() each fill their own components of a ChartPt, and each
// order reads the state left by the order below. A chart of K free values
// keeps its second-order components for the K (K + 1) / 2 pairs k <= l in
// the order of numericals7::tuple_indices(K, 2) (pair_index()); the pair
// (t, s) with t > s is the pair (s, t).
//
// Arithmetic follows R's operation order: x^2 is x * x and any other power
// R_pow(); sum(), rowSums() and colSums() accumulate in long double; %*% is
// the BLAS dgemm, tcrossprod() dsyrk, forwardsolve()/backsolve() dtrsm,
// chol() dpotrf and eigen(only.values = TRUE) dsyevr. exp() takes R's
// vector wrapper (1 + x for |x| <= sqrt(eps)) where R calls it on a vector
// of length two or more, and libm where R calls it on one number
// (SPEC_penalty2 R2', REVIEW_penalty2 findings 2 and 4).
#ifndef CHART7_ENGINE_CHARTS_H
#define CHART7_ENGINE_CHARTS_H

#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <Rmath.h>
#include <vector>
#include <string>
#include <cmath>
#include <cfloat>
#include <algorithm>

#ifndef FCONE
#define FCONE
#endif

namespace chart7 {

typedef long double ldouble;

// ---------------------------------------------------------------- dense
// a column-major matrix (engine_la.h's Dense)
struct Dense {
    int nr = 0, nc = 0;
    std::vector<double> v;
    Dense() {}
    Dense(int r, int c) : nr(r), nc(c), v((std::size_t) r * c, 0.0) {}
    double& operator()(int i, int j) { return v[(std::size_t) j * nr + i]; }
    double operator()(int i, int j) const { return v[(std::size_t) j * nr + i]; }
};

inline Dense identity(int p) {
    Dense I(p, p);
    for (int i = 0; i < p; ++i) I(i, i) = 1.0;
    return I;
}

inline Dense transpose(const Dense& a) {
    Dense t(a.nc, a.nr);
    for (int j = 0; j < a.nc; ++j)
        for (int i = 0; i < a.nr; ++i) t(j, i) = a(i, j);
    return t;
}

// A + B, A - B, -A, elementwise
inline Dense add(const Dense& a, const Dense& b) {
    Dense o(a.nr, a.nc);
    for (std::size_t i = 0; i < o.v.size(); ++i) o.v[i] = a.v[i] + b.v[i];
    return o;
}
inline Dense sub(const Dense& a, const Dense& b) {
    Dense o(a.nr, a.nc);
    for (std::size_t i = 0; i < o.v.size(); ++i) o.v[i] = a.v[i] - b.v[i];
    return o;
}
inline Dense neg(const Dense& a) {
    Dense o(a.nr, a.nc);
    for (std::size_t i = 0; i < o.v.size(); ++i) o.v[i] = -a.v[i];
    return o;
}
// x * A, the scalar on the left as R writes it
inline Dense scal(double x, const Dense& a) {
    Dense o(a.nr, a.nc);
    for (std::size_t i = 0; i < o.v.size(); ++i) o.v[i] = x * a.v[i];
    return o;
}
// A * B elementwise
inline Dense hadamard(const Dense& a, const Dense& b) {
    Dense o(a.nr, a.nc);
    for (std::size_t i = 0; i < o.v.size(); ++i) o.v[i] = a.v[i] * b.v[i];
    return o;
}

// A %*% B: dgemm (R's matprod(); with the reference BLAS its dgemv calls
// for a single row or column give the same sums)
inline Dense mm(const Dense& A, const Dense& B) {
    Dense C(A.nr, B.nc);
    if (!A.nr || !B.nc || !A.nc) return C;
    const double one = 1.0, zero = 0.0;
    int m = A.nr, n = B.nc, k = A.nc;
    F77_CALL(dgemm)("N", "N", &m, &n, &k, &one, A.v.data(), &m, B.v.data(), &k,
                    &zero, C.v.data(), &m FCONE FCONE);
    return C;
}

// tcrossprod(A): dsyrk on the upper triangle, copied to the lower
// (symtcrossprod())
inline Dense tcross(const Dense& A) {
    int n = A.nr, k = A.nc;
    Dense C(n, n);
    if (!n) return C;
    const double one = 1.0, zero = 0.0;
    F77_CALL(dsyrk)("U", "N", &n, &k, &one, A.v.data(), &n, &zero, C.v.data(), &n
                    FCONE FCONE);
    for (int i = 1; i < n; ++i)
        for (int j = 0; j < i; ++j) C(i, j) = C(j, i);
    return C;
}

// backsolve()/forwardsolve() of the k x k triangle of r on B, in place
// (do_backsolve()); false where R signals "singular matrix in 'backsolve'"
inline bool trsolve(const Dense& r, bool upper, Dense& B) {
    const int k = r.nc;
    if (k == 0) return true;
    for (int i = 0; i < k; ++i) if (r(i, i) == 0.0) return false;
    const double one = 1.0;
    int nrr = r.nr, kk = k, nb = B.nc;
    if (nb == 0) return true;
    F77_CALL(dtrsm)("L", upper ? "U" : "L", "N", "N", &kk, &nb, &one,
                    r.v.data(), &nrr, B.v.data(), &kk FCONE FCONE FCONE FCONE);
    return true;
}

// chol(A): dpotrf on the upper triangle with the lower one zeroed (La_chol)
inline bool r_chol(const Dense& A, Dense& U) {
    const int n = A.nr;
    U = A;
    for (int j = 0; j < n; ++j)
        for (int i = j + 1; i < n; ++i) U(i, j) = 0.0;
    int info = 0, m = n;
    F77_CALL(dpotrf)("U", &m, U.v.data(), &m, &info FCONE);
    return info == 0;
}

// eigen(A, symmetric = TRUE, only.values = TRUE)$values: dsyevr on the
// lower triangle, abstol 0, reversed to decreasing order (La_rs)
inline bool r_eigen_values(const Dense& A, std::vector<double>& values) {
    int n = A.nr;
    values.assign(n, 0.0);
    if (n == 0) return true;
    std::vector<double> rx(A.v), rv(n), rz(1);
    std::vector<int> isuppz(2 * (std::size_t) n);
    double vl = 0.0, vu = 0.0, abstol = 0.0, tmp = 0.0;
    int il = 0, iu = 0, m = 0, info = 0, lwork = -1, liwork = -1, itmp = 0;
    F77_CALL(dsyevr)("N", "A", "L", &n, rx.data(), &n, &vl, &vu, &il, &iu,
                     &abstol, &m, rv.data(), rz.data(), &n, isuppz.data(),
                     &tmp, &lwork, &itmp, &liwork, &info FCONE FCONE FCONE);
    if (info != 0) return false;
    lwork = (int) tmp;
    liwork = itmp;
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    F77_CALL(dsyevr)("N", "A", "L", &n, rx.data(), &n, &vl, &vu, &il, &iu,
                     &abstol, &m, rv.data(), rz.data(), &n, isuppz.data(),
                     work.data(), &lwork, iwork.data(), &liwork, &info
                     FCONE FCONE FCONE);
    if (info != 0) return false;
    for (int i = 0; i < n; ++i) values[i] = rv[n - 1 - i];
    return true;
}

// ---------------------------------------------------------------- scalars
// sum(): long double, then R's conversion (rsum())
inline double r_sum(const double* x, std::size_t n) {
    ldouble s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += x[i];
    if (s > DBL_MAX) return R_PosInf;
    if (s < -DBL_MAX) return R_NegInf;
    return (double) s;
}
// colSums() / rowSums() of one line: long double, plain conversion
inline double r_linesum(const double* x, std::size_t n) {
    ldouble s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += x[i];
    return (double) s;
}

// R's x^y: x * x at y = 2, R_pow() otherwise (R_POW() of arithmetic.c;
// R_pow() is R's own, powl() on 64-bit Windows)
inline double rpow(double x, double y) { return y == 2.0 ? x * x : R_pow(x, y); }

// R's vector wrappers of arithmetic.c (d7x::r_exp ... of distributions7)
inline double w_exp(double x) { return (std::fabs(x) <= std::sqrt(DBL_EPSILON)) ? 1. + x : std::exp(x); }
inline double w_expm1(double x) { return (std::fabs(x) <= DBL_EPSILON) ? x : std::expm1(x); }
inline double w_log1p(double x) { return (std::fabs(x) <= DBL_EPSILON) ? x : std::log1p(x); }

// exp_floored()'s floor, (24 / .Machine$double.xmax)^0.25 as R computes it
inline double exp_floor() {
    static const double f = R_pow(24.0 / DBL_MAX, 0.25);
    return f;
}

// ---------------------------------------------------------------- links
// linkfunctions7's registry, resolved once
typedef void (*lf7_inv12p_t)(int, const double*, double, double*, double*, double*);
typedef double (*lf7_clamp_t)(double, double, double);
typedef int (*lf7_class_id_t)(const char*);
typedef double (*lf7_fun1p_t)(int, const double*, double);

struct Lf7 {
    lf7_inv12p_t inv12p = nullptr;
    lf7_clamp_t clamp = nullptr;
    lf7_class_id_t class_id = nullptr;
    lf7_fun1p_t fun1p = nullptr;
};
inline const Lf7& lf7() {
    static Lf7 f;
    if (!f.inv12p) {
        f.inv12p = (lf7_inv12p_t) R_GetCCallable("linkfunctions7", "lf7_inv12p");
        f.clamp = (lf7_clamp_t) R_GetCCallable("linkfunctions7", "lf7_clamp");
        f.class_id = (lf7_class_id_t) R_GetCCallable("linkfunctions7", "lf7_class_id");
        f.fun1p = (lf7_fun1p_t) R_GetCCallable("linkfunctions7", "lf7_fun1p");
    }
    return f;
}

// a link: its registry id (lf7_class_id() of the class name), its own
// parameters (link_scalar_route()$par) and its bounds (link_bounds)
struct LinkSpec {
    int id = -1;
    std::vector<double> par;
    double lwr = R_NegInf, upr = R_PosInf;
};

// linkinv() (clamped), dlinkinv() and d2linkinv() at one value, as R
// computes them on ONE number (libm inside lf7_inv12p())
inline void link_inv3(const LinkSpec& lk, double eta, double& h0, double& h1, double& h2) {
    const Lf7& f = lf7();
    double h = 0.0;
    f.inv12p(lk.id, lk.par.empty() ? nullptr : lk.par.data(), eta, &h, &h1, &h2);
    h0 = f.clamp(h, lk.lwr, lk.upr);
}

// ---------------------------------------------------------------- charts
enum { CH_LC = 1, CH_CS = 2, CH_AR1 = 3, CH_DRP = 4 };

// the static description of a chart, built once in R
struct ChartStatic {
    int code = 0, p = 0, K = 0;
    bool inverted = false;           // a precision chart (kind 4 only)
    // log-Cholesky: the positions of the free values, 0-based (diagonal
    // first, then the strict lower triangle column by column)
    std::vector<int> row, col;
    std::vector<char> ondiag;
    // compound symmetry and AR(1): the scale's and the correlation's link;
    // dr_prod: the scales' link (log only) and the angles' link
    // (bounded_link(0, pi), DoublyBoundedLink only)
    LinkSpec ls, lr;
    // dr_prod: each angle's row and column of the correlation factor,
    // 0-based, and the angle at (row, column) (-1 where none)
    std::vector<int> arow, acol;
    std::vector<int> aidx;           // p x p, column-major
};

inline int n_pairs(int K) { return K * (K + 1) / 2; }
// the position of the pair (k, l), k <= l, in numericals7::tuple_indices(K, 2)
inline int pair_index(int K, int k, int l) {
    if (k > l) std::swap(k, l);
    if (k == l) return k;
    return K + k * (2 * K - k - 1) / 2 + (l - k - 1);
}

// one chart at one eta. The oriented fields (L, sinv, logdet, ai, dld, ai2,
// d2ld) are what a consumer reads; the others are the chart's own state.
struct ChartPt {
    int order = -1;
    bool ok = false;                 // the chart's values (order 0) are finite and
                                     // R signals no error at them
    bool L_ok = false;               // param_factor() succeeded (the eigenvalue test
                                     // of compound symmetry and AR(1), a zero diagonal)
    std::string err;
    // oriented: order 0
    Dense L, sinv;
    double logdet = 0.0;
    // oriented: order 1 (K) and order 2 (K (K + 1) / 2 pairs)
    std::vector<Dense> ai, ai2;
    std::vector<double> dld, d2ld;
    // the chart: Sigma = M (or Omega for a precision chart), its raw inverse
    // si (param_solve()), its log-determinant and its derivatives
    std::vector<double> eta;
    Dense M, si;
    double ld = 0.0;
    std::vector<Dense> a, a2;        // param_d1(), param_d2()
    std::vector<double> rdld, rd2ld; // param_dlogdet(), param_d2logdet()
    std::vector<Dense> rai, rai2;    // param_inv_d1(), param_inv_d2()
    // log-Cholesky: G = L^-1, Lk, B_k = G Lk, C_k = B_k + t(B_k)
    Dense G;
    std::vector<Dense> Lk, B, C;
    // compound symmetry and AR(1)
    double v0 = 0, v1 = 0, v2 = 0, r0 = 0, r1 = 0, r2 = 0, om = 0, pa = 0, pb = 0;
    Dense P0, P1, P2;
    // dr_prod: the scales' sd[0..2] (d = sd0), the angles' tables
    std::vector<double> sd0, sd1, sd2;
    std::vector<double> s0, c0, sn1, sn2, cs1, cs2;
    Dense L0, R0;
};

// -------------------------------------------------- the sandwich of
// matrix_parameter (inv_derivs.R), on the UNSYMMETRIZED si
inline Dense sandwich1(const Dense& si, const Dense& a) {
    return neg(mm(mm(si, a), si));
}
inline Dense sandwich2(const Dense& si, const Dense& ak, const Dense& al, const Dense& a2) {
    Dense t = sub(add(mm(mm(al, si), ak), mm(mm(ak, si), al)), a2);
    return mm(mm(si, t), si);
}

// ======================================================== log-Cholesky
// chol_leibniz_cpp() of parameters7 on one tuple (order 1 or 2): the
// derivative of M = L L' by the Leibniz rule
inline Dense lc_leibniz(const ChartStatic& st, const Dense& L, const int* tup, int order) {
    const int p = st.p;
    Dense m(p, p);
    const int nsub = 1 << order;
    for (int b = 0; b < nsub; ++b) {
        int ks = -1, cntS = 0;
        bool ok = true;
        for (int j = 0; j < order; ++j) {
            if (!(b & (1 << j))) continue;
            int k = tup[j];
            if (cntS == 0) ks = k;
            else if (k != ks) { ok = false; break; }
            ++cntS;
        }
        if (!ok || (cntS > 1 && !st.ondiag[ks])) continue;
        int kc = -1, cntC = 0;
        for (int j = 0; j < order; ++j) {
            if (b & (1 << j)) continue;
            int k = tup[j];
            if (cntC == 0) kc = k;
            else if (k != kc) { ok = false; break; }
            ++cntC;
        }
        if (!ok || (cntC > 1 && !st.ondiag[kc])) continue;
        if (cntS == 0 && cntC > 0) {
            int r = st.row[kc], c = st.col[kc];
            double v = st.ondiag[kc] ? L(r, c) : 1.0;
            for (int ii = 0; ii < p; ++ii) m(ii, r) += v * L(ii, c);
        } else if (cntC == 0 && cntS > 0) {
            int r = st.row[ks], c = st.col[ks];
            double v = st.ondiag[ks] ? L(r, c) : 1.0;
            for (int jj = 0; jj < p; ++jj) m(r, jj) += v * L(jj, c);
        } else if (cntS > 0 && cntC > 0) {
            int ra = st.row[ks], ca = st.col[ks];
            int rb = st.row[kc], cb = st.col[kc];
            if (ca == cb) {
                double va = st.ondiag[ks] ? L(ra, ca) : 1.0;
                double vb = st.ondiag[kc] ? L(rb, cb) : 1.0;
                m(ra, rb) += va * vb;
            }
        }
    }
    return m;
}

inline bool lc_d0(const ChartStatic& st, ChartPt& pt) {
    const int p = st.p, K = st.K;
    const std::vector<double>& eta = pt.eta;
    // chol_assemble(): exp() on the vector of the p diagonal values, R's
    // wrapper for p >= 2 and libm for p = 1
    pt.L = Dense(p, p);
    std::vector<double> dg;
    for (int k = 0; k < K; ++k) {
        double v = eta[k];
        if (st.ondiag[k]) {
            v = (p >= 2) ? w_exp(eta[k]) : std::exp(eta[k]);
            dg.push_back(eta[k]);
        }
        pt.L(st.row[k], st.col[k]) = v;
    }
    pt.L_ok = true;
    for (int i = 0; i < p; ++i) if (pt.L(i, i) == 0.0) pt.L_ok = false;
    pt.ld = 2 * r_sum(dg.data(), dg.size());
    // param_solve(): backsolve(t(L), forwardsolve(L, I))
    pt.si = identity(p);
    if (!trsolve(pt.L, false, pt.si) || !trsolve(transpose(pt.L), true, pt.si)) {
        pt.err = "singular matrix in 'backsolve'";
        return false;
    }
    if (st.inverted) pt.M = tcross(pt.L);
    return true;
}

inline void lc_d1(const ChartStatic& st, ChartPt& pt) {
    const int p = st.p, K = st.K;
    if (st.inverted) {
        pt.a.resize(K);
        for (int k = 0; k < K; ++k) pt.a[k] = lc_leibniz(st, pt.L, &k, 1);
    } else {
        pt.G = identity(p);
        trsolve(pt.L, false, pt.G);
        Dense Gt = transpose(pt.G);
        pt.Lk.resize(K); pt.B.resize(K); pt.C.resize(K); pt.rai.resize(K);
        for (int k = 0; k < K; ++k) {
            Dense lk(p, p);
            lk(st.row[k], st.col[k]) = st.ondiag[k] ? pt.L(st.row[k], st.row[k]) : 1.0;
            pt.Lk[k] = lk;
            pt.B[k] = mm(pt.G, lk);
            pt.C[k] = add(pt.B[k], transpose(pt.B[k]));
            pt.rai[k] = neg(mm(mm(Gt, pt.C[k]), pt.G));
        }
    }
    pt.rdld.assign(K, 0.0);
    for (int k = 0; k < K; ++k) pt.rdld[k] = st.ondiag[k] ? 2.0 : 0.0;
}

inline void lc_d2(const ChartStatic& st, ChartPt& pt) {
    const int K = st.K, np = n_pairs(K);
    pt.rd2ld.assign(np, 0.0);
    if (st.inverted) {
        pt.a2.resize(np);
        for (int k = 0; k < K; ++k)
            for (int l = k; l < K; ++l) {
                int tup[2] = {k, l};
                pt.a2[pair_index(K, k, l)] = lc_leibniz(st, pt.L, tup, 2);
            }
        return;
    }
    Dense Gt = transpose(pt.G);
    pt.rai2.resize(np);
    for (int k = 0; k < K; ++k)
        for (int l = k; l < K; ++l) {
            Dense dbk = mm(neg(pt.B[l]), pt.B[k]);
            if (k == l && st.ondiag[k]) dbk = add(dbk, mm(pt.G, pt.Lk[k]));
            Dense dc = add(dbk, transpose(dbk));
            Dense inner = sub(add(mm(transpose(pt.B[l]), pt.C[k]), mm(pt.C[k], pt.B[l])), dc);
            pt.rai2[pair_index(K, k, l)] = mm(mm(Gt, inner), pt.G);
        }
}

// ======================================== compound symmetry and AR(1)
// chol_pd(): the finite test, the eigenvalue test, then t(chol())
inline bool chol_pd(const Dense& m, Dense& L) {
    for (double x : m.v) if (!R_FINITE(x)) return false;
    std::vector<double> ev;
    if (!r_eigen_values(m, ev)) return false;
    const double mx = ev.front(), mn = ev.back();
    if (mx <= 0 || mn <= 1e-14 * mx) return false;
    Dense U;
    if (!r_chol(m, U)) return false;
    L = transpose(U);
    return true;
}

inline bool econ_d0(const ChartStatic& st, ChartPt& pt) {
    const int p = st.p;
    const bool ar = st.code == CH_AR1;
    link_inv3(st.ls, pt.eta[0], pt.v0, pt.v1, pt.v2);
    link_inv3(st.lr, pt.eta[1], pt.r0, pt.r1, pt.r2);
    const double r0 = pt.r0, r1 = pt.r1, r2 = pt.r2, v0 = pt.v0;
    pt.P0 = Dense(p, p); pt.P1 = Dense(p, p); pt.P2 = Dense(p, p);
    if (!ar) {
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < p; ++i) {
                const double off = (i == j) ? 0.0 : 1.0, id = (i == j) ? 1.0 : 0.0;
                pt.P0(i, j) = id + r0 * off;
                pt.P1(i, j) = r1 * off;
                pt.P2(i, j) = r2 * off;
            }
    } else {
        std::vector<double> w0(p), w1(p), w2(p);
        for (int m = 0; m < p; ++m) {
            w0[m] = rpow(r0, (double) m);
            const double d1 = (m >= 1) ? (double) m * rpow(r0, (double) (m - 1)) : 0.0;
            const double d2 = (m >= 2) ? (double) (m * (m - 1)) * rpow(r0, (double) (m - 2)) : 0.0;
            w1[m] = d1 * r1;
            w2[m] = (d2 * (r1 * r1)) + (d1 * r2);
        }
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < p; ++i) {
                const int lag = std::abs(i - j);
                pt.P0(i, j) = w0[lag]; pt.P1(i, j) = w1[lag]; pt.P2(i, j) = w2[lag];
            }
    }
    pt.M = scal(v0, pt.P0);
    pt.L_ok = chol_pd(pt.M, pt.L);
    if (!pt.L_ok) pt.err = "not positive definite at this eta";
    // param_solve() and param_logdet() in closed form; the product by the
    // identity becomes + 0.0 (REVIEW_penalty 3.1)
    pt.si = Dense(p, p);
    const double e2 = pt.eta[1];
    if (!ar) {
        pt.pa = R::plogis(e2, 0.0, 1.0, 1, 0);
        pt.pb = R::plogis(-e2, 0.0, 1.0, 1, 0);
        const double one_plus = p * pt.pa;
        const double one_minus = ((double) p / (double) (p - 1)) * pt.pb;
        const double rr = r0 / one_plus, den = v0 * one_minus;
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < p; ++i)
                pt.si(i, j) = (((i == j ? 1.0 : 0.0) - rr * 1.0) / den) + 0.0;
        pt.ld = (((p * std::log(v0)) + std::log((double) p)) +
                 R::plogis(e2, 0.0, 1.0, 1, 1)) +
                ((double) (p - 1) * ((std::log((double) p) - std::log((double) (p - 1))) +
                                     R::plogis(-e2, 0.0, 1.0, 1, 1)));
    } else {
        const double z = std::fabs(e2);
        const double az = std::exp((-2) * z);
        pt.om = (4 * az) / ((1 + az) * (1 + az));
        const double den = v0 * pt.om;
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < p; ++i) {
                double t = 0.0;
                if (i == j) t = (i == 0 || i == p - 1) ? 1.0 : 1 + r0 * r0;
                else if (std::abs(i - j) == 1) t = -r0;
                pt.si(i, j) = (t / den) + 0.0;
            }
        const double lsech2 = ((2 * std::log(2.0)) - (2 * z)) -
                              (2 * std::log1p(std::exp((-2) * z)));
        pt.ld = (p * std::log(v0)) + ((double) (p - 1) * lsech2);
    }
    return true;
}

inline void econ_d1(const ChartStatic& st, ChartPt& pt) {
    const int p = st.p;
    pt.a.resize(2);
    pt.a[0] = scal(pt.v1, pt.P0);
    pt.a[1] = scal(pt.v0, pt.P1);
    const double f1 = (1.0 * 1.0) / pt.v0;
    pt.rdld.assign(2, 0.0);
    pt.rdld[0] = p * (f1 * pt.v1);
    if (st.code == CH_CS) pt.rdld[1] = pt.pb - (double) (p - 1) * pt.pa;
    else pt.rdld[1] = (double) (p - 1) * ((-2) * std::tanh(pt.eta[1]));
    if (!st.inverted) {
        pt.rai.resize(2);
        for (int k = 0; k < 2; ++k) pt.rai[k] = sandwich1(pt.si, pt.a[k]);
    }
}

inline void econ_d2(const ChartStatic& st, ChartPt& pt) {
    const int p = st.p;
    pt.a2.resize(3);
    pt.a2[pair_index(2, 0, 0)] = scal(pt.v2, pt.P0);
    pt.a2[pair_index(2, 1, 1)] = scal(pt.v0, pt.P2);
    pt.a2[pair_index(2, 0, 1)] = scal(pt.v1, pt.P1);
    const double f1 = (1.0 * 1.0) / pt.v0;
    const double f2 = ((-1.0) * 1.0) / (pt.v0 * pt.v0);
    pt.rd2ld.assign(3, 0.0);
    pt.rd2ld[pair_index(2, 0, 0)] = p * ((f2 * (pt.v1 * pt.v1)) + (f1 * pt.v2));
    if (st.code == CH_CS) pt.rd2ld[pair_index(2, 1, 1)] = ((double) (-p) * pt.pa) * pt.pb;
    else pt.rd2ld[pair_index(2, 1, 1)] = (double) (p - 1) * ((-2) * pt.om);
    pt.rd2ld[pair_index(2, 0, 1)] = 0.0;
    if (!st.inverted) {
        pt.rai2.resize(3);
        for (int k = 0; k < 2; ++k)
            for (int l = k; l < 2; ++l) {
                const int pi = pair_index(2, k, l);
                pt.rai2[pi] = sandwich2(pt.si, pt.a[k], pt.a[l], pt.a2[pi]);
            }
    }
}

// =============================================== dr_prod / correlation
// corr_dfactor(): the derivative of the correlation's factor in the angles
// ks (0, 1 or 2 of them); false (R's NULL) where they lie in two rows
inline bool drp_D(const ChartStatic& st, const ChartPt& pt, const int* ks, int nks, Dense& m) {
    const int p = st.p;
    for (int i = 1; i < nks; ++i) if (st.arow[ks[i]] != st.arow[ks[0]]) return false;
    m = Dense(p, p);
    m(0, 0) = nks ? 0.0 : 1.0;
    std::vector<int> mult(p, 0);
    for (int i = 0; i < nks; ++i) mult[st.acol[ks[i]]]++;
    int r_lo = 1, r_hi = p - 1;
    if (nks) r_lo = r_hi = st.arow[ks[0]];
    for (int i = r_lo; i <= r_hi; ++i) {
        if (i < 1) continue;
        for (int j = 0; j <= i; ++j) {
            const int nused = (j < i) ? j + 1 : i;
            bool bad = false;
            for (int t = 0; t < nks; ++t)
                if (st.arow[ks[t]] == i && st.acol[ks[t]] >= nused) bad = true;
            if (bad) { m(i, j) = 0.0; continue; }
            double v = 1.0;
            for (int c = 0; c < nused; ++c) {
                const int k = st.aidx[(std::size_t) c * p + i];
                const bool use_cs = (j < i && c == j);
                const int mu = mult[c];
                double tab;
                if (use_cs) tab = (mu == 0) ? pt.c0[k] : (mu == 1) ? pt.cs1[k] : pt.cs2[k];
                else tab = (mu == 0) ? pt.s0[k] : (mu == 1) ? pt.sn1[k] : pt.sn2[k];
                v = v * tab;
            }
            m(i, j) = v;
        }
    }
    return true;
}

// leibniz_gram() of the correlation over a sorted tuple of angles (1 or 2),
// diagonal set to 0
inline Dense drp_corr(const ChartStatic& st, const ChartPt& pt, const std::vector<int>& t) {
    const int p = st.p, n = (int) t.size();
    if (!n) return pt.R0;
    Dense acc(p, p);
    for (int b = 0; b < (1 << n); ++b) {
        std::vector<int> tk, tn;
        for (int i = 0; i < n; ++i) ((b >> i) & 1 ? tk : tn).push_back(t[i]);
        Dense A, B;
        if (!drp_D(st, pt, tk.data(), (int) tk.size(), A)) continue;
        if (!drp_D(st, pt, tn.data(), (int) tn.size(), B)) continue;
        acc = add(acc, mm(A, transpose(B)));
    }
    for (int i = 0; i < p; ++i) acc(i, i) = 0.0;
    return acc;
}

// dr_scale_factor() of the sorted scale indices u (0, 1 or 2 of them)
inline Dense drp_sfac(const ChartStatic& st, const ChartPt& pt, const std::vector<int>& u) {
    const int p = st.p;
    const std::vector<double>& d = pt.sd0;
    Dense out(p, p);
    auto sd = [&](int m, int k) { return m == 0 ? pt.sd0[k] : m == 1 ? pt.sd1[k] : pt.sd2[k]; };
    if (u.empty()) {
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < p; ++i) out(i, j) = d[i] * d[j];
        return out;
    }
    if (u.size() == 1 || u[0] == u[1]) {
        const int k = u[0], m = (int) u.size();
        const double f = sd(m, k);
        for (int j = 0; j < p; ++j) out(k, j) = f * d[j];
        for (int i = 0; i < p; ++i) out(i, k) = f * d[i];
        std::vector<double> terms;
        for (int i = 0; i <= m; ++i) terms.push_back((R::choose((double) m, (double) i) * sd(i, k)) * sd(m - i, k));
        out(k, k) = r_sum(terms.data(), terms.size());
        return out;
    }
    const int i = u[0], j = u[1];
    out(i, j) = pt.sd1[i] * pt.sd1[j];
    out(j, i) = pt.sd1[i] * pt.sd1[j];
    return out;
}

// the derivative of Sigma over a sorted tuple of free values
inline Dense drp_amat(const ChartStatic& st, const ChartPt& pt, const std::vector<int>& tup) {
    std::vector<int> sc, an;
    for (int k : tup) { if (k < st.p) sc.push_back(k); else an.push_back(k - st.p); }
    return hadamard(drp_sfac(st, pt, sc), drp_corr(st, pt, an));
}

inline bool drp_d0(const ChartStatic& st, ChartPt& pt) {
    const int p = st.p, na = st.K - p;
    pt.sd0.resize(p); pt.sd1.resize(p); pt.sd2.resize(p);
    // the log link on the VECTOR of the scales: R's wrapper of exp
    // (REVIEW_penalty2 finding 2), floored; linkinv() clamped
    for (int i = 0; i < p; ++i) {
        double v = w_exp(pt.eta[i]);
        if (!ISNAN(v) && v < exp_floor()) v = exp_floor();
        pt.sd0[i] = lf7().clamp(v, st.ls.lwr, st.ls.upr);
        pt.sd1[i] = v;
        pt.sd2[i] = v;
    }
    pt.s0.resize(na); pt.c0.resize(na); pt.sn1.resize(na); pt.sn2.resize(na);
    pt.cs1.resize(na); pt.cs2.resize(na);
    for (int k = 0; k < na; ++k) {
        double t0, t1, t2;
        link_inv3(st.lr, pt.eta[p + k], t0, t1, t2);
        const double s0 = std::sin(t0), c0 = std::cos(t0);
        pt.s0[k] = s0; pt.c0[k] = c0;
        pt.sn1[k] = c0 * t1;
        pt.sn2[k] = ((-s0) * (t1 * t1)) + (c0 * t2);
        pt.cs1[k] = (-s0) * t1;
        pt.cs2[k] = ((-c0) * (t1 * t1)) + ((-s0) * t2);
    }
    drp_D(st, pt, nullptr, 0, pt.L0);
    pt.R0 = tcross(pt.L0);
    for (int i = 0; i < p; ++i) pt.R0(i, i) = 1.0;
    const std::vector<double>& d = pt.sd0;
    pt.M = Dense(p, p);
    pt.L = Dense(p, p);
    for (int j = 0; j < p; ++j)
        for (int i = 0; i < p; ++i) {
            pt.M(i, j) = (d[i] * d[j]) * pt.R0(i, j);
            pt.L(i, j) = d[i] * pt.L0(i, j);
        }
    pt.L_ok = true;
    for (int i = 0; i < p; ++i) if (pt.L(i, i) == 0.0) pt.L_ok = false;
    // si = backsolve(t(L0), forwardsolve(L0, I / d)) / d
    Dense X(p, p);
    for (int j = 0; j < p; ++j)
        for (int i = 0; i < p; ++i) X(i, j) = (i == j ? 1.0 : 0.0) / d[i];
    if (!trsolve(pt.L0, false, X) || !trsolve(transpose(pt.L0), true, X)) {
        pt.err = "singular matrix in 'backsolve'";
        return false;
    }
    pt.si = Dense(p, p);
    for (int j = 0; j < p; ++j)
        for (int i = 0; i < p; ++i) pt.si(i, j) = X(i, j) / d[i];
    std::vector<double> ld(p), ls(na);
    for (int i = 0; i < p; ++i) ld[i] = std::log(d[i]);
    for (int k = 0; k < na; ++k) ls[k] = std::log(pt.s0[k]);
    pt.ld = (2 * r_sum(ld.data(), p)) + (2 * r_sum(ls.data(), na));
    return true;
}

inline void drp_d1(const ChartStatic& st, ChartPt& pt) {
    const int p = st.p, K = st.K;
    pt.a.resize(K);
    pt.rdld.assign(K, 0.0);
    for (int k = 0; k < K; ++k) {
        pt.a[k] = drp_amat(st, pt, std::vector<int>{k});
        if (k < p) pt.rdld[k] = 2 * (pt.sd1[k] / pt.sd0[k]);
        else {
            const int j = k - p;
            const double g1 = ((2.0 * 1.0) * 1.0) / pt.s0[j];
            pt.rdld[k] = g1 * pt.sn1[j];
        }
    }
    if (!st.inverted) {
        pt.rai.resize(K);
        for (int k = 0; k < K; ++k) pt.rai[k] = sandwich1(pt.si, pt.a[k]);
    }
}

inline void drp_d2(const ChartStatic& st, ChartPt& pt) {
    const int p = st.p, K = st.K, np = n_pairs(K);
    pt.a2.resize(np);
    pt.rd2ld.assign(np, 0.0);
    for (int k = 0; k < K; ++k)
        for (int l = k; l < K; ++l) {
            const int pi = pair_index(K, k, l);
            pt.a2[pi] = drp_amat(st, pt, std::vector<int>{k, l});
            double v = 0.0;
            if (k == l) {
                if (k < p) {
                    const double u1 = pt.sd1[k] / pt.sd0[k], u2 = pt.sd2[k] / pt.sd0[k];
                    v = 2 * (u2 - u1 * u1);
                } else {
                    const int j = k - p;
                    const double g1 = ((2.0 * 1.0) * 1.0) / pt.s0[j];
                    const double g2 = ((2.0 * (-1.0)) * 1.0) / (pt.s0[j] * pt.s0[j]);
                    v = (g2 * (pt.sn1[j] * pt.sn1[j])) + (g1 * pt.sn2[j]);
                }
            }
            pt.rd2ld[pi] = v;
        }
    if (!st.inverted) {
        pt.rai2.resize(np);
        for (int k = 0; k < K; ++k)
            for (int l = k; l < K; ++l) {
                const int pi = pair_index(K, k, l);
                pt.rai2[pi] = sandwich2(pt.si, pt.a[k], pt.a[l], pt.a2[pi]);
            }
    }
}

// ================================================ the three orders
inline bool finite_all(const Dense& m) {
    for (double x : m.v) if (!R_FINITE(x)) return false;
    return true;
}

// order 0: L (param_factor()), si, sinv, logdet; the precision orientation
// reads sinv = M and logdet = -ld (mv_matrix_pieces())
inline bool chart_order0(const ChartStatic& st, const double* eta, ChartPt& pt) {
    pt = ChartPt();
    pt.eta.assign(eta, eta + st.K);
    for (double e : pt.eta)
        if (!R_FINITE(e)) { pt.err = "non-finite free value"; return false; }
    bool ok = false;
    switch (st.code) {
    case CH_LC: ok = lc_d0(st, pt); break;
    case CH_CS: case CH_AR1: ok = econ_d0(st, pt); break;
    case CH_DRP:
        // the scales' wrapper is established for the log link only, the
        // angles' arithmetic for the doubly bounded link only
        if (st.ls.id != 1 || st.lr.id != 13 || st.p < 2) {
            pt.err = "dr_prod links not covered";
            return false;
        }
        ok = drp_d0(st, pt);
        break;
    default: pt.err = "chart not covered"; return false;
    }
    if (!ok) return false;
    const int p = st.p;
    if (st.inverted) {
        pt.sinv = pt.M;
        pt.logdet = -pt.ld;
    } else {
        pt.sinv = Dense(p, p);
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < p; ++i) pt.sinv(i, j) = (pt.si(i, j) + pt.si(j, i)) / 2;
        pt.logdet = pt.ld;
    }
    pt.order = 0;
    pt.ok = finite_all(pt.sinv) && finite_all(pt.si) && R_FINITE(pt.logdet);
    if (!pt.ok && pt.err.empty()) pt.err = "non-finite value";
    return pt.ok && pt.L_ok;
}

// order 1: ai, dld (param_inv_d1() and param_dlogdet(); the precision
// orientation reads param_d1() and -param_dlogdet())
inline bool chart_order1(const ChartStatic& st, ChartPt& pt) {
    if (pt.order < 0 || !pt.ok) return false;
    switch (st.code) {
    case CH_LC: lc_d1(st, pt); break;
    case CH_CS: case CH_AR1: econ_d1(st, pt); break;
    case CH_DRP: drp_d1(st, pt); break;
    default: return false;
    }
    const int K = st.K;
    pt.dld.assign(K, 0.0);
    if (st.inverted) {
        pt.ai = pt.a;
        for (int k = 0; k < K; ++k) pt.dld[k] = -pt.rdld[k];
    } else {
        pt.ai = pt.rai;
        pt.dld = pt.rdld;
    }
    pt.order = 1;
    bool fin = true;
    for (const Dense& m : pt.ai) fin = fin && finite_all(m);
    for (double x : pt.dld) fin = fin && R_FINITE(x);
    if (!fin) pt.err = "non-finite first derivative";
    return fin;
}

// order 2: ai2, d2ld over the pairs k <= l
inline bool chart_order2(const ChartStatic& st, ChartPt& pt) {
    if (pt.order < 1) return false;
    switch (st.code) {
    case CH_LC: lc_d2(st, pt); break;
    case CH_CS: case CH_AR1: econ_d2(st, pt); break;
    case CH_DRP: drp_d2(st, pt); break;
    default: return false;
    }
    const int np = n_pairs(st.K);
    pt.d2ld.assign(np, 0.0);
    if (st.inverted) {
        pt.ai2 = pt.a2;
        for (int i = 0; i < np; ++i) pt.d2ld[i] = -pt.rd2ld[i];
    } else {
        pt.ai2 = pt.rai2;
        pt.d2ld = pt.rd2ld;
    }
    pt.order = 2;
    bool fin = true;
    for (const Dense& m : pt.ai2) fin = fin && finite_all(m);
    for (double x : pt.d2ld) fin = fin && R_FINITE(x);
    if (!fin) pt.err = "non-finite second derivative";
    return fin;
}

}  // namespace chart7

#endif
