// The compiled inner fit (piano_motore.md, phase 2; specification in
// misure_motore_fase2/SPEC_iwls.md). engine_iwls_fit() runs iwls_fit() as
// fit_smooth() calls it on a model of the covered set -- a univariate family
// of the distributions7 registry without a wrapper, links of the
// linkfunctions7 registry, design blocks that do not move, a penalty that is
// quadratic in the coefficients at fixed hyperparameters -- step for step:
// the same information (expected, observed, the switch from the first to the
// second and the backup of the second by the first), the same line search,
// stopping rule, polishing step, damping, fallback and verdict, the same
// square-root and assembled routes with the same factorizations.
//
// The penalty arrives from R as its units (engine_penalty_units()), each
// with what its class's penalty_value() and penalty_gradient() read at the
// call's hyperparameters, and the value and the gradient repeat those methods
// operation by operation (D17); its Hessian S and the square root C =
// penalty_sqrt(S), which iwls_pieces() recomputes at every iteration from the
// same S, arrive as matrices.
//
// Where the R route would signal an error (a parameter outside its domain at
// a point the R route evaluates without a guard, a factorization that fails
// outside a tryCatch) the engine stops and returns status 1, and the caller
// runs iwls_fit() from the same start, which reproduces the error.
//
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <R_ext/Rdynload.h>
#include "engine_la.h"
#include "qr_factor_core.h"
#include <d7_link_scale.h>
#include <string>
#include <chrono>
#include <cfloat>

using namespace sm7e;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Trip;

namespace {

// ---------------------------------------------------------------- registry
typedef int (*d7_id_t)(const char*);
typedef double (*d7_logpdf_t)(int, double, const double*);
typedef void (*d7_parts_t)(int, double, const double*, int, double*);
typedef void (*d7_parts_exp_t)(int, const double*, int, double*);
typedef void (*d7_parts_out_t)(int, const double*, double*);
typedef int (*d7_cov_t)(int);
typedef int (*d7_max_t)(void);
typedef int (*lf7_id_t)(const char*);
typedef void (*lf7_inv12p_t)(int, const double*, double, double*, double*,
                             double*);
typedef double (*lf7_clamp_t)(double, double, double);
typedef double (*lf7_fun1p_t)(int, const double*, double);
typedef int (*d7_expect_t)(int, int, int, int, const double*, const double*,
                           double, const double*, double*);

struct Registry {
    bool ok = false;
    d7_id_t id = nullptr;
    d7_logpdf_t logpdf = nullptr;
    d7_parts_t parts = nullptr;
    d7_parts_exp_t parts_expected = nullptr;
    d7_parts_out_t grad = nullptr, hess = nullptr, info = nullptr;
    d7_cov_t covered = nullptr;
    d7_max_t pmax = nullptr;
    lf7_id_t lid = nullptr;
    lf7_inv12p_t inv12p = nullptr;
    lf7_clamp_t clamp = nullptr;
    lf7_fun1p_t fun = nullptr;
    d7_expect_t expect = nullptr;
};

Registry& registry() {
    static Registry r;
    if (!r.ok) {
        r.id = (d7_id_t) R_GetCCallable("distributions7", "d7_scalar_id");
        r.logpdf = (d7_logpdf_t) R_GetCCallable("distributions7", "d7_logpdf");
        r.parts = (d7_parts_t) R_GetCCallable("distributions7", "d7_parts");
        r.parts_expected = (d7_parts_exp_t)
            R_GetCCallable("distributions7", "d7_parts_expected");
        r.grad = (d7_parts_out_t) R_GetCCallable("distributions7", "d7_parts_grad");
        r.hess = (d7_parts_out_t) R_GetCCallable("distributions7", "d7_parts_hess");
        r.info = (d7_parts_out_t) R_GetCCallable("distributions7", "d7_parts_info");
        r.covered = (d7_cov_t) R_GetCCallable("distributions7", "d7_parts_covered");
        r.pmax = (d7_max_t) R_GetCCallable("distributions7", "d7_parts_max");
        r.lid = (lf7_id_t) R_GetCCallable("linkfunctions7", "lf7_class_id");
        r.inv12p = (lf7_inv12p_t) R_GetCCallable("linkfunctions7", "lf7_inv12p");
        r.clamp = (lf7_clamp_t) R_GetCCallable("linkfunctions7", "lf7_clamp");
        r.fun = (lf7_fun1p_t) R_GetCCallable("linkfunctions7", "lf7_fun1p");
        r.expect = (d7_expect_t) R_GetCCallable("distributions7", "d7_expect");
        r.ok = true;
    }
    return r;
}

// ------------------------------------------------------------------ model
struct Block {
    int m = 0;
    bool sp = false;
    const double* d = nullptr;   // n x m, column-major
    SpMat s;                     // n x m
};

// a matrix as R holds it: a base matrix or a dgCMatrix
struct MatAny {
    bool sp = false;
    Dense d;
    SpMat s;
    int nr = 0, nc = 0;
};

// M %*% x: dgemv on a base matrix, the column-wise sparse product otherwise
inline void mat_vec(const MatAny& M, const double* x, double* y) {
    if (!M.sp) { r_matvec(M.d.v.data(), M.nr, M.nc, x, y); return; }
    for (int i = 0; i < M.nr; ++i) y[i] = 0.0;
    for (int j = 0; j < M.nc; ++j) {
        const double t = x[j];
        for (SpMat::InnerIterator it(M.s, j); it; ++it)
            y[it.row()] = y[it.row()] + it.value() * t;
    }
}

// crossprod(M, x)
inline void mat_tvec(const MatAny& M, const double* x, double* y) {
    if (!M.sp) { r_crossprod_vec(M.d.v.data(), M.nr, M.nc, x, y); return; }
    for (int j = 0; j < M.nc; ++j) {
        double t = 0.0;
        for (SpMat::InnerIterator it(M.s, j); it; ++it) t = t + it.value() * x[it.row()];
        y[j] = t;
    }
}

// one penalized unit, described by engine_penalty_units(): kind 1 a
// QuadraticPenalty, 2 an AdditivePenalty, 3 a DistribPenalty with a gaussian
// prior of fixed mean, 4 one with a multivariate gaussian prior
struct PenUnit {
    int kind = 0;
    std::vector<int> idx;            // 0-based stacked positions
    bool has_map = false;
    MatAny map;
    MatAny P, Mg, S;                 // kind 1: P, penalty_matrix; kind 2: S
    std::vector<MatAny> mats;        // kind 2
    std::vector<double> lams;
    double lam = 0, a = 0, c1 = 0, c2 = 0, c3 = 0;
    double mu = 0, sigma = 1, lj = NA_REAL;  // kind 3 (lj also kind 4)
    int block = 1;                   // kind 4
    std::vector<double> muv;
    Dense L, Sinv;
    bool inverted = false;
    double logdet = 0, p2pi = 0;
};

struct Model {
    int n = 0, K = 0, p = 0;
    std::vector<int> npar, offs;
    std::vector<Block> X;
    bool sparse = false;
    const double* y = nullptr;
    const double* w = nullptr;
    std::vector<const double*> off;
    int fam = -1, nc = 0;
    const double* cons = nullptr;    // n x nc
    std::vector<double> lo, hi;
    std::vector<int> lid;
    std::vector<std::vector<double>> lpar;
    std::vector<double> llo, lhi;    // the links' bounds (link_bounds)
    int pmax = 0;
    // penalty at the call's hyperparameters
    Dense Sd;
    SpMat Ss;
    std::vector<PenUnit> units;
    int gauss_id = -1;
    bool hasC = false, Csp = false;
    Dense Cd;
    SpMat Cs;
    std::vector<double> Ccol2;       // colSums(C^2)
    int threads = 1;
};

// the method, settled in R (iwls_resolve()) and read off once
struct Method {
    std::string how;                 // "qr", "svd", "chol"
    bool has_stol = false;
    double stol = 0.0, floor = 1e-8;
    int maxit = 1000, step_halving = 30;
    double tol = 1e-6, dec_limit = 1e-6, switch_after = 10.0,
           switch_rate = 0.5, tol_def = 1e-8;
    bool kind_expected = false;      // the information of pieces_at
    bool backup = false;             // config D: expected stands in
    bool switching = false;          // config C: switch_at exists
    bool damp_on_reject = true, guard = false;
    bool polish = true;              // iwls_polish()
    int emethod = 0;                 // expected information: 0 closed form, 1 OPG,
                                     // 2 simulation, 3 integration
    bool escore = false;             // 2, 3: the score identity
    double nsim = 0;                 // 2: the draws per setting
    double etol[6] = {0, 0, 0, 0, 0, 0};  // 3: quad_vec and series_vec's limits
    double support[2] = {0, 0};      // the family's support
};

// R's error, which the caller reproduces by running the R route
struct RError {};

// seconds spent per phase, returned with the result
enum { T_POINT, T_PARTS, T_PIECES, T_SPQR, T_DQR, T_ASM, T_FN, T_N };
double g_time[T_N];
struct Timer {
    int k;
    std::chrono::steady_clock::time_point t0;
    explicit Timer(int k_) : k(k_), t0(std::chrono::steady_clock::now()) {}
    ~Timer() {
        g_time[k] += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    }
};

// ------------------------------------------------------------- evaluation
struct Point {
    std::vector<double> beta;
    std::vector<double> eta, th, h1, h2;   // n*K, parameter k at k*n
    bool inb = true;
    bool parts = false, eparts = false;
    std::vector<double> buf, ebuf;         // n * pmax
};

inline bool same_beta(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!(a[i] == b[i]) && !(ISNAN(a[i]) && ISNAN(b[i]))) return false;
    return true;
}

class Engine {
public:
    Model M;
    Method me;
    Registry& reg;
    Point cur;          // the cache for the last point whose pieces were asked
    std::vector<double> oob;   // the first point outside the domain: k, i, eta, theta
    std::vector<double> thbuf;

    Engine() : reg(registry()) {}

    // statmod_eta(): eta_k = X_k b_k (+ offset), theta = h(eta); h' and h''
    // with them
    void point_at(const std::vector<double>& beta, Point& P) {
        Timer tm(T_POINT);
        const int n = M.n, K = M.K;
        P.beta = beta;
        P.eta.assign((std::size_t) n * K, 0.0);
        P.th.resize((std::size_t) n * K);
        P.h1.resize((std::size_t) n * K);
        P.h2.resize((std::size_t) n * K);
        P.parts = P.eparts = false;
        for (int k = 0; k < K; ++k) {
            double* e = P.eta.data() + (std::size_t) k * n;
            const Block& B = M.X[k];
            const double* b = beta.data() + M.offs[k];
            if (B.m > 0) {
                bool allfin = true, allnz = true, anynz = false;
                for (int j = 0; j < B.m; ++j) {
                    if (!R_FINITE(b[j])) allfin = false;
                    if (b[j] != 0.0) anynz = true; else allnz = false;
                }
                if (!B.sp) {
                    if (!allfin || allnz) {
                        r_matvec(B.d, n, B.m, b, e);
                    } else if (anynz) {
                        // X[, nz] %*% b[nz]
                        for (int i = 0; i < n; ++i) e[i] = 0.0;
                        for (int j = 0; j < B.m; ++j) {
                            if (b[j] == 0.0) continue;
                            const double t = b[j];
                            const double* xj = B.d + (std::size_t) j * n;
                            for (int i = 0; i < n; ++i) e[i] = e[i] + t * xj[i];
                        }
                    }
                } else if (anynz || !allfin) {
                    for (int j = 0; j < B.m; ++j) {
                        if (allfin && b[j] == 0.0) continue;
                        const double t = b[j];
                        for (SpMat::InnerIterator it(B.s, j); it; ++it)
                            e[it.row()] = e[it.row()] + it.value() * t;
                    }
                }
            }
            if (M.off[k]) {
                const double* o = M.off[k];
                for (int i = 0; i < n; ++i) e[i] = e[i] + o[i];
            }
            double* t = P.th.data() + (std::size_t) k * n;
            double* a = P.h1.data() + (std::size_t) k * n;
            double* c = P.h2.data() + (std::size_t) k * n;
            const double* par = M.lpar[k].empty() ? nullptr : M.lpar[k].data();
            // linkinv(): the class's inverse, then the generic's clamp to
            // the link's bounds (link_bounds_clamp(), lf7_clamp()); the
            // derivatives of the inverse at linkfun() of that value, as
            // inverse_link_derivs() takes them
            for (int i = 0; i < n; ++i) {
                double h0;
                reg.inv12p(M.lid[k], par, e[i], &h0, &a[i], &c[i]);
                t[i] = reg.clamp(h0, M.llo[k], M.lhi[k]);
                const double e2 = reg.fun(M.lid[k], par, t[i]);
                reg.inv12p(M.lid[k], par, e2, &h0, &a[i], &c[i]);
            }
        }
        P.inb = true;
        for (int k = 0; k < K && P.inb; ++k) {
            const double* t = P.th.data() + (std::size_t) k * n;
            for (int i = 0; i < n; ++i) {
                const double v = t[i];
                if (!(R_FINITE(v) && v > M.lo[k] && v < M.hi[k])) {
                    P.inb = false;
                    if (oob.empty()) oob = {(double) k, (double) i, P.eta[(std::size_t) k * n + i], v};
                    break;
                }
            }
        }
    }

    const double* theta_obs(const Point& P, int i) {
        const int n = M.n, K = M.K;
        for (int k = 0; k < K; ++k) thbuf[k] = P.th[(std::size_t) k * n + i];
        for (int c = 0; c < M.nc; ++c) thbuf[K + c] = M.cons[(std::size_t) c * n + i];
        return thbuf.data();
    }

    // t = map_apply(penalty, b) for one unit
    void unit_t(const PenUnit& u, const std::vector<double>& beta, std::vector<double>& bu,
                std::vector<double>& t) {
        bu.resize(u.idx.size());
        for (std::size_t j = 0; j < u.idx.size(); ++j) bu[j] = beta[u.idx[j]];
        if (!u.has_map) { t = bu; return; }
        t.assign(u.map.nr, 0.0);
        mat_vec(u.map, bu.data(), t.data());
    }

    // penalty_value() of one unit, as its class's method computes it
    double unit_value(const PenUnit& u, const std::vector<double>& beta) {
        std::vector<double> bu, t;
        unit_t(u, beta, bu, t);
        const std::size_t m = t.size();
        if (u.kind == 1) {
            std::vector<double> pt(u.P.nr), tp(m);
            mat_vec(u.P, t.data(), pt.data());
            for (std::size_t i = 0; i < m; ++i) tp[i] = t[i] * pt[i];
            const double q = r_sum(tp.data(), m);
            return ((u.a * q - u.c1) + u.c2) - u.c3;
        }
        if (u.kind == 2) {
            std::vector<double> st(u.S.nr), tp(m);
            mat_vec(u.S, t.data(), st.data());
            for (std::size_t i = 0; i < m; ++i) tp[i] = t[i] * st[i];
            const double q = r_sum(tp.data(), m);
            return (0.5 * q - u.c1) + u.c2;
        }
        double v;
        if (u.kind == 3) {
            std::vector<double> lp(m);
            const double th[2] = {u.mu, u.sigma};
            for (std::size_t i = 0; i < m; ++i) lp[i] = reg.logpdf(M.gauss_id, t[i], th);
            v = -r_sum(lp.data(), m);
        } else {
            const int p = u.block, G = (int) (m / p);
            std::vector<double> ld(G), r(p), w(p), sq(p);
            for (int g = 0; g < G; ++g) {
                for (int j = 0; j < p; ++j) r[j] = t[(std::size_t) g * p + j] - u.muv[j];
                if (u.inverted) {
                    r_crossprod_vec(u.L.v.data(), p, p, r.data(), w.data());
                } else {
                    w = r;
                    if (!r_trisolve(u.L, false, w)) throw RError();
                }
                for (int j = 0; j < p; ++j) sq[j] = w[j] * w[j];
                const double q = r_colsum(sq.data(), p);
                ld[g] = -0.5 * ((u.p2pi + u.logdet) + q);
            }
            v = -r_sum(ld.data(), G);
        }
        if (!u.has_map || ISNAN(u.lj)) return v;
        return v - u.lj;
    }

    // penalty_gradient() of one unit, added into pg at its positions
    void unit_gradient(const PenUnit& u, const std::vector<double>& beta,
                       std::vector<double>& pg) {
        std::vector<double> bu, t;
        unit_t(u, beta, bu, t);
        const std::size_t k = u.idx.size(), m = t.size();
        std::vector<double> g(k);
        if (u.kind == 1) {
            std::vector<double> y(u.Mg.nr);
            mat_vec(u.Mg, bu.data(), y.data());
            for (std::size_t j = 0; j < k; ++j) g[j] = u.lam * y[j];
        } else if (u.kind == 2) {
            std::vector<double> y(k);
            for (std::size_t h = 0; h < u.mats.size(); ++h) {
                mat_vec(u.mats[h], bu.data(), y.data());
                for (std::size_t j = 0; j < k; ++j) {
                    const double v = u.lams[h] * y[j];
                    g[j] = h == 0 ? v : g[j] + v;
                }
            }
        } else {
            // distrib_grad_y() + 0 * a, flattened, then -map_back()
            std::vector<double> gy(m);
            if (u.kind == 3) {
                const double s2 = u.sigma * u.sigma;
                for (std::size_t i = 0; i < m; ++i)
                    gy[i] = (-(t[i] - u.mu)) / s2 + 0.0 * t[i];
            } else {
                const int p = u.block, G = (int) (m / p);
                std::vector<double> r(p);
                for (int gg = 0; gg < G; ++gg) {
                    for (int j = 0; j < p; ++j) r[j] = t[(std::size_t) gg * p + j] - u.muv[j];
                    for (int j = 0; j < p; ++j) {
                        double w = 0.0;
                        for (int l = 0; l < p; ++l) w = w + u.Sinv(l, j) * r[l];
                        const std::size_t e = (std::size_t) gg * p + j;
                        gy[e] = -w + 0.0 * t[e];
                    }
                }
            }
            if (!u.has_map) {
                for (std::size_t j = 0; j < k; ++j) g[j] = -gy[j];
            } else {
                std::vector<double> bb(k);
                mat_tvec(u.map, gy.data(), bb.data());
                for (std::size_t j = 0; j < k; ++j) g[j] = -bb[j];
            }
        }
        for (std::size_t j = 0; j < k; ++j) pg[u.idx[j]] = pg[u.idx[j]] + g[j];
    }

    // statmod_penalty_at(): the value, and the gradient
    double pen_value(const std::vector<double>& b) {
        double value = 0.0;
        for (const PenUnit& u : M.units) value = value + unit_value(u, b);
        return value;
    }
    void pen_sb(const std::vector<double>& b, std::vector<double>& pg) {
        pg.assign(M.p, 0.0);
        for (const PenUnit& u : M.units) unit_gradient(u, b, pg);
    }

    // obj$fn(): -sum(w * log f) + penalty; false where R signals an error
    bool fn(const std::vector<double>& beta, double& value) {
        Point P;
        point_at(beta, P);
        if (!P.inb) return false;
        value = value_at(P);
        return true;
    }
    // D22: two families whose R density tests a condition over the whole
    // vector and switches every observation's formula on it, where the
    // registry decides observation by observation. Where R's switch would
    // differ from the registry the fit is handed to iwls_fit(): weibull3
    // (registry id 21) when any scale exp(log(mean) - lgamma(1 + 1/sigma))
    // is not finite and positive; betabinom2 (id 42) when the criterion of
    // betabinom_log_mass() holds for some observations in the support and
    // not for others.
    void vector_branch_check(Point& P) {
        const int n = M.n;
        if (M.fam == 21) {
            for (int i = 0; i < n; ++i) {
                const double* t = theta_obs(P, i);
                const double b = std::exp(std::log(t[0]) - R::lgammafn(1 + 1 / t[1]));
                if (!(R_FINITE(b) && b > 0)) throw RError();
            }
        } else if (M.fam == 42) {
            int seen = -1;
            for (int i = 0; i < n; ++i) {
                const double* t = theta_obs(P, i);
                const double y = M.y[i], size = t[2];
                if (!(y >= 0 && y <= size && y == std::nearbyint(y))) continue;
                const double s = t[0] + t[1];
                const int c = R_FINITE(s) && R::lgammafn(s + size) * DBL_EPSILON < 1e-8;
                if (seen < 0) seen = c;
                else if (c != seen) throw RError();
            }
        }
    }
    double value_at(Point& P) {
        Timer tm(T_FN);
        const int n = M.n;
        if (M.fam == 21 || M.fam == 42) vector_branch_check(P);
        std::vector<double> v(n);
        for (int i = 0; i < n; ++i)
            v[i] = M.w[i] * reg.logpdf(M.fam, M.y[i], theta_obs(P, i));
        const double ll = r_sum(v.data(), n);
        return -ll + pen_value(P.beta);
    }

    void ensure_parts(Point& P) {
        if (P.parts) return;
        Timer tm(T_PARTS);
        const int n = M.n;
        P.buf.resize((std::size_t) n * M.pmax);
        for (int i = 0; i < n; ++i)
            reg.parts(M.fam, M.y[i], theta_obs(P, i), 2,
                      P.buf.data() + (std::size_t) i * M.pmax);
        P.parts = true;
    }
    void ensure_eparts(Point& P) {
        if (P.eparts) return;
        const int n = M.n;
        P.ebuf.resize((std::size_t) n * M.pmax);
        for (int i = 0; i < n; ++i)
            reg.parts_expected(M.fam, theta_obs(P, i), 2,
                               P.ebuf.data() + (std::size_t) i * M.pmax);
        P.eparts = true;
    }

    Point& point(const std::vector<double>& beta) {
        if (!same_beta(cur.beta, beta) || cur.th.empty()) point_at(beta, cur);
        return cur;
    }

    // obj$gr(): -X'(w * dl/deta) + S b; false where R signals an error
    bool gr(const std::vector<double>& beta, std::vector<double>& g) {
        Point& P = point(beta);
        if (!P.inb) return false;
        ensure_parts(P);
        const int n = M.n, K = M.K;
        std::vector<double> gth(K), gl(K), gv((std::size_t) n * K);
        for (int i = 0; i < n; ++i) {
            reg.grad(M.fam, P.buf.data() + (std::size_t) i * M.pmax, gth.data());
            for (int k = 0; k < K; ++k) {
                const double h1 = P.h1[(std::size_t) k * n + i];
                d7ls::link_grad(1, &gth[k], &h1, &gl[k]);
                gv[(std::size_t) k * n + i] = M.w[i] * gl[k];
            }
        }
        std::vector<double> sb;
        pen_sb(beta, sb);
        g.assign(M.p, 0.0);
        for (int k = 0; k < K; ++k) {
            const Block& B = M.X[k];
            if (B.m == 0) continue;
            std::vector<double> s(B.m);
            const double* v = gv.data() + (std::size_t) k * n;
            if (!B.sp) {
                r_crossprod_vec(B.d, n, B.m, v, s.data());
            } else {
                for (int j = 0; j < B.m; ++j) {
                    double t = 0.0;
                    for (SpMat::InnerIterator it(B.s, j); it; ++it)
                        t = t + it.value() * v[it.row()];
                    s[j] = t;
                }
            }
            for (int j = 0; j < B.m; ++j) {
                const int c = M.offs[k] + j;
                g[c] = -s[j] + sb[c];
            }
        }
        return true;
    }

    // distrib_expected_hessian(method = expect_mc() / expect_integrate()) at
    // every observation, on the parameter scale, through d7_expect(); called
    // at every evaluation, as the R route calls it, since a simulation draws
    // anew each time
    std::vector<double> Eall;
    void expect_all(Point& P) {
        const int n = M.n, K = M.K, ncol = K + M.nc;
        std::vector<double> th((std::size_t) n * ncol);
        for (int i = 0; i < n; ++i) {
            const double* t = theta_obs(P, i);
            for (int c = 0; c < ncol; ++c) th[i + (std::size_t) n * c] = t[c];
        }
        Eall.assign((std::size_t) n * K * K, 0.0);
        const int method = me.emethod == 2 ? (me.escore ? 1 : 0) : (me.escore ? 3 : 2);
        const int st = reg.expect(M.fam, method, n, ncol, th.data(), me.support, me.nsim,
                                  me.etol, Eall.data());
        if (st != 0) throw RError();
    }

    // the K x K link-scale Hessian of the log density at observation i
    void hess_obs(Point& P, int i, bool expected, double* out) {
        const int n = M.n, K = M.K;
        double h1[16], h2[16], gth[16], H[256];
        for (int k = 0; k < K; ++k) {
            h1[k] = P.h1[(std::size_t) k * n + i];
            h2[k] = P.h2[(std::size_t) k * n + i];
        }
        if (expected && me.emethod == 1) {
            // expected_by_opg(): -g_a g_b at the observation, on the
            // parameter scale, then the expected route's link scale
            reg.grad(M.fam, P.buf.data() + (std::size_t) i * M.pmax, gth);
            for (int a = 0; a < K; ++a)
                for (int b = a; b < K; ++b) {
                    const double v = -gth[a] * gth[b];
                    H[a + b * K] = v;
                    H[b + a * K] = v;
                }
            d7ls::link_info(K, H, h1, h2, out);
        } else if (expected && me.emethod >= 2) {
            // the parameter-scale expected Hessian of expect_all()
            d7ls::link_info(K, Eall.data() + (std::size_t) i * K * K, h1, h2, out);
        } else if (expected) {
            reg.info(M.fam, P.ebuf.data() + (std::size_t) i * M.pmax, H);
            d7ls::link_info(K, H, h1, h2, out);
        } else {
            const double* b = P.buf.data() + (std::size_t) i * M.pmax;
            reg.grad(M.fam, b, gth);
            reg.hess(M.fam, b, H);
            d7ls::link_hess(K, gth, H, h1, h2, out);
        }
    }

    // ---------------------------------------------------------- pieces
    struct Pieces {
        bool sqrt = false;
        Dense Rd;     // nK x p
        SpMat Rs;
        Dense Ad;     // p x p
        SpMat As;
        bool asp = false;
        // R'R + C'C of the sparse square-root route, formed at the first solve
        mutable bool haveK = false;
        mutable SpMat K;
    };

    // iwls_pieces(): false where R signals an error
    bool pieces(const std::vector<double>& beta, bool expected, Pieces& pc) {
        Point& P = point(beta);
        if (!P.inb) return false;
        if (expected && me.emethod == 0) ensure_eparts(P); else ensure_parts(P);
        if (expected && me.emethod >= 2) expect_all(P);
        Timer tm(T_PIECES);
        const int n = M.n, K = M.K, p = M.p;
        // H (link scale), per observation K x K; Om = -w H
        std::vector<double> Hall((std::size_t) n * K * K);
        for (int i = 0; i < n; ++i)
            hess_obs(P, i, expected, Hall.data() + (std::size_t) i * K * K);
        pc = Pieces();
        bool sqrt_route = (me.how == "qr" || me.how == "svd") && M.hasC;
        std::vector<double> L;
        if (sqrt_route) {
            // chol_blocks(info_blocks())
            L.assign((std::size_t) n * K * K, 0.0);
            auto Om = [&](int i, int a, int b) {
                return -M.w[i] * Hall[(std::size_t) i * K * K + a + b * K];
            };
            auto Lr = [&](int i, int a, int b) -> double& {
                return L[(std::size_t) i * K * K + a + b * K];
            };
            for (int j = 0; j < K && sqrt_route; ++j) {
                // the pivot of column j at every observation, then the test
                std::vector<double> s(n);
                for (int i = 0; i < n; ++i) {
                    double v = Om(i, j, j);
                    for (int m = 0; m < j; ++m) v = v - Lr(i, j, m) * Lr(i, j, m);
                    s[i] = v;
                }
                for (int i = 0; i < n; ++i)
                    if (ISNAN(s[i]) || s[i] <= 0) { sqrt_route = false; break; }
                if (!sqrt_route) break;
                for (int i = 0; i < n; ++i) Lr(i, j, j) = std::sqrt(s[i]);
                for (int a = j + 1; a < K; ++a) {
                    for (int i = 0; i < n; ++i) {
                        double v = Om(i, a, j);
                        for (int m = 0; m < j; ++m) v = v - Lr(i, a, m) * Lr(i, j, m);
                        Lr(i, a, j) = v / Lr(i, j, j);
                    }
                }
            }
            if (sqrt_route)
                for (double v : L) if (!R_FINITE(v)) { sqrt_route = false; break; }
        }
        if (sqrt_route) {
            // sqrt_design(): row block a holds L[, b, a] * X_b for b >= a
            pc.sqrt = true;
            const std::size_t nr = (std::size_t) n * K;
            auto Lv = [&](int i, int b, int a) {
                return L[(std::size_t) i * K * K + b + a * K];
            };
            if (!M.sparse) {
                pc.Rd = Dense((int) nr, p);
                for (int a = 0; a < K; ++a) {
                    for (int b = a; b < K; ++b) {
                        const Block& B = M.X[b];
                        for (int c = 0; c < B.m; ++c) {
                            double* col = pc.Rd.col(M.offs[b] + c) + (std::size_t) a * n;
                            const double* xc = B.d + (std::size_t) c * n;
                            for (int i = 0; i < n; ++i) col[i] = Lv(i, b, a) * xc[i];
                        }
                    }
                }
            } else {
                // the pattern of R is that of the design blocks, the same at
                // every iteration: built once, then only the values change
                if (!rp.built) build_rpattern();
                std::vector<double> val(rp.sx.size());
                for (std::size_t e = 0; e < val.size(); ++e)
                    val[e] = Lv(rp.si[e], rp.sb[e], rp.sa[e]) * rp.sx[e];
                Eigen::Map<const SpMat> R((int) nr, p, (int) val.size(),
                                          rp.outer.data(), rp.inner.data(), val.data());
                pc.Rs = R;
                // exact zeros (a vanishing block of L) are not stored, as
                // cbind2() does not store them on the R route
                pc.Rs.prune([](int, int, double v) { return v != 0.0; });
            }
            return true;
        }
        // assembled(H): -sum_{blocks} X_a' (w H_ab) X_b + S
        pc.sqrt = false;
        std::vector<double> wv(n);
        if (!M.sparse) {
            pc.Ad = Dense(p, p);
            for (int a = 0; a < K; ++a) {
                const Block& Xa = M.X[a];
                if (Xa.m == 0) continue;
                for (int b = a; b < K; ++b) {
                    const Block& Xb = M.X[b];
                    if (Xb.m == 0) continue;
                    for (int i = 0; i < n; ++i)
                        wv[i] = M.w[i] * Hall[(std::size_t) i * K * K + a + b * K];
                    std::vector<double> xw(n);
                    for (int j = 0; j < Xa.m; ++j) {
                        const double* xj = Xa.d + (std::size_t) j * n;
                        for (int i = 0; i < n; ++i) xw[i] = xj[i] * wv[i];
                        for (int k = 0; k < Xb.m; ++k) {
                            const double* xk = Xb.d + (std::size_t) k * n;
                            double t = 0.0;
                            for (int i = 0; i < n; ++i) t = t + xw[i] * xk[i];
                            pc.Ad(M.offs[a] + j, M.offs[b] + k) = -t;
                            if (b > a) pc.Ad(M.offs[b] + k, M.offs[a] + j) = -t;
                        }
                    }
                }
            }
            for (int j = 0; j < p; ++j)
                for (int i = 0; i < p; ++i) pc.Ad(i, j) = pc.Ad(i, j) + M.Sd(i, j);
        } else {
            pc.asp = true;
            SpMat A(p, p);
            for (int a = 0; a < K; ++a) {
                const Block& Xa = M.X[a];
                if (Xa.m == 0) continue;
                for (int b = a; b < K; ++b) {
                    const Block& Xb = M.X[b];
                    if (Xb.m == 0) continue;
                    for (int i = 0; i < n; ++i)
                        wv[i] = M.w[i] * Hall[(std::size_t) i * K * K + a + b * K];
                    SpMat sa = Xa.sp ? Xa.s : dense_to_sparse(Xa);
                    SpMat sb = Xb.sp ? Xb.s : dense_to_sparse(Xb);
                    Eigen::Map<const Eigen::VectorXd> wvm(wv.data(), n);
                    SpMat blk = SpMat(sa.transpose() * wvm.asDiagonal() * sb);
                    std::vector<Trip> tr;
                    for (int c = 0; c < blk.outerSize(); ++c)
                        for (SpMat::InnerIterator it(blk, c); it; ++it) {
                            tr.emplace_back(M.offs[a] + it.row(), M.offs[b] + c, -it.value());
                            if (b > a)
                                tr.emplace_back(M.offs[b] + c, M.offs[a] + it.row(), -it.value());
                        }
                    SpMat part(p, p);
                    part.setFromTriplets(tr.begin(), tr.end());
                    A = A + part;
                }
            }
            pc.As = A + M.Ss;
        }
        return true;
    }

    // the structural pattern of the sparse square root R: column offs[b] + c
    // holds rows a * n + i, for a <= b and the nonzero rows i of column c of
    // X_b, in ascending order
    struct RPattern {
        bool built = false;
        std::vector<int> outer, inner, sb, sa, si;
        std::vector<double> sx;
    } rp;
    void build_rpattern() {
        const int n = M.n, K = M.K;
        rp.outer.assign(1, 0);
        for (int b = 0; b < K; ++b) {
            const Block& B = M.X[b];
            for (int c = 0; c < B.m; ++c) {
                for (int a = 0; a <= b; ++a) {
                    if (B.sp) {
                        for (SpMat::InnerIterator it(B.s, c); it; ++it) {
                            rp.inner.push_back(a * n + (int) it.row());
                            rp.sb.push_back(b); rp.sa.push_back(a);
                            rp.si.push_back((int) it.row()); rp.sx.push_back(it.value());
                        }
                    } else {
                        const double* xc = B.d + (std::size_t) c * n;
                        for (int i = 0; i < n; ++i) {
                            if (xc[i] == 0.0) continue;
                            rp.inner.push_back(a * n + i);
                            rp.sb.push_back(b); rp.sa.push_back(a);
                            rp.si.push_back(i); rp.sx.push_back(xc[i]);
                        }
                    }
                }
                rp.outer.push_back((int) rp.inner.size());
            }
        }
        rp.built = true;
    }

    SpMat dense_to_sparse(const Block& B) {
        std::vector<Trip> tr;
        for (int c = 0; c < B.m; ++c) {
            const double* xc = B.d + (std::size_t) c * M.n;
            for (int i = 0; i < M.n; ++i) if (xc[i] != 0.0) tr.emplace_back(i, c, xc[i]);
        }
        SpMat s(M.n, B.m);
        s.setFromTriplets(tr.begin(), tr.end());
        return s;
    }

    Dense dense_of(const SpMat& s) {
        Dense d((int) s.rows(), (int) s.cols());
        for (int c = 0; c < s.outerSize(); ++c)
            for (SpMat::InnerIterator it(s, c); it; ++it) d(it.row(), c) = it.value();
        return d;
    }

    // ------------------------------------------------- scale and checks
    double col2_R(const Pieces& pc, int j) {
        if (!M.sparse) {
            const double* c = pc.Rd.col(j);
            const std::size_t nr = (std::size_t) pc.Rd.nr;
            ldouble s = 0.0;
            for (std::size_t i = 0; i < nr; ++i) s += c[i] * c[i];
            return (double) s;
        }
        ldouble s = 0.0;
        for (SpMat::InnerIterator it(pc.Rs, j); it; ++it) s += it.value() * it.value();
        return (double) s;
    }
    double adiag(const Pieces& pc, int j) {
        return pc.asp ? pc.As.coeff(j, j) : pc.Ad(j, j);
    }

    // iwls_scale()
    double scale(const Pieces& pc) {
        double m = R_NegInf;
        bool any = false;
        for (int j = 0; j < M.p; ++j) {
            double v = pc.sqrt ? col2_R(pc, j) + (M.hasC ? M.Ccol2[j] : 0.0)
                               : adiag(pc, j);
            if (!R_FINITE(v)) continue;
            if (!any || v > m) m = v;
            any = true;
        }
        if (!any) return 1.0;
        return (R_FINITE(m) && m > 0) ? m : 1.0;
    }

    // iwls_info_diag()
    std::vector<double> info_diag(const Pieces& pc) {
        std::vector<double> d(M.p);
        for (int j = 0; j < M.p; ++j) d[j] = pc.sqrt ? col2_R(pc, j) : adiag(pc, j);
        return d;
    }

    // pieces_definite()
    bool definite(const Pieces& pc, double rel) {
        if (pc.sqrt) return true;
        Dense A = pc.asp ? dense_of(pc.As) : pc.Ad;
        for (double v : A.v) if (!R_FINITE(v)) return true;
        std::vector<double> ev;
        Dense V;
        if (!r_eigen_sym(A, true, ev, V)) throw RError();
        double mx = 0.0;
        for (double v : ev) mx = std::max(mx, std::fabs(v));
        const double mn = *std::min_element(ev.begin(), ev.end());
        return mn > rel * mx;
    }

    // --------------------------------------------------------- solves
    struct Sol {
        std::vector<double> delta;
        int rank = 0;
        std::string route;
        std::vector<int> held, dropped;   // 0-based, stacked numbering
    };

    // pd_repair(): false for NULL
    bool pd_repair(Dense& A, double rel) {
        if (A.nr == 0) return true;
        for (double v : A.v) if (!R_FINITE(v)) return false;
        Dense U;
        if (r_chol(A, U)) return true;
        std::vector<double> ev;
        Dense V;
        if (!r_eigen_sym(A, false, ev, V)) throw RError();
        double mx = 0.0;
        for (double v : ev) mx = std::max(mx, std::fabs(v));
        double fl = rel * mx;
        if (!R_FINITE(fl) || fl <= 0) fl = rel;
        const int n = A.nr;
        Dense T(n, n);   // pmax(values, floor) * t(vectors)
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) T(i, j) = std::max(ev[i], fl) * V(j, i);
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) A(i, j) = 0.0;
            for (int l = 0; l < n; ++l) {
                const double t = T(l, j);
                for (int i = 0; i < n; ++i) A(i, j) = A(i, j) + t * V(i, l);
            }
        }
        return true;
    }

    // augmented_solve() on the dense augmented matrix
    void aug_dense(const Dense& A, const std::vector<double>& u,
                   std::vector<double>& delta, int& rank,
                   std::vector<int>& dropped, bool& err) {
        const int k = A.nc;
        delta.assign(k, 0.0);
        dropped.clear();
        for (double v : A.v) if (!R_FINITE(v)) { err = true; return; }
        if (me.how == "svd") {
            std::vector<double> d;
            Dense v;
            if (!r_svd(A, d, v)) { err = true; return; }
            double dm = 0.0;
            for (double x : d) dm = std::max(dm, x);
            const double cut = me.has_stol ? me.stol * dm
                : std::max(A.nr, A.nc) * DBL_EPSILON * dm;
            const int np = (int) d.size();
            std::vector<double> dinv(np, 0.0), t(np);
            rank = 0;
            for (int j = 0; j < np; ++j)
                if (d[j] > cut) { dinv[j] = 1.0 / (d[j] * d[j]); ++rank; }
            r_crossprod_vec(v.v.data(), v.nr, np, u.data(), t.data());
            for (int j = 0; j < np; ++j) t[j] = dinv[j] * t[j];
            r_matvec(v.v.data(), v.nr, np, t.data(), delta.data());
            return;
        }
        const double qtol = me.has_stol ? me.stol : 1e-7;
        if (M.threads > 1 && me.how == "qr" &&
            (double) A.nr * A.nc * A.nc >= 5e7) {
            std::vector<double> cn(k);
            bool okc = true;
            for (int j = 0; j < k; ++j) {
                const double* c = A.col(j);
                ldouble s = 0.0;
                for (int i = 0; i < A.nr; ++i) s += c[i] * c[i];
                cn[j] = std::sqrt((double) s);
                if (!(R_FINITE(cn[j]) && cn[j] > 0)) okc = false;
            }
            if (okc) {
                Dense B = A, Rt(k, k);
                sm7::qr_factor_core(B.v.data(), A.nr, k, M.threads, Rt.v.data());
                std::vector<double> eq(k);
                bool fin = true;
                double emx = R_NegInf, emn = R_PosInf;
                for (int j = 0; j < k; ++j) {
                    eq[j] = std::fabs(Rt(j, j)) / cn[j];
                    if (!R_FINITE(eq[j])) fin = false;
                    emx = std::max(emx, eq[j]);
                    emn = std::min(emn, eq[j]);
                }
                if (fin && emx > 0 && emn > qtol * emx) {
                    std::vector<double> d(u);
                    if (r_chol_solve(Rt, d)) {
                        bool f = true;
                        for (double x : d) if (!R_FINITE(x)) f = false;
                        if (f) { delta = d; rank = k; return; }
                    }
                }
            }
        }
        QRResult q;
        Timer tm(T_DQR);
        if (!r_qr(A, qtol, q)) { err = true; return; }
        rank = q.rank;
        // backsolve() rejects k = 0 ("invalid 'k' argument")
        if (rank == 0) { err = true; return; }
        Dense Rf(rank, rank);
        for (int j = 0; j < rank; ++j)
            for (int i = 0; i <= j; ++i) Rf(i, j) = q.qr(i, j);
        std::vector<double> z(rank);
        for (int j = 0; j < rank; ++j) z[j] = u[q.pivot[j] - 1];
        if (!r_chol_solve(Rf, z)) { err = true; return; }
        std::vector<char> in(k, 0);
        for (int j = 0; j < rank; ++j) { delta[q.pivot[j] - 1] = z[j]; in[q.pivot[j] - 1] = 1; }
        for (int j = 0; j < k; ++j) if (!in[j]) dropped.push_back(j);
    }

    // sparse_augmented_solve() by the Cholesky factor of A'A (D16): the
    // equilibrated diagonal L_jj / ||a_j|| is |R_jj| / ||a_j|| of the QR in
    // the same column order, and the rank test is the R route's with a
    // safety margin, below which the dense pivoted QR decides; false for NULL
    static constexpr double kCholSafety = 1e-5;
    bool chol_sparse(const SpMat& K, const std::vector<int>& keep, double damp,
                     const std::vector<double>& u, std::vector<double>& delta) {
        Timer tm(T_SPQR);
        const int k = (int) keep.size();
        const double qtol = me.has_stol ? me.stol : 1e-7;
        std::vector<int> pos(K.cols(), -1);
        for (int j = 0; j < k; ++j) pos[keep[j]] = j;
        std::vector<Trip> tr;
        tr.reserve(K.nonZeros() + k);
        for (int c = 0; c < K.outerSize(); ++c) {
            if (pos[c] < 0) continue;
            for (SpMat::InnerIterator it(K, c); it; ++it)
                if (pos[it.row()] >= 0) tr.emplace_back(pos[it.row()], pos[c], it.value());
        }
        if (damp > 0) for (int j = 0; j < k; ++j) tr.emplace_back(j, j, damp);
        SpMat Kk(k, k);
        Kk.setFromTriplets(tr.begin(), tr.end());
        Eigen::SimplicialLLT<SpMat, Eigen::Lower, Eigen::AMDOrdering<int>> ch(Kk);
        if (ch.info() != Eigen::Success) return false;
        const auto& Lm = ch.matrixL().nestedExpression();
        Eigen::VectorXi pinv = ch.permutationPinv().indices();
        double emx = 0.0, emn = R_PosInf;
        for (int j = 0; j < k; ++j) {
            const double l = std::fabs(Lm.coeff(j, j));
            const double cn = std::sqrt(Kk.coeff(pinv[j], pinv[j]));
            if (!R_FINITE(l) || !R_FINITE(cn) || cn == 0) return false;
            const double e = l / cn;
            emx = std::max(emx, e);
            emn = std::min(emn, e);
        }
        if (!(emx > 0)) return false;
        if (emn <= std::max(qtol, kCholSafety) * emx) return false;
        Eigen::VectorXd uu = Eigen::Map<const Eigen::VectorXd>(u.data(), k);
        Eigen::VectorXd d = ch.solve(uu);
        for (int j = 0; j < k; ++j) if (!R_FINITE(d[j])) return false;
        delta.assign(d.data(), d.data() + k);
        return true;
    }

    SpMat dense_to_sparse_mat(const Dense& D) {
        std::vector<Trip> tr;
        for (int c = 0; c < D.nc; ++c)
            for (int i = 0; i < D.nr; ++i) if (D(i, c) != 0.0) tr.emplace_back(i, c, D(i, c));
        SpMat s(D.nr, D.nc);
        s.setFromTriplets(tr.begin(), tr.end());
        return s;
    }

    // colSums(abs(M)) is not finite
    static bool bad_col_dense(const double* c, std::size_t n) {
        ldouble s = 0.0;
        for (std::size_t i = 0; i < n; ++i) s += std::fabs(c[i]);
        return !R_FINITE((double) s);
    }

    // iwls_solve(); err where R would signal an error
    void solve(const Pieces& pc, const std::vector<double>& u, double damp,
               const std::vector<char>& froz, Sol& sol, bool& err) {
        const int p = M.p;
        sol = Sol();
        sol.delta.assign(p, 0.0);
        if ((me.how == "qr" || me.how == "svd") && pc.sqrt) {
            sol.route = me.how;
            std::vector<char> bad(p, 0);
            for (int j = 0; j < p; ++j) {
                bool b = froz[j] || !R_FINITE(u[j]);
                if (!b) {
                    if (!M.sparse) b = bad_col_dense(pc.Rd.col(j), (std::size_t) pc.Rd.nr);
                    else {
                        ldouble s = 0.0;
                        for (SpMat::InnerIterator it(pc.Rs, j); it; ++it) s += std::fabs(it.value());
                        b = !R_FINITE((double) s);
                    }
                }
                if (!b && M.hasC) {
                    if (!M.Csp) b = bad_col_dense(M.Cd.col(j), (std::size_t) M.Cd.nr);
                    else {
                        ldouble s = 0.0;
                        for (SpMat::InnerIterator it(M.Cs, j); it; ++it) s += std::fabs(it.value());
                        b = !R_FINITE((double) s);
                    }
                }
                bad[j] = b;
            }
            std::vector<int> keep;
            for (int j = 0; j < p; ++j) { if (bad[j]) sol.held.push_back(j); else keep.push_back(j); }
            if (keep.empty()) { sol.rank = 0; return; }
            const int k = (int) keep.size();
            std::vector<double> uk(k);
            for (int j = 0; j < k; ++j) uk[j] = u[keep[j]];
            const int nR = M.sparse ? (int) pc.Rs.rows() : pc.Rd.nr;
            const int nC = M.hasC ? (M.Csp ? (int) M.Cs.rows() : M.Cd.nr) : 0;
            const int nD = damp > 0 ? k : 0;
            const double sd = damp > 0 ? std::sqrt(damp) : 0.0;
            std::vector<double> dk;
            std::vector<int> dropk;
            bool done = false;
            if (M.sparse && me.how == "qr" && nR + nC + nD >= k) {
                // sparse_augmented_solve() by the Cholesky factor of A'A
                // (D16): K = R'R + C'C is formed once per pieces, restricted
                // to the kept columns and damped here
                if (!pc.haveK) {
                    Timer tk(T_SPQR);
                    SpMat K = SpMat(pc.Rs.transpose() * pc.Rs);
                    if (M.hasC) {
                        SpMat Cs = M.Csp ? M.Cs : dense_to_sparse_mat(M.Cd);
                        K = K + SpMat(Cs.transpose() * Cs);
                    }
                    // the exact zeros of R (an information block that
                    // vanishes) are dropped, as cbind2() drops them on the R
                    // route, so that they do not fill the factor
                    K.prune([](int, int, double v) { return v != 0.0; });
                    pc.K = K;
                    pc.haveK = true;
                }
                if (chol_sparse(pc.K, keep, damp, uk, dk)) { done = true; sol.rank = k; }
            }
            if (!done) {
                Dense A(nR + nC + nD, k);
                for (int j = 0; j < k; ++j) {
                    double* c = A.col(j);
                    if (M.sparse) {
                        for (SpMat::InnerIterator it(pc.Rs, keep[j]); it; ++it)
                            c[it.row()] = it.value();
                    } else {
                        const double* r = pc.Rd.col(keep[j]);
                        for (int i = 0; i < nR; ++i) c[i] = r[i];
                    }
                    if (M.hasC) {
                        if (M.Csp) {
                            for (SpMat::InnerIterator it(M.Cs, keep[j]); it; ++it)
                                c[nR + it.row()] = it.value();
                        } else {
                            for (int i = 0; i < nC; ++i) c[nR + i] = M.Cd(i, keep[j]);
                        }
                    }
                    if (nD) c[nR + nC + j] = sd;
                }
                aug_dense(A, uk, dk, sol.rank, dropk, err);
                if (err) return;
            }
            for (int j = 0; j < k; ++j) sol.delta[keep[j]] = dk[j];
            for (int d : dropk) sol.dropped.push_back(keep[d]);
            return;
        }
        // the assembled route
        Timer tm(T_ASM);
        sol.route = (me.how == "qr" || me.how == "svd") ? me.how + "->chol" : me.how;
        std::vector<char> bad(p, 0);
        for (int j = 0; j < p; ++j)
            bad[j] = froz[j] || !R_FINITE(u[j]) || !R_FINITE(adiag(pc, j));
        Dense A0 = pc.asp ? Dense() : pc.Ad;
        for (;;) {
            std::vector<int> keep;
            for (int j = 0; j < p; ++j) if (!bad[j]) keep.push_back(j);
            if (keep.empty()) break;
            bool more = false;
            for (int b : keep) {
                ldouble s = 0.0;
                for (int a : keep) s += std::fabs(pc.asp ? pc.As.coeff(a, b) : pc.Ad(a, b));
                if (!R_FINITE((double) s)) { bad[b] = 1; more = true; }
            }
            if (!more) break;
        }
        std::vector<int> keep;
        for (int j = 0; j < p; ++j) { if (bad[j]) sol.held.push_back(j); else keep.push_back(j); }
        if (!sol.held.empty()) {
            if (keep.empty()) { sol.rank = 0; return; }
            const int k = (int) keep.size();
            Dense A(k, k);
            for (int b = 0; b < k; ++b)
                for (int a = 0; a < k; ++a)
                    A(a, b) = pc.asp ? pc.As.coeff(keep[a], keep[b]) : pc.Ad(keep[a], keep[b]);
            if (!pd_repair(A, me.floor)) { sol.rank = 0; return; }
            if (damp > 0) for (int j = 0; j < k; ++j) A(j, j) = A(j, j) + damp;
            Dense U;
            if (!r_chol(A, U)) { err = true; return; }
            std::vector<double> d(k);
            for (int j = 0; j < k; ++j) d[j] = u[keep[j]];
            if (!r_chol_solve(U, d)) { err = true; return; }
            for (int j = 0; j < k; ++j) sol.delta[keep[j]] = d[j];
            sol.rank = k;
            return;
        }
        if (pc.asp) {
            // a sparse curvature: the probe and the solve by a sparse
            // Cholesky, the repair on the dense matrix
            bool fin = true;
            for (int c = 0; c < pc.As.outerSize() && fin; ++c)
                for (SpMat::InnerIterator it(pc.As, c); it; ++it)
                    if (!R_FINITE(it.value())) { fin = false; break; }
            if (!fin) { sol.rank = 0; return; }
            Eigen::SimplicialLLT<SpMat, Eigen::Upper, Eigen::NaturalOrdering<int>> ch;
            ch.compute(pc.As);
            if (ch.info() == Eigen::Success) {
                SpMat A = pc.As;
                if (damp > 0) {
                    for (int j = 0; j < p; ++j) A.coeffRef(j, j) += damp;
                    ch.compute(A);
                    if (ch.info() != Eigen::Success) { err = true; return; }
                }
                Eigen::VectorXd uu = Eigen::Map<const Eigen::VectorXd>(u.data(), p);
                Eigen::VectorXd d = ch.solve(uu);
                for (int j = 0; j < p; ++j) sol.delta[j] = d[j];
                sol.rank = p;
                return;
            }
            A0 = dense_of(pc.As);
        }
        Dense A = A0;
        if (!pd_repair(A, me.floor)) { sol.rank = 0; return; }
        if (damp > 0) for (int j = 0; j < p; ++j) A(j, j) = A(j, j) + damp;
        Dense U;
        if (!r_chol(A, U)) { err = true; return; }
        std::vector<double> d(u);
        if (!r_chol_solve(U, d)) { err = true; return; }
        sol.delta = d;
        sol.rank = p;
    }

    // ------------------------------------------------------------ loop
    struct LineSearch {
        bool ok = false;
        std::vector<double> cand;
        double vnew = NA_REAL, step = 1.0;
        Sol sol;
        int errors = 0;
    };

    // iwls_line_search(): the solve is outside a tryCatch in R
    void line_search(const std::vector<double>& beta, double value,
                     const std::vector<double>& g, const Pieces& pc,
                     double damp, const std::vector<char>& froz, bool guard,
                     LineSearch& st) {
        const int p = M.p;
        std::vector<double> u(p);
        for (int j = 0; j < p; ++j) u[j] = -g[j];
        bool err = false;
        st = LineSearch();
        solve(pc, u, damp, froz, st.sol, err);
        if (err) throw RError();
        std::vector<double> delta = st.sol.delta;
        for (double& d : delta) if (!R_FINITE(d)) d = 0.0;
        std::vector<double> gdv(p);
        for (int j = 0; j < p; ++j) gdv[j] = g[j] * delta[j];
        const double gd = r_sum(gdv.data(), p);
        double step = 1.0;
        st.cand = beta;
        for (int h = 0; h < me.step_halving + 1; ++h) {
            for (int j = 0; j < p; ++j) st.cand[j] = beta[j] + step * delta[j];
            double vnew;
            if (!fn(st.cand, vnew)) {
                if (!guard) throw RError();
                st.errors += 1;
                vnew = R_PosInf;
            }
            st.vnew = vnew;
            if (R_FINITE(vnew) && vnew <= value + 1e-4 * step * gd) { st.ok = true; break; }
            step = step / 2;
        }
        st.step = step;
    }

    // iwls_decrement(): NA on an error
    double decrement(const Pieces& pc, const std::vector<double>& g,
                     double damp, const std::vector<char>& froz) {
        const int p = M.p;
        std::vector<double> u(p);
        for (int j = 0; j < p; ++j) u[j] = -g[j];
        Sol sol;
        bool err = false;
        try { solve(pc, u, damp, froz, sol, err); } catch (RError&) { err = true; }
        if (err) return NA_REAL;
        std::vector<double> v(p);
        for (int j = 0; j < p; ++j) {
            const double d = R_FINITE(sol.delta[j]) ? sol.delta[j] : 0.0;
            v[j] = (froz[j] ? 0.0 : g[j]) * d;
        }
        return -0.5 * r_sum(v.data(), p);
    }

    static bool slow(const std::vector<double>& s, double rate) {
        const std::size_t k = s.size();
        if (ISNAN(rate) || k < 3) return false;
        const double r1 = s[k - 2] / s[k - 3], r2 = s[k - 1] / s[k - 2];
        return R_FINITE(r1) && R_FINITE(r2) && r1 > rate && r2 > rate;
    }

    double escalate(double damp, const Pieces& pc) {
        return damp > 0 ? damp * 100 : scale(pc) * 1e-8;
    }

    static double max_abs_div(const std::vector<double>& g, double n) {
        double m = 0.0;
        for (double v : g) m = std::max(m, std::fabs(v / n));
        return m;
    }
    static double max_abs(const std::vector<double>& g) {
        double m = 0.0;
        for (double v : g) m = std::max(m, std::fabs(v));
        return m;
    }

    static std::string fmt_step(double x) {
        if (!R_FINITE(x)) return "-";
        char b[64];
        std::snprintf(b, sizeof b, "%.3g", x);
        return b;
    }
};

}  // namespace

// [[Rcpp::export]]
Rcpp::List engine_family_probe(std::string family, Rcpp::CharacterVector links) {
    Registry& r = registry();
    const int id = r.id(family.c_str());
    const int cov = id >= 0 ? r.covered(id) : 0;
    Rcpp::IntegerVector lids(links.size());
    for (int k = 0; k < links.size(); ++k)
        lids[k] = r.lid(Rcpp::as<std::string>(links[k]).c_str());
    return Rcpp::List::create(Rcpp::Named("id") = id, Rcpp::Named("covered") = cov,
                              Rcpp::Named("links") = lids);
}

namespace {

void read_sparse(SEXP x, SpMat& out) {
    Rcpp::S4 s(x);
    Rcpp::IntegerVector dim = s.slot("Dim"), ii = s.slot("i"), pp = s.slot("p");
    Rcpp::NumericVector xx = s.slot("x");
    const int nr = dim[0], nc = dim[1];
    std::vector<Trip> tr;
    tr.reserve(xx.size());
    for (int c = 0; c < nc; ++c)
        for (int k = pp[c]; k < pp[c + 1]; ++k) tr.emplace_back(ii[k], c, xx[k]);
    out = SpMat(nr, nc);
    out.setFromTriplets(tr.begin(), tr.end());
}

bool is_dgc(SEXP x) {
    return Rf_isS4(x) && Rf_inherits(x, "dgCMatrix");
}

void read_any(SEXP x, MatAny& m) {
    if (is_dgc(x)) {
        m.sp = true;
        read_sparse(x, m.s);
        m.nr = (int) m.s.rows();
        m.nc = (int) m.s.cols();
    } else {
        Rcpp::NumericMatrix d(x);
        m.sp = false;
        m.nr = d.nrow();
        m.nc = d.ncol();
        m.d = Dense(m.nr, m.nc);
        std::copy(d.begin(), d.end(), m.d.v.begin());
    }
}

}  // namespace

namespace {
// the model and the penalty of a call, read into an engine
void engine_setup(Engine& E, Rcpp::List model, Rcpp::List pen) {
    Model& M = E.M;
    Rcpp::List X = model["X"];
    M.K = X.size();
    Rcpp::NumericVector y = model["y"], w = model["w"];
    M.n = y.size();
    M.y = y.begin();
    M.w = w.begin();
    M.X.resize(M.K);
    M.npar.resize(M.K);
    M.offs.resize(M.K);
    M.sparse = Rcpp::as<bool>(model["sparse"]);
    int off = 0;
    for (int k = 0; k < M.K; ++k) {
        SEXP xk = X[k];
        Block& B = M.X[k];
        if (Rf_isNull(xk)) { B.m = 0; }
        else if (is_dgc(xk)) { B.sp = true; read_sparse(xk, B.s); B.m = (int) B.s.cols(); }
        else { B.d = REAL(xk); B.m = Rf_ncols(xk); }
        M.npar[k] = B.m;
        M.offs[k] = off;
        off += B.m;
    }
    M.p = off;
    Rcpp::List offl = model["offsets"];
    M.off.assign(M.K, nullptr);
    for (int k = 0; k < M.K; ++k) if (!Rf_isNull(offl[k])) M.off[k] = REAL(offl[k]);
    Registry& reg = registry();
    std::string fam = Rcpp::as<std::string>(model["family"]);
    M.fam = reg.id(fam.c_str());
    Rcpp::NumericMatrix cons = model["constants"];
    M.nc = cons.ncol();
    M.cons = cons.begin();
    Rcpp::NumericVector lo = model["lo"], hi = model["hi"];
    M.lo.assign(lo.begin(), lo.end());
    M.hi.assign(hi.begin(), hi.end());
    Rcpp::List links = model["links"];
    M.lid.resize(M.K);
    M.lpar.resize(M.K);
    for (int k = 0; k < M.K; ++k) {
        Rcpp::List lk = links[k];
        M.lid[k] = reg.lid(Rcpp::as<std::string>(lk["name"]).c_str());
        Rcpp::NumericVector pr = lk["par"];
        M.lpar[k].assign(pr.begin(), pr.end());
        Rcpp::NumericVector lb = lk["bounds"];
        M.llo.push_back(lb[0]);
        M.lhi.push_back(lb[1]);
    }
    M.pmax = reg.pmax();
    M.threads = Rcpp::as<int>(model["threads"]);
    E.thbuf.assign(M.K + M.nc + 1, 0.0);

    SEXP S = pen["S"];
    if (M.sparse) {
        if (is_dgc(S)) read_sparse(S, M.Ss);
        else {
            Rcpp::NumericMatrix Sm(S);
            std::vector<Trip> tr;
            for (int j = 0; j < Sm.ncol(); ++j)
                for (int i = 0; i < Sm.nrow(); ++i)
                    if (Sm(i, j) != 0.0) tr.emplace_back(i, j, Sm(i, j));
            M.Ss = SpMat(M.p, M.p);
            M.Ss.setFromTriplets(tr.begin(), tr.end());
        }
    } else {
        Rcpp::NumericMatrix Sm(S);
        M.Sd = Dense(M.p, M.p);
        std::copy(Sm.begin(), Sm.end(), M.Sd.v.begin());
    }
    M.gauss_id = reg.id("Gaussian1Distrib");
    Rcpp::List units = pen["units"];
    for (int k = 0; k < units.size(); ++k) {
        Rcpp::List ul = units[k];
        PenUnit u;
        u.kind = Rcpp::as<int>(ul["kind"]);
        Rcpp::IntegerVector ix = ul["idx"];
        for (int v : ix) u.idx.push_back(v - 1);
        SEXP mp = ul["map"];
        u.has_map = !Rf_isNull(mp);
        if (u.has_map) read_any(mp, u.map);
        if (u.kind == 1) {
            read_any(ul["P"], u.P);
            read_any(ul["M"], u.Mg);
            u.lam = Rcpp::as<double>(ul["lam"]);
            u.a = Rcpp::as<double>(ul["a"]);
            u.c1 = Rcpp::as<double>(ul["c1"]);
            u.c2 = Rcpp::as<double>(ul["c2"]);
            u.c3 = Rcpp::as<double>(ul["c3"]);
        } else if (u.kind == 2) {
            read_any(ul["S"], u.S);
            u.c1 = Rcpp::as<double>(ul["c1"]);
            u.c2 = Rcpp::as<double>(ul["c2"]);
            Rcpp::List ms = ul["matrices"];
            for (int h = 0; h < ms.size(); ++h) { MatAny mm; read_any(ms[h], mm); u.mats.push_back(mm); }
            Rcpp::NumericVector lm = ul["lams"];
            u.lams.assign(lm.begin(), lm.end());
        } else if (u.kind == 3) {
            u.mu = Rcpp::as<double>(ul["mu"]);
            u.sigma = Rcpp::as<double>(ul["sigma"]);
            u.lj = Rcpp::as<double>(ul["lj"]);
        } else {
            u.block = Rcpp::as<int>(ul["block"]);
            Rcpp::NumericVector mu = ul["mu"];
            u.muv.assign(mu.begin(), mu.end());
            Rcpp::NumericMatrix L = ul["L"], Si = ul["sinv"];
            u.L = Dense(L.nrow(), L.ncol());
            std::copy(L.begin(), L.end(), u.L.v.begin());
            u.Sinv = Dense(Si.nrow(), Si.ncol());
            std::copy(Si.begin(), Si.end(), u.Sinv.v.begin());
            u.inverted = Rcpp::as<bool>(ul["inverted"]);
            u.logdet = Rcpp::as<double>(ul["logdet"]);
            u.p2pi = Rcpp::as<double>(ul["p2pi"]);
            u.lj = Rcpp::as<double>(ul["lj"]);
        }
        M.units.push_back(u);
    }
    SEXP C = pen["C"];
    M.hasC = !Rf_isNull(C);
    M.Ccol2.assign(M.p, 0.0);
    if (M.hasC) {
        if (is_dgc(C)) {
            M.Csp = true;
            read_sparse(C, M.Cs);
            for (int j = 0; j < M.p; ++j) {
                ldouble s = 0.0;
                for (SpMat::InnerIterator it(M.Cs, j); it; ++it) s += it.value() * it.value();
                M.Ccol2[j] = (double) s;
            }
        } else {
            Rcpp::NumericMatrix Cm(C);
            M.Cd = Dense(Cm.nrow(), Cm.ncol());
            std::copy(Cm.begin(), Cm.end(), M.Cd.v.begin());
            for (int j = 0; j < M.p; ++j) {
                ldouble s = 0.0;
                for (int i = 0; i < M.Cd.nr; ++i) s += M.Cd(i, j) * M.Cd(i, j);
                M.Ccol2[j] = (double) s;
            }
        }
    }

}
}  // namespace

// The compiled inner fit. `model`: K design blocks (NULL, a base matrix or a
// dgCMatrix), the response, weights, offsets, the family's registry name,
// constants and parameter bounds, the links' registry names and parameters.
// `pen`: S, C (NULL for none) and the units. `method`: the settled iwls() as a
// list. Returns status 0 with iwls_fit()'s fields (the note as a code), or
// status 1 where the R route signals an error.
// [[Rcpp::export]]
Rcpp::List engine_iwls_fit(Rcpp::List model, Rcpp::List pen,
                           Rcpp::List method, Rcpp::NumericVector start,
                           Rcpp::IntegerVector frozen, Rcpp::List groups,
                           bool verbose) {
    for (int k = 0; k < T_N; ++k) g_time[k] = 0.0;
    Engine E;
    Model& M = E.M;
    engine_setup(E, model, pen);
    Method& me = E.me;
    me.how = Rcpp::as<std::string>(method["how"]);
    double st = Rcpp::as<double>(method["stol"]);
    me.has_stol = !ISNAN(st);
    me.stol = st;
    me.floor = Rcpp::as<double>(method["floor"]);
    me.maxit = Rcpp::as<int>(method["maxit"]);
    me.tol = Rcpp::as<double>(method["tol"]);
    me.step_halving = Rcpp::as<int>(method["step_halving"]);
    me.dec_limit = Rcpp::as<double>(method["decrement_limit"]);
    me.switch_after = Rcpp::as<double>(method["switch_after"]);
    me.switch_rate = Rcpp::as<double>(method["switch_rate"]);
    me.tol_def = Rcpp::as<double>(method["tol_def"]);
    me.kind_expected = Rcpp::as<bool>(method["expected"]);
    me.backup = Rcpp::as<bool>(method["backup"]);
    me.switching = Rcpp::as<bool>(method["switching"]);
    me.damp_on_reject = Rcpp::as<bool>(method["damp_on_reject"]);
    me.guard = Rcpp::as<bool>(method["guard"]);
    me.polish = Rcpp::as<bool>(method["polish"]);
    me.emethod = Rcpp::as<int>(method["emethod"]);
    me.escore = Rcpp::as<bool>(method["escore"]);
    me.nsim = Rcpp::as<double>(method["nsim"]);
    {
        Rcpp::NumericVector et = method["etol"], sp = method["support"];
        for (int k = 0; k < 6; ++k) me.etol[k] = et[k];
        me.support[0] = sp[0];
        me.support[1] = sp[1];
    }

    const int p = M.p;
    const double nn = (double) M.n;
    std::vector<char> froz(p, 0);
    for (int f : frozen) froz[f - 1] = 1;
    std::vector<std::vector<int>> grp;
    for (int k = 0; k < groups.size(); ++k) {
        Rcpp::IntegerVector gk = groups[k];
        std::vector<int> v;
        for (int j : gk) v.push_back(j - 1);
        grp.push_back(v);
    }
    if (grp.empty()) { std::vector<int> v(p); for (int j = 0; j < p; ++j) v[j] = j; grp.push_back(v); }

    // pieces_at, backup_at, switch_at as the kind of information each reads
    bool pa_expected = me.kind_expected;
    bool has_backup = me.backup, backup_expected = true;
    const bool has_switch = me.switching;
    double switch_after = has_switch ? me.switch_after : R_PosInf;

    std::vector<double> beta(start.begin(), start.end());
    const double tol_def = me.tol_def;
    int switched = 0;
    double value = 0.0;
    std::vector<int> aliased;
    bool converged = false;
    int it = 0;
    double score = R_PosInf;
    std::vector<double> scores;
    int note_kind = 0, note_it = 0;
    double dimless = NA_REAL;
    double damp = 0.0;
    int damp_tries = 0;
    int fb_indef = 0, fb_search = 0, fb_error = 0;
    double step_used = NA_REAL;
    // history
    std::vector<int> h_it, h_rank, h_held;
    std::vector<double> h_obj, h_score, h_step, h_damp;
    std::vector<std::string> h_route;

    try {
        if (!E.fn(beta, value)) throw RError();
        if (verbose) Rprintf("  %-5s %14s %12s %10s\n", "iter", "objective", "score/n", "step");
        typedef Engine::Pieces Pieces;
        for (it = 1; it <= me.maxit; ++it) {
            std::vector<double> g;
            if (!E.gr(beta, g)) throw RError();
            std::vector<double> gf(g);
            for (int j = 0; j < p; ++j) if (froz[j]) gf[j] = 0.0;
            bool fin = true;
            for (double v : gf) if (!R_FINITE(v)) fin = false;
            if (!fin) { note_kind = 1; note_it = it; break; }
            score = Engine::max_abs(gf) / nn;
            scores.push_back(score);
            if (!switched && has_switch &&
                ((double) it > switch_after || Engine::slow(scores, me.switch_rate))) {
                has_backup = true;
                backup_expected = pa_expected;
                pa_expected = false;
                switched = it;
            }
            if (verbose) {
                std::string s = it == 1 ? "-" : Engine::fmt_step(step_used);
                Rprintf("  %-5d %14.6f %12.3e %10s\n", it - 1, value, score, s.c_str());
            }
            bool met = Engine::max_abs_div(gf, nn) < me.tol;
            bool have_pc_met = false;
            Pieces pc_met;
            if (met) {
                if (!E.pieces(beta, pa_expected, pc_met)) throw RError();
                have_pc_met = true;
                double dec;
                if (has_backup && !E.definite(pc_met, tol_def)) {
                    Pieces pcd;
                    if (!E.pieces(beta, backup_expected, pcd)) throw RError();
                    dec = E.decrement(pcd, g, damp, froz);
                } else {
                    dec = E.decrement(pc_met, g, damp, froz);
                }
                if (R_FINITE(dec) && dec > me.dec_limit) met = false;
            }
            if (met) {
                converged = true;
                if (me.polish) {
                Pieces pp_own;
                const Pieces* ppp = &pc_met;
                bool have_pp = true;
                if (has_switch && !switched) {
                    if (!E.pieces(beta, false, pp_own)) throw RError();
                    ppp = &pp_own;
                }
                const Pieces& pp = *ppp;
                if (!E.definite(pp, tol_def)) have_pp = false;
                if (have_pp) {
                    std::vector<double> u(p);
                    for (int j = 0; j < p; ++j) u[j] = -g[j];
                    Engine::Sol sol;
                    bool err = false;
                    try { E.solve(pp, u, damp, froz, sol, err); } catch (RError&) { err = true; }
                    bool f = !err;
                    if (f) for (double d : sol.delta) if (!R_FINITE(d)) f = false;
                    if (f) {
                        std::vector<double> cand(p);
                        for (int j = 0; j < p; ++j) cand[j] = beta[j] + sol.delta[j];
                        double vnew;
                        if (!E.fn(cand, vnew)) vnew = R_PosInf;
                        if (R_FINITE(vnew) && vnew <= value + 1e-8 * std::max(1.0, std::fabs(value))) {
                            beta = cand;
                            value = vnew;
                        }
                    }
                }
                }
                break;
            }
            Pieces pc_own;
            Pieces* pcp = &pc_met;
            if (!have_pc_met) {
                if (!E.pieces(beta, pa_expected, pc_own)) throw RError();
                pcp = &pc_own;
            }
            if (has_backup && !E.definite(*pcp, tol_def)) fb_indef += 1;
            Engine::LineSearch st;
            E.line_search(beta, value, g, *pcp, damp, froz, me.guard || has_backup, st);
            fb_error += st.errors;
            if (!st.ok && has_backup) {
                if (!E.pieces(beta, backup_expected, pc_own)) throw RError();
                pcp = &pc_own;
                fb_search += 1;
                E.line_search(beta, value, g, *pcp, damp, froz, true, st);
                fb_error += st.errors;
            }
            aliased = st.sol.dropped;
            step_used = st.step;
            if (!st.ok) {
                if (me.damp_on_reject && damp_tries < 8) {
                    damp = E.escalate(damp, *pcp);
                    damp_tries += 1;
                    continue;
                }
                note_kind = 2;
                note_it = it;
                break;
            }
            h_it.push_back(it);
            h_obj.push_back(st.vnew);
            h_score.push_back(score);
            h_step.push_back(st.step);
            h_rank.push_back(st.sol.rank);
            h_route.push_back(st.sol.route);
            h_held.push_back((int) st.sol.held.size());
            h_damp.push_back(damp);
            const bool stalled = value - st.vnew <= 1e-12 * std::max(1.0, std::fabs(value));
            beta = st.cand;
            value = st.vnew;
            if (stalled) {
                if (st.step < 1e-3 && damp_tries < 8) {
                    damp = E.escalate(damp, *pcp);
                    damp_tries += 1;
                    continue;
                }
                note_kind = 3;
                note_it = it;
                break;
            }
            if (damp > 0) {
                damp = damp / 100;
                if (damp < E.scale(*pcp) * 1e-12) damp = 0;
            }
        }
        if (it > me.maxit) it = me.maxit;
        if (!converged) {
            std::vector<double> gfin;
            if (!E.gr(beta, gfin)) throw RError();
            for (int j = 0; j < p; ++j) if (froz[j]) gfin[j] = 0.0;
            bool fin = true;
            for (double v : gfin) if (!R_FINITE(v)) fin = false;
            if (fin) {
                score = Engine::max_abs(gfin) / nn;
                converged = Engine::max_abs_div(gfin, nn) < me.tol;
                if (!converged) {
                    Pieces pf;
                    if (!E.pieces(beta, pa_expected, pf)) throw RError();
                    std::vector<double> hj = E.info_diag(pf);
                    double out = 0.0;
                    for (const auto& idx : grp) {
                        if (idx.empty()) continue;
                        std::vector<double> h;
                        for (int j : idx) h.push_back(hj[j]);
                        const double m = r_median(h);
                        const double s = (R_FINITE(m) && m > 0) ? std::sqrt(m / nn) : 1.0;
                        double mg = 0.0;
                        for (int j : idx) mg = std::max(mg, std::fabs(gfin[j]));
                        out = std::max(out, mg / (nn * s));
                    }
                    if (out < me.tol) { converged = true; dimless = out; note_kind = 4; }
                }
            } else {
                score = R_PosInf;
            }
        }
    } catch (RError&) {
        return Rcpp::List::create(Rcpp::Named("status") = 1);
    }

    Rcpp::IntegerVector al(aliased.size());
    for (std::size_t j = 0; j < aliased.size(); ++j) al[j] = aliased[j] + 1;
    SEXP hist = R_NilValue;
    if (!h_it.empty()) {
        hist = Rcpp::List::create(
            Rcpp::Named("iteration") = Rcpp::wrap(h_it),
            Rcpp::Named("objective") = Rcpp::wrap(h_obj),
            Rcpp::Named("score") = Rcpp::wrap(h_score),
            Rcpp::Named("step") = Rcpp::wrap(h_step),
            Rcpp::Named("rank") = Rcpp::wrap(h_rank),
            Rcpp::Named("route") = Rcpp::wrap(h_route),
            Rcpp::Named("held") = Rcpp::wrap(h_held),
            Rcpp::Named("damp") = Rcpp::wrap(h_damp));
    }
    return Rcpp::List::create(
        Rcpp::Named("status") = 0,
        Rcpp::Named("par") = Rcpp::wrap(beta),
        Rcpp::Named("value") = value,
        Rcpp::Named("converged") = converged,
        Rcpp::Named("iterations") = it,
        Rcpp::Named("score") = score,
        Rcpp::Named("note_kind") = note_kind,
        Rcpp::Named("note_it") = note_it,
        Rcpp::Named("dimless") = dimless,
        Rcpp::Named("aliased") = al,
        Rcpp::Named("fallback") = Rcpp::IntegerVector::create(
            Rcpp::Named("indefinite") = fb_indef, Rcpp::Named("search") = fb_search,
            Rcpp::Named("error") = fb_error, Rcpp::Named("newton") = 0),
        Rcpp::Named("switched") = switched,
        Rcpp::Named("history") = hist,
        Rcpp::Named("oob") = Rcpp::wrap(E.oob),
        Rcpp::Named("timing") = Rcpp::NumericVector::create(
            Rcpp::Named("point") = g_time[T_POINT], Rcpp::Named("parts") = g_time[T_PARTS],
            Rcpp::Named("pieces") = g_time[T_PIECES], Rcpp::Named("sparse_qr") = g_time[T_SPQR],
            Rcpp::Named("dense_qr") = g_time[T_DQR], Rcpp::Named("assembled") = g_time[T_ASM],
            Rcpp::Named("fn") = g_time[T_FN]));
}

// The engine's objective, gradient and parameters at one point, for the
// twin tests: what obj$fn(), obj$gr() and statmod_eta() return on the R route.
// [[Rcpp::export]]
Rcpp::List engine_point_probe(Rcpp::List model, Rcpp::List pen, Rcpp::NumericVector beta) {
    Engine E;
    engine_setup(E, model, pen);
    std::vector<double> b(beta.begin(), beta.end());
    double value = NA_REAL;
    const bool okv = E.fn(b, value);
    std::vector<double> g;
    const bool okg = E.gr(b, g);
    Point& P = E.point(b);
    return Rcpp::List::create(Rcpp::Named("ok") = okv && okg, Rcpp::Named("value") = value,
                              Rcpp::Named("gradient") = Rcpp::wrap(g),
                              Rcpp::Named("theta") = Rcpp::wrap(P.th),
                              Rcpp::Named("h1") = Rcpp::wrap(P.h1),
                              Rcpp::Named("eta") = Rcpp::wrap(P.eta));
}
