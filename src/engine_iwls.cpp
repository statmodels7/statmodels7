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
    double rank = -1;                // kind 1: the penalty's rank (phase 3)
    int chart = 0;                   // kind 4 (phase 3): 1 log-Cholesky, 2 compound
                                     // symmetry, 3 AR(1); 0 not described
    // a chart that the context does not write (code 9): R's derivatives of
    // Sinv and of log|Sigma| in the free values (mv_matrix_pieces())
    std::vector<Dense> gai, gai2;    // K, K*K (t * K + s)
    std::vector<double> gdld, gd2ld; // K, K*K
    std::vector<double> ceta, ch;    // the chart's free values; CS and AR(1): the
                                     // scale's and the correlation's h0, h1, h2
    int add_m = 0;                   // kind 2 (phase 3): the blocks, 0 not described
    double add_tol = 1e-10;          // its rank tolerance
    std::vector<Dense> add_block;    // one block of each component
    // kinds 5-8 (phase 3b), the penalties with a kink: 5 the lasso (a
    // laplace2 prior at zero), 6 the elastic net (an enet prior at zero),
    // 7 SCAD, 8 MCP; a map is NULL or diagonal, held as its diagonal `dg`
    std::vector<double> dg;          // the map's diagonal, empty without a map
    double klam = 0;                 // lambda
    double kL = 0;                   // 5: log(lambda / 2)
    double ka = 0, kc = 0, klogz = 0;// 6: lambda alpha, lambda (1 - alpha), log Z
    double kpar = 0;                 // 7: a; 8: gamma
    std::vector<double> cv;          // 7, 8: the curvature (ones without one)
    // kind 9 (phase 3b): a DistribPenalty with a univariate Student t prior of
    // held location (fixed() of student_t1_distrib()), mu and sigma above,
    // without a map; its Hessian moves with the coefficients
    double tnu = 0;
    std::vector<int> tfree;          // outer: the hyperparameters' positions in
                                     // (sigma, nu), 1 or 2
    bool bq = false;                 // outer: beta_quadratic() at the point
    // kind 10 (phase 3b): a multivariate Student t prior of held location
    // (fixed() of an MvStudentTDistrib, not inverted, without a map); block,
    // muv, L, Sinv and tnu as above, and the value's constant
    // (cp0 - (p/2) log(pi)) - logdet/2 and (nu + p)/2
    double mc0 = 0, mhalf = 0;
    // outer (kind 10): the penalty's derivatives in its H hyperparameters at
    // the mode, from R (penalty_grad_theta() ... penalty_d2hessian()), the
    // dS blocks group by group, and the chart's derivatives of Sinv with the
    // derivative of N = nu + p, per hyperparameter
    bool mo = false;
    int mH = 0;
    std::vector<double> mgt, mrho2, mc, mc2;     // H, H*H, H*m, H*H*m
    std::vector<double> mdS, md2S;               // H*G*p*p, H*H*G*p*p
    std::vector<Dense> mai;                      // H
    std::vector<double> mNa;                     // H
};

inline bool kinked_kind(int k) { return k >= 5 && k <= 8; }

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
    // a penalty whose Hessian moves with the coefficients (kind 9): S and C
    // are those of the last point update_penalty() was called at
    bool moving = false;
    std::vector<double> Sbeta;       // that point
    // a sparse S: the entries of the units that do not move, kept once
    bool ss0 = false;
    SpMat Ss0;
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
        if (kinked_kind(u.kind)) {
            // a diagonal map: Diagonal %*% b is the elementwise product
            t.resize(bu.size());
            for (std::size_t j = 0; j < bu.size(); ++j) t[j] = u.dg[j] * bu[j];
            return;
        }
        t.assign(u.map.nr, 0.0);
        mat_vec(u.map, bu.data(), t.data());
    }

    // R's sign()
    static double rsign(double x) {
        if (ISNAN(x)) return x;
        return x > 0 ? 1.0 : (x < 0 ? -1.0 : 0.0);
    }

    // penalty_value() of a kinked unit at its mapped coefficients t: the
    // lasso and the elastic net as a DistribPenalty (minus the sum of the log
    // density of their fixed() parent, minus the map's log-Jacobian), SCAD and
    // MCP as their own methods (scad_mcp.R), operation by operation
    double kink_value(const PenUnit& u, const std::vector<double>& t) {
        const std::size_t m = t.size();
        std::vector<double> e(m);
        if (u.kind == 5 || u.kind == 6) {
            for (std::size_t j = 0; j < m; ++j) {
                const double z = t[j] - 0.0;
                if (u.kind == 5) e[j] = u.kL - u.klam * std::fabs(z);
                else e[j] = ((-u.ka) * std::fabs(z) - (u.kc * (z * z)) / 2) - u.klogz;
            }
            const double v = -r_sum(e.data(), m);
            if (!u.has_map || ISNAN(u.lj)) return v;
            return v - u.lj;
        }
        for (std::size_t j = 0; j < m; ++j) {
            const double c = u.cv[j], lam = u.klam / c, a = u.kpar, x = std::fabs(t[j]);
            double val;
            if (u.kind == 7) {
                if (x <= lam) val = lam * x;
                else if (!(x > a * lam)) val = (((2 * a) * lam) * x - x * x - lam * lam) / (2 * (a - 1));
                else val = ((lam * lam) * (a + 1)) / 2;
            } else {
                if (x <= a * lam) val = lam * x - (x * x) / (2 * a);
                else val = (a * (lam * lam)) / 2;
            }
            e[j] = c * val;
        }
        return r_sum(e.data(), m);
    }

    // penalty_gradient() of a kinked unit, on its coefficients (map_back())
    void kink_gradient(const PenUnit& u, const std::vector<double>& t,
                       std::vector<double>& g) {
        const std::size_t m = t.size();
        for (std::size_t j = 0; j < m; ++j) {
            double x;
            if (u.kind == 5 || u.kind == 6) {
                // -map_back(distrib_grad_y(parent, t) + 0 * t)
                const double z = t[j] - 0.0, s = rsign(z);
                double gy = u.kind == 5 ? (-u.klam) * s : (-u.ka) * s - u.kc * z;
                gy = gy + 0.0 * t[j];
                x = u.has_map ? -(u.dg[j] * gy) : -gy;
            } else {
                const double c = u.cv[j], lam = u.klam / c, a = u.kpar,
                             ab = std::fabs(t[j]), s = rsign(t[j]);
                double d1;
                if (u.kind == 7) {
                    if (ab <= lam) d1 = lam;
                    else if (!(ab > a * lam)) d1 = (a * lam - ab) / (a - 1);
                    else d1 = 0.0;
                } else {
                    d1 = ab <= a * lam ? lam - ab / a : 0.0;
                }
                const double y = (s * c) * d1;
                x = u.has_map ? u.dg[j] * y : y;
            }
            g[j] = x;
        }
    }

    // the diagonal of penalty_hessian() of a kinked unit (map_quad() is
    // crossprod(D, D * h), so a mapped entry is d * (d * h))
    void kink_hessian(const PenUnit& u, const std::vector<double>& t,
                      std::vector<double>& h) {
        const std::size_t m = t.size();
        h.resize(m);
        for (std::size_t j = 0; j < m; ++j) {
            double in;
            bool neg = false;
            if (u.kind == 5) { in = 0.0 + 0.0 * t[j]; neg = true; }
            else if (u.kind == 6) { in = -u.kc + 0.0 * t[j]; neg = true; }
            else {
                const double c = u.cv[j], lam = u.klam / c, a = u.kpar, ab = std::fabs(t[j]);
                if (u.kind == 7) {
                    const bool r2 = !(ab <= lam) && !(ab > a * lam);
                    in = c * (r2 ? -1 / (a - 1) : 0.0);
                } else {
                    in = c * (ab <= a * lam ? -1 / a : 0.0);
                }
            }
            double v = u.has_map ? u.dg[j] * (u.dg[j] * in) : in;
            h[j] = neg ? -v : v;
        }
    }

    // mvt_yscaled() of one group's coefficients y (p of them): the residual
    // divided by m = max(1, max |r|), ws = (rt Sinv) / (m st), is, cw
    static void mvt_yscaled(const PenUnit& u, const double* y, std::vector<double>& ws,
                            double& is, double& cw) {
        const int p = u.block;
        std::vector<double> r(p), rt(p), w2(p), sq(p);
        double mx = 1.0;
        for (int j = 0; j < p; ++j) {
            r[j] = y[j] - u.muv[j];
            mx = std::max(mx, std::fabs(r[j]));
        }
        for (int j = 0; j < p; ++j) rt[j] = r[j] / mx;
        for (int j = 0; j < p; ++j) w2[j] = rt[j];
        if (!r_trisolve(u.L, false, w2)) throw RError();
        for (int j = 0; j < p; ++j) sq[j] = w2[j] * w2[j];
        const double qt = r_colsum(sq.data(), p);
        const double im = 1 / mx, im2 = im * im;
        const double st = u.tnu * im2 + qt;
        is = im2 / st;
        cw = (u.tnu + p) * is;
        ws.assign(p, 0.0);
        for (int c = 0; c < p; ++c) {
            double w = 0.0;
            for (int l = 0; l < p; ++l) w = w + rt[l] * u.Sinv(l, c);
            ws[c] = w / (mx * st);
        }
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
        if (kinked_kind(u.kind)) return kink_value(u, t);
        double v;
        if (u.kind == 10) {
            // distrib_pdf() of the multivariate t, group by group
            const int p = u.block, G = (int) (m / p);
            std::vector<double> ld(G), w(p), sq(p);
            for (int g = 0; g < G; ++g) {
                for (int j = 0; j < p; ++j) w[j] = t[(std::size_t) g * p + j] - u.muv[j];
                if (!r_trisolve(u.L, false, w)) throw RError();
                for (int j = 0; j < p; ++j) sq[j] = w[j] * w[j];
                const double q = r_colsum(sq.data(), p);
                ld[g] = u.mc0 - u.mhalf * std::log1p(q / u.tnu);
            }
            return -r_sum(ld.data(), G);
        }
        if (u.kind == 9) {
            // distrib_pdf() of the t: dt((y - mu) / sigma, nu, log) - log(sigma)
            std::vector<double> lp(m);
            for (std::size_t i = 0; i < m; ++i)
                lp[i] = R::dt((t[i] - u.mu) / u.sigma, u.tnu, 1) - std::log(u.sigma);
            return -r_sum(lp.data(), m);
        }
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
        if (kinked_kind(u.kind)) {
            kink_gradient(u, t, g);
        } else if (u.kind == 1) {
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
            if (u.kind == 10) {
                // distrib_grad_y(): -(nu + p) * ws of mvt_yscaled(), + 0 * a
                const int p = u.block, G = (int) (m / p);
                std::vector<double> ws(p);
                double is, cw;
                for (int gg = 0; gg < G; ++gg) {
                    mvt_yscaled(u, t.data() + (std::size_t) gg * p, ws, is, cw);
                    for (int j = 0; j < p; ++j) {
                        const std::size_t e = (std::size_t) gg * p + j;
                        gy[e] = (-(u.tnu + p)) * ws[j] + 0.0 * t[e];
                    }
                }
            } else if (u.kind == 9) {
                // student_t1_dy1_cpp(), + 0 * a
                const double s = u.sigma, ik = 1.0 / (u.tnu - 0.0);
                for (std::size_t i = 0; i < m; ++i) {
                    const double z = (t[i] - u.mu) / s;
                    const double q = z * z * ik;
                    const double tt = 1.0 / (1.0 + q);
                    const double zt = (std::fabs(z) <= 1.0) ? z * tt : 1.0 / (1.0 / z + z * ik);
                    gy[i] = -zt*(ik + 1)/s + 0.0 * t[i];
                }
            } else if (u.kind == 3) {
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

    // ------------------------------------------------------- the subset
    // fit_smooth() on a subset (phase 3b): the loop works on the positions
    // `sub` (0-based, ascending) with every other coordinate held at `base`;
    // empty `sub` is the whole vector
    std::vector<int> sub;
    std::vector<double> base;
    int np() const { return sub.empty() ? M.p : (int) sub.size(); }
    void expand(const std::vector<double>& b, std::vector<double>& v) const {
        if (sub.empty()) { v = b; return; }
        v = base;
        for (std::size_t j = 0; j < sub.size(); ++j) v[sub[j]] = b[j];
    }
    // the full objective at the expanded point, as fit_smooth()'s fn()
    bool fn_s(const std::vector<double>& b, double& value) {
        if (sub.empty()) return fn(b, value);
        std::vector<double> v;
        expand(b, v);
        return fn(v, value);
    }
    // obj$gr_sub(): the score over the subset's columns and the stacked
    // penalty gradient read at them; every entry of obj$gr() is its own column
    // product, so the subset of the full gradient is the same double
    bool gr_s(const std::vector<double>& b, std::vector<double>& g) {
        if (sub.empty()) return gr(b, g);
        std::vector<double> v, gf;
        expand(b, v);
        if (!gr(v, gf)) return false;
        g.resize(sub.size());
        for (std::size_t j = 0; j < sub.size(); ++j) g[j] = gf[sub[j]];
        return true;
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

    // penalty_sqrt(S) of a dense S: the zero matrix gives no rows, a diagonal
    // one diagonal_sqrt(), any other the eigen route on (S + t(S)) / 2; false
    // (NULL) when S is indefinite beyond the tolerance
    bool penalty_sqrt_dense(const Dense& S, Dense& C) {
        const int p = S.nc;
        if (p == 0) { C = Dense(0, 0); return true; }
        bool zero = true, diag = true;
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < p; ++i) {
                const double v = S(i, j);
                if (!(v == 0.0)) { zero = false; if (i != j) diag = false; }
            }
        if (zero) { C = Dense(0, p); return true; }
        std::vector<int> keep;
        if (diag) {
            // diagonal_sqrt()
            std::vector<double> d(p);
            double mx = 0.0, mn = R_PosInf;
            for (int j = 0; j < p; ++j) {
                d[j] = S(j, j);
                mx = std::max(mx, std::fabs(d[j]));
                mn = std::min(mn, d[j]);
            }
            const double tol = p * DBL_EPSILON * mx;
            if (mn < -std::max(tol, 1e-10 * mx)) return false;
            for (int j = 0; j < p; ++j) if (d[j] > tol) keep.push_back(j);
            C = Dense((int) keep.size(), p);
            for (std::size_t r = 0; r < keep.size(); ++r)
                C((int) r, keep[r]) = std::sqrt(d[keep[r]]);
            return true;
        }
        Dense A(p, p);
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < p; ++i) A(i, j) = (S(i, j) + S(j, i)) / 2;
        std::vector<double> ev;
        Dense V;
        if (!r_eigen_sym(A, false, ev, V)) throw RError();
        double mx = 0.0, mn = R_PosInf;
        for (double v : ev) { mx = std::max(mx, std::fabs(v)); mn = std::min(mn, v); }
        const double tol = p * DBL_EPSILON * mx;
        if (mn < -std::max(tol, 1e-10 * mx)) return false;
        for (int k = 0; k < p; ++k) if (ev[k] > tol) keep.push_back(k);
        C = Dense((int) keep.size(), p);
        for (std::size_t r = 0; r < keep.size(); ++r) {
            const double sq = std::sqrt(ev[keep[r]]);
            for (int j = 0; j < p; ++j) C((int) r, j) = sq * V(j, keep[r]);
        }
        return true;
    }

    // statmod_penalty_at(..., "hessian") and penalty_sqrt() at beta, for a
    // penalty that moves: the moving units' diagonal entries, 0 + their
    // block, are written into S (no other unit reaches their positions), and
    // C is formed anew
    void update_penalty(const std::vector<double>& beta) {
        if (!M.moving || same_beta(M.Sbeta, beta)) return;
        // the moving units' entries, 0 + their block; on a dense S written
        // in place, on a sparse one (statmod_penalty_at()'s accumulator
        // follows the design) S is the other units' entries and these, the
        // exact zeros not stored, as the subassignment of the R route stores
        // them
        std::vector<Trip> mv;
        auto put = [&](int i, int j, double v) {
            if (M.sparse) mv.emplace_back(i, j, v);
            else M.Sd(i, j) = v;
        };
        for (const PenUnit& u : M.units) {
            if (u.kind == 10) {
                // -dp_blockdiag(distrib_hess_y()): the block of group g is
                // -((-cw) Sinv + (2 N) ws ws')
                const int p = u.block, G = (int) (u.idx.size() / p);
                const double N = u.tnu + p;
                std::vector<double> y(p), ws(p);
                double is, cw;
                for (int g = 0; g < G; ++g) {
                    for (int j = 0; j < p; ++j) y[j] = beta[u.idx[(std::size_t) g * p + j]];
                    mvt_yscaled(u, y.data(), ws, is, cw);
                    for (int b = 0; b < p; ++b)
                        for (int a = 0; a < p; ++a) {
                            const double h = ((-cw) * u.Sinv(a, b)) + ((2 * N) * (ws[a] * ws[b]));
                            put(u.idx[(std::size_t) g * p + a], u.idx[(std::size_t) g * p + b],
                                0.0 + (-h));
                        }
                }
                continue;
            }
            if (u.kind != 9) continue;
            const double s = u.sigma, ik = 1.0 / (u.tnu - 0.0);
            for (std::size_t j = 0; j < u.idx.size(); ++j) {
                const double b = beta[u.idx[j]];
                const double z = (b - u.mu) / s;
                const double q = z * z * ik;
                const double t = 1.0 / (1.0 + q);
                const double uu = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
                const double h = t*(-ik*t + ik*uu - t + uu)/std::pow(s, 2);
                put(u.idx[j], u.idx[j], 0.0 + (-(h + 0.0 * b)));
            }
        }
        if (M.sparse) {
            if (!M.ss0) {
                std::vector<char> moving_pos(M.p, 0);
                for (const PenUnit& u : M.units)
                    if (u.kind == 9 || u.kind == 10) for (int i : u.idx) moving_pos[i] = 1;
                std::vector<Trip> keep;
                for (int j = 0; j < M.Ss.outerSize(); ++j)
                    for (SpMat::InnerIterator it(M.Ss, j); it; ++it)
                        if (!moving_pos[it.row()] && !moving_pos[j])
                            keep.emplace_back((int) it.row(), j, it.value());
                M.Ss0 = SpMat(M.p, M.p);
                M.Ss0.setFromTriplets(keep.begin(), keep.end());
                M.ss0 = true;
            }
            std::vector<Trip> all;
            for (int j = 0; j < M.Ss0.outerSize(); ++j)
                for (SpMat::InnerIterator it(M.Ss0, j); it; ++it)
                    all.emplace_back((int) it.row(), j, it.value());
            for (const Trip& t : mv) if (t.value() != 0.0) all.push_back(t);
            M.Ss = SpMat(M.p, M.p);
            M.Ss.setFromTriplets(all.begin(), all.end());
        }
        if (me.how == "qr" || me.how == "svd") {
            Dense C;
            M.hasC = penalty_sqrt_dense(M.sparse ? dense_of(M.Ss) : M.Sd, C);
            M.Ccol2.assign(M.p, 0.0);
            M.Csp = false;
            if (M.hasC && M.sparse) {
                // a diagonal S stored sparse gives a sparse factor
                // (diagonal_sqrt()), any other S a dense one
                bool diag = true;
                for (int j = 0; j < M.Ss.outerSize() && diag; ++j)
                    for (SpMat::InnerIterator it(M.Ss, j); it; ++it)
                        if (it.row() != j && it.value() != 0.0) { diag = false; break; }
                if (diag && M.Ss.nonZeros() > 0) {
                    M.Csp = true;
                    M.Cs = dense_to_sparse_mat(C);
                }
            }
            if (M.hasC) {
                M.Cd = C;
                for (int j = 0; j < M.p; ++j) {
                    ldouble sq = 0.0;
                    for (int i = 0; i < M.Cd.nr; ++i) sq += M.Cd(i, j) * M.Cd(i, j);
                    M.Ccol2[j] = (double) sq;
                }
            }
        }
        M.Sbeta = beta;
    }

    // iwls_pieces(): false where R signals an error
    bool pieces(const std::vector<double>& beta, bool expected, Pieces& pc) {
        Point& P = point(beta);
        if (!P.inb) return false;
        update_penalty(beta);
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

    // subset_pieces(): H over the subset's columns (information_subset(),
    // one product per pair of equations and column pair, the diagonal blocks
    // in full) plus S[sub, sub], always on the assembled route
    bool pieces_s(const std::vector<double>& b, bool expected, Pieces& pc) {
        if (sub.empty()) return pieces(b, expected, pc);
        std::vector<double> beta;
        expand(b, beta);
        Point& P = point(beta);
        if (!P.inb) return false;
        update_penalty(beta);
        if (expected && me.emethod == 0) ensure_eparts(P); else ensure_parts(P);
        if (expected && me.emethod >= 2) expect_all(P);
        Timer tm(T_PIECES);
        const int n = M.n, K = M.K, ps = (int) sub.size();
        std::vector<double> Hall((std::size_t) n * K * K);
        for (int i = 0; i < n; ++i)
            hess_obs(P, i, expected, Hall.data() + (std::size_t) i * K * K);
        pc = Pieces();
        pc.sqrt = false;
        std::vector<int> pos(M.p, -1);
        for (int j = 0; j < ps; ++j) pos[sub[j]] = j;
        pc.Ad = Dense(ps, ps);
        std::vector<double> wv(n), xw(n);
        for (int a = 0; a < K; ++a) {
            const Block& Xa = M.X[a];
            if (Xa.m == 0) continue;
            for (int bq = a; bq < K; ++bq) {
                const Block& Xb = M.X[bq];
                if (Xb.m == 0) continue;
                for (int i = 0; i < n; ++i)
                    wv[i] = M.w[i] * Hall[(std::size_t) i * K * K + a + bq * K];
                for (int j = 0; j < Xa.m; ++j) {
                    const int pj = pos[M.offs[a] + j];
                    if (pj < 0) continue;
                    const double* xj = Xa.d + (std::size_t) j * n;
                    for (int i = 0; i < n; ++i) xw[i] = xj[i] * wv[i];
                    for (int k = 0; k < Xb.m; ++k) {
                        const int pk = pos[M.offs[bq] + k];
                        if (pk < 0) continue;
                        const double* xk = Xb.d + (std::size_t) k * n;
                        double t = 0.0;
                        for (int i = 0; i < n; ++i) t = t + xw[i] * xk[i];
                        pc.Ad(pj, pk) = -t;
                        if (bq > a) pc.Ad(pk, pj) = -t;
                    }
                }
            }
        }
        for (int j = 0; j < ps; ++j)
            for (int i = 0; i < ps; ++i) pc.Ad(i, j) = pc.Ad(i, j) + M.Sd(sub[i], sub[j]);
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
        for (int j = 0; j < np(); ++j) {
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
        std::vector<double> d(np());
        for (int j = 0; j < np(); ++j) d[j] = pc.sqrt ? col2_R(pc, j) : adiag(pc, j);
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
        const int p = np();
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
        const int p = np();
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
            if (!fn_s(st.cand, vnew)) {
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
        const int p = np();
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
        if (u.has_map && !kinked_kind(u.kind)) read_any(mp, u.map);
        if (kinked_kind(u.kind)) {
            if (u.has_map) {
                Rcpp::NumericVector d(mp);
                u.dg.assign(d.begin(), d.end());
            }
            u.klam = Rcpp::as<double>(ul["lambda"]);
            u.lj = Rcpp::as<double>(ul["lj"]);
            if (u.kind == 5) u.kL = Rcpp::as<double>(ul["L"]);
            if (u.kind == 6) {
                u.ka = Rcpp::as<double>(ul["a"]);
                u.kc = Rcpp::as<double>(ul["c"]);
                u.klogz = Rcpp::as<double>(ul["logz"]);
            }
            if (u.kind >= 7) {
                u.kpar = Rcpp::as<double>(ul["par"]);
                Rcpp::NumericVector cvv = ul["cv"];
                u.cv.assign(cvv.begin(), cvv.end());
            }
        } else if (u.kind == 1) {
            read_any(ul["P"], u.P);
            read_any(ul["M"], u.Mg);
            u.lam = Rcpp::as<double>(ul["lam"]);
            u.a = Rcpp::as<double>(ul["a"]);
            u.c1 = Rcpp::as<double>(ul["c1"]);
            u.c2 = Rcpp::as<double>(ul["c2"]);
            u.c3 = Rcpp::as<double>(ul["c3"]);
            if (ul.containsElementNamed("rank")) u.rank = Rcpp::as<double>(ul["rank"]);
        } else if (u.kind == 2) {
            read_any(ul["S"], u.S);
            u.c1 = Rcpp::as<double>(ul["c1"]);
            u.c2 = Rcpp::as<double>(ul["c2"]);
            Rcpp::List ms = ul["matrices"];
            for (int h = 0; h < ms.size(); ++h) { MatAny mm; read_any(ms[h], mm); u.mats.push_back(mm); }
            Rcpp::NumericVector lm = ul["lams"];
            u.lams.assign(lm.begin(), lm.end());
            if (ul.containsElementNamed("add")) {
                Rcpp::List al = ul["add"];
                u.add_m = Rcpp::as<int>(al["m"]);
                u.rank = Rcpp::as<double>(al["rank"]);
                u.add_tol = Rcpp::as<double>(al["tol"]);
                Rcpp::List bm = al["block_mats"];
                for (int h = 0; h < bm.size(); ++h) {
                    Rcpp::NumericMatrix B(bm[h]);
                    Dense D(B.nrow(), B.ncol());
                    std::copy(B.begin(), B.end(), D.v.begin());
                    u.add_block.push_back(D);
                }
            }
        } else if (u.kind == 3) {
            u.mu = Rcpp::as<double>(ul["mu"]);
            u.sigma = Rcpp::as<double>(ul["sigma"]);
            u.lj = Rcpp::as<double>(ul["lj"]);
        } else if (u.kind == 10) {
            u.block = Rcpp::as<int>(ul["block"]);
            Rcpp::NumericVector mu = ul["mu"];
            u.muv.assign(mu.begin(), mu.end());
            Rcpp::NumericMatrix L = ul["L"], Si = ul["sinv"];
            u.L = Dense(L.nrow(), L.ncol());
            std::copy(L.begin(), L.end(), u.L.v.begin());
            u.Sinv = Dense(Si.nrow(), Si.ncol());
            std::copy(Si.begin(), Si.end(), u.Sinv.v.begin());
            u.logdet = Rcpp::as<double>(ul["logdet"]);
            u.tnu = Rcpp::as<double>(ul["nu"]);
            u.mc0 = Rcpp::as<double>(ul["c0"]);
            u.mhalf = Rcpp::as<double>(ul["half"]);
            M.moving = true;
        } else if (u.kind == 9) {
            u.mu = Rcpp::as<double>(ul["mu"]);
            u.sigma = Rcpp::as<double>(ul["sigma"]);
            u.tnu = Rcpp::as<double>(ul["nu"]);
            M.moving = true;
            if (ul.containsElementNamed("tfree")) {
                Rcpp::IntegerVector tf = ul["tfree"];
                u.tfree.assign(tf.begin(), tf.end());
                u.bq = Rcpp::as<bool>(ul["bq"]);
            }
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
            if (ul.containsElementNamed("chart")) {
                Rcpp::List cl = ul["chart"];
                u.chart = Rcpp::as<int>(cl["code"]);
                Rcpp::NumericVector ce = cl["eta"], chh = cl["h"];
                u.ceta.assign(ce.begin(), ce.end());
                u.ch.assign(chh.begin(), chh.end());
                if (u.chart == 9) {
                    auto mats = [&](const char* nm, std::vector<Dense>& v) {
                        Rcpp::List l = cl[nm];
                        v.clear();
                        for (int h = 0; h < l.size(); ++h) {
                            Rcpp::NumericMatrix A(l[h]);
                            Dense D(A.nrow(), A.ncol());
                            std::copy(A.begin(), A.end(), D.v.begin());
                            v.push_back(D);
                        }
                    };
                    mats("ai", u.gai);
                    mats("ai2", u.gai2);
                    Rcpp::NumericVector d1 = cl["dld"], d2 = cl["d2ld"];
                    u.gdld.assign(d1.begin(), d1.end());
                    u.gd2ld.assign(d2.begin(), d2.end());
                }
            }
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
// the outer description of a kind 10 unit (engine_outer_mvt())
void read_mvt_outer(PenUnit& u, Rcpp::List ol) {
    u.mo = true;
    u.mH = Rcpp::as<int>(ol["H"]);
    u.bq = Rcpp::as<bool>(ol["bq"]);
    auto rd = [&](const char* nm, std::vector<double>& v) {
        Rcpp::NumericVector x = ol[nm];
        v.assign(x.begin(), x.end());
    };
    rd("gt", u.mgt); rd("rho2", u.mrho2); rd("c", u.mc); rd("c2", u.mc2);
    rd("dS", u.mdS); rd("d2S", u.md2S); rd("Na", u.mNa);
    u.mai.clear();
    Rcpp::List ai = ol["ai"];
    for (int h = 0; h < ai.size(); ++h) {
        Rcpp::NumericMatrix A(ai[h]);
        Dense D(A.nrow(), A.ncol());
        std::copy(A.begin(), A.end(), D.v.begin());
        u.mai.push_back(D);
    }
}
}  // namespace

// The iwls_fit() loop on the engine's problem: the whole vector, or the
// positions E.sub with the others held at E.base (fit_smooth() on a subset).
// False where the R route signals an error.
struct IwlsOut {
    std::vector<double> beta;
    double value = 0.0, score = R_PosInf, dimless = NA_REAL;
    bool converged = false;
    int it = 0, note_kind = 0, note_it = 0, switched = 0;
    int fb_indef = 0, fb_search = 0, fb_error = 0;
    std::vector<int> aliased, h_it, h_rank, h_held;
    std::vector<double> h_obj, h_score, h_step, h_damp;
    std::vector<std::string> h_route;
};

bool iwls_loop(Engine& E, const std::vector<double>& start, const std::vector<char>& froz,
               std::vector<std::vector<int>> grp, bool verbose, IwlsOut& o) {
    Model& M = E.M;
    Method& me = E.me;
    (void) M;
    const int p = E.np();
    const double nn = (double) M.n;
    if (grp.empty()) { std::vector<int> v(p); for (int j = 0; j < p; ++j) v[j] = j; grp.push_back(v); }

    // pieces_at, backup_at, switch_at as the kind of information each reads
    bool pa_expected = me.kind_expected;
    bool has_backup = me.backup, backup_expected = true;
    const bool has_switch = me.switching;
    double switch_after = has_switch ? me.switch_after : R_PosInf;

    std::vector<double> beta(start);
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
        if (!E.fn_s(beta, value)) throw RError();
        if (verbose) Rprintf("  %-5s %14s %12s %10s\n", "iter", "objective", "score/n", "step");
        typedef Engine::Pieces Pieces;
        for (it = 1; it <= me.maxit; ++it) {
            std::vector<double> g;
            if (!E.gr_s(beta, g)) throw RError();
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
                if (!E.pieces_s(beta, pa_expected, pc_met)) throw RError();
                have_pc_met = true;
                double dec;
                if (has_backup && !E.definite(pc_met, tol_def)) {
                    Pieces pcd;
                    if (!E.pieces_s(beta, backup_expected, pcd)) throw RError();
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
                    if (!E.pieces_s(beta, false, pp_own)) throw RError();
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
                        if (!E.fn_s(cand, vnew)) vnew = R_PosInf;
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
                if (!E.pieces_s(beta, pa_expected, pc_own)) throw RError();
                pcp = &pc_own;
            }
            if (has_backup && !E.definite(*pcp, tol_def)) fb_indef += 1;
            Engine::LineSearch st;
            E.line_search(beta, value, g, *pcp, damp, froz, me.guard || has_backup, st);
            fb_error += st.errors;
            if (!st.ok && has_backup) {
                if (!E.pieces_s(beta, backup_expected, pc_own)) throw RError();
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
            if (!E.gr_s(beta, gfin)) throw RError();
            for (int j = 0; j < p; ++j) if (froz[j]) gfin[j] = 0.0;
            bool fin = true;
            for (double v : gfin) if (!R_FINITE(v)) fin = false;
            if (fin) {
                score = Engine::max_abs(gfin) / nn;
                converged = Engine::max_abs_div(gfin, nn) < me.tol;
                if (!converged) {
                    Pieces pf;
                    if (!E.pieces_s(beta, pa_expected, pf)) throw RError();
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
        return false;
    }
    o.beta = beta; o.value = value; o.score = score; o.dimless = dimless;
    o.converged = converged; o.it = it; o.note_kind = note_kind; o.note_it = note_it;
    o.switched = switched; o.fb_indef = fb_indef; o.fb_search = fb_search; o.fb_error = fb_error;
    o.aliased = aliased; o.h_it = h_it; o.h_rank = h_rank; o.h_held = h_held;
    o.h_obj = h_obj; o.h_score = h_score; o.h_step = h_step; o.h_damp = h_damp;
    o.h_route = h_route;
    return true;
}

// the settled iwls() as engine_route() describes it
void read_method(Engine& E, Rcpp::List method) {
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

}

// iwls_fit()'s fields from the loop's result
Rcpp::List loop_result(const IwlsOut& o, const Engine& E) {
    const std::vector<double>& beta = o.beta;
    const double value = o.value, score = o.score, dimless = o.dimless;
    const bool converged = o.converged;
    const int it = o.it, note_kind = o.note_kind, note_it = o.note_it, switched = o.switched;
    const int fb_indef = o.fb_indef, fb_search = o.fb_search, fb_error = o.fb_error;
    const std::vector<int>& aliased = o.aliased;
    const auto& h_it = o.h_it; const auto& h_rank = o.h_rank; const auto& h_held = o.h_held;
    const auto& h_obj = o.h_obj; const auto& h_score = o.h_score; const auto& h_step = o.h_step;
    const auto& h_damp = o.h_damp; const auto& h_route = o.h_route;

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
    read_method(E, method);
    const int p = M.p;
    std::vector<char> froz(p, 0);
    for (int f : frozen) froz[f - 1] = 1;
    std::vector<std::vector<int>> grp;
    for (int k = 0; k < groups.size(); ++k) {
        Rcpp::IntegerVector gk = groups[k];
        std::vector<int> v;
        for (int j : gk) v.push_back(j - 1);
        grp.push_back(v);
    }

    std::vector<double> start_v(start.begin(), start.end());
    IwlsOut o;
    if (!iwls_loop(E, start_v, froz, grp, verbose, o))
        return Rcpp::List::create(Rcpp::Named("status") = 1);
    return loop_result(o, E);
}

// fit_smooth() on a subset (phase 3b): the loop over the positions `sub`
// (1-based, ascending) with the others held at `start`, on the assembled
// route, as iwls_fit() runs on subset_pieces(). `frozen` and `groups` are in
// the subset's numbering, and so are the returned `par` and `aliased`.
// Status 3 where the compiled subset fit does not apply (a sparse design).
// [[Rcpp::export]]
Rcpp::List engine_subset_fit(Rcpp::List model, Rcpp::List pen, Rcpp::List method,
                             Rcpp::NumericVector start, Rcpp::IntegerVector sub,
                             Rcpp::IntegerVector frozen, Rcpp::List groups,
                             bool verbose) {
    for (int k = 0; k < T_N; ++k) g_time[k] = 0.0;
    Engine E;
    Model& M = E.M;
    engine_setup(E, model, pen);
    read_method(E, method);
    if (M.sparse) return Rcpp::List::create(Rcpp::Named("status") = 3);
    E.base.assign(start.begin(), start.end());
    for (int j : sub) E.sub.push_back(j - 1);
    const int ps = (int) E.sub.size();
    std::vector<char> froz(ps, 0);
    for (int f : frozen) froz[f - 1] = 1;
    std::vector<std::vector<int>> grp;
    for (int k = 0; k < groups.size(); ++k) {
        Rcpp::IntegerVector gk = groups[k];
        std::vector<int> v;
        for (int j : gk) v.push_back(j - 1);
        grp.push_back(v);
    }
    std::vector<double> start_v(ps);
    for (int j = 0; j < ps; ++j) start_v[j] = E.base[E.sub[j]];
    IwlsOut o;
    if (!iwls_loop(E, start_v, froz, grp, verbose, o))
        return Rcpp::List::create(Rcpp::Named("status") = 1);
    return loop_result(o, E);
}

#include "engine_kink.h"

// The compiled alternation of a model with a kinked penalty (phase 3b):
// statmod_alternate() on the covered route. `kblocks`: one list per kinked
// block (engine_alternate_blocks()). `smooth`, `held_where` 1-based;
// `smooth_groups`, `smooth_frozen` in the smooth block's numbering. Returns
// status 0 with the point, the history and the curvatures, 1 where the R
// route signals an error, 2 where it takes a branch the engine does not
// reproduce (the proximal fallback of sparse_fit()).
// [[Rcpp::export]]
Rcpp::List engine_alternate(Rcpp::List model, Rcpp::List pen, Rcpp::List method,
                            Rcpp::NumericVector start, Rcpp::List kblocks,
                            Rcpp::IntegerVector smooth, Rcpp::List smooth_groups,
                            Rcpp::IntegerVector smooth_frozen,
                            Rcpp::IntegerVector held_where,
                            Rcpp::NumericVector held_value,
                            Rcpp::LogicalVector kinked, Rcpp::IntegerVector pinned,
                            bool expected, int maxit, double tol) {
    for (int k = 0; k < T_N; ++k) g_time[k] = 0.0;
    Engine E;
    Model& M = E.M;
    engine_setup(E, model, pen);
    read_method(E, method);
    if (M.sparse) return Rcpp::List::create(Rcpp::Named("status") = 3);
    Kink Kn(E);
    Kn.expected = expected;
    for (int j : smooth) Kn.smooth.push_back(j - 1);
    for (int k = 0; k < smooth_groups.size(); ++k) {
        Rcpp::IntegerVector gk = smooth_groups[k];
        std::vector<int> v;
        for (int j : gk) v.push_back(j - 1);
        Kn.smooth_groups.push_back(v);
    }
    Kn.smooth_froz.assign(Kn.smooth.size(), 0);
    for (int f : smooth_frozen) Kn.smooth_froz[f - 1] = 1;
    for (int j : held_where) Kn.held_where.push_back(j - 1);
    Kn.held_value.assign(held_value.begin(), held_value.end());
    Kn.kinked.assign(M.p, 0);
    for (int j = 0; j < M.p && j < kinked.size(); ++j) Kn.kinked[j] = kinked[j] == TRUE;
    for (int j : pinned) Kn.pinned.push_back(j - 1);
    bool curv_scaled = false;
    for (int b = 0; b < kblocks.size(); ++b) {
        Rcpp::List bl = kblocks[b];
        KBlock B;
        B.unit = Rcpp::as<int>(bl["unit"]) - 1;
        B.param = Rcpp::as<int>(bl["param"]) - 1;
        Rcpp::IntegerVector cols = bl["cols"], oth = bl["other"];
        for (int j : cols) B.cols.push_back(j - 1);
        for (int j : oth) B.other.push_back(j - 1);
        B.ctr = Rcpp::as<bool>(bl["ctr"]);
        B.i0 = B.ctr ? Rcpp::as<int>(bl["i0"]) - 1 : -1;
        B.curv_todo = Rcpp::as<bool>(bl["curv_todo"]);
        B.curv_centered = Rcpp::as<bool>(bl["curv_centered"]);
        B.has_curv = Rcpp::as<bool>(bl["has_curv"]);
        SEXP pk = bl["prev_kink"];
        B.has_prev = !Rf_isNull(pk);
        if (B.has_prev) B.prev_kink = Rcpp::as<double>(pk);
        B.label = Rcpp::as<std::string>(bl["label"]);
        B.X = Rcpp::NumericMatrix(Rcpp::as<Rcpp::NumericMatrix>(bl["X"]));
        SEXP xo = bl["Xo"];
        if (!Rf_isNull(xo)) {
            Rcpp::NumericMatrix Xo(xo);
            B.Xo.assign(Xo.begin(), Xo.end());
        }
        if (B.has_curv) curv_scaled = true;
        Kn.blocks.push_back(B);
    }

    std::vector<double> beta(start.begin(), start.end());
    double value = 0.0;
    bool converged = false;
    std::vector<int> hb_pass, hb_it, hb_conv;
    std::vector<std::string> hb_block;
    std::vector<double> hb_obj, hb_change;
    std::vector<int> hi_pass, hi_it, hi_rank, hi_held;
    std::vector<double> hi_obj, hi_score, hi_step, hi_damp;
    std::vector<std::string> hi_route;
    int note_kind = 0, note_it = 0;
    double note_dimless = NA_REAL, note_score = NA_REAL;
    bool any_smooth = false;
    std::set<int> aliased;
    int passes = 0;
    try {
        value = Kn.fn(beta);
        for (int pass = 1; pass <= maxit; ++pass) {
            passes = pass;
            if (curv_scaled) {
                Kn.curvatures(beta);
                value = Kn.fn(beta);
            }
            const double before = value;
            if (!Kn.smooth.empty()) {
                for (std::size_t h = 0; h < Kn.held_where.size(); ++h)
                    beta[Kn.held_where[h]] = Kn.held_value[h];
                E.sub = Kn.smooth;
                E.base = beta;
                std::vector<double> st(E.sub.size());
                for (std::size_t j = 0; j < E.sub.size(); ++j) st[j] = beta[E.sub[j]];
                IwlsOut o;
                const bool ok = iwls_loop(E, st, Kn.smooth_froz, Kn.smooth_groups, false, o);
                E.sub.clear();
                if (!ok) throw RError();
                for (std::size_t j = 0; j < Kn.smooth.size(); ++j) beta[Kn.smooth[j]] = o.beta[j];
                value = Kn.fn(beta);
                any_smooth = true;
                note_kind = o.note_kind; note_it = o.note_it;
                note_dimless = o.dimless; note_score = o.score;
                for (int a : o.aliased) aliased.insert(Kn.smooth[a]);
                hb_pass.push_back(pass); hb_block.push_back("smooth"); hb_obj.push_back(value);
                hb_change.push_back(before - value); hb_it.push_back(o.it);
                hb_conv.push_back(o.converged);
                for (std::size_t r = 0; r < o.h_it.size(); ++r) {
                    hi_pass.push_back(pass); hi_it.push_back(o.h_it[r]);
                    hi_obj.push_back(o.h_obj[r]); hi_score.push_back(o.h_score[r]);
                    hi_step.push_back(o.h_step[r]); hi_rank.push_back(o.h_rank[r]);
                    hi_route.push_back(o.h_route[r]); hi_held.push_back(o.h_held[r]);
                    hi_damp.push_back(o.h_damp[r]);
                }
            }
            for (KBlock& B : Kn.blocks) {
                const double v0 = value;
                Kink::CoordOut co;
                if (!Kn.coord_fit(beta, B, co)) {
                    return Rcpp::List::create(Rcpp::Named("status") = 2);
                }
                value = co.value;
                hb_pass.push_back(pass); hb_block.push_back(B.label); hb_obj.push_back(value);
                hb_change.push_back(v0 - value); hb_it.push_back(co.iterations);
                hb_conv.push_back(co.converged);
            }
            const double rel = std::fabs(before - value) / std::max(1.0, std::fabs(value));
            if (rel < tol) {
                Kn.aliased.assign(aliased.begin(), aliased.end());
                converged = Kn.settled(beta);
                break;
            }
        }
    } catch (RError&) {
        return Rcpp::List::create(Rcpp::Named("status") = 1);
    }
    Rcpp::List cvs(Kn.blocks.size());
    for (std::size_t b = 0; b < Kn.blocks.size(); ++b)
        cvs[b] = Rcpp::wrap(M.units[Kn.blocks[b].unit].cv);
    Rcpp::IntegerVector al(aliased.size());
    {
        int k = 0;
        for (int a : aliased) al[k++] = a + 1;
    }
    return Rcpp::List::create(
        Rcpp::Named("status") = 0, Rcpp::Named("par") = Rcpp::wrap(beta),
        Rcpp::Named("value") = value, Rcpp::Named("converged") = converged,
        Rcpp::Named("passes") = passes, Rcpp::Named("aliased") = al,
        Rcpp::Named("note") = Rcpp::List::create(
            Rcpp::Named("any") = any_smooth, Rcpp::Named("kind") = note_kind,
            Rcpp::Named("it") = note_it, Rcpp::Named("dimless") = note_dimless,
            Rcpp::Named("score") = note_score),
        Rcpp::Named("blocks") = Rcpp::List::create(
            Rcpp::Named("pass") = Rcpp::wrap(hb_pass), Rcpp::Named("block") = Rcpp::wrap(hb_block),
            Rcpp::Named("objective") = Rcpp::wrap(hb_obj), Rcpp::Named("change") = Rcpp::wrap(hb_change),
            Rcpp::Named("iterations") = Rcpp::wrap(hb_it),
            Rcpp::Named("converged") = Rcpp::LogicalVector(hb_conv.begin(), hb_conv.end())),
        Rcpp::Named("inner") = Rcpp::List::create(
            Rcpp::Named("pass") = Rcpp::wrap(hi_pass), Rcpp::Named("iteration") = Rcpp::wrap(hi_it),
            Rcpp::Named("objective") = Rcpp::wrap(hi_obj), Rcpp::Named("score") = Rcpp::wrap(hi_score),
            Rcpp::Named("step") = Rcpp::wrap(hi_step), Rcpp::Named("rank") = Rcpp::wrap(hi_rank),
            Rcpp::Named("route") = Rcpp::wrap(hi_route), Rcpp::Named("held") = Rcpp::wrap(hi_held),
            Rcpp::Named("damp") = Rcpp::wrap(hi_damp)),
        Rcpp::Named("curvature") = cvs,
        Rcpp::Named("timing") = Rcpp::NumericVector::create(
            Rcpp::Named("point") = g_time[T_POINT], Rcpp::Named("parts") = g_time[T_PARTS],
            Rcpp::Named("pieces") = g_time[T_PIECES], Rcpp::Named("assembled") = g_time[T_ASM],
            Rcpp::Named("fn") = g_time[T_FN], Rcpp::Named("coord") = g_time[T_SPQR],
            Rcpp::Named("settled") = g_time[T_DQR]));
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

#include "engine_outer.h"

// The outer context at one point (phase 3): `model` and `pen` as for the
// inner fit, `beta` the mode, `hold` the pinned coordinates (1-based),
// `basis` the integrated basis of ml() or NULL, `order` 2 where the Hessian
// may be asked. Returns list(status, ctx, why): status 0 with the context,
// 2 where the point goes to the R route.
// [[Rcpp::export]]
Rcpp::List engine_outer_new(Rcpp::List model, Rcpp::List pen, Rcpp::NumericVector beta,
                            Rcpp::IntegerVector hold, SEXP basis, int order,
                            SEXP hold_mode = R_NilValue) {
    Rcpp::XPtr<Outer> X(new Outer(), true);
    Outer& O = *X;
    O.model = model;
    O.pen = pen;
    try {
        engine_setup(O.E, model, pen);
        for (int h : hold) O.hold.push_back(h - 1);
        if (!Rf_isNull(hold_mode)) {
            O.two = true;
            Rcpp::IntegerVector hm(hold_mode);
            for (int h : hm) O.hold_mode.push_back(h - 1);
        }
        if (!Rf_isNull(basis)) {
            Rcpp::NumericMatrix Bm(basis);
            O.ml = true;
            O.B = Dense(Bm.nrow(), Bm.ncol());
            std::copy(Bm.begin(), Bm.end(), O.B.v.begin());
        }
        std::vector<double> b(beta.begin(), beta.end());
        O.setup(b, order >= 2);
    } catch (RError&) {
        O.give_back("the R route signals an error");
    }
    return Rcpp::List::create(Rcpp::Named("status") = O.status, Rcpp::Named("ctx") = X,
                              Rcpp::Named("why") = O.why);
}

// [[Rcpp::export]]
Rcpp::List engine_outer_value(SEXP ctx) {
    Rcpp::XPtr<Outer> X(ctx);
    Outer& O = *X;
    if (O.status) return Rcpp::List::create(Rcpp::Named("status") = O.status);
    const double v = O.ll - O.rho + O.q / 2.0 * std::log(2.0 * M_PI) - O.logdet_c / 2.0;
    return Rcpp::List::create(Rcpp::Named("status") = 0, Rcpp::Named("value") = v,
                              Rcpp::Named("loglik") = O.ll, Rcpp::Named("penalty") = O.rho,
                              Rcpp::Named("logdet") = O.logdet_c, Rcpp::Named("q") = O.q);
}

// members: `mrow`, `munit` 1-based; `mchain` dlinkinv at each member's
// value; `gam` the estimated coefficients' stacked positions (1-based)
// [[Rcpp::export]]
Rcpp::List engine_outer_grad(SEXP ctx, Rcpp::IntegerVector mrow, Rcpp::IntegerVector munit,
                             Rcpp::IntegerVector mname, Rcpp::NumericVector mchain, int nh,
                             Rcpp::IntegerVector gam, bool free) {
    Rcpp::XPtr<Outer> X(ctx);
    Outer& O = *X;
    if (O.status) return Rcpp::List::create(Rcpp::Named("status") = O.status);
    std::vector<int> r, u, g, nmv;
    for (int v : mrow) r.push_back(v - 1);
    for (int v : munit) u.push_back(v - 1);
    for (int v : mname) nmv.push_back(v - 1);
    for (int v : gam) g.push_back(v - 1);
    std::vector<double> ch(mchain.begin(), mchain.end()), out;
    if (!O.gradient(r, u, nmv, ch, nh, g, free, out))
        return Rcpp::List::create(Rcpp::Named("status") = 2);
    return Rcpp::List::create(Rcpp::Named("status") = 0,
                              Rcpp::Named("gradient") = Rcpp::wrap(out));
}

// [[Rcpp::export]]
Rcpp::List engine_outer_hess(SEXP ctx, Rcpp::IntegerVector mrow, Rcpp::IntegerVector munit,
                             Rcpp::IntegerVector mname, Rcpp::NumericVector mchain, int nh,
                             Rcpp::IntegerVector gam, Rcpp::NumericVector h1,
                             Rcpp::NumericVector h2) {
    Rcpp::XPtr<Outer> X(ctx);
    Outer& O = *X;
    if (O.status) return Rcpp::List::create(Rcpp::Named("status") = O.status);
    std::vector<int> r, u, g, nmv;
    for (int v : mrow) r.push_back(v - 1);
    for (int v : munit) u.push_back(v - 1);
    for (int v : mname) nmv.push_back(v - 1);
    for (int v : gam) g.push_back(v - 1);
    std::vector<double> ch(mchain.begin(), mchain.end()), a1(h1.begin(), h1.end()),
        a2(h2.begin(), h2.end());
    Dense out;
    if (!O.hessian(r, u, nmv, ch, nh, g, a1, a2, out))
        return Rcpp::List::create(Rcpp::Named("status") = 2);
    Rcpp::NumericMatrix R(out.nr, out.nc);
    std::copy(out.v.begin(), out.v.end(), R.begin());
    return Rcpp::List::create(Rcpp::Named("status") = 0, Rcpp::Named("hessian") = R);
}

// The kind 10 units of an outer context whose derivatives in the
// hyperparameters have not been given yet (1-based), which the criterion does
// not read and the gradient and the Hessian do (engine_outer_complete())
// [[Rcpp::export]]
Rcpp::IntegerVector engine_outer_mvt_need(SEXP ctx) {
    Rcpp::XPtr<Outer> X(ctx);
    Outer& O = *X;
    Rcpp::IntegerVector out;
    for (std::size_t k = 0; k < O.E.M.units.size(); ++k) {
        const PenUnit& u = O.E.M.units[k];
        if (u.kind == 10 && !u.mo) out.push_back((int) k + 1);
    }
    return out;
}

// [[Rcpp::export]]
void engine_outer_mvt_set(SEXP ctx, int k, Rcpp::List ob) {
    Rcpp::XPtr<Outer> X(ctx);
    Outer& O = *X;
    read_mvt_outer(O.E.M.units[k - 1], ob);
    O.have_uder = false;
}

// K^-1 x with K the mode's pinned matrix of an outer context, the inverse
// that inner_mode_error() and criterion_resolution() read (engine_ctx_solve())
// [[Rcpp::export]]
Rcpp::List engine_outer_kinv(SEXP ctx, Rcpp::NumericVector x) {
    Rcpp::XPtr<Outer> X(ctx);
    Outer& O = *X;
    if (O.status || x.size() != O.p || !O.ensure_inverse())
        return Rcpp::List::create(Rcpp::Named("status") = 2);
    std::vector<double> y(O.p);
    dense_matvec(O.Kinv, x.begin(), y.data());
    return Rcpp::List::create(Rcpp::Named("status") = 0, Rcpp::Named("x") = Rcpp::wrap(y));
}
