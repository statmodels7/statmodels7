// The outer criterion in the engine (piano_motore.md, phase 3, stage A):
// reml() and ml() on the observed information, the value, the exact
// gradient and the exact Hessian in the outer coordinates, at one point
// (beta_hat, theta), as statmod_marginal(), statmod_marginal_grad() and
// statmod_marginal_hess() compute them on the covered route
// (misure_motore_fase3/SPEC_outer.md, sections 4-8).
//
// One context per outer point (D24): at the mode the derivatives of orders
// 1-4 on the link scale, H, the pinned K with its factor and inverse, the
// trace matrix and the leverages are computed once, and the criterion, the
// gradient and the Hessian read them.
//
// The context hands a point back to the R route (status 2) wherever the R
// route takes a branch that it does not reproduce: a non-finite curvature or
// a value outside the domain. The factorization follows every route of R's
// pd_factor() (D37, 2026-11-07): Cholesky accepted at R's own threshold, the
// eigenvalues otherwise, and the deficient-coordinate retry; where R's
// criterion or inverse is NULL the entry points answer status 4. Near the
// threshold the two routes may decide differently on rounding (and R's
// sparse route reads another condition estimate), so there the comparison
// with R is at a tolerance scaled on the condition number. Included once, by
// engine_iwls.cpp, after the Engine class.
#ifndef SM7_ENGINE_OUTER_H
#define SM7_ENGINE_OUTER_H

namespace {

typedef void (*lf7_inv1234p_t)(int, const double*, double, double*);
typedef void (*d7_parts_d34_t)(int, const double*, double*);

struct OuterRegistry {
    bool ok = false;
    lf7_inv1234p_t inv4 = nullptr;
    d7_parts_d34_t d3 = nullptr, d4 = nullptr;
};
OuterRegistry& outer_registry() {
    static OuterRegistry r;
    if (!r.ok) {
        r.inv4 = (lf7_inv1234p_t) R_GetCCallable("linkfunctions7", "lf7_inv1234p");
        r.d3 = (d7_parts_d34_t) R_GetCCallable("distributions7", "d7_parts_d3");
        r.d4 = (d7_parts_d34_t) R_GetCCallable("distributions7", "d7_parts_d4");
        r.ok = true;
    }
    return r;
}

// one row of a design block, its nonzero columns and values
struct RowNz {
    std::vector<int> start;    // n + 1
    std::vector<int> col;
    std::vector<double> val;
};

inline void dense_matmul(const Dense& A, const Dense& B, Dense& C) {
    C = Dense(A.nr, B.nc);
    if (!A.nr || !B.nc || !A.nc) return;
    const double one = 1.0, zero = 0.0;
    int m = A.nr, n = B.nc, k = A.nc;
    F77_CALL(dgemm)("N", "N", &m, &n, &k, &one, A.v.data(), &m, B.v.data(), &k,
                    &zero, C.v.data(), &m FCONE FCONE);
}

inline void dense_matvec(const Dense& A, const double* x, double* y) {
    r_matvec(A.v.data(), A.nr, A.nc, x, y);
}

}  // namespace

#include "engine_additive.h"

namespace {

class Outer {
public:
    Engine E;
    Rcpp::List model, pen;          // keep the R objects the engine points into
    OuterRegistry& oreg;
    int n = 0, K = 0, p = 0;
    std::vector<double> beta;
    std::vector<int> hold;          // pinned_coords(), 0-based
    bool ml = false;
    Dense B;                        // ml(): the integrated basis, p x nb
    int status = 0;                 // 0 usable, 2 hand back to the R route
    std::string why;

    // per observation, link scale: l1 K, l2 K*K, l3 K^3, l4 K^4
    std::vector<double> l1, l2, l3, l4;
    bool have4 = false;
    std::vector<RowNz> rows;        // per equation
    Dense H, S, Kp, Kinv, Mt;       // Mt the trace matrix (held rows zeroed)
    std::vector<int> held_at;       // pinned positions of Kp (no boundary here)
    // phase 3b: beside a kinked penalty the criterion's matrix is pinned at
    // the pinned and the kinked coordinates (laplace_pinned()), the mode's at
    // the pinned coordinates and the kinked ones at zero (zero_kinked()); the
    // mode's matrix gives Kinv, the movement of the mode and the Hessian's
    // Jmat, the criterion's the determinant and the trace matrix
    bool two = false;
    std::vector<int> hold_mode, held_mode;
    Dense Kmode, Umode;
    double logdet_c = 0;            // the criterion's log-determinant
    int q = 0;
    double ll = 0, rho = 0;
    // leverages G_ab per observation for a <= b, at index a*K + b
    std::vector<std::vector<double>> G;
    std::vector<double> u;
    bool have_u = false;

    Outer() : oreg(outer_registry()) {}

    static int idx3(int K, int a, int b, int c) { return a + K * (b + K * c); }
    static int idx4(int K, int a, int b, int c, int d) {
        return a + K * (b + K * (c + K * d));
    }

    void give_back(const char* w) { if (status == 0) { status = 2; why = w; } }

    // ---------------------------------------------------- the point
    void setup(const std::vector<double>& b, bool need4) {
        Model& M = E.M;
        n = M.n; K = M.K; p = M.p;
        beta = b;
        Point& P = E.cur;
        E.point_at(beta, P);
        if (!P.inb) { give_back("outside the domain"); return; }
        // where linkinv()'s clamp acted the R route reads a curvature that is
        // not finite and pins the coordinate as a boundary (pin_boundary());
        // the context leaves such a point to it
        if (P.clamped) { give_back("a parameter at its link's clamp"); return; }
        // row-wise nonzeros of each design block
        rows.assign(K, RowNz());
        for (int k = 0; k < K; ++k) {
            const Block& Bk = M.X[k];
            RowNz& R = rows[k];
            R.start.assign(n + 1, 0);
            if (Bk.m == 0) continue;
            if (!Bk.sp) {
                for (int i = 0; i < n; ++i) R.start[i + 1] = R.start[i] + Bk.m;
                R.col.resize((std::size_t) n * Bk.m);
                R.val.resize((std::size_t) n * Bk.m);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < Bk.m; ++j) {
                        R.col[(std::size_t) i * Bk.m + j] = j;
                        R.val[(std::size_t) i * Bk.m + j] = Bk.d[(std::size_t) j * n + i];
                    }
            } else {
                std::vector<int> cnt(n, 0);
                for (int j = 0; j < Bk.m; ++j)
                    for (SpMat::InnerIterator it(Bk.s, j); it; ++it) cnt[it.row()]++;
                for (int i = 0; i < n; ++i) R.start[i + 1] = R.start[i] + cnt[i];
                R.col.resize(R.start[n]);
                R.val.resize(R.start[n]);
                std::vector<int> pos(R.start.begin(), R.start.end() - 1);
                for (int j = 0; j < Bk.m; ++j)
                    for (SpMat::InnerIterator it(Bk.s, j); it; ++it) {
                        const int i = it.row();
                        R.col[pos[i]] = j;
                        R.val[pos[i]] = it.value();
                        pos[i]++;
                    }
            }
        }
        (void) need4;
        derivs(P, 2);
        if (status) return;
        // the log-likelihood and the penalty, as the objective computes them
        {
            std::vector<double> lv(n);
            for (int i = 0; i < n; ++i)
                lv[i] = M.w[i] * E.reg.logpdf(M.fam, M.y[i], E.theta_obs(P, i));
            ll = r_sum(lv.data(), n);
            rho = E.pen_value(beta);
        }
        if (!R_FINITE(ll) || !R_FINITE(rho)) { give_back("non-finite value"); return; }
        information();
        if (status) return;
        factor();
    }

    // the derivatives on the link scale to order `level` (2, 3 or 4): the
    // criterion reads orders 1-2, the gradient 3, the Hessian 4, so the
    // higher orders are computed only when they are asked for
    int have_level = 0;
    bool ensure_level(int level) {
        if (status) return false;
        if (have_level < level) derivs(E.cur, level);
        return status == 0;
    }
    void derivs(Point& P, int level) {
        Model& M = E.M;
        const int K2 = K * K, K3 = K2 * K, K4 = K3 * K;
        const bool need3 = level >= 3, need4 = level >= 4;
        l1.assign((std::size_t) n * K, 0.0);
        l2.assign((std::size_t) n * K2, 0.0);
        if (need3) l3.assign((std::size_t) n * K3, 0.0);
        if (need4) l4.assign((std::size_t) n * K4, 0.0);
        have4 = need4;
        have_level = level;
        std::vector<double> buf(M.pmax), g(K), Hh(K2), T(need3 ? K3 : 1), Q(need4 ? K4 : 1);
        std::vector<double> h1(K), h2(K), h3(K), h4(K);
        for (int i = 0; i < n; ++i) {
            const double* th = E.theta_obs(P, i);
            E.reg.parts(M.fam, M.y[i], th, level, buf.data());
            E.reg.grad(M.fam, buf.data(), g.data());
            E.reg.hess(M.fam, buf.data(), Hh.data());
            if (need3) oreg.d3(M.fam, buf.data(), T.data());
            if (need4) oreg.d4(M.fam, buf.data(), Q.data());
            for (int k = 0; k < K && need3; ++k) {
                const double* par = M.lpar[k].empty() ? nullptr : M.lpar[k].data();
                const double e2 = E.reg.fun(M.lid[k], par, P.th[(std::size_t) k * n + i]);
                double h[5];
                oreg.inv4(M.lid[k], par, e2, h);
                h1[k] = P.h1[(std::size_t) k * n + i];
                h2[k] = P.h2[(std::size_t) k * n + i];
                h3[k] = h[3];
                h4[k] = h[4];
            }
            double* o1 = l1.data() + (std::size_t) i * K;
            double* o2 = l2.data() + (std::size_t) i * K2;
            if (!need3)
                for (int k = 0; k < K; ++k) {
                    h1[k] = P.h1[(std::size_t) k * n + i];
                    h2[k] = P.h2[(std::size_t) k * n + i];
                }
            for (int a = 0; a < K; ++a) o1[a] = g[a] * h1[a];
            for (int a = 0; a < K; ++a)
                for (int b = 0; b < K; ++b)
                    o2[a + K * b] = Hh[a + K * b] * h1[a] * h1[b] +
                                    (a == b ? g[a] * h2[a] : 0.0);
            if (!need3) continue;
            double* o3 = l3.data() + (std::size_t) i * K3;
            for (int a = 0; a < K; ++a)
                for (int b = 0; b < K; ++b)
                    for (int c = 0; c < K; ++c) {
                        double v = T[idx3(K, a, b, c)] * h1[a] * h1[b] * h1[c];
                        if (a == b) v += Hh[a + K * c] * h2[a] * h1[c];
                        if (a == c) v += Hh[a + K * b] * h2[a] * h1[b];
                        if (b == c) v += Hh[a + K * b] * h1[a] * h2[b];
                        if (a == b && b == c) v += g[a] * h3[a];
                        o3[idx3(K, a, b, c)] = v;
                    }
            if (need4) {
                double* o4 = l4.data() + (std::size_t) i * K4;
                for (int a = 0; a < K; ++a)
                    for (int b = 0; b < K; ++b)
                        for (int c = 0; c < K; ++c)
                            for (int d = 0; d < K; ++d) {
                                const double t = T[idx3(K, a, b, c)];
                                double v = Q[idx4(K, a, b, c, d)] * h1[a] * h1[b] * h1[c] * h1[d];
                                if (a == d) v += t * h2[a] * h1[b] * h1[c];
                                if (b == d) v += t * h1[a] * h2[b] * h1[c];
                                if (c == d) v += t * h1[a] * h1[b] * h2[c];
                                if (a == b) {
                                    v += T[idx3(K, a, c, d)] * h1[d] * h2[a] * h1[c];
                                    if (a == d) v += Hh[a + K * c] * h3[a] * h1[c];
                                    if (c == d) v += Hh[a + K * c] * h2[a] * h2[c];
                                }
                                if (a == c) {
                                    v += T[idx3(K, a, b, d)] * h1[d] * h2[a] * h1[b];
                                    if (a == d) v += Hh[a + K * b] * h3[a] * h1[b];
                                    if (b == d) v += Hh[a + K * b] * h2[a] * h2[b];
                                }
                                if (b == c) {
                                    v += T[idx3(K, a, b, d)] * h1[d] * h1[a] * h2[b];
                                    if (a == d) v += Hh[a + K * b] * h2[a] * h2[b];
                                    if (b == d) v += Hh[a + K * b] * h1[a] * h3[b];
                                }
                                if (a == b && b == c) {
                                    v += Hh[a + K * d] * h1[d] * h3[a];
                                    if (a == d) v += g[a] * h4[a];
                                }
                                o4[idx4(K, a, b, c, d)] = v;
                            }
            }
        }
        for (double v : l2) if (!R_FINITE(v)) { give_back("non-finite curvature"); return; }
        if (need3) for (double v : l3) if (!R_FINITE(v)) { give_back("non-finite third derivative"); return; }
        if (need4) for (double v : l4) if (!R_FINITE(v)) { give_back("non-finite fourth derivative"); return; }
    }

    // H = -sum_i w_i l2_ab x_ia x_ib', dense
    void information() {
        Model& M = E.M;
        H = Dense(p, p);
        for (int a = 0; a < K; ++a)
            for (int b = a; b < K; ++b) {
                if (!M.npar[a] || !M.npar[b]) continue;
                const RowNz& Ra = rows[a];
                const RowNz& Rb = rows[b];
                const int oa = M.offs[a], ob = M.offs[b];
                for (int i = 0; i < n; ++i) {
                    const double wv = -M.w[i] * l2[(std::size_t) i * K * K + a + K * b];
                    if (wv == 0.0) continue;
                    for (int s = Ra.start[i]; s < Ra.start[i + 1]; ++s) {
                        const double xa = wv * Ra.val[s];
                        const int ja = oa + Ra.col[s];
                        for (int t = Rb.start[i]; t < Rb.start[i + 1]; ++t)
                            H(ja, ob + Rb.col[t]) += xa * Rb.val[t];
                    }
                }
            }
        // the lower blocks and the off-diagonal halves from the upper ones
        for (int a = 0; a < K; ++a)
            for (int b = a + 1; b < K; ++b) {
                if (!M.npar[a] || !M.npar[b]) continue;
                for (int j = 0; j < M.npar[a]; ++j)
                    for (int k = 0; k < M.npar[b]; ++k)
                        H(M.offs[b] + k, M.offs[a] + j) = H(M.offs[a] + j, M.offs[b] + k);
            }
        S = Dense(p, p);
        if (M.sparse) {
            for (int j = 0; j < M.Ss.outerSize(); ++j)
                for (SpMat::InnerIterator it(M.Ss, j); it; ++it) S(it.row(), j) = it.value();
        } else {
            S = M.Sd;
        }
        for (double v : S.v) if (!R_FINITE(v)) { give_back("non-finite penalty"); return; }
    }

    // a Cholesky factor accepted clearly by the test of pd_factor() on the
    // free block, or false (hand back)
    bool clear_chol(const Dense& A, const std::vector<int>& held, Dense& U) {
        if (!r_chol(A, U)) return false;
        return chol_test(A, held, U, 100.0);
    }

    // the test of pd_logdet_dense() on the free block, lmin = rcond * ||A_ff||_1
    // against sqrt(eps) ||A_ff||_1, with a margin (1 is R's own test)
    bool chol_test(const Dense& A, const std::vector<int>& held, const Dense& U,
                   double margin) {
        const int m = A.nr;
        std::vector<char> isheld(m, 0);
        for (int j : held) isheld[j] = 1;
        std::vector<int> fr;
        for (int j = 0; j < m; ++j) if (!isheld[j]) fr.push_back(j);
        if (fr.empty()) return true;
        const int f = (int) fr.size();
        Dense Uf(f, f);
        double af = 0.0;
        for (int jj = 0; jj < f; ++jj) {
            ldouble s = 0.0;
            for (int ii = 0; ii < f; ++ii) {
                Uf(ii, jj) = U(fr[ii], fr[jj]);
                s += std::fabs(A(fr[ii], fr[jj]));
            }
            af = std::max(af, (double) s);
        }
        std::vector<double> work(3 * (std::size_t) f);
        std::vector<int> iwork(f);
        double rc = 0.0;
        int info = 0, ff = f;
        F77_CALL(dpocon)("U", &ff, Uf.v.data(), &ff, &af, &rc, work.data(),
                         iwork.data(), &info FCONE);
        if (info != 0 || !R_FINITE(rc)) return false;
        const double lmin = rc * af, thr = std::sqrt(DBL_EPSILON) * af;
        return lmin > margin * thr;
    }

    // pd_factor() / pd_logdet_dense() of a pinned dense matrix (D37): the
    // Cholesky factor accepted at R's own threshold, otherwise the
    // eigenvalues, accepted where the smallest exceeds p eps times the
    // largest; the factor is kept wherever chol() succeeds, since
    // ctx_penalized() inverts through it on either route
    struct Pd { bool ok = false; double logdet = 0; bool have_U = false; Dense U; };
    Pd pd_dense(const Dense& A, const std::vector<int>& held) {
        Pd o;
        for (double v : A.v) if (!R_FINITE(v)) return o;
        if (r_chol(A, o.U)) {
            o.have_U = true;
            if (chol_test(A, held, o.U, 1.0)) {
                o.ok = true;
                o.logdet = logdet_of(o.U);
                return o;
            }
        }
        std::vector<double> ev;
        Dense vv;
        if (!r_eigen_sym(A, true, ev, vv)) return o;
        if (ev.empty()) return o;
        double mn = ev[0], mx = ev[0];
        for (double e : ev) {
            if (!R_FINITE(e)) return o;
            mn = std::min(mn, e);
            mx = std::max(mx, e);
        }
        if (mn > (double) A.nr * DBL_EPSILON * mx) {
            std::vector<double> lg(ev.size());
            for (std::size_t j = 0; j < ev.size(); ++j) lg[j] = std::log(ev[j]);
            o.ok = true;
            o.logdet = r_sum(lg.data(), (int) lg.size());
        }
        return o;
    }

    // the matrix H + S pinned at `held` (sorted, unique)
    Dense pinned_sum(const std::vector<int>& held) {
        Dense A(p, p);
        for (std::size_t t = 0; t < A.v.size(); ++t) A.v[t] = H.v[t] + S.v[t];
        pin(A, held);
        return A;
    }

    // pd_factor() with the deficient_coords() retry of ctx_penalized() and
    // statmod_marginal(): where the first factorization fails, the coordinates
    // the matrix does not identify are pinned too and the matrix factorized
    // again
    Pd pd_retry(Dense& A, std::vector<int>& held, bool retry) {
        Pd r = pd_dense(A, held);
        if (r.ok || !retry) return r;
        std::vector<int> fl = Kink::deficient(A, held);
        if (fl.empty()) return r;
        held.insert(held.end(), fl.begin(), fl.end());
        std::sort(held.begin(), held.end());
        held.erase(std::unique(held.begin(), held.end()), held.end());
        A = pinned_sum(held);
        return pd_dense(A, held);
    }

    static double logdet_of(const Dense& U) {
        std::vector<double> lg(U.nr);
        for (int j = 0; j < U.nr; ++j) lg[j] = std::log(U(j, j));
        return 2.0 * r_sum(lg.data(), U.nr);
    }

    static Dense chol2inv(const Dense& U) {
        Dense V = U;
        int m = U.nr, info = 0;
        if (m) F77_CALL(dpotri)("U", &m, V.v.data(), &m, &info FCONE);
        for (int j = 0; j < m; ++j)
            for (int i = j + 1; i < m; ++i) V(i, j) = V(j, i);
        return V;
    }

    void pin(Dense& A, const std::vector<int>& j) {
        for (int c : j) {
            for (int i = 0; i < A.nr; ++i) { A(i, c) = 0.0; A(c, i) = 0.0; }
            A(c, c) = 1.0;
        }
    }

    // D37: every route of R's factorization is taken here; where R's criterion
    // is NULL the context records it (crit_na) and where R's inverse is NULL
    // it records that (inv_ok false), and the entry points answer status 4
    bool crit_na = false, inv_ok = true;
    void factor() {
        held_at = hold;
        std::sort(held_at.begin(), held_at.end());
        held_at.erase(std::unique(held_at.begin(), held_at.end()), held_at.end());
        // a coordinate whose curvature is not finite is pinned by R as a
        // boundary (pin_boundary(), boundary_coords()), which the criterion's
        // gradient reads as a constant row; that branch stays on the R route
        {
            std::vector<char> isheld(p, 0);
            for (int j : held_at) isheld[j] = 1;
            for (int j = 0; j < p; ++j)
                if (!isheld[j] && !R_FINITE(H(j, j) + S(j, j))) {
                    give_back("boundary coordinate");
                    return;
                }
        }
        Kp = pinned_sum(held_at);
        // the criterion's matrix: statmod_marginal() retries only without a
        // kinked penalty and without an integrated basis; ctx_penalized()
        // retries the same matrix where nothing is kinked
        Pd pc = pd_retry(Kp, held_at, !two && !ml);
        if (two) {
            if (ml) { give_back("ml beside a kinked penalty"); return; }
            // the mode's matrix, ctx_penalized(laplace = FALSE): with the retry
            held_mode = hold_mode;
            std::sort(held_mode.begin(), held_mode.end());
            held_mode.erase(std::unique(held_mode.begin(), held_mode.end()), held_mode.end());
            Kmode = pinned_sum(held_mode);
            Pd pm = pd_retry(Kmode, held_mode, true);
            if (!pm.ok || !pm.have_U) inv_ok = false; else Umode = pm.U;
            // the trace matrix inverts the criterion's matrix
            if (!pc.have_U) inv_ok = false; else Ufac = pc.U;
        } else if (ml) {
            // the mode's matrix (the inverse) is K pinned at hold, with the retry
            Kmode = Kp;
            held_mode = held_at;
            Pd pm = pd_retry(Kmode, held_mode, true);
            if (!pm.ok || !pm.have_U) inv_ok = false; else Umode = pm.U;
            Ufac = Umode;
        } else {
            Kmode = Kp;
            held_mode = held_at;
            if (!pc.ok || !pc.have_U) inv_ok = false; else { Umode = pc.U; Ufac = pc.U; }
        }
        if (!ml) {
            if (!pc.ok) { crit_na = true; return; }
            logdet_c = pc.logdet;
            q = p - (int) held_at.size();
        } else {
            const int nb = B.nc;
            q = nb;
            // the criterion's matrix, B'(H + S)B, pinned only where its diagonal
            // is not finite (the context gives such points back before here)
            Dense Kc(p, p);
            for (std::size_t t = 0; t < Kc.v.size(); ++t) Kc.v[t] = H.v[t] + S.v[t];
            Dense KB, BKB, Bt = transpose(B);
            if (nb == 0) { logdet_c = 0; return; }
            dense_matmul(Kc, B, KB);
            dense_matmul(Bt, KB, BKB);
            std::vector<int> none;
            Pd pb = pd_dense(BKB, none);
            if (!pb.ok) { crit_na = true; return; }
            logdet_c = pb.logdet;
        }
    }

    // the inverse of the mode's matrix and the matrix of the traces, which
    // only the gradient and the Hessian read
    Dense Ufac;
    bool have_inv = false;
    bool ensure_inverse() {
        if (status || !inv_ok) return false;
        if (have_inv) return true;
        Kinv = chol2inv(Umode);
        if (!ml) {
            Mt = two ? chol2inv(Ufac) : Kinv;
        } else if (B.nc == 0) {
            Mt = Dense(p, p);
        } else {
            Dense Bt = transpose(B);
            // the trace matrix, B (B'K_pB)^-1 B' with held rows zeroed
            Dense KpB, BKpB;
            dense_matmul(Kp, B, KpB);
            dense_matmul(Bt, KpB, BKpB);
            Dense Up;
            if (!r_chol(BKpB, Up)) { give_back("ml trace factor"); return false; }
            Dense inner = chol2inv(Up), Bi;
            dense_matmul(B, inner, Bi);
            dense_matmul(Bi, Bt, Mt);
        }
        for (int c : held_at)
            for (int i = 0; i < p; ++i) { Mt(i, c) = 0.0; Mt(c, i) = 0.0; }
        have_inv = true;
        return true;
    }

    // ---------------------------------------------------- leverages
    void leverage() {
        if (!G.empty()) return;
        Model& M = E.M;
        G.assign(K * K, std::vector<double>());
        for (int a = 0; a < K; ++a)
            for (int b = a; b < K; ++b) {
                if (!M.npar[a] || !M.npar[b]) continue;
                std::vector<double>& g = G[a * K + b];
                g.assign(n, 0.0);
                const RowNz& Ra = rows[a];
                const RowNz& Rb = rows[b];
                const int oa = M.offs[a], ob = M.offs[b];
                for (int i = 0; i < n; ++i) {
                    double s = 0.0;
                    for (int t = Ra.start[i]; t < Ra.start[i + 1]; ++t) {
                        double r = 0.0;
                        const int ja = oa + Ra.col[t];
                        for (int v = Rb.start[i]; v < Rb.start[i + 1]; ++v)
                            r += Mt(ja, ob + Rb.col[v]) * Rb.val[v];
                        s += Ra.val[t] * r;
                    }
                    g[i] = s;
                }
            }
    }
    const std::vector<double>& Gab(int a, int b) const {
        return a <= b ? G[a * K + b] : G[b * K + a];
    }

    // X_k' v for the stacked equation k
    void xt(int k, const double* v, double* out) {
        Model& M = E.M;
        const RowNz& R = rows[k];
        const int m = M.npar[k];
        std::vector<ldouble> acc(m, 0.0);
        for (int i = 0; i < n; ++i)
            for (int t = R.start[i]; t < R.start[i + 1]; ++t)
                acc[R.col[t]] += (ldouble) R.val[t] * v[i];
        for (int j = 0; j < m; ++j) out[j] = (double) acc[j];
    }
    // X_k b_k
    void xv(int k, const double* b, std::vector<double>& out) {
        const RowNz& R = rows[k];
        out.assign(n, 0.0);
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int t = R.start[i]; t < R.start[i + 1]; ++t) s += R.val[t] * b[R.col[t]];
            out[i] = s;
        }
    }

    // u_c = tr(M dK/dbeta_c) = -sum_i w_i x_ic sum_{a,b} l3_abk G_ab
    void u_vector() {
        if (have_u) return;
        leverage();
        Model& M = E.M;
        u.assign(p, 0.0);
        const int K3 = K * K * K;
        std::vector<double> s(n), w(n);
        for (int k = 0; k < K; ++k) {
            if (!M.npar[k]) continue;
            std::fill(s.begin(), s.end(), 0.0);
            for (int a = 0; a < K; ++a)
                for (int b = 0; b < K; ++b) {
                    if (!M.npar[a] || !M.npar[b]) continue;
                    const std::vector<double>& g = Gab(a, b);
                    for (int i = 0; i < n; ++i)
                        s[i] += l3[(std::size_t) i * K3 + idx3(K, a, b, k)] * g[i];
                }
            for (int i = 0; i < n; ++i) w[i] = M.w[i] * s[i];
            std::vector<double> o(M.npar[k]);
            xt(k, w.data(), o.data());
            for (int j = 0; j < M.npar[k]; ++j) u[M.offs[k] + j] = -o[j];
        }
        have_u = true;
    }

    // ---------------------------------------------------- penalty pieces
    // A unit's derivatives in its H hyperparameters, on the parameter scale:
    // the values (gt, rho2), the cross derivatives (c, c2, over the unit's m
    // positions) and the derivatives of its Hessian (dS, d2S), which are a
    // blk x blk block repeated down the diagonal (blk = m for a quadratic
    // penalty, 1 for a gaussian prior, the dimension for a multivariate one)
    struct UnitDer {
        int H = 0, m = 0, blk = 0;
        std::vector<double> gt;
        std::vector<std::vector<double>> c;       // H x m
        std::vector<Dense> dS;                    // H
        Dense rho2;                               // H x H
        std::vector<std::vector<double>> c2;      // H*H x m
        std::vector<Dense> d2S;                   // H*H
        // a prior whose Hessian moves (kind 9, phase 3b): dS and d2S are
        // diagonal with an entry per coefficient (dSv, d2Sv), and the
        // movement in beta is read from three coefficients per position:
        // dS/dbeta[x] = diag(a3 x), d2S/dbeta2[v, w] = diag(a4 v w),
        // d(dS/dtheta_h)/dbeta[v] = diag(a3h[h] v); `moving` is false where
        // beta_quadratic() holds at the point
        bool vec = false, moving = false;
        std::vector<std::vector<double>> dSv, d2Sv;
        std::vector<double> a3, a4;
        std::vector<std::vector<double>> a3h;
        // a multivariate t prior (kind 10): dS and d2S block by block
        // (dSg[h][g], d2Sg[h * H + k][g]), and per group the pieces that the
        // movement in beta reads: from mvt_yscaled() x = ws and is, from
        // mvt_dpieces() w, s, and per hyperparameter sd and wd
        bool grp = false;
        int G = 0;
        std::vector<std::vector<Dense>> dSg, d2Sg;
        std::vector<std::vector<double>> gx, gw, gsd, gwd;   // G x p; H x G; H x G*p
        std::vector<double> gis, gs;
    };
    std::vector<UnitDer> uder;
    std::vector<char> uder_ok;
    bool have_uder = false;

    // the chart of a multivariate unit: Sigma (or Omega) and its first and
    // second derivatives in the free values
    struct Chart {
        int p = 0, K = 0;
        Dense V;                                  // Sigma, or Omega if inverted
        std::vector<Dense> A;                     // K
        std::vector<Dense> A2;                    // K*K, (t, s) at t * K + s
    };
    bool chart_build(const PenUnit& un, Chart& ch) {
        const int p = un.block, K = (int) un.ceta.size();
        ch.p = p; ch.K = K;
        ch.V = Dense(p, p);
        ch.A.assign(K, Dense(p, p));
        ch.A2.assign((std::size_t) K * K, Dense(p, p));
        if (un.chart == 1) {
            // log-Cholesky: the diagonal positions first, then the strict
            // lower triangle column by column
            std::vector<int> row(K), col(K);
            int k = 0;
            for (int i = 0; i < p && k < K; ++i) { row[k] = i; col[k] = i; ++k; }
            for (int j = 0; j < p; ++j)
                for (int i = j + 1; i < p && k < K; ++i) { row[k] = i; col[k] = j; ++k; }
            if (k != K || K != p * (p + 1) / 2) return false;
            Dense L(p, p);
            std::vector<double> v(K);
            for (int t = 0; t < K; ++t) {
                const double e = un.ceta[t];
                L(row[t], col[t]) = row[t] == col[t] ? std::exp(e) : e;
                v[t] = row[t] == col[t] ? L(row[t], col[t]) : 1.0;
            }
            for (int i = 0; i < p; ++i)
                for (int j = 0; j < p; ++j) {
                    double s = 0.0;
                    for (int a = 0; a < p; ++a) s += L(i, a) * L(j, a);
                    ch.V(i, j) = s;
                }
            // dV/de_t = E_t L' + L E_t', E_t holding v_t at (row, col)
            for (int t = 0; t < K; ++t) {
                const int r = row[t], c = col[t];
                for (int j = 0; j < p; ++j) {
                    ch.A[t](r, j) += v[t] * L(j, c);
                    ch.A[t](j, r) += v[t] * L(j, c);
                }
            }
            for (int t = 0; t < K; ++t)
                for (int s = 0; s < K; ++s) {
                    Dense& D = ch.A2[(std::size_t) t * K + s];
                    // E_t E_s' + E_s E_t'
                    if (col[t] == col[s]) {
                        D(row[t], row[s]) += v[t] * v[s];
                        D(row[s], row[t]) += v[t] * v[s];
                    }
                    // the second derivative of exp on a diagonal entry
                    if (t == s && row[t] == col[t]) {
                        const int r = row[t];
                        for (int j = 0; j < p; ++j) {
                            D(r, j) += v[t] * L(j, r);
                            D(j, r) += v[t] * L(j, r);
                        }
                    }
                }
            return true;
        }
        if (un.chart == 2 || un.chart == 3) {
            // compound symmetry (2) and AR(1) (3): V = s P(r), K = 2, with the
            // scale's and the correlation's inverse links and their first two
            // derivatives at the free values
            if (K != 2 || un.ch.size() != 6) return false;
            const double s0 = un.ch[0], s1 = un.ch[1], s2 = un.ch[2];
            const double r0 = un.ch[3], r1 = un.ch[4], r2 = un.ch[5];
            Dense P0(p, p), P1(p, p), P2(p, p);
            for (int i = 0; i < p; ++i)
                for (int j = 0; j < p; ++j) {
                    const int lag = std::abs(i - j);
                    if (un.chart == 2) {
                        P0(i, j) = i == j ? 1.0 : r0;
                        P1(i, j) = i == j ? 0.0 : r1;
                        P2(i, j) = i == j ? 0.0 : r2;
                    } else {
                        P0(i, j) = lag == 0 ? 1.0 : R_pow_di(r0, lag);
                        const double d1 = lag == 0 ? 0.0 : lag * R_pow_di(r0, lag - 1);
                        const double d2 = lag < 2 ? 0.0 : lag * (lag - 1.0) * R_pow_di(r0, lag - 2);
                        P1(i, j) = d1 * r1;
                        P2(i, j) = d2 * r1 * r1 + d1 * r2;
                    }
                }
            for (int t = 0; t < p * p; ++t) {
                ch.V.v[t] = s0 * P0.v[t];
                ch.A[0].v[t] = s1 * P0.v[t];
                ch.A[1].v[t] = s0 * P1.v[t];
                ch.A2[0].v[t] = s2 * P0.v[t];
                ch.A2[3].v[t] = s0 * P2.v[t];
                ch.A2[1].v[t] = s1 * P1.v[t];
                ch.A2[2].v[t] = s1 * P1.v[t];
            }
            return true;
        }
        return false;
    }
    // sum(A * t(B)), the trace of A B
    static double tr_prod(const Dense& A, const Dense& B) {
        std::vector<double> e((std::size_t) A.nr * A.nc);
        for (int j = 0; j < A.nc; ++j)
            for (int i = 0; i < A.nr; ++i) e[(std::size_t) j * A.nr + i] = A(i, j) * B(j, i);
        return r_sum(e.data(), e.size());
    }

    // a unit's derivatives on the coefficients the penalty reads: for a
    // gaussian prior with a map (kinds 3 and 4) the mapped ones t = A b,
    // carried back to b by map_der() as map_back() and map_quad_full() carry
    // the penalty's own derivatives
    bool unit_der(int k, UnitDer& out) {
        const PenUnit& un = E.M.units[k];
        std::vector<double> bu, t;
        E.unit_t(un, beta, bu, t);
        const bool mapped = un.has_map && (un.kind == 3 || un.kind == 4);
        if (!mapped) return unit_der_core(un, bu, out);
        PenUnit u2 = un;
        u2.has_map = false;
        if (!unit_der_core(u2, t, out)) return false;
        return map_der(un, out);
    }
    bool map_der(const PenUnit& un, UnitDer& out) {
        if (out.vec || out.grp) return false;
        Dense A;
        if (!un.map.sp) A = un.map.d;
        else {
            A = Dense(un.map.nr, un.map.nc);
            for (int j = 0; j < un.map.s.outerSize(); ++j)
                for (SpMat::InnerIterator it(un.map.s, j); it; ++it) A(it.row(), j) = it.value();
        }
        const int mt = A.nr, kk = A.nc;
        if (mt != out.m || kk != (int) un.idx.size()) return false;
        Dense At = transpose(A);
        auto back = [&](std::vector<double>& c) {
            std::vector<double> r(kk);
            dense_matvec(At, c.data(), r.data());
            c = r;
        };
        auto quad = [&](Dense& X) {
            Dense D(mt, mt), DA, R;
            for (int g = 0; g * out.blk < mt; ++g)
                for (int b = 0; b < out.blk; ++b)
                    for (int a = 0; a < out.blk; ++a)
                        D(g * out.blk + a, g * out.blk + b) = X(a, b);
            dense_matmul(D, A, DA);
            dense_matmul(At, DA, R);
            X = R;
        };
        for (auto& c : out.c) back(c);
        for (auto& c : out.c2) back(c);
        for (auto& X : out.dS) quad(X);
        for (auto& X : out.d2S) quad(X);
        out.blk = kk;
        out.m = kk;
        return true;
    }
    bool unit_der_core(const PenUnit& un, const std::vector<double>& bu, UnitDer& out) {
        const int m = (int) bu.size();
        out.m = m;
        if (un.kind == 1) {
            if (un.Mg.nr != m || un.rank < 0) return false;
            out.H = 1; out.blk = m;
            std::vector<double> sb(m), qv(m);
            mat_vec(un.Mg, bu.data(), sb.data());
            for (int j = 0; j < m; ++j) qv[j] = bu[j] * sb[j];
            const double qf = r_sum(qv.data(), m), lam = un.lam, r = un.rank;
            out.gt = {qf / 2.0 - r / (2.0 * lam)};
            out.rho2 = Dense(1, 1);
            out.rho2(0, 0) = r / (2.0 * lam * lam);
            out.c = {sb};
            out.c2 = {std::vector<double>(m, 0.0)};
            Dense D(m, m);
            if (!un.Mg.sp) D = un.Mg.d;
            else
                for (int j = 0; j < un.Mg.s.outerSize(); ++j)
                    for (SpMat::InnerIterator it(un.Mg.s, j); it; ++it) D(it.row(), j) = it.value();
            out.dS = {D};
            out.d2S = {Dense(m, m)};
            return true;
        }
        if (un.kind == 2 && un.add_m > 0) {
            // an AdditivePenalty: S = sum_k lam_k P_k, the log
            // pseudo-determinant through additive_stable() (additive_sum())
            const int Kh = (int) un.mats.size();
            std::vector<Dense> P(Kh), Pb;
            for (int h = 0; h < Kh; ++h) {
                const MatAny& A = un.mats[h];
                if (A.nr != m) return false;
                P[h] = Dense(m, m);
                if (!A.sp) P[h] = A.d;
                else
                    for (int j = 0; j < A.s.outerSize(); ++j)
                        for (SpMat::InnerIterator it(A.s, j); it; ++it) P[h](it.row(), j) = it.value();
            }
            double val = 0.0;
            std::vector<double> dlog;
            Dense d2log;
            const bool one = un.add_m == 1;
            const std::vector<Dense>& src = one ? P : un.add_block;
            const int rk = (int) std::llround(one ? un.rank : un.rank / un.add_m);
            if (!additive_stable_c(src, un.lams, rk, un.add_tol, val, dlog, d2log)) return false;
            if (!one) {
                for (double& z : dlog) z *= un.add_m;
                for (double& z : d2log.v) z *= un.add_m;
            }
            // the rank identity additive_sum() checks (additive_check_tol() = 1e-3)
            std::vector<double> ld(Kh);
            for (int h = 0; h < Kh; ++h) ld[h] = un.lams[h] * dlog[h];
            if (!(std::fabs(r_sum(ld.data(), Kh) - un.rank) <= 1e-3 * un.rank)) return false;
            out.H = Kh; out.blk = m;
            out.gt.assign(Kh, 0.0);
            out.c.assign(Kh, std::vector<double>(m));
            out.dS = P;
            out.rho2 = Dense(Kh, Kh);
            out.c2.assign((std::size_t) Kh * Kh, std::vector<double>(m, 0.0));
            out.d2S.assign((std::size_t) Kh * Kh, Dense(m, m));
            for (int h = 0; h < Kh; ++h) {
                dense_matvec(P[h], bu.data(), out.c[h].data());
                std::vector<double> q(m);
                for (int j = 0; j < m; ++j) q[j] = bu[j] * out.c[h][j];
                out.gt[h] = 0.5 * r_sum(q.data(), m) - 0.5 * dlog[h];
                for (int s = 0; s < Kh; ++s) out.rho2(h, s) = 0.5 * d2log(h, s);
            }
            return true;
        }
        if (un.kind == 3 && !un.has_map) {
            out.H = 1; out.blk = 1;
            const double s = un.sigma, mu = un.mu;
            std::vector<double> d2(m);
            for (int j = 0; j < m; ++j) { const double d = bu[j] - mu; d2[j] = d * d; }
            const double ss = r_sum(d2.data(), m);
            const double s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
            out.gt = {m / s - ss / s3};
            out.rho2 = Dense(1, 1);
            out.rho2(0, 0) = 3.0 * ss / s4 - m / s2;
            out.c.assign(1, std::vector<double>(m));
            out.c2.assign(1, std::vector<double>(m));
            for (int j = 0; j < m; ++j) {
                const double d = bu[j] - mu;
                out.c[0][j] = -2.0 * d / s3;
                out.c2[0][j] = 6.0 * d / s4;
            }
            Dense a(1, 1), b(1, 1);
            a(0, 0) = -2.0 / s3;
            b(0, 0) = 6.0 / s4;
            out.dS = {a};
            out.d2S = {b};
            return true;
        }
        if (un.kind == 9 && !un.tfree.empty()) {
            // the t prior: every derivative from the registry at level 4 in
            // (mu, sigma, nu); a derivative in the response is one in mu with
            // the sign of y - mu
            const int Ht = (int) un.tfree.size();
            const int tid = E.reg.id("StudentT1Distrib");
            if (tid < 0) return false;
            out.H = Ht; out.blk = 1; out.vec = true; out.moving = !un.bq;
            std::vector<double> buf(E.M.pmax), g(3), Hh(9), T(27), Q(81);
            std::vector<std::vector<double>> gv(Ht, std::vector<double>(m));
            std::vector<std::vector<double>> hv((std::size_t) Ht * Ht, std::vector<double>(m));
            out.c.assign(Ht, std::vector<double>(m));
            out.dSv.assign(Ht, std::vector<double>(m));
            out.c2.assign((std::size_t) Ht * Ht, std::vector<double>(m));
            out.d2Sv.assign((std::size_t) Ht * Ht, std::vector<double>(m));
            out.a3.assign(m, 0.0);
            out.a4.assign(m, 0.0);
            out.a3h.assign(Ht, std::vector<double>(m));
            const double th[3] = {un.mu, un.sigma, un.tnu};
            for (int j = 0; j < m; ++j) {
                E.reg.parts(tid, bu[j], th, 4, buf.data());
                E.reg.grad(tid, buf.data(), g.data());
                E.reg.hess(tid, buf.data(), Hh.data());
                oreg.d3(tid, buf.data(), T.data());
                oreg.d4(tid, buf.data(), Q.data());
                for (int h = 0; h < Ht; ++h) {
                    const int a = un.tfree[h];
                    gv[h][j] = g[a];
                    out.c[h][j] = Hh[0 + 3 * a];
                    out.dSv[h][j] = -T[idx3(3, 0, 0, a)];
                    out.a3h[h][j] = Q[idx4(3, 0, 0, 0, a)];
                    for (int k = 0; k < Ht; ++k) {
                        const int b = un.tfree[k];
                        const std::size_t key = (std::size_t) h * Ht + k;
                        hv[key][j] = Hh[a + 3 * b];
                        out.c2[key][j] = T[idx3(3, 0, a, b)];
                        out.d2Sv[key][j] = -Q[idx4(3, 0, 0, a, b)];
                    }
                }
                out.a3[j] = T[0];
                out.a4[j] = -Q[0];
            }
            out.gt.assign(Ht, 0.0);
            out.rho2 = Dense(Ht, Ht);
            for (int h = 0; h < Ht; ++h) {
                out.gt[h] = -r_sum(gv[h].data(), m);
                for (int k = 0; k < Ht; ++k)
                    out.rho2(h, k) = -r_sum(hv[(std::size_t) h * Ht + k].data(), m);
            }
            out.dS.assign(Ht, Dense(1, 1));
            out.d2S.assign((std::size_t) Ht * Ht, Dense(1, 1));
            return true;
        }
        if (un.kind == 10 && un.mo) {
            const int pp = un.block, Hh = un.mH;
            if (pp <= 0 || m % pp) return false;
            const int G = m / pp;
            out.H = Hh; out.blk = pp; out.grp = true; out.G = G; out.moving = !un.bq;
            out.gt = un.mgt;
            out.rho2 = Dense(Hh, Hh);
            std::copy(un.mrho2.begin(), un.mrho2.end(), out.rho2.v.begin());
            out.c.assign(Hh, std::vector<double>(m));
            for (int h = 0; h < Hh; ++h)
                for (int j = 0; j < m; ++j) out.c[h][j] = un.mc[(std::size_t) h * m + j];
            out.c2.assign((std::size_t) Hh * Hh, std::vector<double>(m));
            for (int key = 0; key < Hh * Hh; ++key)
                for (int j = 0; j < m; ++j) out.c2[key][j] = un.mc2[(std::size_t) key * m + j];
            const std::size_t pq = (std::size_t) pp * pp;
            auto blocks = [&](const std::vector<double>& src, int nk,
                              std::vector<std::vector<Dense>>& dst) {
                dst.assign(nk, std::vector<Dense>(G, Dense(pp, pp)));
                for (int k = 0; k < nk; ++k)
                    for (int g = 0; g < G; ++g)
                        std::copy(src.begin() + ((std::size_t) k * G + g) * pq,
                                  src.begin() + ((std::size_t) k * G + g + 1) * pq,
                                  dst[k][g].v.begin());
            };
            blocks(un.mdS, Hh, out.dSg);
            blocks(un.md2S, Hh * Hh, out.d2Sg);
            out.dS.assign(Hh, Dense(1, 1));
            out.d2S.assign((std::size_t) Hh * Hh, Dense(1, 1));
            // the pieces of the movement, group by group
            out.gx.assign(G, std::vector<double>(pp));
            out.gw.assign(G, std::vector<double>(pp));
            out.gis.assign(G, 0.0);
            out.gs.assign(G, 0.0);
            out.gsd.assign(Hh, std::vector<double>(G));
            out.gwd.assign(Hh, std::vector<double>((std::size_t) G * pp));
            std::vector<double> ws, r(pp), w2(pp), sq(pp);
            for (int g = 0; g < G; ++g) {
                const double* y = bu.data() + (std::size_t) g * pp;
                double is, cw;
                Engine::mvt_yscaled(un, y, ws, is, cw);
                out.gx[g] = ws;
                out.gis[g] = is;
                for (int j = 0; j < pp; ++j) r[j] = y[j] - un.muv[j];
                for (int c = 0; c < pp; ++c) {
                    double w = 0.0;
                    for (int l = 0; l < pp; ++l) w += r[l] * un.Sinv(l, c);
                    out.gw[g][c] = w;
                }
                w2 = r;
                if (!r_trisolve(un.L, false, w2)) return false;
                for (int j = 0; j < pp; ++j) sq[j] = w2[j] * w2[j];
                out.gs[g] = un.tnu + r_colsum(sq.data(), pp);
                for (int h = 0; h < Hh; ++h) {
                    const Dense& A = un.mai[h];
                    double sd = 0.0;
                    for (int c = 0; c < pp; ++c) {
                        double wd = 0.0;
                        for (int l = 0; l < pp; ++l) wd += r[l] * A(l, c);
                        out.gwd[h][(std::size_t) g * pp + c] = wd;
                        sd += wd * r[c];
                    }
                    out.gsd[h][g] = un.mNa[h] != 0.0 ? 1.0 : sd;
                }
            }
            return true;
        }
        if (un.kind == 4 && !un.has_map && un.chart == 9) {
            // a chart that the context does not write, its derivatives from R
            const int pp = un.block, KK = (int) un.gai.size();
            if (pp <= 0 || m % pp || (int) un.gai2.size() != KK * KK ||
                (int) un.gdld.size() != KK || (int) un.gd2ld.size() != KK * KK) return false;
            const int G = m / pp;
            std::vector<Dense> ai = un.gai, ai2 = un.gai2;
            std::vector<double> dld = un.gdld;
            Dense d2ld(KK, KK);
            for (int t = 0; t < KK; ++t)
                for (int s2 = 0; s2 < KK; ++s2) d2ld(t, s2) = un.gd2ld[(std::size_t) t * KK + s2];
            // the groups' deviations from the mean, group by group
            std::vector<double> R(m);
            for (int g = 0; g < G; ++g)
                for (int c = 0; c < pp; ++c)
                    R[(std::size_t) g * pp + c] = bu[(std::size_t) g * pp + c] -
                        (un.muv.size() == (std::size_t) pp ? un.muv[c] : 0.0);
            // W = R A by group, and sum_g R_g' A R_g
            auto quad = [&](const Dense& A, std::vector<double>& W) {
                W.assign(m, 0.0);
                std::vector<double> q(G);
                for (int g = 0; g < G; ++g) {
                    std::vector<double> e(pp);
                    for (int c = 0; c < pp; ++c) {
                        double w = 0.0;
                        for (int a = 0; a < pp; ++a) w += R[(std::size_t) g * pp + a] * A(a, c);
                        W[(std::size_t) g * pp + c] = w;
                        e[c] = w * R[(std::size_t) g * pp + c];
                    }
                    q[g] = r_sum(e.data(), pp);
                }
                return r_sum(q.data(), G);
            };
            out.H = KK; out.blk = pp;
            out.gt.assign(KK, 0.0);
            out.c.assign(KK, std::vector<double>());
            out.dS = ai;
            out.rho2 = Dense(KK, KK);
            out.c2.assign((std::size_t) KK * KK, std::vector<double>());
            out.d2S = ai2;
            for (int t = 0; t < KK; ++t)
                out.gt[t] = 0.5 * G * dld[t] + 0.5 * quad(ai[t], out.c[t]);
            for (int t = 0; t < KK; ++t)
                for (int s = 0; s < KK; ++s)
                    out.rho2(t, s) = 0.5 * G * d2ld(t, s) +
                        0.5 * quad(ai2[(std::size_t) t * KK + s], out.c2[(std::size_t) t * KK + s]);
            return true;
        }
        if (un.kind == 4 && !un.has_map && un.chart > 0 && un.chart != 9) {
            Chart ch;
            if (!chart_build(un, ch)) return false;
            const int pp = ch.p, KK = ch.K;
            if (pp <= 0 || m % pp) return false;
            const int G = m / pp;
            Dense U;
            if (!r_chol(ch.V, U)) return false;
            Dense Vi = chol2inv(U);
            // ai, ai2: derivatives of Sigma^-1; dld, d2ld: of log|Sigma|
            std::vector<Dense> ai(KK), ai2((std::size_t) KK * KK), ViA(KK);
            std::vector<double> dld(KK);
            Dense d2ld(KK, KK);
            for (int t = 0; t < KK; ++t) dense_matmul(Vi, ch.A[t], ViA[t]);
            if (un.chart == 1 && !un.inverted) {
                // log-Cholesky of Sigma: through Gm = L^-1, as the R route
                // takes it (chol_inv_pieces(), SPEC_blockpen 4.1), which
                // keeps the digits that inverting L L' loses where a
                // diagonal of L goes to zero
                std::vector<int> row(KK), col(KK);
                int k = 0;
                for (int i = 0; i < pp; ++i) { row[k] = i; col[k] = i; ++k; }
                for (int j = 0; j < pp; ++j)
                    for (int i = j + 1; i < pp; ++i) { row[k] = i; col[k] = j; ++k; }
                Dense L(pp, pp);
                std::vector<double> v(KK);
                for (int t = 0; t < KK; ++t) {
                    L(row[t], col[t]) = row[t] == col[t] ? std::exp(un.ceta[t]) : un.ceta[t];
                    v[t] = row[t] == col[t] ? L(row[t], col[t]) : 1.0;
                }
                Dense Gm(pp, pp);
                for (int j = 0; j < pp; ++j) {
                    std::vector<double> e(pp, 0.0);
                    e[j] = 1.0;
                    r_trisolve(L, false, e);
                    for (int i = 0; i < pp; ++i) Gm(i, j) = e[i];
                }
                Dense Gt = transpose(Gm);
                std::vector<Dense> Bk(KK), Ck(KK);
                for (int t = 0; t < KK; ++t) {
                    Bk[t] = Dense(pp, pp);
                    for (int i = 0; i < pp; ++i) Bk[t](i, col[t]) = v[t] * Gm(i, row[t]);
                    Ck[t] = Dense(pp, pp);
                    for (int j = 0; j < pp; ++j)
                        for (int i = 0; i < pp; ++i) Ck[t](i, j) = Bk[t](i, j) + Bk[t](j, i);
                    Dense X;
                    dense_matmul(Gt, Ck[t], X);
                    dense_matmul(X, Gm, ai[t]);
                    for (double& z : ai[t].v) z = -z;
                    dld[t] = row[t] == col[t] ? 2.0 : 0.0;
                }
                for (int t = 0; t < KK; ++t)
                    for (int s = 0; s < KK; ++s) {
                        // dB = (-B_s) B_t (+ Gm D2), X = (B_s' C_t + C_t B_s) - (dB + dB')
                        Dense mB = Bk[s], dB;
                        for (double& z : mB.v) z = -z;
                        dense_matmul(mB, Bk[t], dB);
                        if (t == s && row[t] == col[t])
                            for (int i = 0; i < pp; ++i)
                                dB(i, col[t]) += Gm(i, row[t]) * L(row[t], col[t]);
                        Dense BsT = transpose(Bk[s]), X1, X2, X(pp, pp), Y;
                        dense_matmul(BsT, Ck[t], X1);
                        dense_matmul(Ck[t], Bk[s], X2);
                        for (int j = 0; j < pp; ++j)
                            for (int i = 0; i < pp; ++i)
                                X(i, j) = (X1(i, j) + X2(i, j)) - (dB(i, j) + dB(j, i));
                        dense_matmul(Gt, X, Y);
                        dense_matmul(Y, Gm, ai2[(std::size_t) t * KK + s]);
                        d2ld(t, s) = 0.0;
                    }
            } else {
            for (int t = 0; t < KK; ++t) {
                if (!un.inverted) {
                    dense_matmul(ViA[t], Vi, ai[t]);
                    for (double& z : ai[t].v) z = -z;
                } else {
                    ai[t] = ch.A[t];
                }
                double tr = 0.0;
                for (int i = 0; i < pp; ++i) tr += ViA[t](i, i);
                dld[t] = un.inverted ? -tr : tr;
            }
            for (int t = 0; t < KK; ++t)
                for (int s = 0; s < KK; ++s) {
                    const Dense& A2 = ch.A2[(std::size_t) t * KK + s];
                    Dense ViA2;
                    dense_matmul(Vi, A2, ViA2);
                    double tr = 0.0;
                    for (int i = 0; i < pp; ++i) tr += ViA2(i, i);
                    const double d2 = tr - tr_prod(ViA[t], ViA[s]);
                    d2ld(t, s) = un.inverted ? -d2 : d2;
                    if (un.inverted) { ai2[(std::size_t) t * KK + s] = A2; continue; }
                    // Vi (A_s Vi A_t + A_t Vi A_s - A_ts) Vi
                    Dense X1, X2, Y(pp, pp), Z;
                    dense_matmul(ch.A[s], ViA[t], X1);
                    dense_matmul(ch.A[t], ViA[s], X2);
                    for (int z = 0; z < pp * pp; ++z) Y.v[z] = X1.v[z] + X2.v[z] - A2.v[z];
                    dense_matmul(Vi, Y, Z);
                    dense_matmul(Z, Vi, ai2[(std::size_t) t * KK + s]);
                }
            }
            // the groups' deviations from the mean, group by group
            std::vector<double> R(m);
            for (int g = 0; g < G; ++g)
                for (int c = 0; c < pp; ++c)
                    R[(std::size_t) g * pp + c] = bu[(std::size_t) g * pp + c] -
                        (un.muv.size() == (std::size_t) pp ? un.muv[c] : 0.0);
            // W = R A by group, and sum_g R_g' A R_g
            auto quad = [&](const Dense& A, std::vector<double>& W) {
                W.assign(m, 0.0);
                std::vector<double> q(G);
                for (int g = 0; g < G; ++g) {
                    std::vector<double> e(pp);
                    for (int c = 0; c < pp; ++c) {
                        double w = 0.0;
                        for (int a = 0; a < pp; ++a) w += R[(std::size_t) g * pp + a] * A(a, c);
                        W[(std::size_t) g * pp + c] = w;
                        e[c] = w * R[(std::size_t) g * pp + c];
                    }
                    q[g] = r_sum(e.data(), pp);
                }
                return r_sum(q.data(), G);
            };
            out.H = KK; out.blk = pp;
            out.gt.assign(KK, 0.0);
            out.c.assign(KK, std::vector<double>());
            out.dS = ai;
            out.rho2 = Dense(KK, KK);
            out.c2.assign((std::size_t) KK * KK, std::vector<double>());
            out.d2S = ai2;
            for (int t = 0; t < KK; ++t)
                out.gt[t] = 0.5 * G * dld[t] + 0.5 * quad(ai[t], out.c[t]);
            for (int t = 0; t < KK; ++t)
                for (int s = 0; s < KK; ++s)
                    out.rho2(t, s) = 0.5 * G * d2ld(t, s) +
                        0.5 * quad(ai2[(std::size_t) t * KK + s], out.c2[(std::size_t) t * KK + s]);
            return true;
        }
        return false;
    }
    void ensure_uder() {
        if (have_uder) return;
        const int nu = (int) E.M.units.size();
        uder.assign(nu, UnitDer());
        uder_ok.assign(nu, 0);
        for (int k = 0; k < nu; ++k) uder_ok[k] = unit_der(k, uder[k]) ? 1 : 0;
        have_uder = true;
    }
    // the unit's block X, repeated down the diagonal, added at its positions
    void place(const PenUnit& un, const UnitDer& ud, const Dense& X, Dense& T) {
        for (int g = 0; g * ud.blk < ud.m; ++g)
            for (int b = 0; b < ud.blk; ++b)
                for (int a = 0; a < ud.blk; ++a)
                    T(un.idx[g * ud.blk + a], un.idx[g * ud.blk + b]) += X(a, b);
    }
    // tr(M X) for the unit's block X repeated down the diagonal
    double trace_unit(const PenUnit& un, const UnitDer& ud, const Dense& X) {
        std::vector<double> e;
        e.reserve((std::size_t) ud.m * ud.blk);
        for (int g = 0; g * ud.blk < ud.m; ++g)
            for (int b = 0; b < ud.blk; ++b)
                for (int a = 0; a < ud.blk; ++a)
                    e.push_back(Mt(un.idx[g * ud.blk + a], un.idx[g * ud.blk + b]) * X(a, b));
        return r_sum(e.data(), e.size());
    }

    // the per-coefficient forms of place() and trace_unit() (kind 9)
    void place_v(const PenUnit& un, const std::vector<double>& d, Dense& T) {
        for (std::size_t j = 0; j < d.size(); ++j) T(un.idx[j], un.idx[j]) += d[j];
    }
    double trace_v(const PenUnit& un, const std::vector<double>& d) {
        std::vector<double> e(d.size());
        for (std::size_t j = 0; j < d.size(); ++j) e[j] = Mt(un.idx[j], un.idx[j]) * d[j];
        return r_sum(e.data(), e.size());
    }
    // the group-block forms of place() and trace_unit() (kind 10)
    void place_g(const PenUnit& un, const std::vector<Dense>& B, Dense& T) {
        const int pp = un.block;
        for (std::size_t g = 0; g < B.size(); ++g)
            for (int b = 0; b < pp; ++b)
                for (int a = 0; a < pp; ++a)
                    T(un.idx[g * pp + a], un.idx[g * pp + b]) += B[g](a, b);
    }
    double trace_g(const PenUnit& un, const std::vector<Dense>& B) {
        const int pp = un.block;
        std::vector<double> e;
        e.reserve(B.size() * pp * pp);
        for (std::size_t g = 0; g < B.size(); ++g)
            for (int b = 0; b < pp; ++b)
                for (int a = 0; a < pp; ++a)
                    e.push_back(Mt(un.idx[g * pp + a], un.idx[g * pp + b]) * B[g](a, b));
        return r_sum(e.data(), e.size());
    }
    // the closed forms of section 6.3 of SPEC_tprior.md, the penalty's (the
    // negated response derivatives), for group g of a kind 10 unit; v and w
    // the group's slices of the directions
    static double dotp(const double* a, const double* b, int p) {
        double s = 0.0;
        for (int i = 0; i < p; ++i) s += a[i] * b[i];
        return s;
    }
    void mvt_dSb(const PenUnit& un, const UnitDer& ud, int g, const double* v, Dense& B) {
        const int p = un.block;
        const double N = un.tnu + p, is = ud.gis[g];
        const double* x = ud.gx[g].data();
        std::vector<double> Ov(p, 0.0);
        for (int a = 0; a < p; ++a)
            for (int c = 0; c < p; ++c) Ov[a] += un.Sinv(a, c) * v[c];
        const double xv = dotp(x, v, p);
        B = Dense(p, p);
        for (int b = 0; b < p; ++b)
            for (int a = 0; a < p; ++a)
                B(a, b) = -(2 * N * is * (un.Sinv(a, b) * xv + Ov[a] * x[b] + x[a] * Ov[b]) -
                            8 * N * xv * x[a] * x[b]);
    }
    void mvt_d2Sb(const PenUnit& un, const UnitDer& ud, int g, const double* v,
                  const double* w, Dense& B) {
        const int p = un.block;
        const double N = un.tnu + p, is = ud.gis[g];
        const double* x = ud.gx[g].data();
        std::vector<double> Ov(p, 0.0), Ow(p, 0.0);
        for (int a = 0; a < p; ++a)
            for (int c = 0; c < p; ++c) { Ov[a] += un.Sinv(a, c) * v[c]; Ow[a] += un.Sinv(a, c) * w[c]; }
        const double xv = dotp(x, v, p), xw = dotp(x, w, p), vOw = dotp(v, Ow.data(), p);
        B = Dense(p, p);
        for (int b = 0; b < p; ++b)
            for (int a = 0; a < p; ++a) {
                const double O = un.Sinv(a, b);
                const double t1 = O * vOw + Ov[a] * Ow[b] + Ow[a] * Ov[b];
                const double t2 = O * xv * xw + Ov[a] * x[b] * xw + Ow[a] * x[b] * xv +
                                  x[a] * Ov[b] * xw + x[a] * Ow[b] * xv + vOw * x[a] * x[b];
                B(a, b) = -(2 * N * is * is * t1 - 8 * N * is * t2 + 48 * N * xv * xw * x[a] * x[b]);
            }
    }
    void mvt_dSbt(const PenUnit& un, const UnitDer& ud, int h, int g, const double* v, Dense& B) {
        const int p = un.block;
        const double N = un.tnu + p, Na = un.mNa[h], s = ud.gs[g], sd = ud.gsd[h][g];
        const double* w = ud.gw[g].data();
        const double* wd = ud.gwd[h].data() + (std::size_t) g * p;
        const Dense& Oa = un.mai[h];
        std::vector<double> Ov(p, 0.0), Oav(p, 0.0);
        for (int a = 0; a < p; ++a)
            for (int c = 0; c < p; ++c) { Ov[a] += un.Sinv(a, c) * v[c]; Oav[a] += Oa(a, c) * v[c]; }
        const double wv = dotp(w, v, p), wdv = dotp(wd, v, p);
        const double s2 = s * s, s3 = s2 * s, s4 = s3 * s;
        const double k1 = 2 * Na / s2 - 4 * N * sd / s3, k2 = 2 * N / s2;
        const double k3 = 8 * Na / s3 - 24 * N * sd / s4, k4 = 8 * N / s3;
        B = Dense(p, p);
        for (int b = 0; b < p; ++b)
            for (int a = 0; a < p; ++a) {
                const double O = un.Sinv(a, b);
                const double T3 = O * wv + Ov[a] * w[b] + w[a] * Ov[b];
                const double dT3 = Oa(a, b) * wv + Oav[a] * w[b] + w[a] * Oav[b] +
                                   O * wdv + Ov[a] * wd[b] + wd[a] * Ov[b];
                const double ww = wv * w[a] * w[b];
                const double dww = wdv * w[a] * w[b] + wv * (wd[a] * w[b] + w[a] * wd[b]);
                B(a, b) = -(k1 * T3 + k2 * dT3 - k3 * ww - k4 * dww);
            }
    }
    // the kind 10 movements of a unit over every group
    void mvt_blocks(int which, int k, int h, const std::vector<double>& v,
                    const std::vector<double>* w, std::vector<Dense>& out) {
        const PenUnit& un = E.M.units[k];
        const UnitDer& ud = uder[k];
        const int p = un.block;
        out.assign(ud.G, Dense());
        std::vector<double> vg(p), wg(p);
        for (int g = 0; g < ud.G; ++g) {
            for (int a = 0; a < p; ++a) {
                vg[a] = v[un.idx[(std::size_t) g * p + a]];
                if (w) wg[a] = (*w)[un.idx[(std::size_t) g * p + a]];
            }
            if (which == 1) mvt_dSb(un, ud, g, vg.data(), out[g]);
            else if (which == 2) mvt_d2Sb(un, ud, g, vg.data(), wg.data(), out[g]);
            else mvt_dSbt(un, ud, h, g, vg.data(), out[g]);
        }
    }

    // trp = sum over the moving units of tr(M dS/dbeta[x])
    double trp(const std::vector<double>& x) {
        double s = 0.0;
        for (std::size_t k = 0; k < uder.size(); ++k) {
            if (!uder_ok[k] || !uder[k].moving) continue;
            const PenUnit& un = E.M.units[k];
            const UnitDer& ud = uder[k];
            if (ud.grp) {
                std::vector<Dense> B;
                mvt_blocks(1, (int) k, 0, x, nullptr, B);
                s = s + trace_g(un, B);
                continue;
            }
            std::vector<double> e(ud.m);
            for (int j = 0; j < ud.m; ++j)
                e[j] = Mt(un.idx[j], un.idx[j]) * (ud.a3[j] * x[un.idx[j]]);
            s = s + r_sum(e.data(), ud.m);
        }
        return s;
    }
    // P[x] = sum over the moving units of dS/dbeta[x], added into T
    void add_P(const std::vector<double>& x, Dense& T) {
        for (std::size_t k = 0; k < uder.size(); ++k) {
            if (!uder_ok[k] || !uder[k].moving) continue;
            const PenUnit& un = E.M.units[k];
            const UnitDer& ud = uder[k];
            if (ud.grp) {
                std::vector<Dense> B;
                mvt_blocks(1, (int) k, 0, x, nullptr, B);
                place_g(un, B, T);
                continue;
            }
            for (int j = 0; j < ud.m; ++j)
                T(un.idx[j], un.idx[j]) += ud.a3[j] * x[un.idx[j]];
        }
    }

    // ---------------------------------------------------- the gradient
    // members: row, unit, name (the hyperparameter's position in its unit),
    // chain (dlinkinv at the member's value); gam: the stacked positions the
    // criterion estimates. Parameter scale where `free` is false.
    bool gradient(const std::vector<int>& mrow, const std::vector<int>& munit,
                  const std::vector<int>& mname, const std::vector<double>& mchain,
                  int nh, const std::vector<int>& gam, bool free,
                  std::vector<double>& out) {
        if (!ensure_level(3) || !ensure_inverse()) return false;
        u_vector();
        ensure_uder();
        const int nm = (int) mrow.size();
        out.assign(nh + gam.size(), 0.0);
        for (int i = 0; i < nm; ++i) {
            const int k = munit[i], h = mname[i];
            if (!uder_ok[k] || h < 0 || h >= uder[k].H) return false;
            const UnitDer& ud = uder[k];
            const PenUnit& un = E.M.units[k];
            std::vector<double> cm(p, 0.0), v(p);
            for (int j = 0; j < ud.m; ++j) cm[un.idx[j]] = ud.c[h][j];
            dense_matvec(Kinv, cm.data(), v.data());
            for (double& x : v) x = -x;
            const double tr = ud.grp ? trace_g(un, ud.dSg[h])
                            : ud.vec ? trace_v(un, ud.dSv[h]) : trace_unit(un, ud, ud.dS[h]);
            std::vector<double> uv(p);
            for (int j = 0; j < p; ++j) uv[j] = u[j] * v[j];
            const double dt = -ud.gt[h] - ((tr + r_sum(uv.data(), p)) + trp(v)) / 2.0;
            out[mrow[i]] += free ? dt * mchain[i] : dt;
        }
        if (gam.empty()) return true;
        // the score at the estimated coefficients, and the mode's movement
        Model& M = E.M;
        std::vector<int> holdv = held_mode;
        for (std::size_t j = 0; j < gam.size(); ++j) {
            const int wj = gam[j];
            int a = -1;
            for (int e = 0; e < K; ++e)
                if (M.npar[e] && wj >= M.offs[e] && wj < M.offs[e] + M.npar[e]) a = e;
            if (a < 0) return false;
            std::vector<double> gv(n);
            for (int i = 0; i < n; ++i) gv[i] = M.w[i] * l1[(std::size_t) i * K + a];
            std::vector<double> sc(M.npar[a]);
            xt(a, gv.data(), sc.data());
            std::vector<double> col(p), v(p);
            for (int r = 0; r < p; ++r) col[r] = H(r, wj) + S(r, wj);
            for (int hh : holdv) col[hh] = 0.0;
            dense_matvec(Kinv, col.data(), v.data());
            for (double& x : v) x = -x;
            for (int hh : holdv) v[hh] = 0.0;
            v[wj] = 1.0;
            std::vector<double> ub(p);
            for (int r = 0; r < p; ++r) ub[r] = u[r] * v[r];
            out[nh + j] = sc[wj - M.offs[a]] - (r_sum(ub.data(), p) + trp(v)) / 2.0;
        }
        return true;
    }
    // ---------------------------------------------------- the Hessian
    // T[v] = -sum_i w_i omega_ab x_ia x_ib', omega_ab = sum_k l3_abk tv_k
    void contract3(const std::vector<std::vector<double>>& tv, Dense& T) {
        Model& M = E.M;
        T = Dense(p, p);
        const int K3 = K * K * K;
        for (int a = 0; a < K; ++a)
            for (int b = a; b < K; ++b) {
                if (!M.npar[a] || !M.npar[b]) continue;
                const RowNz& Ra = rows[a];
                const RowNz& Rb = rows[b];
                const int oa = M.offs[a], ob = M.offs[b];
                for (int i = 0; i < n; ++i) {
                    double om = 0.0;
                    for (int k = 0; k < K; ++k)
                        if (!tv[k].empty())
                            om += l3[(std::size_t) i * K3 + idx3(K, a, b, k)] * tv[k][i];
                    const double wv = -M.w[i] * om;
                    if (wv == 0.0) continue;
                    for (int s = Ra.start[i]; s < Ra.start[i + 1]; ++s) {
                        const double xa = wv * Ra.val[s];
                        const int ja = oa + Ra.col[s];
                        for (int t = Rb.start[i]; t < Rb.start[i + 1]; ++t)
                            T(ja, ob + Rb.col[t]) += xa * Rb.val[t];
                    }
                }
                if (a != b)
                    for (int j = 0; j < M.npar[a]; ++j)
                        for (int k = 0; k < M.npar[b]; ++k)
                            T(ob + k, oa + j) = T(oa + j, ob + k);
            }
    }
    // tr(M T[tv]) or tr(M U[tv, tu])
    double trace_form(const std::vector<std::vector<double>>& tv,
                      const std::vector<std::vector<double>>* tu) {
        Model& M = E.M;
        const int K3 = K * K * K, K4 = K3 * K;
        ldouble tot = 0.0;
        std::vector<double> om(n), prod(n);
        for (int a = 0; a < K; ++a)
            for (int b = a; b < K; ++b) {
                if (!M.npar[a] || !M.npar[b]) continue;
                std::fill(om.begin(), om.end(), 0.0);
                for (int k = 0; k < K; ++k) {
                    if (tv[k].empty()) continue;
                    if (!tu) {
                        for (int i = 0; i < n; ++i)
                            om[i] += l3[(std::size_t) i * K3 + idx3(K, a, b, k)] * tv[k][i];
                    } else {
                        for (int qq = 0; qq < K; ++qq) {
                            if ((*tu)[qq].empty()) continue;
                            for (int i = 0; i < n; ++i)
                                om[i] += l4[(std::size_t) i * K4 + idx4(K, a, b, k, qq)] *
                                         tv[k][i] * (*tu)[qq][i];
                        }
                    }
                }
                const std::vector<double>& g = Gab(a, b);
                for (int i = 0; i < n; ++i) prod[i] = M.w[i] * om[i] * g[i];
                const double s = r_sum(prod.data(), n);
                tot += (a == b ? 1.0 : 2.0) * s;
            }
        return -(double) tot;
    }
    void predictors(const std::vector<double>& v, std::vector<std::vector<double>>& tv) {
        Model& M = E.M;
        tv.assign(K, std::vector<double>());
        for (int k = 0; k < K; ++k)
            if (M.npar[k]) xv(k, v.data() + M.offs[k], tv[k]);
    }

    bool hessian(const std::vector<int>& mrow, const std::vector<int>& munit,
                 const std::vector<int>& mname, const std::vector<double>& mchain,
                 int nh, const std::vector<int>& gam,
                 const std::vector<double>& h1, const std::vector<double>& h2,
                 Dense& out) {
        if (!ensure_level(4) || !ensure_inverse()) return false;
        u_vector();
        ensure_uder();
        const int ng = (int) gam.size(), nv = nh + ng, nm = (int) mrow.size();
        std::vector<int> holdv = held_mode;
        for (int i = 0; i < nm; ++i)
            if (!uder_ok[munit[i]] || mname[i] < 0 || mname[i] >= uder[munit[i]].H) return false;
        // per row: S_m and c_m
        std::vector<Dense> Sm(nv, Dense(p, p));
        std::vector<std::vector<double>> cm(nh, std::vector<double>(p, 0.0));
        for (int i = 0; i < nm; ++i) {
            const UnitDer& ud = uder[munit[i]];
            const PenUnit& un = E.M.units[munit[i]];
            if (ud.grp) place_g(un, ud.dSg[mname[i]], Sm[mrow[i]]);
            else if (ud.vec) place_v(un, ud.dSv[mname[i]], Sm[mrow[i]]);
            else place(un, ud, ud.dS[mname[i]], Sm[mrow[i]]);
            for (int a = 0; a < ud.m; ++a) cm[mrow[i]][un.idx[a]] += ud.c[mname[i]][a];
        }
        // the mode's movements
        std::vector<std::vector<double>> bh(nv, std::vector<double>(p, 0.0));
        for (int m = 0; m < nh; ++m) {
            dense_matvec(Kinv, cm[m].data(), bh[m].data());
            for (double& x : bh[m]) x = -x;
        }
        Dense Hm(p, p);
        for (std::size_t t = 0; t < Hm.v.size(); ++t) Hm.v[t] = H.v[t] + S.v[t];
        for (int j = 0; j < ng; ++j) {
            const int w = gam[j];
            std::vector<double> col(p);
            for (int r = 0; r < p; ++r) col[r] = Hm(r, w);
            for (int hh : holdv) col[hh] = 0.0;
            dense_matvec(Kinv, col.data(), bh[nh + j].data());
            for (double& x : bh[nh + j]) x = -x;
            for (int hh : holdv) bh[nh + j][hh] = 0.0;
            bh[nh + j][w] = 1.0;
        }
        std::vector<std::vector<std::vector<double>>> tv(nv);
        std::vector<Dense> Tm(nv), MK(nv);
        for (int m = 0; m < nv; ++m) {
            predictors(bh[m], tv[m]);
            contract3(tv[m], Tm[m]);
            add_P(bh[m], Tm[m]);
            Dense Km(p, p);
            for (std::size_t t = 0; t < Km.v.size(); ++t) Km.v[t] = Sm[m].v[t] + Tm[m].v[t];
            dense_matmul(Mt, Km, MK[m]);
        }
        out = Dense(nv, nv);
        std::vector<double> tmp(p), rhs(p), bml(p), Jb(p);
        for (int m = 0; m < nv; ++m)
            for (int l = m; l < nv; ++l) {
                const bool both = m < nh && l < nh;
                // the second-order pieces of the pair: every ordered pair of
                // member lines of one unit with rows (m, l)
                double rho2 = 0.0, trS = 0.0;
                std::vector<double> c2(p, 0.0);
                if (both)
                    for (int i = 0; i < nm; ++i) {
                        if (mrow[i] != m) continue;
                        for (int j = 0; j < nm; ++j) {
                            if (mrow[j] != l || munit[j] != munit[i]) continue;
                            const UnitDer& ud = uder[munit[i]];
                            const PenUnit& un = E.M.units[munit[i]];
                            const int key = mname[i] * ud.H + mname[j];
                            rho2 += ud.rho2(mname[i], mname[j]);
                            trS += ud.grp ? trace_g(un, ud.d2Sg[key])
                                 : ud.vec ? trace_v(un, ud.d2Sv[key])
                                          : trace_unit(un, ud, ud.d2S[key]);
                            for (int a = 0; a < ud.m; ++a) c2[un.idx[a]] += ud.c2[key][a];
                        }
                    }
                // rhs = (S_l + T_l) b_m + S_m b_l + c_ml
                Dense SlTl(p, p);
                for (std::size_t t = 0; t < SlTl.v.size(); ++t) SlTl.v[t] = Sm[l].v[t] + Tm[l].v[t];
                dense_matvec(SlTl, bh[m].data(), rhs.data());
                dense_matvec(Sm[m], bh[l].data(), tmp.data());
                for (int r = 0; r < p; ++r) rhs[r] = rhs[r] + tmp[r] + c2[r];
                for (int hh : holdv) rhs[hh] = 0.0;
                dense_matvec(Kinv, rhs.data(), bml.data());
                for (double& x : bml) x = -x;
                for (int hh : holdv) bml[hh] = 0.0;
                std::vector<std::vector<double>> tml;
                predictors(bml, tml);
                double tr_dK = trS + trace_form(tv[l], &tv[m]) + trace_form(tml, nullptr);
                // tr(M E_ml) of the moving units: d2S/dbeta2[b_l, b_m] +
                // dS/dbeta[b_ml] + the movement of the m row's dS along b_l
                // and of the l row's along b_m
                for (std::size_t k = 0; k < uder.size(); ++k) {
                    if (!uder_ok[k] || !uder[k].moving) continue;
                    const PenUnit& un = E.M.units[k];
                    const UnitDer& ud = uder[k];
                    if (ud.grp) {
                        std::vector<Dense> B, X;
                        mvt_blocks(2, (int) k, 0, bh[l], &bh[m], B);
                        mvt_blocks(1, (int) k, 0, bml, nullptr, X);
                        for (int g = 0; g < ud.G; ++g)
                            for (std::size_t t = 0; t < B[g].v.size(); ++t) B[g].v[t] += X[g].v[t];
                        for (int i = 0; i < nm; ++i) {
                            if (munit[i] != (int) k) continue;
                            if (mrow[i] == m) {
                                mvt_blocks(3, (int) k, mname[i], bh[l], nullptr, X);
                                for (int g = 0; g < ud.G; ++g)
                                    for (std::size_t t = 0; t < B[g].v.size(); ++t) B[g].v[t] += X[g].v[t];
                            }
                            if (mrow[i] == l) {
                                mvt_blocks(3, (int) k, mname[i], bh[m], nullptr, X);
                                for (int g = 0; g < ud.G; ++g)
                                    for (std::size_t t = 0; t < B[g].v.size(); ++t) B[g].v[t] += X[g].v[t];
                            }
                        }
                        tr_dK += trace_g(un, B);
                        continue;
                    }
                    std::vector<double> Ed(ud.m);
                    for (int j = 0; j < ud.m; ++j) {
                        const int q = un.idx[j];
                        Ed[j] = ud.a4[j] * bh[l][q] * bh[m][q] + ud.a3[j] * bml[q];
                    }
                    for (int i = 0; i < nm; ++i) {
                        if (munit[i] != (int) k) continue;
                        if (mrow[i] == m)
                            for (int j = 0; j < ud.m; ++j)
                                Ed[j] += ud.a3h[mname[i]][j] * bh[l][un.idx[j]];
                        if (mrow[i] == l)
                            for (int j = 0; j < ud.m; ++j)
                                Ed[j] += ud.a3h[mname[i]][j] * bh[m][un.idx[j]];
                    }
                    tr_dK += trace_v(un, Ed);
                }
                double first;
                if (both) {
                    dense_matvec(Kmode, bh[l].data(), Jb.data());
                    for (int r = 0; r < p; ++r) tmp[r] = bh[m][r] * Jb[r];
                    first = r_sum(tmp.data(), p);
                } else if (m >= nh) {
                    const int w = gam[m - nh];
                    for (int r = 0; r < p; ++r) tmp[r] = Hm(w, r) * bh[l][r];
                    first = -r_sum(tmp.data(), p);
                } else {
                    const int w = gam[l - nh];
                    for (int r = 0; r < p; ++r) tmp[r] = Hm(w, r) * bh[m][r];
                    first = -r_sum(tmp.data(), p);
                }
                const double trKK = tr_prod(MK[l], MK[m]);
                const double v = -rho2 + first + trKK / 2.0 - tr_dK / 2.0;
                out(m, l) = v;
                out(l, m) = v;
            }
        // onto the free scale
        std::vector<double> g;
        if (!gradient(mrow, munit, mname, mchain, nh, gam, false, g)) return false;
        std::vector<double> c1(nv, 1.0), c2v(nv, 0.0);
        for (int r = 0; r < nh; ++r) { c1[r] = h1[r]; c2v[r] = h2[r]; }
        for (int m = 0; m < nv; ++m)
            for (int l = 0; l < nv; ++l) out(m, l) = out(m, l) * (c1[m] * c1[l]);
        for (int m = 0; m < nv; ++m) out(m, m) = out(m, m) + c2v[m] * g[m];
        return true;
    }
};

}  // namespace

#endif
