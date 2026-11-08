// The compiled alternation of a model with a kinked penalty (piano_motore.md,
// phase 3b, D32 and D35; specification in misure_motore_fase3b/SPEC_kink.md
// with its errata). statmod_alternate() on the covered route: per pass the
// curvature of a scaled SCAD or MCP (statmod_curv()), the smooth block by the
// iwls loop on a subset (fit_smooth(), iwls_loop() with E.sub), every kinked
// block by coord_fit() around the coordinate-descent kernel of
// coord_descent.h, the relative change of the objective, and at the break
// alternation_settled(). Every operation follows the R route's order: the
// products of crossprod() and %*% as the BLAS forms them (r_crossprod_vec(),
// r_matvec()), R's sum() in long double (r_sum()).
//
// Included by engine_iwls.cpp after iwls_loop().

#include "coord_descent.h"
#include <set>

namespace {

// one kinked block of statmod_blocks()$sparse
struct KBlock {
    int unit = -1;                 // its unit in M.units
    int param = 0;                 // its equation, 0-based
    std::vector<int> cols;         // its columns within the equation, 0-based
    std::vector<int> other;        // the equation's other columns
    bool ctr = false;              // coord_centers()
    int i0 = -1;                   // the intercept's stacked position (0-based)
    bool curv_todo = false;        // statmod_curv() writes its curvature
    bool curv_centered = false;    // curv_centered()
    bool has_curv = false;         // the penalty carries a curvature
    bool has_prev = false;         // prev_kink (a path)
    double prev_kink = 0.0;
    std::string label;             // "param/term" of the history
    Rcpp::NumericMatrix X;         // X[, cols]
    std::vector<double> Xo;        // X[, other], n x |other|
};

struct Working {
    std::vector<double> w, z;
};

// R's identical() on two numeric vectors of the same length
inline bool same_vec(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (ISNAN(a[i]) || ISNAN(b[i])) {
            if (!(ISNAN(a[i]) && ISNAN(b[i]))) return false;
            continue;
        }
        if (a[i] != b[i]) return false;
    }
    return true;
}

// crossprod(w, X^2) column by column, X n x m column-major
inline void wxsq_cols(const double* X, int n, int m, const std::vector<double>& w,
                      std::vector<double>& out) {
    std::vector<double> x2((std::size_t) n * m);
    for (std::size_t e = 0; e < x2.size(); ++e) x2[e] = X[e] * X[e];
    out.assign(m, 0.0);
    r_crossprod_vec(x2.data(), n, m, w.data(), out.data());
}

class Kink {
public:
    Engine& E;
    std::vector<KBlock> blocks;
    bool expected = false;                       // the alternation's information
    std::vector<int> smooth;                     // 0-based
    std::vector<int> held_where;                 // 0-based
    std::vector<double> held_value;
    std::vector<std::vector<int>> smooth_groups; // subset numbering
    std::vector<char> smooth_froz;
    std::vector<char> kinked;                    // coef_labels()$kinked
    std::vector<int> pinned;                     // pinned_coords(), 0-based
    std::vector<int> aliased;
    int status = 0;                              // 0 ok, 2 the R route takes it

    explicit Kink(Engine& e) : E(e) {}

    double fn(const std::vector<double>& b) {
        double v;
        if (!E.fn(b, v)) throw RError();
        return v;
    }

    // coord_working(): the working weights and response of equation p, or
    // false where the curvature is not positive or the score not finite
    bool working(const std::vector<double>& beta, int p, bool expect, Working& wq) {
        Model& M = E.M;
        Point& P = E.point(beta);
        if (!P.inb) throw RError();
        E.ensure_parts(P);
        if (expect && E.me.emethod == 0) E.ensure_eparts(P);
        const int n = M.n, K = M.K;
        std::vector<double> s(n), h(n), gth(K), Hk((std::size_t) K * K);
        for (int i = 0; i < n; ++i) {
            E.reg.grad(M.fam, P.buf.data() + (std::size_t) i * M.pmax, gth.data());
            const double h1 = P.h1[(std::size_t) p * n + i];
            double gl;
            d7ls::link_grad(1, &gth[p], &h1, &gl);
            s[i] = M.w[i] * gl;
            E.hess_obs(P, i, expect, Hk.data());
            h[i] = (-M.w[i]) * Hk[p + (std::size_t) p * K];
        }
        for (int i = 0; i < n; ++i) {
            if (!(M.w[i] > 0)) continue;
            if (!R_FINITE(h[i]) || h[i] <= 0 || !R_FINITE(s[i])) return false;
        }
        wq.z.resize(n);
        const double* e0 = P.eta.data() + (std::size_t) p * n;
        for (int i = 0; i < n; ++i) {
            if (M.w[i] > 0) wq.z[i] = e0[i] + s[i] / h[i];
            else { wq.z[i] = e0[i]; h[i] = 0.0; }
        }
        wq.w = h;
        return true;
    }

    // coord_colsq() on a dense block: the centered weighted sums of squares,
    // the uncentered one where the centered falls below 1e-10 of it
    void colsq(const KBlock& B, const std::vector<double>& w,
               const std::vector<int>& k, const std::vector<double>& m,
               std::vector<double>& out) {
        const int n = E.M.n, nk = (int) k.size();
        std::vector<double> D2((std::size_t) n * nk), X2((std::size_t) n * nk);
        for (int c = 0; c < nk; ++c) {
            const double* x = B.X.begin() + (std::size_t) k[c] * n;
            const double mk = m[k[c]];
            for (int i = 0; i < n; ++i) {
                const double d = x[i] - mk;
                D2[(std::size_t) c * n + i] = d * d;
                X2[(std::size_t) c * n + i] = x[i] * x[i];
            }
        }
        std::vector<double> s(nk), raw(nk);
        r_crossprod_vec(D2.data(), n, nk, w.data(), s.data());
        r_crossprod_vec(X2.data(), n, nk, w.data(), raw.data());
        for (int c = 0; c < nk; ++c) if (s[c] <= 1e-10 * raw[c]) s[c] = raw[c];
        out = s;
    }

    // crossprod(w, X[, k]^2)
    void curv(const KBlock& B, const std::vector<double>& w, const std::vector<int>& k,
              std::vector<double>& out) {
        const int n = E.M.n, nk = (int) k.size();
        std::vector<double> Xk((std::size_t) n * nk);
        for (int c = 0; c < nk; ++c)
            std::copy(B.X.begin() + (std::size_t) k[c] * n,
                      B.X.begin() + (std::size_t) (k[c] + 1) * n,
                      Xk.begin() + (std::size_t) c * n);
        wxsq_cols(Xk.data(), n, nk, w, out);
    }

    // crossprod(X, v) / sum(v)-style weighted means of every column
    void means(const KBlock& B, const std::vector<double>& w, std::vector<double>& mw) {
        const int n = E.M.n, m = B.X.ncol();
        mw.assign(m, 0.0);
        r_crossprod_vec(B.X.begin(), n, m, w.data(), mw.data());
        const double sw = r_sum(w.data(), n);
        for (int j = 0; j < m; ++j) mw[j] = mw[j] / sw;
    }

    // coord_kinks(): 2 at(eps/2) - at(eps), at(h) the half difference of the
    // penalty's gradient at rep(h) and rep(-h); the first finite kink is 0
    void kinks(const PenUnit& u, std::vector<double>& out) {
        const std::size_t q = u.idx.size();
        auto grad_at = [&](double h, std::vector<double>& g) {
            std::vector<double> t(q);
            for (std::size_t j = 0; j < q; ++j)
                t[j] = u.has_map ? u.dg[j] * (0.0 + h) : 0.0 + h;
            g.assign(q, 0.0);
            E.kink_gradient(u, t, g);
        };
        auto at = [&](double h, std::vector<double>& a) {
            std::vector<double> up, dn;
            grad_at(h, up);
            grad_at(-h, dn);
            a.resize(q);
            for (std::size_t j = 0; j < q; ++j) a[j] = (up[j] - dn[j]) / 2;
        };
        const double eps = 1e-4;
        std::vector<double> a1, a2;
        at(eps / 2, a1);
        at(eps, a2);
        out.resize(q);
        for (std::size_t j = 0; j < q; ++j) out[j] = 2 * a1[j] - a2[j];
    }

    // penalty_prox_spec() of the unit restricted to `keep`, at the steps
    // `step` (one per kept coordinate), with the curvature `tcv` (SCAD, MCP);
    // false for NULL
    bool table(const PenUnit& u, const std::vector<int>& keep, const std::vector<double>& step,
               const std::vector<double>& tcv, Rcpp::NumericMatrix& cut,
               Rcpp::NumericMatrix& slope, Rcpp::NumericMatrix& icept) {
        const int m = (int) keep.size();
        const int k = (u.kind == 5 || u.kind == 6) ? 2 : (u.kind == 7 ? 4 : 3);
        cut = Rcpp::NumericMatrix(m, k);
        slope = Rcpp::NumericMatrix(m, k);
        icept = Rcpp::NumericMatrix(m, k);
        std::vector<double> d(m, 1.0), t(m);
        for (int a = 0; a < m; ++a) {
            if (u.has_map) d[a] = std::fabs(u.dg[keep[a]]);
            t[a] = u.has_map ? step[a] * (d[a] * d[a]) : step[a];
        }
        if (u.kind >= 7) {
            for (int a = 0; a < m; ++a) {
                const double te = t[a] * tcv[a];
                if (u.kind == 7 ? te >= u.kpar - 1 : te >= u.kpar) return false;
            }
        }
        for (int a = 0; a < m; ++a) {
            const double ta = t[a];
            if (u.kind == 5) {
                const double lam = u.klam;
                cut(a, 0) = ta * lam; cut(a, 1) = R_PosInf;
                slope(a, 0) = 0; slope(a, 1) = 1;
                icept(a, 0) = 0; icept(a, 1) = (-ta) * lam;
            } else if (u.kind == 6) {
                const double dd = 1 + ta * u.kc;
                cut(a, 0) = ta * u.ka; cut(a, 1) = R_PosInf;
                slope(a, 0) = 0; slope(a, 1) = 1 / dd;
                icept(a, 0) = 0; icept(a, 1) = ((-ta) * u.ka) / dd;
            } else if (u.kind == 7) {
                const double tt = ta * tcv[a], lam = u.klam / tcv[a], A = u.kpar;
                const double dd = 1 - tt / (A - 1);
                cut(a, 0) = tt * lam; cut(a, 1) = (1 + tt) * lam; cut(a, 2) = A * lam;
                cut(a, 3) = R_PosInf;
                slope(a, 0) = 0; slope(a, 1) = 1; slope(a, 2) = 1 / dd; slope(a, 3) = 1;
                icept(a, 0) = 0; icept(a, 1) = (-tt) * lam;
                icept(a, 2) = (((-tt) * A) * lam) / ((A - 1) * dd); icept(a, 3) = 0;
            } else {
                const double tt = ta * tcv[a], lam = u.klam / tcv[a], G = u.kpar;
                const double dd = 1 - tt / G;
                cut(a, 0) = tt * lam; cut(a, 1) = G * lam; cut(a, 2) = R_PosInf;
                slope(a, 0) = 0; slope(a, 1) = 1 / dd; slope(a, 2) = 1;
                icept(a, 0) = 0; icept(a, 1) = ((-tt) * lam) / dd; icept(a, 2) = 0;
            }
            if (u.has_map) {
                for (int c = 0; c < k; ++c) {
                    cut(a, c) = cut(a, c) / d[a];
                    icept(a, c) = icept(a, c) / d[a];
                }
            }
        }
        return true;
    }

    struct CoordOut {
        double value = 0.0;
        bool converged = false;
        int iterations = 0;
    };

    // coord_fit(): false where it returns NULL (the R route's proximal
    // fallback, which the engine does not reproduce)
    bool coord_fit(std::vector<double>& cur, KBlock& B, CoordOut& out) {
        Timer tm(T_SPQR);
        Model& M = E.M;
        PenUnit& u = M.units[B.unit];
        const int n = M.n, m = B.X.ncol(), p = B.param;
        const double tol = 1e-8;
        const int maxit = 100;
        if (m == 0) return false;
        {
            // the probe at a step of 1e-10
            std::vector<int> all(m);
            for (int j = 0; j < m; ++j) all[j] = j;
            std::vector<double> st(m, 1e-10), tc = u.cv;
            Rcpp::NumericMatrix c1, c2, c3;
            if (u.kind >= 7 && tc.size() != (std::size_t) m) return false;
            if (!table(u, all, st, tc, c1, c2, c3)) return false;
        }
        const bool scaled = u.kind >= 7 && B.has_curv;
        std::vector<double> c_prev;
        if (scaled) c_prev = u.cv;
        double c_gap = 0.0;
        bool done = false, stuck = false, have_prev = false;
        std::vector<double> pw, pz;
        int sweeps = 0;
        const double eps = DBL_EPSILON;
        for (int it = 1; it <= maxit; ++it) {
            Working wq;
            if (!working(cur, p, expected, wq)) {
                if (expected || !working(cur, p, true, wq)) return false;
            }
            std::vector<double> mw;
            if (B.ctr) means(B, wq.w, mw);
            std::vector<double> off(n, 0.0);
            if (!B.other.empty()) {
                std::vector<double> bo(B.other.size());
                for (std::size_t j = 0; j < B.other.size(); ++j)
                    bo[j] = cur[M.offs[p] + B.other[j]];
                r_matvec(B.Xo.data(), n, (int) B.other.size(), bo.data(), off.data());
            }
            std::vector<double> z(n);
            const double* ofs = M.off[p];
            for (int i = 0; i < n; ++i) z[i] = (wq.z[i] - off[i]) - (ofs ? ofs[i] : 0.0);
            std::vector<double> zc(z);
            double sw = 0.0, swz = 0.0;
            {
                std::vector<double> t(n);
                for (int i = 0; i < n; ++i) t[i] = wq.w[i] * z[i];
                swz = r_sum(t.data(), n);
                sw = r_sum(wq.w.data(), n);
            }
            if (B.ctr) {
                const double c = swz / sw;
                for (int i = 0; i < n; ++i) zc[i] = z[i] - c;
            }
            if (have_prev && same_vec(wq.w, pw)) {
                double dz = 0.0, mz = 0.0;
                for (int i = 0; i < n; ++i) {
                    dz = std::max(dz, std::fabs(zc[i] - pz[i]));
                    mz = std::max(mz, std::fabs(z[i]));
                }
                if (dz <= 64 * eps * mz && c_gap <= 1e-10) { done = true; break; }
            }
            pw = wq.w; pz = zc; have_prev = true;
            // the column curvatures: every one finite and positive, or NULL
            std::vector<double> v(m, NA_REAL);
            {
                std::vector<double> vall;
                wxsq_cols(B.X.begin(), n, m, wq.w, vall);
                for (double x : vall) if (!R_FINITE(x) || x <= 0) return false;
                if (!B.ctr) v = vall;
            }
            std::vector<double> b0(m);
            for (int j = 0; j < m; ++j) b0[j] = cur[u.idx[j]];
            std::vector<double> s_now;
            kinks(u, s_now);
            std::vector<double> s_prev;
            const bool use_prev = B.has_prev && s_now[0] > 0;
            if (use_prev) {
                s_prev.resize(m);
                const double r = B.prev_kink / s_now[0];
                for (int j = 0; j < m; ++j) s_prev[j] = s_now[j] * r;
            }
            std::vector<int> keep;
            screen(B, wq.w, z, b0, s_now, use_prev ? &s_prev : nullptr, keep);
            Rcpp::List res;
            for (;;) {
                std::vector<int> need;
                for (int j : keep) if (ISNAN(v[j])) need.push_back(j);
                if (!need.empty()) {
                    std::vector<double> vn;
                    if (B.ctr) colsq(B, wq.w, need, mw, vn); else curv(B, wq.w, need, vn);
                    for (std::size_t a = 0; a < need.size(); ++a) v[need[a]] = vn[a];
                }
                const int mk = (int) keep.size();
                std::vector<double> tcv(mk, 1.0);
                if (u.kind >= 7) {
                    for (int a = 0; a < mk; ++a) tcv[a] = u.cv[keep[a]];
                }
                if (scaled) {
                    std::vector<double> cj(mk), tgt(mk);
                    bool ok = true;
                    for (int a = 0; a < mk; ++a) {
                        double x = v[keep[a]];
                        if (u.has_map) x = x / (u.dg[keep[a]] * u.dg[keep[a]]);
                        tgt[a] = x;
                        cj[a] = std::sqrt(c_prev[keep[a]] * x);
                        if (!(R_FINITE(cj[a]) && cj[a] > 0)) ok = false;
                    }
                    if (ok) tcv = cj;
                    double gmax = 0.0;
                    for (int a = 0; a < mk; ++a) {
                        const double g = std::fabs(std::log(tcv[a] / tgt[a]));
                        if (ISNAN(g) || g > gmax) gmax = g;
                        if (ISNAN(g)) { gmax = g; break; }
                    }
                    c_gap = gmax;
                    for (int a = 0; a < mk; ++a) c_prev[keep[a]] = tcv[a];
                }
                std::vector<double> st(mk);
                for (int a = 0; a < mk; ++a) st[a] = 1 / v[keep[a]];
                Rcpp::NumericMatrix cut, slope, icept;
                if (!table(u, keep, st, tcv, cut, slope, icept)) return false;
                Rcpp::NumericVector zz(z.begin(), z.end()), ww(wq.w.begin(), wq.w.end()),
                    bb(b0.begin(), b0.end()), mm;
                if (B.ctr) mm = Rcpp::NumericVector(mw.begin(), mw.end());
                Rcpp::IntegerVector scr(mk);
                for (int a = 0; a < mk; ++a) scr[a] = keep[a];
                const bool cov = mk <= 32 && n > 8 * mk;
                res = sm7cd::coord_run(sm7cd::DenseCols(B.X), zz, ww, bb, cut, slope,
                                       icept, scr, 500, tol, cov, mm);
                sweeps += Rcpp::as<int>(res["sweeps"]);
                if (mk == m) break;
                Rcpp::NumericVector gr = res["grad"];
                std::vector<int> back;
                std::vector<char> in(m, 0);
                for (int j : keep) in[j] = 1;
                for (int j = 0; j < m; ++j)
                    if (std::fabs(gr[j]) > s_now[j] * (1 + 1e-10) && !in[j]) back.push_back(j);
                if (back.empty()) break;
                for (int j : back) keep.push_back(j);
                std::sort(keep.begin(), keep.end());
            }
            Rcpp::NumericVector nb = res["beta"];
            std::vector<double> nxt(cur);
            for (int j = 0; j < m; ++j) nxt[u.idx[j]] = nb[j];
            if (B.ctr) {
                std::vector<double> t(m);
                for (int j = 0; j < m; ++j) t[j] = mw[j] * nb[j];
                nxt[B.i0] = (cur[B.i0] + swz / sw) - r_sum(t.data(), m);
            }
            if (!scaled) {
                const double f_cur = fn(cur);
                double f_nxt = fn(nxt);
                double step = 1.0;
                std::vector<double> dir(cur.size());
                for (std::size_t j = 0; j < cur.size(); ++j) dir[j] = nxt[j] - cur[j];
                auto ok = [&](double f) {
                    return R_FINITE(f) && f <= f_cur + 64 * eps * std::max(1.0, std::fabs(f_cur));
                };
                while (!ok(f_nxt) && step > std::ldexp(1.0, -30)) {
                    step = step / 2;
                    for (std::size_t j = 0; j < cur.size(); ++j) nxt[j] = cur[j] + step * dir[j];
                    f_nxt = fn(nxt);
                }
                if (!ok(f_nxt)) { stuck = true; break; }
            }
            double moved = 0.0;
            for (int j = 0; j < m; ++j) moved = std::max(moved, std::fabs(nxt[u.idx[j]] - b0[j]));
            cur = nxt;
            if (moved < tol && c_gap <= 1e-10) { done = true; break; }
        }
        out.value = fn(cur);
        out.converged = done || stuck;
        out.iterations = sweeps;
        return true;
    }

    // coord_screen()
    void screen(const KBlock& B, const std::vector<double>& w, const std::vector<double>& z,
                const std::vector<double>& b0, const std::vector<double>& s_now,
                const std::vector<double>* s_prev, std::vector<int>& keep) {
        const int n = E.M.n, m = B.X.ncol();
        keep.clear();
        bool all = s_prev == nullptr;
        if (!all) {
            for (int j = 0; j < m; ++j)
                if (!R_FINITE((*s_prev)[j]) || !R_FINITE(s_now[j]) || s_now[j] <= 0) all = true;
        }
        if (all) { for (int j = 0; j < m; ++j) keep.push_back(j); return; }
        // x_times_b(): the nonzero columns only, unless every one is nonzero
        // or a coefficient is not finite
        std::vector<double> xb(n, 0.0);
        std::vector<int> nz;
        bool fin = true;
        for (int j = 0; j < m; ++j) { if (b0[j] != 0) nz.push_back(j); if (!R_FINITE(b0[j])) fin = false; }
        if ((int) nz.size() == m || !fin) {
            r_matvec(B.X.begin(), n, m, b0.data(), xb.data());
        } else if (!nz.empty()) {
            const int k = (int) nz.size();
            std::vector<double> Xk((std::size_t) n * k), bk(k);
            for (int c = 0; c < k; ++c) {
                std::copy(B.X.begin() + (std::size_t) nz[c] * n,
                          B.X.begin() + (std::size_t) (nz[c] + 1) * n, Xk.begin() + (std::size_t) c * n);
                bk[c] = b0[nz[c]];
            }
            r_matvec(Xk.data(), n, k, bk.data(), xb.data());
        }
        std::vector<double> r(n);
        for (int i = 0; i < n; ++i) r[i] = z[i] - xb[i];
        if (B.ctr) {
            std::vector<double> t(n);
            for (int i = 0; i < n; ++i) t[i] = w[i] * r[i];
            const double c = r_sum(t.data(), n) / r_sum(w.data(), n);
            for (int i = 0; i < n; ++i) r[i] = r[i] - c;
        }
        std::vector<double> wr(n), g(m);
        for (int i = 0; i < n; ++i) wr[i] = w[i] * r[i];
        r_crossprod_vec(B.X.begin(), n, m, wr.data(), g.data());
        for (int j = 0; j < m; ++j) {
            g[j] = std::fabs(g[j]);
            if (g[j] >= 2 * s_now[j] - (*s_prev)[j] || b0[j] != 0) keep.push_back(j);
        }
        if (keep.empty()) {
            int jm = 0;
            for (int j = 1; j < m; ++j) if (g[j] > g[jm]) jm = j;
            keep.push_back(jm);
        }
    }

    // statmod_curv(): the curvature of every scaled SCAD or MCP at `beta`
    void curvatures(const std::vector<double>& beta) {
        for (KBlock& B : blocks) {
            if (!B.curv_todo) continue;
            PenUnit& u = E.M.units[B.unit];
            Working wq;
            if (!working(beta, B.param, expected, wq)) {
                if (expected || !working(beta, B.param, true, wq)) continue;
            }
            const int m = B.X.ncol();
            std::vector<double> cj;
            if (B.curv_centered) {
                std::vector<double> mw;
                means(B, wq.w, mw);
                std::vector<int> all(m);
                for (int j = 0; j < m; ++j) all[j] = j;
                colsq(B, wq.w, all, mw, cj);
            } else {
                wxsq_cols(B.X.begin(), E.M.n, m, wq.w, cj);
            }
            if (u.has_map)
                for (int j = 0; j < m; ++j) cj[j] = cj[j] / (u.dg[j] * u.dg[j]);
            bool ok = cj.size() == u.cv.size();
            for (double x : cj) if (!R_FINITE(x) || x <= 0) ok = false;
            if (!ok) continue;
            u.cv = cj;
        }
    }

    // the information at `beta`, H alone (statmod_information_at(), assembled)
    Dense information(const std::vector<double>& beta) {
        Engine::Pieces pc;
        std::vector<int> keep = E.sub;
        std::vector<double> kb = E.base;
        E.sub.clear();
        // assembled: the full pieces with the square root switched off
        const bool hc = E.M.hasC;
        E.M.hasC = false;
        Dense S0 = E.M.Sd;
        E.M.Sd = Dense(E.M.p, E.M.p);
        bool ok = E.pieces(beta, expected, pc);
        E.M.Sd = S0;
        E.M.hasC = hc;
        E.sub = keep;
        E.base = kb;
        if (!ok) throw RError();
        return pc.Ad;
    }

    // the penalty's Hessian at `beta`: the constant units as R passed them,
    // every kinked unit's diagonal at its current coefficients
    Dense penalty_hessian(const std::vector<double>& beta) {
        Dense S = E.M.Sd;
        for (const PenUnit& u : E.M.units) {
            if (!kinked_kind(u.kind)) continue;
            std::vector<double> bu, t, h;
            E.unit_t(u, beta, bu, t);
            E.kink_hessian(u, t, h);
            for (std::size_t j = 0; j < u.idx.size(); ++j) S(u.idx[j], u.idx[j]) = h[j];
        }
        for (double& x : S.v) if (!R_FINITE(x)) x = 0.0;
        return S;
    }

    // solve(A, b) for a general matrix: dgesv, then the reciprocal condition
    // number against .Machine$double.eps (La_solve())
    static bool solve_general(const Dense& A, std::vector<double>& b) {
        int n = A.nr, nrhs = 1, info = 0, lda = n;
        if (n == 0) return true;
        Dense a = A;
        std::vector<int> ipiv(n);
        F77_CALL(dgesv)(&n, &nrhs, a.v.data(), &lda, ipiv.data(), b.data(), &n, &info);
        if (info != 0) return false;
        // anorm of the original matrix, dlange("1")
        std::vector<double> wk(4 * (std::size_t) n);
        std::vector<int> iw(n);
        const double anorm = F77_CALL(dlange)("1", &n, &n, A.v.data(), &n, wk.data() FCONE);
        double rcond = 0.0;
        F77_CALL(dgecon)("1", &n, a.v.data(), &n, &anorm, &rcond, wk.data(), iw.data(),
                         &info FCONE);
        if (rcond < DBL_EPSILON) return false;
        return true;
    }

    // deficient_coords(K, held) of a symmetric matrix, 0-based positions
    // within it; `held` (0-based) is left out of the test, as setdiff() does
    static std::vector<int> deficient(const Dense& K,
                                      const std::vector<int>& held = std::vector<int>()) {
        std::vector<int> outv;
        const int p = K.nr;
        if (p < 2) return outv;
        std::vector<char> isheld(p, 0);
        for (int j : held) if (j >= 0 && j < p) isheld[j] = 1;
        std::vector<int> rest;
        for (int j = 0; j < p; ++j)
            if (R_FINITE(K(j, j)) && K(j, j) > 0 && !isheld[j]) rest.push_back(j);
        const int r = (int) rest.size();
        if (r < 2) return outv;
        Dense A(r, r);
        for (int b = 0; b < r; ++b)
            for (int a = 0; a < r; ++a) {
                A(a, b) = K(rest[a], rest[b]);
                if (!R_FINITE(A(a, b))) return outv;
            }
        std::vector<double> s(r);
        for (int j = 0; j < r; ++j) s[j] = 1 / std::sqrt(A(j, j));
        Dense Ae(r, r);
        for (int b = 0; b < r; ++b)
            for (int a = 0; a < r; ++a) Ae(a, b) = A(a, b) * (s[a] * s[b]);
        // solve_pd(): the factor and its condition estimate
        Dense U;
        if (r_chol(Ae, U)) {
            ldouble mx = 0.0;
            for (int b = 0; b < r; ++b) {
                ldouble cs = 0.0;
                for (int a = 0; a < r; ++a) cs += std::fabs(Ae(a, b));
                if (b == 0 || cs > mx) mx = cs;
            }
            double anorm = (double) mx;
            if (R_FINITE(anorm) && anorm > 0) {
                std::vector<double> work(3 * (std::size_t) r);
                std::vector<int> iwork(r);
                double rc = 0.0;
                int info = 0, rr = r;
                std::vector<double> uf(U.v);
                F77_CALL(dpocon)("U", &rr, uf.data(), &rr, &anorm, &rc, work.data(),
                                 iwork.data(), &info FCONE);
                const double lmin = info != 0 ? 0.0 : rc * anorm;
                if (R_FINITE(lmin) && lmin > 1e-12 * anorm) return outv;
            }
        }
        QRResult q;
        if (!r_qr(Ae, 1e-7, q)) return outv;
        if (q.rank >= r) return outv;
        std::vector<char> kept(r, 0);
        for (int j = 0; j < q.rank; ++j) kept[q.pivot[j] - 1] = 1;
        for (int j = 0; j < r; ++j) if (!kept[j]) outv.push_back(rest[j]);
        std::sort(outv.begin(), outv.end());
        return outv;
    }

    // alternation_readings() and alternation_settled()
    bool settled(const std::vector<double>& beta) {
        Timer tm(T_DQR);
        Model& M = E.M;
        const int p = M.p;
        double mode = NA_REAL, kkt = NA_REAL;
        std::vector<char> drop(p, 0), zero(p, 0);
        for (int j = 0; j < p; ++j) {
            zero[j] = kinked[j] && beta[j] == 0;
            if (zero[j]) drop[j] = 1;
        }
        for (int j : aliased) drop[j] = 1;
        for (int j : pinned) drop[j] = 1;
        std::vector<double> g;
        bool have = true;
        Dense K;
        try {
            if (!E.gr(beta, g)) throw RError();
            Dense H = information(beta);
            Dense S = penalty_hessian(beta);
            K = Dense(p, p);
            for (std::size_t e = 0; e < K.v.size(); ++e) K.v[e] = H.v[e] + S.v[e];
        } catch (RError&) {
            have = false;
        }
        if (!have) return false;
        std::vector<int> keep;
        for (int j = 0; j < p; ++j)
            if (!drop[j] && R_FINITE(g[j]) && R_FINITE(K(j, j))) keep.push_back(j);
        if (keep.size() > 1) {
            const int k = (int) keep.size();
            Dense A0(k, k);
            for (int b = 0; b < k; ++b)
                for (int a = 0; a < k; ++a)
                    A0(a, b) = (K(keep[a], keep[b]) + K(keep[b], keep[a])) / 2;
            std::vector<int> dep = deficient(A0);
            if (!dep.empty()) {
                std::vector<char> out(k, 0);
                for (int d : dep) out[d] = 1;
                std::vector<int> k2;
                for (int a = 0; a < k; ++a) if (!out[a]) k2.push_back(keep[a]);
                keep = k2;
            }
        }
        if (keep.empty()) {
            mode = 0.0;
        } else {
            const int k = (int) keep.size();
            Dense A(k, k);
            for (int b = 0; b < k; ++b)
                for (int a = 0; a < k; ++a)
                    A(a, b) = (K(keep[a], keep[b]) + K(keep[b], keep[a])) / 2;
            std::vector<double> sol(k);
            for (int a = 0; a < k; ++a) sol[a] = g[keep[a]];
            if (solve_general(A, sol)) {
                bool fin = true;
                for (double x : sol) if (!R_FINITE(x)) fin = false;
                if (fin) {
                    std::vector<double> t(k);
                    for (int a = 0; a < k; ++a) t[a] = g[keep[a]] * sol[a];
                    mode = 0.5 * r_sum(t.data(), k);
                }
            }
        }
        std::vector<int> z;
        for (int j = 0; j < p; ++j) if (zero[j]) z.push_back(j);
        if (!z.empty()) {
            double kk = NA_REAL;
            try {
                std::vector<double> f0, f1;
                E.pen_sb(beta, f0);
                std::vector<double> bp(beta);
                const double tiny = std::sqrt(DBL_EPSILON) * DBL_EPSILON;
                for (int j : z) bp[j] = tiny;
                E.pen_sb(bp, f1);
                Dense H = information(beta);
                std::vector<double> dec;
                for (int j : z) {
                    const double pull = g[j] - f0[j];
                    const double kink = f1[j];
                    const double df = std::fabs(pull) - std::fabs(kink);
                    // pmax(0, NaN) is NaN in R
                    const double x = ISNAN(df) ? df : std::max(0.0, df);
                    const double dv = (x * x) / (2 * H(j, j));
                    if (R_FINITE(dv)) dec.push_back(dv);
                }
                if (!dec.empty()) {
                    kk = dec[0];
                    for (double x : dec) if (x > kk) kk = x;
                }
            } catch (RError&) {
                kk = NA_REAL;
            }
            kkt = ISNAN(kk) ? R_PosInf : kk;
        }
        const double lim = 1e-3;
        return (!ISNAN(mode) && mode <= lim) && (ISNAN(kkt) || kkt <= lim);
    }
};

}  // namespace
