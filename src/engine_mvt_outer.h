// Kind 10 of the compiled engine, a multivariate Student t prior of held
// location (fixed() of an MvStudentTDistrib, not inverted, without a map),
// written to reproduce R's numbers bit for bit (SPEC_penalty 3.7,
// SPEC_penalty2 section 4, REVIEW_penalty 3.3, REVIEW_penalty2 finding 3):
//
//   mvt_order0()          the unit's L, sinv, logdet, nu, c0, half
//   mvt_moving_blocks()   the G diagonal blocks of penalty_hessian() at the
//                         coefficients (mvt_yscaled() with mv_quad())
//   mvt_groups()          the per-group pieces at the mode (4.1)
//   mvt_outer_order1()    gt, c, dS, ai, Na   (4.2)
//   mvt_outer_order2()    rho2, c2, d2S       (4.3)
//   mvt_beta_quadratic()  bq                  (4.5)
//
// and the functions of nu one order each: mvt_cp_d0(), mvt_cp_d1(),
// mvt_cp_d2(), mvt_D(), mvt_h_d2(). The outer fields are laid out as
// engine_outer_mvt() returns them and read_mvt_outer() reads them.
#ifndef CHART7_ENGINE_MVT_OUTER_H
#define CHART7_ENGINE_MVT_OUTER_H

#include "engine_charts.h"

namespace chart7 {

// the static description of a kind 10 unit
struct MvtStatic {
    ChartStatic ch;                  // never inverted
    int p = 0;
    std::vector<double> mu;          // the held location, p values
    std::vector<int> role;           // per hyperparameter, in pen@params order:
                                     // the chart's free value (0-based) or -1 for nu
    // odd p: the Horner coefficients of mvt_cp_d() for nu >= max(20, 4 p),
    // 40 each, computed in R as mvt_cp_d() computes them:
    // cf0[n] = -g[n + 1] * factorial(n - 1), cf1[n] = g[n + 1] * factorial(n),
    // cf2[n] = -g[n + 1] * factorial(n + 1), g = mvt_nu_coef(p)
    std::vector<double> cf0, cf1, cf2;
};

// ------------------------------------------------- the functions of nu
inline bool cp_big(double nu, int p) {
    return nu >= std::max(20.0, 4.0 * p);
}
inline double cp_horner(const std::vector<double>& cf, double x) {
    double s = 0.0;
    for (int i = (int) cf.size(); i >= 1; --i) s = (s + cf[i - 1]) * x;
    return s;
}

// mvt_cp_d(nu, p, 0)
inline double mvt_cp_d0(double nu, int p, const std::vector<double>& cf0) {
    if (p % 2 == 0) {
        double acc = (-((double) p / 2)) * std::log(2.0);
        for (int j = 1; j <= p / 2 - 1; ++j) acc = acc + std::log1p((2.0 * j) / nu);
        return acc;
    }
    if (cp_big(nu, p)) {
        const double s = cp_horner(cf0, 1 / nu);
        return ((-((double) p / 2)) * std::log(2.0)) + s;
    }
    return (R::lgammafn((nu + p) / 2) - R::lgammafn(nu / 2)) - (((double) p / 2) * std::log(nu));
}

// mvt_cp_d(nu, p, 1)
inline double mvt_cp_d1(double nu, int p, const std::vector<double>& cf1) {
    if (p % 2 == 0) {
        double acc = 0.0;
        for (int j = 1; j <= p / 2 - 1; ++j)
            acc = acc + std::expm1((-1.0) * std::log1p((2.0 * j) / nu));
        return ((1.0 * 1.0) * acc) / rpow(nu, 1.0);
    }
    if (cp_big(nu, p)) {
        const double x = 1 / nu;
        return cp_horner(cf1, x) * rpow(x, 1.0);
    }
    return (rpow(0.5, 1.0) * (R::psigamma((nu + p) / 2, 0.0) - R::psigamma(nu / 2, 0.0))) -
           (((((double) p / 2) * 1.0) * 1.0) / rpow(nu, 1.0));
}

// mvt_cp_d(nu, p, 2)
inline double mvt_cp_d2(double nu, int p, const std::vector<double>& cf2) {
    if (p % 2 == 0) {
        double acc = 0.0;
        for (int j = 1; j <= p / 2 - 1; ++j)
            acc = acc + std::expm1((-2.0) * std::log1p((2.0 * j) / nu));
        return (((-1.0) * 1.0) * acc) / rpow(nu, 2.0);
    }
    if (cp_big(nu, p)) {
        const double x = 1 / nu;
        return cp_horner(cf2, x) * rpow(x, 2.0);
    }
    return (rpow(0.5, 2.0) * (R::psigamma((nu + p) / 2, 1.0) - R::psigamma(nu / 2, 1.0))) -
           (((((double) p / 2) * (-1.0)) * 1.0) / rpow(nu, 2.0));
}

// mvt_D(u) = u / (1 + u) - log1p(u), by its series below |u| = 1e-3
inline double mvt_D(double u) {
    if (std::fabs(u) < 1e-3)
        return (u * u) * (-0.5 + u * (2.0 / 3.0 + u * (-0.75 + u * 0.8)));
    return (u / (1 + u)) - std::log1p(u);
}

// mvt_h_d(u, nu, p, 2) on the vector u of the G groups: R calls log1p()
// and expm1() of a_k on that vector, so its wrappers when G >= 2 and libm
// when G = 1 (REVIEW_penalty2 finding 3)
inline void mvt_h_d2(const std::vector<double>& u, double nu, int p, std::vector<double>& out) {
    const int n = (int) u.size();
    const bool wr = n >= 2;
    out.assign(n, 0.0);
    for (int g = 0; g < n; ++g) {
        const double x = u[g];
        const double lp = wr ? w_log1p(x) : std::log1p(x);
        const double ak = wr ? w_expm1((-2.0) * lp) : std::expm1((-2.0) * lp);
        double bk;
        if (std::fabs(x) < 0.1) {
            double s = 0.0;
            for (int m = 22; m >= 2; --m) {
                const double cf = rpow(-1.0, (double) m) * ((double) (m + 1) - 2.0 * 1.0);
                s = (s + cf) * x;
            }
            bk = s * x;
        } else {
            const double l1 = std::log1p(x);
            bk = std::expm1((-2.0) * l1) - 2.0 * std::expm1((-1.0) * l1);
        }
        out[g] = (((0.5 * (-1.0)) * 1.0) / rpow(nu, 1.0)) * (bk + (((double) p / nu) * ak));
    }
}

// ------------------------------------------------- order 0 of the unit
struct MvtUnit0 {
    Dense L, sinv;
    double logdet = 0, nu = 0, c0 = 0, half = 0;
};

// engine_penalty_units()'s kind 10 fields, from the chart at order 0
inline bool mvt_order0(const MvtStatic& ms, const ChartPt& ch, double nu, MvtUnit0& out) {
    if (ch.order < 0 || !ch.ok || !ch.L_ok) return false;
    const int p = ms.p;
    out.L = ch.L;
    out.sinv = ch.sinv;
    out.logdet = ch.logdet;
    out.nu = nu;
    out.c0 = (mvt_cp_d0(nu, p, ms.cf0) - (((double) p / 2) * std::log(M_PI))) - (0.5 * ch.logdet);
    out.half = (nu + p) / 2;
    return R_FINITE(out.c0) && R_FINITE(out.half);
}

// ------------------------------------------- the scaled residuals
// mvt_yscaled() of the rows of Y (n x p): ws (n x p), is and cw, with
// qt = colSums(forwardsolve(L, t(rt))^2) (mv_quad())
inline bool mvt_yscaled(const MvtStatic& ms, const ChartPt& ch, double nu, const Dense& Y,
                        Dense& ws, std::vector<double>& is, std::vector<double>& cw) {
    const int p = ms.p, n = Y.nr;
    Dense rt(n, p);
    std::vector<double> m(n);
    for (int g = 0; g < n; ++g) {
        double mx = 1.0;
        for (int j = 0; j < p; ++j) mx = std::max(mx, std::fabs(Y(g, j) - ms.mu[j]));
        m[g] = mx;
    }
    for (int j = 0; j < p; ++j)
        for (int g = 0; g < n; ++g) rt(g, j) = (Y(g, j) - ms.mu[j]) / m[g];
    Dense wt = mm(rt, ch.sinv);
    Dense X = transpose(rt);
    if (!trsolve(ch.L, false, X)) return false;
    ws = Dense(n, p);
    is.assign(n, 0.0);
    cw.assign(n, 0.0);
    std::vector<double> sq(p);
    for (int g = 0; g < n; ++g) {
        for (int j = 0; j < p; ++j) sq[j] = X(j, g) * X(j, g);
        const double qt = r_linesum(sq.data(), p);
        const double im = 1 / m[g], im2 = im * im;
        const double st = (nu * im2) + qt;
        is[g] = im2 / st;
        cw[g] = (nu + p) * is[g];
        for (int j = 0; j < p; ++j) ws(g, j) = wt(g, j) / (m[g] * st);
    }
    return true;
}

// distrib_hess_y() of one row: ((-cw) Sinv) + ((2 N) tcrossprod(ws))
inline Dense mvt_hess_slice(const ChartPt& ch, double N, const Dense& ws, int g, double cw) {
    const int p = ch.sinv.nr;
    Dense h(p, p);
    for (int b = 0; b < p; ++b)
        for (int a = 0; a < p; ++a)
            h(a, b) = ((-cw) * ch.sinv(a, b)) + ((2 * N) * (0.0 + ws(g, a) * ws(g, b)));
    return h;
}

// penalty_hessian()'s G diagonal blocks at the coefficients b (m = G p of
// them, group g at g p .. g p + p - 1): -(distrib_hess_y()) per group
inline bool mvt_moving_blocks(const MvtStatic& ms, const ChartPt& ch, double nu,
                              const double* b, int m, std::vector<Dense>& blocks) {
    const int p = ms.p, G = m / p;
    Dense Y(G, p);
    for (int g = 0; g < G; ++g)
        for (int j = 0; j < p; ++j) Y(g, j) = b[(std::size_t) g * p + j];
    Dense ws;
    std::vector<double> is, cw;
    if (!mvt_yscaled(ms, ch, nu, Y, ws, is, cw)) return false;
    const double N = nu + p;
    blocks.assign(G, Dense());
    for (int g = 0; g < G; ++g) blocks[g] = neg(mvt_hess_slice(ch, N, ws, g, cw[g]));
    return true;
}

// ------------------------------------------- 4.5 beta_quadratic()
// the four probe rows, distrib_hess_y() on each, and all.equal() of the
// 4 p^2 numbers against the first slice repeated, tolerance 1e-12
inline bool mvt_beta_quadratic(const MvtStatic& ms, const ChartPt& ch, double nu) {
    const int p = ms.p;
    static const double base[4] = {-1.73, -0.29, 0.61, 2.04};
    Dense T(4, p);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < p; ++j) T(i, j) = base[i] + 0.37 * (double) j;
    Dense ws;
    std::vector<double> is, cw;
    if (!mvt_yscaled(ms, ch, nu, T, ws, is, cw)) return false;
    const double N = nu + p;
    const int k = p * p;
    std::vector<double> v(4 * (std::size_t) k);
    for (int i = 0; i < 4; ++i) {
        Dense h = mvt_hess_slice(ch, N, ws, i, cw[i]);
        std::copy(h.v.begin(), h.v.end(), v.begin() + (std::size_t) i * k);
    }
    // all.equal.numeric(target = v, current = rep(v[1:k], 4))
    std::vector<double> tg, cu;
    for (std::size_t i = 0; i < v.size(); ++i) {
        const double t = v[i], c = v[i % k];
        if (ISNAN(t) != ISNAN(c)) return false;
        if (ISNAN(t) || t == c) continue;
        tg.push_back(t);
        cu.push_back(c);
    }
    if (tg.empty()) return true;
    const double Nn = (double) tg.size();
    std::vector<double> at(tg.size());
    for (std::size_t i = 0; i < tg.size(); ++i) at[i] = std::fabs(tg[i]) / Nn;
    double scale = 0 + r_sum(at.data(), at.size());
    if (!(R_FINITE(scale) && scale > 1e-12)) scale = 1;
    std::vector<double> d(tg.size());
    for (std::size_t i = 0; i < tg.size(); ++i) d[i] = std::fabs(tg[i] - cu[i]) / (Nn * scale);
    const double xy = r_sum(d.data(), d.size());
    return !ISNAN(xy) && !(xy > 1e-12);
}

// ------------------------------------------- 4.1 the groups at the mode
struct MvtGroups {
    int G = 0, p = 0, H = 0;
    double nu = 0, N = 0;
    Dense Y, R, w;                   // G x p
    std::vector<double> q, s, cw, dw;
    // per hyperparameter
    std::vector<Dense> wa, Sa;       // G x p, p x p
    std::vector<std::vector<double>> sa, ca, da;
    std::vector<double> Na;
};

inline std::vector<double> rowsums_prod(const Dense& A, const Dense& B) {
    const int G = A.nr, p = A.nc;
    std::vector<double> out(G), t(p);
    for (int g = 0; g < G; ++g) {
        for (int j = 0; j < p; ++j) t[j] = A(g, j) * B(g, j);
        out[g] = r_linesum(t.data(), p);
    }
    return out;
}

// the chart must be at order 1 or above
inline bool mvt_groups(const MvtStatic& ms, const ChartPt& ch, double nu, const double* b,
                       int m, MvtGroups& gr) {
    if (ch.order < 1) return false;
    const int p = ms.p, G = m / p, H = (int) ms.role.size();
    gr = MvtGroups();
    gr.G = G; gr.p = p; gr.H = H; gr.nu = nu;
    const double N = nu + p;
    gr.N = N;
    gr.Y = Dense(G, p);
    gr.R = Dense(G, p);
    for (int g = 0; g < G; ++g)
        for (int j = 0; j < p; ++j) {
            gr.Y(g, j) = b[(std::size_t) g * p + j];
            gr.R(g, j) = gr.Y(g, j) - ms.mu[j];
        }
    gr.w = mm(gr.R, ch.sinv);
    Dense X = transpose(gr.R);
    if (!trsolve(ch.L, false, X)) return false;
    gr.q.assign(G, 0.0); gr.s.assign(G, 0.0); gr.cw.assign(G, 0.0); gr.dw.assign(G, 0.0);
    std::vector<double> sq(p);
    for (int g = 0; g < G; ++g) {
        for (int j = 0; j < p; ++j) sq[j] = X(j, g) * X(j, g);
        gr.q[g] = r_linesum(sq.data(), p);
        gr.s[g] = nu + gr.q[g];
        gr.cw[g] = N / gr.s[g];
        gr.dw[g] = N / (gr.s[g] * gr.s[g]);
    }
    gr.wa.assign(H, Dense()); gr.Sa.assign(H, Dense());
    gr.sa.assign(H, std::vector<double>()); gr.ca.assign(H, std::vector<double>());
    gr.da.assign(H, std::vector<double>()); gr.Na.assign(H, 0.0);
    for (int h = 0; h < H; ++h) {
        const int k = ms.role[h];
        std::vector<double>& sa = gr.sa[h];
        std::vector<double>& ca = gr.ca[h];
        std::vector<double>& da = gr.da[h];
        ca.assign(G, 0.0); da.assign(G, 0.0);
        if (k < 0) {
            gr.Na[h] = 1.0;
            gr.wa[h] = Dense(G, p);
            gr.Sa[h] = Dense(p, p);
            sa.assign(G, 1.0);
            for (int g = 0; g < G; ++g)
                ca[g] = (gr.q[g] - p) / (gr.s[g] * gr.s[g]);
        } else {
            gr.Na[h] = 0.0;
            gr.wa[h] = mm(gr.R, ch.ai[k]);
            gr.Sa[h] = ch.ai[k];
            sa = rowsums_prod(gr.wa[h], gr.R);
            for (int g = 0; g < G; ++g)
                ca[g] = ((-N) * sa[g]) / (gr.s[g] * gr.s[g]);
        }
        for (int g = 0; g < G; ++g) {
            const double s = gr.s[g];
            da[g] = (gr.Na[h] / (s * s)) - (((2 * N) * sa[g]) / rpow(s, 3.0));
        }
    }
    return true;
}

// ------------------------------------------- 4.2 order 1
struct MvtOuter1 {
    std::vector<double> gt, c, dS, Na;   // H, H m, H G p p, H
    std::vector<Dense> ai;               // H
};

inline bool mvt_outer_order1(const MvtStatic& ms, const ChartPt& ch, const MvtGroups& gr,
                             MvtOuter1& out) {
    if (ch.order < 1) return false;
    const int p = gr.p, G = gr.G, H = gr.H, m = G * p;
    const double nu = gr.nu, N = gr.N;
    out.gt.assign(H, 0.0);
    out.c.assign((std::size_t) H * m, 0.0);
    out.dS.assign((std::size_t) H * G * p * p, 0.0);
    out.Na = gr.Na;
    out.ai.assign(H, Dense());
    std::vector<double> v(G);
    for (int h = 0; h < H; ++h) {
        const int k = ms.role[h];
        const std::vector<double>& sa = gr.sa[h];
        const std::vector<double>& ca = gr.ca[h];
        const std::vector<double>& da = gr.da[h];
        const Dense& wa = gr.wa[h];
        const Dense& Sa = gr.Sa[h];
        if (k >= 0) {
            for (int g = 0; g < G; ++g)
                v[g] = ((-0.5) * ch.dld[k]) - ((0.5 * gr.cw[g]) * sa[g]);
            out.ai[h] = ch.ai[k];
        } else {
            const double cp1 = mvt_cp_d1(nu, p, ms.cf1);
            for (int g = 0; g < G; ++g) {
                const double u = gr.q[g] / nu;
                v[g] = 0.5 * (((2 * cp1) + mvt_D(u)) + ((((double) p / nu) * u) / (1 + u)));
            }
            out.ai[h] = Dense(p, p);
        }
        out.gt[h] = -r_sum(v.data(), G);
        // c = -(cross + 0 * Y), read row by row
        for (int g = 0; g < G; ++g)
            for (int j = 0; j < p; ++j) {
                const double cross = ((-ca[g]) * gr.w(g, j)) - (gr.cw[g] * wa(g, j));
                out.c[(std::size_t) h * m + (std::size_t) g * p + j] = -(cross + 0 * gr.Y(g, j));
            }
        // dS: the group's block, column-major
        for (int g = 0; g < G; ++g) {
            const std::size_t off = ((std::size_t) h * G + g) * p * p;
            for (int j = 0; j < p; ++j)
                for (int i = 0; i < p; ++i) {
                    const double t =
                        (((ch.sinv(i, j) * (-ca[g])) + (Sa(i, j) * (-gr.cw[g]))) +
                         ((gr.w(g, i) * gr.w(g, j)) * (2 * da[g]))) +
                        (((wa(g, i) * gr.w(g, j)) + (gr.w(g, i) * wa(g, j))) * (2 * gr.dw[g]));
                    out.dS[off + (std::size_t) j * p + i] = -t;
                }
        }
    }
    bool fin = true;
    for (double x : out.gt) fin = fin && R_FINITE(x);
    for (double x : out.c) fin = fin && R_FINITE(x);
    for (double x : out.dS) fin = fin && R_FINITE(x);
    return fin;
}

// ------------------------------------------- 4.3 order 2
struct MvtOuter2 {
    std::vector<double> rho2, c2, d2S;   // H H, H H m, H H G p p
};

inline bool mvt_outer_order2(const MvtStatic& ms, const ChartPt& ch, const MvtGroups& gr,
                             MvtOuter2& out) {
    if (ch.order < 2) return false;
    const int p = gr.p, G = gr.G, H = gr.H, m = G * p, K = ms.ch.K;
    const double nu = gr.nu, N = gr.N;
    const std::size_t blk = (std::size_t) G * p * p;
    // the unordered pairs a <= b, computed once, then placed at (h, k) and (k, h)
    const int nup = H * (H + 1) / 2;
    std::vector<double> rho(nup);
    std::vector<std::vector<double>> c2u(nup), d2u(nup);
    std::vector<double> v(G), hd;
    const bool any_nu = std::find(ms.role.begin(), ms.role.end(), -1) != ms.role.end();
    double cp2 = 0.0;
    if (any_nu) {
        cp2 = mvt_cp_d2(nu, p, ms.cf2);
        std::vector<double> u(G);
        for (int g = 0; g < G; ++g) u[g] = gr.q[g] / nu;
        mvt_h_d2(u, nu, p, hd);
    }
    auto uidx = [&](int a, int b) { return a * H - a * (a - 1) / 2 + (b - a); };
    for (int a = 0; a < H; ++a)
        for (int b = a; b < H; ++b) {
            const int ka = ms.role[a], kb = ms.role[b];
            const bool nua = ka < 0, nub = kb < 0;
            // pair(): S_ab, w_ab, s_ab, c_ab, d_ab
            Dense S(p, p), wab(G, p);
            std::vector<double> sab(G, 0.0), cab(G), dab(G);
            int pk = -1;
            if (!nua && !nub) {
                pk = pair_index(K, ka, kb);
                S = ch.ai2[pk];
                wab = mm(gr.R, S);
                sab = rowsums_prod(wab, gr.R);
            }
            const std::vector<double>& saa = gr.sa[a];
            const std::vector<double>& sbb = gr.sa[b];
            const double Naa = gr.Na[a], Nbb = gr.Na[b];
            for (int g = 0; g < G; ++g) {
                const double s = gr.s[g];
                const double s3 = rpow(s, 3.0), s4 = rpow(s, 4.0);
                if (nua && nub) cab[g] = (2 * (p - gr.q[g])) / s3;
                else
                    cab[g] = (((((-Naa) * sbb[g]) / (s * s)) - ((Nbb * saa[g]) / (s * s))) -
                              ((N * sab[g]) / (s * s))) +
                             ((((2 * N) * saa[g]) * sbb[g]) / s3);
                dab[g] = (((((((-2) * Naa) * sbb[g]) / s3) - (((2 * Nbb) * saa[g]) / s3)) -
                           (((2 * N) * sab[g]) / s3)) +
                          ((((6 * N) * saa[g]) * sbb[g]) / s4));
            }
            // rho2
            for (int g = 0; g < G; ++g) {
                const double s = gr.s[g];
                if (!nua && !nub)
                    v[g] = (((-0.5) * ch.d2ld[pk]) + (((gr.cw[g] * sbb[g]) * saa[g]) / (2 * s))) -
                           ((0.5 * gr.cw[g]) * sab[g]);
                else if (nua && nub)
                    v[g] = cp2 - hd[g];
                else {
                    const std::vector<double>& sk = nua ? sbb : saa;
                    v[g] = ((-0.5) * sk[g]) * ((gr.q[g] - p) / (s * s));
                }
            }
            const int ui = uidx(a, b);
            rho[ui] = -r_sum(v.data(), G);
            // c2 = -(gyh + 0 * Y), read row by row
            const Dense& waa = gr.wa[a];
            const Dense& wbb = gr.wa[b];
            const std::vector<double>& caa = gr.ca[a];
            const std::vector<double>& cbb = gr.ca[b];
            std::vector<double>& c2 = c2u[ui];
            c2.assign(m, 0.0);
            for (int g = 0; g < G; ++g)
                for (int j = 0; j < p; ++j) {
                    const double gyh = ((((-cab[g]) * gr.w(g, j)) - (caa[g] * wbb(g, j))) -
                                        (cbb[g] * waa(g, j))) - (gr.cw[g] * wab(g, j));
                    c2[(std::size_t) g * p + j] = -(gyh + 0 * gr.Y(g, j));
                }
            // d2S, the eight terms left to right
            const Dense& Saa = gr.Sa[a];
            const Dense& Sbb = gr.Sa[b];
            const std::vector<double>& daa = gr.da[a];
            const std::vector<double>& dbb = gr.da[b];
            std::vector<double>& d2 = d2u[ui];
            d2.assign(blk, 0.0);
            auto sym = [&](const Dense& x, const Dense& y, int g, int i, int j) {
                return (x(g, i) * y(g, j)) + (y(g, i) * x(g, j));
            };
            for (int g = 0; g < G; ++g)
                for (int j = 0; j < p; ++j)
                    for (int i = 0; i < p; ++i) {
                        const double t =
                            (((((((ch.sinv(i, j) * (-cab[g])) + (Sbb(i, j) * (-caa[g]))) +
                                 (Saa(i, j) * (-cbb[g]))) +
                                (S(i, j) * (-gr.cw[g]))) +
                               ((gr.w(g, i) * gr.w(g, j)) * (2 * dab[g]))) +
                              (sym(wbb, gr.w, g, i, j) * (2 * daa[g]))) +
                             (sym(waa, gr.w, g, i, j) * (2 * dbb[g]))) +
                            ((sym(wab, gr.w, g, i, j) + sym(waa, wbb, g, i, j)) * (2 * gr.dw[g]));
                        d2[(std::size_t) g * p * p + (std::size_t) j * p + i] = -t;
                    }
        }
    // the H x H layout of engine_outer_mvt(): pair (h, k) at h H + k
    out.rho2.assign((std::size_t) H * H, 0.0);
    out.c2.assign((std::size_t) H * H * m, 0.0);
    out.d2S.assign((std::size_t) H * H * blk, 0.0);
    for (int h = 0; h < H; ++h)
        for (int k = 0; k < H; ++k) {
            const int ui = uidx(std::min(h, k), std::max(h, k));
            const std::size_t i = (std::size_t) h * H + k;
            out.rho2[i] = rho[ui];
            std::copy(c2u[ui].begin(), c2u[ui].end(), out.c2.begin() + i * m);
            std::copy(d2u[ui].begin(), d2u[ui].end(), out.d2S.begin() + i * blk);
        }
    bool fin = true;
    for (double x : out.rho2) fin = fin && R_FINITE(x);
    for (double x : out.c2) fin = fin && R_FINITE(x);
    for (double x : out.d2S) fin = fin && R_FINITE(x);
    return fin;
}

}  // namespace chart7

#endif
