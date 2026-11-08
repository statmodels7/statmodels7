// The penalty at an outer point, built in compiled code (phase 4, stage 3,
// D46): the list that engine_model() and engine_outer_pen() return in R
// (`S`, `C`, `units`), and the member chain of engine_outer_members(), from
// the outer vector. Kinds 1 (quadratic), 2 (additive), 3 (gaussian prior)
// and 4 (multivariate gaussian prior on a log-Cholesky, compound symmetry or
// AR(1) chart) are written here; a model with any other unit keeps the R
// construction. The arithmetic follows misure_motore_fase4/SPEC_penalty.md
// and its review: R's scalar and vector math (SPEC_penalty2 R2'), the order
// of every sum, the product by the identity as + 0.0.
// Included in engine_iwls.cpp, after the engine and before the compiled point.

#ifndef STATMODELS7_ENGINE_PEN_H
#define STATMODELS7_ENGINE_PEN_H

#include "engine_mvt_outer.h"

namespace pen7 {

// the static data of a chart and of a multivariate t prior for the chart7
// code (misure_motore_fase4/cpp_d48d49), as engine_pen_spec() writes them
inline chart7::LinkSpec read_link7(Rcpp::List s, const char* pre) {
    chart7::LinkSpec lk;
    const std::string nm = std::string(pre) + "_name";
    if (!s.containsElementNamed(nm.c_str())) return lk;
    const std::string cls = Rcpp::as<std::string>(s[nm]);
    lk.id = chart7::lf7().class_id(cls.c_str());
    Rcpp::NumericVector par = s[std::string(pre) + "_par"];
    lk.par.assign(par.begin(), par.end());
    lk.lwr = Rcpp::as<double>(s[std::string(pre) + "_lwr"]);
    lk.upr = Rcpp::as<double>(s[std::string(pre) + "_upr"]);
    return lk;
}
inline chart7::ChartStatic read_chart7(Rcpp::List s) {
    chart7::ChartStatic st;
    st.code = Rcpp::as<int>(s["code"]);
    st.p = Rcpp::as<int>(s["p"]);
    st.K = Rcpp::as<int>(s["K"]);
    st.inverted = Rcpp::as<bool>(s["inverted"]);
    if (st.code == chart7::CH_LC) {
        Rcpp::IntegerVector r = s["row"], c = s["col"];
        Rcpp::LogicalVector d = s["ondiag"];
        st.row.assign(r.begin(), r.end());
        st.col.assign(c.begin(), c.end());
        for (int i = 0; i < d.size(); ++i) st.ondiag.push_back((char) (d[i] ? 1 : 0));
    }
    st.ls = read_link7(s, "ls");
    st.lr = read_link7(s, "lr");
    if (st.code == chart7::CH_DRP) {
        Rcpp::IntegerVector r = s["arow"], c = s["acol"];
        st.arow.assign(r.begin(), r.end());
        st.acol.assign(c.begin(), c.end());
        st.aidx.assign((std::size_t) st.p * st.p, -1);
        for (std::size_t k = 0; k < st.arow.size(); ++k)
            st.aidx[(std::size_t) st.acol[k] * st.p + st.arow[k]] = (int) k;
    }
    return st;
}
inline Rcpp::NumericMatrix as_rmat7(const chart7::Dense& d) {
    Rcpp::NumericMatrix m(d.nr, d.nc);
    std::copy(d.v.begin(), d.v.end(), m.begin());
    return m;
}
inline Rcpp::List as_rlist7(const std::vector<chart7::Dense>& v) {
    Rcpp::List out(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) out[i] = as_rmat7(v[i]);
    return out;
}

// R's vector exp (arithmetic.c): 1 + x for |x| <= sqrt(eps)
inline double r_exp_vec(double x) {
    return (std::fabs(x) <= std::sqrt(DBL_EPSILON)) ? 1. + x : std::exp(x);
}

struct LinkRef {
    int id = -1;
    std::vector<double> par;
    double lwr = R_NegInf, upr = R_PosInf;
};

LinkRef read_link(Rcpp::List l) {
    LinkRef o;
    Registry& reg = registry();
    o.id = reg.lid(Rcpp::as<std::string>(l["name"]).c_str());
    Rcpp::NumericVector p = l["par"], b = l["bounds"];
    o.par.assign(p.begin(), p.end());
    o.lwr = b[0];
    o.upr = b[1];
    return o;
}

// linkinv() with link_bounds_clamp(), and the two derivatives at eta
inline void link_inv(const LinkRef& L, double eta, double& h, double& h1, double& h2) {
    Registry& reg = registry();
    reg.inv12p(L.id, L.par.data(), eta, &h, &h1, &h2);
    h = reg.clamp(h, L.lwr, L.upr);
}
inline double link_fun(const LinkRef& L, double th) {
    return registry().fun(L.id, L.par.data(), th);
}

struct UnitSpec {
    int kind = 0;
    std::vector<int> idx;        // 0-based stacked positions
    Rcpp::List ul, ulo;          // the unit's lists at the start (inner, outer)
    std::vector<double> theta;   // its hyperparameters, in pen@params order
    // kind 1
    double rank = 0;
    MatAny M;
    // kind 2
    std::vector<MatAny> mats;
    int add_m = 1;
    double add_rank = 0, add_tol = 0;
    std::vector<Dense> block_mats;
    // kind 3
    int si = 0;
    // kind 9: the positions of sigma and nu among the hyperparameters (-1
    // where held), their held values, the location
    int tsig = -1, tnu = -1;
    double sig0 = 0, nu0 = 0, mu = 0;
    // kinds 5 to 8: the positions of lambda and of the second parameter
    // (alpha, a or gamma; -1 where held), its held value, the map's
    // diagonal, the curvature
    int tl = -1, tp2 = -1;
    double p20 = 0;
    bool has_map = false;
    std::vector<double> dg, cv;
    // kind 4
    int code = 0, pd = 0, G = 0;
    bool inverted = false;
    std::vector<int> lc_row, lc_col;
    std::vector<char> lc_diag;
    LinkRef ls, lr;
    // kind 4 on dr_prod, and kind 10: the chart7 code
    bool c7 = false;
    chart7::ChartStatic cst;
    chart7::MvtStatic mst;   // nu0 above holds a held nu
};

struct PenSpec {
    int p = 0;
    bool sparse = false, need_C = false;
    std::vector<UnitSpec> units;
    // the members of the outer index (0-based) and the links of its rows
    std::vector<int> mrow, munit, mname;
    std::vector<LinkRef> rows;
    int nh = 0;
};

inline std::vector<double> dvec_(SEXP x) {
    Rcpp::NumericVector v(x);
    return std::vector<double>(v.begin(), v.end());
}

// .enet_logM() of one number
inline double enet_logM(double x) {
    if (!(x > 30)) return R::pnorm(-x, 0.0, 1.0, 1, 1) + (x * x) / 2 + 0.5 * std::log(2 * M_PI);
    const double u = 1 / x, u2 = u * u;
    return -std::log(x) + std::log1p(u2 * (-1 + u2 * (3 + u2 * (-15 + u2 * 105))));
}

PenSpec read_pen_spec(Rcpp::List ps) {
    PenSpec o;
    o.p = Rcpp::as<int>(ps["p"]);
    o.sparse = Rcpp::as<bool>(ps["sparse"]);
    o.need_C = Rcpp::as<bool>(ps["need_C"]);
    Rcpp::List us = ps["units"];
    for (int k = 0; k < us.size(); ++k) {
        Rcpp::List u = us[k];
        UnitSpec s;
        s.kind = Rcpp::as<int>(u["kind"]);
        Rcpp::IntegerVector ix = u["idx"];
        for (int v : ix) s.idx.push_back(v - 1);
        s.ul = u["ul"];
        s.ulo = u["ulo"];
        Rcpp::NumericVector th = u["theta"];
        s.theta.assign(th.begin(), th.end());
        if (s.kind == 1) {
            s.rank = Rcpp::as<double>(s.ul["rank"]);
            read_any(s.ul["M"], s.M);
        } else if (s.kind == 2) {
            Rcpp::List ms = s.ul["matrices"];
            for (int h = 0; h < ms.size(); ++h) { MatAny m; read_any(ms[h], m); s.mats.push_back(m); }
            Rcpp::List ad = s.ulo["add"];
            s.add_m = Rcpp::as<int>(ad["m"]);
            s.add_rank = Rcpp::as<double>(ad["rank"]);
            s.add_tol = Rcpp::as<double>(ad["tol"]);
            Rcpp::List bm = ad["block_mats"];
            for (int h = 0; h < bm.size(); ++h) {
                Rcpp::NumericMatrix B(bm[h]);
                Dense D(B.nrow(), B.ncol());
                std::copy(B.begin(), B.end(), D.v.begin());
                s.block_mats.push_back(D);
            }
        } else if (s.kind == 3) {
            s.si = Rcpp::as<int>(u["si"]) - 1;
        } else if (s.kind >= 5 && s.kind <= 8) {
            s.tl = Rcpp::as<int>(u["tl"]) - 1;
            s.tp2 = Rcpp::as<int>(u["tp2"]) - 1;
            s.p20 = Rcpp::as<double>(u["p20"]);
            SEXP mp = s.ul["map"];
            s.has_map = !Rf_isNull(mp);
            if (s.has_map) s.dg = dvec_(mp);
            if (s.kind >= 7) s.cv = dvec_(s.ul["cv"]);
        } else if (s.kind == 9) {
            s.tsig = Rcpp::as<int>(u["tsig"]) - 1;
            s.tnu = Rcpp::as<int>(u["tnu"]) - 1;
            s.sig0 = Rcpp::as<double>(s.ul["sigma"]);
            s.nu0 = Rcpp::as<double>(s.ul["nu"]);
            s.mu = Rcpp::as<double>(s.ul["mu"]);
        } else if (s.kind == 10) {
            Rcpp::List mv = u["mvt"];
            s.mst.ch = read_chart7(mv["chart"]);
            s.mst.p = s.mst.ch.p;
            s.mst.mu = dvec_(mv["mu"]);
            Rcpp::IntegerVector role = mv["role"];
            s.mst.role.assign(role.begin(), role.end());
            s.mst.cf0 = dvec_(mv["cf0"]);
            s.mst.cf1 = dvec_(mv["cf1"]);
            s.mst.cf2 = dvec_(mv["cf2"]);
            s.nu0 = Rcpp::as<double>(s.ul["nu"]);
        } else if (s.kind == 4 && u.containsElementNamed("chart7")) {
            s.c7 = true;
            s.cst = read_chart7(u["chart7"]);
            s.pd = s.cst.p;
            s.G = (int) s.idx.size() / s.pd;
        } else if (s.kind == 4) {
            s.code = Rcpp::as<int>(u["code"]);
            s.pd = Rcpp::as<int>(u["pd"]);
            s.G = (int) s.idx.size() / s.pd;
            s.inverted = Rcpp::as<bool>(u["inverted"]);
            if (s.code == 1) {
                Rcpp::IntegerVector rr = u["lc_row"], cc = u["lc_col"];
                Rcpp::LogicalVector dg = u["lc_diag"];
                for (int j = 0; j < rr.size(); ++j) {
                    s.lc_row.push_back(rr[j] - 1);
                    s.lc_col.push_back(cc[j] - 1);
                    s.lc_diag.push_back(dg[j] == TRUE);
                }
            } else {
                s.ls = read_link(u["link_scale"]);
                s.lr = read_link(u["link_rho"]);
            }
        }
        o.units.push_back(s);
    }
    Rcpp::List mb = ps["members"];
    Rcpp::IntegerVector mr = mb["mrow"], mu = mb["munit"], mn = mb["mname"];
    for (int i = 0; i < mr.size(); ++i) {
        o.mrow.push_back(mr[i] - 1);
        o.munit.push_back(mu[i] - 1);
        o.mname.push_back(mn[i] - 1);
    }
    Rcpp::List rl = mb["links"];
    for (int r = 0; r < rl.size(); ++r) o.rows.push_back(read_link(rl[r]));
    o.nh = (int) o.rows.size();
    return o;
}

// the chart of a kind 4 unit at its free values: L, the inverse, log|Sigma|
// as mv_matrix_pieces() returns them; false where R signals an error
struct Chart {
    Dense L, sinv;
    double logdet = 0;
    std::vector<double> h;   // the outer chart's h (codes 2 and 3)
};

// chol_pd(): the eigenvalue test, then t(chol(m))
inline bool chol_pd(const Dense& m, Dense& L) {
    for (double x : m.v) if (!R_FINITE(x)) return false;
    std::vector<double> ev;
    Dense V;
    if (!r_eigen_sym(m, true, ev, V)) return false;
    double mx = ev[0], mn = ev[0];
    for (double x : ev) { mx = std::max(mx, x); mn = std::min(mn, x); }
    if (mx <= 0 || mn <= 1e-14 * mx) return false;
    Dense U;
    if (!r_chol(m, U)) return false;
    L = transpose(U);
    return true;
}

inline bool chart_at(const UnitSpec& u, const std::vector<double>& eta, Chart& c) {
    const int p = u.pd;
    Dense m(p, p), si(p, p);
    double ld = 0;
    if (u.code == 1) {
        // chol_assemble(): exp() on the vector of diagonal values
        c.L = Dense(p, p);
        int nd = 0;
        for (char d : u.lc_diag) if (d) ++nd;
        for (std::size_t j = 0; j < eta.size(); ++j) {
            double v = eta[j];
            if (u.lc_diag[j]) v = nd >= 2 ? r_exp_vec(eta[j]) : std::exp(eta[j]);
            if (u.lc_row.at(j) >= p || u.lc_col.at(j) >= p)
                Rcpp::stop("engine_pen: a log-Cholesky position is out of range");
            c.L(u.lc_row[j], u.lc_col[j]) = v;
        }
        std::vector<double> dgv;
        for (std::size_t j = 0; j < eta.size(); ++j) if (u.lc_diag[j]) dgv.push_back(eta[j]);
        ld = 2 * r_sum(dgv.data(), dgv.size());
        if (u.inverted) {
            // tcrossprod(L): dsyrk, the upper triangle copied
            int pp = p;
            const double one = 1.0, zero = 0.0;
            F77_CALL(dsyrk)("U", "N", &pp, &pp, &one, c.L.v.data(), &pp, &zero, m.v.data(),
                            &pp FCONE FCONE);
            for (int j = 0; j < p; ++j)
                for (int i = j + 1; i < p; ++i) m(i, j) = m(j, i);
        } else {
            // backsolve(t(L), forwardsolve(L, I)), column by column
            Dense Lt = transpose(c.L);
            for (int j = 0; j < p; ++j) {
                std::vector<double> x(p, 0.0);
                x[j] = 1.0;
                if (!r_trisolve(c.L, false, x)) return false;
                if (!r_trisolve(Lt, true, x)) return false;
                for (int i = 0; i < p; ++i) si(i, j) = x[i];
            }
        }
    } else {
        double v, v1, v2, r, r1, r2;
        link_inv(u.ls, eta[0], v, v1, v2);
        link_inv(u.lr, eta[1], r, r1, r2);
        c.h = {v, v1, v2, r, r1, r2};
        if (u.code == 2) {
            for (int j = 0; j < p; ++j)
                for (int i = 0; i < p; ++i) m(i, j) = v * (i == j ? (1 + r * 0) : (0 + r));
            const double one_plus = p * R::plogis(eta[1], 0.0, 1.0, 1, 0);
            const double one_minus = ((double) p / (p - 1)) * R::plogis(-eta[1], 0.0, 1.0, 1, 0);
            const double q = r / one_plus, den = v * one_minus;
            for (int j = 0; j < p; ++j)
                for (int i = 0; i < p; ++i) si(i, j) = ((i == j ? 1.0 : 0.0) - q * 1.0) / den + 0.0;
            ld = p * std::log(v) + std::log((double) p) + R::plogis(eta[1], 0.0, 1.0, 1, 1) +
                 (p - 1) * (std::log((double) p) - std::log((double) (p - 1)) +
                            R::plogis(-eta[1], 0.0, 1.0, 1, 1));
        } else {
            // AR(1): r^k by R_pow, r * r at k = 2
            std::vector<double> w(p);
            for (int k = 0; k < p; ++k)
                w[k] = k == 0 ? 1.0 : (k == 1 ? r : (k == 2 ? r * r : R_pow(r, (double) k)));
            for (int j = 0; j < p; ++j)
                for (int i = 0; i < p; ++i) m(i, j) = v * w[std::abs(i - j)];
            // sech2() and log_sech2() on one number: libm's exp and log1p
            const double z = eta[1];
            const double a = std::exp(-2 * std::fabs(z));
            const double om = 4 * a / ((1 + a) * (1 + a));
            Dense tm(p, p);
            for (int i = 0; i < p; ++i) tm(i, i) = (i == 0 || i == p - 1) ? 1.0 : 1 + r * r;
            for (int i = 0; i + 1 < p; ++i) { tm(i, i + 1) = -r; tm(i + 1, i) = -r; }
            const double den = v * om;
            for (int j = 0; j < p; ++j)
                for (int i = 0; i < p; ++i) si(i, j) = tm(i, j) / den + 0.0;
            const double ls2 = 2 * std::log(2.0) - 2 * std::fabs(z) -
                               2 * std::log1p(std::exp(-2 * std::fabs(z)));
            ld = p * std::log(v) + (p - 1) * ls2;
        }
        if (!chol_pd(m, c.L)) return false;
    }
    if (u.inverted) {
        c.sinv = m;
        c.logdet = -ld;
    } else {
        c.sinv = Dense(p, p);
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < p; ++i) c.sinv(i, j) = (si(i, j) + si(j, i)) / 2;
        c.logdet = ld;
    }
    return true;
}

// a list with the same elements, without copying them
Rcpp::List shallow(const Rcpp::List& x) {
    Rcpp::List o(x.size());
    for (R_xlen_t i = 0; i < x.size(); ++i) o[i] = x[i];
    o.attr("names") = x.attr("names");
    return o;
}

Rcpp::NumericMatrix as_rmat(const Dense& D) {
    Rcpp::NumericMatrix R(D.nr, D.nc);
    std::copy(D.v.begin(), D.v.end(), R.begin());
    return R;
}

// a dense matrix as a dgCMatrix of its nonzero entries
// returned protected: a bare SEXP held across allocations is collected
Rcpp::RObject as_dgc(const Dense& D) {
    std::vector<int> ii, pp(1, 0);
    std::vector<double> xx;
    for (int j = 0; j < D.nc; ++j) {
        for (int i = 0; i < D.nr; ++i) {
            const double v = D(i, j);
            if (v != 0.0 || ISNAN(v)) { ii.push_back(i); xx.push_back(v); }
        }
        pp.push_back((int) ii.size());
    }
    Rcpp::S4 s("dgCMatrix");
    s.slot("i") = Rcpp::IntegerVector(ii.begin(), ii.end());
    s.slot("p") = Rcpp::IntegerVector(pp.begin(), pp.end());
    s.slot("x") = Rcpp::NumericVector(xx.begin(), xx.end());
    s.slot("Dim") = Rcpp::IntegerVector::create(D.nr, D.nc);
    return Rcpp::RObject((SEXP) s);
}

// distrib_hess_y() of the t at y (update_penalty()'s transcription of
// student_t1_dy2_cpp())
inline double t_dy2(double y, double mu, double s, double nu) {
    const double ik = 1.0 / (nu - 0.0);
    const double z = (y - mu) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    const double uu = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    return t * (-ik * t + ik * uu - t + uu) / std::pow(s, 2);
}

// beta_quadratic() of a t prior of block 1: whether its distrib_hess_y() is
// the same at four probe values, read by all.equal(tolerance = 1e-12)
inline bool t_beta_quadratic(double mu, double s, double nu) {
    const double pr[4] = {-1.73, -0.29, 0.61, 2.04};
    double v[4];
    for (int i = 0; i < 4; ++i) v[i] = t_dy2(pr[i] + 0.37 * 0.0, mu, s, nu);
    for (double x : v) if (ISNAN(x)) return false;
    std::vector<double> tg, cu;
    for (int i = 0; i < 4; ++i) if (!(v[i] == v[0])) { tg.push_back(v[i]); cu.push_back(v[0]); }
    if (tg.empty()) return true;
    const double N = (double) tg.size();
    long double sc = 0.0L;
    for (double x : tg) sc += std::fabs(x) / N;
    double scale = (double) sc;
    if (!(scale > 1e-12)) scale = 1.0;
    long double xy = 0.0L;
    for (std::size_t i = 0; i < tg.size(); ++i) xy += std::fabs(tg[i] - cu[i]) / (N * scale);
    return (double) xy <= 1e-12;
}

// a unit block added at idx, entry by entry (hess[idx, idx] + blk)
inline void add_block(Dense& S, const std::vector<int>& idx, const Dense& B) {
    const int m = (int) idx.size();
    if (B.nr != m || B.nc != m) Rcpp::stop("engine_pen: a unit's block does not match its positions");
    for (int i : idx) if (i < 0 || i >= S.nr) Rcpp::stop("engine_pen: a unit's position is out of range");
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < m; ++i) S(idx[i], idx[j]) = S(idx[i], idx[j]) + B(i, j);
}
inline Dense dense_of_any(const MatAny& A) {
    if (!A.sp) return A.d;
    Dense D((int) A.s.rows(), (int) A.s.cols());
    for (int j = 0; j < A.s.outerSize(); ++j)
        for (SpMat::InnerIterator it(A.s, j); it; ++it) D((int) it.row(), j) = it.value();
    return D;
}

struct PenOut {
    Rcpp::List pen, pen_outer;
    // per unit, the outer pieces of a multivariate t prior (read_mvt_outer())
    std::vector<Rcpp::RObject> mvt;
    std::vector<double> mchain, h1, h2, vals;
    bool ok = false;
};

// the penalty at the outer vector `v` (its first nh entries) and the
// coefficients `b`, which the gaussian prior's block reads (+ 0 * t)
inline bool build_pen(const PenSpec& ps, const std::vector<double>& v,
                      const std::vector<double>& b, PenOut& out,
                      const std::vector<std::vector<double>>* cvs = nullptr) {
    out = PenOut();
    // eta_to_hyper(): every member writes its row's value into its unit
    std::vector<std::vector<double>> th(ps.units.size());
    for (std::size_t k = 0; k < ps.units.size(); ++k) th[k] = ps.units[k].theta;
    std::vector<double> rv(ps.nh, NA_REAL);
    std::vector<char> seen(ps.nh, 0);
    out.mchain.resize(ps.mrow.size());
    for (std::size_t i = 0; i < ps.mrow.size(); ++i) {
        const int r = ps.mrow[i];
        double h, h1, h2;
        link_inv(ps.rows[r], v[r], h, h1, h2);
        th.at(ps.munit[i]).at(ps.mname[i]) = h;
        if (!seen[r]) { rv[r] = h; seen[r] = 1; }
        // the derivative at linkfun(v), not at eta (SPEC_penalty 4)
        double e = link_fun(ps.rows[r], h), g, g1, g2;
        registry().inv12p(ps.rows[r].id, ps.rows[r].par.data(), e, &g, &g1, &g2);
        out.mchain[i] = g1;
    }
    out.h1.resize(ps.nh);
    out.h2.resize(ps.nh);
    for (int r = 0; r < ps.nh; ++r) {
        double e = link_fun(ps.rows[r], rv[r]), g, g1, g2;
        registry().inv12p(ps.rows[r].id, ps.rows[r].par.data(), e, &g, &g1, &g2);
        out.h1[r] = g1;
        out.h2[r] = g2;
    }
    out.vals = rv;

    Dense S(ps.p, ps.p);
    Rcpp::List units(ps.units.size()), unitso(ps.units.size());
    out.mvt.assign(ps.units.size(), Rcpp::RObject());
    for (std::size_t k = 0; k < ps.units.size(); ++k) {
        const UnitSpec& u = ps.units[k];
        const std::vector<double>& t = th[k];
        Rcpp::List ul = shallow(u.ul), ulo = shallow(u.ulo);
        const int m = (int) u.idx.size();
        if (u.kind == 1) {
            const double lam = t[0];
            const double a = lam / 2, c1 = u.rank / 2 * std::log(lam);
            for (Rcpp::List* l : {&ul, &ulo}) {
                (*l)["lam"] = lam;
                (*l)["a"] = a;
                (*l)["c1"] = c1;
            }
            Dense B = dense_of_any(u.M);
            for (double& x : B.v) x = lam * x;
            add_block(S, u.idx, B);
        } else if (u.kind == 2) {
            const int Kh = (int) u.mats.size();
            Dense B;
            std::vector<Dense> P(Kh);
            for (int h = 0; h < Kh; ++h) {
                P[h] = dense_of_any(u.mats[h]);
                Dense Q = P[h];
                for (double& x : Q.v) x = t[h] * x;
                if (h == 0) B = Q;
                else for (std::size_t z = 0; z < B.v.size(); ++z) B.v[z] = B.v[z] + Q.v[z];
            }
            // additive_sum(): logpdet, or NaN where the rank identity fails
            double val = 0.0;
            std::vector<double> dlog;
            Dense d2log;
            const bool one = u.add_m <= 1;
            double logpdet = R_NaN;
            const int rk = (int) std::llround(one ? u.add_rank : u.add_rank / u.add_m);
            if (additive_stable_c(one ? P : u.block_mats, t, rk, u.add_tol, val, dlog, d2log)) {
                if (!one) {
                    val = u.add_m * val;
                    for (double& z : dlog) z = u.add_m * z;
                }
                std::vector<double> ld(Kh);
                for (int h = 0; h < Kh; ++h) ld[h] = t[h] * dlog[h];
                if (R_FINITE(val) &&
                    std::fabs(r_sum(ld.data(), Kh) - u.add_rank) <= 1e-3 * u.add_rank)
                    logpdet = val;
            }
            // NaN goes to the R route (Q3)
            if (ISNAN(logpdet)) return false;
            const double c1 = 0.5 * logpdet;
            for (Rcpp::List* l : {&ul, &ulo}) {
                if (u.mats[0].sp) (*l)["S"] = as_dgc(B);
                else (*l)["S"] = as_rmat(B);
                (*l)["c1"] = c1;
                (*l)["lams"] = Rcpp::NumericVector(t.begin(), t.end());
            }
            add_block(S, u.idx, B);
        } else if (u.kind == 3) {
            const double sg = t[u.si];
            for (Rcpp::List* l : {&ul, &ulo}) (*l)["sigma"] = sg;
            const double hj = (-1) / (sg * sg);
            Dense B(m, m);
            for (int j = 0; j < m; ++j) B(j, j) = -(hj + 0 * b[u.idx[j]]);
            add_block(S, u.idx, B);
        } else if (u.kind >= 5 && u.kind <= 8) {
            const double lam = t[u.tl];
            const double p2 = u.tp2 >= 0 ? t[u.tp2] : u.p20;
            std::vector<double> cv = u.cv;
            if (u.kind >= 7 && cvs && !(*cvs)[k].empty()) cv = (*cvs)[k];
            for (Rcpp::List* l : {&ul, &ulo}) {
                (*l)["lambda"] = lam;
                if (u.kind == 5) (*l)["L"] = std::log(lam / 2);
                if (u.kind == 6) {
                    const double a = lam * p2, c = lam * (1 - p2);
                    const double x = a / std::sqrt(c);
                    (*l)["a"] = a;
                    (*l)["c"] = c;
                    (*l)["logz"] = (std::log(2.0) - 0.5 * std::log(c)) + enet_logM(x);
                }
                if (u.kind >= 7) {
                    (*l)["par"] = p2;
                    (*l)["cv"] = Rcpp::NumericVector(cv.begin(), cv.end());
                }
            }
            // kink_hessian() at the coefficients
            const double kc = lam * (1 - p2);
            Dense B(m, m);
            for (int j = 0; j < m; ++j) {
                const double bj = b[u.idx[j]];
                const double tj = u.has_map ? u.dg[j] * bj : bj;
                double in;
                bool neg = false;
                if (u.kind == 5) { in = 0.0 + 0.0 * tj; neg = true; }
                else if (u.kind == 6) { in = -kc + 0.0 * tj; neg = true; }
                else {
                    const double c = cv[j], lm = lam / c, a = p2, ab = std::fabs(tj);
                    if (u.kind == 7) {
                        const bool r2 = !(ab <= lm) && !(ab > a * lm);
                        in = c * (r2 ? -1 / (a - 1) : 0.0);
                    } else {
                        in = c * (ab <= a * lm ? -1 / a : 0.0);
                    }
                }
                const double vv = u.has_map ? u.dg[j] * (u.dg[j] * in) : in;
                B(j, j) = neg ? -vv : vv;
            }
            add_block(S, u.idx, B);
        } else if (u.kind == 9) {
            const double sg = u.tsig >= 0 ? t[u.tsig] : u.sig0;
            const double nu = u.tnu >= 0 ? t[u.tnu] : u.nu0;
            for (Rcpp::List* l : {&ul, &ulo}) {
                (*l)["sigma"] = sg;
                (*l)["nu"] = nu;
            }
            ulo["bq"] = t_beta_quadratic(u.mu, sg, nu);
            Dense B(m, m);
            for (int j = 0; j < m; ++j) {
                const double bj = b[u.idx[j]];
                B(j, j) = -(t_dy2(bj, u.mu, sg, nu) + 0.0 * bj);
            }
            add_block(S, u.idx, B);
        } else if (u.kind == 4 && u.c7) {
            // dr_prod: the chart and, for the outer context, its pieces
            // (engine_outer_chart() code 9)
            chart7::ChartPt pt;
            if (!chart7::chart_order0(u.cst, t.data(), pt) || !pt.ok || !pt.L_ok) return false;
            if (!chart7::chart_order1(u.cst, pt) || !chart7::chart_order2(u.cst, pt)) return false;
            for (Rcpp::List* l : {&ul, &ulo}) {
                (*l)["L"] = as_rmat7(pt.L);
                (*l)["sinv"] = as_rmat7(pt.sinv);
                (*l)["logdet"] = pt.logdet;
            }
            const int K = u.cst.K;
            Rcpp::List a2(K * K);
            Rcpp::NumericVector d2(K * K);
            for (int a = 0; a < K; ++a)
                for (int c = 0; c < K; ++c) {
                    a2[a * K + c] = as_rmat7(pt.ai2[chart7::pair_index(K, a, c)]);
                    d2[a * K + c] = pt.d2ld[chart7::pair_index(K, a, c)];
                }
            Rcpp::List ch = shallow(Rcpp::List(ulo["chart"]));
            ch["eta"] = Rcpp::NumericVector(t.begin(), t.end());
            ch["h"] = Rcpp::NumericVector(0);
            ch["ai"] = as_rlist7(pt.ai);
            ch["ai2"] = a2;
            ch["dld"] = Rcpp::NumericVector(pt.dld.begin(), pt.dld.end());
            ch["d2ld"] = d2;
            ulo["chart"] = ch;
            const int p = u.pd;
            Dense B(m, m);
            for (int g = 0; g < u.G; ++g)
                for (int j = 0; j < p; ++j)
                    for (int i = 0; i < p; ++i) B(g * p + i, g * p + j) = pt.sinv.v[(std::size_t) j * p + i];
            add_block(S, u.idx, B);
        } else if (u.kind == 10) {
            // the multivariate t prior: the chart's free values and nu
            const int K = u.mst.ch.K;
            std::vector<double> eta(K, 0.0);
            double nu = u.nu0;
            for (std::size_t h = 0; h < u.mst.role.size(); ++h) {
                if (u.mst.role[h] >= 0) eta[u.mst.role[h]] = t[h];
                else nu = t[h];
            }
            chart7::ChartPt ch;
            if (!chart7::chart_order0(u.mst.ch, eta.data(), ch) || !ch.ok || !ch.L_ok) return false;
            chart7::MvtUnit0 u0;
            if (!chart7::mvt_order0(u.mst, ch, nu, u0)) return false;
            for (Rcpp::List* l : {&ul, &ulo}) {
                (*l)["L"] = as_rmat7(u0.L);
                (*l)["sinv"] = as_rmat7(u0.sinv);
                (*l)["logdet"] = u0.logdet;
                (*l)["nu"] = u0.nu;
                (*l)["c0"] = u0.c0;
                (*l)["half"] = u0.half;
            }
            std::vector<double> bu(m);
            for (int j = 0; j < m; ++j) bu[j] = b[u.idx[j]];
            std::vector<chart7::Dense> blocks;
            if (!chart7::mvt_moving_blocks(u.mst, ch, nu, bu.data(), m, blocks)) return false;
            const int p = u.mst.p;
            Dense B(m, m);
            for (std::size_t g = 0; g < blocks.size(); ++g)
                for (int j = 0; j < p; ++j)
                    for (int i = 0; i < p; ++i)
                        B((int) g * p + i, (int) g * p + j) = blocks[g].v[(std::size_t) j * p + i];
            add_block(S, u.idx, B);
            // the outer pieces at these coefficients (engine_outer_mvt())
            chart7::MvtGroups gr;
            chart7::MvtOuter1 o1;
            chart7::MvtOuter2 o2;
            if (chart7::chart_order1(u.mst.ch, ch) && chart7::chart_order2(u.mst.ch, ch) &&
                chart7::mvt_groups(u.mst, ch, nu, bu.data(), m, gr) &&
                chart7::mvt_outer_order1(u.mst, ch, gr, o1) &&
                chart7::mvt_outer_order2(u.mst, ch, gr, o2)) {
                const bool bq = chart7::mvt_beta_quadratic(u.mst, ch, nu);
                out.mvt[k] = Rcpp::List::create(
                    Rcpp::Named("H") = (int) u.mst.role.size(), Rcpp::Named("bq") = bq,
                    Rcpp::Named("gt") = Rcpp::NumericVector(o1.gt.begin(), o1.gt.end()), Rcpp::Named("rho2") = Rcpp::NumericVector(o2.rho2.begin(), o2.rho2.end()),
                    Rcpp::Named("c") = Rcpp::NumericVector(o1.c.begin(), o1.c.end()), Rcpp::Named("c2") = Rcpp::NumericVector(o2.c2.begin(), o2.c2.end()),
                    Rcpp::Named("dS") = Rcpp::NumericVector(o1.dS.begin(), o1.dS.end()), Rcpp::Named("d2S") = Rcpp::NumericVector(o2.d2S.begin(), o2.d2S.end()),
                    Rcpp::Named("ai") = as_rlist7(o1.ai), Rcpp::Named("Na") = Rcpp::NumericVector(o1.Na.begin(), o1.Na.end()));
            }
        } else if (u.kind == 4) {
            Chart c;
            if (!chart_at(u, t, c)) return false;
            for (Rcpp::List* l : {&ul, &ulo}) {
                (*l)["L"] = as_rmat(c.L);
                (*l)["sinv"] = as_rmat(c.sinv);
                (*l)["logdet"] = c.logdet;
            }
            Rcpp::List ch = shallow(Rcpp::List(ulo["chart"]));
            ch["eta"] = Rcpp::NumericVector(t.begin(), t.end());
            ch["h"] = Rcpp::NumericVector(c.h.begin(), c.h.end());
            ulo["chart"] = ch;
            const int p = u.pd;
            Dense B(m, m);
            for (int g = 0; g < u.G; ++g)
                for (int j = 0; j < p; ++j)
                    for (int i = 0; i < p; ++i) B(g * p + i, g * p + j) = c.sinv(i, j);
            add_block(S, u.idx, B);
        } else {
            return false;
        }
        units[k] = ul;
        unitso[k] = ulo;
    }
    for (double x : S.v) if (!R_FINITE(x)) return false;
    Rcpp::RObject Sx;
    if (ps.sparse) Sx = as_dgc(S); else Sx = as_rmat(S);
    Rcpp::RObject C;
    if (ps.need_C) {
        Engine tmp;
        Dense Cd;
        bool diag = true;
        for (int j = 0; j < ps.p; ++j)
            for (int i = 0; i < ps.p; ++i) if (i != j && S(i, j) != 0.0) diag = false;
        // an indefinite S (a t prior's Hessian away from its mode) has no
        // square root: penalty_sqrt() returns NULL
        if (tmp.penalty_sqrt_dense(S, Cd))
        {
            if (ps.sparse && diag) C = as_dgc(Cd); else C = as_rmat(Cd);
        }
    }
    out.pen = Rcpp::List::create(Rcpp::Named("S") = Sx, Rcpp::Named("C") = C,
                                 Rcpp::Named("units") = units);
    out.pen_outer = Rcpp::List::create(Rcpp::Named("S") = Sx, Rcpp::Named("C") = R_NilValue,
                                       Rcpp::Named("units") = unitso);
    out.ok = true;
    return true;
}

} // namespace pen7

#endif
