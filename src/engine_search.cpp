// The outer search of outer_fit() in compiled code (phase 4): the state of
// the search, evaluate() and derivs() around a point source, the searches
// through optimizers7's descent loop, and the exits (the box, the rounding
// floor, the lbfgs fallback, the restart with an empty memory, the edge
// probe). The point source is the PointSource class below; R's outer_fit()
// is the reference, and every step here repeats a step there, in its order.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <optimizers7/descent_loop.h>
#include <optimizers7/rules.h>
#include <R_ext/Lapack.h>
#include "engine_point.h"
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#ifndef FCONE
#define FCONE
#endif

namespace {

using optimizers7::Objective;

const double kInf = std::numeric_limits<double>::infinity();

// a point as the search keeps it (state$last)
struct SPoint {
    bool ok = false;
    double value = NA_REAL;
    bool have_grad = false, have_hess = false;
    std::vector<double> grad, hess;  // hess column-major nv x nv
    // the incumbent before this point, for derivs()'s no_grad() (D5)
    double pb_v = kInf;
    bool pb_has = false;
    std::vector<double> pb_beta;
};

// the point source in R (closures of outer_fit())
class RPointSource : public PointSource {
public:
    explicit RPointSource(Rcpp::List cb)
        : f_eval(cb["eval"]), f_grad(cb["grad"]), f_hess(cb["hess"]), f_read(cb["reading"]) {}
    void eval(const std::vector<double>& v, const std::vector<double>& warm, bool first,
              PointOut& out) {
        Rcpp::List r = f_eval(Rcpp::NumericVector(v.begin(), v.end()),
                              Rcpp::NumericVector(warm.begin(), warm.end()), first);
        out = PointOut();
        out.ok = Rcpp::as<bool>(r["ok"]);
        SEXP e = r["err"];
        out.has_err = !Rf_isNull(e);
        if (out.has_err) out.err = Rcpp::as<std::string>(e);
        if (out.ok) {
            out.value = Rcpp::as<double>(r["value"]);
            out.loglik = Rcpp::as<double>(r["loglik"]);
            out.penalty = Rcpp::as<double>(r["penalty"]);
            Rcpp::NumericVector p = r["par"], vv = r["vals"];
            out.par.assign(p.begin(), p.end());
            out.vals.assign(vv.begin(), vv.end());
        } else {
            out.hv = Rcpp::as<std::string>(r["hv"]);
            out.cause = Rcpp::as<std::string>(r["cause"]);
        }
    }
    bool grad(std::vector<double>& g) {
        Rcpp::RObject r = f_grad();
        if (Rf_isNull(r)) return false;
        Rcpp::NumericVector x(r);
        g.assign(x.begin(), x.end());
        return true;
    }
    bool hess(std::vector<double>& H) {
        Rcpp::RObject r = f_hess();
        if (Rf_isNull(r)) return false;
        Rcpp::NumericVector x(r);
        H.assign(x.begin(), x.end());
        return true;
    }
    void reading(double& r, bool& has_me, double& me) {
        Rcpp::List x = f_read();
        r = Rcpp::as<double>(x["r"]);
        SEXP m = x["mode_error"];
        has_me = !Rf_isNull(m);
        me = has_me ? Rcpp::as<double>(m) : NA_REAL;
    }
private:
    Rcpp::Function f_eval, f_grad, f_hess, f_read;
};

// the error statmod() reads at the start (R: stop(..., call. = FALSE))
struct StartError {
    std::string msg;
};

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

// the settings of one optimizer, as descent_settings() gives them
struct Settings {
    Rcpp::List method;
    optimizers7::LineSearchSpec ls;
    Rcpp::List rule;
    int maxit = 100, max_eval = INT_MAX, refresh = 10;
    bool verbose = false, keep_trace = false;
    double step = 1.0;
    std::string type;
};

Settings read_settings(Rcpp::List s) {
    Settings o;
    o.method = Rcpp::clone(Rcpp::List(s["method"]));
    o.ls = optimizers7::LineSearchSpec::from_list(s["line_search"]);
    o.rule = s["rule_spec"];
    o.maxit = Rcpp::as<int>(s["maxit"]);
    o.max_eval = Rcpp::as<int>(s["max_eval"]);
    o.verbose = Rcpp::as<bool>(s["verbose"]);
    o.refresh = Rcpp::as<int>(s["refresh"]);
    o.keep_trace = Rcpp::as<bool>(s["keep_trace"]);
    o.step = Rcpp::as<double>(s["step"]);
    o.type = Rcpp::as<std::string>(o.method["type"]);
    return o;
}

// crit_any(a, crit_abs_obj(tol)) as criterion_spec() describes it
Rcpp::List rule_any_abs(Rcpp::List a, double tol) {
    return Rcpp::List::create(
        Rcpp::Named("type") = "any",
        Rcpp::Named("criteria") = Rcpp::List::create(
            a, Rcpp::List::create(Rcpp::Named("type") = "abs_obj",
                                  Rcpp::Named("tol") = tol)));
}

// the result of one search, as outer_fit() reads it
struct SRes {
    std::vector<double> par;
    double value = NA_REAL;
    bool converged = false;
    int iterations = 0;
    std::string stopped_by, message;
    bool resolution_stop = false;  // criterion_met == "no decrease above ..."
};

class Search;

// fn, gr and he of outer_fit() as an optimizers7 Objective
class LoopObjective : public Objective {
public:
    LoopObjective(Search* s, bool g, bool h) : s_(s), g_(g), h_(h) {}
    double value(const arma::vec& x);
    arma::vec gradient(const arma::vec& x);
    arma::mat hessian(const arma::vec& x);
    bool has_gradient() const { return g_; }
    bool has_hessian() const { return h_; }
private:
    Search* s_;
    bool g_, h_;
};

class Search {
public:
    // configuration
    int nh = 0, ng = 0, nv = 0;
    bool exact = false, exact2 = false, chose = false, verbose = false, trace_grad = false;
    double sgn = -1.0;
    std::vector<int> gam_where;  // 0-based
    std::vector<char> angle;     // per hyperparameter row
    PointSource* src = nullptr;
    Rcpp::Function say_ok, say_na, say_me;

    // state
    std::vector<double> beta;
    double best_v = kInf;
    bool has_best = false;
    std::vector<double> best_beta;
    bool has_key = false;
    std::vector<double> key;
    SPoint last;
    bool has_last = false;
    int evals = 0;
    double worst = NA_REAL;
    std::vector<double> res_seen;
    double resolution = NA_REAL;
    bool said_mode_error = false;
    std::vector<std::vector<double>> rows;

    Search(Rcpp::List cb) : say_ok(cb["say_ok"]), say_na(cb["say_na"]), say_me(cb["say_me"]) {}

    static bool same(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) return false;
        for (std::size_t i = 0; i < a.size(); ++i) {
            const bool na = ISNAN(a[i]), nb = ISNAN(b[i]);
            if (na || nb) {
                // R's key: NaN and NaN are one key, NA and NaN two
                if (!(na && nb && R_IsNA(a[i]) == R_IsNA(b[i]))) return false;
                continue;
            }
            if (a[i] != b[i]) return false;
        }
        return true;
    }

    static double resolution_summary(const std::vector<double>& x) {
        bool any = false;
        double m = 0.0;
        for (double v : x) {
            if (R_FINITE(v) && v > 0.0) {
                if (!any || v > m) m = v;
                any = true;
            }
        }
        return any ? m : NA_REAL;
    }

    const SPoint& evaluate(const std::vector<double>& v) {
        if (has_key && same(key, v)) return last;
        std::vector<double> warm = has_best ? best_beta : beta;
        for (int j = 0; j < ng; ++j) warm[gam_where[j]] = v[nh + j];
        PointOut po;
        src->eval(v, warm, evals == 0, po);
        evals += 1;
        if (!po.ok) {
            if (verbose) {
                Rcpp::RObject hv_arg, err_arg;
                if (po.hv_num) hv_arg = Rcpp::NumericVector(po.vals.begin(), po.vals.end());
                else hv_arg = Rcpp::wrap(po.hv);
                if (po.has_err) err_arg = Rcpp::wrap(po.err);
                say_na(evals, hv_arg, err_arg);
            }
            if (evals == 1) {
                StartError e;
                e.msg = po.cause;
                throw e;
            }
            const double w = R_FINITE(worst) ? worst : 1.0;
            const double bar = w + std::fabs(w) + 1.0;
            SPoint out;
            out.value = sgn * bar;
            out.ok = false;
            out.have_grad = true;
            out.grad.assign(nv, 0.0);
            if (exact2) {
                out.have_hess = true;
                out.hess.assign((size_t) nv * nv, 0.0);
                for (int i = 0; i < nv; ++i) out.hess[(size_t) i * nv + i] = -1.0;
            }
            key = v;
            has_key = true;
            last = out;
            has_last = true;
            return last;
        }
        beta = po.par;
        const double vfn = sgn * po.value;
        SPoint out;
        out.pb_v = best_v;
        out.pb_has = has_best;
        out.pb_beta = best_beta;
        if (R_FINITE(vfn) && vfn < best_v) {
            best_v = vfn;
            best_beta = po.par;
            has_best = true;
        }
        if (R_FINITE(vfn)) {
            worst = R_FINITE(worst) ? std::max(worst, vfn) : vfn;
        }
        std::vector<double> row;
        row.push_back((double) evals);
        row.push_back(po.value);
        row.push_back(po.loglik);
        row.push_back(po.penalty);
        row.insert(row.end(), po.vals.begin(), po.vals.end());
        rows.push_back(row);
        out.value = po.value;
        out.ok = true;
        key = v;
        has_key = true;
        last = out;
        has_last = true;
        if (verbose) {
            // held protected until the call: a bare SEXP is collected
            Rcpp::RObject gv;
            if (exact && trace_grad) {
                try {
                    const std::vector<double>& g = derivs(1);
                    gv = Rcpp::NumericVector(g.begin(), g.end());
                } catch (StartError&) {
                    throw;
                } catch (...) {
                    gv = Rcpp::RObject();
                }
            }
            say_ok(evals, po.value, gv, Rcpp::NumericVector(po.vals.begin(), po.vals.end()));
        }
        if (chose) {
            double r = NA_REAL, me = NA_REAL;
            bool has_me = false;
            src->reading(r, has_me, me);
            if (R_FINITE(r) && r > 0.0) {
                res_seen.push_back(r);
                resolution = resolution_summary(res_seen);
            } else if (has_me) {
                if (verbose && !said_mode_error) {
                    said_mode_error = true;
                    say_me(me);
                }
            }
        }
        return last;
    }

    // derivs(order) on the criterion's scale, at the last point
    const std::vector<double>& derivs(int order) {
        SPoint& st = last;
        if (order == 1 && st.have_grad) return st.grad;
        if (order == 2 && st.have_hess) return st.hess;
        if (!st.ok) {
            if (order == 1) {
                st.grad.assign(nv, 0.0);
                st.have_grad = true;
                return st.grad;
            }
            if (!exact2) Rcpp::stop("no Hessian at an unavailable point");
            st.hess.assign((size_t) nv * nv, 0.0);
            for (int i = 0; i < nv; ++i) st.hess[(size_t) i * nv + i] = -1.0;
            st.have_hess = true;
            return st.hess;
        }
        if (order == 1) {
            std::vector<double> g;
            if (src->grad(g)) {
                st.grad = g;
            } else {
                // no_grad(): the incumbent the point displaced is put back (D5)
                best_v = st.pb_v;
                has_best = st.pb_has;
                best_beta = st.pb_beta;
                st.ok = false;
                st.grad.assign(nv, R_NaN);
            }
            st.have_grad = true;
            return st.grad;
        }
        std::vector<double> H;
        if (src->hess(H)) {
            st.hess = H;
        } else {
            st.hess.assign((size_t) nv * nv, 0.0);
            for (int i = 0; i < nv; ++i) st.hess[(size_t) i * nv + i] = -1.0;
        }
        st.have_hess = true;
        return st.hess;
    }

    double fn(const std::vector<double>& v) { return sgn * evaluate(v).value; }

    // the box (section 6 of the specification)
    bool has_box = false, use_box = false;

    // one minimize() of outer_fit()'s mini(); throws where minimize() errors
    SRes run(const Settings& s, const Rcpp::List& rule, std::vector<double> start,
             bool with_gr, bool with_he, bool resolution_fn, double resolution_const,
             bool have_resolution_const) {
        std::vector<optimizers7::Coord> coords;
        if (use_box) {
            // optimizer_bounded() is TRUE for every descent method
            for (int j = 0; j < nh; ++j) {
                if (angle[j] && std::fabs(start[j]) >= 9.0) {
                    start[j] = (start[j] > 0 ? 1.0 : (start[j] < 0 ? -1.0 : 0.0)) * 8.5;
                }
            }
            bool any = false;
            for (int j = 0; j < nh; ++j) if (angle[j]) any = true;
            if (any) {
                coords.resize(nv);
                for (int j = 0; j < nv; ++j) {
                    optimizers7::Coord c;
                    if (j < nh && angle[j]) {
                        c.lwr = -9.0;
                        c.upr = 9.0;
                        c.kind = optimizers7::BOUND_BOTH;
                        if (!(start[j] > -9.0 && start[j] < 9.0))
                            Rcpp::stop("The starting value must lie strictly inside its bounds.");
                    }
                    coords[j] = c;
                }
            }
        }
        for (double x : start) if (ISNAN(x)) Rcpp::stop("'par' must be a numeric vector of starting values.");
        optimizers7::LineSearchSpec ls = s.ls;
        if (resolution_fn) {
            ls.resolution_cb = [this]() { return resolution / 10.0; };
        } else if (have_resolution_const) {
            ls.resolution_cb = nullptr;
            ls.resolution_fn = Rcpp::RObject();
            ls.resolution = resolution_const;
        }
        LoopObjective obj(this, with_gr, with_he);
        std::unique_ptr<optimizers7::Direction> dir(optimizers7::make_direction(s.method, nv));
        optimizers7::CompiledStopRule crule(optimizers7::make_compiled_rule(rule));
        arma::vec par(start);
        optimizers7::DescentOut r = optimizers7::descent_loop(
            obj, par, crule, *dir, ls, s.maxit, s.max_eval, s.verbose, s.refresh,
            s.keep_trace, s.step, coords);
        SRes o;
        o.par.assign(r.par.begin(), r.par.end());
        o.value = r.value;
        o.converged = r.converged;
        o.iterations = r.iterations;
        o.stopped_by = r.stopped_by;
        o.message = r.message;
        o.resolution_stop = r.stopped_by == "resolution";
        return o;
    }

    // the decrement test of 8c, 8d and 8f; false where anything errors
    bool decrement(const std::vector<double>& at, double& dec) {
        try {
            evaluate(at);
            std::vector<double> g = derivs(1);
            std::vector<double> A = derivs(2);
            for (double& x : A) x = sgn * x;
            std::vector<double> As(A.size());
            for (int i = 0; i < nv; ++i)
                for (int j = 0; j < nv; ++j)
                    As[(size_t) j * nv + i] = (A[(size_t) j * nv + i] + A[(size_t) i * nv + j]) / 2.0;
            dec = joint_decrement(g, As);
            return R_FINITE(dec);
        } catch (StartError&) {
            throw;
        } catch (...) {
            return false;
        }
    }
};

double LoopObjective::value(const arma::vec& x) {
    ++n_value;
    return s_->fn(std::vector<double>(x.begin(), x.end()));
}
arma::vec LoopObjective::gradient(const arma::vec& x) {
    ++n_grad;
    s_->evaluate(std::vector<double>(x.begin(), x.end()));
    const std::vector<double>& g = s_->derivs(1);
    arma::vec out(g.size());
    for (std::size_t i = 0; i < g.size(); ++i) out[i] = s_->sgn * g[i];
    return out;
}
arma::mat LoopObjective::hessian(const arma::vec& x) {
    ++n_hess;
    s_->evaluate(std::vector<double>(x.begin(), x.end()));
    const std::vector<double>& H = s_->derivs(2);
    arma::mat out(s_->nv, s_->nv);
    for (std::size_t i = 0; i < H.size(); ++i) out[i] = s_->sgn * H[i];
    return out;
}

Rcpp::List sres_list(const SRes& r) {
    return Rcpp::List::create(
        Rcpp::Named("par") = Rcpp::NumericVector(r.par.begin(), r.par.end()),
        Rcpp::Named("value") = r.value, Rcpp::Named("converged") = r.converged,
        Rcpp::Named("iterations") = r.iterations, Rcpp::Named("stopped_by") = r.stopped_by,
        Rcpp::Named("message") = r.message, Rcpp::Named("resolution_stop") = r.resolution_stop);
}

} // namespace

// The search of outer_fit() from its first evaluation to the edge probe
// (sections 7 and 8 of SPEC_outerloop.md). `cfg`: nh, ng, exact, exact2,
// chose, verbose, trace_grad, gam_where (1-based), angle, beta, eta0, own (the
// main optimizer's settings), alt (lbfgs()'s settings), backtracks. `cb`: the
// point source's closures and the trace's printers. Returns the search's
// result and state, or list(error = message) where outer_fit() stops.
// [[Rcpp::export]]
Rcpp::List engine_outer_search(Rcpp::List cfg, Rcpp::List cb) {
    RPointSource rsrc(cb);
    Search S(cb);
    S.src = &rsrc;
    // the point in compiled code where outer_fit() describes it (stage 2)
    std::unique_ptr<PointSource> csrc;
    if (cfg.containsElementNamed("compiled") && !Rf_isNull(cfg["compiled"])) {
        csrc.reset(make_compiled_point_source(cfg["compiled"], cb, &rsrc));
        S.src = csrc.get();
    }
    S.nh = Rcpp::as<int>(cfg["nh"]);
    S.ng = Rcpp::as<int>(cfg["ng"]);
    S.nv = S.nh + S.ng;
    S.exact = Rcpp::as<bool>(cfg["exact"]);
    S.exact2 = Rcpp::as<bool>(cfg["exact2"]);
    S.chose = Rcpp::as<bool>(cfg["chose"]);
    S.verbose = Rcpp::as<bool>(cfg["verbose"]);
    S.trace_grad = Rcpp::as<bool>(cfg["trace_grad"]);
    {
        Rcpp::IntegerVector gw = cfg["gam_where"];
        for (int j : gw) S.gam_where.push_back(j - 1);
        Rcpp::LogicalVector an = cfg["angle"];
        for (int j = 0; j < an.size(); ++j) S.angle.push_back(an[j] == TRUE);
        Rcpp::NumericVector b = cfg["beta"];
        S.beta.assign(b.begin(), b.end());
    }
    Rcpp::NumericVector e0 = cfg["eta0"];
    const std::vector<double> eta0(e0.begin(), e0.end());
    Settings opt = read_settings(cfg["own"]);
    const Rcpp::List own_rule = opt.rule;
    const bool chose = S.chose, exact = S.exact, exact2 = S.exact2;
    const int nh = S.nh, ng = S.ng;
    const double mel = 1e-3;  // mode_error_limit()
    const double cert_tol = Rcpp::as<double>(cfg["cert_tol"]);
    bool has_ls = opt.type == "newton" || opt.type == "lbfgs" || opt.type == "bfgs" ||
                  opt.type == "gd" || opt.type == "cg" || opt.type == "bb";
    for (int j = 0; j < nh; ++j) if (S.angle[j]) S.has_box = chose;

    double r1 = NA_REAL;
    std::vector<double> typical;
    Rcpp::List rule = own_rule;
    // which optimizer the returned one is: 0 the main one, 1 lbfgs()
    int which_opt = 0;
    try {
        if (chose) {
            S.fn(eta0);
            r1 = S.resolution;
            if (R_FINITE(r1) && r1 > 0.0) rule = rule_any_abs(own_rule, r1);
            opt.ls.max_step = Rcpp::as<int>(cfg["backtracks"]);
            // the typical size (D2)
            Rcpp::NumericVector ts = opt.method.containsElementNamed("typical_size")
                ? Rcpp::NumericVector(opt.method["typical_size"]) : Rcpp::NumericVector(1);
            if (ng && exact2 && opt.method.containsElementNamed("typical_size") &&
                ts.size() == 0) {
                bool okh = true;
                std::vector<double> H0;
                try {
                    S.evaluate(eta0);
                    H0 = S.derivs(2);
                } catch (StartError&) {
                    throw;
                } catch (...) {
                    okh = false;
                }
                if (okh) {
                    std::vector<double> se0(ng);
                    bool any = false;
                    for (int j = 0; j < ng; ++j) {
                        const int i = nh + j;
                        double s = 1.0 / std::sqrt(std::fabs(S.sgn * H0[(size_t) i * S.nv + i]));
                        if (!R_FINITE(s)) s = 1.0;
                        se0[j] = s;
                        if (s > 1.0) any = true;
                    }
                    if (any) {
                        typical.assign(nh, 1.0);
                        for (int j = 0; j < ng; ++j) typical.push_back(std::max(1.0, se0[j]));
                        opt.method["typical_size"] = Rcpp::NumericVector(typical.begin(), typical.end());
                    }
                }
            }
        }
        auto mini = [&](const Settings& s, const Rcpp::List& rl, const std::vector<double>& st,
                        bool gr, bool he, bool own_res, double rc, bool has_rc) {
            return S.run(s, rl, st, gr, he, own_res, rc, has_rc);
        };
        const bool res_fn = chose && has_ls;
        SRes res = mini(opt, rule, eta0, exact, exact2, res_fn, 0.0, false);

        // the box (6)
        if (S.has_box && !S.use_box) {
            bool out = false;
            for (int j = 0; j < nh; ++j)
                if (S.angle[j] && std::fabs(res.par[j]) > 9.0) out = true;
            if (out) {
                S.use_box = true;
                try {
                    SRes rb = mini(opt, rule, res.par, exact, exact2, res_fn, 0.0, false);
                    if (R_FINITE(rb.value)) res = rb;
                } catch (StartError&) {
                    throw;
                } catch (...) {}
            }
        }
        // the rounding floor (8c)
        if (chose && exact2 && !res.converged &&
            res.message == "the line search found no acceptable step") {
            const double floor_res = 1e3 * DBL_EPSILON *
                std::max(1.0, R_FINITE(res.value) ? std::fabs(res.value) : 1.0);
            double dec = NA_REAL;
            if (S.decrement(res.par, dec) && dec <= floor_res) {
                res.converged = true;
                res.resolution_stop = true;
                res.message = "";
                res.stopped_by = "resolution";
            }
        }
        // settled (8d)
        bool settled = res.converged;
        if (!settled && exact2) {
            double dec = NA_REAL;
            settled = S.decrement(res.par, dec) && dec <= cert_tol;
        }
        // the lbfgs fallback (8e)
        Settings cur = opt;
        Rcpp::List cur_rule = rule;
        bool cur_res_fn = res_fn;
        if (chose && !settled && exact && opt.type == "newton") {
            Settings alt = read_settings(cfg["alt"]);
            bool ok2 = true;
            SRes r2;
            try {
                r2 = mini(alt, alt.rule, eta0, true, false, false, 0.0, false);
            } catch (StartError&) {
                throw;
            } catch (...) {
                ok2 = false;
            }
            if (ok2 && R_FINITE(r2.value) &&
                r2.value < res.value - 1e-8 * std::max(1.0, std::fabs(res.value))) {
                res = r2;
                cur = alt;
                cur_rule = alt.rule;
                cur_res_fn = false;
                which_opt = 1;
                settled = res.converged;
                if (!settled && exact2) {
                    double dec = NA_REAL;
                    settled = S.decrement(res.par, dec) && dec <= cert_tol;
                }
            }
        }
        // the restart with an empty memory (8f)
        Settings rst = cur;
        Rcpp::List rst_rule = cur_rule;
        const bool rst_has_ls = chose && (cur.type == "newton" || cur.type == "lbfgs" ||
                                          cur.type == "bfgs" || cur.type == "gd" ||
                                          cur.type == "cg" || cur.type == "bb");
        if (rst_has_ls) rst_rule = rule_any_abs(own_rule, mel);
        for (int attempt = 1; attempt <= 3; ++attempt) {
            if (!(chose && exact2 && res.converged &&
                  (attempt > 1 || res.resolution_stop))) break;
            double dec = NA_REAL;
            const bool rising = S.decrement(res.par, dec) && dec > mel;
            if (!rising) break;
            bool ok2 = true;
            SRes r2;
            try {
                r2 = mini(rst, rst_rule, res.par, true, true, false, mel / 10.0, rst_has_ls);
            } catch (StartError&) {
                throw;
            } catch (...) {
                ok2 = false;
            }
            if (!ok2 || !R_FINITE(r2.value) || !(r2.value < res.value - mel)) {
                try { S.evaluate(res.par); } catch (StartError&) { throw; } catch (...) {}
                break;
            }
            res = r2;
            settled = res.converged;
        }
        // the edge probe (8g)
        if (chose && nh) {
            for (int attempt = 1; attempt <= 3; ++attempt) {
                const double k_best_v = S.best_v;
                const bool k_has_best = S.has_best;
                const std::vector<double> k_best_beta = S.best_beta, k_beta = S.beta,
                                          k_res_seen = S.res_seen;
                const double k_resolution = S.resolution;
                auto restore = [&]() {
                    S.best_v = k_best_v;
                    S.has_best = k_has_best;
                    S.best_beta = k_best_beta;
                    S.beta = k_beta;
                    S.res_seen = k_res_seen;
                    S.resolution = k_resolution;
                    try { S.evaluate(res.par); } catch (StartError&) { throw; } catch (...) {}
                };
                // edge_probe_hyper(res@par, eta0, nh, res@value, fn)
                bool found = false;
                std::vector<double> best;
                {
                    std::vector<int> at_edge;
                    for (int j = 0; j < nh; ++j) if (std::fabs(res.par[j]) > 8.0) at_edge.push_back(j);
                    if (!at_edge.empty() && R_FINITE(res.value)) {
                        double bv = res.value - mel;
                        const double fr[4] = {0.75, 0.5, 0.25, 0.0};
                        for (int j : at_edge) {
                            const double s = std::fabs(eta0[j]) <= 8.0 ? eta0[j] : 0.0;
                            for (double f : fr) {
                                std::vector<double> v = res.par;
                                v[j] = s + f * (res.par[j] - s);
                                double fv = NA_REAL;
                                try { fv = S.fn(v); } catch (StartError&) { throw; } catch (...) { fv = NA_REAL; }
                                if (R_FINITE(fv) && fv < bv) {
                                    bv = fv;
                                    best = v;
                                    found = true;
                                }
                            }
                        }
                    }
                }
                if (!found) {
                    restore();
                    break;
                }
                bool ok2 = true;
                SRes r2;
                try {
                    r2 = mini(cur, cur_rule, best, exact, exact2, cur_res_fn, 0.0, false);
                } catch (StartError&) {
                    throw;
                } catch (...) {
                    ok2 = false;
                }
                if (!ok2 || !R_FINITE(r2.value) || !(r2.value < res.value - mel)) {
                    restore();
                    break;
                }
                res = r2;
                settled = res.converged;
                if (!settled && exact2) {
                    double dec = NA_REAL;
                    settled = S.decrement(res.par, dec) && dec <= cert_tol;
                }
            }
        }
        // the final refit at the reported point, from the incumbent (D40)
        Rcpp::List refit;
        if (Rcpp::as<bool>(cfg["refit"])) {
            std::vector<double> warm = S.has_best ? S.best_beta : S.beta;
            for (int j = 0; j < ng; ++j) warm[S.gam_where[j]] = res.par[nh + j];
            refit = S.src->refit(res.par, warm);
        }
        // the rows of hist_outer, one numeric vector per row
        Rcpp::List rows(S.rows.size());
        for (std::size_t i = 0; i < S.rows.size(); ++i)
            rows[i] = Rcpp::NumericVector(S.rows[i].begin(), S.rows[i].end());
        return Rcpp::List::create(
            Rcpp::Named("res") = sres_list(res), Rcpp::Named("settled") = settled,
            Rcpp::Named("which_opt") = which_opt, Rcpp::Named("r1") = r1,
            Rcpp::Named("resolution") = S.resolution,
            Rcpp::Named("typical") = Rcpp::NumericVector(typical.begin(), typical.end()),
            Rcpp::Named("evals") = S.evals, Rcpp::Named("rows") = rows,
            Rcpp::Named("handed_back") = S.src->handed_back(),
            Rcpp::Named("refit") = refit,
            Rcpp::Named("beta") = Rcpp::NumericVector(S.beta.begin(), S.beta.end()),
            Rcpp::Named("best_beta") = S.has_best
                ? (SEXP) Rcpp::NumericVector(S.best_beta.begin(), S.best_beta.end())
                : R_NilValue);
    } catch (StartError& e) {
        return Rcpp::List::create(Rcpp::Named("error") = e.msg);
    }
}
