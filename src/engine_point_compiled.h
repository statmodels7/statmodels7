// ---------------------------------------------------------------------------
// The compiled point of the outer search (phase 4, stage 2): evaluate()'s
// body for a model whose inner fit is one compiled iwls pass over the whole
// vector, with the outer context of phase 3. The penalty at the point comes
// from R (`pen_at`) until stage 3. A point that the compiled route does not
// compute (the inner fit or the context hands back, the first evaluation
// unavailable, an error whose message R writes) is handed to `fallback`, the
// R closures of stage 1, which compute it as outer_fit() does; a derivative
// or a reading that the context hands back is read from R at that point (D47).
// Included at the end of engine_iwls.cpp, once.

#include "engine_point.h"
#include <memory>

namespace {

std::vector<int> ivec0(SEXP x) {
    std::vector<int> o;
    if (Rf_isNull(x)) return o;
    Rcpp::IntegerVector v(x);
    for (int j : v) o.push_back(j - 1);
    return o;
}
std::vector<double> dvec(SEXP x) {
    Rcpp::NumericVector v(x);
    return std::vector<double>(v.begin(), v.end());
}

class CompiledPointSource : public PointSource {
public:
    CompiledPointSource(Rcpp::List cfg, Rcpp::List cb, PointSource* fb)
        : model(cfg["model"]), method(cfg["method"]), fallback(fb),
          pen_at(cb["pen_at"]), r_crit(cb["r_crit"]), r_grad(cb["r_grad"]),
          r_hess(cb["r_hess"]), r_mvt(cb["r_mvt"]), r_set(cb["r_set"]),
          r_hold(cb["r_hold"]), pen_outer_at(cb["pen_outer_at"]) {
        // the penalty in compiled code (stage 3), where every unit is covered
        if (cfg.containsElementNamed("pen_spec") && !Rf_isNull(cfg["pen_spec"])) {
            ps.reset(new pen7::PenSpec(pen7::read_pen_spec(cfg["pen_spec"])));
            Rcpp::List hl = cfg["held"];
            hw0 = ivec0(hl["where"]);
            Rcpp::IntegerVector src = hl["src"];
            hsrc.assign(src.begin(), src.end());
            hval = dvec(hl["value"]);
            read_hold();
        }
        // the compiled alternation beside a kinked block (statmod_alternate()
        // through engine_alternate_fit())
        if (cfg.containsElementNamed("kink") && !Rf_isNull(cfg["kink"])) {
            kink = cfg["kink"];
            has_kink = true;
            kc = ivec0(kink["kc"]);
            Rcpp::List kb = kink["kblocks"];
            for (int b = 0; b < kb.size(); ++b) {
                Rcpp::List bl = kb[b];
                kb_unit.push_back(Rcpp::as<int>(bl["unit"]) - 1);
                kb_todo.push_back(Rcpp::as<bool>(bl["curv_todo"]));
            }
        }
        Rcpp::List g = cfg["groups"];
        for (int k = 0; k < g.size(); ++k) grp.push_back(ivec0(g[k]));
        gam = ivec0(cfg["gam_where"]);
        held_free = ivec0(cfg["held_free"]);
        nh = Rcpp::as<int>(cfg["nh"]);
        SEXP b = cfg["basis"];
        if (!Rf_isNull(b)) {
            Rcpp::NumericMatrix Bm(b);
            ml = true;
            B = Dense(Bm.nrow(), Bm.ncol());
            std::copy(Bm.begin(), Bm.end(), B.v.begin());
        }
    }

    int handed_back() const { return n_back; }

    // outer_fit()'s final refit, statmod_alternate() and statmod_marginal()
    // at v from warm, with no context carried over. Returns list(kind, raw,
    // value0, value, crit, crit_na) where raw is the compiled inner fit's
    // own list (engine_iwls_fit() or engine_alternate()), or an empty list
    // where a step goes to R.
    Rcpp::List refit(const std::vector<double>& v, const std::vector<double>& warm) {
        last_r = false;
        O.reset();
        cur_v = v;
        r_ready = false;
        Rcpp::List pen;
        std::vector<int> hw;
        std::vector<double> hv;
        if (ps) {
            hw = hw0;
            hv.resize(hw.size());
            for (std::size_t j = 0; j < hw.size(); ++j)
                hv[j] = hsrc[j] > 0 ? v.at(hsrc[j] - 1) : (hsrc[j] < 0 ? warm.at(hw[j]) : hval.at(j));
            std::vector<double> b = warm;
            for (std::size_t j = 0; j < hw.size(); ++j) b.at(hw[j]) = hv[j];
            pen7::PenOut po;
            if (!pen7::build_pen(*ps, v, b, po)) return Rcpp::List();
            pen = po.pen;
            pen_outer = po.pen_outer;
        } else {
            Rcpp::List pa = pen_at(Rcpp::NumericVector(v.begin(), v.end()),
                                   Rcpp::NumericVector(warm.begin(), warm.end()));
            r_ready = true;
            if (Rf_isNull(pa["pen"]) || Rf_isNull(pa["pen_outer"])) return Rcpp::List();
            pen = pa["pen"];
            pen_outer = pa["pen_outer"];
            hold = ivec0(pa["hold"]);
            hw = ivec0(pa["held_where"]);
            hv = dvec(pa["held_value"]);
        }
        Rcpp::List raw;
        std::vector<double> par;
        double value0 = NA_REAL, value = NA_REAL;
        std::vector<std::vector<double>> cvs;
        if (has_kink) {
            Rcpp::IntegerVector hw1(hw.size()), pin1(hold.size());
            for (std::size_t j = 0; j < hw.size(); ++j) hw1[j] = hw[j] + 1;
            for (std::size_t j = 0; j < hold.size(); ++j) pin1[j] = hold[j] + 1;
            raw = engine_alternate(
                model, pen, method, Rcpp::NumericVector(warm.begin(), warm.end()),
                kink["kblocks"], kink["smooth"], kink["sgroups"], kink["sfrozen"], hw1,
                Rcpp::NumericVector(hv.begin(), hv.end()), kink["kinked"], pin1,
                Rcpp::as<bool>(kink["expected"]), Rcpp::as<int>(kink["maxit"]),
                Rcpp::as<double>(kink["tol"]));
            if (!(Rcpp::as<int>(raw["status"]) == 0)) return Rcpp::List();
            par = dvec(raw["par"]);
            Rcpp::List cvl = raw["curvature"];
            cv_r = Rcpp::List(cvl.size());
            cvs.assign(ps->units.size(), std::vector<double>());
            for (int b = 0; b < cvl.size(); ++b) {
                if (!kb_todo[b]) { cv_r[b] = R_NilValue; continue; }
                cv_r[b] = cvl[b];
                cvs[kb_unit[b]] = dvec(cvl[b]);
            }
        } else {
            Engine E;
            try {
                engine_setup(E, model, pen);
                read_method(E, method);
                // statmod_alternate()'s objective at its start, before the
                // holds are written in
                if (!E.fn(warm, value0)) return Rcpp::List();
                std::vector<double> start = warm;
                for (std::size_t j = 0; j < hw.size(); ++j) start.at(hw[j]) = hv[j];
                std::vector<char> froz(E.M.p, 0);
                for (int j : hw) froz[j] = 1;
                IwlsOut o;
                for (int k = 0; k < T_N; ++k) g_time[k] = 0.0;
                if (!iwls_loop(E, start, froz, grp, false, o)) return Rcpp::List();
                raw = loop_result(o, E);
                par = o.beta;
                // fit_smooth()'s obj$fn(out)
                if (!E.fn(par, value)) return Rcpp::List();
            } catch (RError&) {
                return Rcpp::List();
            }
        }
        mode = par;
        // statmod_marginal() with no context: a fresh one at the mode
        double crit = NA_REAL;
        bool crit_na = false, crit_r = false;
        std::vector<int> zk;
        for (int j : kc) if (mode[j] == 0) zk.push_back(j);
        if (moving_pen(pen_outer)) {
            if (ps) {
                pen7::PenOut po;
                if (!pen7::build_pen(*ps, v, mode, po, has_kink ? &cvs : nullptr)) return Rcpp::List();
                pen_outer = po.pen_outer;
            } else {
                ready();
                Rcpp::RObject r = pen_outer_at(Rcpp::NumericVector(mode.begin(), mode.end()));
                if (Rf_isNull(r)) return Rcpp::List();
                pen_outer = r;
            }
        }
        {
            Outer O2;
            O2.model = model;
            O2.pen = pen_outer;
            bool gb = false;
            try {
                engine_setup(O2.E, model, pen_outer);
                O2.hold = has_kink ? lap : hold;
                if (has_kink) {
                    O2.two = true;
                    O2.hold_mode = hold;
                    O2.hold_mode.insert(O2.hold_mode.end(), zk.begin(), zk.end());
                }
                if (ml) { O2.ml = true; O2.B = B; }
                O2.setup(mode, false);
            } catch (RError&) {
                gb = true;
            }
            if (gb || O2.status) crit_r = true;
            else if (O2.crit_na) crit_na = true;
            else crit = O2.ll - O2.rho + O2.q / 2.0 * std::log(2.0 * M_PI) - O2.logdet_c / 2.0;
        }
        return Rcpp::List::create(
            Rcpp::Named("kind") = has_kink ? "kink" : "iwls", Rcpp::Named("raw") = raw,
            Rcpp::Named("value0") = value0, Rcpp::Named("value") = value,
            Rcpp::Named("crit") = crit, Rcpp::Named("crit_na") = crit_na,
            Rcpp::Named("crit_r") = crit_r, Rcpp::Named("cv") = cv_r);
    }

    void eval(const std::vector<double>& v, const std::vector<double>& warm, bool first,
              PointOut& out) {
        last_r = false;
        O.reset();
        if (!compiled_eval(v, warm, out) || (first && !out.ok)) {
            hand_back(v, warm, first, out);
        }
    }

    bool grad(std::vector<double>& g) {
        if (last_r) return fallback->grad(g);
        if (!complete()) return r_deriv(r_grad, g);
        if (!O->inv_ok) return false;
        if (!O->gradient(mrow, munit, mname, mchain, nh, gam, true, g)) return r_deriv(r_grad, g);
        return true;
    }
    bool hess(std::vector<double>& H) {
        if (last_r) return fallback->hess(H);
        if (!complete()) return r_deriv(r_hess, H);
        if (!O->inv_ok) return false;
        Dense D;
        if (!O->hessian(mrow, munit, mname, mchain, nh, gam, h1, h2, D)) return r_deriv(r_hess, H);
        H = D.v;
        return true;
    }

    // criterion_resolution() at the last point
    void reading(double& r, bool& has_me, double& me) {
        if (last_r) { fallback->reading(r, has_me, me); return; }
        r = NA_REAL;
        has_me = false;
        me = NA_REAL;
        if (!sc_ok) return;
        for (double x : sc) if (!R_FINITE(x)) return;
        // engine_ctx_solve(); where the context has no inverse the R route's
        // ctx_penalized() has none either
        if (O->status || !O->ensure_inverse()) return;
        const int p = O->p;
        std::vector<double> db(p);
        dense_matvec(O->Kinv, sc.data(), db.data());
        for (double x : db) if (!R_FINITE(x)) return;
        std::vector<double> pr(p);
        for (int i = 0; i < p; ++i) pr[i] = sc[i] * db[i];
        const double quad = 0.5 * r_sum(pr.data(), p);
        if (!R_FINITE(quad) || quad > 1e-3) {
            has_me = true;
            me = quad;
            return;
        }
        std::vector<double> par2(p);
        for (int i = 0; i < p; ++i) par2[i] = mode[i] - db[i];
        double m2 = NA_REAL;
        bool have = false;
        {
            // the criterion at par2 reads the penalty there, as R builds a
            // fresh context at the displaced coefficients
            Rcpp::List pen2 = pen_outer;
            bool gb = false;
            if (moving_pen(pen_outer)) {
                if (ps) {
                    pen7::PenOut po;
                    if (pen7::build_pen(*ps, cur_v, par2, po, has_kink ? &cvs_last : nullptr))
                        pen2 = po.pen_outer;
                    else gb = true;
                } else {
                    ready();
                    Rcpp::RObject x = pen_outer_at(Rcpp::NumericVector(par2.begin(), par2.end()));
                    if (Rf_isNull(x)) gb = true; else pen2 = x;
                }
            }
            Outer O2;
            O2.model = model;
            O2.pen = pen2;
            if (!gb) try {
                engine_setup(O2.E, model, pen2);
                O2.hold = has_kink ? lap : hold;
                if (has_kink) {
                    O2.two = true;
                    O2.hold_mode = hold;
                    for (int j : kc) if (par2[j] == 0) O2.hold_mode.push_back(j);
                }
                if (ml) { O2.ml = true; O2.B = B; }
                O2.setup(par2, false);
            } catch (RError&) {
                gb = true;
            }
            if (gb || O2.status) {
                ready();
                Rcpp::RObject x = r_crit(Rcpp::NumericVector(par2.begin(), par2.end()));
                if (!Rf_isNull(x)) { m2 = Rcpp::as<double>(x); have = true; }
            } else if (!O2.crit_na) {
                m2 = O2.ll - O2.rho + O2.q / 2.0 * std::log(2.0 * M_PI) - O2.logdet_c / 2.0;
                have = true;
            }
        }
        if (!have || !R_FINITE(m2)) return;
        const double d = std::fabs(m2 - value);
        if (R_FINITE(d) && d > 0) r = d;
    }

private:
    Rcpp::List model, method;
    PointSource* fallback;
    Rcpp::Function pen_at, r_crit, r_grad, r_hess, r_mvt, r_set, r_hold, pen_outer_at;

    static bool moving_pen(const Rcpp::List& pen) {
        Rcpp::List us = pen["units"];
        for (int k = 0; k < us.size(); ++k) {
            const int kd = Rcpp::as<int>(Rcpp::List(us[k])["kind"]);
            if (kd >= 5 && kd <= 10) return true;
        }
        return false;
    }
    std::unique_ptr<pen7::PenSpec> ps;
    std::vector<int> hw0, hsrc, pinned;
    std::vector<double> hval, cur_v;
    bool r_ready = false;

    // pinned_coords() and attr(design, "pinned"), which the first evaluation
    // may set (its R branch)
    void read_hold() {
        Rcpp::List h = r_hold();
        hold = ivec0(h["hold"]);
        pinned = ivec0(h["pinned"]);
        lap = ivec0(h["lap"]);
    }
    Rcpp::List kink;
    bool has_kink = false;
    std::vector<int> kc, lap, kb_unit;
    std::vector<char> kb_todo;
    Rcpp::List cv_r;   // the curvatures the alternation settled, per block
    std::vector<std::vector<double>> cvs_last;  // the same, per unit
    std::vector<Rcpp::RObject> mvt_last;        // kind 10's outer pieces at the mode
    // the R closures' point (pt$sp, pt$hy) for a fallback at this point
    void ready() {
        if (!r_ready) {
            r_set(Rcpp::NumericVector(cur_v.begin(), cur_v.end()), cv_r);
            r_ready = true;
        }
    }
    std::vector<std::vector<int>> grp;
    std::vector<int> gam, held_free;
    int nh = 0;
    bool ml = false;
    Dense B;
    int n_back = 0;
    // the last point
    bool last_r = false;
    std::unique_ptr<Outer> O;
    Rcpp::List pen_outer;
    std::vector<int> hold, mrow, munit, mname;
    std::vector<double> mchain, h1, h2, mode, sc;
    bool have_mem = false, sc_ok = false, mvt_done = false;
    double value = NA_REAL;

    void hand_back(const std::vector<double>& v, const std::vector<double>& warm, bool first,
                   PointOut& out) {
        ++n_back;
        last_r = true;
        O.reset();
        fallback->eval(v, warm, first, out);
        if (ps) read_hold();
    }

    template <class T>
    bool r_deriv(Rcpp::Function& f, T& x) {
        ready();
        Rcpp::RObject r = f(Rcpp::NumericVector(mode.begin(), mode.end()));
        if (Rf_isNull(r)) return false;
        Rcpp::NumericVector y(r);
        x.assign(y.begin(), y.end());
        return true;
    }

    // engine_outer_complete() and engine_outer_members()
    bool complete() {
        if (!have_mem || O->status) return false;
        if (!mvt_done) {
            std::vector<int> need;
            for (std::size_t k = 0; k < O->E.M.units.size(); ++k) {
                const PenUnit& u = O->E.M.units[k];
                if (u.kind == 10 && !u.mo) need.push_back((int) k);
            }
            if (!need.empty() && ps) {
                // built in C++ with the penalty at the mode
                for (int k : need) {
                    if ((std::size_t) k >= mvt_last.size() || mvt_last[k].isNULL()) return false;
                    read_mvt_outer(O->E.M.units[k], Rcpp::List(mvt_last[k]));
                }
                O->have_uder = false;
            } else if (!need.empty()) {
                Rcpp::IntegerVector nk(need.size());
                for (std::size_t j = 0; j < need.size(); ++j) nk[j] = need[j] + 1;
                ready();
                Rcpp::RObject r = r_mvt(Rcpp::NumericVector(mode.begin(), mode.end()), nk);
                if (Rf_isNull(r)) return false;
                Rcpp::List obs(r);
                for (std::size_t j = 0; j < need.size(); ++j)
                    read_mvt_outer(O->E.M.units[need[j]], obs[j]);
                O->have_uder = false;
            }
            mvt_done = true;
        }
        return true;
    }

    // false where the point goes to the R closures
    bool compiled_eval(const std::vector<double>& v, const std::vector<double>& warm,
                       PointOut& out) {
        out = PointOut();
        have_mem = sc_ok = mvt_done = false;
        mvt_last.clear();
        cur_v = v;
        r_ready = false;
        Rcpp::List pen;
        std::vector<int> hw;
        std::vector<double> hv;
        if (ps) {
            // the held values: the criterion's coefficients from v, the
            // others as they were (held_positions())
            hw = hw0;
            hv.resize(hw.size());
            for (std::size_t j = 0; j < hw.size(); ++j)
                hv[j] = hsrc[j] > 0 ? v.at(hsrc[j] - 1) : (hsrc[j] < 0 ? warm.at(hw[j]) : hval.at(j));
            std::vector<double> b = warm;
            for (std::size_t j = 0; j < hw.size(); ++j) b.at(hw[j]) = hv[j];
            pen7::PenOut po;
            if (!pen7::build_pen(*ps, v, b, po)) return false;
            pen = po.pen;
            pen_outer = po.pen_outer;
            mrow = ps->mrow;
            munit = ps->munit;
            mname = ps->mname;
            mchain = po.mchain;
            h1 = po.h1;
            h2 = po.h2;
            have_mem = true;
            out.vals = po.vals;
            for (std::size_t j = nh; j < v.size(); ++j) out.vals.push_back(v[j]);
            out.hv_num = true;
        } else {
        Rcpp::List pa = pen_at(Rcpp::NumericVector(v.begin(), v.end()),
                               Rcpp::NumericVector(warm.begin(), warm.end()));
        r_ready = true;
        if (Rf_isNull(pa["pen"]) || Rf_isNull(pa["pen_outer"])) return false;
        pen = pa["pen"];
        pen_outer = pa["pen_outer"];
        hold = ivec0(pa["hold"]);
        pinned = ivec0(pa["pinned"]);
        out.vals = dvec(pa["vals"]);
        out.hv = Rcpp::as<std::string>(pa["hv"]);
        hw = ivec0(pa["held_where"]);
        hv = dvec(pa["held_value"]);
        SEXP mem = pa["mem"];
        if (!Rf_isNull(mem)) {
            Rcpp::List m(mem);
            mrow = ivec0(m["mrow"]);
            munit = ivec0(m["munit"]);
            mname = ivec0(m["mname"]);
            mchain = dvec(m["mchain"]);
            h1 = dvec(m["h1"]);
            h2 = dvec(m["h2"]);
            have_mem = true;
        }
        }
        // fit_smooth(): the held values written in, the compiled pass
        std::vector<double> start = warm;
        for (std::size_t j = 0; j < hw.size(); ++j) start.at(hw[j]) = hv[j];
        IwlsOut o;
        std::vector<std::vector<double>> cvs;
        cv_r = Rcpp::List();
        if (has_kink) {
            // engine_alternate_fit(): from the warm start, the holds passed
            Rcpp::IntegerVector hw1(hw.size()), pin1(hold.size());
            for (std::size_t j = 0; j < hw.size(); ++j) hw1[j] = hw[j] + 1;
            for (std::size_t j = 0; j < hold.size(); ++j) pin1[j] = hold[j] + 1;
            Rcpp::List r = engine_alternate(
                model, pen, method, Rcpp::NumericVector(warm.begin(), warm.end()),
                kink["kblocks"], kink["smooth"], kink["sgroups"], kink["sfrozen"], hw1,
                Rcpp::NumericVector(hv.begin(), hv.end()), kink["kinked"], pin1,
                Rcpp::as<bool>(kink["expected"]), Rcpp::as<int>(kink["maxit"]),
                Rcpp::as<double>(kink["tol"]));
            if (!(Rcpp::as<int>(r["status"]) == 0)) return false;
            o.beta = dvec(r["par"]);
            o.converged = Rcpp::as<bool>(r["converged"]);
            // the curvatures written into the specification (curv_todo)
            Rcpp::List cvl = r["curvature"];
            cv_r = Rcpp::List(cvl.size());
            cvs.assign(ps->units.size(), std::vector<double>());
            for (int b = 0; b < cvl.size(); ++b) {
                if (!kb_todo[b]) { cv_r[b] = R_NilValue; continue; }
                cv_r[b] = cvl[b];
                cvs[kb_unit[b]] = dvec(cvl[b]);
            }
        } else {
            Engine E;
            try {
                engine_setup(E, model, pen);
                read_method(E, method);
                std::vector<char> froz(E.M.p, 0);
                for (int j : hw) froz[j] = 1;
                if (!iwls_loop(E, start, froz, grp, false, o)) return false;
            } catch (RError&) {
                return false;
            }
        }
        mode = o.beta;
        cvs_last = cvs;
        // the kinked coordinates held at zero (zero_kinked())
        std::vector<int> zk;
        for (int j : kc) if (mode[j] == 0) zk.push_back(j);
        // a penalty whose Hessian moves with the coefficients (a t prior, a
        // kink) enters the outer context as it stands at the mode
        if (moving_pen(pen_outer)) {
            if (ps) {
                pen7::PenOut po;
                if (!pen7::build_pen(*ps, v, mode, po, has_kink ? &cvs : nullptr)) return false;
                pen_outer = po.pen_outer;
                mvt_last = po.mvt;
            } else {
                ready();
                Rcpp::RObject r = pen_outer_at(Rcpp::NumericVector(mode.begin(), mode.end()));
                if (Rf_isNull(r)) return false;
                pen_outer = r;
            }
        }
        // the outer context at the mode
        O.reset(new Outer());
        O->model = model;
        O->pen = pen_outer;
        try {
            engine_setup(O->E, model, pen_outer);
            // beside a kink the criterion's matrix is pinned at every kinked
            // coordinate and the mode's at the ones at zero
            O->hold = has_kink ? lap : hold;
            if (has_kink) {
                O->two = true;
                O->hold_mode = hold;
                O->hold_mode.insert(O->hold_mode.end(), zk.begin(), zk.end());
            }
            if (ml) { O->ml = true; O->B = B; }
            O->setup(mode, true);
        } catch (RError&) {
            O->give_back("the R route signals an error");
        }
        if (O->status) return false;
        const int p = O->p;
        // the inner score, free of the held and estimated coordinates
        std::vector<double> g;
        if (!O->E.gr(mode, g)) return false;
        for (int j : gam) g.at(j) = 0.0;
        for (int j : held_free) g.at(j) = 0.0;
        for (int j : zk) g[j] = 0.0;
        if (!o.converged) {
            // inner_mode_error() on the context's inverse
            bool usable = false;
            if (O->ensure_inverse()) {
                std::vector<double> db(p), pr(p);
                dense_matvec(O->Kinv, g.data(), db.data());
                bool fin = true;
                for (double x : db) if (!R_FINITE(x)) fin = false;
                if (fin) {
                    for (int i = 0; i < p; ++i) pr[i] = g[i] * db[i];
                    const double q = 0.5 * r_sum(pr.data(), p);
                    usable = R_FINITE(q) && !(q > 1e-3);
                }
            }
            if (!usable) {
                out.ok = false;
                return true;
            }
        }
        if (O->crit_na) {
            out.ok = false;
            return true;
        }
        value = O->ll - O->rho + O->q / 2.0 * std::log(2.0 * M_PI) - O->logdet_c / 2.0;
        // the score for the resolution, the pinned coordinates zeroed (D4)
        sc = g;
        for (int j : pinned) sc.at(j) = 0.0;
        sc_ok = true;
        out.ok = true;
        out.value = value;
        out.loglik = O->ll;
        out.penalty = O->rho;
        out.par = mode;
        return true;
    }
};

} // namespace

PointSource* make_compiled_point_source(Rcpp::List cfg, Rcpp::List cb, PointSource* fallback) {
    return new CompiledPointSource(cfg, cb, fallback);
}
