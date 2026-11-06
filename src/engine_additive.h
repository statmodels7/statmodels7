// additive_stable() of penalties7 (additive_penalty.R), the stable
// log pseudo-determinant of S = sum_k lam_k P_k with its first and second
// derivatives in lam (Wood 2011, appendix B), for the compiled outer
// criterion (piano_motore.md, D28). The recursion is the R function's,
// level by level: split the active components into a dominant group, peel
// the range of its sum, reduce the others onto the kernel, recurse. The
// closures tr1() of the R function, which read tr(S^+ X) through the levels,
// are the stack of levels kept here.
#ifndef SM7_ENGINE_ADDITIVE_H
#define SM7_ENGINE_ADDITIVE_H

namespace {

inline Dense ad_mul(const Dense& A, const Dense& B) { Dense C; dense_matmul(A, B, C); return C; }
inline Dense ad_tmul(const Dense& A, const Dense& B) {  // crossprod(A, B)
    Dense At = transpose(A), C;
    dense_matmul(At, B, C);
    return C;
}
inline Dense ad_sym(const Dense& Z) {
    Dense S(Z.nr, Z.nc);
    for (int j = 0; j < Z.nc; ++j)
        for (int i = 0; i < Z.nr; ++i) S(i, j) = (Z(i, j) + Z(j, i)) / 2.0;
    return S;
}
inline double ad_maxabs(const Dense& P) {
    double m = 0.0;
    for (double v : P.v) m = std::max(m, std::fabs(v));
    return m;
}
inline double ad_dot(const Dense& A, const Dense& B) {  // sum(A * B)
    std::vector<double> e(A.v.size());
    for (std::size_t t = 0; t < e.size(); ++t) e[t] = A.v[t] * B.v[t];
    return r_sum(e.data(), e.size());
}
inline Dense ad_cols(const Dense& V, const std::vector<int>& k) {
    Dense O(V.nr, (int) k.size());
    for (std::size_t j = 0; j < k.size(); ++j)
        for (int i = 0; i < V.nr; ++i) O(i, (int) j) = V(i, k[j]);
    return O;
}

inline Dense ad_chol2inv(const Dense& U) {
    Dense V = U;
    int m = U.nr, info = 0;
    if (m) F77_CALL(dpotri)("U", &m, V.v.data(), &m, &info FCONE);
    for (int j = 0; j < m; ++j)
        for (int i = j + 1; i < m; ++i) V(i, j) = V(j, i);
    return V;
}

struct AdLevel {
    bool base = false;
    Dense Sp;                  // base: S^+
    Dense Up, U0, Fm, Ainv;    // a peeled level
};
struct AdResult {
    double value = NAN;
    std::vector<double> g;
    Dense T2;
    std::vector<AdLevel> lv;   // from this level down
};

// tr(S^+ X) through the levels from `from` down
inline double ad_tr1(const std::vector<AdLevel>& lv, std::size_t from, const Dense& X) {
    const AdLevel& L = lv[from];
    if (L.base) return ad_dot(L.Sp, X);
    Dense Xa = ad_tmul(L.Up, ad_mul(X, L.Up));
    Dense Xc = ad_tmul(L.Up, ad_mul(X, L.U0));
    Dense Xb = ad_tmul(L.U0, ad_mul(X, L.U0));
    Dense A1 = ad_tmul(Xc, L.Fm), A2 = ad_tmul(L.Fm, Xc), A3 = ad_tmul(L.Fm, ad_mul(Xa, L.Fm));
    Dense Y(Xb.nr, Xb.nc);
    for (std::size_t t = 0; t < Y.v.size(); ++t) Y.v[t] = Xb.v[t] - A1.v[t] - A2.v[t] + A3.v[t];
    return ad_dot(L.Ainv, Xa) + ad_tr1(lv, from + 1, ad_sym(Y));
}

inline int ad_grp_rank(const std::vector<Dense>& R, const std::vector<int>& dom, double rtol) {
    const int n = R[dom[0]].nr;
    Dense St(n, n);
    for (int k : dom) {
        const double m = ad_maxabs(R[k]);
        for (std::size_t t = 0; t < St.v.size(); ++t) St.v[t] += R[k].v[t] / m;
    }
    std::vector<double> ev;
    Dense vv;
    if (!r_eigen_sym(St, true, ev, vv)) return -1;
    const double mx = ev.empty() ? 0.0 : ev[0];
    int c = 0;
    for (double v : ev) if (v > rtol * mx) ++c;
    return c;
}

inline bool ad_rec(const std::vector<Dense>& R, const std::vector<double>& lam,
                   std::vector<char> peeled, const Dense* E, int r,
                   double dtol, double rtol, AdResult& out) {
    const int K = (int) R.size(), n = R[0].nr;
    std::vector<int> act;
    for (int k = 0; k < K; ++k) if (!peeled[k]) act.push_back(k);
    Dense cur(n, n);
    for (int k : act)
        for (std::size_t t = 0; t < cur.v.size(); ++t) cur.v[t] += lam[k] * R[k].v[t];
    if (E) for (std::size_t t = 0; t < cur.v.size(); ++t) cur.v[t] += E->v[t];
    cur = ad_sym(cur);
    for (double v : cur.v) if (!R_FINITE(v)) return false;
    std::vector<double> scal;
    double smax = 0.0;
    for (int k : act) { scal.push_back(lam[k] * ad_maxabs(R[k])); smax = std::max(smax, scal.back()); }
    std::vector<int> dom;
    for (std::size_t i = 0; i < act.size(); ++i) if (scal[i] >= dtol * smax) dom.push_back(act[i]);
    int ra = r;
    if (dom.size() != act.size()) {
        const int gr = ad_grp_rank(R, dom, rtol);
        if (gr < 0) return false;
        ra = std::min(gr, r);
    }
    out.g.assign(K, 0.0);
    out.T2 = Dense(K, K);
    if (dom.size() == act.size() || ra == 0 || ra >= r) {
        std::vector<double> ev;
        Dense V;
        if (!r_eigen_sym(cur, false, ev, V)) return false;
        if (r < 1 || r > n) return false;
        if (ev[r - 1] <= 0) return false;
        std::vector<int> k(r);
        for (int i = 0; i < r; ++i) k[i] = i;
        Dense Vk = ad_cols(V, k);
        Dense Vt = transpose(Vk);
        for (int j = 0; j < Vt.nc; ++j)
            for (int i = 0; i < r; ++i) Vt(i, j) = Vt(i, j) / ev[i];
        AdLevel L;
        L.base = true;
        L.Sp = ad_mul(Vk, Vt);
        std::vector<Dense> SpR(K);
        for (int i = 0; i < K; ++i) SpR[i] = ad_mul(L.Sp, R[i]);
        for (int i = 0; i < K; ++i)
            for (int j = i; j < K; ++j) {
                Dense Ti = transpose(SpR[i]);
                const double v = ad_dot(Ti, SpR[j]);
                out.T2(i, j) = v; out.T2(j, i) = v;
            }
        std::vector<double> lg(r);
        for (int i = 0; i < r; ++i) lg[i] = std::log(ev[i]);
        out.value = r_sum(lg.data(), r);
        for (int i = 0; i < K; ++i) out.g[i] = ad_dot(L.Sp, R[i]);
        out.lv.clear();
        out.lv.push_back(L);
        return true;
    }
    Dense Sa(n, n);
    for (int k : dom)
        for (std::size_t t = 0; t < Sa.v.size(); ++t) Sa.v[t] += lam[k] * R[k].v[t];
    std::vector<double> ea;
    Dense EV;
    if (!r_eigen_sym(ad_sym(Sa), false, ea, EV)) return false;
    std::vector<int> ip, i0;
    for (int i = 0; i < n; ++i) (i < ra ? ip : i0).push_back(i);
    AdLevel L;
    L.Up = ad_cols(EV, ip);
    L.U0 = ad_cols(EV, i0);
    std::vector<int> sub;
    for (int k : act) if (std::find(dom.begin(), dom.end(), k) == dom.end()) sub.push_back(k);
    Dense subm(n, n);
    for (int k : sub)
        for (std::size_t t = 0; t < subm.v.size(); ++t) subm.v[t] += lam[k] * R[k].v[t];
    if (E) for (std::size_t t = 0; t < subm.v.size(); ++t) subm.v[t] += E->v[t];
    Dense A = ad_sym(ad_tmul(L.Up, ad_mul(cur, L.Up)));
    Dense Cm = ad_tmul(L.Up, ad_mul(subm, L.U0));
    Dense cA;
    if (!r_chol(A, cA)) return false;
    L.Ainv = ad_chol2inv(cA);
    L.Fm = ad_mul(L.Ainv, Cm);
    std::vector<Dense> ak(K), red(K), JX(K);
    for (int k = 0; k < K; ++k) ak[k] = ad_tmul(L.Up, ad_mul(R[k], L.Up));
    auto isdom = [&](int k) { return std::find(dom.begin(), dom.end(), k) != dom.end(); };
    for (int k = 0; k < K; ++k) {
        Dense FaF = ad_tmul(L.Fm, ad_mul(ak[k], L.Fm));
        if (isdom(k)) { red[k] = ad_sym(FaF); continue; }
        Dense ck = ad_tmul(L.Up, ad_mul(R[k], L.U0));
        Dense bk = ad_tmul(L.U0, ad_mul(R[k], L.U0));
        Dense A1 = ad_tmul(ck, L.Fm), A2 = ad_tmul(L.Fm, ck);
        Dense Y(bk.nr, bk.nc);
        for (std::size_t t = 0; t < Y.v.size(); ++t) Y.v[t] = bk.v[t] - A1.v[t] - A2.v[t] + FaF.v[t];
        red[k] = ad_sym(Y);
    }
    const int n0 = (int) i0.size();
    Dense E2(n0, n0);
    for (int k : dom)
        for (std::size_t t = 0; t < E2.v.size(); ++t) E2.v[t] += lam[k] * red[k].v[t];
    if (E) {
        Dense Ea = ad_tmul(L.Up, ad_mul(*E, L.Up));
        Dense Ec = ad_tmul(L.Up, ad_mul(*E, L.U0));
        Dense Eb = ad_tmul(L.U0, ad_mul(*E, L.U0));
        Dense A1 = ad_tmul(Ec, L.Fm), A2 = ad_tmul(L.Fm, Ec), A3 = ad_tmul(L.Fm, ad_mul(Ea, L.Fm));
        Dense Y(n0, n0);
        for (std::size_t t = 0; t < Y.v.size(); ++t) Y.v[t] = Eb.v[t] - A1.v[t] - A2.v[t] + A3.v[t];
        Dense Ys = ad_sym(Y);
        for (std::size_t t = 0; t < E2.v.size(); ++t) E2.v[t] += Ys.v[t];
    }
    std::vector<char> peel2 = peeled;
    for (int k : dom) peel2[k] = 1;
    AdResult nx;
    if (!ad_rec(red, lam, peel2, &E2, r - ra, dtol, rtol, nx)) return false;
    if (!R_FINITE(nx.value)) return false;
    for (int k = 0; k < K; ++k) {
        Dense aF = ad_mul(ak[k], L.Fm);
        if (isdom(k)) {
            JX[k] = aF;
            for (double& z : JX[k].v) z = -z;
        } else {
            Dense c = ad_tmul(L.Up, ad_mul(R[k], L.U0));
            JX[k] = Dense(c.nr, c.nc);
            for (std::size_t t = 0; t < c.v.size(); ++t) JX[k].v[t] = c.v[t] - aF.v[t];
        }
    }
    std::vector<Dense> AiA(K);
    for (int k = 0; k < K; ++k) AiA[k] = ad_mul(L.Ainv, ak[k]);
    std::vector<AdLevel> lv;
    lv.push_back(L);
    lv.insert(lv.end(), nx.lv.begin(), nx.lv.end());
    for (int i = 0; i < K; ++i)
        for (int j = i; j < K; ++j) {
            Dense mid = ad_tmul(JX[j], ad_mul(L.Ainv, JX[i]));
            Dense mm(mid.nr, mid.nc);
            for (int c = 0; c < mid.nc; ++c)
                for (int a = 0; a < mid.nr; ++a) mm(a, c) = mid(a, c) + mid(c, a);
            Dense Ti = transpose(AiA[i]);
            const double v = ad_dot(Ti, AiA[j]) + ad_tr1(lv, 1, mm) + nx.T2(i, j);
            out.T2(i, j) = v; out.T2(j, i) = v;
        }
    std::vector<double> lg(cA.nr);
    for (int i = 0; i < cA.nr; ++i) lg[i] = std::log(cA(i, i));
    out.value = 2.0 * r_sum(lg.data(), cA.nr) + nx.value;
    for (int k = 0; k < K; ++k) out.g[k] = ad_dot(L.Ainv, ak[k]) + nx.g[k];
    out.lv = lv;
    return true;
}

// value, dlog and d2log; false where the R function returns NaN
inline bool additive_stable_c(const std::vector<Dense>& mats, const std::vector<double>& lam,
                              int rank, double rtol, double& value,
                              std::vector<double>& dlog, Dense& d2log) {
    const double dtol = std::pow(DBL_EPSILON, 0.4);
    AdResult res;
    std::vector<char> peeled(mats.size(), 0);
    if (!ad_rec(mats, lam, peeled, nullptr, rank, dtol, rtol, res)) return false;
    if (!R_FINITE(res.value)) return false;
    value = res.value;
    dlog = res.g;
    d2log = res.T2;
    return true;
}

}  // namespace

#endif
