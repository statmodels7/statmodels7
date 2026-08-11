#include <Rcpp.h>
using namespace Rcpp;

// Cyclic coordinate descent on a weighted quadratic with a separable penalty.
//
// The loss is (1/2) sum_i w_i (z_i - x_i'beta)^2 and the penalty is applied
// one coordinate at a time. The update of coordinate j minimizes
//
//     (v_j/2) (b - u_j)^2 + rho_j(b),   v_j = sum_i w_i x_ij^2,
//     u_j = beta_j + (sum_i w_i x_ij r_i) / v_j,
//
// which is the penalty's own proximal operator at the step 1/v_j. That
// operator arrives as a piecewise linear table rather than as a callback: it
// is applied once per coordinate per sweep at a point that moves every time,
// so calling back into R for it would spend the whole gain on the calls.
// The kernel names no family; it evaluates whatever map the table describes.
//
// The table is odd, so it acts on |u| and the sign is restored: on the k-th
// piece, |u| <= cut(j, k), the value is slope(j, k) * |u| + icept(j, k).

static inline double apply_table(double u, int j, const NumericMatrix& cut,
                                 const NumericMatrix& slope,
                                 const NumericMatrix& icept) {
  double a = std::abs(u);
  int k = cut.ncol() - 1;
  for (int q = 0; q < cut.ncol(); q++) {
    if (a <= cut(j, q)) { k = q; break; }
  }
  double val = slope(j, k) * a + icept(j, k);
  return u < 0 ? -val : val;
}

// [[Rcpp::export]]
List coord_descent(NumericMatrix X, NumericVector z, NumericVector w,
                   NumericVector beta0, NumericMatrix cut, NumericMatrix slope,
                   NumericMatrix icept, int maxit, double tol) {
  int n = X.nrow(), p = X.ncol();
  NumericVector beta = clone(beta0);
  std::vector<double> v(p), r(n), wr(n);
  std::vector<int> active;

  for (int j = 0; j < p; j++) {
    double s = 0.0;
    if (n > 0) {
      const double* xj = &X(0, j);
      for (int i = 0; i < n; i++) s += w[i] * xj[i] * xj[i];
    }
    v[j] = s;
  }
  for (int i = 0; i < n; i++) {
    double e = 0.0;
    for (int j = 0; j < p; j++) e += X(i, j) * beta[j];
    r[i] = z[i] - e;
  }

  int sweeps = 0;
  bool full = true;
  for (int it = 1; it <= maxit; it++) {
    sweeps = it;
    double delta = 0.0;
    // a full sweep visits every coordinate; the passes between visit only the
    // ones that are away from the kink, which is where the answer moves
    int m = full ? p : (int) active.size();
    for (int q = 0; q < m; q++) {
      int j = full ? q : active[q];
      if (v[j] <= 0.0) continue;
      const double* xj = &X(0, j);
      double g = 0.0;
      for (int i = 0; i < n; i++) g += w[i] * xj[i] * r[i];
      double u = beta[j] + g / v[j];
      double nb = apply_table(u, j, cut, slope, icept);
      double d = nb - beta[j];
      if (d != 0.0) {
        for (int i = 0; i < n; i++) r[i] -= xj[i] * d;
        beta[j] = nb;
        double ad = std::abs(d);
        if (ad > delta) delta = ad;
      }
    }
    if (full) {
      active.clear();
      for (int j = 0; j < p; j++) if (beta[j] != 0.0) active.push_back(j);
    }
    if (delta < tol) {
      // a converged inner pass has to be confirmed by a full sweep: a
      // coordinate at zero can still want to move, and only a full sweep asks
      if (full) break;
      full = true;
    } else {
      full = active.empty();
    }
  }

  double rss = 0.0;
  for (int i = 0; i < n; i++) rss += w[i] * r[i] * r[i];
  return List::create(_["beta"] = beta, _["sweeps"] = sweeps,
                      _["rss"] = 0.5 * rss);
}
