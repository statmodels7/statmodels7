#include <Rcpp.h>
using namespace Rcpp;

// Cyclic coordinate descent on a weighted quadratic with a separable penalty.
//
// The loss is (1/2) sum_i w_i (z_i - x_i'beta)^2 and the penalty is applied
// one coordinate at a time. The update of coordinate j minimizes
//
//     (v_j/2) (b - u_j)^2 + rho_j(b),   v_j = sum_i w_i x_ij^2,
//     u_j = beta_j + g_j / v_j,         g_j = sum_i w_i x_ij r_i,
//
// which is the penalty's own proximal operator at the step 1/v_j. That
// operator arrives as a piecewise linear table rather than as a callback: it
// is applied once per coordinate per sweep at a point that moves every time,
// so calling back into R for it would spend the whole gain on the calls.
// The kernel names no family; it evaluates whatever map the table describes.
//
// Two ways to keep g_j, and which one is cheaper depends on the shape of the
// problem rather than on the penalty:
//
//   naive       hold the residual r and read g_j = sum_i w_i x_ij r_i, which
//               is O(n) per visit and O(n) again to update r after a change;
//   covariance  hold g itself, using
//                   g_j = (X'Wz)_j - sum_k (X'WX)_jk beta_k,
//               so that a change of beta_k by d costs one pass over the
//               screened set, O(m), and the column (X'WX)_.k is computed once
//               the first time coordinate k moves off zero. It wins when n is
//               large next to the number of live coordinates, which is the
//               regime a kinked penalty is used in, and pays for it in the
//               memory the cached columns take.

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
                   NumericMatrix icept, IntegerVector screen, int maxit,
                   double tol, bool covariance) {
  int n = X.nrow(), p = X.ncol(), m = screen.size();
  NumericVector beta = clone(beta0);
  std::vector<double> v(p, 0.0);

  for (int a = 0; a < m; a++) {
    int j = screen[a];
    double s = 0.0;
    if (n > 0) {
      const double* xj = &X(0, j);
      for (int i = 0; i < n; i++) s += w[i] * xj[i] * xj[i];
    }
    v[j] = s;
  }

  std::vector<double> r(n), g(m, 0.0);
  std::vector< std::vector<double> > gram(m);   // lazily filled, over screen
  if (covariance) {
    for (int a = 0; a < m; a++) {
      const double* xa = &X(0, screen[a]);
      double s = 0.0;
      for (int i = 0; i < n; i++) s += w[i] * xa[i] * z[i];
      g[a] = s;
    }
    for (int b = 0; b < m; b++) {
      double bk = beta[screen[b]];
      if (bk == 0.0) continue;
      if (gram[b].empty()) {
        gram[b].resize(m);
        const double* xb = &X(0, screen[b]);
        for (int a = 0; a < m; a++) {
          const double* xa = &X(0, screen[a]);
          double s = 0.0;
          for (int i = 0; i < n; i++) s += w[i] * xa[i] * xb[i];
          gram[b][a] = s;
        }
      }
      for (int a = 0; a < m; a++) g[a] -= gram[b][a] * bk;
    }
  } else {
    for (int i = 0; i < n; i++) {
      double e = 0.0;
      for (int j = 0; j < p; j++) e += X(i, j) * beta[j];
      r[i] = z[i] - e;
    }
  }

  std::vector<int> active;
  int sweeps = 0;
  bool full = true;
  for (int it = 1; it <= maxit; it++) {
    sweeps = it;
    double delta = 0.0;
    int q_end = full ? m : (int) active.size();
    for (int q = 0; q < q_end; q++) {
      int a = full ? q : active[q];
      int j = screen[a];
      if (v[j] <= 0.0) continue;
      double gj;
      if (covariance) {
        gj = g[a];
      } else {
        const double* xj = &X(0, j);
        gj = 0.0;
        for (int i = 0; i < n; i++) gj += w[i] * xj[i] * r[i];
      }
      double u = beta[j] + gj / v[j];
      double nb = apply_table(u, a, cut, slope, icept);
      double d = nb - beta[j];
      if (d != 0.0) {
        if (covariance) {
          if (gram[a].empty()) {
            gram[a].resize(m);
            const double* xa = &X(0, j);
            for (int b = 0; b < m; b++) {
              const double* xb = &X(0, screen[b]);
              double s = 0.0;
              for (int i = 0; i < n; i++) s += w[i] * xb[i] * xa[i];
              gram[a][b] = s;
            }
          }
          for (int b = 0; b < m; b++) g[b] -= gram[a][b] * d;
        } else {
          const double* xj = &X(0, j);
          for (int i = 0; i < n; i++) r[i] -= xj[i] * d;
        }
        beta[j] = nb;
        double ad = std::abs(d);
        if (ad > delta) delta = ad;
      }
    }
    if (full) {
      active.clear();
      for (int a = 0; a < m; a++) if (beta[screen[a]] != 0.0) active.push_back(a);
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

  // The gradient over EVERY column, at the point reached: what the screening
  // has to be checked against, since a strong rule is a heuristic and only
  // this comparison makes the answer exact. Where nothing was screened out
  // there is nothing to check, and the pass is skipped -- computing it
  // unconditionally cost more than the screening saved.
  NumericVector grad(p);
  if (m < p) {
    std::vector<double> res(n);
    for (int i = 0; i < n; i++) {
      double e = 0.0;
      for (int j = 0; j < p; j++) e += X(i, j) * beta[j];
      res[i] = z[i] - e;
    }
    for (int j = 0; j < p; j++) {
      const double* xj = &X(0, j);
      double s = 0.0;
      for (int i = 0; i < n; i++) s += w[i] * xj[i] * res[i];
      grad[j] = s;
    }
  }

  return List::create(_["beta"] = beta, _["sweeps"] = sweeps,
                      _["grad"] = grad);
}
