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
//
// DENSE AND SPARSE SHARE THIS ALGORITHM. Every access to the design is a walk
// over one column, so the two differ only in what that walk visits: every row
// for a dense block, the stored nonzeros for a compressed-column one. The
// algorithm is written once against an accessor and instantiated twice, which
// is what keeps the two from drifting apart. A coordinate update on a sparse
// column touching only its nonzeros is not an optimization bolted on: it is
// the algorithm this method has always been, read on the storage it is given.
//
// The two paths agree BIT FOR BIT and not merely to a tolerance, which is a
// property of the arithmetic rather than luck: skipping a structural zero
// omits `s += w[i] * 0.0 * y[i]`, and adding zero to a running sum is exact.
// The tests assert identity, so a reordering would fail them.

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

// --- the two ways to walk a column ------------------------------------------

struct DenseCols {
  const double* X;
  int n, p;
  DenseCols(const NumericMatrix& M)
    : X(M.begin()), n(M.nrow()), p(M.ncol()) {}
  template <class F> inline void each(int j, F f) const {
    const double* xj = X + (std::size_t) j * (std::size_t) n;
    for (int i = 0; i < n; i++) f(i, xj[i]);
  }
};

// A dgCMatrix is compressed by column, so column j's stored entries are
// Ax[Ap[j] .. Ap[j+1]-1] at rows Ai[..] -- exactly the walk this kernel
// wants, and the reason no conversion to a dense buffer is needed anywhere.
struct SparseCols {
  const int* Ai;
  const int* Ap;
  const double* Ax;
  int n, p;
  SparseCols(const IntegerVector& i, const IntegerVector& pp,
             const NumericVector& x, int nrow, int ncol)
    : Ai(i.begin()), Ap(pp.begin()), Ax(x.begin()), n(nrow), p(ncol) {}
  template <class F> inline void each(int j, F f) const {
    for (int k = Ap[j]; k < Ap[j + 1]; k++) f(Ai[k], Ax[k]);
  }
};

// --- the algorithm, written once --------------------------------------------

template <class ACC>
List coord_run(const ACC& A, NumericVector z, NumericVector w,
               NumericVector beta0, NumericMatrix cut, NumericMatrix slope,
               NumericMatrix icept, IntegerVector screen, int maxit,
               double tol, bool covariance) {
  const int n = A.n, p = A.p, m = screen.size();
  NumericVector beta = clone(beta0);
  std::vector<double> v(p, 0.0);

  for (int a = 0; a < m; a++) {
    int j = screen[a];
    double s = 0.0;
    A.each(j, [&](int i, double x) { s += w[i] * x * x; });
    v[j] = s;
  }

  std::vector<double> r(n), g(m, 0.0);
  std::vector< std::vector<double> > gram(m);   // lazily filled, over screen
  // a scatter buffer for the gram entries: one column is written into it and
  // the other walked against it, so a sparse pair costs its own nonzeros
  // rather than n. Cleared through the same walk that filled it.
  std::vector<double> work(covariance ? n : 0, 0.0);

  if (covariance) {
    for (int a = 0; a < m; a++) {
      double s = 0.0;
      A.each(screen[a], [&](int i, double x) { s += w[i] * x * z[i]; });
      g[a] = s;
    }
    for (int b = 0; b < m; b++) {
      double bk = beta[screen[b]];
      if (bk == 0.0) continue;
      if (gram[b].empty()) {
        gram[b].resize(m);
        A.each(screen[b], [&](int i, double x) { work[i] = x; });
        for (int a = 0; a < m; a++) {
          double s = 0.0;
          A.each(screen[a], [&](int i, double x) { s += w[i] * x * work[i]; });
          gram[b][a] = s;
        }
        A.each(screen[b], [&](int i, double) { work[i] = 0.0; });
      }
      for (int a = 0; a < m; a++) g[a] -= gram[b][a] * bk;
    }
  } else {
    for (int i = 0; i < n; i++) r[i] = z[i];
    for (int j = 0; j < p; j++) {
      double bj = beta[j];
      if (bj == 0.0) continue;
      A.each(j, [&](int i, double x) { r[i] -= x * bj; });
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
        gj = 0.0;
        A.each(j, [&](int i, double x) { gj += w[i] * x * r[i]; });
      }
      double u = beta[j] + gj / v[j];
      double nb = apply_table(u, a, cut, slope, icept);
      double d = nb - beta[j];
      if (d != 0.0) {
        if (covariance) {
          if (gram[a].empty()) {
            gram[a].resize(m);
            A.each(j, [&](int i, double x) { work[i] = x; });
            for (int b = 0; b < m; b++) {
              double s = 0.0;
              A.each(screen[b], [&](int i, double x) {
                s += w[i] * x * work[i];
              });
              gram[a][b] = s;
            }
            A.each(j, [&](int i, double) { work[i] = 0.0; });
          }
          for (int b = 0; b < m; b++) g[b] -= gram[a][b] * d;
        } else {
          A.each(j, [&](int i, double x) { r[i] -= x * d; });
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
    for (int i = 0; i < n; i++) res[i] = z[i];
    for (int j = 0; j < p; j++) {
      double bj = beta[j];
      if (bj == 0.0) continue;
      A.each(j, [&](int i, double x) { res[i] -= x * bj; });
    }
    for (int j = 0; j < p; j++) {
      double s = 0.0;
      A.each(j, [&](int i, double x) { s += w[i] * x * res[i]; });
      grad[j] = s;
    }
  }

  return List::create(_["beta"] = beta, _["sweeps"] = sweeps,
                      _["grad"] = grad);
}

// [[Rcpp::export]]
List coord_descent(NumericMatrix X, NumericVector z, NumericVector w,
                   NumericVector beta0, NumericMatrix cut, NumericMatrix slope,
                   NumericMatrix icept, IntegerVector screen, int maxit,
                   double tol, bool covariance) {
  return coord_run(DenseCols(X), z, w, beta0, cut, slope, icept, screen,
                   maxit, tol, covariance);
}

// The slots of a dgCMatrix, passed as they are stored. Taking the S4 object
// apart in R rather than here keeps this file free of any dependency on the
// Matrix package's C API.
// [[Rcpp::export]]
List coord_descent_sparse(IntegerVector Ai, IntegerVector Ap,
                          NumericVector Ax, int nrow, int ncol,
                          NumericVector z, NumericVector w,
                          NumericVector beta0, NumericMatrix cut,
                          NumericMatrix slope, NumericMatrix icept,
                          IntegerVector screen, int maxit, double tol,
                          bool covariance) {
  return coord_run(SparseCols(Ai, Ap, Ax, nrow, ncol), z, w, beta0, cut,
                   slope, icept, screen, maxit, tol, covariance);
}
