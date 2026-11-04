// The derivatives of the scaling matrix S = I^-d of a joint score-driven
// filter (d = 1 or 1/2, I the expected information of the filtered
// parameters on the link scale), and of the scaled driving vector u = S g,
// in the predictors. One function per order: jsc_s<m>_cpp() returns S and
// M (M = I for d = 1, the symmetric root for d = 1/2) for every multiset of
// m predictor indices, jsc_u<m>_cpp() the driving vector's derivative for
// every multiset, and each lists its Leibniz terms explicitly.
//
// Storage. A matrix component of order k is an array n x K x K x C_k (C_k
// the number of multisets of size k over the P predictors, in the order of
// joint_scale_multisets()), an entry at observation i, row a, column c and
// multiset b being x[i + n (a + K c) + n K K b]. A vector component is
// n x K x C_k. The table `tab` of order m has one row per multiset B of size
// m and one column per subset A of B's positions (a bit mask, position j
// being bit j): the index, among the multisets of size |A|, of the
// sub-multiset B[A]. Mask 0 is the empty multiset and the full mask is B.
//
// For d = 1, S_B = -S sum_{A nonempty} I_A S_{B\A} (Leibniz on I S = 1).
// For d = 1/2, with X the root, X X_B + X_B X = I_B - sum_{A proper} X_A
// X_{B\A} is solved in the eigenbasis of I, dividing by sqrt(l_a) + sqrt(l_b),
// and S_B = -S sum_{A nonempty} X_A S_{B\A}. u_B = sum_{A} S_{B\A} g_A with
// g_A the family's derivative of order 1 + |A|.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>
#include <vector>

#ifndef FCONE
# define FCONE
#endif

using namespace Rcpp;

namespace {

// out += X Y, X and Y stored with stride n at observation i, out K x K
inline void mm_add(double* out, const double* X, const double* Y, int i,
                   int n, int K) {
  for (int c = 0; c < K; ++c) {
    for (int a = 0; a < K; ++a) {
      double s = 0.0;
      for (int k = 0; k < K; ++k) {
        s += X[i + (size_t) n * (a + K * k)] * Y[i + (size_t) n * (k + K * c)];
      }
      out[a + K * c] += s;
    }
  }
}

// out -= X Y
inline void mm_sub(double* out, const double* X, const double* Y, int i,
                   int n, int K) {
  for (int c = 0; c < K; ++c) {
    for (int a = 0; a < K; ++a) {
      double s = 0.0;
      for (int k = 0; k < K; ++k) {
        s += X[i + (size_t) n * (a + K * k)] * Y[i + (size_t) n * (k + K * c)];
      }
      out[a + K * c] -= s;
    }
  }
}

// u += S g, S stored with stride n (K x K), g with stride n (K)
inline void mv_add(double* u, const double* S, const double* g, int i, int n,
                   int K) {
  for (int l = 0; l < K; ++l) {
    double s = 0.0;
    for (int c = 0; c < K; ++c) {
      s += S[i + (size_t) n * (l + K * c)] * g[i + (size_t) n * c];
    }
    u[l] += s;
  }
}

// out (stride n, offset at observation i) = -S0 acc
inline void neg_left(double* out, const double* S0, const double* acc, int i,
                     int n, int K) {
  for (int c = 0; c < K; ++c) {
    for (int a = 0; a < K; ++a) {
      double s = 0.0;
      for (int k = 0; k < K; ++k) {
        s += S0[i + (size_t) n * (a + K * k)] * acc[k + K * c];
      }
      out[i + (size_t) n * (a + K * c)] = -s;
    }
  }
}

// out (stride n) = Q ((Q' R Q) / (r_a + r_b)) Q', the solution X of
// X0 X + X X0 = R with X0 = Q diag(r) Q'
inline void sylvester(double* out, const double* Q, const double* r,
                      const double* R, int i, int n, int K,
                      std::vector<double>& T1, std::vector<double>& T2) {
  // T1 = Q' R
  for (int c = 0; c < K; ++c) {
    for (int a = 0; a < K; ++a) {
      double s = 0.0;
      for (int k = 0; k < K; ++k) s += Q[i + (size_t) n * (k + K * a)] * R[k + K * c];
      T1[a + K * c] = s;
    }
  }
  // T2 = (T1 Q) / den
  for (int c = 0; c < K; ++c) {
    for (int a = 0; a < K; ++a) {
      double s = 0.0;
      for (int k = 0; k < K; ++k) s += T1[a + K * k] * Q[i + (size_t) n * (k + K * c)];
      T2[a + K * c] = s / (r[i + (size_t) n * a] + r[i + (size_t) n * c]);
    }
  }
  // T1 = Q T2
  for (int c = 0; c < K; ++c) {
    for (int a = 0; a < K; ++a) {
      double s = 0.0;
      for (int k = 0; k < K; ++k) s += Q[i + (size_t) n * (a + K * k)] * T2[k + K * c];
      T1[a + K * c] = s;
    }
  }
  // out = T1 Q'
  for (int c = 0; c < K; ++c) {
    for (int a = 0; a < K; ++a) {
      double s = 0.0;
      for (int k = 0; k < K; ++k) s += T1[a + K * k] * Q[i + (size_t) n * (c + K * k)];
      out[i + (size_t) n * (a + K * c)] = s;
    }
  }
}

inline int popcount(int m) {
  int c = 0;
  while (m) { c += m & 1; m >>= 1; }
  return c;
}

} // namespace


// The Scaling Matrix and Its Eigenbasis at Order Zero
//
// I0 The expected information of the filtered parameters, n x K x K.
// n,K The observations and the filtered parameters.
// d The power, 1/2 or 1.
// Returns A list with `S` (I^-d), `Q` (eigenvectors, n x K x K), `r`
//   (eigenvalues, or their square roots for d = 1/2, n x K) and `ok`
//   (`FALSE` where an eigenvalue is not positive).
// [[Rcpp::export]]
List jsc_s0_cpp(NumericVector I0, int n, int K, double d) {
  NumericVector S(I0.size()), Q(I0.size()), r((size_t) n * K);
  std::vector<double> A((size_t) K * K), w(K);
  int lwork = -1, info = 0;
  double wq = 0.0;
  F77_CALL(dsyev)("V", "L", &K, A.data(), &K, w.data(), &wq, &lwork, &info
                  FCONE FCONE);
  lwork = (int) wq;
  if (lwork < 3 * K) lwork = 3 * K;
  std::vector<double> work(lwork);
  bool ok = true;
  for (int i = 0; i < n; ++i) {
    for (int c = 0; c < K; ++c) {
      for (int a = 0; a < K; ++a) A[a + K * c] = I0[i + (size_t) n * (a + K * c)];
    }
    F77_CALL(dsyev)("V", "L", &K, A.data(), &K, w.data(), work.data(), &lwork,
                    &info FCONE FCONE);
    if (info != 0) stop("dsyev failed at observation %d (info %d).", i + 1, info);
    for (int k = 0; k < K; ++k) {
      if (!(w[k] > 0)) ok = false;
      r[i + (size_t) n * k] = (d == 0.5) ? std::sqrt(w[k]) : w[k];
    }
    for (int c = 0; c < K; ++c) {
      for (int a = 0; a < K; ++a) Q[i + (size_t) n * (a + K * c)] = A[a + K * c];
    }
    for (int c = 0; c < K; ++c) {
      for (int a = 0; a < K; ++a) {
        double s = 0.0;
        for (int k = 0; k < K; ++k) {
          s += A[a + K * k] * A[c + K * k] / r[i + (size_t) n * k];
        }
        S[i + (size_t) n * (a + K * c)] = s;
      }
    }
  }
  return List::create(_["S"] = S, _["Q"] = Q, _["r"] = r, _["ok"] = ok);
}


// The Scaled Driving Vector at Order Zero
//
// S0 The scaling matrix, n x K x K.
// g0 The link-scale gradient of the filtered parameters, n x K.
// n,K The observations and the filtered parameters.
// Returns u = S g, n x K.
// [[Rcpp::export]]
NumericVector jsc_u0_cpp(NumericVector S0, NumericVector g0, int n, int K) {
  NumericVector out((size_t) n * K);
  std::vector<double> u(K);
  for (int i = 0; i < n; ++i) {
    for (int l = 0; l < K; ++l) u[l] = 0.0;
    mv_add(u.data(), S0.begin(), g0.begin(), i, n, K);
    for (int l = 0; l < K; ++l) out[i + (size_t) n * l] = u[l];
  }
  return out;
}


// The Scaling Matrix's Derivatives of Order 1
//
// Im The information's derivatives of order 1, n x K x K x C_1.
// M The components of M of orders 1 to 0 (a list; `NULL` at
//   order 1).
// S The components of S of orders 0 to 0.
// tab The table of sub-multisets of order 1.
// Q,r The eigenbasis of the information, from [jsc_s0_cpp()].
// n,K The observations and the filtered parameters.
// d The power, 1/2 or 1.
// Returns A list with `M` and `S` of order 1.
// [[Rcpp::export]]
List jsc_s1_cpp(NumericVector Im, List M, List S, IntegerMatrix tab,
                NumericVector Q, NumericVector r, int n, int K, double d) {
  const int m = 1, F = 1;
  const int C = tab.nrow();
  const size_t blk = (size_t) n * K * K;
  std::vector<const double*> Mlow(m), Slow(m);
  for (int k = 1; k < m; ++k) Mlow[k] = REAL(M[k - 1]);
  for (int k = 0; k < m; ++k) Slow[k] = REAL(S[k]);
  NumericVector Mout = (d == 1.0) ? Im : NumericVector(Im.size());
  NumericVector Sout(Im.size());
  std::vector<double> R((size_t) K * K), acc((size_t) K * K), T1((size_t) K * K),
    T2((size_t) K * K);
  const double* mp[F + 1];
  const double* sp[F + 1];
  for (int b = 0; b < C; ++b) {
    for (int A = 1; A < F; ++A) {
      mp[A] = Mlow[popcount(A)] + blk * tab(b, A);
      sp[A] = Slow[popcount(A)] + blk * tab(b, A);
    }
    sp[0] = Slow[0];
    mp[F] = Mout.begin() + blk * b;
    if (d != 1.0) {
      const double* Ib = Im.begin() + blk * b;
      for (int i = 0; i < n; ++i) {
        for (int c = 0; c < K; ++c) {
          for (int a = 0; a < K; ++a) R[a + K * c] = Ib[i + (size_t) n * (a + K * c)];
        }
        // no proper split at order 1
        sylvester(Mout.begin() + blk * b, Q.begin(), r.begin(), R.data(), i, n, K,
                  T1, T2);
      }
    }
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < K * K; ++k) acc[k] = 0.0;
      mm_add(acc.data(), mp[1], sp[0], i, n, K);
      neg_left(Sout.begin() + blk * b, Slow[0], acc.data(), i, n, K);
    }
  }
  return List::create(_["M"] = Mout, _["S"] = Sout);
}


// The Scaled Driving Vector's Derivatives of Order 1
//
// S The components of S of orders 0 to 1.
// g The family's derivatives, by order 0 to 1 of the multiset,
//   each n x K x C_k.
// tab The table of sub-multisets of order 1.
// n,K The observations and the filtered parameters.
// Returns The derivative for every multiset of order 1, n x K x C_1.
// [[Rcpp::export]]
NumericVector jsc_u1_cpp(List S, List g, IntegerMatrix tab, int n, int K) {
  const int m = 1, F = 1;
  const int C = tab.nrow();
  const size_t blk = (size_t) n * K * K, vblk = (size_t) n * K;
  std::vector<const double*> Sk(m + 1), gk(m + 1);
  for (int k = 0; k <= m; ++k) {
    Sk[k] = REAL(S[k]);
    gk[k] = REAL(g[k]);
  }
  NumericVector out(vblk * C);
  std::vector<double> u(K);
  const double* sp[F + 1];
  const double* gp[F + 1];
  for (int b = 0; b < C; ++b) {
    for (int A = 0; A <= F; ++A) {
      sp[A] = Sk[popcount(A)] + blk * tab(b, A);
      gp[A] = gk[popcount(A)] + vblk * tab(b, A);
    }
    for (int i = 0; i < n; ++i) {
      for (int l = 0; l < K; ++l) u[l] = 0.0;
      mv_add(u.data(), sp[1], gp[0], i, n, K);
      mv_add(u.data(), sp[0], gp[1], i, n, K);
      for (int l = 0; l < K; ++l) out[i + (size_t) n * l + vblk * b] = u[l];
    }
  }
  return out;
}


// The Scaling Matrix's Derivatives of Order 2
//
// Im The information's derivatives of order 2, n x K x K x C_2.
// M The components of M of orders 1 to 1 (a list; `NULL` at
//   order 1).
// S The components of S of orders 0 to 1.
// tab The table of sub-multisets of order 2.
// Q,r The eigenbasis of the information, from [jsc_s0_cpp()].
// n,K The observations and the filtered parameters.
// d The power, 1/2 or 1.
// Returns A list with `M` and `S` of order 2.
// [[Rcpp::export]]
List jsc_s2_cpp(NumericVector Im, List M, List S, IntegerMatrix tab,
                NumericVector Q, NumericVector r, int n, int K, double d) {
  const int m = 2, F = 3;
  const int C = tab.nrow();
  const size_t blk = (size_t) n * K * K;
  std::vector<const double*> Mlow(m), Slow(m);
  for (int k = 1; k < m; ++k) Mlow[k] = REAL(M[k - 1]);
  for (int k = 0; k < m; ++k) Slow[k] = REAL(S[k]);
  NumericVector Mout = (d == 1.0) ? Im : NumericVector(Im.size());
  NumericVector Sout(Im.size());
  std::vector<double> R((size_t) K * K), acc((size_t) K * K), T1((size_t) K * K),
    T2((size_t) K * K);
  const double* mp[F + 1];
  const double* sp[F + 1];
  for (int b = 0; b < C; ++b) {
    for (int A = 1; A < F; ++A) {
      mp[A] = Mlow[popcount(A)] + blk * tab(b, A);
      sp[A] = Slow[popcount(A)] + blk * tab(b, A);
    }
    sp[0] = Slow[0];
    mp[F] = Mout.begin() + blk * b;
    if (d != 1.0) {
      const double* Ib = Im.begin() + blk * b;
      for (int i = 0; i < n; ++i) {
        for (int c = 0; c < K; ++c) {
          for (int a = 0; a < K; ++a) R[a + K * c] = Ib[i + (size_t) n * (a + K * c)];
        }
        mm_sub(R.data(), mp[1], mp[2], i, n, K);
        mm_sub(R.data(), mp[2], mp[1], i, n, K);
        sylvester(Mout.begin() + blk * b, Q.begin(), r.begin(), R.data(), i, n, K,
                  T1, T2);
      }
    }
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < K * K; ++k) acc[k] = 0.0;
      mm_add(acc.data(), mp[1], sp[2], i, n, K);
      mm_add(acc.data(), mp[2], sp[1], i, n, K);
      mm_add(acc.data(), mp[3], sp[0], i, n, K);
      neg_left(Sout.begin() + blk * b, Slow[0], acc.data(), i, n, K);
    }
  }
  return List::create(_["M"] = Mout, _["S"] = Sout);
}


// The Scaled Driving Vector's Derivatives of Order 2
//
// S The components of S of orders 0 to 2.
// g The family's derivatives, by order 0 to 2 of the multiset,
//   each n x K x C_k.
// tab The table of sub-multisets of order 2.
// n,K The observations and the filtered parameters.
// Returns The derivative for every multiset of order 2, n x K x C_2.
// [[Rcpp::export]]
NumericVector jsc_u2_cpp(List S, List g, IntegerMatrix tab, int n, int K) {
  const int m = 2, F = 3;
  const int C = tab.nrow();
  const size_t blk = (size_t) n * K * K, vblk = (size_t) n * K;
  std::vector<const double*> Sk(m + 1), gk(m + 1);
  for (int k = 0; k <= m; ++k) {
    Sk[k] = REAL(S[k]);
    gk[k] = REAL(g[k]);
  }
  NumericVector out(vblk * C);
  std::vector<double> u(K);
  const double* sp[F + 1];
  const double* gp[F + 1];
  for (int b = 0; b < C; ++b) {
    for (int A = 0; A <= F; ++A) {
      sp[A] = Sk[popcount(A)] + blk * tab(b, A);
      gp[A] = gk[popcount(A)] + vblk * tab(b, A);
    }
    for (int i = 0; i < n; ++i) {
      for (int l = 0; l < K; ++l) u[l] = 0.0;
      mv_add(u.data(), sp[3], gp[0], i, n, K);
      mv_add(u.data(), sp[2], gp[1], i, n, K);
      mv_add(u.data(), sp[1], gp[2], i, n, K);
      mv_add(u.data(), sp[0], gp[3], i, n, K);
      for (int l = 0; l < K; ++l) out[i + (size_t) n * l + vblk * b] = u[l];
    }
  }
  return out;
}


// The Scaling Matrix's Derivatives of Order 3
//
// Im The information's derivatives of order 3, n x K x K x C_3.
// M The components of M of orders 1 to 2 (a list; `NULL` at
//   order 1).
// S The components of S of orders 0 to 2.
// tab The table of sub-multisets of order 3.
// Q,r The eigenbasis of the information, from [jsc_s0_cpp()].
// n,K The observations and the filtered parameters.
// d The power, 1/2 or 1.
// Returns A list with `M` and `S` of order 3.
// [[Rcpp::export]]
List jsc_s3_cpp(NumericVector Im, List M, List S, IntegerMatrix tab,
                NumericVector Q, NumericVector r, int n, int K, double d) {
  const int m = 3, F = 7;
  const int C = tab.nrow();
  const size_t blk = (size_t) n * K * K;
  std::vector<const double*> Mlow(m), Slow(m);
  for (int k = 1; k < m; ++k) Mlow[k] = REAL(M[k - 1]);
  for (int k = 0; k < m; ++k) Slow[k] = REAL(S[k]);
  NumericVector Mout = (d == 1.0) ? Im : NumericVector(Im.size());
  NumericVector Sout(Im.size());
  std::vector<double> R((size_t) K * K), acc((size_t) K * K), T1((size_t) K * K),
    T2((size_t) K * K);
  const double* mp[F + 1];
  const double* sp[F + 1];
  for (int b = 0; b < C; ++b) {
    for (int A = 1; A < F; ++A) {
      mp[A] = Mlow[popcount(A)] + blk * tab(b, A);
      sp[A] = Slow[popcount(A)] + blk * tab(b, A);
    }
    sp[0] = Slow[0];
    mp[F] = Mout.begin() + blk * b;
    if (d != 1.0) {
      const double* Ib = Im.begin() + blk * b;
      for (int i = 0; i < n; ++i) {
        for (int c = 0; c < K; ++c) {
          for (int a = 0; a < K; ++a) R[a + K * c] = Ib[i + (size_t) n * (a + K * c)];
        }
        mm_sub(R.data(), mp[1], mp[6], i, n, K);
        mm_sub(R.data(), mp[2], mp[5], i, n, K);
        mm_sub(R.data(), mp[3], mp[4], i, n, K);
        mm_sub(R.data(), mp[4], mp[3], i, n, K);
        mm_sub(R.data(), mp[5], mp[2], i, n, K);
        mm_sub(R.data(), mp[6], mp[1], i, n, K);
        sylvester(Mout.begin() + blk * b, Q.begin(), r.begin(), R.data(), i, n, K,
                  T1, T2);
      }
    }
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < K * K; ++k) acc[k] = 0.0;
      mm_add(acc.data(), mp[1], sp[6], i, n, K);
      mm_add(acc.data(), mp[2], sp[5], i, n, K);
      mm_add(acc.data(), mp[3], sp[4], i, n, K);
      mm_add(acc.data(), mp[4], sp[3], i, n, K);
      mm_add(acc.data(), mp[5], sp[2], i, n, K);
      mm_add(acc.data(), mp[6], sp[1], i, n, K);
      mm_add(acc.data(), mp[7], sp[0], i, n, K);
      neg_left(Sout.begin() + blk * b, Slow[0], acc.data(), i, n, K);
    }
  }
  return List::create(_["M"] = Mout, _["S"] = Sout);
}


// The Scaled Driving Vector's Derivatives of Order 3
//
// S The components of S of orders 0 to 3.
// g The family's derivatives, by order 0 to 3 of the multiset,
//   each n x K x C_k.
// tab The table of sub-multisets of order 3.
// n,K The observations and the filtered parameters.
// Returns The derivative for every multiset of order 3, n x K x C_3.
// [[Rcpp::export]]
NumericVector jsc_u3_cpp(List S, List g, IntegerMatrix tab, int n, int K) {
  const int m = 3, F = 7;
  const int C = tab.nrow();
  const size_t blk = (size_t) n * K * K, vblk = (size_t) n * K;
  std::vector<const double*> Sk(m + 1), gk(m + 1);
  for (int k = 0; k <= m; ++k) {
    Sk[k] = REAL(S[k]);
    gk[k] = REAL(g[k]);
  }
  NumericVector out(vblk * C);
  std::vector<double> u(K);
  const double* sp[F + 1];
  const double* gp[F + 1];
  for (int b = 0; b < C; ++b) {
    for (int A = 0; A <= F; ++A) {
      sp[A] = Sk[popcount(A)] + blk * tab(b, A);
      gp[A] = gk[popcount(A)] + vblk * tab(b, A);
    }
    for (int i = 0; i < n; ++i) {
      for (int l = 0; l < K; ++l) u[l] = 0.0;
      mv_add(u.data(), sp[7], gp[0], i, n, K);
      mv_add(u.data(), sp[6], gp[1], i, n, K);
      mv_add(u.data(), sp[5], gp[2], i, n, K);
      mv_add(u.data(), sp[4], gp[3], i, n, K);
      mv_add(u.data(), sp[3], gp[4], i, n, K);
      mv_add(u.data(), sp[2], gp[5], i, n, K);
      mv_add(u.data(), sp[1], gp[6], i, n, K);
      mv_add(u.data(), sp[0], gp[7], i, n, K);
      for (int l = 0; l < K; ++l) out[i + (size_t) n * l + vblk * b] = u[l];
    }
  }
  return out;
}


// The Scaling Matrix's Derivatives of Order 4
//
// Im The information's derivatives of order 4, n x K x K x C_4.
// M The components of M of orders 1 to 3 (a list; `NULL` at
//   order 1).
// S The components of S of orders 0 to 3.
// tab The table of sub-multisets of order 4.
// Q,r The eigenbasis of the information, from [jsc_s0_cpp()].
// n,K The observations and the filtered parameters.
// d The power, 1/2 or 1.
// Returns A list with `M` and `S` of order 4.
// [[Rcpp::export]]
List jsc_s4_cpp(NumericVector Im, List M, List S, IntegerMatrix tab,
                NumericVector Q, NumericVector r, int n, int K, double d) {
  const int m = 4, F = 15;
  const int C = tab.nrow();
  const size_t blk = (size_t) n * K * K;
  std::vector<const double*> Mlow(m), Slow(m);
  for (int k = 1; k < m; ++k) Mlow[k] = REAL(M[k - 1]);
  for (int k = 0; k < m; ++k) Slow[k] = REAL(S[k]);
  NumericVector Mout = (d == 1.0) ? Im : NumericVector(Im.size());
  NumericVector Sout(Im.size());
  std::vector<double> R((size_t) K * K), acc((size_t) K * K), T1((size_t) K * K),
    T2((size_t) K * K);
  const double* mp[F + 1];
  const double* sp[F + 1];
  for (int b = 0; b < C; ++b) {
    for (int A = 1; A < F; ++A) {
      mp[A] = Mlow[popcount(A)] + blk * tab(b, A);
      sp[A] = Slow[popcount(A)] + blk * tab(b, A);
    }
    sp[0] = Slow[0];
    mp[F] = Mout.begin() + blk * b;
    if (d != 1.0) {
      const double* Ib = Im.begin() + blk * b;
      for (int i = 0; i < n; ++i) {
        for (int c = 0; c < K; ++c) {
          for (int a = 0; a < K; ++a) R[a + K * c] = Ib[i + (size_t) n * (a + K * c)];
        }
        mm_sub(R.data(), mp[1], mp[14], i, n, K);
        mm_sub(R.data(), mp[2], mp[13], i, n, K);
        mm_sub(R.data(), mp[3], mp[12], i, n, K);
        mm_sub(R.data(), mp[4], mp[11], i, n, K);
        mm_sub(R.data(), mp[5], mp[10], i, n, K);
        mm_sub(R.data(), mp[6], mp[9], i, n, K);
        mm_sub(R.data(), mp[7], mp[8], i, n, K);
        mm_sub(R.data(), mp[8], mp[7], i, n, K);
        mm_sub(R.data(), mp[9], mp[6], i, n, K);
        mm_sub(R.data(), mp[10], mp[5], i, n, K);
        mm_sub(R.data(), mp[11], mp[4], i, n, K);
        mm_sub(R.data(), mp[12], mp[3], i, n, K);
        mm_sub(R.data(), mp[13], mp[2], i, n, K);
        mm_sub(R.data(), mp[14], mp[1], i, n, K);
        sylvester(Mout.begin() + blk * b, Q.begin(), r.begin(), R.data(), i, n, K,
                  T1, T2);
      }
    }
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < K * K; ++k) acc[k] = 0.0;
      mm_add(acc.data(), mp[1], sp[14], i, n, K);
      mm_add(acc.data(), mp[2], sp[13], i, n, K);
      mm_add(acc.data(), mp[3], sp[12], i, n, K);
      mm_add(acc.data(), mp[4], sp[11], i, n, K);
      mm_add(acc.data(), mp[5], sp[10], i, n, K);
      mm_add(acc.data(), mp[6], sp[9], i, n, K);
      mm_add(acc.data(), mp[7], sp[8], i, n, K);
      mm_add(acc.data(), mp[8], sp[7], i, n, K);
      mm_add(acc.data(), mp[9], sp[6], i, n, K);
      mm_add(acc.data(), mp[10], sp[5], i, n, K);
      mm_add(acc.data(), mp[11], sp[4], i, n, K);
      mm_add(acc.data(), mp[12], sp[3], i, n, K);
      mm_add(acc.data(), mp[13], sp[2], i, n, K);
      mm_add(acc.data(), mp[14], sp[1], i, n, K);
      mm_add(acc.data(), mp[15], sp[0], i, n, K);
      neg_left(Sout.begin() + blk * b, Slow[0], acc.data(), i, n, K);
    }
  }
  return List::create(_["M"] = Mout, _["S"] = Sout);
}


// The Scaled Driving Vector's Derivatives of Order 4
//
// S The components of S of orders 0 to 4.
// g The family's derivatives, by order 0 to 4 of the multiset,
//   each n x K x C_k.
// tab The table of sub-multisets of order 4.
// n,K The observations and the filtered parameters.
// Returns The derivative for every multiset of order 4, n x K x C_4.
// [[Rcpp::export]]
NumericVector jsc_u4_cpp(List S, List g, IntegerMatrix tab, int n, int K) {
  const int m = 4, F = 15;
  const int C = tab.nrow();
  const size_t blk = (size_t) n * K * K, vblk = (size_t) n * K;
  std::vector<const double*> Sk(m + 1), gk(m + 1);
  for (int k = 0; k <= m; ++k) {
    Sk[k] = REAL(S[k]);
    gk[k] = REAL(g[k]);
  }
  NumericVector out(vblk * C);
  std::vector<double> u(K);
  const double* sp[F + 1];
  const double* gp[F + 1];
  for (int b = 0; b < C; ++b) {
    for (int A = 0; A <= F; ++A) {
      sp[A] = Sk[popcount(A)] + blk * tab(b, A);
      gp[A] = gk[popcount(A)] + vblk * tab(b, A);
    }
    for (int i = 0; i < n; ++i) {
      for (int l = 0; l < K; ++l) u[l] = 0.0;
      mv_add(u.data(), sp[15], gp[0], i, n, K);
      mv_add(u.data(), sp[14], gp[1], i, n, K);
      mv_add(u.data(), sp[13], gp[2], i, n, K);
      mv_add(u.data(), sp[12], gp[3], i, n, K);
      mv_add(u.data(), sp[11], gp[4], i, n, K);
      mv_add(u.data(), sp[10], gp[5], i, n, K);
      mv_add(u.data(), sp[9], gp[6], i, n, K);
      mv_add(u.data(), sp[8], gp[7], i, n, K);
      mv_add(u.data(), sp[7], gp[8], i, n, K);
      mv_add(u.data(), sp[6], gp[9], i, n, K);
      mv_add(u.data(), sp[5], gp[10], i, n, K);
      mv_add(u.data(), sp[4], gp[11], i, n, K);
      mv_add(u.data(), sp[3], gp[12], i, n, K);
      mv_add(u.data(), sp[2], gp[13], i, n, K);
      mv_add(u.data(), sp[1], gp[14], i, n, K);
      mv_add(u.data(), sp[0], gp[15], i, n, K);
      for (int l = 0; l < K; ++l) out[i + (size_t) n * l + vblk * b] = u[l];
    }
  }
  return out;
}


// The Scaling Matrix's Derivatives of Order 5
//
// Im The information's derivatives of order 5, n x K x K x C_5.
// M The components of M of orders 1 to 4 (a list; `NULL` at
//   order 1).
// S The components of S of orders 0 to 4.
// tab The table of sub-multisets of order 5.
// Q,r The eigenbasis of the information, from [jsc_s0_cpp()].
// n,K The observations and the filtered parameters.
// d The power, 1/2 or 1.
// Returns A list with `M` and `S` of order 5.
// [[Rcpp::export]]
List jsc_s5_cpp(NumericVector Im, List M, List S, IntegerMatrix tab,
                NumericVector Q, NumericVector r, int n, int K, double d) {
  const int m = 5, F = 31;
  const int C = tab.nrow();
  const size_t blk = (size_t) n * K * K;
  std::vector<const double*> Mlow(m), Slow(m);
  for (int k = 1; k < m; ++k) Mlow[k] = REAL(M[k - 1]);
  for (int k = 0; k < m; ++k) Slow[k] = REAL(S[k]);
  NumericVector Mout = (d == 1.0) ? Im : NumericVector(Im.size());
  NumericVector Sout(Im.size());
  std::vector<double> R((size_t) K * K), acc((size_t) K * K), T1((size_t) K * K),
    T2((size_t) K * K);
  const double* mp[F + 1];
  const double* sp[F + 1];
  for (int b = 0; b < C; ++b) {
    for (int A = 1; A < F; ++A) {
      mp[A] = Mlow[popcount(A)] + blk * tab(b, A);
      sp[A] = Slow[popcount(A)] + blk * tab(b, A);
    }
    sp[0] = Slow[0];
    mp[F] = Mout.begin() + blk * b;
    if (d != 1.0) {
      const double* Ib = Im.begin() + blk * b;
      for (int i = 0; i < n; ++i) {
        for (int c = 0; c < K; ++c) {
          for (int a = 0; a < K; ++a) R[a + K * c] = Ib[i + (size_t) n * (a + K * c)];
        }
        mm_sub(R.data(), mp[1], mp[30], i, n, K);
        mm_sub(R.data(), mp[2], mp[29], i, n, K);
        mm_sub(R.data(), mp[3], mp[28], i, n, K);
        mm_sub(R.data(), mp[4], mp[27], i, n, K);
        mm_sub(R.data(), mp[5], mp[26], i, n, K);
        mm_sub(R.data(), mp[6], mp[25], i, n, K);
        mm_sub(R.data(), mp[7], mp[24], i, n, K);
        mm_sub(R.data(), mp[8], mp[23], i, n, K);
        mm_sub(R.data(), mp[9], mp[22], i, n, K);
        mm_sub(R.data(), mp[10], mp[21], i, n, K);
        mm_sub(R.data(), mp[11], mp[20], i, n, K);
        mm_sub(R.data(), mp[12], mp[19], i, n, K);
        mm_sub(R.data(), mp[13], mp[18], i, n, K);
        mm_sub(R.data(), mp[14], mp[17], i, n, K);
        mm_sub(R.data(), mp[15], mp[16], i, n, K);
        mm_sub(R.data(), mp[16], mp[15], i, n, K);
        mm_sub(R.data(), mp[17], mp[14], i, n, K);
        mm_sub(R.data(), mp[18], mp[13], i, n, K);
        mm_sub(R.data(), mp[19], mp[12], i, n, K);
        mm_sub(R.data(), mp[20], mp[11], i, n, K);
        mm_sub(R.data(), mp[21], mp[10], i, n, K);
        mm_sub(R.data(), mp[22], mp[9], i, n, K);
        mm_sub(R.data(), mp[23], mp[8], i, n, K);
        mm_sub(R.data(), mp[24], mp[7], i, n, K);
        mm_sub(R.data(), mp[25], mp[6], i, n, K);
        mm_sub(R.data(), mp[26], mp[5], i, n, K);
        mm_sub(R.data(), mp[27], mp[4], i, n, K);
        mm_sub(R.data(), mp[28], mp[3], i, n, K);
        mm_sub(R.data(), mp[29], mp[2], i, n, K);
        mm_sub(R.data(), mp[30], mp[1], i, n, K);
        sylvester(Mout.begin() + blk * b, Q.begin(), r.begin(), R.data(), i, n, K,
                  T1, T2);
      }
    }
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < K * K; ++k) acc[k] = 0.0;
      mm_add(acc.data(), mp[1], sp[30], i, n, K);
      mm_add(acc.data(), mp[2], sp[29], i, n, K);
      mm_add(acc.data(), mp[3], sp[28], i, n, K);
      mm_add(acc.data(), mp[4], sp[27], i, n, K);
      mm_add(acc.data(), mp[5], sp[26], i, n, K);
      mm_add(acc.data(), mp[6], sp[25], i, n, K);
      mm_add(acc.data(), mp[7], sp[24], i, n, K);
      mm_add(acc.data(), mp[8], sp[23], i, n, K);
      mm_add(acc.data(), mp[9], sp[22], i, n, K);
      mm_add(acc.data(), mp[10], sp[21], i, n, K);
      mm_add(acc.data(), mp[11], sp[20], i, n, K);
      mm_add(acc.data(), mp[12], sp[19], i, n, K);
      mm_add(acc.data(), mp[13], sp[18], i, n, K);
      mm_add(acc.data(), mp[14], sp[17], i, n, K);
      mm_add(acc.data(), mp[15], sp[16], i, n, K);
      mm_add(acc.data(), mp[16], sp[15], i, n, K);
      mm_add(acc.data(), mp[17], sp[14], i, n, K);
      mm_add(acc.data(), mp[18], sp[13], i, n, K);
      mm_add(acc.data(), mp[19], sp[12], i, n, K);
      mm_add(acc.data(), mp[20], sp[11], i, n, K);
      mm_add(acc.data(), mp[21], sp[10], i, n, K);
      mm_add(acc.data(), mp[22], sp[9], i, n, K);
      mm_add(acc.data(), mp[23], sp[8], i, n, K);
      mm_add(acc.data(), mp[24], sp[7], i, n, K);
      mm_add(acc.data(), mp[25], sp[6], i, n, K);
      mm_add(acc.data(), mp[26], sp[5], i, n, K);
      mm_add(acc.data(), mp[27], sp[4], i, n, K);
      mm_add(acc.data(), mp[28], sp[3], i, n, K);
      mm_add(acc.data(), mp[29], sp[2], i, n, K);
      mm_add(acc.data(), mp[30], sp[1], i, n, K);
      mm_add(acc.data(), mp[31], sp[0], i, n, K);
      neg_left(Sout.begin() + blk * b, Slow[0], acc.data(), i, n, K);
    }
  }
  return List::create(_["M"] = Mout, _["S"] = Sout);
}


// The Scaled Driving Vector's Derivatives of Order 5
//
// S The components of S of orders 0 to 5.
// g The family's derivatives, by order 0 to 5 of the multiset,
//   each n x K x C_k.
// tab The table of sub-multisets of order 5.
// n,K The observations and the filtered parameters.
// Returns The derivative for every multiset of order 5, n x K x C_5.
// [[Rcpp::export]]
NumericVector jsc_u5_cpp(List S, List g, IntegerMatrix tab, int n, int K) {
  const int m = 5, F = 31;
  const int C = tab.nrow();
  const size_t blk = (size_t) n * K * K, vblk = (size_t) n * K;
  std::vector<const double*> Sk(m + 1), gk(m + 1);
  for (int k = 0; k <= m; ++k) {
    Sk[k] = REAL(S[k]);
    gk[k] = REAL(g[k]);
  }
  NumericVector out(vblk * C);
  std::vector<double> u(K);
  const double* sp[F + 1];
  const double* gp[F + 1];
  for (int b = 0; b < C; ++b) {
    for (int A = 0; A <= F; ++A) {
      sp[A] = Sk[popcount(A)] + blk * tab(b, A);
      gp[A] = gk[popcount(A)] + vblk * tab(b, A);
    }
    for (int i = 0; i < n; ++i) {
      for (int l = 0; l < K; ++l) u[l] = 0.0;
      mv_add(u.data(), sp[31], gp[0], i, n, K);
      mv_add(u.data(), sp[30], gp[1], i, n, K);
      mv_add(u.data(), sp[29], gp[2], i, n, K);
      mv_add(u.data(), sp[28], gp[3], i, n, K);
      mv_add(u.data(), sp[27], gp[4], i, n, K);
      mv_add(u.data(), sp[26], gp[5], i, n, K);
      mv_add(u.data(), sp[25], gp[6], i, n, K);
      mv_add(u.data(), sp[24], gp[7], i, n, K);
      mv_add(u.data(), sp[23], gp[8], i, n, K);
      mv_add(u.data(), sp[22], gp[9], i, n, K);
      mv_add(u.data(), sp[21], gp[10], i, n, K);
      mv_add(u.data(), sp[20], gp[11], i, n, K);
      mv_add(u.data(), sp[19], gp[12], i, n, K);
      mv_add(u.data(), sp[18], gp[13], i, n, K);
      mv_add(u.data(), sp[17], gp[14], i, n, K);
      mv_add(u.data(), sp[16], gp[15], i, n, K);
      mv_add(u.data(), sp[15], gp[16], i, n, K);
      mv_add(u.data(), sp[14], gp[17], i, n, K);
      mv_add(u.data(), sp[13], gp[18], i, n, K);
      mv_add(u.data(), sp[12], gp[19], i, n, K);
      mv_add(u.data(), sp[11], gp[20], i, n, K);
      mv_add(u.data(), sp[10], gp[21], i, n, K);
      mv_add(u.data(), sp[9], gp[22], i, n, K);
      mv_add(u.data(), sp[8], gp[23], i, n, K);
      mv_add(u.data(), sp[7], gp[24], i, n, K);
      mv_add(u.data(), sp[6], gp[25], i, n, K);
      mv_add(u.data(), sp[5], gp[26], i, n, K);
      mv_add(u.data(), sp[4], gp[27], i, n, K);
      mv_add(u.data(), sp[3], gp[28], i, n, K);
      mv_add(u.data(), sp[2], gp[29], i, n, K);
      mv_add(u.data(), sp[1], gp[30], i, n, K);
      mv_add(u.data(), sp[0], gp[31], i, n, K);
      for (int l = 0; l < K; ++l) out[i + (size_t) n * l + vblk * b] = u[l];
    }
  }
  return out;
}
