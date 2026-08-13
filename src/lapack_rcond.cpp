// The reciprocal condition number of a symmetric positive definite matrix,
// estimated from a Cholesky factor that has already been computed.
//
// This is LAPACK's dpocon, and what it buys is the difference between
// O(p^2) and O(p^3): asking whether the smallest eigenvalue of a penalized
// information is real curvature or a flat direction used to cost a full
// eigendecomposition -- 1.18 s at p = 1022 against 0.25 s for the Cholesky
// the inverse needs anyway. The estimator returns 1/(||A||_1 ||A^-1||_1),
// so the caller reads the smallest eigenvalue off it as rcond * ||A||_1;
// see solve_pd() for which way that estimate errs and why the slack does
// not reach the two cases it has to keep apart.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>

#ifndef FCONE
# define FCONE
#endif

using namespace Rcpp;

// [[Rcpp::export]]
double chol_rcond_cpp(NumericMatrix R, double anorm) {
    int n = R.nrow();
    if (n <= 0 || !R_FINITE(anorm) || anorm <= 0) return NA_REAL;
    // dpocon takes the factor as dpotrf leaves it, which is what chol()
    // returns for uplo = "U": A = R'R with R upper triangular
    std::vector<double> a(R.begin(), R.end());
    std::vector<double> work(3 * (size_t) n);
    std::vector<int> iwork((size_t) n);
    double rcond = 0.0;
    int info = 0;
    F77_CALL(dpocon)("U", &n, a.data(), &n, &anorm, &rcond, work.data(),
                     iwork.data(), &info FCONE);
    if (info != 0) return NA_REAL;
    return rcond;
}
