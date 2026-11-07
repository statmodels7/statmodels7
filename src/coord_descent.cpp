#include <Rcpp.h>
#include "coord_descent.h"
using namespace Rcpp;

// [[Rcpp::export]]
List coord_descent(NumericMatrix X, NumericVector z, NumericVector w,
                   NumericVector beta0, NumericMatrix cut, NumericMatrix slope,
                   NumericMatrix icept, IntegerVector screen, int maxit,
                   double tol, bool covariance,
                   NumericVector means = NumericVector::create()) {
  return sm7cd::coord_run(sm7cd::DenseCols(X), z, w, beta0, cut, slope, icept, screen,
                   maxit, tol, covariance, means);
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
                          bool covariance,
                          NumericVector means = NumericVector::create()) {
  return sm7cd::coord_run(sm7cd::SparseCols(Ai, Ap, Ax, nrow, ncol), z, w, beta0, cut,
                   slope, icept, screen, maxit, tol, covariance, means);
}
