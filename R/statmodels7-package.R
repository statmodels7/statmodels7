#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib statmodels7, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

# A design carrying a grouping indicator is SPARSE, and the operations that
# assemble a fit -- the transpose, the crossproduct, the solve, the Cholesky
# -- are S4 generics in Matrix rather than the base functions of the same
# name. A package that only Imports Matrix does not get that dispatch: base's
# `t()` reaches `t.default()` and reports that its argument is not a matrix,
# from inside the assembly and naming nothing a caller wrote. Importing the
# generics is what Matrix documents for this, and it changes nothing for a
# dense design, each generic falling back to the base function there.
#' @importFrom Matrix t crossprod tcrossprod solve chol diag rowSums colSums
NULL

# S7 methods registered on a BASE generic -- print, summary, coef, vcov,
# confint, predict, fitted, logLik -- reach S3's dispatch table only when this
# runs. Under pkgload they are registered as a side effect of evaluating the
# package's code, so a suite run with test_local() is green either way; from an
# INSTALLED package they are simply absent, and a fit prints as its raw
# property dump while confint() fails inside the method it never reached.
#' @noRd
.onLoad <- function(libname, pkgname) {
  S7::methods_register()
}

# The members sit in Imports rather than Depends. Depends would attach them
# through R's own mechanism, in the order the field happens to list them and
# with no message, which leaves a caller no way to see what was attached or at
# which version. Imports pulls them all at installation, which is what makes
# one install command install the toolkit, and the attaching is done here so
# that it can be reported.
#
# The names are read out of the package's own DESCRIPTION rather than written
# out a second time here: a list repeated in two places is a list that will
# disagree with itself.
#
# R CMD check reports "All declared Imports should be used" for every member,
# and the note is correct as far as its heuristic can see: nothing here calls
# any member by name. It is left standing rather than silenced. The Imports
# field is what makes installing this package install the toolkit, which is
# the package's whole purpose, and the alternatives are worse -- importing an
# arbitrary symbol from each member to satisfy the check would put five
# unused bindings in the namespace and say something false about what the
# package uses, while moving the members to Depends would give up control of
# the attaching that this package exists to exercise. The note is a NOTE, and
# the CI treats only warnings as failures.
