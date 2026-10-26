#' @include information.R
NULL

#' @title S7 Classes for the Solver of a Scoring Step
#' @name IwlsSolver
#' @aliases SolveQR SolveSVD SolveChol
#'
#' @description
#' `IwlsSolver` is the abstract parent of the objects that select how the linear
#' system of one step of [iwls()] is solved. Its three subclasses are `SolveQR`,
#' `SolveSVD` and `SolveChol`. They are built by [solve_qr()], [solve_svd()] and
#' [solve_chol()], which are the functions to call; the classes are exported so
#' that `S7::S7_inherits()` can test an object. [iwls_solvers] compares the
#' three. The three classes have no properties.
#'
#' @return An S7 object of the class named, which inherits from `IwlsSolver`.
#'
#' @examples
#' S7::S7_inherits(solve_svd(), IwlsSolver)
#' S7::S7_inherits(solve_svd(), SolveSVD)
#'
#' @seealso [iwls_solvers] for the comparison of the three; [solve_qr()],
#'   [solve_svd()] and [solve_chol()] for the constructors.
#' @export
IwlsSolver <- S7::new_class("IwlsSolver", abstract = TRUE)

#' @rdname IwlsSolver
#' @export
SolveQR <- S7::new_class("SolveQR", parent = IwlsSolver)

#' @rdname IwlsSolver
#' @export
SolveSVD <- S7::new_class("SolveSVD", parent = IwlsSolver)

#' @rdname IwlsSolver
#' @export
SolveChol <- S7::new_class("SolveChol", parent = IwlsSolver)


#' @title Solvers Available to a Scoring Step
#' @name iwls_solvers
#'
#' @description
#' [iwls()] takes an argument `solver`, an object that selects the
#' factorization with which each scoring step is solved. Three objects are
#' available, built by [solve_qr()], [solve_svd()] and [solve_chol()]; all
#' three inherit from [IwlsSolver]. A character string is not accepted.
#'
#' @details
#' The step solves \eqn{(X^\top W X + S)\,\delta = u}, where \eqn{X} is the
#' design, \eqn{W} the diagonal matrix of the working weights (or, more
#' generally, the per-observation curvature), \eqn{S} the Hessian of the
#' penalty, and \eqn{u} the penalized score. Let
#' \eqn{A = [\sqrt{W}X;\ \mathrm{chol}(S)]}, the design augmented by a factor of
#' the penalty, so that \eqn{A^\top A = X^\top W X + S}.
#'
#' \describe{
#'   \item{[solve_qr()]}{Factors \eqn{A} by a QR decomposition, which does not
#'     form \eqn{X^\top W X}, so the condition number of the system is not
#'     squared. Dense matrices use a QR with column pivoting, which drops a
#'     column that is numerically dependent on the others (the coefficient is
#'     reported as aliased) and records the rank. Sparse matrices use a sparse
#'     QR, which declines a rank-deficient system, and the step is then solved
#'     by the dense route. This is the default.}
#'   \item{[solve_svd()]}{Factors \eqn{A} by a singular value decomposition and
#'     applies the pseudo-inverse. A singular value below
#'     \eqn{\max(\mathrm{dim}(A))\,\epsilon\,\sigma_{\max}} (with \eqn{\epsilon}
#'     the machine precision) is discarded, so a deficient system gives a
#'     minimum-norm step and a numerical rank, and no coordinate is reported as
#'     aliased because a discarded direction is a combination of columns. It is
#'     the dearest of the three and has no sparse counterpart.}
#'   \item{[solve_chol()]}{Factors the assembled information
#'     \eqn{X^\top W X + S} by a Cholesky decomposition, after flooring the
#'     eigenvalues of a matrix that is not positive definite. The matrix is
#'     formed explicitly, so the condition number is squared; the design of a
#'     break-point term is bounded only to \eqn{\epsilon^{-1/2}}, so the
#'     squaring can exhaust double precision.}
#' }
#'
#' [solve_qr()] and [solve_svd()] need a square root of the per-observation
#' curvature and of the penalty. Where either does not exist (an observed
#' curvature far from the optimum, a non-convex penalty), the step is solved
#' from the assembled information, as [solve_chol()] does.
#'
#' Four tolerances are fixed inside the solvers and are not arguments of the
#' constructors: the rank tolerance of the pivoted QR, the singular-value
#' tolerance of the SVD above, the relative floor \eqn{10^{-8}} of the repair of
#' an indefinite matrix, and the equilibrated-diagonal test of the sparse QR.
#'
#' @return No value. The topic describes the three objects, of classes
#'   [SolveQR][IwlsSolver], [SolveSVD][IwlsSolver] and [SolveChol][IwlsSolver],
#'   that [solve_qr()], [solve_svd()] and [solve_chol()] return.
#'
#' @seealso [solve_qr()], [solve_svd()] and [solve_chol()] for the
#'   constructors; [iwls()] for the function that takes the object;
#'   [information_objects] for the matrix that the solver factors.
#'
#' @examples
#' solve_qr()
#' solve_svd()
#' solve_chol()
#'
#' # The three reach the same fit on a design of full rank.
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' f_qr <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                 inner_optimizer = iwls(solver = solve_qr()),
#'                 outer_criterion = NULL)
#' f_chol <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                   inner_optimizer = iwls(solver = solve_chol()),
#'                   outer_criterion = NULL)
#' max(abs(unlist(coef(f_qr)) - unlist(coef(f_chol))))
NULL


#' @title The QR Solver
#' @name solve_qr
#'
#' @description
#' Selects a QR decomposition of the augmented design
#' \eqn{A = [\sqrt{W}X;\ \mathrm{chol}(S)]} for the scoring step of [iwls()].
#' The cross-product \eqn{X^\top W X} is never formed, so the condition number
#' of the system is not squared. A dense matrix is factored with column
#' pivoting: a column that is numerically dependent on the others is dropped
#' from the step, the coefficient is reported as aliased, and the rank of the
#' system is recorded. A sparse matrix is factored by a sparse QR, which
#' declines a rank-deficient system and leaves the step to the dense route.
#'
#' It is the default of [iwls()]. Choose it for ordinary designs and for a
#' break-point term, whose design is ill-conditioned. [iwls_solvers] compares it
#' with the other two solvers.
#'
#' @return An S7 object of class `SolveQR`, which inherits from [IwlsSolver],
#'   with no properties.
#'
#' @examples
#' solve_qr()
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_optimizer = iwls(solver = solve_qr()))
#' fit@methods$smooth
#'
#' @seealso [solve_svd()] and [solve_chol()] for the other solvers;
#'   [iwls_solvers] for the comparison; [iwls()] for the function that takes the
#'   object.
#' @export
solve_qr <- function() SolveQR()


#' @title The SVD Solver
#' @name solve_svd
#'
#' @description
#' Selects a singular value decomposition of the augmented design
#' \eqn{A = [\sqrt{W}X;\ \mathrm{chol}(S)]} for the scoring step of [iwls()].
#' The step applies the pseudo-inverse of \eqn{A^\top A}: a singular value of
#' \eqn{A} below \eqn{\max(\mathrm{dim}(A))\,\epsilon\,\sigma_{\max}}, with
#' \eqn{\epsilon} the machine precision, is discarded. A rank-deficient system
#' gives a minimum-norm step and a numerical rank, where a Cholesky
#' decomposition would fail. No coordinate is reported as aliased, because a
#' discarded direction is a combination of columns. The solver has no sparse
#' counterpart, and it costs more than [solve_qr()] and [solve_chol()].
#'
#' [iwls_solvers] compares it with the other two solvers.
#'
#' @return An S7 object of class `SolveSVD`, which inherits from [IwlsSolver],
#'   with no properties.
#'
#' @examples
#' solve_svd()
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_optimizer = iwls(solver = solve_svd()))
#' fit@methods$smooth
#'
#' @seealso [solve_qr()] and [solve_chol()] for the other solvers;
#'   [iwls_solvers] for the comparison; [iwls()] for the function that takes the
#'   object.
#' @export
solve_svd <- function() SolveSVD()


#' @title The Cholesky Solver
#' @name solve_chol
#'
#' @description
#' Selects a Cholesky decomposition of the assembled penalized information
#' \eqn{X^\top W X + S} for the scoring step of [iwls()]. The information is
#' formed explicitly, so the condition number of the system is squared; a
#' matrix that is not positive definite has its eigenvalues floored, at a
#' relative \eqn{10^{-8}}, before it is factored. The decomposition names no
#' coordinate, so none is reported as aliased. It is the route that
#' [solve_qr()] and [solve_svd()] follow where the square root of the curvature
#' does not exist.
#'
#' [iwls_solvers] compares it with the other two solvers.
#'
#' @return An S7 object of class `SolveChol`, which inherits from [IwlsSolver],
#'   with no properties.
#'
#' @examples
#' solve_chol()
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_optimizer = iwls(solver = solve_chol()))
#' fit@methods$smooth
#'
#' @seealso [solve_qr()] and [solve_svd()] for the other solvers;
#'   [iwls_solvers] for the comparison; [iwls()] for the function that takes the
#'   object.
#' @export
solve_chol <- function() SolveChol()


#' The Name of a Solver
#'
#' @description
#' Internal generic that returns the word with which the numerical routines of
#' the step (`iwls_solve()`, `augmented_solve()`) identify a solver: `"qr"`,
#' `"svd"` or `"chol"`. It dispatches on the class of `solver`.
#'
#' @param solver An object inheriting from [IwlsSolver].
#'
#' @return A single string.
#'
#' @seealso [iwls_solvers]
#' @keywords internal
solver_how <- S7::new_generic("solver_how", "solver",
  function(solver) S7::S7_dispatch())

S7::method(solver_how, SolveQR) <- function(solver) "qr"
S7::method(solver_how, SolveSVD) <- function(solver) "svd"
S7::method(solver_how, SolveChol) <- function(solver) "chol"


#' Check That an Argument Is a Solver
#'
#' @description
#' Returns `x` invisibly when it inherits from [IwlsSolver] and signals an
#' error otherwise. A character string of the former interface (`"qr"`,
#' `"svd"` or `"chol"`) produces a message that names the constructor that
#' replaces it.
#'
#' @param x The value to check.
#'
#' @return `x`, invisibly.
#'
#' @seealso [iwls_solvers]
#' @keywords internal
check_solver <- function(x) {
  if (S7::S7_inherits(x, IwlsSolver)) return(invisible(x))
  repl <- c(qr = "solve_qr()", svd = "solve_svd()", chol = "solve_chol()")
  if (is.character(x) && length(x) >= 1L && !anyNA(x) && x[[1L]] %in% names(repl)) {
    stop(sprintf(paste0("'solver' must be a solver object such as ",
                        "solve_qr(); the string \"%s\" is no longer ",
                        "accepted. Use %s."), x[[1L]], repl[[x[[1L]]]]),
         call. = FALSE)
  }
  stop("'solver' must be a solver object: solve_qr(), solve_svd() or solve_chol().",
       call. = FALSE)
}


#' Print a Solver
#'
#' @description
#' Writes one line naming the factorization.
#'
#' @param x An object inheriting from [IwlsSolver].
#' @param ... Unused.
#'
#' @return `x`, invisibly. The function is called for the output that it writes.
#'
#' @examples
#' solve_qr()
#'
#' @seealso [iwls_solvers]
#' @name print.IwlsSolver
#' @keywords internal
print.IwlsSolver <- function(x, ...) {
  cat(switch(solver_how(x), qr = "QR solver", svd = "SVD solver",
             chol = "Cholesky solver"), "\n", sep = "")
  invisible(x)
}
S7::method(print, IwlsSolver) <- print.IwlsSolver
