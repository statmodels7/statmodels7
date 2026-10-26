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
#' three.
#'
#' `SolveQR` has the properties `tol` and `floor`, `SolveSVD` has `tol` and
#' `floor`, and `SolveChol` has `floor`. The raw constructors do not validate
#' them; the functions [solve_qr()], [solve_svd()] and [solve_chol()] do.
#'
#' @param tol For `SolveQR`, the tolerance of the pivoted QR, a single positive
#'   number. For `SolveSVD`, `NULL` or the relative cut of the singular values,
#'   a single positive number.
#' @param floor The relative floor of the eigenvalues of a penalized matrix
#'   that is not positive definite, a single positive number.
#'
#' @return An S7 object of the class named, which inherits from `IwlsSolver`.
#'
#' @examples
#' S7::S7_inherits(solve_svd(), IwlsSolver)
#' S7::S7_inherits(solve_svd(), SolveSVD)
#' solve_svd()@tol
#'
#' @seealso [iwls_solvers] for the comparison of the three; [solve_qr()],
#'   [solve_svd()] and [solve_chol()] for the constructors.
#' @export
IwlsSolver <- S7::new_class("IwlsSolver", abstract = TRUE)

#' @rdname IwlsSolver
#' @export
SolveQR <- S7::new_class("SolveQR", parent = IwlsSolver,
  properties = list(tol = S7::class_numeric, floor = S7::class_numeric)
)

#' @rdname IwlsSolver
#' @export
SolveSVD <- S7::new_class("SolveSVD", parent = IwlsSolver,
  properties = list(
    tol = S7::new_property(S7::new_union(NULL, S7::class_numeric),
                           default = NULL),
    floor = S7::class_numeric
  )
)

#' @rdname IwlsSolver
#' @export
SolveChol <- S7::new_class("SolveChol", parent = IwlsSolver,
  properties = list(floor = S7::class_numeric)
)


#' @title Solvers Available to a Scoring Step
#' @name iwls_solvers
#'
#' @description
#' [iwls()] takes an argument `solver`, an object that selects the
#' factorization with which each scoring step is solved. Three objects are
#' available, built by [solve_qr()], [solve_svd()] and [solve_chol()]; all
#' three inherit from [IwlsSolver]. A character string is not accepted. Each
#' constructor has the numerical tolerances of its factorization as
#' arguments.
#'
#' @details
#' The step solves \eqn{(X^\top W X + S)\,\delta = u}, where \eqn{X} is the
#' design, \eqn{W} the diagonal matrix of the working weights (or, more
#' generally, the per-observation curvature), \eqn{S} the Hessian of the
#' penalty, and \eqn{u} the penalized score. The curvature is the matrix that
#' the argument `inner_information` of [statmod()] selects (see
#' [information_objects]). Let
#' \eqn{A = [\sqrt{W}X;\ \mathrm{chol}(S)]}, the design augmented by a factor of
#' the penalty, so that \eqn{A^\top A = X^\top W X + S}.
#'
#' \describe{
#'   \item{[solve_qr()]}{Factors \eqn{A} by a QR decomposition, which does not
#'     form \eqn{X^\top W X}, so the condition number of the system is not
#'     squared. Dense matrices use a QR with column pivoting, which drops a
#'     column that is numerically dependent on the others (the coefficient is
#'     reported as aliased) and records the rank; the argument `tol` is the
#'     tolerance of that decision. Sparse matrices use a sparse QR, which
#'     declines a rank-deficient system, and the step is then solved by the
#'     dense route. This is the default.}
#'   \item{[solve_svd()]}{Factors \eqn{A} by a singular value decomposition and
#'     applies the pseudo-inverse. A singular value below the cut `tol`
#'     (by default \eqn{\max(\mathrm{dim}(A))\,\epsilon\,\sigma_{\max}}, with
#'     \eqn{\epsilon} the machine precision) is discarded, so a deficient
#'     system gives a minimum-norm step and a numerical rank, and no
#'     coordinate is reported as aliased because a discarded direction is a
#'     combination of columns. It is the dearest of the three and has no
#'     sparse counterpart.}
#'   \item{[solve_chol()]}{Factors the assembled information
#'     \eqn{X^\top W X + S} by a Cholesky decomposition, after flooring the
#'     eigenvalues of a matrix that is not positive definite at the relative
#'     level `floor`. The matrix is formed explicitly, so the condition number
#'     is squared; the design of a break-point term is bounded only to
#'     \eqn{\epsilon^{-1/2}}, so the squaring can exhaust double precision.}
#' }
#'
#' [solve_qr()] and [solve_svd()] need a square root of the per-observation
#' curvature and of the penalty. Where either does not exist (an observed
#' curvature far from the optimum, a non-convex penalty), the step is solved
#' from the assembled information, as [solve_chol()] does, and the eigenvalues
#' of a matrix that is not positive definite are floored at the `floor` of
#' the solver. For this reason the three constructors have the argument
#' `floor`.
#'
#' The tolerances and their defaults are:
#'
#' \tabular{lll}{
#'   Argument \tab Meaning \tab Default \cr
#'   `tol` of [solve_qr()] \tab pivoted QR rank tolerance \tab `1e-7` \cr
#'   `tol` of [solve_svd()] \tab relative cut of the singular values \tab
#'     `NULL`, that is \eqn{\max(\mathrm{dim}(A))\,\epsilon} \cr
#'   `floor` of the three \tab relative floor of the eigenvalues \tab `1e-8` \cr
#' }
#'
#' The test of the sparse QR, which compares the equilibrated diagonal of the
#' triangular factor with \eqn{\mathrm{ncol}(A)\,\epsilon}, is not an
#' argument.
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
#' solve_svd(tol = 1e-12)
#' solve_chol(floor = 1e-10)
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
#' @param tol The tolerance of the pivoted QR (the argument `tol` of
#'   [base::qr()]), a single positive number. A column is declared dependent
#'   on the previous ones when its norm after the projection is under `tol`
#'   times its original norm. The same value is the threshold of the
#'   equilibrated diagonal under which the threaded triangular factor, used
#'   for large designs with several threads, defers to the pivoted route. The
#'   default is `1e-7`.
#' @param floor The relative floor of the eigenvalues, a single positive
#'   number, applied when the step is solved from the assembled matrix
#'   because a square root of the curvature or of the penalty does not exist.
#'   The eigenvalues of a matrix that is not positive definite are raised to
#'   `floor` times the largest absolute eigenvalue before the factorization.
#'   The default is `1e-8`.
#'
#' @return An S7 object of class `SolveQR`, which inherits from [IwlsSolver],
#'   with the properties `tol` and `floor`.
#'
#' @examples
#' solve_qr()
#' solve_qr(tol = 1e-9)
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_optimizer = iwls(solver = solve_qr(tol = 1e-9)))
#' fit@methods$smooth
#'
#' @seealso [solve_svd()] and [solve_chol()] for the other solvers;
#'   [iwls_solvers] for the comparison; [iwls()] for the function that takes the
#'   object.
#' @export
solve_qr <- function(tol = 1e-7, floor = 1e-8) {
  check_positive_number(tol, "tol")
  check_positive_number(floor, "floor")
  SolveQR(tol = as.numeric(tol), floor = as.numeric(floor))
}


#' @title The SVD Solver
#' @name solve_svd
#'
#' @description
#' Selects a singular value decomposition of the augmented design
#' \eqn{A = [\sqrt{W}X;\ \mathrm{chol}(S)]} for the scoring step of [iwls()].
#' The step applies the pseudo-inverse of \eqn{A^\top A}: a singular value of
#' \eqn{A} below the cut `tol` times the largest singular value is discarded.
#' A rank-deficient system gives a minimum-norm step and a numerical rank,
#' where a Cholesky decomposition would fail. No coordinate is reported as
#' aliased, because a discarded direction is a combination of columns. The
#' solver has no sparse counterpart, and it costs more than [solve_qr()] and
#' [solve_chol()].
#'
#' [iwls_solvers] compares it with the other two solvers.
#'
#' @param tol The relative cut of the singular values: `NULL` or a single
#'   positive number. With a number, a singular value \eqn{\sigma} is discarded
#'   when \eqn{\sigma \le \mathrm{tol}\,\sigma_{\max}}. With `NULL`, the
#'   default, the cut is \eqn{\max(\mathrm{dim}(A))\,\epsilon}, with
#'   \eqn{\epsilon} the machine precision (`.Machine$double.eps`) and
#'   \eqn{\mathrm{dim}(A)} the dimensions of the augmented design, which
#'   depend on the model.
#' @param floor The relative floor of the eigenvalues, as in [solve_qr()]. It
#'   is applied when the step is solved from the assembled matrix because a
#'   square root of the curvature or of the penalty does not exist. The
#'   default is `1e-8`.
#'
#' @return An S7 object of class `SolveSVD`, which inherits from [IwlsSolver],
#'   with the properties `tol` and `floor`.
#'
#' @examples
#' solve_svd()
#' solve_svd(tol = 1e-12)
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_optimizer = iwls(solver = solve_svd(tol = 1e-12)))
#' fit@methods$smooth
#'
#' @seealso [solve_qr()] and [solve_chol()] for the other solvers;
#'   [iwls_solvers] for the comparison; [iwls()] for the function that takes the
#'   object.
#' @export
solve_svd <- function(tol = NULL, floor = 1e-8) {
  if (!is.null(tol)) check_positive_number(tol, "tol")
  check_positive_number(floor, "floor")
  SolveSVD(tol = if (is.null(tol)) NULL else as.numeric(tol),
           floor = as.numeric(floor))
}


#' @title The Cholesky Solver
#' @name solve_chol
#'
#' @description
#' Selects a Cholesky decomposition of the assembled penalized information
#' \eqn{X^\top W X + S} for the scoring step of [iwls()]. The information is
#' formed explicitly, so the condition number of the system is squared; a
#' matrix that is not positive definite has its eigenvalues floored, at the
#' relative level `floor`, before it is factored. The decomposition names no
#' coordinate, so none is reported as aliased. It is the route that
#' [solve_qr()] and [solve_svd()] follow where the square root of the curvature
#' does not exist.
#'
#' [iwls_solvers] compares it with the other two solvers.
#'
#' @param floor The relative floor of the eigenvalues, a single positive
#'   number. A penalized matrix that is not positive definite is rebuilt from
#'   its eigendecomposition with every eigenvalue below `floor` times the
#'   largest absolute eigenvalue raised to that level. The default is `1e-8`.
#'
#' @return An S7 object of class `SolveChol`, which inherits from [IwlsSolver],
#'   with the property `floor`.
#'
#' @examples
#' solve_chol()
#' solve_chol(floor = 1e-10)
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(120, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(120, sd = 0.3)
#' G <- distributions7::gaussian1_distrib()
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 8)), G, dd,
#'                inner_optimizer = iwls(solver = solve_chol(floor = 1e-10)))
#' fit@methods$smooth
#'
#' @seealso [solve_qr()] and [solve_svd()] for the other solvers;
#'   [iwls_solvers] for the comparison; [iwls()] for the function that takes the
#'   object.
#' @export
solve_chol <- function(floor = 1e-8) {
  check_positive_number(floor, "floor")
  SolveChol(floor = as.numeric(floor))
}


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


#' The Tolerances of a Solver
#'
#' @description
#' `solver_tol()` returns the `tol` of a [solve_qr()] or [solve_svd()] object
#' and `NULL` for a [solve_chol()] object, which has none. `solver_floor()`
#' returns the `floor` of any of the three.
#'
#' @param solver An object inheriting from [IwlsSolver].
#'
#' @return A single number, or `NULL` for `solver_tol()` of a [solve_svd()]
#'   object left at its default and of a [solve_chol()] object.
#'
#' @seealso [iwls_solvers]
#' @keywords internal
solver_tol <- function(solver) {
  if (S7::S7_inherits(solver, SolveChol)) NULL else solver@tol
}

#' @rdname solver_tol
#' @keywords internal
solver_floor <- function(solver) solver@floor


#' Check That an Argument Is a Solver
#'
#' @description
#' Returns `x` invisibly when it inherits from [IwlsSolver] and signals an
#' error otherwise.
#'
#' @param x The value to check.
#'
#' @return `x`, invisibly.
#'
#' @seealso [iwls_solvers]
#' @keywords internal
check_solver <- function(x) {
  if (S7::S7_inherits(x, IwlsSolver)) return(invisible(x))
  stop("'solver' must be a solver object: solve_qr(), solve_svd() or solve_chol().",
       call. = FALSE)
}


#' Print a Solver
#'
#' @description
#' Writes one line naming the factorization and its tolerances.
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
             chol = "Cholesky solver"),
      if (S7::S7_inherits(x, SolveChol)) {
        sprintf(" (floor %g)", x@floor)
      } else {
        sprintf(" (tol %s, floor %g)",
                if (is.null(x@tol)) "default" else format(x@tol), x@floor)
      }, "\n", sep = "")
  invisible(x)
}
S7::method(print, IwlsSolver) <- print.IwlsSolver
