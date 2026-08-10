#' @include augmented.R
NULL

#' The Iterated Weighted Least Squares Method
#'
#' @description
#' The S7 class of the scoring step \code{\link{statmod}} uses for the smooth
#' block, carrying its own choice of curvature and of decomposition.
#'
#' @param hessian Either \code{"expected"} -- Fisher scoring -- or
#'   \code{"observed"}, which is Newton.
#' @param approx The approximation of the expected information, read only
#'   where the family has no closed form.
#' @param decomposition How the step is solved: \code{"qr"}, \code{"svd"},
#'   \code{"chol"} or \code{"chol_crossprod"}.
#' @param maxit The iteration budget.
#' @param tol The stopping tolerance on the scaled score.
#' @param step_halving The number of halvings allowed before a step is
#'   abandoned.
#'
#' @return An object of class \code{Iwls}.
#'
#' @seealso \code{\link{iwls}}
#'
#' @examples
#' S7::S7_inherits(iwls(), Iwls)
#'
#' @name Iwls-class
#' @aliases Iwls
#' @keywords internal
#' @export
Iwls <- S7::new_class("Iwls",
  properties = list(
    hessian = S7::class_character,
    approx = S7::class_character,
    decomposition = S7::class_character,
    maxit = S7::class_numeric,
    tol = S7::class_numeric,
    step_halving = S7::class_numeric
  )
)


#' Iterated Weighted Least Squares
#'
#' @description
#' The scoring step, written out, with the curvature and the decomposition
#' left to the caller.
#'
#' @details
#' The name is \code{iwls} and not \code{irls}: the re-weighting is already
#' implied by the iteration, so the second \code{r} says nothing that
#' \emph{iterated} has not said.
#'
#' \strong{The curvature.} \code{hessian = "expected"} inverts the expected
#' information, which is Fisher scoring; \code{"observed"} inverts the
#' observed Hessian, which is Newton. \code{approx} is handed to
#' \pkg{distributions7} and is read only where the family has no closed
#' expected information; asking for one where it would be ignored is refused,
#' since an argument accepted and ignored is worse than one that errors.
#'
#' \strong{The decomposition.} The step solves
#' \eqn{(X'WX + S)\delta = X'g - S\beta}, and how depends on this argument.
#'
#' \describe{
#'   \item{\code{"qr"}}{a QR of the augmented design
#'     \eqn{[\sqrt{W}X;\ \mathrm{chol}(S)]}. The default, and not as a matter
#'     of taste: forming \eqn{X'X} squares a conditioning that the
#'     break-point terms bound to \eqn{\epsilon^{-1/2}} by construction, and
#'     it is why a closed-form ridge is beaten at \eqn{p = 200} by a method
#'     that never forms it.}
#'   \item{\code{"svd"}}{a singular value decomposition of the same augmented
#'     design, which reports the numerical rank rather than failing on a
#'     deficient block.}
#'   \item{\code{"chol"}}{a Cholesky factor of the penalized information,
#'     assembled without the augmentation.}
#'   \item{\code{"chol_crossprod"}}{the same, forming \eqn{X'WX} explicitly.
#'     The fastest per iteration and the worst conditioned; offered because
#'     the choice is the user's.}
#' }
#'
#' @inheritParams Iwls-class
#'
#' @return An object of class \code{\link{Iwls}}.
#'
#' @seealso \code{\link{statmod}}
#'
#' @examples
#' iwls()
#' iwls(hessian = "observed", decomposition = "svd")
#'
#' @export
iwls <- function(hessian = c("expected", "observed"),
                 approx = c("bartlett", "integrate", "mc", "opg"),
                 decomposition = c("qr", "svd", "chol", "chol_crossprod"),
                 maxit = 100L, tol = 1e-6, step_halving = 30L) {
  hessian <- match.arg(hessian)
  approx <- match.arg(approx)
  decomposition <- match.arg(decomposition)
  if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1) {
    stop("'maxit' must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || tol <= 0) {
    stop("'tol' must be a single positive number.", call. = FALSE)
  }
  Iwls(hessian = hessian, approx = approx, decomposition = decomposition,
       maxit = as.numeric(maxit), tol = tol,
       step_halving = as.numeric(step_halving))
}

#' @export
#' @param x An \code{\link{Iwls}} object.
#' @param ... Unused.
#' @rdname iwls
print.Iwls <- function(x, ...) {
  cat(sprintf("iwls: %s information, %s\n", x@hessian, x@decomposition))
  cat(sprintf("  maxit %d, tol %g\n", as.integer(x@maxit), x@tol))
  invisible(x)
}
S7::method(print, Iwls) <- print.Iwls


#' Solve One Weighted Least Squares Step
#'
#' @description
#' Returns the increment \eqn{\delta} solving \eqn{(H + S)\delta = u} by the
#' requested decomposition.
#'
#' @details
#' \code{"qr"} and \code{"svd"} decompose the augmented matrix
#' \eqn{[R;\ C]}, whose cross-product IS the penalized information, so
#' \eqn{X'X} is never formed and the conditioning is never squared. That route
#' needs the per-observation curvature to be positive definite and the penalty
#' positive semidefinite; where either fails -- an observed Hessian far from
#' the optimum, a non-convex penalty -- there is no square root to take, and
#' the step falls back to the assembled matrix with its eigenvalues floored,
#' reporting in \code{route} that it did.
#'
#' \code{"chol"} and \code{"chol_crossprod"} factor the assembled information
#' directly: faster per iteration and worse conditioned, which is the trade
#' the caller is choosing.
#'
#' @param pieces A list with \code{R}, \code{C} and \code{A}, as
#'   \code{\link{iwls_pieces}} builds it.
#' @param u The right-hand side.
#' @param how The decomposition.
#'
#' @return A list with \code{delta}, \code{rank} and \code{route}.
#'
#' @keywords internal
iwls_solve <- function(pieces, u, how) {
  if (how %in% c("qr", "svd") && !is.null(pieces$R) && !is.null(pieces$C)) {
    out <- augmented_solve(pieces$R, pieces$C, u, how)
    out$route <- how
    return(out)
  }
  A <- pd_repair(pieces$A)
  route <- if (how %in% c("qr", "svd")) paste0(how, "->chol") else how
  ch <- chol(A)
  list(delta = as.numeric(backsolve(ch, forwardsolve(t(ch), u))),
       rank = ncol(A), route = route)
}


#' The Pieces One Scoring Step Needs
#'
#' @description
#' Assembles, at the current coefficients, whichever of the square-root
#' design, the penalty's factor and the penalized information the requested
#' decomposition will actually use.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients, per parameter.
#' @param hyper The hyperparameters.
#' @param method An \code{\link{Iwls}} object.
#'
#' @return A list with \code{R}, \code{C} and \code{A}.
#'
#' @keywords internal
iwls_pieces <- function(spec, design, coef, hyper, method) {
  expected <- identical(method@hessian, "expected")
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  assembled <- function() {
    list(R = NULL, C = NULL,
         A = statmod_information_at(spec, coef, design, expected,
                                    method@approx) + S)
  }
  if (!(method@decomposition %in% c("qr", "svd"))) return(assembled())
  th <- statmod_eta(spec, design, coef)$theta
  L <- chol_blocks(info_blocks(spec, th, expected, method@approx))
  R <- sqrt_design(design, L)
  C <- penalty_sqrt(S)
  if (is.null(R) || is.null(C)) return(assembled())
  list(R = R, C = C, A = NULL)
}


#' Fit the Smooth Block by Iterated Weighted Least Squares
#'
#' @description
#' Runs the scoring iteration on the objective of
#' \code{\link{statmod_objective}} at fixed hyperparameters, returning the
#' coefficients and the history of the run.
#'
#' @details
#' Two things the loop does that a plain Newton iteration does not, both
#' recorded elsewhere in the toolkit as the reason a run reported failure at
#' the answer. A non-positive-definite curvature is repaired by flooring its
#' eigenvalues rather than abandoning the start, since \code{solve()} would
#' otherwise force one. And the stopping rule is read at the ITERATE, on a
#' score scaled by the sample size, so that a threshold means the same thing
#' at \eqn{n = 10} and at \eqn{n = 10^7} while the objective itself stays
#' unaveraged.
#'
#' @param obj The objective, from \code{\link{statmod_objective}}.
#' @param start The starting coefficients, stacked.
#' @param method An \code{\link{Iwls}} object.
#' @param n The number of observations, for the scaled stopping rule.
#' @param pieces_at A function of the stacked coefficients returning what
#'   \code{\link{iwls_solve}} needs, as \code{\link{iwls_pieces}} builds it.
#' @param verbose Whether to print a line per iteration.
#'
#' @return A list with \code{par}, \code{value}, \code{converged},
#'   \code{iterations}, \code{score} and \code{history}.
#'
#' @seealso \code{\link{iwls}}
#'
#' @keywords internal
iwls_fit <- function(obj, start, method, n, pieces_at, verbose = FALSE) {
  beta <- start
  value <- obj$fn(beta)
  hist <- list()
  converged <- FALSE
  it <- 0L
  score <- Inf

  if (verbose) {
    cat(sprintf("  %-5s %14s %12s %10s\n", "iter", "objective", "score/n",
                "step"))
  }
  for (it in seq_len(as.integer(method@maxit))) {
    g <- obj$gr(beta)
    score <- max(abs(g)) / n
    if (verbose) {
      cat(sprintf("  %-5d %14.6f %12.3e %10s\n", it - 1L, value, score,
                  if (it == 1L) "-" else fmt_step(step_used)))
    }
    if (score < method@tol) {
      converged <- TRUE
      break
    }
    sol <- iwls_solve(pieces_at(beta), -g, method@decomposition)
    delta <- sol$delta
    delta[!is.finite(delta)] <- 0

    # sufficient decrease, not mere non-increase
    step_used <- 1
    ok <- FALSE
    for (h in seq_len(as.integer(method@step_halving))) {
      cand <- beta + step_used * delta
      vnew <- obj$fn(cand)
      if (is.finite(vnew) && vnew <= value - 1e-4 * step_used * sum(g * delta)) {
        ok <- TRUE
        break
      }
      step_used <- step_used / 2
    }
    if (!ok) break
    hist[[length(hist) + 1L]] <- data.frame(
      iteration = it, objective = vnew, score = score, step = step_used,
      rank = sol$rank, route = sol$route
    )
    beta <- cand
    value <- vnew
  }
  # the rule is asked once more at the point reached, so that a run which
  # could not move is still converged when the POINT says so
  if (!converged) {
    score <- max(abs(obj$gr(beta))) / n
    converged <- score < method@tol
  }
  if (verbose) {
    cat(sprintf("  %-5s %14.6f %12.3e %10s\n", "end", value, score,
                if (converged) "converged" else "stopped"))
  }
  list(par = beta, value = value, converged = converged,
       iterations = it, score = score,
       history = if (length(hist)) do.call(rbind, hist) else NULL)
}


#' Floor the Eigenvalues of a Curvature Matrix
#'
#' @description
#' Returns the matrix unchanged when it is positive definite, and otherwise
#' its eigendecomposition with the eigenvalues floored.
#'
#' @details
#' Abandoning a start because the curvature is indefinite is what
#' \code{solve()} forces and what a repair avoids; the floor is relative to
#' the largest eigenvalue, since an absolute one means nothing across scales.
#'
#' @param A A symmetric matrix.
#' @param rel The floor, relative to the largest eigenvalue.
#'
#' @return A symmetric positive definite matrix.
#'
#' @keywords internal
pd_repair <- function(A, rel = 1e-8) {
  if (nrow(A) == 0L) return(A)
  ok <- tryCatch({
    chol(A)
    TRUE
  }, error = function(e) FALSE)
  if (ok) return(A)
  e <- eigen(A, symmetric = TRUE)
  floor_ <- rel * max(abs(e$values))
  if (!is.finite(floor_) || floor_ <= 0) floor_ <- rel
  e$vectors %*% (pmax(e$values, floor_) * t(e$vectors))
}


#' Format a Step Length
#'
#' @param x A number.
#'
#' @return A single string.
#'
#' @keywords internal
fmt_step <- function(x) {
  if (is.null(x) || !is.finite(x)) return("-")
  formatC(x, format = "g", digits = 3)
}
