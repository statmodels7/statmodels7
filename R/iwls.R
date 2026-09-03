#' @include augmented.R
NULL

#' The Iterated Weighted Least Squares Method
#'
#' @description
#' The S7 class of the scoring step [statmod()] uses for the jointly fitted
#' smooth block. An object of this class carries the whole configuration of
#' that step: which curvature to invert, how to solve the system, how long to
#' run and when to stop. Build one with [iwls()], which validates every
#' argument; this raw constructor does not.
#'
#' @param hessian `"expected"` for Fisher scoring or `"observed"` for
#'   Newton, a single string.
#' @param approx How the expected information is approximated for a family
#'   with no closed form: `"opg"`, `"bartlett"`, `"integrate"` or `"mc"`.
#' @param decomposition How the step is solved: `"qr"`, `"svd"`, `"chol"` or
#'   `"chol_crossprod"`.
#' @param maxit The iteration budget, a single positive number.
#' @param tol The stopping tolerance on the score per observation.
#' @param criterion An \pkg{optimizers7} `criterion` object driving the loop
#'   in place of `tol`, or `NULL` for the built-in rule.
#' @param step_halving How many halvings are allowed before a step is
#'   abandoned, a single non-negative number.
#'
#' @return An object of class `Iwls` with one property per argument above,
#'   each holding what was passed.
#'
#' @seealso [iwls()], the constructor to use.
#'
#' @examples
#' S7::S7_inherits(iwls(), Iwls)
#' iwls()@decomposition
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
    criterion = S7::class_any,
    step_halving = S7::class_numeric
  )
)


#' Iterated Weighted Least Squares
#'
#' @description
#' Configures the scoring step `statmod()` fits the smooth block with, and
#' returns it as an object to pass as `inner_optimizer`. This is the default,
#' so `statmod(...)` and `statmod(..., inner_optimizer = iwls())` fit the
#' same model; call it by name to change the curvature, the decomposition,
#' the budget or the stopping rule.
#'
#' Every argument has a default, and the object carries no data: one `iwls()`
#' can be passed to any number of fits.
#'
#' @details
#' # The name
#'
#' `iwls`, not `irls`. The iteration already implies the re-weighting, so the
#' second `r` says nothing *iterated* has not said.
#'
#' # The curvature
#'
#' `hessian = "expected"` inverts the expected information, which is Fisher
#' scoring, and `"observed"` inverts the observed Hessian, which is Newton.
#' The expected information is positive definite wherever the family is
#' regular, so Fisher scoring takes usable steps from a poor start where
#' Newton's matrix may be indefinite; Newton converges faster near the
#' optimum and is what the exact outer gradient needs.
#'
#' `approx` reaches \pkg{distributions7} and is read only where the family
#' has no closed expected information; elsewhere the family's own method
#' answers and the argument does not enter. It is accepted rather than
#' refused there because an `iwls()` object is built before it meets a
#' distribution and can be passed to any number of fits, so the question it
#' would be refused on cannot be asked at construction.
#' [distributions7::expected_hessian_exact()] is the predicate that says
#' which case a family is in.
#'
#' Its default is `"opg"`, the outer product of the observed scores, which
#' costs one gradient. The alternative readings evaluate the expectation
#' itself -- a sum over the support for a discrete family, a quadrature for a
#' continuous one -- and are orders of magnitude dearer inside a loop that
#' rebuilds the curvature at every iteration: measured on a Poisson-inverse
#' gaussian regression at \eqn{n = 500}, `statmod()` takes 0.64 s under
#' `"opg"` and 89.06 s under `"bartlett"`, reaching the same log-likelihood
#' and coefficients agreeing to \eqn{10^{-6}}. What the cheap matrix gives up
#' is variance, not the answer: the score is exact, so a scoring step driven
#' by any positive definite matrix reaches the same stationary point, and the
#' outer product is positive semidefinite by construction where the observed
#' Hessian need not be. A standard error is a different question and is read
#' off the observed information instead -- see [vcov.StatmodFit()].
#'
#' # The decomposition
#'
#' The step solves \eqn{(X'WX + S)\delta = X'g - S\beta}, and this argument
#' says how.
#'
#' \describe{
#'   \item{`"qr"`}{a QR of the augmented design
#'     \eqn{[\sqrt{W}X;\ \mathrm{chol}(S)]}, and the default. Forming
#'     \eqn{X'X} squares the condition number, and a break-point term's
#'     design is bounded only to \eqn{\epsilon^{-1/2}} by construction, so
#'     squaring it exhausts double precision. This is also why a closed-form
#'     ridge is beaten at \eqn{p = 200} by a method that never forms
#'     \eqn{X'X}.}
#'   \item{`"svd"`}{a singular value decomposition of the same augmented
#'     design. It reports a numerical rank where a deficient block would make
#'     a Cholesky fail.}
#'   \item{`"chol"`}{a Cholesky factor of the penalized information,
#'     assembled without the augmentation.}
#'   \item{`"chol_crossprod"`}{the same, forming \eqn{X'WX} explicitly.
#'     The fastest per iteration and the worst conditioned; offered because
#'     the choice is the user's.}
#' }
#'
#' # The stopping rule
#'
#' A scoring step is not an optimizer and carries its own loop, but the rule
#' that ends the loop belongs to the caller.
#'
#' With `criterion = NULL`, the default, the built-in rule applies: the score
#' per observation \eqn{\max_j \lvert g_j \rvert / n} against `tol`, with the
#' dimensionless reading of [iwls_score()] arbitrating the final verdict. Any
#' \pkg{optimizers7} criterion may drive the loop instead, and then `tol` is
#' never read, so passing both is an error.
#'
#' What a criterion is shown is the state `optimizers7::crit_met()`
#' documents, with three things worth knowing about it here:
#'
#' - Its `gradient` is the score **per observation**, the quantity the
#'   built-in rule compares. So `criterion = optimizers7::crit_grad(t)` is
#'   exactly `tol = t`, and either threshold means the same at \eqn{n = 10}
#'   and at \eqn{n = 10^7}.
#' - Its objective is the penalized log-likelihood **unaveraged**, which is
#'   the scale the penalty is added on. A rule reading the objective's
#'   absolute value therefore carries the sample size with it; one reading
#'   its relative change does not.
#' - On the first iteration `f_old` and `x_old` are `NULL`, there being no
#'   previous point, so a rule reading a change returns `FALSE` there by
#'   construction.
#'
#' A criterion needing something the step does not compute is rejected at
#' construction. A stationarity measure is the case: it belongs to the
#' derivative-free methods, and a rule asking for one would otherwise sit in
#' the loop and never fire.
#'
#' @inheritParams Iwls-class
#'
#' @return An [Iwls()] object, holding the seven settings above and no data.
#'
#' @seealso [statmod()] for where it is passed, [iwls_score()] for the
#'   built-in rule's final verdict, [optimizers7::crit_grad()] for the
#'   criteria that can replace it.
#'
#' @examples
#' iwls()
#'
#' # Newton instead of Fisher scoring, solved through an SVD, which reports a
#' # numerical rank where a deficient block would make a Cholesky fail.
#' iwls(hessian = "observed", decomposition = "svd")
#'
#' # An optimizers7 stopping rule in place of tol. crit_grad(t) here is
#' # exactly tol = t, the state's gradient being the score per observation.
#' iwls(criterion = optimizers7::crit_any(optimizers7::crit_grad(1e-8),
#'                                        optimizers7::crit_rel_obj(1e-12)))
#'
#' # Both at once is an error: tol would be read by nobody.
#' try(iwls(tol = 1e-8, criterion = optimizers7::crit_grad(1e-8)))
#'
#' @export
iwls <- function(hessian = c("expected", "observed"),
                 approx = c("opg", "bartlett", "integrate", "mc"),
                 decomposition = c("qr", "svd", "chol", "chol_crossprod"),
                 maxit = 100L, tol = 1e-6, criterion = NULL,
                 step_halving = 30L) {
  hessian <- match.arg(hessian)
  approx <- match.arg(approx)
  decomposition <- match.arg(decomposition)
  if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1) {
    stop("'maxit' must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || tol <= 0) {
    stop("'tol' must be a single positive number.", call. = FALSE)
  }
  if (!is.null(criterion)) {
    if (!S7::S7_inherits(criterion, optimizers7::criterion)) {
      stop("'criterion' must be an optimizers7 criterion, or NULL.",
           call. = FALSE)
    }
    # the step computes a gradient and nothing else; a rule asking for more
    # would never fire, which is what optimizers7's own check_criterion()
    # refuses at construction for the same reason
    missing_ <- setdiff(optimizers7::crit_needs(criterion), "gradient")
    if (length(missing_)) {
      stop(sprintf(paste0("The stopping rule needs %s, which a scoring step ",
                          "does not compute.\n  Choose a rule it can ",
                          "evaluate, or leave 'criterion' NULL."),
                   paste(missing_, collapse = ", ")), call. = FALSE)
    }
    if (!missing(tol)) {
      stop(paste0("'tol' and 'criterion' both say when to stop: pass one.\n",
                  "  The rule reads the score per observation, so ",
                  "crit_grad(tol) is what 'tol' means."), call. = FALSE)
    }
  }
  Iwls(hessian = hessian, approx = approx, decomposition = decomposition,
       maxit = as.numeric(maxit), tol = tol, criterion = criterion,
       step_halving = as.numeric(step_halving))
}

#' @export
#' @param x An [Iwls()] object.
#' @param ... Unused.
#' @rdname iwls
print.Iwls <- function(x, ...) {
  cat(sprintf("iwls: %s information, %s\n", x@hessian, x@decomposition))
  cat(sprintf("  maxit %d, %s\n", as.integer(x@maxit),
              if (is.null(x@criterion)) sprintf("tol %g", x@tol)
              else x@criterion@label))
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
#' # The two families of route
#'
#' `"qr"` and `"svd"` decompose the augmented matrix \eqn{[R;\, C]}, whose
#' cross-product is the penalized information, so \eqn{X'X} is never formed
#' and the condition number is never squared. That route needs the
#' per-observation curvature to be positive definite and the penalty positive
#' semidefinite. Where either fails, as with an observed Hessian far from the
#' optimum or a non-convex penalty, there is no square root to take: the step
#' falls back to the assembled matrix with its eigenvalues floored, and says
#' so in `route`.
#'
#' `"chol"` and `"chol_crossprod"` factor the assembled information directly.
#' Faster per iteration, worse conditioned, and the trade the caller chose.
#'
#' # The damping is Levenberg's, and is zero unless a step has failed
#'
#' `damp` adds \eqn{\lambda I} to the system, which shortens a coordinate in
#' proportion to how little curvature it has. A coordinate whose diagonal is
#' 0.24 beside neighbors at 2328 is shortened by
#' \eqn{(0.24 + \lambda)/0.24} while theirs move by
#' \eqn{(2328 + \lambda)/2328}: a differential shrink, which a scalar step
#' length cannot give.
#'
#' It goes in as further rows of the augmented design, so the QR route keeps
#' its conditioning. The matrix added is the identity. A damping
#' proportional to \eqn{\mathrm{diag}(K)} would shrink every coordinate
#' alike and leave the disparity exactly where it was.
#'
#' # A coordinate at a boundary is held
#'
#' Where a parameter has reached the clamp its link keeps it strictly inside,
#' the family's curvature there is zero or `NaN`, and no step can move that
#' coordinate. Those coordinates are dropped from the system and the rest is
#' solved, which is the active set the boundary defines; the step stays a
#' descent step for the reduced problem.
#'
#' Holding them keeps one boundary coordinate from stopping the others. With
#' a single `NaN` in the curvature the whole step came back zero, and a
#' Student t whose \eqn{\nu} had reached `double.xmax` stopped with a score
#' of -617.6 in \eqn{\sigma}, an ordinary coordinate nowhere near stationary,
#' while \eqn{\nu}'s own score was exactly 0.
#'
#' The test is on the **diagonal** of the curvature. A boundary coordinate
#' makes its entire row non-finite, cross terms included, so a test over
#' whole columns marks its neighbors too and holds coordinates that were
#' perfectly movable.
#'
#' @param pieces A list with `R`, `C` and `A`, as [iwls_pieces()] builds it.
#'   Which of the three are present depends on `how`.
#' @param u The right-hand side \eqn{X'g - S\beta}, a numeric vector of
#'   length `p`.
#' @param how The decomposition: `"qr"`, `"svd"`, `"chol"` or
#'   `"chol_crossprod"`.
#' @param damp The Levenberg damping \eqn{\lambda}, a single non-negative
#'   number. `0`, the default, is the plain scoring step.
#'
#' @return A list of four:
#'   \describe{
#'     \item{`delta`}{the increment, a numeric vector of length `p`, zero in
#'       every held coordinate.}
#'     \item{`rank`}{the rank of the system solved, an integer.}
#'     \item{`route`}{which decomposition was used, including `"floor"` when
#'       the augmented route fell back to the floored assembled matrix.}
#'     \item{`held`}{the positions dropped from the system, an integer
#'       vector, empty when nothing was held.}
#'   }
#'
#' @seealso [iwls_pieces()] for the input, [iwls()] for the `decomposition`
#'   argument this serves.
#'
#' @keywords internal
iwls_solve <- function(pieces, u, how, damp = 0) {
  p <- length(u)
  # sqrt(lambda) I as further rows of the augmented design, which is how a
  # ridge is added there: [R; C; sqrt(lambda) I] has cross-product
  # K + lambda I and the route keeps its conditioning.
  damp_rows <- function(C, k) {
    if (!(damp > 0)) return(C)
    D <- if (isS4(C)) Matrix::Diagonal(k, sqrt(damp)) else diag(sqrt(damp), k)
    rbind(C, D)
  }
  # A COORDINATE WHOSE CURVATURE IS NOT FINITE IS ONE COORDINATE, NOT THE
  # WHOLE POINT. It happens where a parameter has reached the clamp its link
  # keeps it strictly inside: the family's derivatives there are 0 or NaN,
  # and the coordinate cannot move whatever is solved. Until 0.82.0 a single
  # such entry made the whole system unusable and the step came back zero for
  # EVERY coordinate, so a run stopped with the others far from stationary --
  # measured on a Student t whose nu reached double.xmax, the score in sigma
  # was -617.6 there while nu's was exactly 0. Those coordinates are held and
  # the rest is solved, which is the active set a boundary defines: the step
  # is a descent step for the reduced problem, and the held coordinate's own
  # score being zero is what says it is a stationary point of the constrained
  # problem rather than a place the fit was stopped at.
  # `colSums(abs(.))` and not `is.finite()` over the whole matrix: a column
  # sum keeps a sparse design sparse where an elementwise test would build a
  # dense logical of its size, and `abs` cannot overflow into a false
  # positive where a square could. A NaN or an infinity anywhere in the
  # column reaches the sum.
  bad_cols <- function(M) !is.finite(as.numeric(Matrix::colSums(abs(M))))
  held <- integer(0)
  if (how %in% c("qr", "svd") && !is.null(pieces$R) && !is.null(pieces$C)) {
    bad <- !is.finite(u) | bad_cols(pieces$R) | bad_cols(pieces$C)
    if (any(bad)) {
      held <- which(bad)
      keep <- which(!bad)
      if (!length(keep)) {
        return(list(delta = numeric(p), rank = 0L, route = how, held = held))
      }
      out <- augmented_solve(pieces$R[, keep, drop = FALSE],
                             damp_rows(pieces$C[, keep, drop = FALSE],
                                       length(keep)),
                             u[keep], how,
                             threads = if (is.null(pieces$threads)) 1L
                                       else pieces$threads)
      d <- numeric(p)
      d[keep] <- out$delta
      return(list(delta = d, rank = out$rank, route = how, held = held))
    }
    out <- augmented_solve(pieces$R, damp_rows(pieces$C, p), u, how,
                           threads = if (is.null(pieces$threads)) 1L
                                     else pieces$threads)
    out$route <- how
    out$held <- held
    return(out)
  }
  route <- if (how %in% c("qr", "svd")) paste0(how, "->chol") else how
  A0 <- pieces$A
  if (!is.null(A0)) {
    # BY THE DIAGONAL AND NOT BY THE COLUMN. A boundary coordinate makes its
    # whole ROW non-finite, cross terms included, so a column test marks its
    # neighbors too: measured on the Student t, testing columns held sigma
    # along with nu and left the fit exactly where it had been. The diagonal
    # names the coordinate whose own curvature is gone; dropping its row and
    # column then clears the cross terms, and the loop repeats in case the
    # reduced matrix still carries one from another source.
    bad <- !is.finite(u) | !is.finite(as.numeric(Matrix::diag(A0)))
    repeat {
      keep <- which(!bad)
      if (!length(keep)) break
      more <- bad_cols(A0[keep, keep, drop = FALSE])
      if (!any(more)) break
      bad[keep[more]] <- TRUE
    }
    if (any(bad)) {
      held <- which(bad)
      keep <- which(!bad)
      if (!length(keep)) {
        return(list(delta = numeric(p), rank = 0L, route = route, held = held))
      }
      A <- pd_repair(A0[keep, keep, drop = FALSE])
      if (is.null(A)) {
        return(list(delta = numeric(p), rank = 0L, route = route, held = held))
      }
      if (damp > 0) Matrix::diag(A) <- Matrix::diag(A) + damp
      ch <- chol(A)
      d <- numeric(p)
      d[keep] <- as.numeric(backsolve(ch, forwardsolve(t(ch), u[keep])))
      return(list(delta = d, rank = ncol(A), route = route, held = held))
    }
  }
  A <- pd_repair(A0)
  # nothing can be solved at a point whose curvature is not finite ANYWHERE;
  # the caller reads a zero step as "this iterate is unusable" and keeps the
  # last good one
  if (is.null(A)) {
    return(list(delta = numeric(p), rank = 0L, route = route, held = held))
  }
  if (damp > 0) Matrix::diag(A) <- Matrix::diag(A) + damp
  ch <- chol(A)
  list(delta = as.numeric(backsolve(ch, forwardsolve(t(ch), u))),
       rank = ncol(A), route = route, held = held)
}


#' The Pieces One Scoring Step Needs
#'
#' @description
#' Assembles, at the current coefficients, whichever of the square-root
#' design, the penalty's factor and the penalized information the requested
#' decomposition will actually use.
#'
#' @details
#' Only what the decomposition will read is built. `"qr"` and `"svd"` need
#' the square-root design `R` and the penalty's factor `C` and never form the
#' information; `"chol"` and `"chol_crossprod"` need the assembled `A` and
#' neither square root. Building all three unconditionally would pay for the
#' route not taken at every iteration of the loop.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, refreshed at `coef` if any term needs it.
#' @param coef A named list of coefficient vectors, one per distribution
#'   parameter.
#' @param hyper The hyperparameters, per penalized term.
#' @param method An [Iwls()] object, read for `decomposition`, `hessian` and
#'   `approx`.
#'
#' @return A list of three, of which the ones the route does not need are
#'   `NULL`:
#'   \describe{
#'     \item{`R`}{the square-root design, as [sqrt_design()] returns it, or
#'       `NULL`.}
#'     \item{`C`}{the penalty's factor, as [penalty_sqrt()] returns it, or
#'       `NULL`.}
#'     \item{`A`}{the assembled penalized information, or `NULL`.}
#'   }
#'   `R` is also `NULL` when the per-observation curvature had no Cholesky
#'   factor, and [iwls_solve()] falls back to the assembled route there.
#'
#' @seealso [iwls_solve()], which consumes this,
#'   [sqrt_design()] and [penalty_sqrt()] for the two factors.
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
  R <- sqrt_design(statmod_design_at(spec, coef, design), L)
  C <- penalty_sqrt(S)
  if (is.null(R) || is.null(C)) return(assembled())
  list(R = R, C = C, A = NULL, threads = spec@threads)
}


#' Fit the Smooth Block by Iterated Weighted Least Squares
#'
#' @description
#' Runs the scoring iteration on the objective of
#' [statmod_objective()] at fixed hyperparameters, returning the
#' coefficients and the history of the run.
#'
#' @details
#' # Two departures from a plain Newton iteration
#'
#' A curvature that is not positive definite is repaired by flooring its
#' eigenvalues, and the run continues. `solve()` would signal an error there
#' and abandon a start that is merely far from the optimum.
#'
#' The stopping rule is read at the iterate, on a score divided by the sample
#' size, so a threshold means the same at \eqn{n = 10} and at
#' \eqn{n = 10^7}. The objective itself stays unaveraged, that being the
#' scale the penalty is added on.
#'
#' # The dimensionless final verdict
#'
#' The absolute score of a location equation carries the units \eqn{1/y}. On
#' a response three decades small its rounding floor sits above any fixed
#' threshold, and a run stalled exactly at the optimum reads as a failure.
#'
#' The verdict therefore adds a second, dimensionless reading,
#' \deqn{\max_j \lvert g_j \rvert / (n\,s_p), \qquad
#'       s_p = \sqrt{\mathrm{median}_j H_{jj} / n},}
#' with \eqn{s_p} taken over the equation each coordinate belongs to.
#'
#' It can relabel a run that has already stopped as converged, and it can
#' never stop a run that is still moving. Driving the loop with it was tried
#' and made the tolerance unreachable at the other end of the scale, where
#' the objective's magnitude grows with \eqn{\log y} and the stall guard
#' fires first.
#'
#' @param obj The objective, as [statmod_objective()] returns it.
#' @param start The starting coefficients, stacked into one numeric vector.
#' @param method An [Iwls()] object, carrying the curvature, the
#'   decomposition, the budget and the stopping rule.
#' @param n The number of observations, the divisor of the scaled stopping
#'   rule.
#' @param pieces_at A function of the stacked coefficients returning what
#'   [iwls_solve()] needs, as [iwls_pieces()] builds it.
#' @param verbose `TRUE` to print one line per iteration: the objective, the
#'   score, the step length and the route taken.
#'
#' @return A list of six:
#'   \describe{
#'     \item{`par`}{the stacked coefficients reached, a numeric vector.}
#'     \item{`value`}{the penalized objective there, unaveraged.}
#'     \item{`converged`}{a single logical: whether a stopping rule fired,
#'       the dimensionless verdict included.}
#'     \item{`iterations`}{how many steps were taken. `1` means the rule was
#'       met at the first check with no step taken, as a warm
#'       start already at the mode reports.}
#'     \item{`score`}{the score per observation at the point reached.}
#'     \item{`history`}{a data frame with one row per iteration.}
#'   }
#'
#' @seealso [iwls()] for the settings, [iwls_solve()] for one step,
#'   [iwls_score()] for the dimensionless reading.
#'
#' @keywords internal
iwls_fit <- function(obj, start, method, n, pieces_at, verbose = FALSE,
                     groups = NULL) {
  beta <- start
  value <- obj$fn(beta)
  hist <- list()
  converged <- FALSE
  it <- 0L
  score <- Inf
  note <- NULL
  # NULL until a step has been taken: a rule reading a change in the
  # objective has nothing to read at the starting point, and says FALSE
  # there rather than comparing a number with itself
  f_old <- NULL
  x_old <- NULL
  # zero unless a step has failed, so the plain scoring iteration is what
  # every run that does not stall performs
  damp <- 0
  damp_tries <- 0L
  # the equations' coordinate ranges, for the dimensionless reading of the
  # final verdict; the caller says them where the coefficients it hands in
  # are a subset, since the objective's own split maps the full vector
  if (is.null(groups)) groups <- list(seq_along(beta))

  if (verbose) {
    cat(sprintf("  %-5s %14s %12s %10s\n", "iter", "objective", "score/n",
                "step"))
  }
  for (it in seq_len(as.integer(method@maxit))) {
    g <- obj$gr(beta)
    # a step whose OBJECTIVE is finite can still land where the score is not:
    # the line search below checks the value and nothing checked this, so a
    # non-finite gradient reached `if (score < tol)` and stopped the run with
    # "missing value where TRUE/FALSE needed", naming neither the iteration
    # nor the cause. The point is unusable and the last good one is kept.
    if (!all(is.finite(g))) {
      note <- paste0("the score is not finite at iteration ", it,
                     ", so the run stopped at the last usable point")
      break
    }
    score <- max(abs(g)) / n
    if (verbose) {
      cat(sprintf("  %-5d %14.6f %12.3e %10s\n", it - 1L, value, score,
                  if (it == 1L) "-" else fmt_step(step_used)))
    }
    if (iwls_met(method, list(iter = it - 1L, f_new = value, f_old = f_old,
                              x_new = beta, x_old = x_old, gradient = g / n,
                              stationarity = NULL))) {
      converged <- TRUE
      break
    }
    pc <- pieces_at(beta)
    sol <- iwls_solve(pc, -g, method@decomposition, damp)
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
      rank = sol$rank, route = sol$route, held = length(sol$held),
      damp = damp
    )
    # A step accepted for a decrease below the rounding of the objective has
    # not moved it, and the sufficient-decrease test cannot tell the two
    # apart: with a step of 1e-9 the decrease it asks for is of that order,
    # which any point satisfies. The run then spends its budget standing
    # still, which is what a term whose gradient belongs to a working model
    # rather than to the objective produces at its own fixed point.
    stalled <- value - vnew <= 1e-12 * max(1, abs(value))
    x_old <- beta
    f_old <- value
    beta <- cand
    value <- vnew
    # A STALL WITH THE SCORE STILL LARGE IS ONE COORDINATE HOLDING THE OTHERS,
    # and a shorter step cannot cure it because the step length is scalar.
    # Where a coordinate's curvature is orders below its neighbors' -- a
    # shape approaching the clamp its link keeps it inside, whose information
    # falls as nu^-4 -- the scoring step in it is astronomically long, and the
    # line search shortens the WHOLE step to keep it admissible: measured on a
    # Student t at nu = 9.1e+12, the step was -1790 in nu against at most 3.9
    # in every mean coordinate, the search ran 1, 0.125, 1.5e-05, 1.5e-08 and
    # the run stalled at a score of 2.9e-02 with the mean nowhere near
    # stationary.
    #
    # Levenberg's lambda shortens that coordinate in proportion to how little
    # curvature it has and leaves the others where they were, which is the
    # differential shrink the step length cannot give. It is raised only
    # HERE, where the loop used to give up, so a run that never stalls never
    # sees it and is unchanged; it decays again on every step that moves.
    #
    # ⚠️ ONLY WHERE THE STEP WAS SHRUNK TO NOTHING, which is the signature of
    # the deadlock and not of every stall. A term whose block is a working
    # linearization stalls AT ITS OWN FIXED POINT with the full step taken --
    # the gradient there belongs to the working model rather than to the
    # objective -- and damping past it costs iterations and ends a run that
    # had arrived: measured, escalating on every stall turned a converged
    # break-point fit into a non-converged one and moved three exact-gradient
    # checks from machine precision to 1e-3, the mode being left worse
    # located. The deadlock is the case where the line search had to shrink
    # the whole step by orders to keep one coordinate admissible.
    if (stalled) {
      if (step_used < 1e-3 && damp_tries < 8L) {
        damp <- iwls_escalate(damp, pc)
        damp_tries <- damp_tries + 1L
        next
      }
      note <- paste0("the objective stopped moving at iteration ", it)
      break
    }
    if (damp > 0) {
      damp <- damp / 100
      if (damp < iwls_scale(pc) * 1e-12) damp <- 0
    }
  }
  # The rule is asked once more at the point reached, so that a run which
  # could not move is still converged when the POINT says so. Two readings,
  # in order. The absolute score first, which is the rule that drove the
  # loop. Where it fails, the DIMENSIONLESS score arbitrates: the absolute
  # rule carries the units of the response (a location equation's score is
  # 1/y), so on a response three decades small its rounding floor sits
  # above the threshold and a run that stalled AT the optimum reads as a
  # failure -- measured, 1.37e-6 against 1e-6 at y ~ 1e-3, the objective no
  # longer moving and the fit right. The dimensionless reading, each
  # equation's score against its own information scale, survives the units.
  # It only ever RELABELS a run that has already stopped; it never stops
  # one, so every trajectory is the absolute rule's -- the lesson of
  # fit_distrib()'s dropped crit_rel_obj, and of a first version of this
  # rule that drove the loop dimensionless and made the tolerance
  # unreachable at y * 1e4, where the stall guard on the objective (whose
  # magnitude grows with log y) fires before a rule 1/s_p stricter can.
  # A caller's rule gets the same second reading and NOT the dimensionless
  # relabel, which is the built-in rule's own repair of its own units.
  if (!converged) {
    gfin <- obj$gr(beta)
    if (all(is.finite(gfin))) {
      score <- max(abs(gfin)) / n
      converged <- iwls_met(method, list(iter = it, f_new = value,
                                         f_old = f_old, x_new = beta,
                                         x_old = x_old, gradient = gfin / n,
                                         stationarity = NULL))
      if (!converged && is.null(method@criterion)) {
        hj <- iwls_info_diag(pieces_at(beta))
        dimless <- iwls_score(gfin, hj, groups, n)
        if (dimless < method@tol) {
          converged <- TRUE
          note <- paste0("converged by the dimensionless rule (score ",
                         format(dimless, digits = 3), " against ",
                         format(score, digits = 3),
                         " absolute, the response's scale being the gap)")
        }
      }
    } else {
      score <- Inf
    }
  }
  if (verbose) {
    cat(sprintf("  %-5s %14.6f %12.3e %10s\n", "end", value, score,
                if (converged) "converged" else "stopped"))
    if (!is.null(note)) cat("  ", note, "\n", sep = "")
  }
  list(par = beta, value = value, converged = converged,
       iterations = it, score = score, note = note,
       history = if (length(hist)) do.call(rbind, hist) else NULL)
}


#' The Curvature's Own Scale at One Point
#'
#' @description
#' The largest diagonal entry of the penalized information, read off the
#' pieces without assembling it.
#'
#' @details
#' The Levenberg damping is measured against this, so [iwls_escalate()]
#' carries no constant with units of its own.
#'
#' On the augmented route the diagonal of \eqn{K = R'R + C'C} is the column
#' sums of the squares of \eqn{R} and \eqn{C}, which reads a sparse design
#' without densifying it. On the assembled route it is the matrix's own
#' diagonal.
#'
#' @param pieces What [iwls_pieces()] built, carrying `R` and `C`, or `A`.
#'
#' @return A single positive number: the largest diagonal entry. `1` where
#'   the pieces carry nothing usable, which makes the damping schedule
#'   absolute there instead of relative.
#'
#' @seealso [iwls_escalate()], which divides by this,
#'   [iwls_pieces()] for the input.
#'
#' @keywords internal
iwls_scale <- function(pieces) {
  v <- if (!is.null(pieces$R)) {
    as.numeric(Matrix::colSums(pieces$R^2)) +
      if (is.null(pieces$C)) 0 else as.numeric(Matrix::colSums(pieces$C^2))
  } else if (!is.null(pieces$A)) {
    as.numeric(Matrix::diag(pieces$A))
  } else {
    return(1)
  }
  v <- v[is.finite(v)]
  if (!length(v)) return(1)
  m <- max(v)
  if (is.finite(m) && m > 0) m else 1
}


#' Raise the Levenberg Damping
#'
#' @description
#' The next \eqn{\lambda} to try after a step that failed or moved nothing.
#'
#' @details
#' The first value is \eqn{10^{-8}} of the curvature's own scale, as
#' [iwls_scale()] measures it. That is negligible beside a well-curved
#' coordinate and already large beside one whose information has fallen by
#' eight orders, which is the disparity the damping exists to correct.
#'
#' Each further escalation multiplies by a hundred, so eight of them span
#' sixteen orders and reach a damping that dominates the largest diagonal.
#' Measured on the case this was built for, the offending coordinate needed
#' \eqn{\lambda \approx 100} against a largest diagonal of 2328, which the
#' sixth escalation passes.
#'
#' @param damp The current damping, a single non-negative number. `0` starts
#'   the schedule.
#' @param pieces What [iwls_pieces()] built, read only for the scale.
#'
#' @return The next damping to try, a single positive number.
#'
#' @seealso [iwls_scale()] for the scale, [iwls_solve()] for what the damping
#'   does to the system.
#'
#' @keywords internal
iwls_escalate <- function(damp, pieces) {
  if (damp > 0) return(damp * 100)
  iwls_scale(pieces) * 1e-8
}


#' Has the Step's Stopping Rule Been Met?
#'
#' @description
#' Applies whichever stopping rule the step was configured with: the built-in
#' one, the score per observation against `tol`, or the caller's
#' \pkg{optimizers7} criterion read on the same state.
#'
#' @details
#' Both routes live here, in one place, although the loop asks the question
#' at two. That is what keeps the state a criterion is shown identical in
#' both.
#'
#' Two properties of that state matter to a caller writing a rule.
#' `state$gradient` is already the score per observation, so `crit_grad(t)`
#' and `tol = t` are the same rule. `state$objective` is the penalized
#' log-likelihood unaveraged, that being the scale the penalty is added
#' on.
#'
#' @param method An [Iwls()] object.
#' @param state The iteration state, as `optimizers7::crit_met`
#'   documents it.
#'
#' @return A single logical.
#'
#' @keywords internal
iwls_met <- function(method, state) {
  if (is.null(method@criterion)) {
    return(max(abs(state$gradient)) < method@tol)
  }
  isTRUE(optimizers7::crit_met(method@criterion, state))
}

#' The Dimensionless Reading of the Stopping Rule
#'
#' @description
#' The absolute \eqn{\max\lvert g\rvert / n} read per equation against that
#' equation's own information scale, \eqn{\max_j \lvert g_j\rvert /
#' (n\,s_p)} with \eqn{s_p = \sqrt{\mathrm{median}_j H_{jj} / n}} over the
#' equation's coordinates.
#'
#' @details
#' # Why the division is per equation
#'
#' The score of a location equation carries the units \eqn{1/y} and its
#' curvature \eqn{1/y^2}, so dividing by \eqn{s_p} survives any rescaling of
#' the response. Taking \eqn{s_p} over the whole equation changes nothing
#' within it: a stiff, heavily penalized coordinate is held to the same rule
#' as its neighbors, and that is what the envelope identities the outer
#' gradient rests on require of the mode.
#'
#' This reading arbitrates the final verdict and never drives the loop.
#'
#' # Two designs refused before this one
#'
#' A per-coordinate normalization by \eqn{(H+S)_{jj}} is dimensionless too
#' and is wrong: at extreme shrinkage the penalized coordinates converge
#' loosely, the envelope identities lose their footing, and every outer
#' trajectory moved.
#'
#' Driving the loop with the per-equation form is right at small \eqn{y} and
#' unreachable at large. At \eqn{y \times 10^4} the objective's magnitude
#' grows with \eqn{\log y}, so the stall guard on it fires before a rule
#' \eqn{1/s_p} stricter can. The inner fit then reported failure across the
#' whole corridor of smoothing parameters between the criterion's plateau and
#' its optimum, and the outer search, which reads a non-converged point as
#' unavailable, never crossed it: 1482 evaluations and a fit at cor 0.82,
#' where the absolute rule reaches 0.998.
#'
#' @param g The gradient at the point, a numeric vector over the stacked
#'   coefficients.
#' @param hj The information's diagonal at the same point, as
#'   [iwls_info_diag()] returns it, the same length as `g`.
#' @param groups A list of integer index vectors, one per equation, giving
#'   the coordinates that equation owns.
#' @param n The number of observations.
#'
#' @return A single non-negative number: the largest per-equation reading.
#'   An equation whose median curvature is zero or not finite contributes
#'   nothing.
#'
#' @seealso [fit_smooth()], whose verdict this arbitrates,
#'   [iwls()] for the absolute rule that drives the loop.
#'
#' @keywords internal
iwls_score <- function(g, hj, groups, n) {
  out <- 0
  for (idx in groups) {
    if (!length(idx)) next
    s <- sqrt(stats::median(hj[idx]) / n)
    if (!is.finite(s) || s <= 0) s <- 1
    out <- max(out, max(abs(g[idx])) / (n * s))
  }
  out
}

#' The Diagonal of the Information a Step Uses
#'
#' @description
#' Returns the diagonal of the information [iwls_score()] normalizes by, one
#' entry per stacked coefficient.
#'
#' @details
#' On the augmented route it is the column sums of the squares of the
#' square-root design, whose cross-product is the unpenalized information, so
#' the penalty does not enter.
#'
#' On the assembled route it is the diagonal of the matrix itself, which
#' folds the penalty in. That is a different normalizer, and it is accepted
#' because the assembled route is already the fallback: it is reached only
#' where no square root exists, and the reading it feeds arbitrates a verdict
#' rather than driving the loop.
#'
#' @param pieces The pieces, as [iwls_pieces()] builds them, carrying `R` or
#'   `A`.
#'
#' @return An unnamed numeric vector, one entry per stacked coefficient.
#'   [iwls_pieces()] always builds one of the two matrices, so there is no
#'   third branch.
#'
#' @seealso [iwls_score()], the only caller.
#'
#' @keywords internal
iwls_info_diag <- function(pieces) {
  if (!is.null(pieces$R)) {
    return(as.numeric(Matrix::colSums(pieces$R^2)))
  }
  as.numeric(diag(as_dense(pieces$A)))
}


#' Floor the Eigenvalues of a Curvature Matrix
#'
#' @description
#' Returns the matrix unchanged when it is positive definite, and otherwise
#' its eigendecomposition with the eigenvalues floored.
#'
#' @details
#' Abandoning a start because the curvature is indefinite is what
#' `solve()` forces and what a repair avoids; the floor is relative to
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
  # a curvature with a non-finite entry is not a matrix to repair: eigen()
  # signals "infinite or missing values in 'x'" from three frames down, which
  # names the arithmetic and not the iterate that produced it
  if (!all(is.finite(A))) return(NULL)
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
