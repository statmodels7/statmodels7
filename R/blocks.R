#' @include iwls.R
NULL

#' Split a Specification's Terms Into the Smooth Block and the Rest
#'
#' @description
#' Returns which coefficients are fitted jointly and which carry a penalty
#' with a kink and are therefore fitted by a method of their own.
#'
#' @details
#' The property that decides is one each term already reports: a penalty whose
#' \code{\link[penalties7]{penalty_kinks}} is non-empty is not twice
#' differentiable in its coefficients, so its block cannot enter a system
#' solved by a curvature. Everything else -- an unpenalized block, a ridge, a
#' spline, a random effect, a structured or additive penalty -- goes into one
#' system and is estimated all together, because their joint curvature exists
#' and using it is what makes a fit converge in a handful of iterations.
#'
#' A term whose penalty gains a smooth approximation would answer
#' \code{penalty_kinks()} differently and move into the smooth block with no
#' change here.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design, as \code{\link{statmod_design}} returns it.
#'
#' @return A list with \code{smooth} (the stacked column indices fitted
#'   jointly) and \code{sparse} (a list, one entry per non-smooth term, each
#'   with the parameter, the term's name, its stacked columns and its
#'   penalty).
#'
#' @seealso \code{\link{statmod}}
#'
#' @keywords internal
statmod_blocks <- function(spec, design) {
  npar <- vapply(design, function(d) d$npar, integer(1))
  sparse <- list()
  taken <- integer(0)
  for (u in statmod_penalized(spec, design)) {
    if (!penalty_has_kink(u$penalty,
                          sprintf("The penalty of %s", u$key))) next
    sparse[[length(sparse) + 1L]] <- list(
      param = u$param, term = u$key, cols = u$cols, index = u$index,
      penalty = u$penalty)
    taken <- c(taken, u$index)
  }
  list(smooth = setdiff(seq_len(sum(npar)), taken), sparse = sparse)
}


#' Which Coefficients Are Not Sitting at a Kink
#'
#' @description
#' A logical vector over the stacked coefficients, \code{FALSE} where one lies
#' at a point its penalty is not differentiable at.
#'
#' @details
#' The kink locations come from \code{\link[penalties7]{penalty_kinks}} read at
#' the hyperparameters in force, and a coefficient is inactive when it sits at
#' one of them. Everything outside a kinked block is active, having no kink to
#' sit at.
#'
#' The map of a kinked penalty is the identity here: a separable penalty under
#' a general map is the generalized-lasso problem, which
#' \code{\link[penalties7]{penalty_prox}} rejects, so a block that reaches this
#' function penalizes its coefficients one at a time and a kink of the penalty
#' is a kink in a coefficient.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param blocks The blocks, as \code{\link{statmod_blocks}} returns them.
#' @param beta The stacked coefficients.
#' @param hyper The hyperparameters.
#' @param tol How close to a kink counts as at it.
#'
#' @return A logical vector as long as \code{beta}.
#'
#' @seealso \code{\link{outer_tau}}
#'
#' @keywords internal
statmod_active <- function(spec, blocks, beta, hyper, tol = 1e-8) {
  active <- rep(TRUE, length(beta))
  for (b in blocks$sparse) {
    th <- as.list(hyper[[b$param]][[b$term]])
    k <- penalties7::penalty_kinks(b$penalty, th)
    k <- k[is.finite(k)]
    if (!length(k)) next
    v <- beta[b$index]
    at <- vapply(v, function(x) any(abs(x - k) <= tol), logical(1))
    active[b$index] <- !at
  }
  active
}


#' Does a Penalty Have a Kink?
#'
#' @description
#' \code{TRUE} when the penalty is not differentiable somewhere a coefficient
#' can be, which is what puts its block outside the jointly fitted system.
#'
#' @details
#' \code{\link[penalties7]{penalty_kinks}} is read at a probe value of the
#' hyperparameters -- the midpoint of their bounds, the rule
#' \pkg{modelterms7} already uses -- because whether a kink exists is a
#' property of the family and not of a point.
#'
#' A penalty that stops when asked is reported rather than treated as smooth.
#' Reading the failure as an answer sends the term to the scheme for the
#' opposite property: \code{scad()} and \code{mcp()} were fitted by the
#' curvature of a function that has none, reporting an effective 19.00
#' degrees of freedom out of 20 on a design of pure noise, which is no
#' selection at all.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param what How to name the penalty if it cannot answer.
#'
#' @return A single logical.
#'
#' @keywords internal
penalty_has_kink <- function(pen, what = "a penalty") {
  th <- as.list(penalty_theta_start(pen))
  k <- tryCatch(penalties7::penalty_kinks(pen, th), error = function(e) e)
  if (inherits(k, "error")) {
    stop(sprintf(paste0("%s cannot say whether it is differentiable.",
                        "\n  penalty_kinks() stopped with: %s",
                        "\n  The two fitting schemes differ by exactly that",
                        " property, so a penalty\n  that does not answer",
                        " cannot be assigned to one of them."),
                 what, conditionMessage(k)), call. = FALSE)
  }
  length(k) > 0L && any(is.finite(unlist(k)))
}


#' Fit One Non-Smooth Block, the Others Held Fixed
#'
#' @description
#' Runs a proximal gradient iteration on the block's own coefficients, the
#' smooth part of the objective supplying the gradient and the penalty its
#' proximal operator.
#'
#' @details
#' The route is \code{\link[optimizers7]{prox_grad}} with
#' \code{\link[penalties7]{penalty_prox}}, accelerated: measured, at a
#' condition number of 3 the plain iteration wins narrowly, at 55 it is 4153
#' iterations against 126, and at 480 the plain one does not converge in
#' 50000. A coordinate descent on the same objective is 1.1 to 5.3 times
#' faster again and is the next thing to write; it needs the columns of the
#' design and the running residual, so it belongs here and not behind an
#' optimizer's interface.
#'
#' @param obj The full objective.
#' @param beta The current stacked coefficients.
#' @param block One entry of \code{statmod_blocks()$sparse}.
#' @param hyper The hyperparameters.
#' @param maxit The iteration budget.
#' @param tol The stopping tolerance.
#' @param verbose Whether the optimizer prints its own trace.
#'
#' @return A list with \code{par} (the whole vector, updated in this block),
#'   \code{value}, \code{converged} and \code{iterations}.
#'
#' @seealso \code{\link{statmod}}
#'
#' @keywords internal
sparse_fit <- function(obj, beta, block, hyper, maxit = 500, tol = 1e-8,
                       verbose = FALSE, spec = NULL, design = NULL,
                       expected = TRUE, approx = "bartlett") {
  # a coordinate descent reads the block's own columns and the running
  # residual, which is the model rather than the objective, so it is not an
  # optimizer and lives here. It applies where the penalty can describe its
  # operator as a table; everything else takes the proximal route below.
  if (!is.null(spec) && !is.null(design)) {
    cd <- coord_fit(obj, beta, block, hyper, spec, design, expected, approx,
                    tol = tol, prev_kink = block$prev_kink)
    if (!is.null(cd)) return(cd)
  }
  idx <- block$index
  th <- as.list(hyper[[block$param]][[block$term]])
  pen <- block$penalty

  # the smooth part of the objective seen as a function of this block alone:
  # the full objective minus this penalty, which the operator applies instead
  smooth_fn <- function(b) {
    v <- beta
    v[idx] <- b
    obj$fn(v) - penalties7::penalty_value(pen, b, th)
  }
  smooth_gr <- function(b) {
    v <- beta
    v[idx] <- b
    obj$gr(v)[idx] - penalties7::penalty_gradient(pen, b, th)
  }
  prox <- function(v, step) penalties7::penalty_prox(pen, v, step, th)
  # `g` is the VALUE of the non-smooth part, not its gradient: prox_grad adds
  # it to the smooth objective to report the total it is minimizing
  gval <- function(b) penalties7::penalty_value(pen, b, th)

  opt <- optimizers7::prox_grad(prox = prox, g = gval,
                                criterion = optimizers7::crit_grad(tol),
                                maxit = maxit, verbose = verbose)
  res <- optimizers7::minimize(opt, smooth_fn, beta[idx], gr = smooth_gr)
  out <- beta
  out[idx] <- res@par
  list(par = out, value = obj$fn(out), converged = res@converged,
       iterations = res@iterations)
}


#' Every Penalized Unit of a Specification
#'
#' @description
#' One entry per penalty in the model, whatever term it belongs to and whether
#' or not that term has more than one.
#'
#' @details
#' Twelve places used to run the same loop -- over the distribution parameters,
#' over each one's terms, asking each term for its penalty -- and each of them
#' assumed a term carries at most one. A term may carry several, over different
#' subsets of its own parameters, which is what a panel model with a population
#' value and a shrunk deviation per group needs. Enumerating once is both the
#' generalization and the removal of eleven copies.
#'
#' \strong{The key} is the term's name in the formula, and the entry's own name
#' appended after \code{::} when the term carries more than one. Two
#' \code{ridge()} terms are two terms with two keys and two hyperparameters,
#' which they already were; a term with one penalty over the whole of itself
#' keys exactly as before, so nothing that reads a hyperparameter by term name
#' changes.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design, as \code{\link{statmod_design}} returns it.
#'
#' @return A list of entries, each with \code{param}, \code{term} (the name in
#'   the formula), \code{key}, \code{cols} (positions within the parameter's
#'   coefficients), \code{index} (positions in the stacked vector) and
#'   \code{penalty}.
#'
#' @seealso \code{\link{statmod_blocks}},
#'   \code{\link[modelterms7]{term_penalties}}
#'
#' @keywords internal
statmod_penalized <- function(spec, design) {
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  out <- list()
  for (a in seq_along(params)) {
    p <- params[a]
    for (nm in names(spec@terms[[p]])) {
      ent <- modelterms7::term_penalties(spec@terms[[p]][[nm]])
      if (!length(ent)) next
      block <- design[[p]]$blocks[[nm]]
      for (e in ent) {
        key <- if (length(ent) > 1L && nzchar(e$name))
          paste0(nm, "::", e$name) else nm
        cols <- block[e$index]
        out[[length(out) + 1L]] <- list(
          param = p, term = nm, key = key, cols = cols,
          index = offs[a] + cols, penalty = e$penalty)
      }
    }
  }
  out
}
