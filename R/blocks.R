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
#' [penalties7::penalty_kinks()] is non-empty is not twice
#' differentiable in its coefficients, so its block cannot enter a system
#' solved by a curvature. Everything else -- an unpenalized block, a ridge, a
#' spline, a random effect, a structured or additive penalty -- goes into one
#' system and is estimated all together, because their joint curvature exists
#' and using it is what makes a fit converge in a handful of iterations.
#'
#' A term whose penalty gains a smooth approximation would answer
#' `penalty_kinks()` differently and move into the smooth block with no
#' change here.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, as [statmod_design()] returns it.
#'
#' @return A list with `smooth` (the stacked column indices fitted
#'   jointly) and `sparse` (a list, one entry per non-smooth term, each
#'   with the parameter, the term's name, its stacked columns and its
#'   penalty).
#'
#' @seealso [statmod()]
#'
#' @keywords internal
statmod_blocks <- function(spec, design) {
  npar <- vapply(design, function(d) d$npar, integer(1))
  sparse <- list()
  taken <- integer(0)
  # A kinked penalty on a block that MOVES with its coefficients was refused
  # here until 0.60.0, because coord_fit() read the block as it arrived and
  # solved against the linear predictor X beta, where a term registering
  # term_refresh() has neither property: its block is the Jacobian at the
  # current coefficients and what it contributes is X beta + adj. Measured on
  # nl(~ a * exp(-r * x), a ~ 0 + lasso(~grp)) at a held lambda small enough
  # that neither a lasso nor a ridge shrinks, that solved a different model
  # and reported convergence -- a log-likelihood of -339.74 against the ridge
  # control's +155.45 and a rate of 0.22 against a truth of 0.70.
  #
  # coord_fit() now reads the block at the current coefficients and subtracts
  # `adj` from the working response, so the refusal is gone. The two agree:
  # 155.4618 against 155.4548, with the rate at 0.7360 against 0.7377.
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
#' A logical vector over the stacked coefficients, `FALSE` where one lies
#' at a point its penalty is not differentiable at.
#'
#' @details
#' The kink locations come from [penalties7::penalty_kinks()] read at
#' the hyperparameters in force, and a coefficient is inactive when it sits at
#' one of them. Everything outside a kinked block is active, having no kink to
#' sit at.
#'
#' The map of a kinked penalty is the identity here: a separable penalty under
#' a general map is the generalized-lasso problem, which
#' [penalties7::penalty_prox()] rejects, so a block that reaches this
#' function penalizes its coefficients one at a time and a kink of the penalty
#' is a kink in a coefficient.
#'
#' @param spec A [StatmodSpec()].
#' @param blocks The blocks, as [statmod_blocks()] returns them.
#' @param beta The stacked coefficients.
#' @param hyper The hyperparameters.
#' @param tol How close to a kink counts as at it.
#'
#' @return A logical vector as long as `beta`.
#'
#' @seealso [outer_tau()]
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
#' `TRUE` when the penalty is not differentiable somewhere a coefficient
#' can be, which is what puts its block outside the jointly fitted system.
#'
#' @details
#' [penalties7::penalty_kinks()] is read at a probe value of the
#' hyperparameters -- the midpoint of their bounds, the rule
#' \pkg{modelterms7} already uses -- because whether a kink exists is a
#' property of the family and not of a point.
#'
#' A penalty that stops when asked is reported rather than treated as smooth.
#' Reading the failure as an answer sends the term to the scheme for the
#' opposite property: `scad()` and `mcp()` were fitted by the
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
#' The route is [optimizers7::prox_grad()] with
#' [penalties7::penalty_prox()], accelerated: measured, at a
#' condition number of 3 the plain iteration wins narrowly, at 55 it is 4153
#' iterations against 126, and at 480 the plain one does not converge in
#' 50000. A coordinate descent on the same objective is 1.1 to 5.3 times
#' faster again and is the next thing to write; it needs the columns of the
#' design and the running residual, so it belongs here and not behind an
#' optimizer's interface.
#'
#' @param obj The full objective.
#' @param beta The current stacked coefficients.
#' @param block One entry of `statmod_blocks()$sparse`.
#' @param hyper The hyperparameters.
#' @param maxit The iteration budget.
#' @param tol The stopping tolerance.
#' @param verbose Whether the optimizer prints its own trace.
#'
#' @return A list with `par` (the whole vector, updated in this block),
#'   `value`, `converged` and `iterations`.
#'
#' @seealso [statmod()]
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
       method = "proximal gradient",
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
#' **The key** is the term's name in the formula, and the entry's own name
#' appended after `::` when the term carries more than one. Two
#' `ridge()` terms are two terms with two keys and two hyperparameters,
#' which they already were; a term with one penalty over the whole of itself
#' keys exactly as before, so nothing that reads a hyperparameter by term name
#' changes.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, as [statmod_design()] returns it.
#'
#' @return A list of entries, each with `param`, `term` (the name in
#'   the formula), `key`, `cols` (positions within the parameter's
#'   coefficients), `index` (positions in the stacked vector) and
#'   `penalty`.
#'
#' @seealso [statmod_blocks()],
#'   [modelterms7::term_penalties()]
#'
#' @keywords internal
statmod_penalized <- function(spec, design) {
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  lapply(statmod_penalty_keys(spec), function(u) {
    # A structural term contributes no design columns, so its penalty covers
    # positions among the TERM'S OWN parameters and there is nothing to look
    # up in the design. Reading `blocks[[term]]` for it returned NULL, and the
    # positions then indexed the equation's coefficients -- 25 of them where
    # the equation has one -- so the penalty was evaluated at NA and every
    # quantity built on it was not finite.
    if (S7::S7_inherits(spec@terms[[u$param]][[u$term]],
                        modelterms7::structural_term)) {
      return(c(u, list(structural = TRUE, cols = u$within, index = NULL)))
    }
    a <- match(u$param, params)
    cols <- design[[u$param]]$blocks[[u$term]][u$within]
    c(u, list(structural = FALSE, cols = cols, index = offs[a] + cols))
  })
}


#' Every Penalty in a Model, Without the Design
#'
#' @description
#' The same enumeration as [statmod_penalized()] minus the column
#' positions, for the callers that need to know which penalties exist before a
#' design has been built.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A list of entries with `param`, `term`, `key`,
#'   `within` (positions among the term's own parameters) and
#'   `penalty`.
#'
#' @seealso [statmod_penalized()]
#'
#' @keywords internal
statmod_penalty_keys <- function(spec) {
  out <- list()
  for (p in spec@distrib@params) {
    for (nm in names(spec@terms[[p]])) {
      ent <- modelterms7::term_penalties(spec@terms[[p]][[nm]])
      if (!length(ent)) next
      for (e in ent) {
        out[[length(out) + 1L]] <- list(
          param = p, term = nm, key = statmod_entry_key(nm, ent, e),
          within = e$index, penalty = e$penalty,
          # what the TERM holds and how finely it wants each swept: the
          # entry carries both, so a penalty reached through a sub-term of
          # a structural one carries them too
          fixed = if (is.null(e$fixed)) list() else e$fixed,
          n_values = if (is.null(e$n_values)) list() else e$n_values,
          values = if (is.null(e$values)) list() else e$values,
          min_ratio = if (is.null(e$min_ratio)) numeric(0) else
            as.numeric(e$min_ratio),
          search = if (is.null(e$search)) character(0) else
            as.character(e$search))
      }
    }
  }
  out
}


#' Which Hyperparameters the Terms Hold
#'
#' @description
#' One key per hyperparameter a term fixed in its constructor, as
#' `parameter`, the penalty's key and its name joined by a carriage
#' return.
#'
#' @details
#' WHICH hyperparameters are estimated is a property of the terms, not of the
#' criterion: the term is where the penalty is named, and a criterion argument
#' saying otherwise was read by nothing when the two disagreed. Everything
#' here consults this one enumeration -- the outer index, the path, and the
#' summary's account of what was estimated and what was given.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#'
#' @return A character vector, possibly empty.
#'
#' @seealso [modelterms7::term_hyper()], [outer_hyper_index()]
#'
#' @keywords internal
statmod_held <- function(spec, design = NULL) {
  if (is.null(design)) design <- statmod_design(spec)
  out <- character(0)
  for (u in statmod_penalized(spec, design)) {
    for (h in names(u$fixed)) {
      out <- c(out, paste(u$param, u$key, h, sep = "\r"))
    }
  }
  out
}


#' The Key of One of a Term's Penalties
#'
#' @description
#' What a hyperparameter row is filed under: the term's name where the term
#' carries one penalty, and `term::entry` where it carries several.
#'
#' @details
#' The composition is written once because two callers reading a
#' hyperparameter by a key they each compose would agree only by accident. A
#' term carrying one penalty keys as it always did, so a formula of ordinary
#' terms is unaffected by the entry names the terms now supply.
#'
#' @param term The term's name in the formula.
#' @param entries The term's entries, as
#'   [modelterms7::term_penalties()] returns them.
#' @param entry One of them.
#'
#' @return A single string.
#'
#' @keywords internal
statmod_entry_key <- function(term, entries, entry) {
  if (length(entries) > 1L && nzchar(entry$name)) {
    paste0(term, "::", entry$name)
  } else {
    term
  }
}


#' One Penalized Unit, by Parameter and Key
#'
#' @description
#' The entry of [statmod_penalized()] a hyperparameter row names, or
#' `NULL` where there is none.
#'
#' @details
#' The places that read a penalty from a `(parameter, term)` pair used to
#' fetch it with `term_penalty()` and take the term's whole block, which
#' assumes one penalty per term. Looking it up here answers the same question
#' where that holds and the right question where it does not.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param param The distribution parameter.
#' @param key The key, as `statmod_penalized()` composes it.
#'
#' @return One entry, or `NULL`.
#'
#' @keywords internal
statmod_unit <- function(spec, design, param, key) {
  for (u in statmod_penalized(spec, design)) {
    if (identical(u$param, param) && identical(u$key, key)) return(u)
  }
  NULL
}


#' How Many Values a Path Visits for One Hyperparameter
#'
#' @description
#' The grid size the TERM asked for, or the layer's fallback where it asked
#' for nothing.
#'
#' @details
#' How finely a hyperparameter is swept is a property of the term for the
#' same reason as whether it is swept at all: a penalized block of four
#' columns and one of four hundred want different grids, and a criterion
#' applies to every hyperparameter of the model at once -- the smooth ones
#' included, which are read at the mode and not swept -- so it cannot know
#' which it is looking at. [modelterms7::term_grid()] is where a
#' term says so, and the value travels with the penalty's entry, so one
#' reached through a sub-term of a structural term carries it too.
#'
#' @param spec A [StatmodSpec()].
#' @param row One row of [path_rows()]'s index.
#' @param default [path_fallbacks()]'s, for a term that named none.
#'
#' @return A single integer.
#'
#' @seealso [statmod_held()], [path_values()]
#'
#' @keywords internal
statmod_grid_size <- function(spec, row, default) {
  statmod_path_setting(spec, row, "n_values", default, row$name)
}


#' How Far Down the Path Reaches for One Term
#'
#' @description
#' The depth the TERM asked for, or the layer's fallback where it asked for
#' nothing.
#'
#' @details
#' One number per term rather than one per hyperparameter, because only the
#' sweep by kink size uses it: a bounded hyperparameter is swept over its own
#' interval and a shape that does not move the kink over a geometric grid
#' above its lower bound.
#'
#' @param spec A [StatmodSpec()].
#' @param row One row of [path_rows()]'s index.
#' @param default [path_fallbacks()]'s, for a term that named none.
#'
#' @return A single number.
#'
#' @seealso [statmod_grid_size()], [path_values()]
#'
#' @keywords internal
statmod_min_ratio <- function(spec, row, default) {
  statmod_path_setting(spec, row, "min_ratio", default, NULL)
}


#' How a Term Covers Its Own Hyperparameters
#'
#' @description
#' `"grid"` for every combination of the term's kinked hyperparameters,
#' `"cyclic"` for one at a time, or the default where the term named
#' neither.
#'
#' @details
#' It belongs to the TERM and not to the criterion, for the reason the whole
#' enumeration does: a criterion is asked of every hyperparameter of the
#' model, and a smooth one is read at the mode rather than swept, so most of
#' what it is asked about could not answer. A penalty with a kink is fitted
#' by a scheme of its own, and how that scheme covers the term's own
#' hyperparameters is part of the scheme.
#'
#' Being per term is also what keeps one term's choice off another's:
#' `y ~ lasso(X) + enet(R, search = "cyclic")` sweeps the elastic net
#' one coordinate at a time and leaves the lasso alone.
#'
#' @param spec A [StatmodSpec()].
#' @param row One row of [path_rows()]'s index.
#' @param default What to use where the term named nothing.
#'
#' @return `"grid"` or `"cyclic"`.
#'
#' @seealso [modelterms7::term_search()], [statmod_path()]
#'
#' @keywords internal
statmod_search <- function(spec, row, default = "grid") {
  for (u in statmod_penalty_keys(spec)) {
    if (!identical(u$param, row$parameter) ||
        !identical(u$key, row$term)) next
    v <- u$search
    if (!is.null(v) && length(v) && nzchar(v)) return(as.character(v)[[1L]])
  }
  default
}


#' The Values a Term Wrote Out for One Hyperparameter
#'
#' @description
#' The grid the TERM named, or `NULL` where it left the path to build
#' one.
#'
#' @details
#' The third state of a hyperparameter's argument, beside holding it at one
#' number and leaving it to be estimated over a grid the path constructs. It
#' travels with the penalty's entry like the grid size and the depth, so a
#' penalty reached through a sub-term of a structural term carries it too.
#'
#' @param spec A [StatmodSpec()].
#' @param row One row of [path_rows()]'s index.
#'
#' @return A numeric vector, or `NULL`.
#'
#' @seealso [modelterms7::term_values()], [path_forced()]
#'
#' @keywords internal
statmod_values <- function(spec, row) {
  for (u in statmod_penalty_keys(spec)) {
    if (!identical(u$param, row$parameter) ||
        !identical(u$key, row$term)) next
    v <- u$values[[row$name]]
    if (!is.null(v) && length(v)) return(as.numeric(v))
  }
  NULL
}


#' One Setting of the Path, Read From the Term
#'
#' @param spec A [StatmodSpec()].
#' @param row One row of [path_rows()]'s index.
#' @param field Which field of the penalty's entry to read.
#' @param default What the criterion asks for.
#' @param name The hyperparameter's name, or `NULL` where the setting
#'   is one per term.
#'
#' @return A single number.
#'
#' @keywords internal
statmod_path_setting <- function(spec, row, field, default, name = NULL) {
  for (u in statmod_penalty_keys(spec)) {
    if (!identical(u$param, row$parameter) ||
        !identical(u$key, row$term)) next
    v <- u[[field]]
    n <- if (is.null(name)) v else v[[name]]
    if (!is.null(n) && length(n) && is.finite(n)) return(n)
  }
  default
}
