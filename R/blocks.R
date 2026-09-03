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
#' solved by a curvature. Everything else goes into one system and is
#' estimated all together: an unpenalized block, a ridge, a spline, a random
#' effect, a structured or additive penalty. Their joint curvature exists,
#' and using it closes a fit in a handful of iterations.
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
      penalty = u$penalty,
      # the sharing labels travel to the path the way the held values travel
      # to the outer index: it is the path that has to sweep one axis for a
      # group rather than one per member
      ids = u$ids)
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
#' `TRUE` when the penalty is not differentiable somewhere a coefficient can
#' be. A block whose penalty answers `TRUE` is estimated outside the jointly
#' fitted system, by a coordinate descent of its own.
#'
#' @details
#' [penalties7::penalty_kinks()] is read at a probe value of the
#' hyperparameters, the midpoint of their bounds, which is the rule
#' \pkg{modelterms7} already uses. Whether a kink exists is a property of the
#' family, not of a point, so any admissible probe answers.
#'
#' A penalty that stops when asked is reported, never treated as smooth.
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
#' faster again and is the next thing to write. It needs the columns of the
#' design and the running residual, which is the model itself, so it belongs
#' here and could not live behind an optimizer's black-box interface.
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
#' Twelve places used to run the same loop: over the distribution parameters,
#' over each one's terms, asking each term for its penalty. Every one of them
#' assumed a term carries at most one.
#'
#' A term may carry several, over different subsets of its own parameters,
#' as a panel model with a population value and a shrunk deviation
#' per group needs. Enumerating once is both the generalization and the
#' removal of eleven copies.
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
    if (!is.null(u$class)) {
      return(c(u, list(structural = FALSE, cols = NULL,
                       index = class_index(u$class, design, params, offs),
                       pieces = class_pieces(u$class, design, params, offs))))
    }
    if (S7::S7_inherits(spec@terms[[u$param]][[u$term]],
                        modelterms7::structural_term)) {
      return(c(u, list(structural = TRUE, cols = u$within, index = NULL)))
    }
    a <- match(u$param, params)
    cols <- design[[u$param]]$blocks[[u$term]][u$within]
    c(u, list(structural = FALSE, cols = cols, index = offs[a] + cols))
  })
}


#' Where a Covariance Class's Members Sit in the Stacked Vector
#'
#' @description
#' `class_pieces()` gives one entry per member with its parameter, its columns
#' in that parameter's coefficients and their positions in the stacked vector;
#' `class_index()` interleaves those positions group by group, which is the
#' order the class's penalty reads.
#'
#' @details
#' A blockwise penalty reads its argument in consecutive chunks of the prior's
#' dimension, one per group: \pkg{penalties7} reshapes the vector by row. The
#' class's prior is over the \eqn{d} columns one group carries across every
#' member, so the index must list, for each group in turn, that group's columns
#' from each member.
#'
#' Each member's own block is already group-major, so member \eqn{k}'s columns
#' for group \eqn{i} are its positions \eqn{(i-1)d_k + 1} to \eqn{id_k}. The
#' union across members is **not** contiguous -- one member's columns all
#' precede the next member's in the stacked vector, and the two may be in
#' different equations -- so what comes out is a permutation rather than a
#' range. Nothing downstream minds: reading and writing a matrix at
#' `[index, index]` is correct for any index, provided the penalty's own
#' output is in the same order, which is what this ordering arranges.
#'
#' @param cl One class, from [statmod_classes()].
#' @param design The design.
#' @param params The distribution's parameters, in order.
#' @param offs Where each parameter's coefficients start in the stacked vector.
#'
#' @return `class_pieces()` a list of lists with `param`, `term`, `cols` and
#'   `index`; `class_index()` an integer vector of `m * dim` positions.
#'
#' @seealso [statmod_penalized()], their caller.
#'
#' @keywords internal
class_pieces <- function(cl, design, params, offs) {
  lapply(cl$pieces, function(pc) {
    a <- match(pc$param, params)
    cols <- design[[pc$param]]$blocks[[pc$term]]
    # a piece written in a SUBFORMULA is part of its parent's block, and
    # `within` says which part; one written as a term of the equation is the
    # whole of it
    if (!is.null(pc$within)) cols <- cols[pc$within]
    c(pc, list(cols = cols, index = offs[a] + cols))
  })
}

#' @rdname class_pieces
#' @keywords internal
class_index <- function(cl, design, params, offs) {
  pcs <- class_pieces(cl, design, params, offs)
  out <- integer(0)
  for (i in seq_len(cl$m)) {
    for (pc in pcs) {
      d <- as.integer(pc$dim)
      out <- c(out, pc$index[(i - 1L) * d + seq_len(d)])
    }
  }
  out
}


#' The Coefficients a Penalized Unit Covers
#'
#' @description
#' The values the unit's penalty is read at, in the unit's own order.
#'
#' @details
#' Written once because a unit's coefficients are addressed differently
#' depending on what it is, and every caller wants the same vector. An ordinary
#' unit sits in one equation and its columns are positions in that parameter's
#' coefficients; a covariance class spans several and is addressed only in the
#' stacked vector. Stacking and indexing answers both, and for an ordinary unit
#' it returns exactly what `coef[[u$param]][u$cols]` returned, the stacked
#' index being that parameter's offset plus those columns.
#'
#' A structural unit has no position in the stacked vector at all: its penalty
#' covers the term's own parameters, which contribute no design column. It is
#' read from the design's structural state instead and never reaches here.
#'
#' @param u One unit, from [statmod_penalized()].
#' @param coef The coefficients, a named list by distribution parameter.
#' @param params The distribution's parameters, in order.
#'
#' @return A numeric vector as long as the unit's index.
#'
#' @seealso [statmod_penalty_at()], [statmod_marginal_grad()],
#'   [statmod_edf_correction()].
#'
#' @keywords internal
unit_beta <- function(u, coef, params) {
  unlist(coef[params], use.names = FALSE)[u$index]
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
      # a LABELLED term declares no penalty of its own: its coefficients are
      # covered by the class's, appended below
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
            as.character(e$search),
          # WHICH of its hyperparameters this entry shares with others, and
          # under what label. It travels with the entry for the reason the
          # held values do, so a penalty reached through a sub-term of a
          # structural one carries it too.
          ids = if (is.null(e$ids)) character(0) else e$ids)
      }
    }
  }
  # one unit per covariance class, over the stacked columns of its members.
  # `param` is the FIRST piece's, which is what the hyperparameter store is
  # keyed by and is a convention rather than a fact; `params` carries all of
  # them, and a reader reporting to a user reads that one.
  for (cl in statmod_classes(spec@terms)) {
    out[[length(out) + 1L]] <- list(
      param = cl$pieces[[1L]]$param, term = cl$key, key = cl$key,
      within = NULL, penalty = cl$penalty,
      fixed = list(), n_values = list(), values = list(),
      min_ratio = numeric(0), search = character(0), ids = character(0),
      class = cl,
      params = vapply(cl$pieces, function(z) z$param, ""))
  }
  out
}


#' The Covariance Classes of a Specification
#'
#' @description
#' The groups of terms that share a covariance block: those carrying the same
#' label and the same grouping variable, collected across the equations.
#'
#' @details
#' # What identifies a class
#'
#' A label and a grouping, both. The label comes from
#' [modelterms7::term_tag()] and the grouping from
#' [modelterms7::term_group()], which returns the expression, the levels and
#' the within-group column count. Two terms belong together only if the levels
#' agree as well as the expression: `droplevels(id)` and `id` are different
#' expressions for one grouping, and `id` under two subsets is one expression
#' for two groupings.
#'
#' Correlating effects on different groupings means nothing -- they are indexed
#' by different things and there is no block to estimate -- so a label used on
#' two groupings is an error naming both.
#'
#' # Whose prior it is
#'
#' The joint prior belongs to the class and not to any of its members. At most
#' one term may name a `distrib`, and its dimension must be the class's total;
#' where none does, the default is a centered multivariate Gaussian on
#' [parameters7::dr_prod()], whose coordinates are the log standard deviations
#' and the correlations' angles, so a printed hyperparameter is the quantity it
#' names. At a total of one column there is no correlation and the default is
#' the centered univariate Gaussian a single unlabelled term would have built.
#'
#' A class whose members carry priors of **different families** -- a Gaussian
#' intercept correlated with a Student t slope -- is a copula and not an
#' elliptical family; \pkg{modelterms7} rejects a univariate `distrib` on a
#' labelled term for that reason, and this function never sees the case.
#'
#' @param terms The built terms, a named list of named lists: one element
#'   per distribution parameter, in the family's order.
#'
#' @return A list, one element per class, each with `tag`, `group` (the
#'   deparsed expression), `levels`, `m`, `dim` (the class total), `penalty`
#'   and `pieces` -- one per member, with its parameter, its term's name and
#'   its within-group column count, in the order the equations were walked.
#'   Empty where no term carries a label.
#'
#' @seealso [statmod_penalty_keys()], which turns each into one penalized
#'   unit; [statmod_penalized()] for the interleaved index it is read at.
#'
#' @keywords internal
statmod_classes <- function(terms) {
  found <- list()
  for (p in names(terms)) {
    for (nm in names(terms[[p]])) {
      for (pc in label_pieces(terms[[p]][[nm]], p, nm)) {
        gr <- pc$group
        tg <- pc$tag
        key <- sprintf("%s | %s", tg, deparse1(gr$expr))
        if (is.null(found[[key]])) {
          found[[key]] <- list(tag = tg, group = deparse1(gr$expr),
                               levels = gr$levels, m = length(gr$levels),
                               pieces = list(pc))
        } else {
          if (!identical(found[[key]]$levels, gr$levels)) {
            stop(sprintf(paste0(
              "the covariance label '%s' is used on two different groupings:\n",
              "  '%s' in '%s' and '%s' in '%s' write the same expression and\n",
              "  take different levels. Effects correlated with each other are\n",
              "  indexed by one grouping."),
              tg, found[[key]]$pieces[[1L]]$term,
              found[[key]]$pieces[[1L]]$param, nm, p), call. = FALSE)
          }
          found[[key]]$pieces <- c(found[[key]]$pieces, list(pc))
        }
      }
    }
  }
  # a label written on two groupings gives two keys, which is a modelling
  # mistake and not two classes: the effects cannot be correlated
  tags <- vapply(found, function(z) z$tag, "")
  dup <- tags[duplicated(tags)]
  if (length(dup)) {
    which_ <- names(found)[tags == dup[1L]]
    stop(sprintf(paste0(
      "the covariance label '%s' is used on more than one grouping: %s.\n",
      "  Effects correlated with each other are indexed by one grouping, so\n",
      "  a label belongs to one of them."),
      dup[1L], paste(sprintf("'%s'", which_), collapse = " and ")),
      call. = FALSE)
  }
  lapply(stats::setNames(names(found), names(found)), function(key) {
    cl <- found[[key]]
    cl$key <- key
    cl$dim <- sum(vapply(cl$pieces, function(z) as.integer(z$dim), integer(1)))
    cl$penalty <- class_penalty(cl)
    cl
  })
}


#' Every Labelled Effect a Term Carries, Its Sub-Terms Included
#'
#' @description
#' One piece per labelled random effect reachable from a term: the term itself
#' where it carries a label, and the sub-terms developing its own parameters
#' otherwise, walked to any depth.
#'
#' @details
#' # Why a sub-term is reachable at all
#'
#' `seg(x, psi ~ random(~ 1 | u | id))` develops a break-point over a labelled
#' random effect. The labelled term is not one of the equation's terms -- the
#' equation carries one `SegTerm`, whose own label is absent -- and its
#' coefficients are columns of that term's block. Measured, they are exactly
#' the ones [modelterms7::term_components()] reports as that component's
#' `sub_index`, so a piece records them as `within`, positions in the parent's
#' block, and the design turns them into positions in the stacked vector the
#' same way it does for any other term.
#'
#' That is the whole of what a subformula costs here, and it is why the case is
#' covered: a labelled effect written in a subformula of an **additive** term
#' lives in the same vector as one written in an equation. A **structural**
#' parent is different -- its coefficients are its own parameters and it
#' contributes no design column -- and is rejected before reaching this, by
#' [unfittable_reason()].
#'
#' # Depth
#'
#' A sub-term is an ordinary term and may develop parameters of its own, so the
#' walk recurses and composes `within` on the way down: a depth-two sub-term's
#' columns are positions in its parent's block, which are themselves positions
#' in the equation-level term's.
#'
#' A labelled term's own sub-terms are not walked. Its columns are the class's
#' already, and anything inside it belongs to that block rather than to another.
#'
#' @param term One built term.
#' @param param The distribution parameter its equation belongs to.
#' @param nm The equation-level term's name, which is what the design is keyed
#'   by.
#' @param within The piece's columns in that term's block, or `NULL` at the top
#'   level, where the piece is the whole of it.
#'
#' @return A list of pieces, each with `param`, `term`, `within`, `dim`, `tag`,
#'   `group` (as [modelterms7::term_group()] returns it) and `distrib`. Empty
#'   where nothing under the term is labelled.
#'
#' @seealso [statmod_classes()], its caller; [class_pieces()] for the mapping
#'   of `within` onto the stacked vector.
#'
#' @keywords internal
label_pieces <- function(term, param, nm, within = NULL) {
  tg <- tryCatch(modelterms7::term_tag(term), error = function(e) NA_character_)
  if (length(tg) == 1L && !is.na(tg)) {
    gr <- modelterms7::term_group(term)
    if (is.null(gr)) {
      stop(sprintf(paste0(
        "'%s' in '%s' carries the covariance label '%s' but reports no\n",
        "  grouping variable, and effects are correlated within a grouping."),
        nm, param, tg), call. = FALSE)
    }
    return(list(list(param = param, term = nm, within = within, dim = gr$dim,
                     tag = tg, group = gr,
                     distrib = tryCatch(term@distrib, error = function(e) NULL))))
  }
  out <- list()
  comp <- tryCatch(modelterms7::term_components(term), error = function(e) list())
  for (cp in comp) {
    for (k in seq_along(cp$subs)) {
      w <- cp$sub_index[[k]]
      if (!is.null(within)) w <- within[w]
      out <- c(out, label_pieces(cp$subs[[k]], param, nm, w))
    }
  }
  out
}


#' The Joint Prior of a Covariance Class
#'
#' @description
#' The penalty over a class's stacked coefficients: the distribution one of its
#' members named, or the centered multivariate Gaussian that is the default.
#'
#' @details
#' The prior describes the effects of **one group** over every column the label
#' collects, so its dimension is the class's total and the penalty covers
#' \eqn{m} such blocks. Which chart the covariance rides is a modeling choice
#' and is [parameters7::dr_prod()] by default, where a coordinate is the
#' logarithm of a standard deviation exactly; the alternative, log-Cholesky,
#' is the same family of matrices written so that only the first coordinate
#' reads as one.
#'
#' Naming a `distrib` on more than one member is an error rather than a
#' precedence rule: the prior is one object and there is nothing to say which
#' of two should win.
#'
#' @param cl One class, as [statmod_classes()] assembles it before the penalty.
#'
#' @return A \pkg{penalties7} penalty over `m * dim` coefficients.
#'
#' @seealso [statmod_classes()], its only caller.
#'
#' @keywords internal
class_penalty <- function(cl) {
  named <- Filter(function(z) !is.null(z$distrib), cl$pieces)
  if (length(named) > 1L) {
    stop(sprintf(paste0(
      "the covariance label '%s' has a 'distrib' on %d of its terms: %s.\n",
      "  The joint prior belongs to the class and is one object, so it is\n",
      "  named on one term or on none."),
      cl$tag, length(named),
      paste(sprintf("'%s'", vapply(named, function(z) z$term, "")),
            collapse = ", ")), call. = FALSE)
  }
  d <- cl$dim
  if (length(named)) {
    pr <- named[[1L]]$distrib
    if (!identical(pr@n_dim, d)) {
      stop(sprintf(paste0(
        "the prior on '%s' is %d-variate and the label '%s' collects %d\n",
        "  columns per group. A class's prior describes the effects of one\n",
        "  group over every column the label collects."),
        named[[1L]]$term, pr@n_dim, cl$tag, d), call. = FALSE)
    }
  } else if (d == 1L) {
    # no correlation to carry at one column: the same object a single
    # unlabelled term would have built
    pr <- distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0)
  } else {
    mv <- distributions7::mvgaussian1_distrib(d, sigma = parameters7::dr_prod(d))
    pr <- do.call(distributions7::fixed,
                  c(list(mv), stats::setNames(as.list(rep(0, d)),
                                              paste0("mu", seq_len(d)))))
  }
  penalties7::distrib_penalty(pr, n_coef = cl$m * d)
}




#' Which Hyperparameters the Terms Hold
#'
#' @description
#' One key per hyperparameter a term fixed in its constructor, as
#' `parameter`, the penalty's key and its name joined by a carriage
#' return.
#'
#' @details
#' which hyperparameters are estimated is a property of the terms, not of the
#' criterion: the term is where the penalty is named, and a criterion argument
#' saying otherwise was read by nothing when the two disagreed. Everything
#' here consults this one enumeration: the outer index, the path, and the
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
#' The grid size the term asked for, or the layer's fallback where it asked
#' for nothing.
#'
#' @details
#' How finely a hyperparameter is swept is a property of the term for the
#' same reason as whether it is swept at all: a penalized block of four
#' columns and one of four hundred want different grids, and a criterion
#' applies to every hyperparameter of the model at once, the smooth ones
#' included, and those are read at the mode instead of being swept, so a
#' criterion cannot know which kind it is looking at. [modelterms7::term_grid()] is where a
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
#' The depth the term asked for, or the layer's fallback where it asked for
#' nothing.
#'
#' @details
#' One number per term, not one per hyperparameter, because only the sweep by
#' kink size reads it. A bounded hyperparameter is swept over its own
#' interval, and a shape that does not move the kink over a geometric grid
#' above its lower bound; neither takes its length from here.
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
#' It belongs to the term, for the reason the whole enumeration does: a
#' criterion is asked about every hyperparameter of the model, and a smooth
#' one is read at the mode instead of being swept, so most of what it is
#' asked about could not answer. A penalty with a kink is fitted
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
#' The grid the term named, or `NULL` where it left the path to build
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
