#' @include iwls.R
NULL

#' Split a Specification's Terms Into the Smooth Block and the Rest
#'
#' @description
#' Returns which coefficients are fitted jointly and which carry a penalty
#' with a kink and are therefore fitted by a method of their own.
#'
#' @details
#' The split depends on the penalty of each penalized unit (see
#' [statmod_penalized()]). A penalty for which [penalties7::penalty_kinks()]
#' returns a finite point is not twice differentiable in its coefficients, so
#' its block cannot enter a system solved by a curvature. Everything else goes
#' into one system and is estimated all together: an unpenalized block, a
#' ridge, a spline, a random effect, a structured or additive penalty. Their
#' joint curvature exists, so a scoring or Newton step applies to them.
#'
#' The classification is made by [penalty_has_kink()], so a penalty that
#' becomes smooth moves into the smooth block with no change here.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, as [statmod_design()] returns it.
#'
#' @return A list with `smooth` (the stacked column indices fitted
#'   jointly) and `sparse` (a list with one entry per penalized unit that has
#'   a kink). Each entry holds `param`, `term` (the penalty's key), `cols`
#'   (positions within the parameter's coefficients), `index` (positions in the
#'   stacked vector), `penalty`, `structural`, `zterm` (the term's own name)
#'   and `ids` (the labels under which hyperparameters are shared).
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
      # A STRUCTURAL UNIT IS ADDRESSED IN ANOTHER VECTOR, and both routes
      # below need to know it: `cols` are positions among the TERM's own
      # parameters rather than columns of any design, and `index` is empty
      # because such a penalty covers no coefficient. `zterm` is the term's
      # own name, which the design's structural state is keyed by, where
      # `term` is the penalty's key and carries `::entry` where the term
      # declares more than one.
      structural = isTRUE(u$structural), zterm = u$term,
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
#' Returns a logical vector over the stacked coefficients, with `FALSE` where a
#' coefficient lies at a point at which its penalty is not differentiable.
#'
#' @details
#' The kink locations come from [penalties7::penalty_kinks()] read at
#' the hyperparameters in force, and a coefficient is inactive when it sits at
#' one of them. Everything outside a kinked block is active, because it has no
#' kink at which to sit.
#'
#' A kinked penalty reaches this function with no map or with a diagonal map;
#' `standardize = TRUE` writes a diagonal map. A separable penalty under a
#' general map is the generalized-lasso problem, which
#' [penalties7::penalty_prox()] rejects. A block that reaches this function
#' therefore penalizes its coefficients one at a time, and a kink of the
#' penalty is a kink in a coefficient.
#'
#' @param spec A [StatmodSpec()].
#' @param blocks The blocks, as [statmod_blocks()] returns them.
#' @param beta The stacked coefficients.
#' @param hyper The hyperparameters.
#' @param tol The distance to a kink within which a coefficient is taken to be
#'   at the kink.
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
#' Returns `TRUE` when the penalty is not differentiable at some point that a
#' coefficient can take. A block whose penalty gives `TRUE` is estimated outside
#' the jointly fitted system by [sparse_fit()]: a compiled coordinate descent
#' where the penalty describes its proximal operator as a table, and an
#' accelerated proximal gradient iteration otherwise.
#'
#' @details
#' [penalties7::penalty_kinks()] is read at a probe value of the
#' hyperparameters, which [penalty_theta_start()] gives. The probe is the
#' midpoint of the interval where both bounds are finite, one unit inside the
#' bound where only one is finite, and zero where both bounds are infinite.
#' \pkg{modelterms7} uses the same rule. Whether a kink exists is a property
#' of the family and not of a point, so any admissible probe gives the same
#' answer.
#'
#' If [penalties7::penalty_kinks()] signals an error, `penalty_has_kink()`
#' signals an error that names the penalty, and the penalty is not treated as
#' smooth. If the failure were read as an answer, a kinked penalty, such as
#' `scad()` or `mcp()`, would be sent to the scheme for smooth penalties, which
#' fits it by the curvature of a function that has none.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param what How to name the penalty in the error message.
#'
#' @return A single logical.
#'
#' @keywords internal
penalty_has_kink <- function(pen, what = "a penalty") {
  th <- as.list(penalty_theta_start(pen))
  k <- tryCatch(penalties7::penalty_kinks(pen, th), error = function(e) e)
  if (inherits(k, "error")) {
    stop(sprintf(paste0("It cannot be determined whether %s is differentiable.",
                        "\n  penalty_kinks() failed with: %s",
                        "\n  The two fitting schemes differ by exactly that",
                        " property, so a penalty\n  whose kinks cannot be",
                        " evaluated cannot be assigned to either of them."),
                 what, conditionMessage(k)), call. = FALSE)
  }
  length(k) > 0L && any(is.finite(unlist(k)))
}


#' Fit One Non-Smooth Block While the Others Are Held Fixed
#'
#' @description
#' Estimates the coefficients of one block with a kinked penalty while the other
#' coefficients are held fixed. A penalty over a structural term's own
#' parameters is fitted by [sparse_fit_structural()]. Otherwise a compiled
#' coordinate descent, [coord_fit()], is used where the penalty describes its
#' proximal operator as a table, and an accelerated proximal gradient iteration
#' runs on the block where it does not, in which the smooth part of the
#' objective supplies the gradient and the penalty supplies its proximal
#' operator.
#'
#' @details
#' The coordinate descent of [coord_fit()] reads the columns of the design and
#' the running residual, which belong to the model, so it is called here and
#' not through an optimizer's black-box interface. Where [coord_fit()] returns
#' `NULL` (the penalty has no proximal table at the steps of this block, or the
#' working weights are unusable), or where `spec` or `design` is `NULL`, the
#' route is [optimizers7::prox_grad()] with [penalties7::penalty_prox()],
#' accelerated.
#'
#' @param obj The full objective.
#' @param beta The current stacked coefficients.
#' @param block One entry of `statmod_blocks()$sparse`.
#' @param hyper The hyperparameters.
#' @param maxit The iteration budget of the proximal gradient route and of the
#'   structural route. The coordinate descent route has a budget of its own,
#'   set in [coord_fit()] and [coord_call()].
#' @param tol The stopping tolerance.
#' @param verbose Whether the optimizer of the proximal routes prints its own
#'   trace. The coordinate descent route does not use it.
#' @param spec A [StatmodSpec()]. If `spec` or `design` is `NULL`, the
#'   coordinate descent is skipped.
#' @param design The design, as [statmod_design()] returns it.
#' @param expected `TRUE` for the expected information in the working weights
#'   of the coordinate descent, `FALSE` for the observed one.
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#'
#' @return A list with `par` (the whole vector, updated in this block),
#'   `value`, `converged`, `iterations` and `method` (`"coordinate descent"`,
#'   `"proximal gradient"` or `"proximal gradient (structural)"`).
#'
#' @seealso [statmod()]
#'
#' @keywords internal
sparse_fit <- function(obj, beta, block, hyper, maxit = 500, tol = 1e-8,
                       verbose = FALSE, spec = NULL, design = NULL,
                       expected = TRUE, emethod = NULL) {
  # A PENALTY OVER A STRUCTURAL TERM'S OWN PARAMETERS COVERS NO COLUMN AND NO
  # ENTRY OF `beta`, so neither route below can address it: the coordinate
  # descent reads `cols` as columns of a design, which those positions are
  # not, and the proximal branch reads `beta[block$index]`, which is empty.
  # Measured, a lasso held there died on "subscript out of bounds" inside
  # coord_block(), in the published 0.124.0 as well.
  if (isTRUE(block$structural)) {
    return(sparse_fit_structural(obj, beta, block, hyper, spec, design,
                                 maxit = maxit, tol = tol, verbose = verbose))
  }
  # a coordinate descent reads the block's own columns and the running
  # residual, which is the model rather than the objective, so it is not an
  # optimizer and lives here. It applies where the penalty can describe its
  # operator as a table; everything else takes the proximal route below.
  if (!is.null(spec) && !is.null(design)) {
    cd <- coord_fit(obj, beta, block, hyper, spec, design, expected, emethod,
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
#' Returns one entry per penalty in the model, regardless of the term to which
#' the penalty belongs and of whether that term has more than one penalty.
#'
#' @details
#' The enumeration runs over the distribution parameters and over the terms of
#' each, and reads the penalties of each term. A term may carry several
#' penalties, over different subsets of its own parameters, as a panel model
#' with a population value and a shrunk deviation per group needs. A
#' covariance class (see [statmod_classes()]) is one unit over the stacked
#' columns of its members. The callers that need the penalties of a model
#' read this one enumeration.
#'
#' The key is the term's name in the formula. When the term carries more than
#' one penalty and the entry is named, the entry's own name is appended after
#' `::`. Two `ridge()` terms are two terms with two keys and two
#' hyperparameters. A term with one penalty has only its name as key, so a
#' hyperparameter is read by the term's name.
#'
#' `statmod_penalized_units()` builds the list. `statmod_penalized()`
#' returns the list from the design's `eta_memo` environment when the design
#' has one and the terms and the family's parameters are the ones for which the
#' design was built. An objective evaluated many times over one design reads
#' the list from there.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, as [statmod_design()] returns it.
#'
#' @return A list of entries, each with `param`, `term` (the name in
#'   the formula), `key`, `within`, `penalty`, `fixed`, `n_values`, `values`,
#'   `min_ratio`, `search`, `ids`, `structural`, `cols` and `index`. For an
#'   ordinary unit `cols` are positions within the parameter's coefficients and
#'   `index` positions in the stacked vector. For a structural unit `cols` are
#'   positions among the term's own parameters and `index` is `NULL`. The entry
#'   of a covariance class also carries `class`, `params`, `pieces`, `mixed`
#'   and, for a mixed class, `joint` and `beta_index`.
#'
#' @seealso [statmod_blocks()],
#'   [modelterms7::term_penalties()]
#'
#' @keywords internal
statmod_penalized <- function(spec, design) {
  params <- spec@distrib@params
  # a function of the terms, the family's parameters and the design alone,
  # and asked at every evaluation of the objective: where the design keeps a
  # memo the enumeration is kept there, keyed on what it reads
  mm <- attr(design, "eta_memo")
  if (!is.null(mm) && !is.null(mm$pen_value) &&
      identical(mm$pen_terms, spec@terms) &&
      identical(mm$pen_params, params)) {
    return(mm$pen_value)
  }
  out <- statmod_penalized_units(spec, design)
  if (!is.null(mm)) {
    mm$pen_terms <- spec@terms
    mm$pen_params <- params
    mm$pen_value <- out
  }
  out
}

#' @rdname statmod_penalized
#' @keywords internal
statmod_penalized_units <- function(spec, design) {
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
      cl <- u$class
      # a class inside a structural term is addressed exactly as that term's
      # own penalties are: positions among its parameters, read from the
      # design's structural state, with nothing to look up in the design and
      # no place in the stacked vector
      if (identical(cl$space, "zeta")) {
        return(c(u, list(structural = TRUE, mixed = FALSE,
                         cols = class_zeta_cols(cl), index = NULL,
                         pieces = lapply(cl$pieces, function(pc)
                           c(pc, list(cols = pc$within, index = NULL))))))
      }
      # a MIXED class is in neither vector on its own. It has no `index` and
      # no `cols`, so every reader written for one of the two skips it, and
      # what it does carry is `joint`: positions in [beta ; zeta_free], which
      # is the vector the inner step, the criterion and the variance build.
      if (identical(cl$space, "mixed")) {
        jp <- class_joint_pieces(cl, design, params, offs)
        return(c(u, list(structural = FALSE, mixed = TRUE, cols = NULL,
                         index = NULL,
                         joint = class_interleave(
                           cl, lapply(jp, function(z) z$joint)),
                         beta_index = unlist(lapply(jp, function(z) z$index),
                                             use.names = FALSE),
                         pieces = jp)))
      }
      return(c(u, list(structural = FALSE, mixed = FALSE, cols = NULL,
                       index = class_index(cl, design, params, offs),
                       pieces = class_pieces(cl, design, params, offs))))
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
#' `class_index()` interleaves those positions group by group, in the order in
#' which the class's penalty reads them.
#'
#' @details
#' A blockwise penalty reads its argument in consecutive chunks of the prior's
#' dimension, one per group: \pkg{penalties7} reshapes the vector by row. The
#' class's prior is over the \eqn{d} columns that one group carries across every
#' member, so the index must list, for each group in turn, that group's columns
#' from each member.
#'
#' Each member's own block is already group-major, so member \eqn{k}'s columns
#' for group \eqn{i} are its positions \eqn{(i-1)d_k + 1} to \eqn{id_k}. The
#' union across members is not contiguous (one member's columns all
#' precede the next member's in the stacked vector, and the two may be in
#' different equations), so the index is a permutation and not a range.
#' Reading and writing a matrix at `[index, index]` is correct for any index,
#' provided the penalty's own output is in the same order, and this ordering
#' arranges that.
#'
#' A class whose members are all inside a structural term is addressed in that
#' term's own parameters instead. That vector is the one at which its penalty
#' is read, and the design has no column in it. `class_zeta_cols()`
#' interleaves those positions by the same rule, and `class_pieces()` is not
#' consulted at all, because there is no design block in which to look up the
#' columns. A class split
#' between coefficients and a structural term's parameters is addressed in the
#' joint vector by [class_joint_pieces()].
#'
#' @param cl One class, from [statmod_classes()].
#' @param design The design.
#' @param params The distribution's parameters, in order.
#' @param offs Where each parameter's coefficients start in the stacked vector.
#'
#' @return `class_pieces()` a list with one element per member: the member's
#'   piece (see [label_pieces()]), to which `cols` and `index` are added;
#'   `class_index()` and `class_zeta_cols()` an integer vector of `m * dim`
#'   positions.
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
#' @param pos One integer vector per member, each group-major, giving that
#'   member's positions in the vector in which the class is addressed.
#' @keywords internal
class_interleave <- function(cl, pos) {
  out <- integer(0)
  for (i in seq_len(cl$m)) {
    for (k in seq_along(pos)) {
      d <- as.integer(cl$pieces[[k]]$dim)
      out <- c(out, pos[[k]][(i - 1L) * d + seq_len(d)])
    }
  }
  out
}

#' @rdname class_pieces
#' @keywords internal
class_index <- function(cl, design, params, offs) {
  class_interleave(cl, lapply(class_pieces(cl, design, params, offs),
                              function(pc) pc$index))
}

#' @rdname class_pieces
#' @keywords internal
class_zeta_cols <- function(cl) {
  class_interleave(cl, lapply(cl$pieces, function(pc) pc$within))
}


#' Where a Mixed Covariance Class's Members Sit in the Joint Vector
#'
#' @description
#' Returns one entry per member of a class split between the design's
#' coefficients and a filter's own parameters, each with its positions in the
#' joint vector \eqn{[\beta; \zeta_{\mathrm{free}}]} that the inner step, the
#' marginal criterion and the variance all assemble.
#'
#' @details
#' The joint order is the one that those three already write and is not
#' invented here: the stacked coefficients first, then the free parameters of
#' the one structural term of the filter shape, in the order in which the term
#' holds them, less whichever one a linear intercept already carries.
#'
#' A coefficient member's positions are its parameter's offset plus its
#' columns, exactly as [class_pieces()] gives them. A structural member's
#' positions are `nb` plus its place among the free parameters, which is not its
#' place among all of them, because a held parameter is not in the joint vector
#' at all.
#'
#' # A held coordinate cannot be in a class
#'
#' The prior's dimension is fixed when the class is assembled, from the members'
#' widths. If one of the coordinates that it collects is then held (the level
#' of a filter is, wherever a linear intercept in the same equation already
#' carries the constant), the class has fewer coordinates than its prior
#' describes and there is no covariance matrix of the prior's dimension to
#' estimate. It is rejected here, where the held set is visible, naming the
#' coordinate and the two ways out.
#'
#' @param cl One class, from [statmod_classes()], whose `space` is `"mixed"`.
#' @param design The design.
#' @param params The distribution's parameters, in order.
#' @param offs Where each parameter's coefficients start in the stacked vector.
#'
#' @return A list of lists, one per member; each is the piece to which `joint`
#'   is added.
#'   A coefficient member also carries `cols` and `index`. A structural member
#'   carries `zcols` (positions among the term's own parameters), and its
#'   `cols` and `index` are `NULL`.
#'
#' @seealso [statmod_penalized()], its caller; [joint_penalty_at()], which
#'   reads the positions.
#'
#' @keywords internal
class_joint_pieces <- function(cl, design, params, offs) {
  sst <- statmod_structural_state(design)
  nb <- sum(vapply(design[params], function(d) d$npar, integer(1)))
  key <- cl$sterm
  all_nm <- names(sst$zeta[[key]])
  free <- setdiff(all_nm, sst$held[[key]])
  lapply(cl$pieces, function(pc) {
    if (isTRUE(pc$structural)) {
      nm <- all_nm[pc$within]
      j <- match(nm, free)
      if (anyNA(j)) {
        stop(sprintf(paste0(
          "the covariance label '%s' collects '%s', which is held: a linear\n",
          "  intercept in the same equation already carries it, so it is not\n",
          "  estimated and the block would have one coordinate fewer than the\n",
          "  prior describes. Write the equation with `0 +` so the term\n",
          "  carries its own level, or drop the middle bar on that effect."),
          cl$tag, nm[is.na(j)][[1L]]), call. = FALSE)
      }
      # `cols` stays NULL so that a reader marking DESIGN columns from the
      # pieces skips this one rather than marking the columns its parameter
      # positions happen to number
      return(c(pc, list(cols = NULL, index = NULL, zcols = pc$within,
                        joint = nb + j)))
    }
    a <- match(pc$param, params)
    cols <- design[[pc$param]]$blocks[[pc$term]]
    if (!is.null(pc$within)) cols <- cols[pc$within]
    c(pc, list(cols = cols, index = offs[a] + cols, joint = offs[a] + cols))
  })
}


#' The Coefficients That a Penalized Unit Covers
#'
#' @description
#' Returns the values at which the unit's penalty is read, in the unit's own
#' order.
#'
#' @details
#' This function is written once because a unit's coefficients are addressed
#' differently depending on what the unit is, and every caller needs the same
#' vector. An ordinary unit sits in one equation and its columns are positions
#' in that parameter's coefficients; a covariance class may span several
#' equations and is addressed only in the stacked vector. Stacking and indexing
#' serves both, and for an ordinary unit it returns exactly what
#' `coef[[u$param]][u$cols]` returned, because the stacked index is that
#' parameter's offset plus those columns.
#'
#' A structural unit has no position in the stacked vector at all: its penalty
#' covers the term's own parameters, which contribute no design column. It is
#' read from the design's structural state instead. A mixed covariance class
#' is addressed in the joint vector. [unit_joint_beta()] reads both of them.
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


#' Where a Penalized Unit's Coordinates Are
#'
#' @description
#' Returns the positions that the unit's penalty covers, in the vector in which
#' the unit is addressed: the stacked coefficients for an ordinary unit, the
#' structural term's own parameters for one whose coefficients are a filter's,
#' and the joint vector of coefficients and free structural parameters for a
#' mixed class.
#'
#' @details
#' A reader that only needs to know how many coordinates a penalty covers, or
#' to which of a term's components they belong, needs the same answer for every
#' kind of unit, and the kinds keep their positions in different fields. An
#' ordinary unit and a covariance class among the coefficients carry them in
#' `index`. A structural unit has no `index` and carries them in `cols`, among
#' the term's own parameters. A mixed class has neither and carries them in
#' `joint`.
#'
#' @param u One unit, from [statmod_penalized()], or a class carrying the same
#'   fields.
#'
#' @return An integer vector.
#'
#' @seealso [unit_beta()] for the values at those positions.
#'
#' @keywords internal
unit_positions <- function(u) {
  # a mixed class spans coefficients and a filter's parameters, and only the
  # joint vector addresses both
  if (isTRUE(u$mixed)) return(as.integer(u$joint))
  if (isTRUE(u$structural)) u$cols else u$index
}


#' Where a Penalized Unit's Coordinates Are in the Joint Vector
#'
#' @description
#' Returns the same positions as [unit_positions()], carried onto the one vector
#' \eqn{[\beta; \zeta_{\mathrm{free}}]} that the inner step, the marginal
#' criterion and the variance all assemble, so that a unit of any of the three
#' kinds can be read off one diagonal.
#'
#' @details
#' The three kinds are addressed in three different vectors, and a caller that
#' forms a trace needs one. A mixed class already records the answer, because
#' [class_joint_pieces()] builds its `joint` field. A structural unit
#' records positions among all of a term's parameters, where the joint vector
#' holds only the free ones, so a held coordinate has to be mapped out rather
#' than offset past. An ordinary unit's stacked index is already a position in
#' the joint vector, because the coefficients come first and every one of them
#' is a free coordinate.
#'
#' This function is written once so that the convention is composed in one
#' place and not in three.
#'
#' @param u One unit, from [statmod_penalized()].
#' @param spec A [StatmodSpec()].
#' @param design The design.
#'
#' @return An integer vector of positions, empty where the unit names a
#'   coordinate that the joint vector does not carry.
#'
#' @seealso [unit_positions()] for the same coordinates in the unit's own
#'   vector; [class_joint_pieces()], which builds the mixed case.
#'
#' @keywords internal
unit_joint_positions <- function(u, spec, design) {
  if (isTRUE(u$mixed)) return(as.integer(u$joint))
  if (!isTRUE(u$structural)) return(as.integer(u$index))
  sst <- statmod_structural_state(design)
  if (is.null(sst) || is.null(sst$zeta[[u$term]])) return(integer(0))
  nb <- sum(vapply(design[spec@distrib@params], function(d) d$npar,
                   integer(1)))
  all_nm <- names(sst$zeta[[u$term]])
  j <- match(all_nm[u$cols], setdiff(all_nm, sst$held[[u$term]]))
  if (anyNA(j)) integer(0) else as.integer(nb + j)
}


#' A Unit's Own Coordinates, in the Joint Vector's Order
#'
#' @description
#' Returns the values whose positions [unit_joint_positions()] gives, in the
#' order in which that function returns those positions.
#'
#' @details
#' The three kinds are read from two different places. An ordinary unit's
#' coordinates are coefficients and come from `coef`. A structural unit's
#' coordinates are a filter's own parameters and come from the design's
#' structural state, on the unconstrained scale, which is the scale on which the
#' recursion is differentiated and which the penalty covers. A mixed class holds
#' some of each and interleaves them group by group, so its pieces are read
#' separately and put back in the class's own order, which
#' [class_joint_pieces()] records.
#'
#' This function is written beside [unit_joint_positions()] because a value and
#' its position must follow the same convention. A caller that composed them
#' separately could place one unit's values at another unit's positions.
#'
#' @param u One unit, from [statmod_penalized()].
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients, one entry per distribution parameter.
#'
#' @return A numeric vector as long as [unit_joint_positions()]'s answer.
#'
#' @seealso [unit_beta()], the coefficient-only case that this function
#'   generalizes.
#'
#' @keywords internal
unit_joint_beta <- function(u, spec, design, coef) {
  params <- spec@distrib@params
  if (!isTRUE(u$structural) && !isTRUE(u$mixed))
    return(unit_beta(u, coef, params))
  sst <- statmod_structural_state(design)
  if (is.null(sst)) return(numeric(0))
  if (isTRUE(u$structural)) {
    z <- sst$zeta[[u$term]]
    return(if (is.null(z)) numeric(0) else as.numeric(z[u$cols]))
  }
  z <- sst$zeta[[u$class$sterm]]
  flat <- unlist(coef[params], use.names = FALSE)
  bt <- numeric(length(u$joint))
  for (pc in u$pieces) {
    bt[match(pc$joint, u$joint)] <- if (isTRUE(pc$structural)) {
      as.numeric(z[pc$zcols])
    } else {
      flat[pc$index]
    }
  }
  bt
}


#' Every Penalty in a Model, Without the Design
#'
#' @description
#' Returns the same enumeration as [statmod_penalized()] without the column
#' positions, for the callers that need to know which penalties exist before
#' a design has been built.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A list of entries with `param`, `term`, `key`,
#'   `within` (positions among the term's own parameters), `penalty`, `fixed`,
#'   `n_values`, `values`, `min_ratio`, `search` and `ids`. The entry of a
#'   covariance class also carries `class` and `params`.
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
      # `term` says which of the design's structural terms the class is
      # addressed in, which is what every structural consumer looks its zeta
      # up by; for a class among the coefficients there is no such term and
      # the key stands in its place, as it did before there were any
      param = cl$pieces[[1L]]$param,
      term = if (cl$structural) cl$sterm else cl$key, key = cl$key,
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
#' Returns the groups of terms that share a covariance block: those carrying the
#' same label and the same grouping variable, collected across the equations.
#'
#' @details
#' # What identifies a class
#'
#' A class is identified by both a label and a grouping. The label comes from
#' [modelterms7::term_tag()] and the grouping from
#' [modelterms7::term_group()], which returns the expression, the levels and
#' the within-group column count. Two terms belong together only if the levels
#' agree as well as the expression: `droplevels(id)` and `id` are different
#' expressions for one grouping, and `id` under two subsets is one expression
#' for two groupings.
#'
#' Correlating effects on different groupings means nothing, since they are
#' indexed by different things and there is no block to estimate, so a label
#' used on two groupings is an error naming both.
#'
#' # Whose prior it is
#'
#' The joint prior belongs to the class and not to any of its members. At most
#' one term may name a `distrib`, and its dimension must be the class's total;
#' if no term names one, the default is a centered multivariate Gaussian on
#' [parameters7::dr_prod()], whose coordinates are the log standard deviations
#' and the correlations' angles, so a printed hyperparameter is the quantity
#' that it names. At a total of one column there is no correlation and the
#' default is the centered univariate Gaussian that a single unlabeled term
#' would have built.
#'
#' A class whose members carry priors of different families (a Gaussian
#' intercept correlated with a Student t slope) is a copula and not an
#' elliptical family. \pkg{modelterms7} rejects a univariate `distrib` on a
#' labeled term for that reason, so the case does not reach this function.
#'
#' @param terms The built terms, a named list of named lists: one element
#'   per distribution parameter, in the family's order.
#'
#' @return A list named by the class key (`<label> | <grouping>`), one element
#'   per class, each with `tag`, `group` (the deparsed expression), `levels`,
#'   `m`, `key`, `dim` (the class total), `space` (see [class_space()]),
#'   `structural`, `sterm` (the name of the structural term, present only
#'   where the class touches one), `penalty` and `pieces`. The pieces are one
#'   per member, as [label_pieces()] returns them, in the order in which the
#'   equations were walked. The list is empty if the terms carry no label.
#'
#' @seealso [statmod_penalty_keys()], which turns each into one penalized
#'   unit; [statmod_penalized()] for the interleaved index at which it is read.
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
    cl$space <- class_space(cl)
    cl$structural <- identical(cl$space, "zeta")
    # a class touching a structural term is addressed in that term's own
    # parameters, and there is at most one such term in a model, so every
    # structural piece of it names the same one
    st <- Filter(function(z) isTRUE(z$structural), cl$pieces)
    if (length(st)) cl$sterm <- st[[1L]]$term
    cl$penalty <- class_penalty(cl)
    cl
  })
}


#' The Vector in Which a Covariance Class Is Addressed
#'
#' @description
#' Returns which of the three vectors addresses the coordinates of a class.
#'
#' @details
#' A member written in an equation, or in the subformula of an additive term,
#' has columns in the stacked coefficient vector. A member inside a structural
#' term has none: its coordinates are that term's own parameters, which the
#' design carries in its structural state.
#'
#' A class whose members are all of one kind is read in that vector. A class
#' split between them is read in the joint vector that the fit already builds,
#' \eqn{[\beta; \zeta_{\mathrm{free}}]}: the inner step, the marginal
#' criterion and the variance all assemble it by placing the coefficients first
#' and a filter's free parameters after them. Such a class adds to those three
#' the cross block, which the penalty produces on its own as one Hessian over
#' the class's stacked vector, computed without reference to the split between
#' coefficients and filter parameters.
#'
#' @param cl One class, as [statmod_classes()] assembles it.
#'
#' @return One of `"coef"`, `"zeta"` or `"mixed"`.
#'
#' @seealso [statmod_classes()], its only caller; [class_joint_pieces()] for
#'   the positions at which a mixed class is read.
#'
#' @keywords internal
class_space <- function(cl) {
  st <- vapply(cl$pieces, function(z) isTRUE(z$structural), TRUE)
  if (all(st)) return("zeta")
  if (!any(st)) return("coef")
  "mixed"
}


#' Every Labeled Effect in a Term, Including Its Sub-Terms
#'
#' @description
#' Returns one piece per labeled random effect reachable from a term: the term
#' itself where it carries a label, and the sub-terms developing its own
#' parameters otherwise, walked to any depth.
#'
#' @details
#' # Why a sub-term is reachable at all
#'
#' `seg(x, psi ~ random(~ 1 | u | id))` develops a break-point over a labeled
#' random effect. The labeled term is not one of the equation's terms (the
#' equation carries one `SegTerm`, whose own label is absent), and its
#' coefficients are columns of that term's block. They are the columns that
#' [modelterms7::term_components()] reports as that component's `sub_index`,
#' so a piece records them as `within`, positions in the parent's
#' block, and the design turns them into positions in the stacked vector as it
#' does for any other term.
#'
#' A labeled effect written in a subformula of an additive term lives in
#' the same vector as one written in an equation, so the same positions apply.
#'
#' A structural parent addresses its coefficients differently: they are the
#' term's own parameters, which contribute no design column and live in the
#' design's structural state. The walk is the same and the positions that it
#' composes are the same numbers (for a labeled effect inside
#' `gas(kappa1 ~ 1 + random(~ 1 | u | g))` they are the `cols` at which the
#' unlabeled sub-term's own penalty is read), so a piece records which vector
#' the positions index and not a different arithmetic.
#'
#' # Depth
#'
#' A sub-term is an ordinary term and may develop parameters of its own, so the
#' walk recurses and composes `within` on the way down: a depth-two sub-term's
#' columns are positions in its parent's block, which are themselves positions
#' in the equation-level term's.
#'
#' A labeled term's own sub-terms are not walked. Its columns are the class's
#' already, and anything inside it belongs to that block rather than to another.
#'
#' @param term One built term.
#' @param param The distribution parameter to which its equation belongs.
#' @param nm The equation-level term's name, by which the design is keyed.
#' @param within The piece's columns in that term's block, or `NULL` at the top
#'   level, where the piece is the whole of it.
#' @param path Which of the parent's parameters were developed to reach the
#'   piece, outermost first, empty at the top level. It identifies the quantity
#'   on which the effect acts: a labeled random intercept written inside
#'   `seg(x, psi ~ ...)` is an effect on the break-point and not on the mean,
#'   and a report naming only the equation would be read as the second.
#' @param structural Whether `within` indexes a structural term's own
#'   parameters rather than columns of the design. `NULL` at the top level,
#'   where it is read from the term; a recursive call passes the value that it
#'   received.
#'
#' @return A list of pieces, each with `param`, `term`, `within`, `path`,
#'   `dim`, `tag`, `group` (as [modelterms7::term_group()] returns it),
#'   `structural` and `distrib`. The list is empty if nothing under the term
#'   is labeled.
#'
#' @seealso [statmod_classes()], its caller; [class_pieces()] for the mapping
#'   of `within` onto the stacked vector.
#'
#' @keywords internal
label_pieces <- function(term, param, nm, within = NULL,
                        path = character(0), structural = NULL) {
  # WHICH VECTOR THE PIECE'S POSITIONS ARE POSITIONS IN, which is a property
  # of the equation-level term and not of the sub-term carrying the label.
  # It is read once at the top of the walk and carried down: a labelled
  # effect inside a filter is addressed among that term's own parameters,
  # where the design has no column at all.
  if (is.null(structural)) {
    structural <- S7::S7_inherits(term, modelterms7::structural_term)
  }
  tg <- tryCatch(modelterms7::term_tag(term), error = function(e) NA_character_)
  if (length(tg) == 1L && !is.na(tg)) {
    gr <- modelterms7::term_group(term)
    if (is.null(gr)) {
      stop(sprintf(paste0(
        "'%s' in '%s' carries the covariance label '%s' but reports no\n",
        "  grouping variable, and effects are correlated within a grouping."),
        nm, param, tg), call. = FALSE)
    }
    return(list(list(param = param, term = nm, within = within, path = path,
                     dim = gr$dim, tag = tg, group = gr,
                     structural = structural,
                     distrib = tryCatch(term@distrib, error = function(e) NULL))))
  }
  out <- list()
  comp <- tryCatch(modelterms7::term_components(term), error = function(e) list())
  for (cp in comp) {
    for (k in seq_along(cp$subs)) {
      w <- cp$sub_index[[k]]
      if (!is.null(within)) w <- within[w]
      out <- c(out, label_pieces(cp$subs[[k]], param, nm, w,
                                 c(path, cp$name), structural))
    }
  }
  out
}


#' The Joint Prior of a Covariance Class
#'
#' @description
#' Returns the penalty over a class's stacked coefficients: the distribution
#' that one of its members named, or the centered multivariate Gaussian that is
#' the default.
#'
#' @details
#' The prior describes the effects of one group over every column that the label
#' collects, so its dimension is the class's total and the penalty covers
#' \eqn{m} such blocks. The default prior is a centered multivariate Gaussian
#' whose covariance is written by [parameters7::dr_prod()], where a coordinate
#' is the logarithm of a standard deviation exactly.
#'
#' Naming a `distrib` on more than one member is an error and not a precedence
#' rule: the prior is one object, so two of them cannot both be used.
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
        "  group over every column that the label collects."),
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
#' Returns one key per hyperparameter that a term fixed in its constructor,
#' composed of `parameter`, the penalty's key and the hyperparameter's name
#' joined by a carriage return.
#'
#' @details
#' Which hyperparameters are estimated is a property of the terms and not of
#' the criterion, because the term is where the penalty is named.
#' [path_rows()] and the summary's account of what was estimated and what was
#' held read this enumeration. The outer index reads the same `fixed` entries
#' of [statmod_penalized()].
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
#' Returns the key under which a hyperparameter row is filed: the term's name
#' where the term carries one penalty or the entry has no name, and
#' `term::entry` where the term carries several penalties and the entry is
#' named.
#'
#' @details
#' The composition is written once so that every caller reading a
#' hyperparameter by key uses the same key. A term carrying one penalty is
#' keyed by its name without the entry's name, whatever the entry is called.
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
#' Returns the entry of [statmod_penalized()] that a hyperparameter row names,
#' or `NULL` where there is none.
#'
#' @details
#' A reader that starts from a `(parameter, term)` pair looks the penalty up
#' here and not from the term's whole block, which assumes one penalty per
#' term. The lookup returns that penalty where the term has one, and the
#' penalty with the given key where it has several.
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
#' Returns the grid size that the term specifies, or the fallback of
#' [path_fallbacks()] where the term specifies none.
#'
#' @details
#' How finely a hyperparameter is swept is a property of the term for the
#' same reason as whether it is swept at all: a penalized block of four
#' columns and one of four hundred require different grids. A criterion applies
#' to every hyperparameter of the model at once, including the smooth ones, and
#' those are read at the mode instead of being swept, so the grid size cannot
#' be a property of the criterion. A term gives its grid size through
#' [modelterms7::term_grid()], and the value travels with the penalty's entry,
#' so one reached through a sub-term of a structural term carries it too.
#'
#' @param spec A [StatmodSpec()].
#' @param row One row of [path_rows()]'s index.
#' @param default The fallback of [path_fallbacks()], used for a term that
#'   specifies none.
#'
#' @return A single integer.
#'
#' @seealso [statmod_held()], [path_values()]
#'
#' @keywords internal
statmod_grid_size <- function(spec, row, default) {
  statmod_path_setting(spec, row, "n_values", default, row$name)
}


#' The Depth of a Path for One Term
#'
#' @description
#' Returns the depth that the term specifies, or the fallback of
#' [path_fallbacks()] where the term specifies none.
#'
#' @details
#' The depth is one number per term and not one per hyperparameter, because
#' only the sweep by kink size reads it. A bounded hyperparameter is swept over
#' its own interval, and a shape that does not move the kink over a grid whose
#' distances above the smallest admissible shape are geometric; in both sweeps
#' the length does not come from this function.
#'
#' @param spec A [StatmodSpec()].
#' @param row One row of [path_rows()]'s index.
#' @param default The fallback of [path_fallbacks()], used for a term that
#'   specifies none.
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
#' Returns `"grid"` for every combination of the term's kinked hyperparameters,
#' `"cyclic"` for one at a time, or the default where the term specifies
#' neither.
#'
#' @details
#' The search setting belongs to the term, for the reason that the whole
#' enumeration does: a
#' criterion applies to every hyperparameter of the model, and a smooth one is
#' read at the mode instead of being swept, so a setting that concerns only
#' the kinked hyperparameters cannot belong to the criterion. A penalty with a
#' kink is fitted by a scheme of its own, and how that scheme covers the
#' term's own hyperparameters is part of the scheme.
#'
#' Because the setting is per term, one term's choice does not change
#' another's: `y ~ lasso(X) + enet(R, search = "cyclic")` sweeps the elastic
#' net one coordinate at a time and does not change how the lasso is swept.
#'
#' @param spec A [StatmodSpec()].
#' @param row One row of [path_rows()]'s index.
#' @param default What to use where the term specifies nothing.
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


#' The Values That a Term Wrote Out for One Hyperparameter
#'
#' @description
#' Returns the grid that the term specifies, or `NULL` where the path builds
#' one.
#'
#' @details
#' This is the third state of a hyperparameter's argument, beside holding it at
#' one number and leaving it to be estimated over a grid that the path
#' constructs. It
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
#' @description
#' Reads one setting of the path from the penalty's entry in the term, and
#' returns `default` where the term does not set it.
#'
#' @param spec A [StatmodSpec()].
#' @param row One row of [path_rows()]'s index.
#' @param field Which field of the penalty's entry to read.
#' @param default The fallback of [path_fallbacks()] for this setting, used
#'   where the term specifies none.
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
