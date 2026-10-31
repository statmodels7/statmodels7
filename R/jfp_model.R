#' @include jfp_engine.R
NULL

#' The Constant Direction and the Power of Each Jeffreys Unit
#'
#' @description
#' Completes the units of [statmod_penalized()] that carry a
#' [modelterms7::JeffreysPenalty()]: each receives `power`, the power of its
#' determinant, and `const`, whether the constant direction of its equation
#' enters the determinant.
#'
#' @details
#' The constant direction follows the rule of [jfp_constant()] with the
#' setting `intercept` of the term. The level's columns are the columns of
#' the equation that no penalized unit covers, read from the design at which
#' the units are enumerated.
#'
#' @param units The units, as [statmod_penalized_units()] builds them.
#' @param spec A [StatmodSpec()].
#' @param design The design.
#'
#' @return `units`, with `power` and `const` set on each Jeffreys unit.
#'
#' @seealso [jfp_blocks()]
#'
#' @keywords internal
jfp_resolve_units <- function(units, spec, design) {
  jf <- which(vapply(units, function(u) isTRUE(u$jeffreys), logical(1)))
  if (!length(jf)) return(units)
  for (k in jf) {
    u <- units[[k]]
    p <- u$param
    covered <- integer(0)
    for (v in units) {
      if (!is.null(v$pieces)) {
        for (pc in v$pieces) {
          if (identical(pc$param, p)) covered <- c(covered, pc$cols)
        }
      } else if (identical(v$param, p) && !isTRUE(v$structural)) {
        covered <- c(covered, v$cols)
      }
    }
    units[[k]]$power <- u$penalty@power
    if (!S7::S7_inherits(spec@terms[[p]][[u$term]],
                         modelterms7::JeffreysTerm)) {
      units[[k]] <- c(units[[k]], jfp_subformula_unit(u, spec, design))
      next
    }
    units[[k]]$kind <- "static"
    X <- design[[p]]$X
    level <- setdiff(seq_len(ncol(X)), covered)
    units[[k]]$const <- jfp_constant(
      u$penalty@intercept,
      if (length(level)) X[, level, drop = FALSE] else NULL,
      X[, u$cols, drop = FALSE])
  }
  units
}


#' A Jeffreys Block in the Subformula of a Refreshed Term
#'
#' @description
#' Resolves a Jeffreys unit that a term such as [modelterms7::nl()] or a
#' smoothed [modelterms7::seg()] carries up from the subformula of one of its
#' parameters: the term's columns from which the block's rows are combined,
#' and whether the constant direction of the subformula enters them.
#'
#' @details
#' The level of the subformula is its unpenalized sub-terms, read from
#' [modelterms7::term_components()]. The constant direction is the derivative
#' of the predictor along a constant shift of the parameter, which is the
#' combination \eqn{X_{\rm level}\gamma} of the level's refreshed columns
#' with \eqn{Z_{\rm level}\gamma = 1}; the span is tested on the
#' subformula's own design, by [jfp_constant()].
#'
#' A term whose block is a working linearization (a sharp
#' [modelterms7::jump()] or [modelterms7::jseg()]) is rejected, and so is the
#' break-point subformula of a sharp [modelterms7::seg()], whose rows jump
#' when the break-point crosses an observation. A parameter whose own
#' development moves is rejected as well.
#'
#' @param u The unit.
#' @param spec A [StatmodSpec()].
#' @param design The design.
#'
#' @return A list with `kind` (`"refresh"`), `const`, `basis` and `C`, as
#'   [jfp_block()] takes them.
#'
#' @keywords internal
jfp_subformula_unit <- function(u, spec, design) {
  p <- u$param
  tm <- spec@terms[[p]][[u$term]]
  rf <- Filter(function(r) identical(r$param, p) && identical(r$term, u$term),
               attr(design, "refresh"))
  no <- function(why) {
    stop(sprintf("jfp() inside the subformula of '%s' is not available: %s.",
                 u$term, why), call. = FALSE)
  }
  if (!length(rf)) no("the term's block is not recomputed at the coefficients")
  if (isTRUE(rf[[1L]]$frozen)) {
    no("the term's block is a working linearization with a frozen weight")
  }
  comp <- modelterms7::term_components(tm)
  hit <- NULL
  for (cp in comp) {
    for (j in seq_along(cp$sub_index)) {
      if (identical(as.integer(cp$sub_index[[j]]), as.integer(u$within))) {
        hit <- list(cp = cp, j = j)
      }
    }
  }
  if (is.null(hit)) no("its columns are not one sub-term of one parameter")
  cp <- hit$cp
  if (isTRUE(tm@build_info$mobile[[cp$name]])) {
    no(sprintf("the development of '%s' moves with the coefficients",
               cp$name))
  }
  if (S7::S7_inherits(tm, modelterms7::SegTerm) &&
      startsWith(cp$name, "psi") && is.null(tm@spec$smoothed)) {
    no("the rows of a sharp break-point's subformula jump at the observations")
  }
  lev <- which(vapply(seq_along(cp$subs), function(j) {
    j != hit$j && !length(modelterms7::term_penalties(cp$subs[[j]]))
  }, logical(1)))
  ZA <- as_dense(modelterms7::term_matrix(cp$subs[[hit$j]]))
  Zl <- if (length(lev)) {
    do.call(cbind, lapply(cp$subs[lev], function(s)
      as_dense(modelterms7::term_matrix(s))))
  } else NULL
  const <- jfp_constant(u$penalty@intercept, Zl, ZA)
  tcols <- design[[p]]$blocks[[u$term]]
  pA <- length(u$within)
  if (!const) {
    return(list(kind = "refresh", const = FALSE,
                basis = tcols[u$within], C = diag(pA)))
  }
  li <- unlist(cp$sub_index[lev], use.names = FALSE)
  gam <- qr.coef(qr(Zl), rep(1, nrow(Zl)))
  gam[is.na(gam)] <- 0
  C <- matrix(0, length(li) + pA, pA + 1L)
  C[seq_along(li), 1L] <- gam
  C[length(li) + seq_len(pA), 1L + seq_len(pA)] <- diag(pA)
  list(kind = "refresh", const = TRUE, basis = tcols[c(li, u$within)],
       C = C)
}


#' The Jeffreys Blocks of a Model
#'
#' @description
#' Returns one [jfp_block()] per Jeffreys unit of the model, in the order of
#' [statmod_penalized()].
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param keys `NULL` for every Jeffreys unit, or the keys of the units to
#'   keep.
#'
#' @return A list of [jfp_block()] objects, empty where the model carries no
#'   Jeffreys term.
#'
#' @keywords internal
jfp_blocks <- function(spec, design, keys = NULL) {
  out <- list()
  for (u in statmod_penalized(spec, design)) {
    if (!isTRUE(u$jeffreys)) next
    if (!is.null(keys) && !u$key %in% keys) next
    out[[length(out) + 1L]] <- if (identical(u$kind, "refresh")) {
      jfp_block(u$param, u$cols, power = u$power, const = u$const,
                kind = "refresh", term = u$term, basis = u$basis, C = u$C)
    } else {
      jfp_block(u$param, u$cols, power = u$power, const = u$const)
    }
  }
  out
}


#' The Method of the Expected Information That a Jeffreys Penalty Reads
#'
#' @description
#' Returns `NULL` where the family has the closed form of its expected
#' information and of its derivatives, and
#' [distributions7::expect_integrate()] where it has only its density.
#'
#' @details
#' The penalty is the Jeffreys prior only on the expected information, so the
#' method of the fit is not read. The outer product of the scores at the data,
#' which [distributions7::expect_opg()] computes, is not the expected
#' information, and a family that reaches it is rejected; so is a family on
#' which [distributions7::expect_integrate()] gives no first and second
#' derivative of the matrix.
#'
#' @param distrib A [distributions7::distrib] object.
#'
#' @return `NULL` or an [distributions7::expect_integrate()] object.
#'
#' @keywords internal
jfp_emethod <- function(distrib) {
  m <- resolve_emethod(distrib, NULL)
  if (is.null(m)) return(NULL)
  ei <- distributions7::expect_integrate()
  if (isTRUE(expected_deriv_ok(distrib, ei)) &&
      isTRUE(expected_deriv2_ok(distrib, ei))) {
    return(ei)
  }
  stop(paste0("jfp() needs the expected information of the family and its ",
              "first two derivatives,\n  in closed form or by ",
              "distributions7::expect_integrate(), and this family has ",
              "neither."), call. = FALSE)
}


#' The Jeffreys State at Given Coefficients, Kept Between Calls
#'
#' @description
#' Returns [jfp_state()] for the Jeffreys blocks of the model at `coef`, or
#' `NULL` where the model carries no Jeffreys term. The objective, its
#' gradient and its curvature read the same point in turn, so the last state
#' is kept on the design's memo.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef A named list of coefficient vectors.
#' @param keys As in [jfp_blocks()]; a state over some of the units is not
#'   kept.
#'
#' @return An environment, or `NULL`.
#'
#' @keywords internal
jfp_state_at <- function(spec, design, coef, keys = NULL) {
  blocks <- jfp_blocks(spec, design, keys)
  if (!length(blocks)) return(NULL)
  mm <- if (is.null(keys)) attr(design, "eta_memo") else NULL
  key <- list(unlist(coef[spec@distrib@params], use.names = FALSE),
              statmod_structural_state(design)$zeta,
              spec@offsets, spec@weights, spec@n_obs)
  if (!is.null(mm) && !is.null(mm$jfp_key) && identical(mm$jfp_key, key)) {
    return(mm$jfp_state)
  }
  st <- jfp_state(spec, design, coef, blocks,
                  emethod = jfp_emethod(spec@distrib))
  if (!is.null(mm)) {
    mm$jfp_key <- key
    mm$jfp_state <- st
  }
  st
}


#' The Jeffreys Penalty of a Model, With the Sign of a Penalty
#'
#' @description
#' Computes \eqn{-\sum_b a_b\log|\tilde I_b|}, the quantity that the
#' objective adds to the negative log-likelihood, or its gradient or Hessian
#' over the stacked coefficients.
#'
#' @details
#' Where a block's information is not positive definite the value is
#' `Inf`, so a line search rejects the point, and the gradient and the
#' Hessian signal an error.
#'
#' `"parts"` returns the Hessian split as [jfp_hessian()] gives it, with the
#' sign of the penalty: the objective's curvature gains
#' [statmod_information_at()] at `H = -weights`, plus `-refresh`, plus
#' `tcrossprod(lowrank)`.
#'
#' Beside a filter the state spans the joint vector of the coefficients and
#' the filter's free parameters. With `joint = FALSE` the result is then
#' restricted to the coefficients, which is what the objective over the
#' coefficients reads; [jfp_joint_extra()] returns the rest.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef A named list of coefficient vectors.
#' @param what `"value"`, `"gradient"`, `"hessian"` or `"parts"`.
#' @param keys As in [jfp_blocks()].
#' @param joint Whether to return the whole joint vector beside a filter.
#'
#' @return `0`, a zero vector, `NULL` or `NULL` where the model carries no
#'   Jeffreys term; otherwise a number, a numeric vector over the stacked
#'   coefficients, a base matrix, or the list of [jfp_hessian()].
#'
#' @keywords internal
jfp_penalty <- function(spec, design, coef,
                        what = c("value", "gradient", "hessian", "parts"),
                        keys = NULL, joint = FALSE) {
  what <- match.arg(what)
  st <- jfp_state_at(spec, design, coef, keys)
  if (is.null(st)) {
    return(switch(what, value = 0,
                  gradient = numeric(sum(vapply(design, function(d) d$npar,
                                                integer(1)))),
                  NULL))
  }
  if (what == "value") {
    v <- jfp_value(st)
    return(if (is.finite(v)) -v else Inf)
  }
  ib <- if (is.null(st$joint) || joint) seq_len(st$total) else
    seq_len(st$nb)
  if (what == "gradient") return(-jfp_gradient(st)[ib])
  h <- jfp_hessian(st, assemble = what == "hessian")
  if (what == "hessian") return(-h$hessian[ib, ib, drop = FALSE])
  if (!is.null(st$joint) && !joint) {
    h$lowrank <- h$lowrank[ib, , drop = FALSE]
    if (!is.null(h$refresh)) h$refresh <- h$refresh[ib, ib, drop = FALSE]
  }
  h
}


#' The Jeffreys Penalty Over a Filter's Own Parameters
#'
#' @description
#' Returns the part of the gradient or of the Hessian of [jfp_penalty()] that
#' the objective over the coefficients does not carry: the entries in the
#' filter's free parameters, and for the Hessian the cross block, with the
#' block over the coefficients set to zero.
#'
#' @details
#' A joint assembly that places [statmod_penalty_at()] on the coefficients
#' adds this to complete the Jeffreys penalty over the joint vector. It is
#' zero, of the requested size, where the model carries no Jeffreys term or
#' no filter.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef A named list of coefficient vectors.
#' @param what `"gradient"` or `"hessian"`.
#' @param n The length of the joint vector.
#'
#' @return A numeric vector of length `n`, or an `n x n` matrix.
#'
#' @keywords internal
jfp_joint_extra <- function(spec, design, coef, what = c("gradient", "hessian"),
                            n) {
  what <- match.arg(what)
  st <- jfp_state_at(spec, design, coef)
  out <- if (what == "gradient") numeric(n) else matrix(0, n, n)
  if (is.null(st) || is.null(st$joint)) return(out)
  if (st$total != n) {
    stop(sprintf(paste0("The Jeffreys penalty spans %d coordinates and the ",
                        "joint vector %d."), st$total, n), call. = FALSE)
  }
  ib <- seq_len(st$nb)
  if (what == "gradient") {
    out <- jfp_penalty(spec, design, coef, "gradient", joint = TRUE)
    out[ib] <- 0
  } else {
    out <- jfp_penalty(spec, design, coef, "hessian", joint = TRUE)
    out[ib, ib] <- 0
  }
  out
}


#' How the Jeffreys Curvature Moves With the Coefficients
#'
#' @description
#' Returns the derivative of the Jeffreys penalty's Hessian along a
#' direction \eqn{v}, or its second derivative along \eqn{v} and \eqn{w},
#' with the sign of a penalty: \eqn{-\sum_b a_b D_3(v)} and
#' \eqn{-\sum_b a_b D_4(v, w)} of [jfp_deriv3_dir()]. With `M` the trace
#' against it is returned instead of the matrix.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef A named list of coefficient vectors.
#' @param v,w Directions over the coefficients; `w = NULL` for the first
#'   derivative.
#' @param M A symmetric matrix over the coefficients, or `NULL`.
#'
#' @return A base matrix, or a single number when `M` is given; `NULL`
#'   (or `0` when `M` is given) where the model carries no Jeffreys term.
#'
#' @seealso [statmod_penalty_dbeta()], [statmod_penalty_second()]
#'
#' @keywords internal
jfp_penalty_move <- function(spec, design, coef, v, w = NULL, M = NULL) {
  st <- jfp_state_at(spec, design, coef)
  if (is.null(st)) return(if (is.null(M)) NULL else 0)
  if (is.null(w)) {
    if (is.null(M)) -jfp_deriv3_dir(st, v) else -jfp_deriv3(st, M, v)
  } else {
    if (is.null(M)) -jfp_deriv4_dir(st, v, w) else -jfp_deriv4(st, M, v, w)
  }
}


#' Whether the Outer Derivatives Can Read a Jeffreys Penalty
#'
#' @description
#' Returns `TRUE` where the exact outer gradient (`order = 1`) or Hessian
#' (`order = 2`) can carry the movement of a Jeffreys penalty.
#'
#' @details
#' The gradient reads the third derivative of the penalty, so the third
#' derivative of the family's expected information, and the Hessian the
#' fourth; each is probed once at a random point with the method that
#' [jfp_emethod()] selects. The derivatives are written for a model without a
#' filter and without a term whose block is recomputed at the coefficients,
#' and the function returns `FALSE` there.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param order `1` or `2`.
#'
#' @return A single logical.
#'
#' @keywords internal
jfp_outer_ok <- function(spec, design, order = 1L) {
  if (length(attr(design, "structural"))) return(FALSE)
  rf <- attr(design, "refresh")
  if (!is.null(rf) && length(rf)) return(FALSE)
  d <- spec@distrib
  mem <- .jfp_outer_memory
  if (!is.null(mem$distrib) && identical(mem$distrib, d) &&
      !is.null(mem$value[[as.character(order)]])) {
    return(mem$value[[as.character(order)]])
  }
  if (!identical(mem$distrib, d)) {
    mem$distrib <- d
    mem$value <- list()
  }
  # the probe leaves the caller's random stream as it was
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(20261031L)
  ans <- jfp_outer_probe(d, order)
  mem$value[[as.character(order)]] <- ans
  ans
}

# The one-slot memory of jfp_outer_ok(): the family last probed and the
# answer per order.
.jfp_outer_memory <- new.env(parent = emptyenv())

#' @rdname jfp_outer_ok
#' @param d A distribution.
#' @keywords internal
jfp_outer_probe <- function(d, order) {
  em <- tryCatch(jfp_emethod(d), error = function(e) FALSE)
  if (isFALSE(em)) return(FALSE)
  th <- tryCatch(distributions7::generate_random_theta(d),
                 error = function(e) NULL)
  if (is.null(th)) return(FALSE)
  y <- tryCatch(distributions7::distrib_rng(d, 1L, th),
                error = function(e) NULL)
  if (is.null(y)) return(FALSE)
  isTRUE(tryCatch({
    distributions7::distrib_d3expected_hessian(d, y, th, scale = "link",
                                               method = em)
    if (order >= 2L) {
      distributions7::distrib_d4expected_hessian(d, y, th, scale = "link",
                                                 method = em)
    }
    TRUE
  }, error = function(e) FALSE))
}
