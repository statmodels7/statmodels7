#' @include predict_random.R
NULL

#' The Random-Effect Terms Written Inside a Subformula
#'
#' @description
#' Lists the [modelterms7::random()] sub-terms that develop a parameter of
#' another term, such as `Asym ~ 1 + random(~ 1 | Tree)` inside
#' [modelterms7::nl()] or `psi ~ random(~ 1 | id)` inside
#' [modelterms7::seg()].
#'
#' @details
#' The terms are read through [modelterms7::term_components()], which every
#' term carrying developed parameters answers. The key of a row is the outer
#' term's key, the parameter and the sub-term joined by `"::"`, which is the
#' key [hyper()] reports where a term carries more than one penalty.
#'
#' @param spec The fit's specification.
#'
#' @return A data frame with columns `param`, `term`, `comp`, `s` (the
#'   sub-term's place in the parameter's subformula) and `key`, one row per
#'   nested random-effect sub-term, or zero rows.
#'
#' @keywords internal
nested_random_terms <- function(spec) {
  rows <- list()
  for (p in names(spec@terms)) {
    for (k in names(spec@terms[[p]])) {
      cmp <- tryCatch(modelterms7::term_components(spec@terms[[p]][[k]]),
                      error = function(e) list())
      for (cn in names(cmp)) {
        subs <- cmp[[cn]]$subs
        for (s in seq_along(subs)) {
          if (S7::S7_inherits(subs[[s]], modelterms7::RandomTerm)) {
            rows[[length(rows) + 1L]] <- data.frame(
              param = p, term = k, comp = cn, s = s,
              key = paste(k, cn, names(subs)[s], sep = "::"),
              stringsAsFactors = FALSE)
          }
        }
      }
    }
  }
  if (length(rows)) do.call(rbind, rows) else
    data.frame(param = character(0), term = character(0), comp = character(0),
               s = integer(0), key = character(0), stringsAsFactors = FALSE)
}


#' A Term With Its Nested Random Effects Read in One Group's Columns
#'
#' @description
#' Returns a copy of a term carrying developed parameters in which every
#' row reads the within-group design of the listed random-effect sub-terms
#' in the columns of their first level, and nothing in the columns of the
#' other levels.
#'
#' @details
#' A random effect set aside contributes \eqn{z_i^\top u} at every row, with
#' \eqn{z_i} the row's within-group design and \eqn{u} the value it is read
#' at: zero for `random = "zero"`, a node of the prior for `"marginal"`.
#' Inside a nonlinear term the effect enters the parameter, and the
#' contribution of the term is not linear in it, so the columns cannot
#' simply be set to zero as they are for a term written in an equation. With
#' every row in the first level, writing \eqn{u} into that level's
#' coefficients gives every row the effect \eqn{z_i^\top u}, and the term is
#' evaluated at those coefficients by its own methods. A group the fit never
#' saw is read the same way, so it needs no effect of its own.
#'
#' Two places are rewritten. At the fitting rows the term reads the stored
#' design of the developed parameter, `blueprint$Z`, whose columns for the
#' sub-term become the within-group rows (the sum of the group's columns, a
#' row being non-zero in its own group's alone). At new rows the term
#' reapplies each sub-term through [modelterms7::term_predict()], so the
#' sub-term's grouping expression becomes the first level repeated.
#'
#' @param tm A built term carrying developed parameters.
#' @param rows The rows of [nested_random_terms()] for this term.
#' @param coef The term's coefficients to record as the ones it was last
#'   committed at, so that its block at new rows is the Jacobian there.
#'
#' @return The modified term.
#'
#' @keywords internal
pool_nested_random <- function(tm, rows, coef) {
  bp <- tm@blueprint
  cmp <- modelterms7::term_components(tm)
  for (i in seq_len(nrow(rows))) {
    cn <- rows$comp[i]
    s <- rows$s[i]
    sub <- bp$subs[[cn]][[s]]
    gr <- modelterms7::term_group(sub)
    m <- length(gr$levels)
    d <- gr$dim
    pos <- match(cmp[[cn]]$sub_index[[s]], cmp[[cn]]$index)
    Zc <- bp$Z[[cn]]
    if (!is.null(Zc)) {
      W <- as.matrix(Zc[, pos, drop = FALSE] %*%
                       kronecker(rep(1, m), diag(d)))
      P <- matrix(0, nrow(W), length(pos))
      P[, seq_len(d)] <- W
      Zc <- as.matrix(Zc)
      Zc[, pos] <- P
      bp$Z[[cn]] <- Zc
    }
    sb <- sub@blueprint
    v1 <- all.vars(sb$gexpr)[1L]
    sb$gexpr <- bquote(rep(.(as.character(sb$glevels[1L])),
                           length(.(as.name(v1)))))
    sub@blueprint <- sb
    bp$subs[[cn]][[s]] <- sub
  }
  bp$coef <- as.numeric(coef)
  tm@blueprint <- bp
  tm
}


#' The Coordinates of the Nested Random Effects Set Aside
#'
#' @description
#' For each row of [nested_random_terms()], the positions of the sub-term's
#' coefficients in its equation's design, and those of the first level.
#'
#' @param spec,design The specification and its design.
#' @param rows Rows of [nested_random_terms()].
#'
#' @return A list, one entry per row, with `param`, `cols` and `first`.
#'
#' @keywords internal
nested_columns <- function(spec, design, rows) {
  lapply(seq_len(nrow(rows)), function(i) {
    p <- rows$param[i]
    tm <- spec@terms[[p]][[rows$term[i]]]
    cmp <- modelterms7::term_components(tm)
    idx <- design[[p]]$blocks[[rows$term[i]]][
      cmp[[rows$comp[i]]]$sub_index[[rows$s[i]]]]
    d <- modelterms7::term_group(
      cmp[[rows$comp[i]]]$subs[[rows$s[i]]])$dim
    list(param = p, cols = idx, first = idx[seq_len(d)])
  })
}


#' The Prior of One Nested Random-Effect Term
#'
#' @description
#' The Gaussian prior of a random-effect sub-term, as [random_prior()]
#' returns it for a term written in an equation.
#'
#' @details
#' The penalty unit is the outer term's, keyed by the outer term where it
#' carries one penalty and by the row's key where it carries several. Only a
#' Gaussian prior is integrated: inside a nonlinear term the average is taken
#' over the nodes of a product grid at which the term is evaluated, and a
#' prior of another family signals an error naming `random = "zero"`.
#'
#' @param spec,design The specification and its design.
#' @param fit The fitted model, for the hyperparameters.
#' @param row One row of [nested_random_terms()].
#'
#' @return A list with `gaussian = TRUE`, `dim` and `chol`.
#'
#' @keywords internal
nested_prior <- function(spec, design, fit, row) {
  squash <- function(x) gsub("[[:space:]]", "", x)
  units <- Filter(function(u) identical(u$param, row$param) &&
                    identical(u$term, row$term),
                  statmod_penalized(spec, design))
  if (length(units) > 1L) {
    units <- Filter(function(u) identical(squash(u$key), squash(row$key)),
                    units)
  }
  if (length(units) != 1L) {
    stop(sprintf("'%s' has no single prior to integrate over.", row$key),
         call. = FALSE)
  }
  u <- units[[1L]]
  sub <- modelterms7::term_components(
    spec@terms[[row$param]][[row$term]])[[row$comp]]$subs[[row$s]]
  d <- modelterms7::term_group(sub)$dim
  th <- as.list(fit@hyper[[row$param]][[u$key]])
  pen <- u$penalty
  if (!isTRUE(penalties7::beta_quadratic(pen, th))) {
    stop(sprintf(paste0("The prior of '%s' is not Gaussian, and inside a ",
                        "term's subformula only a\n  Gaussian prior is ",
                        "averaged over. random = \"zero\" gives the typical ",
                        "group."), row$key), call. = FALSE)
  }
  H <- as.matrix(penalties7::penalty_hessian(pen, rep(0, pen@n_coef), th))
  Om <- H[seq_len(d), seq_len(d), drop = FALSE]
  list(gaussian = TRUE, dim = d, chol = t(chol(solve(Om))))
}


#' The Design of the Terms Holding Nested Effects, at Other Coefficients
#'
#' @description
#' Recomputes the columns of every term holding a nested random effect at
#' the coefficients given, which is the Jacobian there for a term whose
#' block is one.
#'
#' @param spec,design The specification and its design.
#' @param coef The coefficients, a named list.
#' @param terms A data frame with columns `param` and `term`.
#'
#' @return The design with those columns replaced.
#'
#' @keywords internal
nested_design_at <- function(spec, design, coef, terms) {
  for (i in seq_len(nrow(terms))) {
    p <- terms$param[i]
    k <- terms$term[i]
    cols <- design[[p]]$blocks[[k]]
    tm <- spec@terms[[p]][[k]]
    b <- coef[[p]][cols]
    Xt <- if (is.null(spec@newdata)) {
      modelterms7::term_matrix(modelterms7::term_refresh(tm, b))
    } else {
      bp <- tm@blueprint
      bp$coef <- as.numeric(b)
      tm@blueprint <- bp
      modelterms7::term_predict(tm, newdata = spec@newdata)
    }
    X <- design[[p]]$X
    if (isS4(X)) X <- as.matrix(X)
    X[, cols] <- as.matrix(Xt)
    design[[p]]$X <- X
  }
  design
}


#' A Prediction With Nested Random Effects Set Aside
#'
#' @description
#' Prepares the specification, the design and the coefficients of a
#' prediction in which random effects written inside a subformula are read
#' at zero or averaged over their prior.
#'
#' @details
#' Every term holding such an effect is replaced by the copy
#' [pool_nested_random()] returns, its effects' coefficients are set to zero
#' and its columns are recomputed there. The columns of the effects
#' themselves are then set to zero, so that the delta method of the
#' prediction carries the uncertainty of the other coefficients alone, as it
#' does for a term written in an equation.
#'
#' @param fit The fitted model.
#' @param spec The specification at the prediction's rows.
#' @param nest The rows of [nested_random_terms()] set aside, with `mode`.
#' @param unseen As in [statmod_design()].
#'
#' @return A list with `spec`, `design`, `coef` and `cols`, the latter from
#'   [nested_columns()].
#'
#' @keywords internal
nested_prepare <- function(fit, spec, nest, unseen) {
  coef <- fit@coefficients
  # the positions of the coefficients, read at the fitting rows, where a
  # group the new rows add cannot stop the build
  design0 <- statmod_design(fit@spec)
  cols <- nested_columns(spec, design0, nest)
  for (cl in cols) coef[[cl$param]][cl$cols] <- 0
  terms <- unique(nest[, c("param", "term")])
  for (i in seq_len(nrow(terms))) {
    p <- terms$param[i]
    k <- terms$term[i]
    rows <- nest[nest$param == p & nest$term == k, , drop = FALSE]
    spec@terms[[p]][[k]] <- pool_nested_random(
      spec@terms[[p]][[k]], rows, coef[[p]][design0[[p]]$blocks[[k]]])
  }
  design <- statmod_design(spec, unseen)
  design <- nested_design_at(spec, design, coef, terms)
  for (cl in cols) {
    X <- design[[cl$param]]$X
    X[, cl$cols] <- 0
    design[[cl$param]]$X <- X
  }
  list(spec = spec, design = design, coef = coef, cols = cols,
       terms = terms)
}


#' The Standard Error of a Marginal Parameter With Nested Effects
#'
#' @description
#' The delta method of \eqn{\bar\theta = \sum_k w_k h^{-1}(\eta(\beta,
#' u_k))}, whose gradient in the coefficients is
#' \eqn{\sum_k w_k\, (h^{-1})'(\eta_k)\, J_k} with \eqn{J_k} the design at
#' the coefficients of node \eqn{k}.
#'
#' @param fit The fitted model.
#' @param pr What [nested_prepare()] returns.
#' @param coef_at A function of the node index returning the coefficients
#'   there.
#' @param eta_at A function of the node index returning the predictors there.
#' @param w The nodes' weights.
#' @param p The parameter.
#' @param fitv The marginal parameter.
#' @param level The confidence level.
#' @param ... Passed to [vcov.StatmodFit()].
#'
#' @return A data frame with `fit`, `se`, `lower` and `upper`; the interval
#'   is the delta method's on the link scale, carried through the link.
#'
#' @keywords internal
nested_marginal_se <- function(fit, pr, coef_at, eta_at, w, p, fitv,
                               level = 0.95, ...) {
  V <- vcov(fit, readable = FALSE, ...)
  spec <- pr$spec
  d <- pr$design[[p]]
  n <- spec@n_obs
  key <- paste(p, d$coef_names, sep = ":")
  G <- matrix(0, n, length(key))
  lk <- spec@distrib@link_params[[p]]
  zero_cols <- unlist(lapply(pr$cols, function(cl)
    if (identical(cl$param, p)) cl$cols else integer(0)))
  for (k in seq_along(w)) {
    dk <- nested_design_at(spec, pr$design, coef_at(k), pr$terms)
    Xk <- as.matrix(dk[[p]]$X)
    Xk[, zero_cols] <- 0
    G <- G + w[k] * as.numeric(linkfunctions7::dlinkinv(lk, eta_at(k)[[p]])) * Xk
  }
  v <- rep(NA_real_, n)
  if (all(key %in% rownames(V))) {
    v <- row_quad(G, as.matrix(V[key, key, drop = FALSE]), G)
    v[!is.na(v) & v < 0] <- 0
  }
  se <- sqrt(v)
  z <- stats::qnorm((1 + level) / 2)
  em <- linkfunctions7::linkfun(lk, fitv)
  se_eta <- se / abs(linkfunctions7::dlinkinv(lk, em))
  ends <- cbind(linkfunctions7::linkinv(lk, em - z * se_eta),
                linkfunctions7::linkinv(lk, em + z * se_eta))
  data.frame(fit = fitv, se = se, lower = pmin(ends[, 1L], ends[, 2L]),
             upper = pmax(ends[, 1L], ends[, 2L]))
}
