#' A Block of the Jeffreys Penalty
#'
#' @description
#' Describes one block of coefficients on which the Jeffreys penalty
#' \eqn{a\log|\tilde I|} is computed by [jfp_state()] and its derivatives.
#'
#' @details
#' The block's rows \eqn{\tilde V_{q,t}} are the derivative of the predictor
#' of equation \eqn{q} at observation \eqn{t} with respect to the block's
#' coefficients, preceded by the constant direction of its formula level when
#' `const` is `TRUE`. For a block of `kind = "static"` the rows are the
#' block's columns of the equation's design, which do not depend on the
#' coefficients, and the constant direction is the column of ones.
#'
#' The other kinds describe a block inside the subformula of a term:
#' `"refresh"` for a term whose block is recomputed at the coefficients
#' (`nl()`, a smoothed `seg()`), and `"filter"` for a parameter of `gas()`.
#' Their rows move with the coefficients, and [jfp_state()] signals an error
#' for them: they are reserved for a later phase, and the formulas of the
#' engine do not change when they are added.
#'
#' @param param The name of the distribution parameter whose predictor the
#'   block enters.
#' @param cols Integer positions of the block's columns in that equation's
#'   design.
#' @param power The power \eqn{a > 0} of the determinant.
#' @param const Whether the constant direction is part of the block, as
#'   decided by [jfp_constant()].
#' @param kind One of `"static"`, `"refresh"` or `"filter"`.
#' @param term,slot For the kinds `"refresh"` and `"filter"`, the term and the
#'   parameter of the term whose subformula holds the block.
#'
#' @return A list of class `jfp_block` with the arguments as fields.
#'
#' @seealso [jfp_state()], [jfp_constant()]
#'
#' @keywords internal
jfp_block <- function(param, cols, power = 0.5, const = FALSE,
                      kind = c("static", "refresh", "filter"), term = NULL,
                      slot = NULL) {
  kind <- match.arg(kind)
  if (!is.character(param) || length(param) != 1L) {
    stop("'param' must be a single parameter name.", call. = FALSE)
  }
  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) ||
      power <= 0) {
    stop("'power' must be a single positive number.", call. = FALSE)
  }
  structure(list(param = param, cols = as.integer(cols), power = power,
                 const = isTRUE(const), kind = kind, term = term,
                 slot = slot),
            class = "jfp_block")
}


#' Whether the Constant Direction Enters a Jeffreys Block
#'
#' @description
#' Applies the rule of the `intercept` argument of `jfp()` to the columns of
#' a formula level and to the columns of the block.
#'
#' @details
#' With `intercept = FALSE` the direction is never included. Where the
#' column of ones lies in the span of the block's own columns the direction
#' is not included either, whatever `intercept` says, because the determinant
#' would be zero and the penalty is already invariant to centering a column.
#' Otherwise `TRUE` includes it, and `NULL` includes it when the column of
#' ones lies in the span of `X_level`, that is when the model is invariant to
#' centering a column of the block.
#'
#' The span is tested by [ones_in_span()].
#'
#' @param intercept `NULL`, `TRUE` or `FALSE`.
#' @param X_level The columns of the formula level outside the block, or
#'   `NULL`.
#' @param X_A The block's columns.
#' @param tol The tolerance of the span test.
#'
#' @return A single logical.
#'
#' @seealso [jfp_block()]
#'
#' @keywords internal
jfp_constant <- function(intercept, X_level, X_A, tol = 1e-7) {
  if (!is.null(intercept) &&
      (!is.logical(intercept) || length(intercept) != 1L || is.na(intercept))) {
    stop("'intercept' must be NULL, TRUE or FALSE.", call. = FALSE)
  }
  if (isFALSE(intercept)) return(FALSE)
  if (ones_in_span(X_A, tol)) return(FALSE)
  if (isTRUE(intercept)) return(TRUE)
  ones_in_span(X_level, tol)
}


#' Whether the Column of Ones Lies in the Span of a Matrix
#'
#' @description
#' Compares the rank of the matrix with the rank of the matrix with a column
#' of ones appended, both from the pivoted QR decomposition with the tolerance
#' with which the toolkit detects aliased coefficients.
#'
#' @param X A matrix, a \pkg{Matrix} object, or `NULL`.
#' @param tol The tolerance of [qr()].
#'
#' @return A single logical.
#'
#' @keywords internal
ones_in_span <- function(X, tol = 1e-7) {
  if (is.null(X) || !ncol(X)) return(FALSE)
  X <- as_dense(X)
  qr(cbind(X, 1), tol = tol)$rank == qr(X, tol = tol)$rank
}


#' The State of the Jeffreys Penalty at One Point
#'
#' @description
#' Computes, at given coefficients, everything that the value, the gradient
#' and the Hessian of the Jeffreys penalty \eqn{\sum_b a_b\log|\tilde I_b|}
#' read: for each block its rows, its expected information, the Cholesky
#' factor of that information, the whitened rows and the leverages.
#'
#' @details
#' For a block \eqn{b} with rows \eqn{\tilde V_{q,t}},
#' \deqn{\tilde I_b = \sum_t w_t\sum_{q,r}\mathcal I_{qr}(\eta_t)\,
#'   \tilde V_{q,t}^\top\tilde V_{r,t},}
#' where \eqn{\mathcal I_{qr}} is the family's expected information on the
#' link scale, minus [distributions7::distrib_expected_hessian()]. With
#' \eqn{\tilde I_b = U^\top U}, the whitened rows are
#' \eqn{Z_q = \tilde V_q U^{-1}} and the leverages are
#' \eqn{h_{qr,t} = Z_{q,t}Z_{r,t}^\top}.
#'
#' The point is the joint vector \eqn{u} of [statmod_full_information()]. For
#' a model without a structural term it is the stacked coefficients. For a
#' model carrying a filter it is the coefficients followed by the filter's
#' free parameters, and the derivative of every predictor in \eqn{u} is read
#' from the forward recursion ([jfp_joint_rows()]); a block in an equation
#' that the filter drives moves with the coefficients and is rejected in
#' this phase. A model carrying a mixture over latent states is always
#' rejected, because its expected information does not exist.
#'
#' The derivatives of the expected information are computed at the first
#' request and kept in the state, so that the value, the gradient and the
#' Hessian at one point call the family once per order.
#'
#' Where a block's information is not positive definite its entry carries
#' `ok = FALSE` and a log-determinant of `-Inf`, and the gradient and the
#' Hessian signal an error.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef A named list of coefficient vectors.
#' @param blocks A list of [jfp_block()] objects.
#' @param emethod The `method` of [distributions7::distrib_expected_hessian()],
#'   or `NULL` for the family's own closed form.
#'
#' @return An environment with the fields `spec`, `design` (refreshed at
#'   `coef`), `coef`, `ev`, `params`, `npar`, `offs`, `nb` (the number of
#'   coefficients), `total` (the length of \eqn{u}), `joint` (`NULL`, or what
#'   [jfp_joint_rows()] returns), `E`, `emethod` and `blocks`, each entry of
#'   the last carrying `block`,
#'   `index` (positions in \eqn{u}), `rows`, `info`, `U`, `Z`, `G`,
#'   `logdet` and `ok`.
#'
#' @seealso [jfp_value()], [jfp_gradient()], [jfp_hessian()]
#'
#' @keywords internal
jfp_state <- function(spec, design, coef, blocks, emethod = NULL) {
  if (inherits(blocks, "jfp_block")) blocks <- list(blocks)
  params <- spec@distrib@params
  ev <- statmod_eta(spec, design, coef)
  if (length(ev$regimes)) {
    stop("The Jeffreys penalty needs the expected information, which a ",
         "model carrying a\n  mixture over latent states does not have.",
         call. = FALSE)
  }
  design0 <- design
  design <- statmod_design_at(spec, coef, design)
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  st <- new.env(parent = emptyenv())
  st$spec <- spec
  st$design <- design
  st$coef <- coef
  st$ev <- ev
  st$params <- params
  st$npar <- npar
  st$offs <- offs
  st$nb <- sum(npar)
  st$total <- sum(npar)
  st$joint <- NULL
  st$emethod <- emethod
  st$E <- statmod_family_hessian(spec, ev$theta, TRUE, emethod)
  if (length(ev$filters)) {
    st$joint <- jfp_joint_rows(spec, design0, coef, ev)
    st$total <- length(st$joint$keep)
    driven <- st$joint$driven
    for (b in blocks) {
      if (b$param %in% driven) {
        stop(sprintf(paste0("A Jeffreys block in '%s', an equation that a ",
                            "filter drives, moves with\n  the coefficients, ",
                            "which is not implemented yet."), b$param),
             call. = FALSE)
      }
    }
  }
  st$blocks <- lapply(blocks, function(b) jfp_block_state(st, b))
  st
}


#' The Derivative of Every Predictor Beside a Filter
#'
#' @description
#' Returns, for a model carrying a filter, the derivative of each equation's
#' predictor in the joint vector of the coefficients and the filter's free
#' parameters, with what the second derivative of a filtered predictor needs.
#'
#' @details
#' The rows are those of [joint_design_rows()], with the equation that the
#' filter drives replaced by the forward Jacobian of the recursion
#' ([filter_curvature()], or [joint_filter_curvature()] for a joint filter),
#' and a level that an intercept holds dropped as
#' [statmod_full_information()] drops it.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design as it arrived, carrying the structural state.
#' @param coef The coefficients.
#' @param ev What [statmod_eta()] returns at `coef`.
#'
#' @return A list with `V` (one dense matrix per distribution parameter over
#'   the kept coordinates), `keep`, `driven` (the equations whose predictor
#'   the filter writes), and `curv`, a function of the per-observation weights
#'   (a named list keyed by parameter) returning
#'   \eqn{\sum_t\sum_c g_{c,t}\,\partial^2\eta_{c,t}/\partial u\,\partial u}
#'   over the kept coordinates.
#'
#' @keywords internal
jfp_joint_rows <- function(spec, design, coef, ev) {
  jd <- joint_design_rows(spec, design, coef)
  params <- spec@distrib@params
  f <- jd$f
  Vs <- jd$V
  gl <- distributions7::distrib_gradient(spec@distrib, spec@response,
                                         ev$theta, scale = "link",
                                         threads = spec@threads)
  H <- distributions7::distrib_hessian(spec@distrib, spec@response, ev$theta,
                                       scale = "link", threads = spec@threads)
  D3 <- distributions7::distrib_deriv3(spec@distrib, spec@response, ev$theta,
                                       scale = "link", threads = spec@threads)
  joint <- isTRUE(f$joint)
  driven <- if (joint) f$params else f$param
  cv <- if (joint) {
    joint_filter_curvature(spec, design, f, Vs, gl, H, D3, ev$theta)
  } else {
    filter_curvature(spec, design, f, match(f$param, params), Vs, gl, H, D3,
                     ev$theta)
  }
  V <- Vs
  if (joint) {
    for (k in seq_along(f$params)) {
      V[[match(f$params[k], params)]] <- cv$jacobian[[k]]
    }
  } else {
    V[[match(f$param, params)]] <- cv$jacobian
  }
  keep <- jd$keep
  V <- lapply(V, function(M) M[, keep, drop = FALSE])
  w <- spec@weights
  n <- spec@n_obs
  curv <- function(g) {
    out <- if (joint) {
      gg <- stats::setNames(lapply(f$params, function(p) w * g[[p]]),
                            f$params)
      joint_filter_curvature(spec, design, f, Vs, gl, H, D3, ev$theta,
                             g = gg)$curvature
    } else {
      filter_curvature(spec, design, f, match(f$param, params), Vs, gl, H,
                       D3, ev$theta, g = w * rep_len(g[[f$param]], n))$curvature
    }
    as_dense(out)[keep, keep, drop = FALSE]
  }
  list(V = V, keep = keep, driven = driven, curv = curv, nb = jd$nb)
}


#' The Rows, Information and Factor of One Jeffreys Block
#'
#' @param st The state under construction, from [jfp_state()].
#' @param b A [jfp_block()].
#'
#' @return The block's entry of the state.
#'
#' @keywords internal
jfp_block_state <- function(st, b) {
  if (!inherits(b, "jfp_block")) {
    stop("Every block must be built by jfp_block().", call. = FALSE)
  }
  rows <- jfp_rows(st, b)
  params <- st$params
  spec <- st$spec
  n <- spec@n_obs
  w <- spec@weights
  live <- which(!vapply(rows, is.null, logical(1)))
  p <- ncol(rows[[live[1L]]])
  info <- matrix(0, p, p)
  for (q in live) {
    for (r in live) {
      wv <- -w * rep_len(st$E[[hess_key(params, q, r)]], n)
      info <- info + wcrossprod(rows[[q]], wv, rows[[r]], spec@threads)
    }
  }
  info <- (info + t(info)) / 2
  U <- if (all(is.finite(info))) {
    tryCatch(chol(info), error = function(e) NULL)
  } else NULL
  out <- list(block = b, index = st$offs[match(b$param, params)] + b$cols,
              rows = rows, info = info, U = U, ok = !is.null(U))
  if (!out$ok) {
    out$logdet <- -Inf
    return(out)
  }
  Ui <- backsolve(U, diag(p))
  Z <- lapply(rows, function(R) if (is.null(R)) NULL else R %*% Ui)
  # the leverages h_qr = Z_q Z_r', symmetric in (q, r) and stored both ways
  # round, which is the layout u_vector() reads its G in
  G <- vector("list", length(params))
  for (q in seq_along(params)) G[[q]] <- vector("list", length(params))
  for (q in live) {
    for (r in live) {
      if (r < q) next
      G[[q]][[r]] <- rowSums(Z[[q]] * Z[[r]])
      if (r != q) G[[r]][[q]] <- G[[q]][[r]]
    }
  }
  out$Z <- Z
  out$G <- G
  out$logdet <- 2 * sum(log(diag(U)))
  out
}


#' The Rows of a Jeffreys Block
#'
#' @description
#' Returns \eqn{\tilde V_q}, one matrix per distribution parameter, `NULL`
#' where the block does not enter that parameter's predictor.
#'
#' @details
#' For a block of kind `"static"` the matrix of its equation is the block's
#' columns of the design, preceded by a column of ones when the constant
#' direction is included. A sparse design is densified over the block's
#' columns only.
#'
#' @param st The state under construction, from [jfp_state()].
#' @param b A [jfp_block()].
#'
#' @return A list as long as the distribution's parameters.
#'
#' @keywords internal
jfp_rows <- function(st, b) {
  if (b$kind != "static") {
    stop(sprintf(paste0("A Jeffreys block of kind '%s' moves with the ",
                        "coefficients, which is not\n  implemented yet."),
                 b$kind), call. = FALSE)
  }
  params <- st$params
  a <- match(b$param, params)
  if (is.na(a)) {
    stop(sprintf("The distribution has no parameter '%s'.", b$param),
         call. = FALSE)
  }
  if (!length(b$cols) || any(b$cols < 1L | b$cols > st$npar[a])) {
    stop(sprintf("'cols' must be positions among the %d columns of '%s'.",
                 st$npar[a], b$param), call. = FALSE)
  }
  X <- as_dense(st$design[[b$param]]$X[, b$cols, drop = FALSE])
  if (b$const) X <- cbind(1, X)
  dimnames(X) <- NULL
  rows <- vector("list", length(params))
  rows[[a]] <- X
  rows
}


#' The Value of the Jeffreys Penalty
#'
#' @param st A state, from [jfp_state()].
#'
#' @return \eqn{\sum_b a_b\log|\tilde I_b|}, a single number, `-Inf` where a
#'   block's information is not positive definite.
#'
#' @seealso [jfp_gradient()], [jfp_hessian()]
#'
#' @keywords internal
jfp_value <- function(st) {
  sum(vapply(st$blocks, function(bs) bs$block$power * bs$logdet, numeric(1)))
}


#' The Gradient of the Jeffreys Penalty in the Joint Vector
#'
#' @description
#' Computes \eqn{\partial_k\sum_b a_b\log|\tilde I_b|} over every coordinate
#' of \eqn{u}.
#'
#' @details
#' For a block whose rows do not move,
#' \deqn{\partial_k\log|\tilde I| = \mathrm{tr}(\tilde I^{-1}\partial_k\tilde I)
#'   = \sum_t w_t\sum_c V_{c,t,k}\,s_{c,t},\qquad
#'   s_{c,t} = -\sum_{q,r}h_{qr,t}\,\partial_{\eta_c}E_{qr,t},}
#' with \eqn{V_{c,t,k}} the derivative of the predictor of equation \eqn{c} in
#' \eqn{u_k} and \eqn{\partial_{\eta_c}E_{qr}} from
#' [distributions7::distrib_dexpected_hessian()]. This is the contraction of
#' [u_vector()] with the block's leverages in place of the diagonal of a
#' matrix over \eqn{u}, which cannot hold the constant direction.
#'
#' @param st A state, from [jfp_state()].
#'
#' @return A numeric vector over \eqn{u}.
#'
#' @seealso [jfp_value()], [jfp_hessian()]
#'
#' @keywords internal
jfp_gradient <- function(st) {
  jfp_check_ok(st)
  dE <- jfp_dexpected(st)
  params <- st$params
  key <- function(a, b, k) distributions7::dexpected_key(params, a, b, k)
  out <- numeric(st$total)
  for (bs in st$blocks) {
    out <- out + bs$block$power *
      u_vector(st$spec, st$design, st$coef, NULL, params, st$npar, st$offs,
               st$total, d3 = dE, G = bs$G, key = key,
               rows = if (is.null(st$joint)) NULL else st$joint$V)
  }
  out
}


#' The Hessian of the Jeffreys Penalty in the Joint Vector
#'
#' @description
#' Computes \eqn{\partial_{kl}\sum_b a_b\log|\tilde I_b|} as a matrix with the
#' sparsity of the likelihood's Hessian minus a product of low rank.
#'
#' @details
#' For a block whose rows do not move,
#' \deqn{\partial_{kl}\log|\tilde I| = \mathrm{tr}(\tilde I^{-1}
#'   \partial_{kl}\tilde I) - \mathrm{tr}(\tilde I^{-1}\partial_k\tilde I\,
#'   \tilde I^{-1}\partial_l\tilde I).}
#' The first term is
#' \eqn{\sum_{c,d}X_c^\top\mathrm{diag}(w\,W_{cd})X_d} with
#' \eqn{W_{cd} = -\sum_{q,r}h_{qr}\,\partial^2_{\eta_c\eta_d}E_{qr}}, from
#' [distributions7::distrib_d2expected_hessian()], assembled by
#' [statmod_information_at()]. Where an equation carries a term whose block is
#' recomputed at the coefficients, the second derivative of its predictor
#' adds \eqn{\sum_t w_t\sum_c s_{c,t}\,\partial^2\eta_{c,t}/\partial u_k
#' \partial u_l}, with \eqn{s_c} the weights of [jfp_gradient()], read through
#' [mode_curvature()].
#'
#' The second term is \eqn{BB^\top}: row \eqn{k} of \eqn{B} holds the upper
#' triangle of \eqn{U^{-\top}\partial_k\tilde I\,U^{-1}}, the off-diagonal
#' entries multiplied by \eqn{\sqrt2}, so that its rank is at most
#' \eqn{\tilde p(\tilde p+1)/2}.
#'
#' @param st A state, from [jfp_state()].
#' @param assemble Whether to return the sum as a dense matrix as well.
#'
#' @return A list with `sparse` (the first term times the powers, in the
#'   storage type of the information), `lowrank` (\eqn{B} times the square
#'   roots of the powers, one block of columns per Jeffreys block),
#'   `weights` (the per-observation weights of the first term, keyed as
#'   [statmod_information_at()] reads its `H`, so that `sparse` is that
#'   function's result at `H = weights` plus `refresh`), `refresh` (the part
#'   that a refreshed block contributes, or `NULL`) and, when `assemble` is
#'   `TRUE`, `hessian`, equal to `sparse - tcrossprod(lowrank)` as a base
#'   matrix.
#'
#' @seealso [jfp_value()], [jfp_gradient()]
#'
#' @keywords internal
jfp_hessian <- function(st, assemble = TRUE) {
  jfp_check_ok(st)
  dE <- jfp_dexpected(st)
  d2E <- jfp_d2expected(st)
  spec <- st$spec
  params <- st$params
  np <- length(params)
  n <- spec@n_obs
  w <- spec@weights
  npar <- st$npar
  offs <- st$offs
  rf <- attr(st$design, "refresh")
  weights <- list()
  refresh <- NULL
  lowrank <- vector("list", length(st$blocks))
  for (i in seq_along(st$blocks)) {
    bs <- st$blocks[[i]]
    a <- bs$block$power
    live <- which(!vapply(bs$rows, is.null, logical(1)))
    # the weights of the first term, one per pair of equations (c <= d), keyed
    # as statmod_information_at() reads its H: it forms -X_c' diag(w H) X_d,
    # so H_cd = sum_qr h_qr d2E_qr,cd gives the trace of the second derivative
    for (cc in seq_len(np)) {
      for (dd in cc:np) {
        v <- numeric(n)
        for (q in live) {
          for (r in live) {
            v <- v + bs$G[[q]][[r]] * rep_len(d2E[[
              distributions7::d2expected_key(params, q, r, cc, dd)]], n)
          }
        }
        k <- hess_key(params, cc, dd)
        weights[[k]] <- if (is.null(weights[[k]])) a * v else
          weights[[k]] + a * v
      }
    }
    if ((!is.null(rf) && length(rf)) || !is.null(st$joint)) {
      s <- jfp_first_weights(st, bs, dE)
      D <- matrix(0, st$total, st$total)
      if (!is.null(rf) && length(rf)) {
        nb <- st$nb
        D[seq_len(nb), seq_len(nb)] <- -a * mode_curvature(
          spec, st$design, st$coef, params, npar, offs, nb, gl = s)
      }
      if (!is.null(st$joint)) D <- D + a * st$joint$curv(s)
      refresh <- if (is.null(refresh)) D else refresh + D
    }

    # B: one column per entry i <= j of the whitened derivative
    p <- ncol(bs$Z[[live[1L]]])
    ij <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
    gam <- ifelse(ij[, 1L] == ij[, 2L], 1, sqrt(2))
    B <- matrix(0, st$total, nrow(ij))
    for (cc in seq_len(np)) {
      if (npar[cc] == 0L && is.null(st$joint)) next
      Fc <- matrix(0, n, nrow(ij))
      for (q in live) {
        for (r in live) {
          dv <- -rep_len(dE[[distributions7::dexpected_key(params, q, r,
                                                            cc)]], n)
          Fc <- Fc + dv * bs$Z[[q]][, ij[, 1L], drop = FALSE] *
            bs$Z[[r]][, ij[, 2L], drop = FALSE]
        }
      }
      Fc <- sweep(Fc, 2L, gam, `*`)
      if (is.null(st$joint)) {
        B[offs[cc] + seq_len(npar[cc]), ] <-
          as.matrix(crossprod(st$design[[params[cc]]]$X, w * Fc))
      } else {
        B <- B + crossprod(st$joint$V[[cc]], w * Fc)
      }
    }
    lowrank[[i]] <- sqrt(a) * B
  }
  lowrank <- do.call(cbind, lowrank)
  sparse <- if (is.null(st$joint)) {
    statmod_information_at(spec, st$coef, st$design, expected = TRUE,
                           emethod = st$emethod, H = weights)
  } else {
    # the same first term on the joint rows: -sum V_c' diag(w H_cd) V_d over
    # the ordered pairs, each unordered pair's weights read once
    out <- matrix(0, st$total, st$total)
    for (cc in seq_len(np)) {
      for (dd in seq_len(np)) {
        wv <- w * rep_len(weights[[hess_key(params, cc, dd)]], n)
        out <- out - wcrossprod(st$joint$V[[cc]], wv, st$joint$V[[dd]],
                                spec@threads)
      }
    }
    (out + t(out)) / 2
  }
  if (!is.null(refresh)) sparse <- sparse + refresh
  out <- list(sparse = sparse, lowrank = lowrank, weights = weights,
              refresh = refresh)
  if (assemble) {
    H <- as_dense(sparse) - tcrossprod(lowrank)
    out$hessian <- (H + t(H)) / 2
  }
  out
}


#' The Weights of the First Derivative of a Jeffreys Block
#'
#' @description
#' Returns \eqn{s_{c,t} = -\sum_{q,r}h_{qr,t}\,\partial_{\eta_c}E_{qr,t}} for
#' every equation \eqn{c}, the per-observation weights that [jfp_gradient()]
#' contracts against each equation's design and that [jfp_hessian()] contracts
#' against the second derivative of a refreshed block.
#'
#' @param st A state, from [jfp_state()].
#' @param bs A block's entry of the state.
#' @param dE The derivative of the expected information, from
#'   [jfp_dexpected()].
#'
#' @return A named list of numeric vectors, one per distribution parameter.
#'
#' @keywords internal
jfp_first_weights <- function(st, bs, dE) {
  params <- st$params
  n <- st$spec@n_obs
  live <- which(!vapply(bs$rows, is.null, logical(1)))
  out <- stats::setNames(vector("list", length(params)), params)
  for (cc in seq_along(params)) {
    v <- numeric(n)
    for (q in live) {
      for (r in live) {
        v <- v - bs$G[[q]][[r]] * rep_len(dE[[
          distributions7::dexpected_key(params, q, r, cc)]], n)
      }
    }
    out[[cc]] <- v
  }
  out
}


#' The Derivatives of the Expected Information Kept in a Jeffreys State
#'
#' @description
#' `jfp_dexpected()` returns [distributions7::distrib_dexpected_hessian()] and
#' `jfp_d2expected()` returns [distributions7::distrib_d2expected_hessian()],
#' both on the link scale at the state's predictors, computed at the first
#' request and kept in the state.
#'
#' @param st A state, from [jfp_state()].
#'
#' @return A named list of per-observation components.
#'
#' @keywords internal
jfp_dexpected <- function(st) {
  if (is.null(st$dE)) {
    st$dE <- distributions7::distrib_dexpected_hessian(
      st$spec@distrib, st$spec@response, st$ev$theta, scale = "link",
      method = st$emethod)
  }
  st$dE
}

#' @rdname jfp_dexpected
#' @keywords internal
jfp_d2expected <- function(st) {
  if (is.null(st$d2E)) {
    st$d2E <- distributions7::distrib_d2expected_hessian(
      st$spec@distrib, st$spec@response, st$ev$theta, scale = "link",
      method = st$emethod, threads = st$spec@threads)
  }
  st$d2E
}


#' Signals an Error Where a Jeffreys Block Has No Usable Information
#'
#' @param st A state, from [jfp_state()].
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
jfp_check_ok <- function(st) {
  bad <- which(!vapply(st$blocks, function(bs) bs$ok, logical(1)))
  if (length(bad)) {
    stop(sprintf(paste0("The expected information of Jeffreys block %d is ",
                        "not positive definite\n  at these coefficients."),
                 bad[1L]), call. = FALSE)
  }
  invisible(NULL)
}
