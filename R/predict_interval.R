#' @include predict_random.R
NULL

#' The Within-Group Rows of a Random-Effect Term
#'
#' @description
#' The row of the within-group design each observation carries, at the new
#' data where the specification has some and at the fitting rows otherwise.
#'
#' @param spec A [StatmodSpec()].
#' @param param The distribution parameter.
#' @param key The term's key.
#'
#' @return A numeric matrix with one row per observation and one column per
#'   within-group coordinate.
#'
#' @keywords internal
random_within <- function(spec, param, key) {
  tm <- spec@terms[[param]][[key]]
  if (!is.null(spec@newdata)) {
    return(as.matrix(modelterms7::term_within(tm, spec@newdata)))
  }
  # at the fitting rows the block itself says it: a row is non-zero in its
  # own group's d columns alone, so summing the groups' columns returns the
  # within-group row
  gr <- modelterms7::term_group(tm)
  as.matrix(modelterms7::term_matrix(tm) %*%
              kronecker(rep(1, length(gr$levels)), diag(gr$dim)))
}


#' The Gaussian Prior of the Effects a Prediction Sets Aside
#'
#' @description
#' The covariance of one new group's effects, block by block: a term with a
#' prior of its own is one block, and the terms a label ties together are one
#' block over all of them.
#'
#' @details
#' A new group's effect is drawn from the prior the fit estimated, so its
#' variance is that prior's, at the hyperparameters the fit reached; the
#' uncertainty of those hyperparameters is not propagated, which the page of
#' [predict.StatmodFit()] states. A prior that is not Gaussian has no
#' covariance to read here and is rejected by name.
#'
#' Where a label ties terms together, the prior is one multivariate Gaussian
#' over the coordinates of all of them, and every member has to be set aside:
#' a member read with its own estimated effect would condition the others on
#' it, which is a different question. A member written inside another term's
#' subformula is rejected as well, its coordinates not being a row of a
#' random-effect design.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#' @param fit The [StatmodFit()].
#' @param aside The rows of [random_modes()] that are not `"conditional"`.
#'
#' @return A list of blocks, each a list with `members` (a data frame of
#'   `param`, `key`, `dim`) and `Sigma`, the covariance over the members'
#'   coordinates in that order.
#'
#' @keywords internal
random_blocks <- function(spec, design, fit, aside) {
  out <- list()
  done <- character(0)
  units <- statmod_penalized(spec, design)
  for (i in seq_len(nrow(aside))) {
    p <- aside$param[i]
    k <- aside$key[i]
    if (paste(p, k) %in% done) next
    tm <- spec@terms[[p]][[k]]
    tag <- modelterms7::term_tag(tm)
    if (is.na(tag)) {
      pr <- random_prior(spec, design, fit, p, k)
      if (!isTRUE(pr$gaussian)) {
        stop(sprintf(paste0("The prior of '%s' is not Gaussian, so a new ",
                            "group's interval has no\n  covariance to read. ",
                            "interval = \"confidence\" is available."), k),
             call. = FALSE)
      }
      out[[length(out) + 1L]] <- list(
        members = data.frame(param = p, key = k, dim = pr$dim,
                             stringsAsFactors = FALSE),
        Sigma = pr$chol %*% t(pr$chol))
      done <- c(done, paste(p, k))
      next
    }
    cu <- Filter(function(u) !is.null(u$class) && any(vapply(
      u$class$pieces, function(pc) identical(pc$param, p) &&
        identical(pc$term, k), logical(1))), units)
    if (length(cu) != 1L) {
      stop(sprintf("'%s' has no single covariance class to read.", k),
           call. = FALSE)
    }
    u <- cu[[1L]]
    cl <- u$class
    mem <- do.call(rbind, lapply(cl$pieces, function(pc) {
      if (!is.null(pc$within)) {
        stop(sprintf(paste0("'%s' shares a covariance with an effect written ",
                            "inside a subformula,\n  whose coordinates a new ",
                            "group's interval cannot read."), k),
             call. = FALSE)
      }
      data.frame(param = pc$param, key = pc$term, dim = as.integer(pc$dim),
                 stringsAsFactors = FALSE)
    }))
    miss <- !paste(mem$param, mem$key) %in% paste(aside$param, aside$key)
    if (any(miss)) {
      stop(sprintf(paste0("'%s' shares a covariance with '%s', which is read ",
                          "with its own\n  estimated effects. Set both aside, ",
                          "or neither."), k, mem$key[miss][1L]),
           call. = FALSE)
    }
    th <- as.list(fit@hyper[[u$param]][[u$key]])
    pen <- u$penalty
    if (!isTRUE(penalties7::beta_quadratic(pen, th))) {
      stop(sprintf(paste0("The covariance shared by '%s' is not Gaussian, so a ",
                          "new group's\n  interval has no covariance to read."),
                   k), call. = FALSE)
    }
    D <- sum(mem$dim)
    H <- as.matrix(penalties7::penalty_hessian(pen, rep(0, pen@n_coef), th))
    out[[length(out) + 1L]] <- list(members = mem,
                                    Sigma = solve(H[seq_len(D), seq_len(D),
                                                    drop = FALSE]))
    done <- c(done, paste(mem$param, mem$key))
  }
  out
}


#' The Covariance of Every Predictor at One Row, for a New Group
#'
#' @description
#' For each observation, the covariance over the distribution's parameters of
#' their linear predictors: the estimation uncertainty of the fixed part, from
#' the variance matrix, plus the variance a new group's effects add.
#'
#' @details
#' The fixed part is \eqn{x_{ip}^\top V_{pq} x_{iq}}, read with the columns
#' of the effects set aside at zero, so it is the uncertainty of the typical
#' group's predictor; the part the effects add is
#' \eqn{z_{ip}^\top \Sigma_{pq} z_{iq}}, summed over the blocks of
#' [random_blocks()].
#'
#' @param object The [StatmodFit()].
#' @param spec,design The specification and design at the rows predicted.
#' @param blocks The blocks of [random_blocks()].
#' @param fixed `FALSE` to skip the fixed part, which is then zero.
#' @param ... Passed to [vcov.StatmodFit()].
#'
#' @return A list with `fixed` and `random`, arrays of dimension
#'   `P x P x n` over the parameters.
#'
#' @keywords internal
predictive_cov <- function(object, spec, design, blocks, fixed = TRUE, ...) {
  params <- spec@distrib@params
  P <- length(params)
  n <- spec@n_obs
  fx <- array(0, c(P, P, n), list(params, params, NULL))
  V <- if (fixed) vcov(object, readable = FALSE, ...) else NULL
  for (a in if (fixed) seq_len(P) else integer(0)) {
    for (b in a:P) {
      da <- design[[params[a]]]
      db <- design[[params[b]]]
      if (!da$npar || !db$npar) next
      ka <- paste(params[a], da$coef_names, sep = ":")
      kb <- paste(params[b], db$coef_names, sep = ":")
      v <- rep(NA_real_, n)
      if (all(c(ka, kb) %in% rownames(V))) {
        Vab <- as.matrix(V[ka, kb, drop = FALSE])
        if (all(is.finite(Vab))) {
          v <- rowSums((as.matrix(da$X) %*% Vab) * as.matrix(db$X))
        }
      }
      fx[a, b, ] <- v
      fx[b, a, ] <- v
    }
  }
  random <- array(0, c(P, P, n), list(params, params, NULL))
  for (bl in blocks) {
    mem <- bl$members
    Zs <- lapply(seq_len(nrow(mem)), function(j)
      random_within(spec, mem$param[j], mem$key[j]))
    at <- cumsum(c(0L, mem$dim))
    for (j in seq_len(nrow(mem))) {
      for (l in seq_len(nrow(mem))) {
        S <- bl$Sigma[at[j] + seq_len(mem$dim[j]), at[l] + seq_len(mem$dim[l]),
                      drop = FALSE]
        a <- match(mem$param[j], params)
        b <- match(mem$param[l], params)
        random[a, b, ] <- random[a, b, ] +
          rowSums((Zs[[j]] %*% S) * Zs[[l]])
      }
    }
  }
  list(fixed = fx, random = random)
}


#' The Predictive Distribution of the Response
#'
#' @description
#' The median, the ends of an interval and the standard deviation of the
#' response at each row, under a mixture of the family over the predictors.
#'
#' @details
#' Each component \eqn{c} gives, at row \eqn{i}, the predictors' mean
#' \eqn{\eta_{ci}} and their covariance \eqn{C_{ci}}, the predictors being
#' taken jointly Gaussian within it. The distribution of \eqn{Y} is the
#' equally weighted mixture of the family over the components and, within
#' each, over the predictors,
#' \deqn{G_i(y) = \frac{1}{B}\sum_{c=1}^{B} E_{\eta \sim
#'   \mathrm{N}(\eta_{ci}, C_{ci})}[F(y \mid h^{-1}(\eta))],}
#' each expectation evaluated on a Gauss-Hermite product grid of at most 400
#' nodes, or at the single node \eqn{\eta_{ci}} where \eqn{C_{ci}} is zero
#' at every row. One component carries the default interval, two
#' sources of variation being folded into its covariance; the parametric
#' bootstrap gives one component per replica. A quantile
#' of \eqn{G_i} lies between the smallest and the largest quantile of the
#' conditional laws at the nodes, a mixture's distribution function being an
#' average of theirs, and is found there by bisection -- on the integers for
#' a discrete family, where it is the smallest value at which \eqn{G_i}
#' reaches the level. Only the family's distribution function and quantile
#' are read, so no location parameter is needed. The standard deviation is
#' \eqn{\sqrt{E[\mathrm{Var}(Y\mid\eta)] + \mathrm{Var}(E[Y\mid\eta])}},
#' `NA` where the family's mean or variance does not exist.
#'
#' @param spec The specification at the rows predicted.
#' @param comps A list of components, each a list with `eta`, a named list of
#'   the predictors' means, and `C`, their covariance as an array
#'   `P x P x n`.
#' @param level The interval's level.
#'
#' @return A data frame with `fit` (the median), `se`, `lower` and `upper`.
#'
#' @keywords internal
predictive_response <- function(spec, comps, level) {
  d <- spec@distrib
  params <- d@params
  links <- d@link_params
  P <- length(params)
  n <- spec@n_obs
  kk <- max(2L, min(20L, floor(400^(1 / P))))
  q <- gauss_hermite(kk)
  grid <- as.matrix(expand.grid(rep(list(seq_len(kk)), P)))
  x_gh <- matrix(q$x[grid], nrow(grid), P) * sqrt(2)
  w_gh <- apply(matrix(q$w[grid], nrow(grid), P), 1L, prod) / pi^(P / 2)
  bad <- rep(FALSE, n)
  pieces <- vector("list", length(comps))
  w <- numeric(0)
  for (cc in seq_along(comps)) {
    C <- comps[[cc]]$C
    eta0 <- comps[[cc]]$eta
    # a component with no covariance at any row is its mean alone, so one
    # node carries it rather than a grid of coincident ones
    point <- all(is.finite(C)) && !any(C != 0)
    x <- if (point) matrix(0, 1L, P) else x_gh
    wc <- if (point) 1 else w_gh
    K <- nrow(x)
    M <- matrix(vapply(params, function(p) rep_len(as.numeric(eta0[[p]]), n),
                       numeric(n)), n, P)
    e <- array(NA_real_, c(n, K, P))
    for (i in seq_len(n)) {
      Ci <- C[, , i]
      if (!all(is.finite(Ci)) || !all(is.finite(M[i, ]))) {
        bad[i] <- TRUE
        next
      }
      if (point) {
        e[i, 1L, ] <- M[i, ]
        next
      }
      ev <- eigen((Ci + t(Ci)) / 2, symmetric = TRUE)
      L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), P)
      e[i, , ] <- sweep(x %*% t(L), 2L, M[i, ], "+")
    }
    pieces[[cc]] <- e
    w <- c(w, wc / length(comps))
  }
  K <- length(w)
  eta <- array(NA_real_, c(n, K, P))
  at <- 0L
  for (e in pieces) {
    eta[, at + seq_len(dim(e)[2L]), ] <- e
    at <- at + dim(e)[2L]
  }
  theta <- stats::setNames(lapply(seq_len(P), function(j)
    as.numeric(linkfunctions7::linkinv(links[[params[j]]],
                                       as.vector(eta[, , j])))), params)
  ok <- !bad
  if (any(bad)) theta <- lapply(theta, function(v) {
    v[rep(bad, K)] <- theta_fill(v)
    v
  })
  cdf <- function(y) {
    Fy <- distributions7::distrib_cdf(d, rep(y, K), theta)
    as.numeric(matrix(Fy, n, K) %*% w)
  }
  cond_q <- function(pr) {
    Q <- matrix(distributions7::distrib_quantile(d, rep(pr, n * K), theta),
                n, K)
    cbind(apply(Q, 1L, min), apply(Q, 1L, max))
  }
  discrete <- S7::S7_inherits(d, distributions7::discrete_distrib)
  quant <- function(pr) {
    br <- cond_q(pr)
    lo <- br[, 1L]
    hi <- br[, 2L]
    if (discrete) {
      lo <- lo - 1
      for (it in seq_len(200L)) {
        open <- hi - lo > 1
        if (!any(open, na.rm = TRUE)) break
        mid <- floor((lo + hi) / 2)
        g <- cdf(ifelse(open, mid, hi))
        up <- open & g >= pr
        hi[up] <- mid[up]
        lo[open & !up] <- mid[open & !up]
      }
      return(hi)
    }
    for (it in seq_len(100L)) {
      mid <- (lo + hi) / 2
      if (all(hi - lo <= 1e-10 * pmax(1, abs(mid)), na.rm = TRUE)) break
      g <- cdf(mid)
      up <- g >= pr
      hi[up] <- mid[up]
      lo[!up] <- mid[!up]
    }
    (lo + hi) / 2
  }
  a <- (1 - level) / 2
  med <- quant(0.5)
  lower <- quant(a)
  upper <- quant(1 - a)
  sdv <- tryCatch({
    mu <- matrix(mean(d, theta), n, K)
    vv <- matrix(distributions7::variance(d, theta), n, K)
    m1 <- as.numeric(mu %*% w)
    s2 <- as.numeric(vv %*% w) + as.numeric(mu^2 %*% w) - m1^2
    ifelse(is.finite(s2) & s2 >= 0, sqrt(s2), NA_real_)
  }, error = function(e) rep(NA_real_, n))
  out <- data.frame(fit = med, se = sdv, lower = lower, upper = upper)
  out[bad, ] <- NA_real_
  out
}


#' A Value Standing in for a Row That Cannot Be Read
#'
#' @param v A vector of parameter values.
#'
#' @return The first finite value of `v`, so a row with no covariance is
#'   carried through the vectorized arithmetic and blanked afterwards.
#'
#' @keywords internal
theta_fill <- function(v) {
  f <- v[is.finite(v)]
  if (length(f)) f[[1L]] else 1
}
