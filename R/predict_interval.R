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
#' covariance to read here and is rejected by name; a group interval and a
#' prediction interval reach such a prior through [predictive_mixture()]
#' instead.
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
      # random_prior() is written for the marginal average and its refusal of
      # a heavy-tailed prior under an unbounded link names random = "zero"
      # as the remedy, which is what a caller here has already asked for
      pr <- tryCatch(random_prior(spec, design, fit, p, k), error = function(e) {
        if (grepl("is not Gaussian", conditionMessage(e), fixed = TRUE)) {
          return(list(gaussian = FALSE))
        }
        stop(e)
      })
      if (!isTRUE(pr$gaussian)) {
        stop(sprintf(paste0("The prior of '%s' is not Gaussian, so a new ",
                            "group's interval has no\n  covariance to read. ",
                            "interval = \"confidence\" is available, and a ",
                            "group the fit\n  saw is predicted with random = ",
                            "\"conditional\"."), k),
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
        v <- row_quad(as.matrix(da$X), Vab, as.matrix(db$X))
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
#' each expectation evaluated at the single node \eqn{\eta_{ci}} where
#' \eqn{C_{ci}} is zero at every row. Otherwise \eqn{C_{ci} = LL^\top} is
#' split by its eigenvectors: the directions after the first are a
#' Gauss-Hermite product grid of at most 64 nodes, and the first, the one of
#' largest variance, a 20-node Gauss-Hermite rule. Where that variance is
#' large against the spread of the conditional law, the law's distribution
#' function turns from 0 to 1 between two adjacent nodes, which the rule
#' integrates badly: with \eqn{r} the ratio of the two standard deviations of
#' a Gaussian law, its error in probability is 2e-11 at \eqn{r = 1}, 1.5e-03
#' at \eqn{r = 3} and 1.8e-02 at \eqn{r = 5}. A direction whose values at two
#' adjacent nodes differ by more than 0.3 is therefore integrated by
#' [numericals7::quad_vec()] instead, at its default tolerances. A
#' direction cell is left to Gauss-Hermite whatever its values when it is
#' among the lightest, whose weights add up to at most 1e-8. One
#' component carries the default interval, two sources of variation being
#' folded into its covariance; the parametric bootstrap gives one component
#' per replica. A quantile of \eqn{G_i} lies between the smallest and the
#' largest quantile of the conditional laws at the nodes, a mixture's
#' distribution function being an average of theirs. For a continuous family
#' it is found by Newton's method on \eqn{G_i(y) = p} with the mixture's
#' density as the derivative, started at the weighted median of those
#' quantiles and kept inside a bracket every evaluation shrinks; for a
#' discrete one by bisection on the integers, as the smallest value at which
#' \eqn{G_i} reaches the level. Only the family's distribution function,
#' density and quantile are read, so no location parameter is needed. The
#' standard deviation is
#' \eqn{\sqrt{E[\mathrm{Var}(Y\mid\eta)] + \mathrm{Var}(E[Y\mid\eta])}},
#' `NA` where the family's mean or variance does not exist.
#'
#' @param spec The specification at the rows predicted.
#' @param comps A list of components, each a list with `eta`, a named list of
#'   the predictors' means, `C`, their covariance as an array `P x P x n` or
#'   `NULL` for none, and optionally `weight`, equal weights otherwise. An
#'   attribute `infinite_variance` set to `TRUE` reports the standard
#'   deviation as `NA`.
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
  discrete <- S7::S7_inherits(d, distributions7::discrete_distrib)
  bad <- rep(FALSE, n)
  # every component becomes cells, each a base point B, a leading direction
  # A and a weight; a point component is one cell with no direction
  pts <- list()
  dirs <- list()
  for (cc in seq_along(comps)) {
    C <- comps[[cc]]$C
    eta0 <- comps[[cc]]$eta
    wt <- comps[[cc]]$weight
    if (is.null(wt)) wt <- 1 / length(comps)
    M <- matrix(vapply(params, function(p) rep_len(as.numeric(eta0[[p]]), n),
                       numeric(n)), n, P)
    bad <- bad | !apply(is.finite(M), 1L, all)
    if (is.null(C) || (all(is.finite(C)) && !any(C != 0))) {
      pts[[length(pts) + 1L]] <- list(B = M, w = wt)
      next
    }
    L <- array(0, c(n, P, P))
    for (i in seq_len(n)) {
      Ci <- C[, , i]
      if (!all(is.finite(Ci))) {
        bad[i] <- TRUE
        next
      }
      ev <- eigen((Ci + t(Ci)) / 2, symmetric = TRUE)
      L[i, , ] <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), P)
    }
    # the directions that carry variance at some row: the first, the largest,
    # is integrated adaptively where it needs it, the others on a
    # Gauss-Hermite grid
    sz <- apply(L^2, 3L, sum)
    r <- max(1L, sum(sz > 1e-14 * max(sz)))
    kr <- if (r > 1L) max(2L, min(8L, floor(64^(1 / (r - 1L))))) else 1L
    qg <- gauss_hermite(max(kr, 2L))
    grid <- if (r > 1L) {
      as.matrix(expand.grid(rep(list(seq_len(kr)), r - 1L)))
    } else matrix(0L, 1L, 0L)
    for (g in seq_len(nrow(grid))) {
      B <- M
      wg <- 1
      for (j in seq_len(r - 1L)) {
        B <- B + L[, , j + 1L] * (sqrt(2) * qg$x[grid[g, j]])
        wg <- wg * qg$w[grid[g, j]] / sqrt(pi)
      }
      dirs[[length(dirs) + 1L]] <- list(B = B, A = matrix(L[, , 1L], n, P),
                                        w = wt * wg)
    }
  }
  stack <- function(z, what) {
    array(vapply(z, function(e) e[[what]], numeric(n * P)),
          c(n, P, length(z)))
  }
  NP <- length(pts)
  ND <- length(dirs)
  BP <- stack(pts, "B")
  WP <- vapply(pts, function(e) e$w, numeric(1))
  BD <- stack(dirs, "B")
  AD <- stack(dirs, "A")
  WD <- vapply(dirs, function(e) e$w, numeric(1))
  # the direction cells left to Gauss-Hermite whatever their values: the
  # lightest, whose weights add up to at most 1e-8, so the error they can
  # carry into the mixture's distribution function is at most that
  ow <- order(WD)
  light <- rep(FALSE, ND)
  light[ow[cumsum(WD[ow]) <= 1e-8 * sum(c(WP, WD))]] <- TRUE
  gh <- gauss_hermite(20L)
  ol <- order(gh$x)
  xl <- sqrt(2) * gh$x[ol]
  wl <- gh$w[ol] / sqrt(pi)
  kl <- length(xl)
  theta_at <- function(E) {
    th <- stats::setNames(lapply(seq_len(P), function(j) as.numeric(
      linkfunctions7::linkinv(links[[params[j]]], E[[j]]))), params)
    lapply(th, function(v) {
      v[!is.finite(v)] <- theta_fill(v)
      v
    })
  }
  # the parameters at every row and point cell, then at every row,
  # direction cell and Gauss-Hermite node: row fastest, then cell, then node
  thP <- theta_at(lapply(seq_len(P), function(j) as.vector(BP[, j, ])))
  thD <- theta_at(lapply(seq_len(P), function(j) {
    as.vector(array(BD[, j, ], c(n, ND, kl)) +
                outer(matrix(AD[, j, ], n, ND), xl))
  }))
  cdf_fun <- function(y, th) distributions7::distrib_cdf(d, y, th)
  pdf_fun <- function(y, th) distributions7::distrib_pdf(d, y, th)
  # quad_vec() over the leading direction of the direction cells in `mark`
  adapt <- function(yv, mark, fun) {
    idx <- which(mark)
    ri <- (idx - 1L) %% n + 1L
    ci <- (idx - 1L) %/% n + 1L
    Bm <- matrix(vapply(seq_len(P), function(j) BD[cbind(ri, j, ci)],
                        numeric(length(idx))), length(idx), P)
    Am <- matrix(vapply(seq_len(P), function(j) AD[cbind(ri, j, ci)],
                        numeric(length(idx))), length(idx), P)
    yy <- yv[ri]
    f <- function(x, i) {
      xv <- as.vector(x)
      ii <- rep(i, length.out = length(xv))
      th <- theta_at(lapply(seq_len(P), function(j)
        Bm[ii, j] + Am[ii, j] * xv))
      stats::dnorm(xv) * fun(yy[ii], th)
    }
    numericals7::quad_vec(f, -Inf, rep(Inf, length(idx)))
  }
  # the mixture's distribution function at y, one value a row, and its
  # density where `dens`; a direction cell whose conditional law turns by
  # more than 0.3 between two Gauss-Hermite nodes is integrated by quad_vec()
  mix <- function(y, dens = FALSE) {
    yv <- rep(y, length.out = n)
    out <- list(F = numeric(n), f = numeric(n))
    if (NP) {
      yr <- rep(yv, NP)
      out$F <- as.numeric(matrix(cdf_fun(yr, thP), n, NP) %*% WP)
      if (dens) out$f <- as.numeric(matrix(pdf_fun(yr, thP), n, NP) %*% WP)
    }
    if (ND) {
      yr <- rep(yv, ND * kl)
      Fa <- array(cdf_fun(yr, thD), c(n, ND, kl))
      Fm <- matrix(0, n, ND)
      jump <- matrix(0, n, ND)
      for (k in seq_len(kl)) {
        Fm <- Fm + wl[k] * Fa[, , k]
        if (k > 1L) jump <- pmax(jump, abs(Fa[, , k] - Fa[, , k - 1L]))
      }
      mark <- jump > 0.3
      mark[bad, ] <- FALSE
      mark[, light] <- FALSE
      if (any(mark)) Fm[mark] <- adapt(yv, mark, cdf_fun)
      out$F <- out$F + as.numeric(Fm %*% WD)
      if (dens) {
        fa <- array(pdf_fun(yr, thD), c(n, ND, kl))
        fm <- matrix(0, n, ND)
        for (k in seq_len(kl)) fm <- fm + wl[k] * fa[, , k]
        if (any(mark)) fm[mark] <- adapt(yv, mark, pdf_fun)
        out$f <- out$f + as.numeric(fm %*% WD)
      }
    }
    out
  }
  # a bracket from the conditional quantiles, which the mixture's lies
  # between, and a start at their weighted median
  cond_q <- function(pr) {
    Q <- cbind(
      if (NP) matrix(distributions7::distrib_quantile(
        d, rep(pr, n * NP), thP), n, NP),
      if (ND) matrix(distributions7::distrib_quantile(
        d, rep(pr, n * ND * kl), thD), n, ND * kl))
    w <- c(WP, as.vector(outer(WD, wl)))
    list(lo = apply(Q, 1L, min), hi = apply(Q, 1L, max),
         start = vapply(seq_len(n), function(i) {
           o <- order(Q[i, ])
           Q[i, o][which(cumsum(w[o]) >= 0.5 * sum(w))[1L]]
         }, numeric(1)))
  }
  quant <- function(pr) {
    br <- cond_q(pr)
    lo <- br$lo
    hi <- br$hi
    if (discrete) {
      lo <- lo - 1
      for (it in seq_len(200L)) {
        open <- hi - lo > 1
        if (!any(open, na.rm = TRUE)) break
        mid <- floor((lo + hi) / 2)
        g <- mix(ifelse(open, mid, hi))$F
        up <- which(open & g >= pr)
        dn <- which(open & g < pr)
        hi[up] <- mid[up]
        lo[dn] <- mid[dn]
      }
      return(hi)
    }
    # Newton on G(y) = pr inside a bracket every evaluation shrinks; a step
    # that leaves the bracket is replaced by its midpoint
    y <- pmin(pmax(br$start, lo), hi)
    for (it in seq_len(100L)) {
      v <- mix(y, dens = TRUE)
      g <- v$F - pr
      up <- which(g >= 0)
      dn <- which(g < 0)
      hi[up] <- y[up]
      lo[dn] <- y[dn]
      done <- !is.finite(g) | abs(g) <= 1e-12 |
        hi - lo <= 1e-10 * pmax(1, abs(y))
      if (all(done)) break
      step <- y - g / v$f
      inside <- is.finite(step) & step > lo & step < hi
      y <- ifelse(done, y, ifelse(inside, step, (lo + hi) / 2))
    }
    y
  }
  a <- (1 - level) / 2
  med <- quant(0.5)
  lower <- quant(a)
  upper <- quant(1 - a)
  # the moments are smooth in the predictors, so the Gauss-Hermite nodes
  # serve
  sdv <- tryCatch({
    m1 <- m2 <- numeric(n)
    if (NP) {
      mu <- matrix(mean(d, thP), n, NP)
      vv <- matrix(distributions7::variance(d, thP), n, NP)
      m1 <- m1 + as.numeric(mu %*% WP)
      m2 <- m2 + as.numeric((vv + mu^2) %*% WP)
    }
    if (ND) {
      w <- as.vector(outer(WD, wl))
      mu <- matrix(mean(d, thD), n, ND * kl)
      vv <- matrix(distributions7::variance(d, thD), n, ND * kl)
      m1 <- m1 + as.numeric(mu %*% w)
      m2 <- m2 + as.numeric((vv + mu^2) %*% w)
    }
    s2 <- m2 - m1^2
    ifelse(is.finite(s2) & s2 >= 0, sqrt(s2), NA_real_)
  }, error = function(e) rep(NA_real_, n))
  # a Student t prior with nu <= 2 has no variance, and the mixture over its
  # quadrature nodes or its draws would report a finite one
  if (isTRUE(attr(comps, "infinite_variance"))) sdv <- rep(NA_real_, n)
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


#' A Quadratic Form Row by Row, Where the Variance Has Missing Entries
#'
#' @description
#' \eqn{x_{ia}^\top V x_{ib}} for every row \eqn{i}, with the entries of
#' \eqn{V} that are not finite read as missing: a row is `NA` only where it
#' reaches one of them, that is where both its coefficient in \eqn{x_a} and
#' its coefficient in \eqn{x_b} are non-zero.
#'
#' @details
#' A coefficient a kinked prior holds at its kink, such as a random effect
#' under a Laplace prior estimated at exactly zero, has no variance, and
#' [vcov.StatmodFit()] reports `NA` there. A prediction whose design does not
#' reach that coefficient, the typical group's under `random = "zero"` for
#' one, has a standard error all the same; requiring every entry of the block
#' to be finite made it `NA` as well.
#'
#' @param Xa,Xb Matrices with a row per observation and a column per
#'   coefficient of the two blocks.
#' @param V The block of the variance matrix, rows for `Xa` and columns for
#'   `Xb`.
#'
#' @return A numeric vector, one value a row.
#'
#' @keywords internal
row_quad <- function(Xa, V, Xb) {
  bad <- !is.finite(V)
  V0 <- V
  V0[bad] <- 0
  v <- rowSums((Xa %*% V0) * Xb)
  if (any(bad)) {
    hit <- rowSums(((Xa != 0) * 1) %*% (bad * 1) * ((Xb != 0) * 1)) > 0
    v[hit] <- NA_real_
  }
  v
}
