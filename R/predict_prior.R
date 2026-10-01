#' @include predict_interval.R
NULL

#' The Predictive Mixture Over a New Group's Effects
#'
#' @description
#' The components [predictive_response()] reads for one set of estimates:
#' the predictors' mean `eta` and the covariance `C0` of their estimation
#' error, widened by the effects of a new group for every term `random` sets
#' aside.
#'
#' @details
#' A Gaussian prior adds \eqn{z^\top\Sigma_b z} to the covariance, so one
#' component carries it. A Student t prior, univariate on one coefficient a
#' group or multivariate, is a scale mixture of Gaussians,
#' \deqn{b \mid w \sim \mathrm{N}(0, \Sigma/w), \qquad
#'   w \sim \mathrm{Gamma}(\nu/2, \nu/2),}
#' so it gives one component per node of a trapezoidal rule in
#' \eqn{\log w} ([gamma_nodes()]), each with the Gaussian covariance \eqn{z^\top\Sigma z/w} and the
#' node's weight; two such terms give the product of the two rules. Any other
#' prior (a Laplace, a prior shared by a label that is not Gaussian), or more
#' than two Student t terms, is averaged by Monte Carlo: `n_draw` draws of a
#' new group's effects, from [penalties7::penalty_draw()] or from the scale
#' mixture, each with a draw of the estimation error, so the interval
#' depends on the random seed. The quantiles of the mixture exist whatever the
#' prior; its standard deviation does not under a Student t with
#' \eqn{\nu \le 2}, and the result carries the attribute
#' `infinite_variance` to say so.
#'
#' @param fit The [StatmodFit()] whose hyperparameters give the priors.
#' @param spec,design The specification and design at the rows predicted.
#' @param aside The rows of [random_modes()] that are not `"conditional"`.
#' @param eta A named list of the predictors' means.
#' @param C0 The covariance of their estimation error, an array
#'   `P x P x n`, or `NULL` for none.
#' @param weight The total weight of the components returned.
#' @param n_draw The number of Monte Carlo draws, where they are needed.
#'
#' @return A list of components, each a list with `eta`, `C` (an array or
#'   `NULL`) and `weight`, with the attribute `infinite_variance`.
#'
#' @seealso [predictive_response()], which reads the components,
#'   [predict.StatmodFit()], the caller.
#'
#' @examples
#' set.seed(8)
#' gg <- data.frame(g = factor(rep(1:12, each = 5)), x = rnorm(60))
#' gg$y <- 1 + gg$x + 0.8 * rt(12, df = 3)[gg$g] + rnorm(60, sd = 0.4)
#' tp <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
#' fit <- statmod(y ~ x + random(~ 1 | g, distrib = tp),
#'                distributions7::gaussian1_distrib(), gg)
#' predict(fit, "response", data.frame(x = 0, g = "new"), random = "zero",
#'         interval = "prediction")
#'
#' @keywords internal
predictive_mixture <- function(fit, spec, design, aside, eta, C0,
                               weight = 1, n_draw = 2000L) {
  params <- spec@distrib@params
  P <- length(params)
  n <- spec@n_obs
  parts <- prior_parts(spec, design, fit, aside)
  base <- predictive_cov(fit, spec, design, parts$gaussian,
                         fixed = FALSE)$random
  if (!is.null(C0)) base <- base + C0
  inf_var <- any(vapply(parts$t, function(tb) tb$nu <= 2, logical(1)))
  if (length(parts$other) || length(parts$t) > 2L) {
    out <- prior_mc(spec, eta, base, c(parts$t, parts$other), weight,
                    as.integer(n_draw))
    attr(out, "infinite_variance") <- inf_var
    return(out)
  }
  out <- list(list(eta = eta, C = base, weight = weight))
  for (tb in parts$t) {
    q <- gamma_nodes(tb$nu / 2)
    nxt <- list()
    for (cc in out) {
      for (k in seq_along(q$w)) {
        bl <- list(list(members = tb$members, Sigma = tb$Sigma / q$w[k]))
        add <- predictive_cov(fit, spec, design, bl, fixed = FALSE)$random
        nxt[[length(nxt) + 1L]] <- list(eta = cc$eta, C = cc$C + add,
                                        weight = cc$weight * q$a[k])
      }
    }
    out <- nxt
  }
  attr(out, "infinite_variance") <- inf_var
  out
}


#' The Priors of the Effects a Prediction Sets Aside, by Kind
#'
#' @description
#' Sorts the terms `random` sets aside into Gaussian blocks, Student t
#' scale mixtures and other priors. See [predictive_mixture()].
#'
#' @param spec,design The specification and design.
#' @param fit The [StatmodFit()].
#' @param aside The rows of [random_modes()] that are not `"conditional"`.
#'
#' @return A list with `gaussian` (blocks in the shape of [random_blocks()]),
#'   `t` (each with `members`, `Sigma` the scale matrix and `nu`) and `other`
#'   (each with `members`, `penalty` and `theta`).
#'
#' @examples
#' set.seed(3)
#' gg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
#' gg$y <- 1 + gg$x + rnorm(8)[gg$g] + rnorm(40)
#' fit <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
#'                gg)
#' rm <- statmodels7:::random_modes(fit@spec, "zero")
#' lengths(statmodels7:::prior_parts(fit@spec, statmod_design(fit@spec),
#'                                   fit, rm))
#'
#' @keywords internal
prior_parts <- function(spec, design, fit, aside) {
  out <- list(gaussian = list(), t = list(), other = list())
  done <- character(0)
  units <- statmod_penalized(spec, design)
  for (i in seq_len(nrow(aside))) {
    p <- aside$param[i]
    k <- aside$key[i]
    if (paste(p, k) %in% done) next
    tm <- spec@terms[[p]][[k]]
    if (!is.na(modelterms7::term_tag(tm))) {
      cu <- Filter(function(u) !is.null(u$class) && any(vapply(
        u$class$pieces, function(pc) identical(pc$param, p) &&
          identical(pc$term, k), logical(1))), units)
      if (length(cu) != 1L) {
        stop(sprintf("'%s' has no single covariance class to read.", k),
             call. = FALSE)
      }
      u <- cu[[1L]]
      mem <- do.call(rbind, lapply(u$class$pieces, function(pc)
        data.frame(param = pc$param, key = pc$term, dim = as.integer(pc$dim),
                   stringsAsFactors = FALSE)))
      th <- as.list(fit@hyper[[u$param]][[u$key]])
      if (prior_is_gaussian(u$penalty, th)) {
        sub <- aside[paste(aside$param, aside$key) %in%
                       paste(mem$param, mem$key), , drop = FALSE]
        out$gaussian <- c(out$gaussian, random_blocks(spec, design, fit, sub))
      } else {
        out$other[[length(out$other) + 1L]] <- list(
          members = mem, penalty = u$penalty, theta = th)
      }
      done <- c(done, paste(mem$param, mem$key))
      next
    }
    uu <- Filter(function(u) identical(u$param, p) && identical(u$term, k),
                 units)
    if (length(uu) != 1L) {
      stop(sprintf("'%s' has no single prior to read.", k), call. = FALSE)
    }
    u <- uu[[1L]]
    th <- as.list(fit@hyper[[p]][[u$key]])
    d <- modelterms7::term_group(tm)$dim
    mem <- data.frame(param = p, key = k, dim = as.integer(d),
                      stringsAsFactors = FALSE)
    done <- c(done, paste(p, k))
    if (prior_is_gaussian(u$penalty, th)) {
      out$gaussian <- c(out$gaussian,
                        random_blocks(spec, design, fit, aside[i, , drop = FALSE]))
      next
    }
    tt <- prior_student(u$penalty, th, d)
    if (!is.null(tt)) {
      out$t[[length(out$t) + 1L]] <- c(list(members = mem), tt)
    } else {
      out$other[[length(out$other) + 1L]] <- list(
        members = mem, penalty = u$penalty, theta = th)
    }
  }
  out
}


#' Whether a Prior Is Gaussian
#'
#' @description
#' A penalty whose Hessian in the coefficients does not move with them and
#' has no kink is a Gaussian prior. The second condition matters: a Laplace
#' prior also has a constant Hessian, zero away from its kink.
#'
#' @param pen A penalty.
#' @param th Its hyperparameters.
#'
#' @return `TRUE` or `FALSE`.
#'
#' @examples
#' lp <- distributions7::fixed(distributions7::laplace_distrib(), mu = 0)
#' statmodels7:::prior_is_gaussian(penalties7::distrib_penalty(lp, n_coef = 3),
#'                                 list(sigma = 1))
#'
#' @keywords internal
prior_is_gaussian <- function(pen, th) {
  isTRUE(penalties7::beta_quadratic(pen, th)) && !penalty_has_kink(pen)
}


#' The Scale Matrix and Degrees of Freedom of a Student t Prior
#'
#' @param pen A penalty.
#' @param th Its hyperparameters, a named list.
#' @param d The number of coefficients a group carries.
#'
#' @return A list with `Sigma` and `nu` where the penalty is a Student t prior
#'   with one mixing variable a group (a multivariate t, or a univariate one
#'   on a single coefficient), and `NULL` otherwise.
#'
#' @examples
#' tp <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
#' pen <- penalties7::distrib_penalty(tp, n_coef = 4)
#' statmodels7:::prior_student(pen, list(sigma = 2, nu = 3), 1)
#'
#' @keywords internal
prior_student <- function(pen, th, d) {
  if (!"parent" %in% S7::prop_names(pen)) return(NULL)
  par <- pen@parent
  inner <- if ("parent_distrib" %in% S7::prop_names(par)) {
    par@parent_distrib
  } else par
  if (S7::S7_inherits(inner, distributions7::MvStudentTDistrib) &&
      inner@n_dim == d) {
    S <- distributions7::mv_sigma(par, th)
    return(list(Sigma = unname(as.matrix(S)), nu = th[["nu"]]))
  }
  # the univariate class is not exported, so it is asked for by name
  if (inherits(inner, "distributions7::StudentT1Distrib") && d == 1L &&
      all(c("sigma", "nu") %in% names(th))) {
    return(list(Sigma = matrix(th[["sigma"]]^2), nu = th[["nu"]]))
  }
  NULL
}


#' Nodes and Weights for a Gamma(a, a) Mixing Variable
#'
#' @description
#' A trapezoidal rule in \eqn{v = \log w} for the \eqn{\mathrm{Gamma}(a, a)}
#' density, over the range between its quantiles \eqn{10^{-12}} and
#' \eqn{1 - 10^{-12}}, the weights normalized to sum to one.
#'
#' @details
#' A prediction's distribution function, averaged over \eqn{w}, depends on
#' \eqn{w} through \eqn{\sqrt{w/(w + c)}}, which is not smooth at zero, and
#' its tails come from small \eqn{w}. A Gauss-Laguerre rule in \eqn{w}
#' converges slowly there; in \eqn{\log w} the integrand decays
#' exponentially at both ends and the trapezoidal rule converges
#' exponentially. Measured on the probability that a Gaussian measurement
#' around a Student t effect falls below the 2.5 per cent quantile, against
#' adaptive quadrature: at \eqn{\nu = 1} the error is 1e-04 with 40 nodes and
#' 2e-08 with 80, against 1.2e-02 for Gauss-Laguerre with 96; at
#' \eqn{\nu = 2.54} it is 1e-07 and 4e-13, against 5e-05.
#'
#' @param a The shape, \eqn{\nu/2}.
#' @param k The number of nodes.
#'
#' @return A list with `w`, the nodes, and `a`, the weights, which sum to one.
#'
#' @examples
#' q <- statmodels7:::gamma_nodes(1.5)
#' c(sum(q$a), sum(q$a * q$w))
#'
#' @keywords internal
gamma_nodes <- function(a, k = 64L) {
  lo <- log(stats::qgamma(1e-12, shape = a, rate = a))
  hi <- log(stats::qgamma(1 - 1e-12, shape = a, rate = a))
  w <- exp(seq(lo, hi, length.out = k))
  d <- stats::dgamma(w, shape = a, rate = a) * w
  list(w = w, a = d / sum(d))
}


#' The Predictive Mixture by Monte Carlo
#'
#' @description
#' `n_draw` point components: each a draw of the estimation error from
#' \eqn{\mathrm{N}(0, C)} at every row and a draw of a new group's effects
#' for every term in `heavy`. See [predictive_mixture()].
#'
#' @param spec The specification at the rows predicted.
#' @param eta A named list of the predictors' means.
#' @param base The covariance of the estimation error plus the Gaussian
#'   effects, an array `P x P x n`.
#' @param heavy The Student t and other entries of [prior_parts()].
#' @param weight The total weight.
#' @param n_draw The number of draws.
#'
#' @return A list of `n_draw` components, each with `C = NULL`.
#'
#' @examples
#' set.seed(8)
#' gg <- data.frame(g = factor(rep(1:12, each = 5)), x = rnorm(60))
#' gg$y <- 1 + gg$x + rlogis(12)[gg$g] + rnorm(60, sd = 0.4)
#' lp <- distributions7::fixed(distributions7::laplace_distrib(), mu = 0)
#' fit <- statmod(y ~ x + random(~ 1 | g, distrib = lp),
#'                distributions7::gaussian1_distrib(), gg)
#' predict(fit, "response", data.frame(x = 0, g = "new"), random = "zero",
#'         interval = "prediction")
#'
#' @keywords internal
prior_mc <- function(spec, eta, base, heavy, weight, n_draw) {
  params <- spec@distrib@params
  P <- length(params)
  n <- spec@n_obs
  M <- matrix(vapply(params, function(p) rep_len(as.numeric(eta[[p]]), n),
                     numeric(n)), n, P, dimnames = list(NULL, params))
  Z <- matrix(stats::rnorm(P * n_draw), P, n_draw)
  # one n x n_draw matrix of predictors for every parameter
  E <- stats::setNames(lapply(seq_len(P), function(j)
    matrix(M[, j], n, n_draw)), params)
  for (i in seq_len(n)) {
    Ci <- base[, , i]
    if (!all(is.finite(Ci))) {
      for (j in seq_len(P)) E[[j]][i, ] <- NA_real_
      next
    }
    ev <- eigen((Ci + t(Ci)) / 2, symmetric = TRUE)
    L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), P)
    dz <- L %*% Z
    for (j in seq_len(P)) E[[j]][i, ] <- E[[j]][i, ] + dz[j, ]
  }
  for (h in heavy) {
    mem <- h$members
    D <- sum(mem$dim)
    B <- if (!is.null(h$nu)) {
      w <- stats::rgamma(n_draw, shape = h$nu / 2, rate = h$nu / 2)
      ev <- eigen(h$Sigma, symmetric = TRUE)
      R <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), D)
      t(R %*% matrix(stats::rnorm(D * n_draw), D, n_draw)) / sqrt(w)
    } else {
      prior_draw_groups(h$penalty, h$theta, D, n_draw)
    }
    at <- cumsum(c(0L, mem$dim))
    for (j in seq_len(nrow(mem))) {
      Zj <- random_within(spec, mem$param[j], mem$key[j])
      Bj <- B[, at[j] + seq_len(mem$dim[j]), drop = FALSE]
      E[[mem$param[j]]] <- E[[mem$param[j]]] + Zj %*% t(Bj)
    }
  }
  lapply(seq_len(n_draw), function(k)
    list(eta = lapply(E, function(m) m[, k]), C = NULL,
         weight = weight / n_draw))
}


#' Draws of One Group's Effects From a Prior
#'
#' @param pen The penalty read as a prior.
#' @param th Its hyperparameters.
#' @param D The coordinates one group carries.
#' @param n_draw The number of draws.
#'
#' @return A matrix with `n_draw` rows and `D` columns, a group's coordinates
#'   side by side, each row an independent draw.
#'
#' @examples
#' lp <- distributions7::fixed(distributions7::laplace_distrib(), mu = 0)
#' pen <- penalties7::distrib_penalty(lp, n_coef = 5)
#' dim(statmodels7:::prior_draw_groups(pen, list(sigma = 1), 1, 12))
#'
#' @keywords internal
prior_draw_groups <- function(pen, th, D, n_draw) {
  rows <- list()
  got <- 0L
  while (got < n_draw) {
    v <- penalties7::penalty_draw(pen, th)
    if (anyNA(v)) {
      stop("A random effect's prior did not give a draw for every effect.",
           call. = FALSE)
    }
    # the coefficients are group-major: a group's coordinates together
    m <- matrix(v, ncol = D, byrow = TRUE)
    rows[[length(rows) + 1L]] <- m
    got <- got + nrow(m)
  }
  do.call(rbind, rows)[seq_len(n_draw), , drop = FALSE]
}
