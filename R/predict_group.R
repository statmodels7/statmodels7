#' @include predict_prior.R
NULL

#' The Interval and Standard Deviation of a New Group's Parameter
#'
#' @description
#' For `interval = "group"`: the distribution of a new group's parameter
#' \eqn{\theta^* = h^{-1}(\eta^*)}, with
#' \eqn{\eta^* = \eta_0 + e + z^\top b}, where \eqn{e} is the estimation
#' error of the fixed part and \eqn{b} the effects of a new group, drawn from
#' the prior the fit estimated. Returns the ends of its interval and its
#' standard deviation, on the predictor's scale and on the parameter's.
#'
#' @details
#' The inverse link is monotone, so a quantile of \eqn{\theta^*} is the
#' inverse link of the same quantile of \eqn{\eta^*}, and the interval is
#' exact whatever the link. Under Gaussian priors \eqn{\eta^*} is Gaussian
#' with variance \eqn{s^2 = \mathrm{se}_0^2 + z^\top\Sigma_b z}, and the ends
#' are \eqn{h^{-1}(\eta_0 \pm z_{(1+\ell)/2}\, s)}. Under any other prior
#' \eqn{\eta^*} is a mixture of Gaussians: a scale mixture over the nodes
#' of [gamma_nodes()] for a Student t, the prior's own quantiles for one
#' univariate prior of another family ([group_quantile_nodes()]), and the
#' Monte Carlo draws of [predictive_mixture()] for anything else. Each end is
#' found by bisection on the mixture's distribution function.
#'
#' The standard deviation is that of \eqn{\theta^*}, not the delta method,
#' which under a log link and \eqn{s = 0.85} is 42 per cent too small. Under
#' a Gaussian \eqn{\eta^*} it is \eqn{s} for the identity link,
#' \eqn{e^{\eta_0 + s^2/2}\sqrt{e^{s^2} - 1}} for the log link and a
#' 40-node Gauss-Hermite rule for any other link. It is reported only where
#' it is known to exist:
#' \itemize{
#'   \item never where the predictor's domain is not the whole line (the
#'     square root and inverse links), since a Gaussian \eqn{\eta^*} leaves
#'     it;
#'   \item always for an inverse link bounded on both sides (logit, probit,
#'     cloglog, a bounded link), from the mixture's nodes;
#'   \item for Gaussian priors under every other link;
#'   \item under a prior that is not Gaussian and an unbounded link, only for
#'     the identity link, where the variance is
#'     \eqn{\mathrm{se}_0^2 + z^\top\mathrm{Var}(b)\,z}: a Student t gives
#'     \eqn{\Sigma\nu/(\nu - 2)}, `NA` for \eqn{\nu \le 2}, and a univariate
#'     prior of another family its own variance. Under a log link a Student t
#'     has no moment generating function, so \eqn{E[\theta^*]} is infinite,
#'     and the standard deviation is `NA`.
#' }
#'
#' @param object The [StatmodFit()].
#' @param spec,design The specification and design at the rows predicted.
#' @param aside The rows of [random_modes()] that are not `"conditional"`.
#' @param ep The predictors at the fixed part, as [statmod_eta()] returns
#'   them.
#' @param su The standard errors of [predict_se()].
#' @param level The interval's level.
#' @param ... Passed to [vcov.StatmodFit()].
#'
#' @return `su`, with `se_eta`, `eta_lower`, `eta_upper`, `se`, `lower` and
#'   `upper` replaced for every parameter a new group's effects reach.
#'
#' @keywords internal
group_interval <- function(object, spec, design, aside, ep, su, level, ...) {
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  n <- spec@n_obs
  parts <- prior_parts(spec, design, object, aside)
  gaussian <- !length(parts$t) && !length(parts$other)
  fixed <- predictive_cov(object, spec, design, list(), ...)$fixed
  zq <- stats::qnorm((1 + level) / 2)
  if (gaussian) {
    add <- predictive_cov(object, spec, design, parts$gaussian,
                          fixed = FALSE)$random
  } else {
    rand <- predictive_cov(object, spec, design, parts$gaussian,
                           fixed = FALSE)$random
    comps <- group_quantile_nodes(spec, ep$eta, fixed + rand, parts)
    if (is.null(comps)) {
      comps <- predictive_mixture(object, spec, design, aside, ep$eta, fixed)
    }
    heavy <- heavy_variance(spec, parts)
  }
  for (p in params) {
    g <- links[[p]]
    eta0 <- rep_len(as.numeric(ep$eta[[p]]), n)
    if (gaussian) {
      a <- add[p, p, ]
      if (!any(a != 0)) next
      s <- sqrt(su[[p]]$se_eta^2 + a)
      lo <- eta0 - zq * s
      hi <- eta0 + zq * s
      sd_theta <- group_sd_gaussian(g, eta0, s)
    } else {
      M <- vapply(comps, function(cc) rep_len(as.numeric(cc$eta[[p]]), n),
                  numeric(n))
      S <- vapply(comps, function(cc) {
        if (is.null(cc$C)) rep(0, n) else sqrt(pmax(cc$C[p, p, ], 0))
      }, numeric(n))
      M <- matrix(M, n)
      S <- matrix(S, n)
      W <- vapply(comps, function(cc) cc$weight, numeric(1))
      hv <- heavy[[p]]
      if (!any(rand[p, p, ] != 0) && !any(is.na(hv) | hv != 0)) next
      lo <- mixture_quantile(M, S, W, (1 - level) / 2)
      hi <- mixture_quantile(M, S, W, (1 + level) / 2)
      # the variance of eta*: the estimation error and the Gaussian effects,
      # plus the variance of the heavier priors where it exists
      v_eta <- fixed[p, p, ] + rand[p, p, ] + heavy[[p]]
      s <- sqrt(v_eta)
      sd_theta <- group_sd_mixture(g, M, S, W, s)
    }
    su[[p]]$se_eta <- s
    su[[p]]$eta_lower <- lo
    su[[p]]$eta_upper <- hi
    ends <- cbind(linkfunctions7::linkinv(g, lo),
                  linkfunctions7::linkinv(g, hi))
    su[[p]]$lower <- pmin(ends[, 1L], ends[, 2L])
    su[[p]]$upper <- pmax(ends[, 1L], ends[, 2L])
    su[[p]]$se <- sd_theta
  }
  su
}


#' The Kind of Inverse Link, for the Moments of a New Group's Parameter
#'
#' @param g A link.
#'
#' @return `"identity"`, `"log"`, `"bounded"` (an inverse link bounded on
#'   both sides), `"restricted"` (a predictor whose domain is not the whole
#'   line) or `"other"`.
#'
#' @keywords internal
link_kind <- function(g) {
  eb <- linkfunctions7::eta_bounds(g)
  if (any(is.finite(eb))) return("restricted")
  if (identical(g@link_name, "identity")) return("identity")
  if (identical(g@link_name, "log")) return("log")
  if (all(is.finite(g@link_bounds))) return("bounded")
  "other"
}


#' The Standard Deviation of h^-1(eta) for a Gaussian eta
#'
#' @param g The link.
#' @param m,s The mean and standard deviation of \eqn{\eta}, one a row.
#'
#' @return The standard deviation of \eqn{h^{-1}(\eta)}, `NA` where the
#'   predictor's domain is not the whole line.
#'
#' @keywords internal
group_sd_gaussian <- function(g, m, s) {
  switch(link_kind(g),
         restricted = rep(NA_real_, length(m)),
         identity = s,
         log = exp(m + s^2 / 2) * sqrt(expm1(s^2)),
         {
           mo <- gh_moments(g, matrix(m), matrix(s), 1)
           sqrt(pmax(mo$m2 - mo$m1^2, 0))
         })
}


#' The Standard Deviation of a New Group's Parameter Under a Mixture
#'
#' @param g The link.
#' @param M,S The components' means and standard deviations of \eqn{\eta},
#'   matrices with a row per observation and a column per component.
#' @param W The components' weights.
#' @param s The standard deviation of \eqn{\eta^*} itself, `NA` where it
#'   does not exist.
#'
#' @return The standard deviation of \eqn{h^{-1}(\eta^*)}, or `NA` where it
#'   is not known to exist. See [group_interval()].
#'
#' @keywords internal
group_sd_mixture <- function(g, M, S, W, s) {
  switch(link_kind(g),
         identity = ifelse(is.finite(s), s, NA_real_),
         bounded = {
           mo <- gh_moments(g, M, S, W)
           sqrt(pmax(mo$m2 - mo$m1^2, 0))
         },
         rep(NA_real_, nrow(M)))
}


#' The First Two Moments of h^-1(eta) Under a Mixture of Gaussians
#'
#' @description
#' Each component is integrated by a 40-node Gauss-Hermite rule, or read at
#' its mean where its standard deviation is zero (a Monte Carlo draw).
#'
#' @param g The link.
#' @param M,S,W As in [group_sd_mixture()].
#'
#' @return A list with `m1` and `m2`, one value a row.
#'
#' @keywords internal
gh_moments <- function(g, M, S, W) {
  q <- gauss_hermite(40L)
  x <- sqrt(2) * q$x
  w <- q$w / sqrt(pi)
  n <- nrow(M)
  m1 <- m2 <- numeric(n)
  for (c in seq_len(ncol(M))) {
    for (k in seq_along(x)) {
      th <- as.numeric(linkfunctions7::linkinv(g, M[, c] + S[, c] * x[k]))
      m1 <- m1 + W[c] * w[k] * th
      m2 <- m2 + W[c] * w[k] * th^2
    }
  }
  list(m1 = m1, m2 = m2)
}


#' The Variance the Priors That Are Not Gaussian Add to Each Predictor
#'
#' @description
#' \eqn{z^\top \mathrm{Var}(b)\, z} for every Student t and every other
#' prior [prior_parts()] lists, in each parameter's equation. A Student t
#' with scale \eqn{\Sigma} has variance \eqn{\Sigma\nu/(\nu - 2)}, infinite
#' for \eqn{\nu \le 2}; a univariate prior of another family has the
#' variance its family reports; a prior over several coordinates that is not
#' a Student t has none read here.
#'
#' @param spec The specification at the rows predicted.
#' @param parts What [prior_parts()] returns.
#'
#' @return A named list of numeric vectors, one value a row, `NA` where the
#'   variance is infinite or not read.
#'
#' @keywords internal
heavy_variance <- function(spec, parts) {
  params <- spec@distrib@params
  n <- spec@n_obs
  out <- stats::setNames(lapply(params, function(p) numeric(n)), params)
  addv <- function(mem, V) {
    at <- cumsum(c(0L, mem$dim))
    for (j in seq_len(nrow(mem))) {
      for (l in seq_len(nrow(mem))) {
        if (!identical(mem$param[j], mem$param[l])) next
        Zj <- random_within(spec, mem$param[j], mem$key[j])
        Zl <- random_within(spec, mem$param[l], mem$key[l])
        Vjl <- V[at[j] + seq_len(mem$dim[j]), at[l] + seq_len(mem$dim[l]),
                 drop = FALSE]
        p <- mem$param[j]
        out[[p]] <<- out[[p]] + rowSums((Zj %*% Vjl) * Zl)
      }
    }
  }
  for (tb in parts$t) {
    V <- if (tb$nu > 2) tb$Sigma * tb$nu / (tb$nu - 2) else
      matrix(NA_real_, nrow(tb$Sigma), ncol(tb$Sigma))
    addv(tb$members, V)
  }
  for (ob in parts$other) {
    D <- sum(ob$members$dim)
    v <- NA_real_
    if (D == 1L && "parent" %in% S7::prop_names(ob$penalty)) {
      v <- tryCatch(as.numeric(distributions7::variance(ob$penalty@parent,
                                                        ob$theta))[1L],
                    error = function(e) NA_real_)
      if (!is.finite(v)) v <- NA_real_
    }
    addv(ob$members, matrix(v, D, D))
  }
  out
}


#' A Quantile of a Mixture of Gaussians and Point Masses
#'
#' @description
#' The value at which \eqn{\sum_c w_c \Phi((e - m_c)/s_c)} reaches `p`, one
#' a row, a component with \eqn{s_c = 0} being a point mass at \eqn{m_c}.
#' Found by bisection, which needs no derivative and is exact to the
#' tolerance for a step function as well.
#'
#' @param M,S Matrices with a row per observation and a column per
#'   component.
#' @param W The components' weights, summing to one.
#' @param p The probability.
#'
#' @return One value a row.
#'
#' @keywords internal
mixture_quantile <- function(M, S, W, p) {
  W <- W / sum(W)
  lo <- apply(M - 12 * S, 1L, min) - 1e-8
  hi <- apply(M + 12 * S, 1L, max) + 1e-8
  cdf <- function(e) {
    Z <- (e - M) / S
    Fc <- ifelse(S > 0, stats::pnorm(Z), as.numeric(M <= e))
    as.numeric(Fc %*% W)
  }
  for (it in seq_len(200L)) {
    mid <- (lo + hi) / 2
    up <- cdf(mid) >= p
    hi[up] <- mid[up]
    lo[!up] <- mid[!up]
    if (all(hi - lo <= 1e-12 * pmax(1, abs(mid)))) break
  }
  hi
}


#' Quantile Nodes for One Univariate Prior That Is Not Gaussian
#'
#' @description
#' Where the only prior a new group's interval sets aside that is not
#' Gaussian is one univariate prior of a family other than the Student t,
#' the mixture over its effect is built from the prior's own quantiles,
#' \eqn{b_k = Q((k - 1/2)/K)} with \eqn{K = 1000} and equal weights, each
#' node a Gaussian component carrying the estimation error and the Gaussian
#' effects. The result does not depend on the random seed, where the Monte
#' Carlo draws of [predictive_mixture()] give the ends of a group interval
#' with a standard deviation of 0.030 to 0.037 over 20 seeds on a logistic
#' prior.
#'
#' @param spec The specification at the rows predicted.
#' @param eta A named list of the predictors' means.
#' @param base The covariance of the estimation error plus the Gaussian
#'   effects, an array `P x P x n`.
#' @param parts What [prior_parts()] returns.
#' @param K The number of nodes.
#'
#' @return A list of components in the shape [predictive_response()] reads,
#'   or `NULL` where the priors are of another kind.
#'
#' @keywords internal
group_quantile_nodes <- function(spec, eta, base, parts, K = 1000L) {
  if (length(parts$t) || length(parts$other) != 1L) return(NULL)
  ob <- parts$other[[1L]]
  mem <- ob$members
  if (nrow(mem) != 1L || mem$dim[1L] != 1L ||
      !"parent" %in% S7::prop_names(ob$penalty)) {
    return(NULL)
  }
  u <- (seq_len(K) - 0.5) / K
  bq <- tryCatch(as.numeric(distributions7::distrib_quantile(
    ob$penalty@parent, u, ob$theta)), error = function(e) NULL)
  if (is.null(bq) || length(bq) != K || !all(is.finite(bq))) return(NULL)
  z <- as.numeric(random_within(spec, mem$param[1L], mem$key[1L]))
  p <- mem$param[1L]
  lapply(seq_len(K), function(k) {
    e <- eta
    e[[p]] <- as.numeric(e[[p]]) + z * bq[k]
    list(eta = e, C = base, weight = 1 / K)
  })
}
