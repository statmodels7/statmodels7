#' @include sandwich.R
NULL

#' @title S7 Class for the Meat of a HAC Variance
#'
#' @description
#' The class that [hac_sandwich()] instantiates. It carries the kernel, the
#' bandwidth, and for a model without a structural term the time and the
#' group variables. Use the constructor, which sets the properties.
#'
#' @param type The string `"HAC"`.
#' @param kernel A single string, `"bartlett"`, `"parzen"` or `"truncated"`.
#' @param bandwidth `NULL` or a non-negative whole number.
#' @param time,cluster `NULL` or a one-sided formula with one variable.
#' @param data `NULL` or the data frame that holds those variables.
#'
#' @return An S7 object with the properties `type`, `kernel`, `bandwidth`,
#'   `time`, `cluster` and `data`.
#'
#' @seealso [hac_sandwich()], [frequentist_vcov()]
#'
#' @examples
#' hac_sandwich(bandwidth = 3)@kernel
#'
#' @name HacSandwich
#' @keywords internal
HacSandwich <- S7::new_class(
  "HacSandwich",
  properties = list(
    type = S7::class_character,
    kernel = S7::class_character,
    bandwidth = S7::new_union(NULL, S7::class_numeric),
    time = S7::new_union(NULL, S7::class_formula),
    cluster = S7::new_union(NULL, S7::class_formula),
    data = S7::new_union(NULL, S7::class_data.frame)
  ),
  validator = function(self) {
    if (length(self@kernel) != 1L ||
        !self@kernel %in% c("bartlett", "parzen", "truncated")) {
      return("@kernel must be \"bartlett\", \"parzen\" or \"truncated\".")
    }
    b <- self@bandwidth
    if (!is.null(b) && (length(b) != 1L || !is.finite(b) || b < 0 ||
                        b != round(b))) {
      return("@bandwidth must be NULL or one non-negative whole number.")
    }
    for (f in list(self@time, self@cluster)) {
      if (!is.null(f) && (length(f) != 2L || !is.name(f[[2L]]))) {
        return(paste0("must be a one-sided formula naming one variable, ",
                      "such as ~ t."))
      }
    }
    NULL
  }
)

#' @title HAC Sandwich Variance of the Coefficients
#'
#' @description
#' Select a heteroskedasticity and autocorrelation consistent meat for the
#' frequentist variance \eqn{V = V_b B V_b}, with \eqn{V_b = (H + S)^{-1}}
#' the inverse of the penalized information. Pass the object as the
#' `sandwich` argument of [frequentist_vcov()]. It is the counterpart of
#' [hc0_sandwich()] for dependent observations, and the estimator of Newey
#' and West (1987) with the Bartlett kernel.
#'
#' @details
#' The score \eqn{u_t} of observation \eqn{t} is the derivative of its
#' log-likelihood contribution in the whole vector of unknowns: the
#' coefficients of every distribution parameter, and the own parameters of a
#' structural term on the unconstrained scale. For a [modelterms7::gas()]
#' term the contribution of \eqn{t} depends on the parameters through the
#' filtered level, and \eqn{u_t} carries that dependence: it is the sum over
#' the equations of the derivative of the log-density in the predictor, on
#' the link scale, times the derivative of the predictor in the unknowns,
#' which the forward recursion of the filter supplies. For a
#' [modelterms7::regime()] term the contribution is the increment
#' \eqn{\log p(y_t \mid y_1, \dots, y_{t-1})} of the forward recursion, and
#' \eqn{u_t} its derivative, computed in forward mode through the recursion.
#' In both cases \eqn{\sum_t u_t} is the gradient of the log-likelihood.
#'
#' With a well specified filter the scores are martingale differences and
#' \eqn{B = \sum_t u_t u_t'} already estimates their variance. The HAC meat
#' adds the autocovariances,
#' \deqn{B = \Gamma_0 + \sum_{l=1}^{L} k\!\left(\frac{l}{L + 1}\right)
#'   (\Gamma_l + \Gamma_l'), \qquad \Gamma_l = \sum_t u_t u_{t-l}',}
#' where \eqn{k} is the kernel and \eqn{L} the bandwidth, the largest lag.
#' The Bartlett kernel is \eqn{k(x) = 1 - x}, the Parzen kernel is
#' \eqn{1 - 6x^2 + 6x^3} for \eqn{x \le 1/2} and \eqn{2(1 - x)^3} up to one,
#' and the truncated kernel is one. With `bandwidth = 0` the meat is
#' \eqn{\sum_t u_t u_t'}, the meat of [hc0_sandwich()]. No finite-sample
#' correction is applied.
#'
#' The default bandwidth is \eqn{L = \lfloor 4 (T / 100)^{2/9} \rfloor}, with
#' \eqn{T} the number of observations of a group (the mean number, when the
#' groups differ in length). A bandwidth larger than the longest group minus
#' one is reduced to that.
#'
#' With a `by` in the structural term the series are panels. A lag runs
#' within a group only, so the product of two scores of different groups
#' enters no \eqn{\Gamma_l}, and the groups are independent. Observations
#' are ordered within a group as the structural term orders them (its
#' `time`), and a lag is a step in that order, whatever the gaps between
#' times.
#'
#' A model without a structural term has no order of its own. `time`
#' names the variable that orders the observations and `cluster` the
#' variable that separates independent series; both are read from `data`,
#' or from the data of the call to [statmod()], which is evaluated again in
#' the environment of the model formula. For a model with a structural term
#' they are not accepted, since the term fixes the order and the groups.
#'
#' The prior weights of [statmod()] multiply the log-likelihood contribution
#' of an observation, and so its score; they are not otherwise interpreted.
#' Models with a term whose block moves with its coefficients (the
#' break-point terms and [modelterms7::nl()]) are rejected.
#'
#' @param kernel `"bartlett"` (the default), `"parzen"` or `"truncated"`.
#' @param bandwidth `NULL` (the default), which takes
#'   \eqn{\lfloor 4 (T/100)^{2/9} \rfloor}, or a non-negative whole number,
#'   the largest lag \eqn{L}.
#' @param time `NULL` (the default) or a one-sided formula with one variable,
#'   such as `~ t`, for a model without a structural term.
#' @param cluster `NULL` (the default) or a one-sided formula with one
#'   variable, such as `~ id`, that separates the series of a model without a
#'   structural term.
#' @param data `NULL` (the default) or the data frame that holds `time` and
#'   `cluster`, with one row for each observation of the fit.
#'
#' @return A `HacSandwich` object with the properties `type` (`"HAC"`),
#'   `kernel`, `bandwidth`, `time`, `cluster` and `data`.
#'
#' @references
#' Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite,
#' heteroskedasticity and autocorrelation consistent covariance matrix.
#' *Econometrica*, 55(3), 703--708.
#'
#' @seealso [frequentist_vcov()], which takes the object, and
#'   [hc0_sandwich()] for independent observations.
#'
#' @examples
#' set.seed(1)
#' n <- 120
#' dd <- data.frame(t = seq_len(n), x = rnorm(n))
#' e <- as.numeric(stats::filter(rnorm(n), 0.6, method = "recursive"))
#' dd$y <- 1 + 0.5 * dd$x + e
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' sqrt(diag(vcov(fit, vcov_type = frequentist_vcov(
#'   sandwich = hc0_sandwich()))))[1:2]
#' sqrt(diag(vcov(fit, vcov_type = frequentist_vcov(
#'   sandwich = hac_sandwich(time = ~ t)))))[1:2]
#' @export
hac_sandwich <- function(kernel = c("bartlett", "parzen", "truncated"),
                         bandwidth = NULL, time = NULL, cluster = NULL,
                         data = NULL) {
  kernel <- match.arg(kernel)
  HacSandwich(type = "HAC", kernel = kernel,
              bandwidth = if (is.null(bandwidth)) NULL else
                as.numeric(bandwidth),
              time = time, cluster = cluster, data = data)
}

S7::method(print, HacSandwich) <- function(x, ...) {
  cat("<sandwich> HAC, ", x@kernel, " kernel, bandwidth ",
      if (is.null(x@bandwidth)) "floor(4 (T/100)^(2/9))" else x@bandwidth,
      "\n", sep = "")
  invisible(x)
}


#' Reject a Fit That the HAC Variance Does Not Cover
#'
#' @description
#' Signals an error for a fit with a term whose block moves with its
#' coefficients, and for a model with no structural term and no `time`.
#'
#' @param fit The [StatmodFit()].
#' @param design The design of the fit.
#' @param sw The [HacSandwich] object.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
hac_check <- function(fit, design, sw) {
  if (length(attr(design, "refresh")) > 0L) {
    stop(paste0("A HAC variance is not available for a model with a term ",
                "whose block moves with its coefficients (seg(), jump(), ",
                "jseg() or nl())."), call. = FALSE)
  }
  structural <- length(attr(design, "structural")) > 0L
  if (structural && (!is.null(sw@time) || !is.null(sw@cluster))) {
    stop(paste0("hac_sandwich(time = , cluster = ) is for a model without a ",
                "structural term. The gas() or regime() term fixes the ",
                "order of the observations and the series."), call. = FALSE)
  }
  if (!structural && is.null(sw@time)) {
    stop(paste0("A model without a gas() or regime() term has no order of ",
                "its own. Give the variable that orders the observations ",
                "as hac_sandwich(time = ~ t)."), call. = FALSE)
  }
  invisible(NULL)
}


#' The Scores of the Observations in the Joint Vector of Unknowns
#'
#' @description
#' Returns the matrix whose row \eqn{t} is the derivative of the
#' log-likelihood contribution of observation \eqn{t}, weighted by its prior
#' weight, in the coefficients of every equation followed by the own
#' parameters of a structural term on the unconstrained scale.
#'
#' @details
#' Without a structural term the row is the family's score on the link scale
#' times the design row of each equation. With a filter it is the sum over the
#' equations of the score times the derivative of the predictor in the
#' unknowns, and the derivative of the filtered equation's predictor is the
#' forward Jacobian that [modelterms7::term_curvature()] returns. With a
#' [modelterms7::regime()] term it is the derivative of the increment of the
#' forward recursion, see [regime_scores()]. The columns of a level that an
#' intercept holds are dropped, as they are from the information.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#'
#' @return A matrix with `spec@n_obs` rows and a column for each unknown of
#'   the joint information.
#'
#' @keywords internal
joint_scores <- function(spec, coef, design = statmod_design(spec)) {
  params <- spec@distrib@params
  n <- spec@n_obs
  w <- spec@weights
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  nb <- sum(npar)
  ev <- statmod_eta(spec, design, coef)
  if (!length(ev$filters) && !length(ev$regimes)) {
    g <- statmod_score_obs(spec, coef, design)
    blocks <- list()
    for (q in params) {
      d <- design[[q]]
      if (d$npar == 0L) next
      blocks[[length(blocks) + 1L]] <- as_dense(d$X) * g[[q]]
    }
    return(do.call(cbind, blocks))
  }
  if (length(ev$regimes)) return(regime_scores(spec, coef, design, ev))

  f <- ev$filters[[1L]]
  ap <- match(f$param, params)
  np <- length(modelterms7::term_params(f$tm))
  m <- nb + np
  Vs <- lapply(seq_along(params), function(a) {
    M <- matrix(0, n, m)
    if (npar[a] > 0L) {
      M[, offs[a] + seq_len(npar[a])] <- as_dense(design[[params[a]]]$X)
    }
    M
  })
  gl <- distributions7::distrib_gradient(spec@distrib, spec@response,
                                         ev$theta, scale = "link",
                                         threads = spec@threads)
  H <- distributions7::distrib_hessian(spec@distrib, spec@response, ev$theta,
                                       scale = "link", threads = spec@threads)
  D3 <- distributions7::distrib_deriv3(spec@distrib, spec@response, ev$theta,
                                       scale = "link", threads = spec@threads)
  if (isTRUE(f$joint)) {
    cv <- joint_filter_curvature(spec, design, f, Vs, gl, H, D3, ev$theta)
    for (k in seq_along(f$params)) {
      Vs[[match(f$params[k], params)]] <- cv$jacobian[[k]]
    }
  } else {
    cv <- filter_curvature(spec, design, f, ap, Vs, gl, H, D3, ev$theta)
    Vs[[ap]] <- cv$jacobian
  }
  U <- matrix(0, n, m)
  for (a in seq_along(params)) {
    U <- U + (w * rep_len(as.numeric(gl[[params[a]]]), n)) * Vs[[a]]
  }
  held <- statmod_structural_state(design)$held[[f$term]]
  if (length(held)) {
    zn <- modelterms7::term_params(f$tm)
    U <- U[, -(nb + match(held, zn)), drop = FALSE]
  }
  U
}


#' The Scores of the Increments of a Regime Likelihood
#'
#' @description
#' Differentiates, in forward mode, the increment
#' \eqn{\log p(y_t \mid y_1, \dots, y_{t-1})} of the scaled forward recursion
#' of a [modelterms7::regime()] term, in the unknowns of [joint_scores()].
#'
#' @details
#' It is the first-order part of the recursion of
#' [modelterms7::term_hessian()], which accumulates the same quantities into
#' their sum and returns only the sum. The state \eqn{a_t} and its derivative
#' are carried, the increment is \eqn{\log c_t} with \eqn{c_t} the
#' normalizing constant of the step, and its derivative is
#' \eqn{\dot c_t / c_t}. The chart of the chain, the levels and their
#' derivatives are those of the term, from the internal function that
#' [modelterms7::term_hessian()] calls. The groups run in the order of the
#' term's blueprint. The sum of the rows, weighted, is the `gradient` of
#' [modelterms7::term_hessian()].
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#' @param ev The evaluated predictors, from [statmod_eta()].
#'
#' @return A matrix as [joint_scores()] returns it.
#'
#' @keywords internal
regime_scores <- function(spec, coef, design, ev) {
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  n <- spec@n_obs
  w <- spec@weights
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  nb <- sum(npar)
  r <- ev$regimes[[1L]]
  ap <- match(r$param, params)
  zn <- modelterms7::term_params(r$tm)
  np <- length(zn)
  m <- nb + np
  cols <- nb + seq_len(np)
  Vs <- lapply(seq_along(params), function(a) {
    M <- matrix(0, n, m)
    if (npar[a] > 0L) {
      M[, offs[a] + seq_len(npar[a])] <- as.matrix(design[[params[a]]]$X)
    }
    M
  })
  theta_at <- function(e, i) {
    stats::setNames(lapply(seq_along(params), function(q) {
      x <- if (q == ap) e else rep_len(ev$eta_static[[params[q]]], n)[i]
      linkfunctions7::linkinv(links[[params[q]]], x)
    }), params)
  }
  yof <- function(i) {
    if (is.matrix(spec@response)) spec@response[i, , drop = FALSE]
    else spec@response[i]
  }
  tm <- r$tm
  bp <- tm@blueprint
  v <- unlist(r$psi[zn])
  k <- tm@k
  chart <- utils::getFromNamespace(".regime_chart", "modelterms7")
  ch <- chart(v, k, zn, tm@chain)
  idx <- seq_len(n)
  LF <- matrix(0, n, k)
  GM <- vector("list", k)
  for (j in seq_len(k)) {
    e <- rep_len(r$eta_static, n) + ch$mu[j]
    LF[, j] <- as.numeric(distributions7::distrib_pdf(
      spec@distrib, yof(idx), theta_at(e, idx), log = TRUE,
      threads = spec@threads))
    g <- distributions7::distrib_gradient(spec@distrib, yof(idx),
                                          theta_at(e, idx), scale = "link",
                                          threads = spec@threads)
    gm <- matrix(0, n, m)
    for (q in seq_along(params)) {
      gm <- gm + rep_len(as.numeric(g[[params[q]]]), n) * Vs[[q]]
    }
    gm[, cols] <- gm[, cols] + outer(rep_len(as.numeric(g[[r$param]]), n),
                                     ch$dmu[j, ])
    GM[[j]] <- gm
  }
  P <- ch$P
  U <- matrix(0, n, m)
  for (rows in bp$order) {
    a <- ch$delta
    da <- matrix(0, k, m)
    da[, cols] <- t(do.call(rbind, ch$ddelta))
    for (tt in seq_along(rows)) {
      row <- rows[tt]
      lf <- LF[row, ]
      wv <- exp(lf - max(lf))
      if (tt == 1L) {
        pred <- a
        dpred <- da
      } else {
        pred <- as.numeric(a %*% P)
        dpred <- crossprod(P, da)
        for (l in seq_len(k)) {
          Dl <- matrix(ch$dPz[, l, ], k, ch$np)
          dpred[l, cols] <- dpred[l, cols] + as.numeric(a %*% Dl)
        }
      }
      atil <- wv * pred
      datil <- matrix(0, k, m)
      for (l in seq_len(k)) {
        datil[l, ] <- wv[l] * (GM[[l]][row, ] * pred[l] + dpred[l, ])
      }
      ct <- sum(atil)
      dct <- colSums(datil)
      U[row, ] <- w[row] * dct / ct
      a <- atil / ct
      da <- (datil - outer(atil, dct) / ct) / ct
    }
  }
  held <- statmod_structural_state(design)$held[[r$term]]
  if (length(held)) U <- U[, -(nb + match(held, zn)), drop = FALSE]
  U
}


#' The Series of a HAC Variance
#'
#' @description
#' Returns the rows of each independent series, in the order of the series.
#'
#' @details
#' A structural term gives them in its blueprint. For a model without one
#' they are the rows of each level of `cluster`, ordered by `time`, with ties
#' kept in the order of the data.
#'
#' @param fit The [StatmodFit()].
#' @param design The design of the fit.
#' @param sw The [HacSandwich] object.
#'
#' @return A list of integer vectors of row numbers.
#'
#' @keywords internal
hac_series <- function(fit, design, sw) {
  spec <- fit@spec
  n <- spec@n_obs
  su <- attr(design, "structural")
  if (length(su)) {
    ev <- statmod_eta(spec, design, fit@coefficients)
    tm <- c(ev$filters, ev$regimes)[[1L]]$tm
    return(unname(lapply(tm@blueprint$order, as.integer)))
  }
  data <- sw@data
  if (is.null(data)) {
    data <- tryCatch({
      env <- environment(spec@formula)
      if (is.null(env)) env <- globalenv()
      eval(fit@call$data, env)
    }, error = function(e) NULL)
  }
  if (!is.data.frame(data) || nrow(data) != n) {
    stop(paste0("The data frame with the time variable could not be ",
                "recovered from the call of the fit. Pass it as ",
                "hac_sandwich(time = , data = ), with one row for each of ",
                "the ", n, " observations."), call. = FALSE)
  }
  read <- function(f) {
    v <- all.vars(f)
    if (!v %in% names(data)) {
      stop(sprintf("The variable '%s' is not a column of the data.", v),
           call. = FALSE)
    }
    x <- data[[v]]
    if (anyNA(x)) {
      stop(sprintf("The variable '%s' has missing values.", v), call. = FALSE)
    }
    x
  }
  tt <- read(sw@time)
  grp <- if (is.null(sw@cluster)) rep(1L, n) else
    as.integer(factor(read(sw@cluster)))
  lapply(split(seq_len(n), grp), function(i) i[order(tt[i])])
}


#' The Kernel of a HAC Variance
#'
#' @param kernel The name of the kernel.
#' @param x The lags divided by `L + 1`, in `[0, 1]`.
#'
#' @return The weights.
#'
#' @keywords internal
hac_kernel <- function(kernel, x) {
  switch(kernel,
    bartlett = 1 - x,
    parzen = ifelse(x <= 0.5, 1 - 6 * x^2 + 6 * x^3, 2 * (1 - x)^3),
    truncated = rep(1, length(x)))
}


#' The HAC Meat
#'
#' @description
#' Computes \eqn{\Gamma_0 + \sum_{l=1}^{L} k(l/(L+1))(\Gamma_l + \Gamma_l')}
#' from the scores and the series.
#'
#' @param U The matrix of scores.
#' @param series The rows of each independent series, in time order.
#' @param kernel The name of the kernel.
#' @param bandwidth The largest lag, or `NULL` for the default.
#'
#' @return A list with the symmetric `meat` and the `bandwidth` used.
#'
#' @keywords internal
hac_meat <- function(U, series, kernel = "bartlett", bandwidth = NULL) {
  len <- vapply(series, length, integer(1))
  L <- if (is.null(bandwidth)) {
    floor(4 * (mean(len) / 100)^(2 / 9))
  } else bandwidth
  L <- min(L, max(len) - 1L)
  m <- ncol(U)
  B <- matrix(0, m, m)
  for (l in 0:L) {
    kw <- if (l == 0L) 1 else hac_kernel(kernel, l / (L + 1))
    if (kw == 0) next
    G <- matrix(0, m, m)
    for (ix in series) {
      T <- length(ix)
      if (T <= l) next
      Ut <- U[ix, , drop = FALSE]
      G <- G + crossprod(Ut[(l + 1L):T, , drop = FALSE],
                         Ut[seq_len(T - l), , drop = FALSE])
    }
    B <- B + if (l == 0L) G else kw * (G + t(G))
  }
  list(meat = B, bandwidth = L)
}


#' The HAC Variance of the Coefficients
#'
#' @description
#' Computes \eqn{V_b B V_b} for the coordinates in `cols`, with \eqn{V_b}
#' the inverse penalized information that [vcov.StatmodFit()] has built and
#' \eqn{B} the HAC meat of the scores.
#'
#' @param fit The [StatmodFit()].
#' @param design The design of the fit.
#' @param Vb The inverse penalized information over the coordinates `cols`.
#' @param cols The positions, among the unknowns of the joint information, of
#'   the coordinates that `Vb` covers.
#' @param sw The [HacSandwich] object.
#'
#' @return A symmetric matrix of the dimension of `Vb`.
#'
#' @seealso [hac_sandwich()], [joint_scores()]
#'
#' @keywords internal
hac_vcov <- function(fit, design, Vb, cols, sw) {
  spec <- fit@spec
  U <- joint_scores(spec, fit@coefficients, design)
  series <- hac_series(fit, design, sw)
  if (length(unlist(series)) != spec@n_obs) {
    stop("The series do not cover the observations of the fit.", call. = FALSE)
  }
  U <- U[, cols, drop = FALSE]
  B <- hac_meat(U, series, sw@kernel, sw@bandwidth)$meat
  V <- Vb %*% B %*% Vb
  (V + t(V)) / 2
}


#' The HAC Variance Over the Joint Unknowns
#'
#' @description
#' Computes \eqn{V_b B V_b} over every unknown of the joint information, for
#' the table of a structural term's own parameters.
#'
#' @param fit The [StatmodFit()].
#' @param design The design of the fit.
#' @param Vb The inverse penalized information over the unknowns of the joint
#'   information.
#' @param sw The [HacSandwich] object.
#'
#' @return A matrix of the dimension of `Vb`.
#'
#' @keywords internal
hac_joint <- function(fit, design, Vb, sw) {
  hac_vcov(fit, design, Vb, seq_len(nrow(Vb)), sw)
}
