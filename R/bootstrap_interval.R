#' @include bootstrap_fit.R inference_objects.R
NULL

# Bootstrap intervals for confint(): the percentile, basic and BCa intervals of
# the stored replicates, and the reuse of the stored replicates by the
# predictive bootstrap of pred_bootstrap().


#' @title Bootstrap Intervals for the Coefficients of a Fit
#' @name bootstrap_interval
#'
#' @description
#' Three intervals are read from the replicates that [add_bootstrap()] stores
#' in a fit: the percentile interval, the basic interval and the BCa interval.
#' Each is an object built by [percentile_interval()], [basic_interval()] or
#' [bca_interval()], and the object is passed as the argument `interval` of
#' [confint.StatmodFit()].
#'
#' `bootstrap_interval` is the abstract parent of the three classes. It cannot
#' be instantiated, and serves to test membership with `S7::S7_inherits()`.
#'
#' @details
#' Write \eqn{\hat\theta} for an estimate, \eqn{\theta^*_{(\alpha)}} for the
#' \eqn{\alpha} quantile of its bootstrap replicates and \eqn{\alpha = 1 -}
#' `level`. The percentile interval is
#' \deqn{[\theta^*_{(\alpha/2)},\; \theta^*_{(1 - \alpha/2)}],}
#' and the basic interval is
#' \deqn{[2\hat\theta - \theta^*_{(1 - \alpha/2)},\;
#'   2\hat\theta - \theta^*_{(\alpha/2)}],}
#' the reflection of the percentile interval about the estimate (Davison and
#' Hinkley 1997, chapter 5). The quantiles are those of [stats::quantile()]
#' with its default type.
#'
#' The BCa interval (bias-corrected and accelerated) takes the quantiles of
#' the replicates at adjusted levels,
#' \deqn{\Phi\!\left(z_0 + \frac{z_0 + z_{\alpha'}}{1 - a\,(z_0 + z_{\alpha'})}
#'   \right), \quad \alpha' \in \{\alpha/2,\; 1 - \alpha/2\},}
#' where \eqn{\Phi} is the standard normal distribution function and
#' \eqn{z_{\alpha'}} its quantile. The bias correction
#' \eqn{z_0 = \Phi^{-1}(\hat p)} uses the proportion \eqn{\hat p} of replicates
#' below the estimate, with the replicates equal to it counted by halves. The
#' acceleration is
#' \deqn{a = \frac{\sum_i (\bar\theta_{(\cdot)} - \hat\theta_{(-i)})^3}
#'   {6\,[\sum_i (\bar\theta_{(\cdot)} - \hat\theta_{(-i)})^2]^{3/2}},}
#' with \eqn{\hat\theta_{(-i)}} the estimate of the fit refitted without
#' observation \eqn{i} (the weight of the observation set to zero) and
#' \eqn{\bar\theta_{(\cdot)}} the mean of these estimates (Efron 1987). With
#' [cluster_resample()] the unit deleted is a cluster. The refits are made once
#' and stored in the entry, one per observation or cluster, with the settings
#' of the entry (`refit`), so that the first call of [confint.StatmodFit()]
#' costs as many refits as there are units, and later calls read them. The
#' interval is available for [case_resample()], [cluster_resample()] and
#' [bayesian_resample()]. [parametric_resample()] is rejected, because the
#' acceleration needs the influence of the observations of the data, and so is
#' [cluster_resample()] beside a [modelterms7::random()] term. A coordinate
#' whose estimate lies outside the range of its replicates, or whose
#' adjusted level is undefined, has `NA` limits.
#'
#' The quantity is the one that [confint.StatmodFit()] reports. With
#' `readable = TRUE` the readable quantities are computed on every replicate
#' and the quantiles are taken afterwards, so the percentile interval of a
#' transformed quantity is the transformed percentile interval. A coordinate
#' that is missing in some replicate (the individual random effects under
#' [cluster_resample()]) or in the fit (an aliased coefficient) has `NA` limits.
#'
#' The replicates are those of the stored entry that the arguments select, in
#' the way [bootstrap_vcov()] selects it: a `NULL` argument takes the value of
#' the single entry that matches the others, and the defaults when none
#' matches. An entry with fewer than `n` replicates is extended, and one with
#' more is read up to its first `n`. The default when nothing is stored is
#' `n = 1000`, because a tail of a 95 per cent interval holds only one
#' replicate in forty. A warning is signaled when a tail holds fewer than ten
#' replicates, that is, when \eqn{n \alpha / 2 < 10}.
#'
#' @param n A single whole number, the total number of replicates, or `NULL`.
#' @param resample A [resample_method] object, or `NULL`.
#' @param refit `"coefficients"`, `"full"` or `NULL`.
#' @param name A single string, the short name of the interval.
#' @param label A single string, the name under which an object prints.
#'
#' @return `percentile_interval()`, `basic_interval()` and `bca_interval()`
#'   return an object that inherits from `bootstrap_interval`, with the
#'   properties `name`, `label`, `n`, `resample` and `refit`.
#'   `bootstrap_interval` returns nothing: its constructor signals an error.
#'
#' @references Davison, A. C. and Hinkley, D. V. (1997). *Bootstrap Methods and
#'   Their Application*. Cambridge University Press.
#'
#'   Efron, B. (1987). Better bootstrap confidence intervals. *Journal of the
#'   American Statistical Association*, 82(397), 171--185.
#'
#' @seealso [confint.StatmodFit()], [add_bootstrap()], [bootstrap_vcov()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' fit <- add_bootstrap(fit, n = 100, resample = case_resample())
#' confint(fit, interval = percentile_interval(n = 100))
#' confint(fit, interval = basic_interval(n = 100))
#' confint(fit, interval = bca_interval(n = 100))
#' @export
bootstrap_interval <- S7::new_class(
  "bootstrap_interval",
  properties = list(name = S7::class_character, label = S7::class_character,
                    n = S7::class_any, resample = S7::class_any,
                    refit = S7::class_any),
  abstract = TRUE
)

#' @rdname bootstrap_interval
PercentileInterval <- S7::new_class("PercentileInterval",
                                    parent = bootstrap_interval)

#' @rdname bootstrap_interval
BasicInterval <- S7::new_class("BasicInterval", parent = bootstrap_interval)

#' @rdname bootstrap_interval
BcaInterval <- S7::new_class("BcaInterval", parent = bootstrap_interval)

#' @rdname bootstrap_interval
#' @export
percentile_interval <- function(n = NULL, resample = NULL, refit = NULL) {
  boot_interval_new(PercentileInterval, "percentile", n, resample, refit)
}

#' @rdname bootstrap_interval
#' @export
basic_interval <- function(n = NULL, resample = NULL, refit = NULL) {
  boot_interval_new(BasicInterval, "basic", n, resample, refit)
}

#' @rdname bootstrap_interval
#' @export
bca_interval <- function(n = NULL, resample = NULL, refit = NULL) {
  boot_interval_new(BcaInterval, "bca", n, resample, refit)
}

#' Build a Bootstrap Interval Object
#'
#' @description
#' Checks the arguments of [percentile_interval()], [basic_interval()] and
#' [bca_interval()] and builds the object.
#'
#' @param class The S7 class.
#' @param name The short name of the interval.
#' @param n,resample,refit As in [bootstrap_vcov()].
#'
#' @return An object that inherits from [bootstrap_interval].
#'
#' @keywords internal
boot_interval_new <- function(class, name, n, resample, refit) {
  if (!is.null(n)) n <- check_boot_n(n)
  if (!is.null(resample)) resample <- as_resample_method(resample)
  if (!is.null(refit)) refit <- check_boot_refit(refit)
  class(name = name, label = paste0(name, " bootstrap interval"), n = n,
        resample = resample, refit = refit)
}

#' @title Print a Bootstrap Interval
#' @name print.bootstrap_interval
#' @description Prints one line, the label of the interval and the settings of
#'   the bootstrap that it reads.
#' @param x A [bootstrap_interval] object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' print(percentile_interval(n = 500))
#' @keywords internal
S7::method(print, bootstrap_interval) <- function(x, ...) {
  cat("<interval> ", x@label, " (",
      if (is.null(x@n)) "n as stored" else paste(x@n, "replicates"), ", ",
      if (is.null(x@resample)) "resample as stored" else
        resample_label(x@resample), ", ",
      if (is.null(x@refit)) "refit as stored" else paste0("refit ", x@refit),
      ")\n", sep = "")
  invisible(x)
}

#' Check the Argument That Selects a Bootstrap Interval
#'
#' @description
#' Returns `x` if it is a [bootstrap_interval] object and signals an error that
#' lists the constructors otherwise.
#'
#' @param x The value of an `interval` argument.
#'
#' @return `x`.
#'
#' @keywords internal
as_boot_interval <- function(x) {
  if (S7::S7_inherits(x, bootstrap_interval)) return(x)
  stop("'interval' must be NULL or a bootstrap_interval object: ",
       "percentile_interval(), basic_interval() or bca_interval().",
       call. = FALSE)
}

#' The Bootstrap Intervals of a Table of Coefficients
#'
#' @description
#' Computes the table of [confint.StatmodFit()] for an `interval` argument that
#' is not `NULL`.
#'
#' @details
#' The table is that of the Wald route with the bootstrap standard error, which
#' supplies the rows, the estimates and the standard errors; its limits are
#' replaced by the quantiles of the replicates. The variance convention
#' `vcov_type` and a test other than Wald's do not apply to a bootstrap
#' interval, and giving one is an error.
#'
#' @param object The [StatmodFit()].
#' @param parm,level,readable As in [confint.StatmodFit()].
#' @param test The `test` argument of [confint.StatmodFit()].
#' @param interval A [bootstrap_interval] object.
#' @param vcov_given Whether the caller gave `vcov_type`.
#' @param ... Passed to [vcov.StatmodFit()].
#'
#' @return The data frame of [confint.StatmodFit()].
#'
#' @keywords internal
confint_interval <- function(object, parm, level, readable, test, interval,
                             vcov_given, ...) {
  interval <- as_boot_interval(interval)
  test <- as_test_statistic(test)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (vcov_given) {
    stop(paste0("'interval' reads the quantiles of the bootstrap replicates, ",
                "so 'vcov_type' does not apply.\n  The standard error ",
                "reported is the bootstrap one."), call. = FALSE)
  }
  if (!S7::S7_inherits(test, WaldTest)) {
    stop("'interval' and a test other than wald_test() cannot be combined: ",
         "both set the limits.", call. = FALSE)
  }
  bca <- S7::S7_inherits(interval, BcaInterval)
  if (bca) boot_bca_check(object, interval@resample)
  vt <- boot_interval_vcov(object, interval)
  if (bca) boot_bca_check(object, vt@resample)
  tb <- confint.StatmodFit(object, parm = parm, level = level,
                           vcov_type = vt, readable = readable,
                           test = wald_test(), ...)
  e <- boot_match(boot_entries(object), vt@resample, vt@refit)[[1L]]
  spec <- object@spec
  design <- statmod_design(spec)
  nm <- rownames(coef_labels(spec, design))
  B <- e@coefficients[seq_len(vt@n), , drop = FALSE]
  B[, nm %in% object@aliased] <- NA_real_
  Q <- if (isTRUE(readable)) boot_readable_values(object, spec, design, B) else B
  alpha <- 1 - level
  if (vt@n * alpha / 2 < 10 - 1e-8) {
    warning(sprintf(paste0("The %d replicates leave %.1f in each tail of a ",
                           "%g per cent interval.\n  The quantiles of a ",
                           "tail are unstable below 10 replicates: use n of ",
                           "at least %d."), vt@n, vt@n * alpha / 2,
                    100 * level, as.integer(ceiling(round(20 / alpha, 8)))),
            call. = FALSE)
  }
  basic <- S7::S7_inherits(interval, BasicInterval)
  if (bca) {
    JK <- boot_jackknife_get(object, vt, interval)
    JQ <- if (isTRUE(readable)) boot_readable_values(object, spec, design, JK) else JK
  }
  extreme <- FALSE
  for (i in seq_len(nrow(tb))) {
    j <- match(rownames(tb)[[i]], colnames(Q))
    est <- tb$estimate[[i]]
    if (is.na(j) || is.na(est) || anyNA(Q[, j])) {
      tb$lower[[i]] <- tb$upper[[i]] <- NA_real_
      next
    }
    q <- if (bca) {
      lim <- boot_bca_limits(Q[, j], est, JQ[, j], alpha)
      if (isTRUE(attr(lim, "extreme"))) extreme <- TRUE
      lim
    } else stats::quantile(Q[, j], c(alpha / 2, 1 - alpha / 2), names = FALSE)
    if (basic) q <- 2 * est - rev(q)
    tb$lower[[i]] <- q[[1L]]
    tb$upper[[i]] <- q[[2L]]
  }
  if (extreme) {
    warning(paste0("The BCa levels of some coordinate fall in the first or ",
                   "last replicate, so an extreme order statistic is a limit.",
                   "\n  Use more replicates."), call. = FALSE)
  }
  tb
}

#' The Bootstrap Variance That Supplies a Bootstrap Interval
#'
#' @description
#' Returns the resolved [bootstrap_vcov()] object whose entry the interval
#' reads, computing the replicates that are missing.
#'
#' @details
#' A `NULL` number of replicates takes that of the single entry that matches
#' the other settings, and 1000 when none matches.
#'
#' @param object The [StatmodFit()].
#' @param interval A [bootstrap_interval] object.
#'
#' @return A `BootstrapVcov` object with no `NULL` property.
#'
#' @keywords internal
boot_interval_vcov <- function(object, interval) {
  vt <- BootstrapVcov(name = "bootstrap", n = interval@n,
                      resample = interval@resample, refit = interval@refit)
  none <- !length(boot_match(boot_entries(object), vt@resample, vt@refit))
  if (is.null(vt@n) && none) vt@n <- 1000L
  # the BCa interval needs a scheme that resamples the data
  if (none && is.null(vt@resample) &&
      S7::S7_inherits(interval, BcaInterval)) vt@resample <- case_resample()
  boot_resolve(object, vt)
}

#' The Readable Quantities of Every Replicate
#'
#' @description
#' Returns the matrix whose columns are the quantities that
#' [confint.StatmodFit()] reports with `readable = TRUE`, evaluated at each
#' replicate. It is the matrix whose covariance is taken by
#' [boot_readable_cov()].
#'
#' @param object The [StatmodFit()].
#' @param spec,design Its specification and design.
#' @param B The matrix of replicates, over the joint coordinates.
#'
#' @return A matrix with a row for each replicate and the names of the readable
#'   quantities, or `B` when every readable quantity is a coordinate.
#'
#' @keywords internal
boot_readable_values <- function(object, spec, design, B) {
  rj <- readable_joint(spec, design, object)
  J <- rj$jacobian
  if (nrow(J) == 0L || ncol(J) != ncol(B)) return(B)
  nmq <- paste(rj$parameter, rj$name, sep = ":")
  plain <- all(rowSums(J != 0) == 1L) && all(J[J != 0] == 1)
  Q <- if (plain) {
    B[, apply(J != 0, 1L, which), drop = FALSE]
  } else {
    params <- spec@distrib@params
    lens <- vapply(object@coefficients[params], length, integer(1))
    grp <- rep(params, lens)
    t(vapply(seq_len(nrow(B)), function(b) {
      cf <- split(B[b, ], factor(grp, levels = params))
      cf <- Map(function(v, p) stats::setNames(v, names(object@coefficients[[p]])),
                cf, params)
      readable_joint(spec, design,
                     S7::set_props(object, coefficients = cf))$value
    }, numeric(nrow(J))))
  }
  colnames(Q) <- nmq
  Q
}

#' The Stored Replicates That a Predictive Bootstrap Can Reuse
#'
#' @description
#' Returns the first `n` replicates of the stored parametric entry as a list
#' of refits in the shape that the engine hands to the hook of
#' [predictive_bootstrap()], or `NULL` when the stored bootstrap cannot stand
#' for the one that the prediction would make.
#'
#' @details
#' The entry is used when its scheme is [parametric_resample()] and its refit
#' is `refit`. An entry with fewer than `n` replicates is extended first, from
#' the stream saved with it, and the extension is stored in the fit; the
#' stream of the session is restored afterwards.
#'
#' The stored replicates of a random-effect coordinate are those of the
#' prediction error (\eqn{\hat b^* - b^* + \hat b}), which differ from the
#' refitted coefficients that the predictive bootstrap reads. A model with
#' random-effect terms is therefore reused only when the entry holds the
#' refitted random effects (the property `raw` of [StatmodBootstrap]) and the
#' prediction sets aside every random-effect term, which is how the entry
#' redrew them. An entry without `raw`, or a prediction that keeps some term
#' conditional, is not reused.
#'
#' @param object The [StatmodFit()].
#' @param aside The rows of [random_modes()] that the prediction sets aside.
#' @param refit `"coefficients"` or `"full"`.
#' @param n The number of replicates wanted.
#'
#' @return `NULL`, or a list of `n` lists with `coefficients` and `hyper`.
#'
#' @keywords internal
boot_stored_replicas <- function(object, aside, refit, n) {
  spec0 <- object@spec
  rm0 <- random_modes(spec0, "zero")
  has_random <- nrow(rm0) > 0L
  if (has_random &&
      !setequal(paste(rm0$param, rm0$key), paste(aside$param, aside$key))) {
    return(NULL)
  }
  if (!has_random && nrow(aside)) return(NULL)
  ents <- boot_match(boot_entries(object), parametric_resample(), refit)
  if (length(ents) != 1L) return(NULL)
  e <- ents[[1L]]
  if (has_random && (!identical(e@random, "redrawn") || is.null(e@raw) ||
                      is.null(e@raw_cols))) {
    return(NULL)
  }
  if (e@n < n) {
    # the entry is extended from its saved stream, as bootstrap_vcov() does;
    # the stream of the session is left where it was
    seed <- if (exists(".Random.seed", envir = globalenv(),
                       inherits = FALSE)) get(".Random.seed", envir = globalenv())
    e <- boot_fill(object, n, parametric_resample(), refit, say = TRUE)
    if (!is.null(seed)) assign(".Random.seed", seed, envir = globalenv())
    if (has_random && (is.null(e@raw) || is.null(e@raw_cols))) return(NULL)
  }
  B <- e@coefficients
  if (has_random) B[, e@raw_cols] <- e@raw
  params <- spec0@distrib@params
  lens <- vapply(object@coefficients[params], length, integer(1))
  grp <- factor(rep(params, lens), levels = params)
  lapply(seq_len(n), function(b) {
    cf <- split(B[b, ], grp)
    cf <- Map(function(v, p) stats::setNames(v, names(object@coefficients[[p]])),
              cf, params)
    hy <- if (identical(refit, "full") && !is.null(e@hyper)) {
      utils::relist(e@hyper[b, ], skeleton = object@hyper)
    } else object@hyper
    list(coefficients = cf, hyper = hy)
  })
}
