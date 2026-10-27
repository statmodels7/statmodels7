#' @include bootstrap_interval.R bootstrap_structural.R
NULL

# The BCa interval of confint(): the refits without one observation or one
# cluster, stored in the entry, and the limits from the adjusted levels.


#' Reject a Resampling Scheme That the BCa Interval Cannot Use
#'
#' @description
#' Signals an error for a scheme whose entry has no influence of the
#' observations to accelerate by.
#'
#' @details
#' The acceleration is computed from the refits of the fit without one
#' observation or one cluster. The parametric scheme draws new responses from
#' the model, so it does not resample the data. The cluster scheme beside a
#' [modelterms7::random()] term relabels the groups drawn more than once, and
#' deleting a cluster there has no counterpart in the stored coordinates.
#'
#' @param object The [StatmodFit()].
#' @param resample A [resample_method] object, or `NULL`, which passes.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
boot_bca_check <- function(object, resample) {
  if (is.null(resample)) return(invisible(NULL))
  if (S7::S7_inherits(resample, ParametricResample)) {
    stop(paste0("bca_interval() needs a scheme that resamples the data: ",
                "case_resample(), cluster_resample() or bayesian_resample().",
                "\n  The acceleration is the influence of the observations, ",
                "which a parametric draw does not have.\n  Use ",
                "percentile_interval() or basic_interval() with ",
                "parametric_resample()."), call. = FALSE)
  }
  if (S7::S7_inherits(resample, ClusterResample) &&
      nrow(random_modes(object@spec, "zero")) > 0L) {
    stop(paste0("bca_interval() is not available for cluster_resample() ",
                "beside a random() term."), call. = FALSE)
  }
  invisible(NULL)
}


#' The BCa Limits of One Quantity
#'
#' @description
#' Computes the two limits of the BCa interval of one coordinate from its
#' replicates and from the refits without one unit.
#'
#' @details
#' The bias correction is the normal quantile of the proportion of replicates
#' below the estimate, with the replicates equal to it counted by halves. The
#' acceleration is the skewness form of the jackknife values (Efron 1987). A
#' jackknife value that is not finite is left out. The adjusted levels are
#' undefined where \eqn{1 - a (z_0 + z) \le 0} or where the bias correction is
#' infinite, and the limits are then `NA`. The attribute `extreme` of the
#' result is `TRUE` when an adjusted level selects the first or the last
#' replicate.
#'
#' @param q The replicates of the quantity.
#' @param est Its estimate.
#' @param jk Its values in the refits without one unit.
#' @param alpha One minus the confidence level.
#'
#' @return A numeric vector of length two, with the attribute `extreme`.
#'
#' @references
#' Efron, B. (1987). Better bootstrap confidence intervals. *Journal of the
#' American Statistical Association*, 82(397), 171--185.
#'
#' @keywords internal
boot_bca_limits <- function(q, est, jk, alpha) {
  out <- c(NA_real_, NA_real_)
  z0 <- stats::qnorm(mean(q < est) + 0.5 * mean(q == est))
  jk <- jk[is.finite(jk)]
  if (!is.finite(z0) || length(jk) < 2L) return(out)
  u <- mean(jk) - jk
  den <- sum(u^2)
  a <- if (den > 0) sum(u^3) / (6 * den^1.5) else 0
  z <- stats::qnorm(c(alpha / 2, 1 - alpha / 2))
  d <- 1 - a * (z0 + z)
  if (any(d <= 0)) return(out)
  p <- stats::pnorm(z0 + (z0 + z) / d)
  out <- stats::quantile(q, p, names = FALSE)
  attr(out, "extreme") <- any(p * length(q) < 1 | (1 - p) * length(q) < 1)
  out
}


#' The Plan of a Delete-One Refit
#'
#' @description
#' Returns the function of the position of a unit that gives the draw of the
#' refit without that unit, in the form that [boot_body()] runs.
#'
#' @details
#' The function is built here so that its environment holds only the weights
#' and the units, which is what a worker process receives. The weights of the
#' rows of the unit are set to zero. For a model with a block that moves with
#' its coefficients the draw carries the position as its seed, so that the
#' restarts of the refit do not depend on the number of workers.
#'
#' @param w0 The prior weights of the fit.
#' @param units A list with the rows of each unit.
#' @param refresh Whether the model has a block that moves with its
#'   coefficients.
#'
#' @return A function of the position.
#'
#' @keywords internal
jackknife_plan <- function(w0, units, refresh) {
  force(w0); force(units); force(refresh)
  function(i) {
    w <- w0
    w[units[[i]]] <- 0
    pl <- list(w = w)
    if (refresh) pl$seed <- as.integer(i)
    pl
  }
}


#' The Refits Without One Observation or One Cluster
#'
#' @description
#' Refits the model once for each unit with the weights of its rows set to
#' zero, and returns the coefficients of the refits.
#'
#' @details
#' The units are the observations with a positive weight, or the clusters of
#' [cluster_resample()] that contain one. The refit is the one of the stored
#' entry (`refit`), and the replicate is checked as the engine checks it. The
#' random number stream of the session is not moved. A refit that fails leaves
#' a row of `NA`, and the attribute `failed` holds the number of such rows.
#'
#' @param object The [StatmodFit()].
#' @param resample The [resample_method] of the entry.
#' @param refit `"coefficients"` or `"full"`.
#'
#' @return A matrix with a row for each unit and a column for each coordinate
#'   of the replicates of the entry, with the attribute `failed`.
#'
#' @keywords internal
bootstrap_jackknife <- function(object, resample, refit) {
  spec0 <- object@spec
  design0 <- statmod_design(spec0)
  bootstrap_check(object, resample, refit, spec0, design0)
  light <- boot_light(object)
  w0 <- spec0@weights
  units <- if (S7::S7_inherits(resample, ClusterResample)) {
    id <- boot_cluster_ids(object, resample, spec0)
    unname(split(seq_along(id), match(id, unique(id))))
  } else as.list(seq_along(w0))
  units <- Filter(function(u) any(w0[u] > 0), units)
  draw <- boot_plan_wrap(
    bootstrap_plan(light, resample, spec0, design0,
                   random_modes(spec0, "zero"), TRUE), design0)
  method <- light@methods$smooth
  cfg <- inner_settings(method)
  beta0 <- unlist(light@coefficients[spec0@distrib@params], use.names = FALSE)
  ref <- list(def = integer(0), singular = FALSE)
  if (is.null(boot_mode(design0))) {
    ref <- boot_read_mode(light, spec0, design0, light@coefficients,
                          light@hyper, cfg)
  }
  plans <- jackknife_plan(w0, units, identical(boot_mode(design0), "refresh"))
  body <- boot_body(light, spec0, draw, refit, method, cfg, beta0, plans, ref,
                    TRUE, TRUE, FALSE)
  res <- worker_map(spec0, length(units), body, what = "jackknife refits")
  p <- length(beta0)
  J <- matrix(NA_real_, length(units), p)
  for (i in seq_along(units)) {
    if (res[[i]]$status == 0L && length(res[[i]]$beta) == p) {
      J[i, ] <- res[[i]]$beta
    }
  }
  attr(J, "failed") <- sum(is.na(J[, 1L]))
  J
}


#' The Stored Refits Without One Unit
#'
#' @description
#' Returns the matrix of [bootstrap_jackknife()] for the entry that a
#' [bootstrap_vcov()] object selects, computing it the first time and storing
#' it in the entry.
#'
#' @details
#' The refits do not depend on the number of replicates of the entry, so the
#' stored matrix stays when the entry is extended. A message says how many
#' refits were made. A warning is signaled when some refits failed; they are
#' left out of the acceleration.
#'
#' @param object The [StatmodFit()].
#' @param vt The resolved `BootstrapVcov` object.
#' @param interval The [bootstrap_interval], whose scheme may carry the data
#'   frame of a cluster scheme, which the stored entry does not.
#'
#' @return The matrix of [bootstrap_jackknife()], with the columns of the
#'   replicates of the entry.
#'
#' @keywords internal
boot_jackknife_get <- function(object, vt, interval) {
  env <- boot_env(object)
  ent <- boot_match(boot_entries(object), vt@resample, vt@refit)[[1L]]
  if (!is.null(ent@jackknife)) return(ent@jackknife)
  rs <- ent@resample
  if (S7::S7_inherits(rs, ClusterResample) && is.null(rs@data) &&
      S7::S7_inherits(interval@resample, ClusterResample)) {
    rs@data <- interval@resample@data
  }
  J <- bootstrap_jackknife(object, rs, ent@refit)
  message(sprintf(paste0("bootstrap: computed %d refits without one %s ",
                         "(%s, refit %s)."), nrow(J),
                  if (S7::S7_inherits(rs, ClusterResample)) "cluster" else
                    "observation", resample_label(rs), ent@refit))
  nf <- attr(J, "failed")
  if (nf > 0L) {
    warning(sprintf(paste0("%d of the %d refits without one unit failed; ",
                           "they are left out of the acceleration."), nf,
                    nrow(J)), call. = FALSE)
  }
  key <- boot_key(object, ent@resample, ent@refit)
  assign(key, S7::set_props(ent, jackknife = J), envir = env)
  J
}
