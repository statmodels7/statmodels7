#' @include bootstrap_engine.R
NULL

#' @title Compute and Store the Bootstrap Replicates of a Fit
#'
#' @description
#' Refits the model on resampled data sets and stores the coefficients of the
#' replicates in the fit. The replicates are read by [bootstrap_vcov()] and
#' returned by [statmod_bootstrap()].
#'
#' @details
#' The fit keeps a collection of entries, one for each combination of
#' `resample` and `refit`. A call with a combination that is not in the
#' collection adds an entry, and a call never overwrites one. The collection is
#' kept in an environment that belongs to the fit object, so the function
#' changes the object in place and returns it, and a copy of the fit made with
#' [S7::set_props()] shares the environment. An entry is used only by a fit
#' with the coefficients and hyperparameters that it was computed for.
#'
#' `n` is the total number of replicates that is wanted. If the entry holds
#' \eqn{m} replicates, \eqn{n > m} computes the \eqn{n - m} missing ones,
#' \eqn{n = m} computes nothing, and \eqn{n < m} computes nothing and leaves
#' the \eqn{m} stored (the first \eqn{n} are the ones that
#' [bootstrap_vcov()] uses). Repeating a call does not change the result.
#'
#' The draws come from the random number stream of the session, one replicate
#' after another in the calling process, and the refits are distributed over
#' the `workers` of the fit (see [statmod()]) with the number of workers that
#' the fit carries. The result therefore depends on the state of
#' `.Random.seed` at the call, and not on the number of workers. The state at
#' the start of the stream and after the last replicate kept is stored in the
#' entry, and an extension continues from the second: 500 replicates and then
#' 300 more equal the 800 replicates of a single call from the same seed.
#' After the call `.Random.seed` is the state after the last replicate kept.
#'
#' A replicate that fails (an error in the refit, a mode error above
#' [mode_error_limit()], or a coefficient that is not identified in the
#' replicate and is in the fit) is replaced by the next draw of the stream. The
#' entry records the failures by cause. The function stops with an error when
#' the failures outnumber the successes after ten draws. Replacing a failed
#' replicate conditions the bootstrap on success, so it understates the
#' variability when failures are frequent.
#'
#' The refits use the `inner_information` of the fit. A fit with a response that
#' is not one number per observation, or (with `refit = "full"`) a [cv()]
#' criterion is rejected, and so is `case_resample()` beside a
#' [modelterms7::random()] term.
#'
#' A fit with a score-driven or a regime term ([modelterms7::gas()],
#' [modelterms7::regime()]) takes `parametric_resample()` only. Each replicate
#' is a new path of the series, drawn through [modelterms7::term_simulate()]
#' at the fitted coefficients and parameters of the term, so the recursion of
#' the filter runs on the draws; it is not the fitted filtered path with new
#' noise. A resample of rows or of groups breaks the recursion at every join, and
#' `case_resample()` and `cluster_resample()` signal an error that names
#' `parametric_resample()`. The effects of a [modelterms7::random()] term keep
#' their estimates in the draws. The replicates hold the free parameters of the
#' term after the coefficients, on the unconstrained scale, in the order of
#' [vcov()]. A replicate is kept when the Newton decrement of the refit over
#' the coefficients and the term's parameters is below [mode_error_limit()]
#' (a score-driven term), or when the refit is finite (a regime term).
#'
#' A fit with a break-point term ([modelterms7::seg()], [modelterms7::jump()],
#' [modelterms7::jseg()]) or a nonlinear term ([modelterms7::nl()]) is refitted
#' on every replicate from the break-points and coefficients of the fit, and the
#' refit then runs the restart cycle of [statmod()] (bootstrap restarting with
#' the budget `n_boot` of the terms, and the settling of the sharp
#' break-points), because a refit that stopped at the local optimum nearest to
#' the fit's positions would be another estimator. A position is a coefficient
#' of the term, so its variance is that of the replicates. The slot of a sharp
#' [modelterms7::jump()] or [modelterms7::jseg()], which the fit holds at a
#' working value, stores the position of the break-point in the replicates.
#' The generator of
#' the restart cycle is seeded from the draw of the replicate, so the result
#' does not depend on the number of workers. A replicate is kept when its
#' coefficients are finite.
#'
#' @param fit A [StatmodFit()].
#' @param n A single whole number, the total number of replicates wanted.
#' @param resample A [resample_method] object: [parametric_resample()] (the
#'   default), [case_resample()] or [cluster_resample()].
#' @param refit `"coefficients"` (the default) re-estimates the coefficients at
#'   the hyperparameters of the fit, and `"full"` selects the hyperparameters
#'   again on every replicate.
#'
#' @return `fit`, invisibly.
#'
#' @seealso [statmod_bootstrap()] to read the entries, [bootstrap_vcov()] for
#'   the variance, and [resample_method] for the schemes.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#'
#' set.seed(2)
#' fit <- add_bootstrap(fit, n = 20, resample = case_resample())
#' # extend the entry to 30 replicates: only the last 10 are computed
#' fit <- add_bootstrap(fit, n = 30, resample = case_resample())
#' statmod_bootstrap(fit)
#' @export
add_bootstrap <- function(fit, n = 200L, resample = parametric_resample(),
                          refit = c("coefficients", "full")) {
  if (!S7::S7_inherits(fit, StatmodFit)) {
    stop("'fit' must be a statmod fit.", call. = FALSE)
  }
  n <- check_boot_n(n)
  resample <- as_resample_method(resample)
  refit <- check_boot_refit(match.arg(refit))
  boot_fill(fit, n, resample, refit)
  invisible(fit)
}


#' @title Read the Stored Bootstrap of a Fit
#'
#' @description
#' Returns the table of the bootstrap entries that [add_bootstrap()] stored in
#' a fit, or one entry.
#'
#' @details
#' Without `resample` and `refit` the result is a data frame with one row for
#' each entry: the resampling scheme, the refit, the number of replicates, and
#' the number of replaced replicates in all and by cause. With them the result
#' is the entry that has those settings, a [StatmodBootstrap] object, and an
#' argument that is `NULL` matches any value. If no entry matches, or if more
#' than one does, an error is signaled; the second one prints the table.
#' [coef()] of an entry is the matrix of replicates.
#'
#' @param fit A [StatmodFit()].
#' @param resample A [resample_method] object, or `NULL`.
#' @param refit `"coefficients"`, `"full"` or `NULL`.
#'
#' @return A data frame, or a [StatmodBootstrap] object.
#'
#' @seealso [add_bootstrap()], [bootstrap_vcov()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' statmod_bootstrap(fit)
#' fit <- add_bootstrap(fit, n = 20, resample = case_resample())
#' statmod_bootstrap(fit)
#' head(coef(statmod_bootstrap(fit, resample = case_resample())), 3)
#' @export
statmod_bootstrap <- function(fit, resample = NULL, refit = NULL) {
  if (!S7::S7_inherits(fit, StatmodFit)) {
    stop("'fit' must be a statmod fit.", call. = FALSE)
  }
  if (!is.null(resample)) resample <- as_resample_method(resample)
  if (!is.null(refit)) refit <- check_boot_refit(refit)
  ents <- boot_entries(fit)
  if (is.null(resample) && is.null(refit)) return(boot_table(ents))
  cand <- boot_match(ents, resample, refit)
  if (!length(cand)) {
    stop("No stored bootstrap has these settings.\n",
         boot_table_text(ents), call. = FALSE)
  }
  if (length(cand) > 1L) boot_ambiguous(ents)
  cand[[1L]]
}


#' The Storage of the Bootstrap Entries of a Fit
#'
#' @description
#' `boot_env()` returns the environment in which a fit keeps its entries, or
#' `NULL` for a fit that has none (an object saved before the property
#' existed). `boot_fingerprint()` returns the numbers that identify the fit an
#' entry belongs to, and `boot_entries()` the entries that belong to the fit,
#' in alphabetical order of their keys.
#'
#' @param fit A [StatmodFit()].
#'
#' @return The environment, a numeric vector, or a list of
#'   [StatmodBootstrap] objects.
#'
#' @keywords internal
boot_env <- function(fit) {
  e <- tryCatch(S7::prop(fit, "bootstrap"), error = function(e) NULL)
  if (is.environment(e)) e else NULL
}

#' @rdname boot_env
#' @keywords internal
boot_fingerprint <- function(fit) {
  as.numeric(c(unlist(fit@coefficients, use.names = FALSE),
               unlist(fit@hyper, use.names = FALSE)))
}

#' @rdname boot_env
#' @keywords internal
boot_entries <- function(fit) {
  env <- boot_env(fit)
  if (is.null(env)) return(list())
  keys <- sort(ls(env, all.names = TRUE))
  ents <- mget(keys, envir = env)
  fp <- boot_fingerprint(fit)
  Filter(function(e) identical(e@fingerprint, fp), ents)
}

#' The Key of an Entry in the Collection
#'
#' @description
#' Combines the scheme, the refit and a summary of the fingerprint, so that a
#' copy of a fit with other estimates, which shares the environment, writes
#' beside the entries of the original and not over them.
#'
#' @param fit A [StatmodFit()].
#' @param resample,refit The settings.
#'
#' @return A string.
#'
#' @keywords internal
boot_key <- function(fit, resample, refit) {
  fp <- boot_fingerprint(fit)
  paste(resample_label(resample), refit, length(fp),
        sprintf("%.15g", sum(fp * cos(seq_along(fp)))), sep = "|")
}

#' @rdname boot_env
#' @param ents A list of entries.
#' @param resample,refit The settings to match, `NULL` for any.
#' @keywords internal
boot_match <- function(ents, resample, refit) {
  Filter(function(e) {
    (is.null(resample) ||
       identical(resample_label(e@resample), resample_label(resample))) &&
      (is.null(refit) || identical(e@refit, refit))
  }, unname(ents))
}

#' The Table of Bootstrap Entries
#'
#' @description
#' `boot_table()` returns the data frame that [statmod_bootstrap()] prints, and
#' `boot_table_text()` the same as text for an error message.
#'
#' @param ents A list of entries.
#'
#' @return A data frame, or a string.
#'
#' @keywords internal
boot_table <- function(ents) {
  ents <- unname(ents)
  fl <- function(e) boot_failures(e)
  data.frame(
    resample = vapply(ents, function(e) resample_label(e@resample),
                      character(1)),
    refit = vapply(ents, function(e) e@refit, character(1)),
    replicates = vapply(ents, function(e) e@n, integer(1)),
    replaced = vapply(ents, function(e) sum(fl(e)), numeric(1)),
    error = vapply(ents, function(e) fl(e)[["error"]], numeric(1)),
    mode_error = vapply(ents, function(e) fl(e)[["mode_error"]], numeric(1)),
    unidentified = vapply(ents, function(e) fl(e)[["unidentified"]],
                          numeric(1)),
    stringsAsFactors = FALSE)
}

#' @rdname boot_table
#' @keywords internal
boot_table_text <- function(ents) {
  if (!length(ents)) return("The fit has no stored bootstrap.")
  paste(c("The stored bootstraps are:",
          utils::capture.output(print(boot_table(ents)))), collapse = "\n")
}

#' @rdname boot_table
#' @keywords internal
boot_ambiguous <- function(ents) {
  stop("More than one stored bootstrap matches. Give 'resample' and ",
       "'refit'.\n", boot_table_text(ents), call. = FALSE)
}


#' Resolve a Bootstrap Variance Against the Stored Entries
#'
#' @description
#' Returns the [bootstrap_vcov()] object with `n`, `resample` and `refit` all
#' set, after making sure that the fit holds an entry with at least `n`
#' replicates for those settings.
#'
#' @details
#' A `NULL` argument takes the value of the single entry that matches the
#' arguments given, and the default (`n = 200`, [parametric_resample()],
#' `"coefficients"`) when none matches. Several matching entries are an error.
#' A missing entry, or one with fewer replicates than `n`, is computed with a
#' message, which is the only place where the engine runs on behalf of a
#' variance. An object that is already resolved finds its entry and computes
#' nothing, so [summary.StatmodFit()] and [confint.StatmodFit()] resolve their
#' `vcov_type` once and pass the result on.
#'
#' @param fit A [StatmodFit()].
#' @param vt A `BootstrapVcov` object.
#'
#' @return A `BootstrapVcov` object with no `NULL` property.
#'
#' @keywords internal
boot_resolve <- function(fit, vt) {
  ents <- boot_entries(fit)
  cand <- boot_match(ents, vt@resample, vt@refit)
  if (length(cand) > 1L) boot_ambiguous(ents)
  if (length(cand) == 1L) {
    e <- cand[[1L]]
    resample <- e@resample
    refit <- e@refit
    n <- if (is.null(vt@n)) e@n else vt@n
  } else {
    resample <- if (is.null(vt@resample)) parametric_resample() else vt@resample
    refit <- if (is.null(vt@refit)) "coefficients" else vt@refit
    n <- if (is.null(vt@n)) 200L else vt@n
  }
  boot_fill(fit, n, resample, refit, say = TRUE)
  BootstrapVcov(name = "bootstrap", n = as.integer(n), resample = resample,
                refit = refit)
}

#' Bring an Entry to a Number of Replicates
#'
#' @description
#' Looks the entry for `resample` and `refit` up in the collection of the fit,
#' computes the replicates that are missing for `n` and stores the result.
#'
#' @details
#' A new entry starts the stream at the current state of `.Random.seed`. An
#' existing entry continues from the state stored after its last replicate.
#' The write to the environment of the fit is made here, in the calling
#' process, after the engine has collected the refits.
#'
#' @param fit A [StatmodFit()].
#' @param n The total number of replicates wanted.
#' @param resample,refit The settings.
#' @param say Whether to report the replicates computed with a message.
#'
#' @return The entry.
#'
#' @keywords internal
boot_fill <- function(fit, n, resample, refit, say = FALSE) {
  env <- boot_env(fit)
  if (is.null(env)) {
    stop(paste0("This fit has no place to store a bootstrap: it was created ",
                "by an older version.\n  Fit the model again."), call. = FALSE)
  }
  key <- boot_key(fit, resample, refit)
  old <- if (exists(key, envir = env, inherits = FALSE)) {
    get(key, envir = env, inherits = FALSE)
  } else NULL
  m <- if (is.null(old)) 0L else old@n
  if (m >= n) return(old)
  eng <- bootstrap_engine(fit, n - m, resample, refit,
                          seed = if (is.null(old)) NULL else old@seed$end)
  nw <- nrow(eng$coefficients)
  if (say) {
    message(sprintf(paste0("bootstrap: computed %d replicates (%s, refit %s",
                           "%s)."), nw, resample_label(resample), refit,
                    if (m) sprintf(", %d already stored", m) else ""))
  }
  entry <- StatmodBootstrap(
    coefficients = if (is.null(old)) eng$coefficients else
      rbind(old@coefficients, eng$coefficients),
    hyper = if (is.null(old)) eng$hyper else if (is.null(eng$hyper)) NULL else
      rbind(old@hyper, eng$hyper),
    n = as.integer(m + nw), resample = boot_strip(resample), refit = refit,
    random = eng$random,
    status = c(if (!is.null(old)) old@status, eng$status),
    seed = list(start = if (is.null(old)) eng$seed_start else old@seed$start,
                end = eng$seed_end),
    fingerprint = boot_fingerprint(fit))
  assign(key, entry, envir = env)
  entry
}


#' The Covariance of the Stored Replicates
#'
#' @description
#' Returns the variance matrix that [vcov.StatmodFit()] reports for
#' [bootstrap_vcov()]: the covariance of the first `n` replicates of the entry,
#' on the coordinates of the other conventions.
#'
#' @details
#' An aliased coordinate is missing in its row and column. A coefficient under a
#' kinked penalty that the fit holds at zero is included, with the variance of
#' its replicates. With `readable = TRUE` the quantities that the model reports
#' are computed for each replicate and their covariance is taken: a term that
#' reports a quantity other than its coefficients ([modelterms7::term_readable()])
#' is evaluated at every replicate, which is exact and not a delta method. A
#' readable quantity that reads a missing coordinate is missing.
#'
#' @param object The [StatmodFit()].
#' @param vt A `BootstrapVcov` object.
#' @param readable,param As in [vcov.StatmodFit()].
#'
#' @return A square matrix.
#'
#' @keywords internal
boot_vcov_matrix <- function(object, vt, readable, param) {
  vt <- boot_resolve(object, vt)
  e <- boot_entries(object)
  e <- boot_match(e, vt@resample, vt@refit)[[1L]]
  B <- e@coefficients[seq_len(vt@n), , drop = FALSE]
  spec <- object@spec
  design <- statmod_design(spec)
  nm <- rownames(coef_labels(spec, design))
  if (identical(boot_mode(design), "structural")) {
    out <- boot_vcov_structural(object, B, readable, spec, design, nm)
    if (!is.null(param)) out <- select_parameter(out, param, spec)
    return(out)
  }
  B[, nm %in% object@aliased] <- NA_real_
  out <- if (isTRUE(readable)) {
    boot_readable_cov(object, spec, design, B)
  } else {
    stats::cov(B)
  }
  if (!is.null(param)) out <- select_parameter(out, param, spec)
  out
}

#' The Covariance of the Readable Quantities of the Replicates
#'
#' @description
#' Evaluates the readable quantities at every replicate and returns their
#' covariance.
#'
#' @details
#' When every readable quantity is one coefficient, the replicates are
#' reordered. Otherwise [readable_joint()] is evaluated at each replicate.
#'
#' @param object The [StatmodFit()].
#' @param spec,design Its specification and design.
#' @param B The matrix of replicates, over the joint coordinates.
#'
#' @return A square matrix with the names of [readable_vcov()].
#'
#' @keywords internal
boot_readable_cov <- function(object, spec, design, B) {
  rj <- readable_joint(spec, design, object)
  J <- rj$jacobian
  nmq <- paste(rj$parameter, rj$name, sep = ":")
  if (nrow(J) == 0L || ncol(J) != ncol(B)) {
    out <- stats::cov(B)
    return(out)
  }
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
  out <- stats::cov(Q)
  dimnames(out) <- list(nmq, nmq)
  out
}

#' What the Summary Reports About a Bootstrap Variance
#'
#' @description
#' Returns the number of replicates, how many were replaced, the relative Monte
#' Carlo error of a standard error in per cent, and the frequency with which
#' each coefficient under a kinked penalty was not zero.
#'
#' @param object The [StatmodFit()].
#' @param vt A resolved `BootstrapVcov` object.
#'
#' @return A list.
#'
#' @keywords internal
boot_summary_info <- function(object, vt) {
  e <- boot_match(boot_entries(object), vt@resample, vt@refit)[[1L]]
  lab <- coef_labels(object@spec, statmod_design(object@spec))
  B <- e@coefficients[seq_len(vt@n), , drop = FALSE]
  kin <- which(lab$kinked)
  list(n = vt@n, replaced = sum(boot_failures(e, vt@n)),
       mc = 100 / sqrt(2 * vt@n), resample = resample_label(e@resample),
       refit = e@refit,
       selection = if (length(kin)) {
         stats::setNames(colMeans(B[, kin, drop = FALSE] != 0),
                         rownames(lab)[kin])
       } else numeric(0))
}

#' Resolve the Variance Convention Once
#'
#' @description
#' Returns `vt` unchanged unless it is a [bootstrap_vcov()] object, which is
#' resolved by [boot_resolve()]. [summary.StatmodFit()] and
#' [confint.StatmodFit()] call it on entry so that the replicates are computed
#' once, however many times they read the variance.
#'
#' @param fit A [StatmodFit()].
#' @param vt A [vcov_method] object.
#'
#' @return A [vcov_method] object.
#'
#' @keywords internal
resolve_vcov_method <- function(fit, vt) {
  if (S7::S7_inherits(vt, BootstrapVcov)) boot_resolve(fit, vt) else vt
}

#' The Scheme as It Is Stored
#'
#' @description
#' Returns the scheme without the data frame of a cluster scheme, so that an
#' entry does not hold a copy of the data.
#'
#' @param x A [resample_method] object.
#'
#' @return A [resample_method] object.
#'
#' @keywords internal
boot_strip <- function(x) {
  if (S7::S7_inherits(x, ClusterResample)) x@data <- NULL
  x
}
