#' The Quantity That a Profile Holds
#'
#' @description
#' The class of one label of a profile, which records what the label is and
#' selects the graph that [autoplot.StatmodProfile()] draws for it. The five
#' subclasses are
#' \describe{
#'   \item{`ProfiledCoefficient`}{a coefficient without a penalty, of a linear
#'     predictor or of the subformula of a term;}
#'   \item{`ProfiledPenalizedCoefficient`}{a coefficient under a quadratic
#'     penalty, of [modelterms7::s()], [modelterms7::random()] or
#'     [modelterms7::ridge()];}
#'   \item{`ProfiledBreakpoint`}{the position of a break-point of
#'     [modelterms7::seg()] without a subformula, with the values of its
#'     covariate in `covariate_values`;}
#'   \item{`ProfiledHyperparameter`}{a hyperparameter that [reml()] or [ml()]
#'     estimates, on the scale of its link `link`;}
#'   \item{`ProfiledPathHyperparameter`}{the hyperparameter of a kinked
#'     penalty that a path chooses, on the scale of its link `link`.}
#' }
#'
#' @param label The label in `which`.
#' @param estimate The estimate on the scale of the profile.
#' @param standard_error The standard error of the estimate on that scale,
#'   `NA` where it is not available.
#' @param covariate_values The values of the covariate of a break-point.
#' @param piecewise_constant Whether the profile of a break-point is constant
#'   between two consecutive values of the covariate, which is the case of
#'   the position of a sharp [modelterms7::jump()] or [modelterms7::jseg()].
#' @param link The link of a hyperparameter, a [linkfunctions7::link] object.
#' @param criterion_name The name of the criterion that chose the
#'   hyperparameter of a path (`"BIC"`, `"AIC"` or `"CV"`).
#' @param path_record The points that the path scored, a data frame with the
#'   columns `value` (on the scale of the link), `criterion` and `se` (the
#'   standard error of a cross-validated criterion, `NA` otherwise).
#'
#' @return An object of the class.
#'
#' @keywords internal
profiled_quantity <- S7::new_class(
  "profiled_quantity",
  properties = list(label = S7::class_character,
                    estimate = S7::class_numeric,
                    standard_error = S7::class_numeric),
  abstract = TRUE
)

#' @rdname profiled_quantity
ProfiledCoefficient <- S7::new_class("ProfiledCoefficient",
                                     parent = profiled_quantity)

#' @rdname profiled_quantity
ProfiledPenalizedCoefficient <- S7::new_class("ProfiledPenalizedCoefficient",
                                              parent = profiled_quantity)

#' @rdname profiled_quantity
ProfiledBreakpoint <- S7::new_class("ProfiledBreakpoint",
  parent = profiled_quantity,
  properties = list(covariate_values = S7::class_numeric,
                    piecewise_constant = S7::class_logical))

#' @rdname profiled_quantity
ProfiledHyperparameter <- S7::new_class("ProfiledHyperparameter",
  parent = profiled_quantity,
  properties = list(link = S7::class_any))

#' @rdname profiled_quantity
ProfiledPathHyperparameter <- S7::new_class("ProfiledPathHyperparameter",
  parent = profiled_quantity,
  properties = list(link = S7::class_any,
                    criterion_name = S7::class_character,
                    path_record = S7::class_data.frame))


#' The Profile Likelihood Returned by statmod_profile()
#'
#' @description
#' The object that [statmod_profile()] returns.
#'
#' @details
#' With `joint = TRUE` the table `points` has one column for each label,
#' named by the label, and one row for each combination of values held. With
#' `joint = FALSE` each label is profiled by itself, with the other labels
#' free, and `points` has the columns `profiled_label`, the label held in the
#' row, and `value`, its value. The columns that follow are, in both cases,
#' `profile_loglik`, `profile_penalty`, `profile_pen_loglik`,
#' `outer_criterion` and `converged`, described at [statmod_profile()].
#' `unrestricted_point` has the same columns at the fit itself, with one row
#' for each label when `joint = FALSE`.
#'
#' @param points The data frame of the points.
#' @param profiled_labels The labels of `which`.
#' @param profiled_quantities A list of [profiled_quantity] objects, one for
#'   each label.
#' @param profile_column The name of the column of `points` that is the
#'   profile of each label, `"outer_criterion"` or `"profile_pen_loglik"`, a
#'   character vector named by the labels.
#' @param unrestricted_point The data frame of the unrestricted fit.
#' @param profile_method The [profile_method] object.
#' @param hyper_handling `"held"` or `"refit"`, the argument `hyper` of
#'   [statmod_profile()].
#' @param outer_method The outer criterion object of the fit, or `NULL`.
#' @param is_likelihood_ratio Whether the profile of each label gives a
#'   likelihood ratio statistic, so that confint() inverts it, a logical vector
#'   named by the labels.
#' @param joint Whether the labels were held jointly.
#'
#' @return A `StatmodProfile` object.
#'
#' @seealso [statmod_profile()], [confint.StatmodProfile()],
#'   [autoplot.StatmodProfile()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' prof <- statmod_profile(fit, which = "mu:x",
#'                         values = profile_range(width = 2, n = 9))
#' prof@profile_column
#' head(as.data.frame(prof), n = 3)
#' @keywords internal
#' @export
StatmodProfile <- S7::new_class("StatmodProfile",
  properties = list(
    points = S7::class_data.frame,
    profiled_labels = S7::class_character,
    profiled_quantities = S7::class_list,
    profile_column = S7::class_character,
    unrestricted_point = S7::class_data.frame,
    profile_method = profile_method,
    hyper_handling = S7::class_character,
    outer_method = S7::class_any,
    is_likelihood_ratio = S7::class_logical,
    joint = S7::class_logical
  )
)


#' The Kind of Each Label of a Profile
#'
#' @description
#' Builds the [profiled_quantity] object of each label of `which`, from the
#' context of the call.
#'
#' @param ctx The context of the call.
#'
#' @return A list of [profiled_quantity] objects, named by the labels.
#'
#' @keywords internal
profile_quantities <- function(ctx) {
  fit <- ctx$fit
  spec <- ctx$spec
  design <- statmod_design(spec)
  pen <- integer(0)
  for (u in statmod_penalized(spec, design)) {
    if (isTRUE(u$structural)) next
    pen <- c(pen, as.integer(u$index), as.integer(u$beta_index))
  }
  bps <- ctx$bps
  se <- if (length(ctx$cidx)) {
    tryCatch(sqrt(diag(stats::vcov(fit))), error = function(e) NULL)
  }
  links <- attr(ctx$hyp, "links")
  out <- list()
  for (i in seq_along(ctx$which)) {
    lab <- ctx$which[i]
    est <- ctx$estimate[i]
    if (i %in% ctx$hidx) {
      r <- ctx$hyp[match(i, ctx$hidx), , drop = FALSE]
      link <- links[[r$row]]
      out[[lab]] <- if (isTRUE(r$path)) {
        ProfiledPathHyperparameter(
          label = lab, estimate = est, standard_error = NA_real_, link = link,
          criterion_name = if (is.null(ctx$path_method)) "" else
            toupper(ctx$path_method@kind),
          path_record = profile_path_record(fit, r, link))
      } else {
        ProfiledHyperparameter(label = lab, estimate = est,
                               standard_error = NA_real_, link = link)
      }
      next
    }
    if (i %in% ctx$xidx) {
      r <- ctx$xpos[match(i, ctx$xidx), ]
      out[[lab]] <- ProfiledBreakpoint(label = lab, estimate = est,
                                       standard_error = NA_real_,
                                       covariate_values = r$covariate[[1L]],
                                       piecewise_constant = TRUE)
      next
    }
    pos <- ctx$at[match(i, ctx$cidx)]
    s <- if (is.null(se)) NA_real_ else unname(se[lab])
    if (!length(s) || !is.finite(s)) s <- NA_real_
    r <- bps[bps$at == pos & !bps$sharp, , drop = FALSE]
    out[[lab]] <- if (nrow(r)) {
      ProfiledBreakpoint(label = lab, estimate = est, standard_error = s,
                         covariate_values = r$covariate[[1L]],
                         piecewise_constant = FALSE)
    } else if (pos %in% pen) {
      ProfiledPenalizedCoefficient(label = lab, estimate = est, standard_error = s)
    } else {
      ProfiledCoefficient(label = lab, estimate = est, standard_error = s)
    }
  }
  out
}


#' Whether the Drop of a Profile Is a Likelihood Ratio
#'
#' @description
#' Twice the drop of the profile from its maximum is a likelihood ratio
#' statistic except in two cases: when a label is the hyperparameter of a
#' kinked penalty that a path chooses, whose profile includes the penalty at
#' the value held, and when the profile is the criterion [reml()] and a label
#' is a coefficient without a penalty, whose hold changes the space of the
#' fixed effects over which the criterion integrates.
#'
#' @param ctx The context of the call.
#' @param quantities The [profiled_quantity] objects of the labels.
#' @param column The column of the profile.
#'
#' @return A single logical value.
#'
#' @keywords internal
profile_is_lr <- function(ctx, quantities, column) {
  is <- function(cls) vapply(quantities, function(q) S7::S7_inherits(q, cls), TRUE)
  if (any(is(ProfiledPathHyperparameter))) return(FALSE)
  unpen <- is(ProfiledCoefficient) | is(ProfiledBreakpoint)
  !(identical(column, "outer_criterion") && !is.null(ctx$outer) &&
      identical(ctx$outer@kind, "reml") && any(unpen))
}


S7::method(print, StatmodProfile) <- function(x, ...) {
  k <- length(x@profiled_labels)
  pts <- x@points
  cat("<profile> ", k, if (k == 1L) " label" else " labels",
      if (k > 1L) (if (x@joint) ", held jointly" else ", held one at a time"),
      ", method ", x@profile_method@name, ", hyperparameters ", x@hyper_handling,
      if (!is.null(x@outer_method)) paste0(", outer criterion ",
                                           x@outer_method@kind, "()"),
      "\n", sep = "")
  for (q in x@profiled_quantities) {
    cat(sprintf("  %-28s %s, estimate %s, profile %s\n", q@label,
                profile_kind_name(q), format(signif(q@estimate, 6)),
                x@profile_column[[q@label]]))
  }
  nc <- sum(!pts$converged, na.rm = TRUE)
  cat(sprintf("%d points, %d not converged\n", nrow(pts), nc))
  if (!all(x@is_likelihood_ratio)) {
    cat("no likelihood ratio, and no interval, for ",
        paste(names(x@is_likelihood_ratio)[!x@is_likelihood_ratio],
              collapse = ", "), "\n", sep = "")
  }
  print(utils::head(pts, n = 10L))
  if (nrow(pts) > 10L) cat("... ", nrow(pts) - 10L, " more rows\n", sep = "")
  invisible(x)
}


#' The Name of the Kind of a Profiled Quantity
#'
#' @param q A [profiled_quantity] object.
#'
#' @return A string.
#'
#' @keywords internal
profile_kind_name <- function(q) {
  if (S7::S7_inherits(q, ProfiledPenalizedCoefficient)) "penalized coefficient"
  else if (S7::S7_inherits(q, ProfiledCoefficient)) "coefficient"
  else if (S7::S7_inherits(q, ProfiledBreakpoint)) "break-point position"
  else if (S7::S7_inherits(q, ProfiledPathHyperparameter)) "hyperparameter chosen by a path"
  else "hyperparameter"
}


S7::method(as.data.frame, StatmodProfile) <- function(x, ...) {
  x@points
}


#' The One-Dimensional Profile of One Label
#'
#' @description
#' Returns the profile of one label as a function of its value alone. With
#' `joint = FALSE` these are the rows in which the label is held. With a joint
#' profile the value of the profile at each value of the label is the maximum
#' over the rows with that value, which is the profile over the other labels
#' as far as the grid resolves it. The rows without a value of the profile
#' are dropped, and the estimate is added with the value of the fit.
#'
#' @param x A [StatmodProfile] object.
#' @param label The label.
#' @param column The column of `points` that is the profile.
#'
#' @return A data frame with the columns `value`, `profile`, `converged` and
#'   `drop`, the likelihood ratio statistic \eqn{2(\hat\ell - \ell_p)}, sorted
#'   by `value`. Its attribute `maximum` is \eqn{\hat\ell}.
#'
#' @keywords internal
profile_curve <- function(x, label, column = x@profile_column[[label]]) {
  col <- column
  pts <- x@points
  q <- x@profiled_quantities[[label]]
  if (x@joint) {
    v <- pts[[label]]
    up <- x@unrestricted_point[[col]][1L]
  } else {
    keep <- pts$profiled_label == label
    pts <- pts[keep, , drop = FALSE]
    v <- pts$value
    up <- x@unrestricted_point[[col]][x@unrestricted_point$profiled_label == label][1L]
  }
  pr <- pts[[col]]
  ok <- is.finite(pr)
  d <- data.frame(value = v[ok], profile = pr[ok], converged = pts$converged[ok])
  if (nrow(d)) {
    agg <- stats::aggregate(profile ~ value, data = d, FUN = max)
    cv <- stats::aggregate(converged ~ value, data = d, FUN = all)
    d <- merge(agg, cv, by = "value")
  }
  if (is.finite(up) && !any(abs(d$value - q@estimate) <= 1e-12 * max(1, abs(q@estimate)))) {
    d <- rbind(d, data.frame(value = q@estimate, profile = up, converged = TRUE))
  }
  d <- d[order(d$value), , drop = FALSE]
  mx <- max(c(up, d$profile), na.rm = TRUE)
  d$drop <- pmax(2 * (mx - d$profile), 0)
  rownames(d) <- NULL
  attr(d, "maximum") <- mx
  d
}


#' Profile Likelihood Intervals
#'
#' @description
#' Inverts the profile of each label. Writing \eqn{\hat\ell} for the maximum
#' of the criterion and \eqn{\ell_p(\theta)} for the profile at the value
#' \eqn{\theta}, the interval is the set of values at which the likelihood
#' ratio statistic \eqn{2\{\hat\ell - \ell_p(\theta)\}} does not exceed the
#' quantile of a chi-squared distribution.
#'
#' @details
#' Between two points of the profile the limit is found by a cubic spline
#' through the square root of the statistic, taken with the sign of
#' \eqn{\theta - \hat\theta}, which is close to linear in \eqn{\theta} where
#' the profile is close to quadratic. A limit outside the values of the
#' profile is `NA`, and the function warns.
#'
#' With `type = "marginal"` the quantile has one degree of freedom, and each
#' interval covers its label with probability `level`. With
#' `type = "simultaneous"` it has \eqn{k} degrees of freedom, \eqn{k} being
#' the number of labels of the profile. Each interval is then the projection
#' on its label of the joint region in which the statistic of the \eqn{k}
#' labels does not exceed \eqn{\chi^2_k(\mathrm{level})}, because the profile
#' of one label is the maximum of the criterion over the other labels, and
#' the \eqn{k} intervals cover the labels jointly with probability at least
#' `level`.
#'
#' The profile of the position of a sharp [modelterms7::jump()] or
#' [modelterms7::jseg()] is constant between two consecutive values of the
#' covariate. Its interval runs from the value of the covariate below the
#' first interval of values in which the statistic does not exceed the
#' quantile to the value above the last one.
#'
#' The one-dimensional profile of a label is exact when the labels are
#' profiled one at a time (`joint = FALSE` in [statmod_profile()]). In a
#' joint profile it is the maximum over the grid of the other labels, so its
#' accuracy is limited by the spacing of that grid.
#'
#' The function signals an error when the profile does not give a likelihood
#' ratio statistic (the property `is_likelihood_ratio`): for the hyperparameter of a
#' kinked penalty that a path chooses, and for a coefficient without a
#' penalty held under [reml()], where [ml()] gives the likelihood ratio.
#'
#' @param object A [StatmodProfile] object.
#' @param parm The labels, a subset of `object@profiled_labels`; by default
#'   all of them.
#' @param level The coverage.
#' @param type `"marginal"` or `"simultaneous"`.
#' @param ... Unused.
#'
#' @return A data frame with one row for each label, named by it, and the
#'   columns `estimate`, `lower` and `upper`, on the scale of the profile.
#'
#' @seealso [statmod_profile()], [confint.StatmodFit()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' prof <- statmod_profile(fit, which = "mu:x",
#'                         values = profile_range(width = 3, n = 13))
#' confint(prof)
#' @keywords internal
confint.StatmodProfile <- function(object, parm = NULL, level = 0.95,
                                   type = c("marginal", "simultaneous"), ...) {
  check_dots("confint", ...)
  type <- match.arg(type)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)
  }
  labs <- object@profiled_labels
  if (is.null(parm)) parm <- labs
  bad <- setdiff(parm, labs)
  if (length(bad)) {
    stop(sprintf("'parm' names %s, which the profile does not hold. Its labels are %s.",
                 paste(sprintf("'%s'", bad), collapse = ", "),
                 paste(sprintf("'%s'", labs), collapse = ", ")), call. = FALSE)
  }
  notlr <- parm[!object@is_likelihood_ratio[parm]]
  if (length(notlr)) {
    stop(paste0("The profile of ",
                paste(sprintf("'%s'", notlr), collapse = ", "),
                " does not give a likelihood ratio statistic, so it gives no ",
                "interval.\n  That is the case for the ",
                "hyperparameter of a kinked penalty that a path chooses, and ",
                "for a\n  coefficient without a penalty held under reml(); ",
                "with ml() the drop is a likelihood ratio."),
         call. = FALSE)
  }
  df <- if (type == "marginal") 1 else length(labs)
  cut <- sqrt(stats::qchisq(level, df = df))
  out <- data.frame(estimate = numeric(0), lower = numeric(0), upper = numeric(0))
  outside <- character(0)
  for (lab in parm) {
    q <- object@profiled_quantities[[lab]]
    d <- profile_curve(object, lab)
    if (S7::S7_inherits(q, ProfiledBreakpoint) && isTRUE(q@piecewise_constant)) {
      lim <- profile_step_limits(d, q@covariate_values, cut^2)
      if (anyNA(lim)) outside <- c(outside, lab)
      out[lab, ] <- c(q@estimate, lim)
      next
    }
    r <- sign(d$value - q@estimate) * sqrt(d$drop)
    lim <- c(NA_real_, NA_real_)
    if (nrow(d) >= 2L) {
      f <- stats::splinefun(d$value, r, method = "fmm")
      for (s in c(-1, 1)) {
        side <- if (s < 0) rev(which(d$value <= q@estimate)) else
          which(d$value >= q@estimate)
        hit <- side[s * r[side] >= cut]
        if (!length(hit)) next
        j <- hit[1L]
        prev <- side[match(j, side) - 1L]
        lim[(s + 3) / 2] <- if (is.na(prev)) d$value[j] else {
          a <- d$value[prev]
          b <- d$value[j]
          g <- function(t) s * f(t) - cut
          # the spline can overshoot between two points; the linear
          # interpolation of the signed square root is used where it does
          if (g(a) * g(b) < 0) {
            stats::uniroot(g, sort(c(a, b)), tol = 1e-10)$root
          } else {
            a + (cut - s * r[prev]) / (s * r[j] - s * r[prev]) * (b - a)
          }
        }
      }
    }
    if (anyNA(lim)) outside <- c(outside, lab)
    out[lab, ] <- c(q@estimate, lim)
  }
  if (length(outside)) {
    warning(sprintf(paste0("The profile of %s does not fall to the critical ",
                           "value inside its values, and the limit is NA. ",
                           "Widen the values of statmod_profile()."),
                    paste(sprintf("'%s'", outside), collapse = ", ")),
            call. = FALSE)
  }
  out
}
S7::method(confint, StatmodProfile) <- confint.StatmodProfile


#' Graphs of a Profile
#'
#' @description
#' Draws the profile: the criterion of each point against the value held.
#' The vertical axis is the criterion that the profile maximizes, the
#' log-likelihood, the penalized log-likelihood, or the criterion [reml()] or
#' [ml()] when a hyperparameter is held, and the column of `points` that it
#' reads is the property `profile_column`.
#'
#' @details
#' A profile of one label is drawn as a line through its points, with a
#' dashed vertical line at the estimate and, where the drop of the profile
#' is a likelihood ratio, a dotted horizontal line at the maximum minus
#' \eqn{\chi^2_1(0.95)/2}: the values above it form the interval of
#' [confint.StatmodProfile()] at level 0.95. The graph then depends on the
#' kind of the label:
#' \describe{
#'   \item{a coefficient}{the value on the horizontal axis;}
#'   \item{the position of a break-point}{a rug of the values of the
#'     covariate, at which the profile can have a kink;}
#'   \item{a hyperparameter}{the value on the scale of its link, with the
#'     natural scale on a second axis;}
#'   \item{the hyperparameter of a path}{the criterion that chose it, [bic()]
#'     or [aic()], computed at each point (the column `path_criterion`), with
#'     the points that the path scored as open circles. A hyperparameter
#'     chosen by [cv()] is drawn from the path's record, the cross-validated
#'     criterion with its standard error, which is not recomputed at the
#'     points of the profile.}
#' }
#'
#' Labels profiled one at a time (`joint = FALSE`) give one panel each. A
#' joint profile of two labels is drawn as the contours of the criterion at
#' the maximum minus \eqn{\chi^2_2(\mathrm{level})/2} for the levels 0.5,
#' 0.8, 0.9, 0.95 and 0.99, and a joint profile of three labels or more as one such panel for
#' each pair, where the criterion over the other labels is the maximum over
#' their grid. Where the drop is not a likelihood ratio the contours are
#' drawn at regular levels. Points whose refit did not converge are marked by
#' a cross.
#'
#' @param object,x A [StatmodProfile] object.
#' @param ... Unused.
#'
#' @return `autoplot()` returns a ggplot2 object, or a patchwork object for
#'   several panels; `plot()` returns the same object invisibly, after
#'   printing it.
#'
#' @seealso [statmod_profile()], [confint.StatmodProfile()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' prof <- statmod_profile(fit, which = "mu:x",
#'                         values = profile_range(width = 3, n = 13))
#' ggplot2::autoplot(prof)
#' @name autoplot.StatmodProfile
#' @export
autoplot.StatmodProfile <- function(object, ...) {
  check_dots("autoplot", ...)
  labs <- object@profiled_labels
  k <- length(labs)
  if (!object@joint || k == 1L) {
    ps <- lapply(labs, function(l) profile_panel(object@profiled_quantities[[l]],
                                                 object))
    if (length(ps) == 1L) return(ps[[1L]])
    return(patchwork::wrap_plots(ps))
  }
  prs <- utils::combn(labs, 2L, simplify = FALSE)
  ps <- lapply(prs, function(p) profile_contour(object, p[1L], p[2L]))
  if (length(ps) == 1L) ps[[1L]] else patchwork::wrap_plots(ps)
}

#' @rdname autoplot.StatmodProfile
#' @export
plot.StatmodProfile <- function(x, ...) {
  p <- autoplot.StatmodProfile(x, ...)
  print(p)
  invisible(p)
}

S7::method(autoplot, StatmodProfile) <- autoplot.StatmodProfile
S7::method(plot, StatmodProfile) <- plot.StatmodProfile


#' The Name of the Criterion of a Profile
#'
#' @param x A [StatmodProfile] object.
#' @param label A label of the profile.
#' @param column The column of `points` that is drawn.
#'
#' @return A string for the axis of a graph.
#'
#' @keywords internal
profile_criterion_name <- function(x, label, column = x@profile_column[[label]]) {
  if (identical(column, "path_criterion")) {
    return(x@profiled_quantities[[label]]@criterion_name)
  }
  if (identical(column, "outer_criterion")) {
    return(sprintf("%s criterion", toupper(x@outer_method@kind)))
  }
  if (identical(column, "profile_loglik")) return("log-likelihood")
  pen <- x@points$profile_penalty
  if (!x@joint) pen <- pen[x@points$profiled_label == label]
  if (any(is.finite(pen) & pen != 0)) "penalized log-likelihood" else "log-likelihood"
}


#' The Title and the Axis of a Profiled Quantity
#'
#' @param q A [profiled_quantity] object.
#'
#' @return A list with `title`, the label written as "parameter: name", and
#'   `x`, the name of the quantity on the horizontal axis.
#'
#' @keywords internal
profile_axis <- function(q) {
  lab <- q@label
  if (grepl("/", lab, fixed = TRUE)) {
    parts <- strsplit(lab, "/", fixed = TRUE)[[1L]]
    name <- parts[length(parts)]
    link <- tryCatch(q@link@link_name, error = function(e) "")
    return(list(title = sprintf("%s: %s %s", parts[1L],
                                paste(parts[-c(1L, length(parts))], collapse = "/"),
                                name),
                x = trimws(paste(if (!identical(link, "identity")) link, name))))
  }
  i <- regexpr(":", lab, fixed = TRUE)
  list(title = paste0(substr(lab, 1L, i - 1L), ": ", substr(lab, i + 1L, nchar(lab))),
       x = substr(lab, i + 1L, nchar(lab)))
}


#' The Panel of One Profiled Quantity
#'
#' @description
#' Draws the one-dimensional profile of one label with the graph of its kind,
#' described at [autoplot.StatmodProfile()].
#'
#' @param q A [profiled_quantity] object.
#' @param x The [StatmodProfile] object.
#'
#' @return A ggplot2 object.
#'
#' @keywords internal
profile_panel <- S7::new_generic("profile_panel", "q", function(q, x) S7::S7_dispatch())


#' The Common Layers of a Profile Panel
#'
#' @description
#' Draws the profile of one label against its value, with the estimate, the
#' limit of the interval at level 0.95 where the profile gives a likelihood
#' ratio statistic, and the points whose refit did not converge. Every
#' element is mapped to an aesthetic with a descriptive label, so that the
#' legend names it; [profile_finish()] adds the legend.
#'
#' @param q A [profiled_quantity] object.
#' @param x The [StatmodProfile] object.
#' @param column The column of `points` that is drawn.
#' @param cutoff Whether to draw the line of the interval at level 0.95.
#' @param steps Whether the profile is constant between two consecutive
#'   values of the covariate, and is drawn as one segment over each interval.
#' @param line Whether the points of the profile are joined by a line.
#'
#' @return A list with `plot`, a ggplot2 object, `data`, the data frame
#'   drawn, and `entries`, the legend entries of the colour scale, a list of
#'   lists with `label`, `colour`, `shape` and `linetype`.
#'
#' @keywords internal
profile_base <- function(q, x, column = x@profile_column[[q@label]],
                         cutoff = x@is_likelihood_ratio[[q@label]],
                         steps = FALSE, line = TRUE) {
  d <- profile_curve(x, q@label, column = column)
  ax <- profile_axis(q)
  crit <- profile_criterion_name(x, q@label, column)
  lab <- sprintf("%s at each value held", crit)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$profile)) +
    ggplot2::theme_bw()
  if (cutoff) {
    p <- p + ggplot2::geom_hline(
      ggplot2::aes(yintercept = attr(d, "maximum") - stats::qchisq(0.95, 1) / 2,
                   linetype = "limit of the 95% interval"),
      colour = "#b2182b")
  }
  p <- p + ggplot2::geom_vline(ggplot2::aes(xintercept = q@estimate,
                                            linetype = "estimate"),
                               colour = "#b2182b", linewidth = 0.3)
  if (steps) {
    # constant between two consecutive values of the covariate: one segment
    # over the interval that holds each value of the grid
    u <- sort(unique(q@covariate_values))
    j <- findInterval(d$value, u)
    ok <- j >= 1L & j < length(u)
    seg <- data.frame(x0 = u[j[ok]], x1 = u[j[ok] + 1L], y = d$profile[ok])
    lab <- sprintf("%s, constant between observations", crit)
    p <- p + ggplot2::geom_segment(
      data = seg, ggplot2::aes(x = .data$x0, xend = .data$x1, y = .data$y,
                               yend = .data$y, colour = lab),
      inherit.aes = FALSE)
    entries <- list(list(label = lab, colour = "#1f5fa8", shape = NA,
                         linetype = "solid"))
  } else {
    if (line) p <- p + ggplot2::geom_line(ggplot2::aes(colour = lab))
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = lab), size = 1)
    entries <- list(list(label = lab, colour = "#1f5fa8", shape = 16,
                         linetype = if (line) "solid" else "blank"))
  }
  bad <- d[!d$converged, , drop = FALSE]
  if (nrow(bad)) {
    p <- p + ggplot2::geom_point(data = bad,
                                 ggplot2::aes(shape = "refit did not converge"),
                                 size = 3, colour = "#b2182b")
  }
  p <- p + ggplot2::labs(title = ax$title, x = ax$x, y = crit)
  list(plot = p, data = d, entries = entries)
}


#' The Legend of a Profile Panel
#'
#' @description
#' Adds the scales that name every element of a panel: the colour scale for
#' the drawn series, with the key of each entry drawn as its line or its
#' point, the line type of the estimate and of the limit of the interval,
#' and the shape of the points whose refit did not converge.
#'
#' @param p A ggplot2 object.
#' @param entries The entries of the colour scale, as [profile_base()]
#'   returns them.
#'
#' @return `p` with the scales and the legend at the bottom.
#'
#' @keywords internal
profile_finish <- function(p, entries) {
  labs <- vapply(entries, function(e) e$label, "")
  cols <- vapply(entries, function(e) e$colour, "")
  shp <- vapply(entries, function(e) as.numeric(e$shape), numeric(1))
  lty <- vapply(entries, function(e) e$linetype, "")
  p + ggplot2::scale_colour_manual(
        values = stats::setNames(cols, labs), breaks = labs, name = NULL,
        guide = ggplot2::guide_legend(
          order = 1, override.aes = list(shape = shp, linetype = lty))) +
    ggplot2::scale_linetype_manual(
      values = c("estimate" = "dashed", "limit of the 95% interval" = "dotted",
                 "plus and minus one standard error" = "solid"),
      name = NULL, guide = ggplot2::guide_legend(order = 2)) +
    (if (any(vapply(p$layers, function(l) !is.null(l$mapping$shape), TRUE)))
      ggplot2::scale_shape_manual(values = c("refit did not converge" = 4),
                                  name = NULL,
                                  guide = ggplot2::guide_legend(order = 3))) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                   legend.margin = ggplot2::margin(0, 0, 0, 0))
}

S7::method(profile_panel, profiled_quantity) <- function(q, x) {
  b <- profile_base(q, x)
  profile_finish(b$plot, b$entries)
}

S7::method(profile_panel, ProfiledBreakpoint) <- function(q, x) {
  b <- profile_base(q, x, steps = isTRUE(q@piecewise_constant))
  d <- b$data
  rug <- data.frame(value = q@covariate_values)
  rug <- rug[rug$value >= min(d$value) & rug$value <= max(d$value), , drop = FALSE]
  lab <- "observed values of the covariate"
  p <- b$plot + ggplot2::geom_rug(data = rug,
                                  ggplot2::aes(x = .data$value, colour = lab),
                                  inherit.aes = FALSE, sides = "b", alpha = 0.5)
  profile_finish(p, c(b$entries, list(list(label = lab, colour = "grey30",
                                           shape = NA, linetype = "solid"))))
}

S7::method(profile_panel, ProfiledHyperparameter) <- function(q, x) {
  link <- q@link
  name <- sub(".*/", "", q@label)
  b <- profile_base(q, x)
  p <- b$plot +
    ggplot2::scale_x_continuous(
      sec.axis = ggplot2::sec_axis(
        transform = function(t) linkfunctions7::linkinv(link, t),
        name = sprintf("%s (natural scale)", name)))
  profile_finish(p, b$entries)
}

S7::method(profile_panel, ProfiledPathHyperparameter) <- function(q, x) {
  rec <- q@path_record
  col <- x@profile_column[[q@label]]
  ax <- profile_axis(q)
  if (identical(col, "path_criterion")) {
    b <- profile_base(q, x, cutoff = FALSE, line = FALSE)
    entries <- b$entries
    p <- b$plot
    if (nrow(rec)) {
      lab <- sprintf("%s at the values of the path", q@criterion_name)
      p <- p + ggplot2::geom_point(data = rec,
                                   ggplot2::aes(x = .data$value, y = .data$criterion,
                                                colour = lab),
                                   inherit.aes = FALSE, shape = 1)
      entries <- c(entries, list(list(label = lab, colour = "grey30", shape = 1,
                                      linetype = "blank")))
    }
    return(profile_finish(p, entries))
  }
  if (nrow(rec) && any(is.finite(rec$se))) {
    # a cross-validated criterion is not recomputed at the points of the
    # profile: the graph is the path's own record, with its standard errors
    lab <- "cross-validated criterion at the values of the path"
    p <- ggplot2::ggplot(rec, ggplot2::aes(x = .data$value, y = .data$criterion)) +
      ggplot2::theme_bw() +
      ggplot2::geom_vline(ggplot2::aes(xintercept = q@estimate,
                                       linetype = "estimate"),
                          colour = "#b2182b", linewidth = 0.3) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$criterion - .data$se,
                                          ymax = .data$criterion + .data$se,
                                          linetype = "plus and minus one standard error"),
                             colour = "grey60", width = 0) +
      ggplot2::geom_line(ggplot2::aes(colour = lab)) +
      ggplot2::geom_point(ggplot2::aes(colour = lab), size = 1) +
      ggplot2::labs(title = ax$title, x = ax$x, y = "cross-validated criterion")
    return(profile_finish(p, list(list(label = lab, colour = "#1f5fa8",
                                       shape = 16, linetype = "solid"))))
  }
  b <- profile_base(q, x, column = "profile_loglik", cutoff = FALSE)
  profile_finish(b$plot, b$entries)
}


#' The Contours of a Joint Profile Over Two Labels
#'
#' @description
#' Draws the contours of a joint profile over two of its labels, the profile
#' over any other label being the maximum over its grid.
#'
#' @param x A [StatmodProfile] object with `joint = TRUE`.
#' @param a,b The two labels.
#'
#' @return A ggplot2 object.
#'
#' @keywords internal
profile_contour <- function(x, a, b) {
  col <- x@profile_column[[a]]
  pts <- x@points
  d <- data.frame(a = pts[[a]], b = pts[[b]], profile = pts[[col]],
                  converged = pts$converged)
  d <- d[is.finite(d$profile), , drop = FALSE]
  agg <- stats::aggregate(profile ~ a + b, data = d, FUN = max)
  cv <- stats::aggregate(converged ~ a + b, data = d, FUN = all)
  d <- merge(agg, cv, by = c("a", "b"))
  mx <- max(c(x@unrestricted_point[[col]][1L], d$profile), na.rm = TRUE)
  bad <- d[!d$converged, , drop = FALSE]
  d <- profile_refine(d)
  qa <- x@profiled_quantities[[a]]
  qb <- x@profiled_quantities[[b]]
  crit <- profile_criterion_name(x, a)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$a, y = .data$b)) +
    ggplot2::theme_bw()
  if (x@is_likelihood_ratio[[a]]) {
    cover <- c(0.99, 0.95, 0.9, 0.8, 0.5)
    brk <- mx - stats::qchisq(cover, 2) / 2
    p <- p + ggplot2::geom_contour(
      ggplot2::aes(z = .data$profile, colour = ggplot2::after_stat(factor(.data$level))),
      breaks = brk) +
      ggplot2::scale_colour_manual(
        values = grDevices::colorRampPalette(c("#0b2f5e", "#9ecae1"))(length(cover)),
        labels = format(cover), name = "joint confidence region, level")
  } else {
    p <- p + ggplot2::geom_contour(
      ggplot2::aes(z = .data$profile, colour = sprintf("contours of the %s", crit)),
      bins = 12) +
      ggplot2::scale_colour_manual(values = "#1f5fa8", name = NULL)
  }
  p <- p +
    ggplot2::geom_point(data = data.frame(a = qa@estimate, b = qb@estimate),
                        ggplot2::aes(shape = "estimate"), size = 3,
                        colour = "#b2182b")
  if (nrow(bad)) {
    p <- p + ggplot2::geom_point(data = bad,
                                 ggplot2::aes(shape = "refit did not converge"),
                                 size = 3, colour = "#b2182b")
  }
  p + ggplot2::scale_shape_manual(
        values = c("estimate" = 3, "refit did not converge" = 4), name = NULL) +
    ggplot2::labs(x = profile_axis(qa)$title, y = profile_axis(qb)$title,
                  title = sprintf("profile of the %s", crit)) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
}


#' A Profile Over Two Labels on a Finer Grid
#'
#' @description
#' Interpolates the criterion of a profile over two labels, given on a
#' complete rectangular grid, to a grid of 60 by 60 values by cubic splines
#' along each axis in turn, so that its contours are smooth curves and not
#' the polygons of the coarse grid. A grid that is not complete, or that has
#' fewer than four values of a label, is returned as it is.
#'
#' @param d A data frame with the columns `a`, `b` and `profile`.
#'
#' @return A data frame with the same columns.
#'
#' @keywords internal
profile_refine <- function(d) {
  ua <- sort(unique(d$a))
  ub <- sort(unique(d$b))
  if (length(ua) < 4L || length(ub) < 4L || nrow(d) != length(ua) * length(ub)) {
    return(d[, c("a", "b", "profile")])
  }
  z <- matrix(NA_real_, length(ua), length(ub))
  z[cbind(match(d$a, ua), match(d$b, ub))] <- d$profile
  ga <- seq(ua[1L], ua[length(ua)], length.out = 60L)
  gb <- seq(ub[1L], ub[length(ub)], length.out = 60L)
  za <- apply(z, 2L, function(col) stats::spline(ua, col, xout = ga, method = "natural")$y)
  zz <- t(apply(za, 1L, function(row) stats::spline(ub, row, xout = gb, method = "natural")$y))
  data.frame(a = rep(ga, times = length(gb)), b = rep(gb, each = length(ga)),
             profile = as.vector(zz))
}


#' The Interval of a Profile That Is Constant Between Observations
#'
#' @param d The profile, as [profile_curve()] returns it.
#' @param covariate The values of the covariate.
#' @param cut The quantile that the statistic must not exceed.
#'
#' @return The two limits: the value of the covariate below the first
#'   interval in which the statistic does not exceed `cut`, and the value
#'   above the last one. A limit is `NA` where that interval is the first or
#'   the last of the profile, the region possibly extending beyond it.
#'
#' @keywords internal
profile_step_limits <- function(d, covariate, cut) {
  u <- sort(unique(covariate))
  acc <- which(d$drop <= cut)
  if (!length(acc)) return(c(NA_real_, NA_real_))
  lo <- if (min(acc) == 1L) NA_real_ else {
    j <- findInterval(d$value[min(acc)], u)
    if (j >= 1L) u[j] else NA_real_
  }
  hi <- if (max(acc) == nrow(d)) NA_real_ else {
    j <- findInterval(d$value[max(acc)], u)
    if (j < length(u)) u[j + 1L] else NA_real_
  }
  c(lo, hi)
}


#' The Points That a Path Scored
#'
#' @param fit A [StatmodFit()].
#' @param r The row of the hyperparameter in the table of
#'   [profile_hyper_table()].
#' @param link The link of the hyperparameter.
#'
#' @return A data frame with the columns `value`, on the scale of the link,
#'   `criterion` and `se`, sorted by `value`; it has no rows where the fit
#'   kept no record of the path.
#'
#' @keywords internal
profile_path_record <- function(fit, r, link) {
  h <- fit@history$outer
  empty <- data.frame(value = numeric(0), criterion = numeric(0), se = numeric(0))
  if (is.null(h) || !all(c("parameter", "term", "name", "value", "criterion") %in% names(h))) {
    return(empty)
  }
  h <- h[h$parameter == r$parameter & h$term == r$term & h$name == r$name &
           is.finite(h$criterion), , drop = FALSE]
  if (!nrow(h)) return(empty)
  out <- data.frame(value = linkfunctions7::linkfun(link, h$value),
                    criterion = h$criterion,
                    se = if (is.null(h$se)) NA_real_ else h$se)
  out[order(out$value), , drop = FALSE]
}
