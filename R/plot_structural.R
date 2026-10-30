#' @include plot_core.R plot_effect.R
NULL

# Graphs of the structural terms: the break-point terms and nl(), whose block
# moves with their coefficients, and gas() and regime(), which carry a state.
# The core's plot_effect() method reads a contribution as the block times the
# coefficients, which is wrong for a block that is a Jacobian; the methods
# here read modelterms7::term_value() for the contribution and the block for
# its derivative.


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Standard Errors of a Linear Map of Coefficients with Missing Variances
#'
#' @description
#' Returns \eqn{\sqrt{x_i' V x_i}} for each row \eqn{x_i} of `X`. A coefficient
#' with no variance (an entry of `V` that is `NA`) is treated as absent in the
#' rows that do not read it, and the standard error is `NA` in the rows that do.
#'
#' @param X A matrix, one row per evaluation.
#' @param V The variance matrix of the coefficients that the columns of `X`
#'   multiply.
#'
#' @return A numeric vector with one value per row of `X`.
#'
#' @keywords internal
moving_se <- function(X, V) {
  X <- as.matrix(X)
  na <- is.na(diag(V))
  V0 <- V
  V0[is.na(V0)] <- 0
  se <- sqrt(pmax(rowSums((X %*% V0) * X), 0))
  if (any(na)) {
    se[rowSums(abs(X[, na, drop = FALSE]) > 0) > 0] <- NA_real_
  }
  se
}


#' The Variables That a Term with a Moving Block Reads
#'
#' @description
#' Returns the data variables of a break-point term or of a nonlinear term:
#' the variable of the break-points and those of the subformulas for a
#' break-point term, and the data variables recorded when the term was built
#' for a nonlinear one.
#'
#' @param tm A [modelterms7::SegTerm] or [modelterms7::NlTerm].
#' @param name The name of the term in its equation.
#'
#' @return A character vector.
#'
#' @keywords internal
moving_term_vars <- function(tm, name) {
  if (S7::S7_inherits(tm, modelterms7::NlTerm)) {
    return(names(tm@build_info$data_vars))
  }
  unique(c(all.vars(tm@var), term_variables(tm, name)))
}


#' The Working Weights of an Equation
#'
#' @description
#' Returns the information of the log-likelihood in the predictor of an
#' equation, \eqn{w_i}, expected or observed according to `information`, and
#' the working residuals \eqn{s_i / w_i} of [working_residuals()], from one
#' evaluation of the derivatives.
#'
#' @param fit A [StatmodFit()].
#' @param param The equation.
#' @param information An object from [distributions7::expected_info()] or
#'   [distributions7::observed_info()].
#'
#' @return A list with `wr` and `w`, numeric vectors with one value per
#'   observation; `wr` is `NA` and `w` is zero where the information is not
#'   positive.
#'
#' @keywords internal
working_parts <- function(fit, param, information) {
  wr <- working_residuals(fit, param, information)
  spec <- fit@spec
  th <- fit@fitted
  key <- paste(param, param, sep = "_")
  h <- if (S7::S7_inherits(information, distributions7::ExpectedInfo)) {
    distributions7::distrib_expected_hessian(spec@distrib, spec@response, th,
                                             scale = "link",
                                             method = information@method)
  } else {
    distributions7::distrib_hessian(spec@distrib, spec@response, th,
                                    scale = "link")
  }
  w <- -h[[key]]
  w[!is.finite(w) | w <= 0 | !is.finite(wr)] <- 0
  list(wr = wr, w = w)
}


#' The Data of the Effect of a Term with a Moving Block
#'
#' @description
#' The counterpart of [effect_data()] for [modelterms7::SegTerm] and
#' [modelterms7::NlTerm]. The contribution on the grid is
#' [modelterms7::term_value()], not the block times the coefficients, and the
#' band is the delta method with the block of [modelterms7::term_predict()],
#' which is the derivative of the contribution in the coefficients.
#'
#' @param type A `PlotEffect` object.
#' @param tm The term.
#' @param fit A [StatmodFit()].
#' @param ctx The context described at [plot_data()].
#'
#' @return As [effect_data()].
#'
#' @keywords internal
moving_effect_data <- function(type, tm, fit, ctx) {
  vars_all <- moving_term_vars(tm, ctx$name)
  data <- plot_needs_data(ctx, fit, type, vars = vars_all)
  plotted <- intersect(ctx$vars, vars_all)
  if (!length(plotted)) plotted <- vars_all
  usable <- Filter(function(v) {
    vv <- data[[v]]
    !is.matrix(vv) && (is.numeric(vv) || is.factor(vv) || is.character(vv))
  }, plotted)
  if (!length(usable)) {
    stop(sprintf(paste0("the term '%s' reads no numeric or factor variable ",
                        "to put on an axis."), ctx$name), call. = FALSE)
  }
  num <- Filter(function(v) is.numeric(data[[v]]), usable)
  if (length(num) > 1L) {
    stop(sprintf(paste0("the term '%s' reads the numeric variables %s; ",
                        "plot_surface() draws a function of two."), ctx$name,
                 paste(num, collapse = " and ")), call. = FALSE)
  }
  var <- if (length(num)) num[[1L]] else usable[[1L]]
  kind <- if (is.numeric(data[[var]])) "numeric" else "factor"
  grid <- effect_grid(data, var, vars_all, type@n_grid)
  group <- attr(grid, "group")
  z <- stats::qnorm((1 + type@level) / 2)
  empty <- data.frame(x = numeric(0), y = numeric(0))
  points <- bins <- empty

  if (type@scale == "response") {
    nd <- response_newdata(fit, data, grid, type@at)
    pr <- predict(fit, what = ctx$param, newdata = nd, se = TRUE,
                  level = type@level, random = "zero",
                  vcov_type = type@vcov_type)
    curve <- data.frame(x = grid[[var]], fit = pr$fit, lower = pr$lower,
                        upper = pr$upper)
    if (!is.null(group)) curve$group <- grid[[group]]
    ylab <- ctx$param
  } else {
    cb <- term_coef_block(fit, ctx$param, ctx$name, NULL, type@vcov_type)
    val <- modelterms7::term_value(tm, newdata = grid)
    se <- moving_se(modelterms7::term_predict(tm, newdata = grid), cb$vcov)
    curve <- data.frame(x = grid[[var]], fit = val, lower = val - z * se,
                        upper = val + z * se)
    if (!is.null(group)) curve$group <- grid[[group]]
    if (type@residuals || !is.null(type@binned)) {
      ofit <- modelterms7::term_value(tm)
      wr <- working_residuals(fit, ctx$param, type@information)
      keep <- fit@spec@weights > 0 & is.finite(wr)
      pts <- data.frame(x = data[[var]], y = ofit + wr)[keep, ]
      if (!is.null(group)) pts$group <- as.factor(data[[group]])[keep]
      if (type@residuals) points <- pts
      if (!is.null(type@binned) && kind == "numeric") {
        br <- unique(stats::quantile(pts$x, probs = seq(0, 1,
                                     length.out = type@binned + 1)))
        cl <- cut(pts$x, breaks = br, include.lowest = TRUE)
        bins <- data.frame(x = as.numeric(tapply(pts$x, cl, mean)),
                           y = as.numeric(tapply(pts$y, cl, mean)))
        bins <- bins[stats::complete.cases(bins), , drop = FALSE]
      }
    }
    ylab <- sprintf("effect on the link scale of %s", ctx$param)
  }
  rug <- if (type@rug && kind == "numeric") data.frame(x = data[[var]])
    else data.frame(x = numeric(0))
  out <- list(curve = curve, points = points, bins = bins, rug = rug)
  attr(out, "kind") <- kind
  attr(out, "labels") <- list(title = sprintf("%s: %s", ctx$param,
                                              if (is.null(ctx$title)) ctx$name
                                              else ctx$title),
                              x = var, y = ylab, colour = group, fill = group)
  out
}


S7::method(plot_data, list(PlotEffect, modelterms7::SegTerm)) <-
  function(type, x, fit, ctx, ...) moving_effect_data(type, x, fit, ctx)

S7::method(plot_data, list(PlotEffect, modelterms7::NlTerm)) <-
  function(type, x, fit, ctx, ...) moving_effect_data(type, x, fit, ctx)


# ---------------------------------------------------------------------------
# plot_breakpoint()
# ---------------------------------------------------------------------------

#' @title S7 Class of the Break-Point Graph
#'
#' @description
#' The class that [plot_breakpoint()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,legend,test,level,vcov_type,information,n_grid,n_profile
#'   See [plot_breakpoint()].
#'
#' @return A `PlotBreakpoint` object, inheriting from [term_plot].
#'
#' @name PlotBreakpoint
#' @keywords internal
PlotBreakpoint <- S7::new_class("PlotBreakpoint", parent = term_plot,
  properties = list(
    test = S7::class_any,
    level = S7::class_numeric,
    vcov_type = S7::class_any,
    information = S7::class_any,
    n_grid = S7::class_numeric,
    n_profile = S7::class_numeric
  ))


#' Graph of the Break-Points of a Term
#'
#' @description
#' Selects the graph of a break-point term ([modelterms7::seg()],
#' [modelterms7::jump()] or [modelterms7::jseg()]): the contribution of the
#' term as a function of its variable, with the partial residuals, each
#' break-point as a dashed vertical line, and its confidence interval as a bar
#' at the bottom of the graph.
#'
#' @details
#' The contribution and its band are those of [plot_effect()] for a term whose
#' block moves with its coefficients: \eqn{f(x) = }
#' [modelterms7::term_value()] on a grid of the variable, and the band is
#' \eqn{f(x) \pm z_{(1+\text{level})/2}\,\sqrt{J V J'}}, where \eqn{J} is the
#' block of [modelterms7::term_predict()] (the derivative of \eqn{f} in the
#' coefficients) and \eqn{V} the variance of the coefficients of the term.
#'
#' The interval of a break-point of a [modelterms7::seg()] term is
#' `confint(fit, parm, test = test, readable = FALSE)` on its coefficient, so
#' that `test = lr_test()` inverts the likelihood ratio test of the position
#' and `test = wald_test()` gives the Wald interval. A limit that the data do
#' not bound is infinite, and the bar is drawn to the end of the range of the
#' variable. The position of a
#' [modelterms7::jump()] or [modelterms7::jseg()] term is read off other
#' coefficients and has no coordinate to invert a test on; its interval is
#' the set of positions, connected to the estimate, at which the profile
#' described below stays under the cutoff \eqn{\chi^2_{1}(\text{level})}.
#'
#' The profile is computed on the working response of the equation,
#' \eqn{y_i^* = f(x_i) + s_i / w_i}, where \eqn{s_i} is the score of the
#' log-likelihood in the predictor and \eqn{w_i} the information selected by
#' `information`. For a position \eqn{\psi_k} on a grid, with the other
#' positions and the other terms of the equation at their fitted values,
#' \eqn{R(\psi_k)} is the weighted residual sum of squares of \eqn{y^*} on an
#' intercept, the linear column and the columns of the term, as computed by
#' [modelterms7::seg_profile_rss()], and \eqn{R(\psi_k) - \min R} is the
#' increase of the criterion in the units of the log-likelihood. For a
#' gaussian mean it is the likelihood ratio statistic; for other families it
#' is its quadratic approximation at the fit.
#'
#' A term whose break-points are developed over covariates has one position per
#' observation and no single profile, and a marginal break-point term has a
#' distribution of positions per group; both signal an error.
#'
#' @param test A test object ([wald_test()], [lr_test()] and others), the test
#'   inverted for the interval of a break-point.
#' @param level The confidence level of the band and of the intervals.
#' @param vcov_type A [vcov_method] object, the variance of the band.
#' @param information An object from [distributions7::expected_info()] (the
#'   default) or [distributions7::observed_info()], the weight of the working
#'   response.
#' @param n_grid The number of points of the grid of the variable.
#' @param n_profile The number of positions at which the profile of a
#'   [modelterms7::jump()] or [modelterms7::jseg()] position is computed.
#'
#' @param legend Whether to draw the legend that names the elements of
#'   the graph.
#'
#' @return A `PlotBreakpoint` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_effect()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(120, 0, 10)))
#' dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 6, 0) + rnorm(120, sd = 0.4)
#' fit <- statmod(y ~ seg(x, npsi = 1), distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, term = "seg(x, npsi = 1)", data = dd)
#' ggplot2::autoplot(fit, term = "seg(x, npsi = 1)", data = dd,
#'                   plot_type = plot_breakpoint(test = wald_test()))
#' @export
plot_breakpoint <- function(test = lr_test(), level = 0.95,
                            vcov_type = bayesian_vcov(),
                            information = distributions7::expected_info(),
                            n_grid = 200, n_profile = 100, legend = TRUE) {
  test <- as_test_statistic(test)
  check_level_grid(level, n_grid)
  if (!is.numeric(n_profile) || length(n_profile) != 1L ||
      !is.finite(n_profile) || n_profile < 5 || n_profile != round(n_profile)) {
    stop("'n_profile' must be a single whole number of at least 5.",
         call. = FALSE)
  }
  vcov_type <- as_vcov_method(vcov_type)
  if (!S7::S7_inherits(information, distributions7::ExpectedInfo) &&
      !S7::S7_inherits(information, distributions7::ObservedInfo)) {
    stop("'information' must be expected_info() or observed_info().",
         call. = FALSE)
  }
  PlotBreakpoint(label = "plot_breakpoint()", test = test,
                 level = level, vcov_type = vcov_type,
                 information = information, n_grid = n_grid,
                 n_profile = n_profile, legend = check_legend(legend))
}


#' The Interval of a Position from Its Profile
#'
#' @description
#' Returns the smallest and largest grid positions of the run, connected to
#' `psi`, over which the increase of the criterion stays at or under
#' `cutoff`.
#'
#' @param x The grid of positions, increasing.
#' @param delta The increase of the criterion at each position.
#' @param psi The estimated position.
#' @param cutoff The cutoff of the increase.
#'
#' @return A numeric vector of length two; `NA` where the profile is empty.
#'
#' @keywords internal
profile_interval <- function(x, delta, psi, cutoff) {
  if (!length(x)) return(c(NA_real_, NA_real_))
  i0 <- which.min(abs(x - psi))
  inside <- delta <= cutoff
  lo <- i0
  while (lo > 1L && inside[lo - 1L]) lo <- lo - 1L
  hi <- i0
  while (hi < length(x) && inside[hi + 1L]) hi <- hi + 1L
  c(x[lo], x[hi])
}


S7::method(plot_data, list(PlotBreakpoint, modelterms7::SegTerm)) <-
  function(type, x, fit, ctx, ...) {
    if (length(x@subformulas)) {
      stop(sprintf(paste0("the break-points of '%s' are developed over ",
                          "covariates and have one position per observation; ",
                          "plot_effect() draws the term."), ctx$name),
           call. = FALSE)
    }
    var <- all.vars(x@var)[1L]
    plot_needs_data(ctx, fit, type, vars = var)
    pe <- plot_effect(level = type@level, vcov_type = type@vcov_type,
                      information = type@information, n_grid = type@n_grid,
                      rug = FALSE)
    eff <- moving_effect_data(pe, x, fit, ctx)
    psi <- as.numeric(modelterms7::seg_psi(x))
    K <- length(psi)

    # the interval of each break-point, from the test the caller selected
    lower <- upper <- rep(NA_real_, K)
    des <- statmod_design(fit@spec)[[ctx$param]]
    cn <- des$coef_names[des$blocks[[ctx$name]]]
    if (identical(x@kind, "seg")) {
      pn <- paste0(ctx$param, ":", cn[grepl("psi[0-9]+$", cn)])
      ci <- tryCatch(
        confint(fit, parm = pn, level = type@level,
                vcov_type = type@vcov_type, readable = FALSE,
                test = type@test),
        error = function(e) {
          warning(sprintf("the intervals of the break-points could not be computed: %s",
                          conditionMessage(e)), call. = FALSE)
          NULL
        })
      if (!is.null(ci)) {
        lower <- ci[pn, "lower"]
        upper <- ci[pn, "upper"]
      }
    }

    # the bars sit under the lowest of the data and the band
    ys <- c(eff$points$y, eff$curve$lower, eff$curve$upper)
    ys <- ys[is.finite(ys)]
    yr <- range(ys)
    bar_y <- yr[1L] - 0.04 * diff(yr)
    breaks <- data.frame(k = seq_len(K), psi = psi, lower = lower,
                         upper = upper, y = bar_y)

    profile <- data.frame(x = numeric(0), delta = numeric(0),
                          k = integer(0))
    cutoff <- data.frame(cutoff = stats::qchisq(type@level, 1))
    # a jump() or jseg() position has no coefficient of its own, so its
    # interval is read from the profile
    from_profile <- !identical(x@kind, "seg")
    if (from_profile) {
      wp <- working_parts(fit, ctx$param, type@information)
      yn <- modelterms7::term_value(x) + wp$wr
      w <- fit@spec@weights * wp$w
      yn[!is.finite(yn)] <- 0
      lim <- x@build_info$lim
      pieces <- lapply(seq_len(K), function(k) {
        g <- sort(unique(c(seq(lim[1L], lim[2L], length.out = type@n_profile),
                           psi[k])))
        rss <- vapply(g, function(v) {
          p <- psi
          p[k] <- v
          tryCatch(modelterms7::seg_profile_rss(
            modelterms7::seg_relocate(x, p), yn, weights = w),
            error = function(e) Inf)
        }, numeric(1))
        ok <- is.finite(rss)
        data.frame(x = g[ok], delta = rss[ok] - min(rss[ok]), k = k)
      })
      profile <- do.call(rbind, pieces)
    }
    if (from_profile && nrow(profile)) {
      for (k in seq_len(K)) {
        pk <- profile[profile$k == k, , drop = FALSE]
        iv <- profile_interval(pk$x, pk$delta, psi[k], cutoff$cutoff)
        breaks$lower[k] <- iv[1L]
        breaks$upper[k] <- iv[2L]
      }
    }
    # a limit the data do not bound is infinite, and the bar is drawn to the
    # edge of the variable's range
    rng <- range(x@build_info$xv)
    breaks$from <- pmax(breaks$lower, rng[1L])
    breaks$to <- pmin(breaks$upper, rng[2L])
    out <- list(curve = eff$curve, points = eff$points, breaks = breaks,
                profile = profile, cutoff = cutoff)
    attr(out, "labels") <- list(title = sprintf("%s: %s", ctx$param,
                                                if (is.null(ctx$title))
                                                  ctx$name else ctx$title),
                                x = var,
                                y = sprintf("effect on the link scale of %s",
                                            ctx$param))
    out
  }

S7::method(plot_data, list(PlotBreakpoint, modelterms7::MarginalBreakTerm)) <-
  function(type, x, fit, ctx, ...) {
    stop(sprintf(paste0("the marginal break-point '%s' has a distribution of ",
                        "positions over the groups and no single position to ",
                        "draw; statmod_latent() gives the position of each ",
                        "group."), ctx$name), call. = FALSE)
  }


S7::method(plot_layers, PlotBreakpoint) <- function(type, data, ...) {
  k_pts <- "partial residuals"
  k_fit <- est_label("estimated effect", type@level)
  k_psi <- est_label("estimated break-point", type@level)
  out <- list()
  if (nrow(data$points)) {
    out <- c(out, list(ggplot2::geom_point(
      data = data$points,
      mapping = ggplot2::aes(x = .data$x, y = .data$y, colour = k_pts),
      alpha = 0.35, size = 0.9)))
  }
  out <- c(out, list(
    ggplot2::geom_ribbon(data = data$curve,
                         mapping = ggplot2::aes(x = .data$x,
                                                ymin = .data$lower,
                                                ymax = .data$upper),
                         fill = "#1f5fa8", alpha = 0.18),
    ggplot2::geom_line(data = data$curve,
                       mapping = ggplot2::aes(x = .data$x, y = .data$fit,
                                              colour = k_fit),
                       linewidth = 0.9, key_glyph = draw_key_band),
    ggplot2::geom_vline(data = data$breaks,
                        mapping = ggplot2::aes(xintercept = .data$psi,
                                               colour = k_psi),
                        linetype = "dashed")))
  has <- is.finite(data$breaks$from) & is.finite(data$breaks$to)
  if (any(has)) {
    out <- c(out, list(ggplot2::geom_segment(
      data = data$breaks[has, , drop = FALSE],
      mapping = ggplot2::aes(x = .data$from, xend = .data$to,
                             y = .data$y, yend = .data$y, colour = k_psi),
      linewidth = 1.6)))
  }
  c(out, plot_key(colour = stats::setNames(c("#1f5fa8", "grey30", "#b2182b"),
                                           c(k_fit, k_pts, k_psi))))
}


# ---------------------------------------------------------------------------
# plot_filter()
# ---------------------------------------------------------------------------

#' @title S7 Class of the Filter Graph
#'
#' @description
#' The class that [plot_filter()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,level,scale,response,groups,vcov_type See [plot_filter()].
#'
#' @return A `PlotFilter` object, inheriting from [term_plot].
#'
#' @name PlotFilter
#' @keywords internal
PlotFilter <- S7::new_class("PlotFilter", parent = term_plot,
  properties = list(
    level = S7::class_numeric,
    scale = S7::class_character,
    response = S7::class_logical,
    groups = S7::class_any,
    vcov_type = S7::class_any
  ))


#' Graph of a Filtered Parameter
#'
#' @description
#' Selects the graph of the path that a score-driven term ([modelterms7::gas()])
#' gives to the parameter of its equation: the filtered parameter over the
#' order of the series, with a pointwise confidence band, and optionally the
#' observed response beside it. A model in which several equations carry a
#' `gas()` term joins them in one filter, and the graph then has a panel for
#' each of the filtered equations.
#'
#' @details
#' The path is the fitted parameter \eqn{\theta_t = g^{-1}(\eta_t)}, where
#' \eqn{\eta_t} is the static predictor plus the level of the filter, and the
#' band is the interval of [predict.StatmodFit()] with `se = TRUE`, which
#' differentiates the filter in the coefficients and in the parameters of the
#' term (the standard error of a filtered predictor). The horizontal axis is the
#' variable named by `time` in the term, or the row order when it has none.
#'
#' For a panel (`gas(by = id)`) the graph has one panel for each group, at most
#' twelve; with more groups, `groups` names the ones to draw.
#'
#' @param level The confidence level of the band.
#' @param scale `"response"` (the default), the parameter, or `"link"`, the
#'   predictor.
#' @param response Whether to draw the observed response. It is drawn on the
#'   scale of the parameter, so it requires `scale = "response"`, a response
#'   that is a vector, and a term in one equation.
#' @param groups `NULL` for every group, or the labels of the groups to draw
#'   (the levels of the `by` variable as text).
#' @param vcov_type A [vcov_method] object, the variance of the band.
#'
#' @param legend Whether to draw the legend that names the elements of
#'   the graph.
#'
#' @return A `PlotFilter` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_news_impact()]
#'
#' @examples
#' set.seed(4)
#' n <- 100
#' x <- rnorm(n)
#' f <- 0
#' y <- integer(n)
#' for (t in seq_len(n)) {
#'   lam <- exp(0.7 + 0.3 * x[t] + f)
#'   y[t] <- rpois(1, lam)
#'   f <- 0.15 * (y[t] - lam) * lam^(-0.5) + 0.8 * f
#' }
#' dd <- data.frame(t = seq_len(n), x = x, y = y)
#' fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t, scaling = 0.5),
#'                distributions7::poisson_distrib(), dd)
#' ggplot2::autoplot(fit, term = "gas(p = 1, q = 1, time = t, scaling = 0.5)",
#'                   plot_type = plot_filter(response = TRUE))
#' @export
plot_filter <- function(level = 0.95, scale = "response", response = FALSE,
                        groups = NULL, vcov_type = bayesian_vcov(),
                        legend = TRUE) {
  check_level_grid(level, 2)
  scale <- match.arg(scale, c("response", "link"))
  if (!is.logical(response) || length(response) != 1L || is.na(response)) {
    stop("'response' must be TRUE or FALSE.", call. = FALSE)
  }
  if (response && scale != "response") {
    stop("'response = TRUE' draws the response beside the parameter, so it needs scale = \"response\".",
         call. = FALSE)
  }
  if (!is.null(groups) && (!length(groups) || anyNA(groups))) {
    stop("'groups' must be NULL or the labels of the groups to draw.",
         call. = FALSE)
  }
  if (!is.null(groups)) groups <- as.character(groups)
  vcov_type <- as_vcov_method(vcov_type)
  PlotFilter(label = "plot_filter()", level = level, scale = scale,
             response = response, groups = groups, vcov_type = vcov_type,
             legend = check_legend(legend))
}


#' The Data of a Filter Graph
#'
#' @description
#' Computes the path of the filtered parameter of each equation in `params`
#' with its band, by group, and the observed response.
#'
#' @param type A `PlotFilter` object.
#' @param bp The `build_info` of the term, or of its first component.
#' @param params The equations whose paths are drawn.
#' @param fit A [StatmodFit()].
#' @param ctx The context described at [plot_data()].
#'
#' @return A list with the data frames `path` and `obs`, and the attribute
#'   `"labels"`.
#'
#' @keywords internal
filter_data <- function(type, bp, params, fit, ctx) {
  spec <- fit@spec
  if (type@response) {
    if (length(params) > 1L) {
      stop(paste0("response = TRUE draws the response beside one parameter; ",
                  "the filter has several, so draw the response with a ",
                  "model of one gas() term."), call. = FALSE)
    }
    if (is.matrix(spec@response)) {
      stop("'response = TRUE' needs a response that is a vector.",
           call. = FALSE)
    }
    # the response shares the axis of the parameter only where the parameter
    # is the mean of the family, which is checked on the fitted values and
    # not assumed from the position of the parameter
    m <- tryCatch(as.numeric(predict(fit, what = "mean")),
                  error = function(e) NULL)
    th <- as.numeric(fit@fitted[[params]])
    if (is.null(m) || type@scale == "link" ||
        !isTRUE(all.equal(m, rep_len(th, length(m)), tolerance = 1e-8))) {
      stop(sprintf(paste0("response = TRUE draws the response on the axis of ",
                          "the filtered parameter, and '%s' is not the mean ",
                          "of the family."), params), call. = FALSE)
    }
  }
  n <- spec@n_obs
  lev <- bp$levels
  grp <- bp$group
  if (is.null(lev)) {
    lev <- "1"
    grp <- rep(1L, n)
  }
  glab <- lev[grp]
  keep <- if (is.null(type@groups)) rep(TRUE, n) else glab %in% type@groups
  if (!is.null(type@groups)) {
    bad <- setdiff(type@groups, lev)
    if (length(bad)) {
      stop(sprintf("'groups' names %s, which the term does not have. The groups are: %s.",
                   paste0("'", bad, "'", collapse = ", "),
                   paste(lev, collapse = ", ")), call. = FALSE)
    }
  } else if (length(lev) > 12L) {
    stop(sprintf(paste0("the term has %d groups and the graph draws at most ",
                        "twelve; give 'groups' with the labels to draw."),
                 length(lev)), call. = FALSE)
  }
  tm <- bp$times
  if (is.null(tm) || !(is.numeric(tm) || inherits(tm, c("Date", "POSIXct")))) {
    tm <- unlist(lapply(split(seq_len(n), grp), seq_along))[order(
      unlist(split(seq_len(n), grp)))]
  }
  paths <- lapply(params, function(p) {
    what <- if (type@scale == "link") paste0("link:", p) else p
    pr <- predict(fit, what = what, se = TRUE, level = type@level,
                  vcov_type = type@vcov_type)
    data.frame(x = tm, group = factor(glab, levels = lev), param = p,
               fit = pr$fit, lower = pr$lower, upper = pr$upper)[keep, ]
  })
  path <- do.call(rbind, paths)
  obs <- data.frame(x = numeric(0), group = factor(character(0)),
                    param = character(0), y = numeric(0))
  if (type@response) {
    obs <- data.frame(x = tm, group = factor(glab, levels = lev),
                      param = params, y = spec@response)[keep, ]
  }
  path$param <- factor(path$param, levels = params)
  out <- list(path = path, obs = obs)
  attr(out, "by") <- if (!is.null(bp$by)) paste(deparse(bp$by), collapse = "")
  attr(out, "labels") <- list(
    title = sprintf("%s: filtered %s", paste(params, collapse = ", "),
                    if (type@scale == "link") "predictor" else "parameter"),
    x = if (is.null(bp$time)) "observation" else
      paste(deparse(bp$time), collapse = ""),
    y = if (type@scale == "link") "predictor" else
      paste(params, collapse = ", "))
  out
}


S7::method(plot_data, list(PlotFilter, modelterms7::GasTerm)) <-
  function(type, x, fit, ctx, ...) {
    filter_data(type, x@build_info, ctx$param, fit, ctx)
  }

S7::method(plot_data, list(PlotFilter, modelterms7::GasJoint)) <-
  function(type, x, fit, ctx, ...) {
    comp <- x@terms
    drv <- vapply(comp, function(t) S7::S7_inherits(t, modelterms7::GasDriver),
                  logical(1))
    params <- x@params[!drv[x@params]]
    filter_data(type, comp[[params[1L]]]@build_info, params, fit, ctx)
  }


S7::method(plot_layers, PlotFilter) <- function(type, data, ...) {
  path <- data$path
  k_obs <- "observed response"
  k_fit <- est_label(if (type@scale == "link") "filtered predictor"
                     else "filtered parameter", type@level)
  out <- list()
  if (nrow(data$obs)) {
    out <- c(out, list(ggplot2::geom_point(
      data = data$obs, mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                              colour = k_obs),
      alpha = 0.5, size = 0.9)))
  }
  if (any(is.finite(path$lower) & is.finite(path$upper))) {
    out <- c(out, list(ggplot2::geom_ribbon(
      data = path, mapping = ggplot2::aes(x = .data$x, ymin = .data$lower,
                                          ymax = .data$upper),
      fill = "#1f5fa8", alpha = 0.18)))
  }
  out <- c(out, list(ggplot2::geom_line(
    data = path, mapping = ggplot2::aes(x = .data$x, y = .data$fit,
                                        colour = k_fit),
    linewidth = 0.8, key_glyph = draw_key_band),
    plot_key(colour = stats::setNames(c("#1f5fa8", "grey30"),
                                      c(k_fit, k_obs)))))
  ng <- nlevels(droplevels(path$group))
  np <- nlevels(droplevels(path$param))
  by <- strip_labeller(attr(data, "by"))
  if (ng > 1L && np > 1L) {
    out <- c(out, list(ggplot2::facet_grid(
      param ~ group, scales = "free",
      labeller = ggplot2::labeller(group = by))))
  } else if (ng > 1L) {
    out <- c(out, list(ggplot2::facet_wrap(ggplot2::vars(.data$group),
                                           scales = "free_x",
                                           labeller = by)))
  } else if (np > 1L) {
    out <- c(out, list(ggplot2::facet_wrap(ggplot2::vars(.data$param),
                                           ncol = 1L, scales = "free_y")))
  }
  out
}


# ---------------------------------------------------------------------------
# plot_news_impact()
# ---------------------------------------------------------------------------

#' @title S7 Class of the News-Impact Graph
#'
#' @description
#' The class that [plot_news_impact()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,range,at,n_points See [plot_news_impact()].
#'
#' @return A `PlotNewsImpact` object, inheriting from [term_plot].
#'
#' @name PlotNewsImpact
#' @keywords internal
PlotNewsImpact <- S7::new_class("PlotNewsImpact", parent = term_plot,
  properties = list(
    range = S7::class_numeric,
    at = S7::class_any,
    n_points = S7::class_numeric
  ))


#' Graph of the News Impact of a Score-Driven Term
#'
#' @description
#' Selects the graph of the one-step update that a score-driven term
#' ([modelterms7::gas()]) applies to its predictor as a function of the
#' standardized observation: the loading of the first score lag times the
#' (scaled) score at the observation.
#'
#' @details
#' The level of the filter is updated by \eqn{\kappa_1 u_t}, where
#' \eqn{u_t = s_t \mathcal{I}_t^{-d}} is the score \eqn{s_t} of the
#' log-likelihood in the predictor, scaled by the expected information
#' \eqn{\mathcal{I}_t} to the power \eqn{d} of the term, and \eqn{\kappa_1} is
#' the loading of the first lag. The update depends on the observation and on
#' the parameters of the density at which the score is read. To compare
#' families on one axis, the observation is standardized through its
#' distribution: for \eqn{z} on a grid, the observation is
#' \eqn{y(z) = F^{-1}(\Phi(z);\, \theta_\mathrm{ref})}, with \eqn{F} the
#' distribution function of the family, \eqn{\Phi} the standard normal one and
#' \eqn{\theta_\mathrm{ref}} the fitted parameters at a reference observation,
#' so that \eqn{z} is the quantile residual of the observation. For a gaussian
#' family \eqn{z} is the usual standardized residual.
#'
#' The reference observation is `at`, or by default the observation whose
#' filtered predictor is the median of the series. For a discrete family
#' \eqn{y(z)} takes integer values and the graph is a step function. With
#' loadings for lags beyond the first (\eqn{p > 1}) the graph shows the
#' immediate update only.
#'
#' A term whose parameters have subformulas has a loading that varies with the
#' observation, a family whose response is a matrix has no quantile function
#' on one axis, and both signal an error.
#'
#' @param range A numeric vector of length 2, the range of \eqn{z}.
#' @param at The row of the fitting data that supplies the parameters of the
#'   density, or `NULL` for the observation with the median filtered
#'   predictor.
#' @param n_points The number of points of the grid of \eqn{z}.
#'
#' @param legend Whether to draw the legend that names the elements of
#'   the graph.
#'
#' @return A `PlotNewsImpact` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_filter()]
#'
#' @examples
#' set.seed(4)
#' n <- 100
#' x <- rnorm(n)
#' f <- 0
#' y <- integer(n)
#' for (t in seq_len(n)) {
#'   lam <- exp(0.7 + 0.3 * x[t] + f)
#'   y[t] <- rpois(1, lam)
#'   f <- 0.15 * (y[t] - lam) * lam^(-0.5) + 0.8 * f
#' }
#' dd <- data.frame(t = seq_len(n), x = x, y = y)
#' fit <- statmod(y ~ x + gas(p = 1, q = 1, time = t, scaling = 0.5),
#'                distributions7::poisson_distrib(), dd)
#' ggplot2::autoplot(fit, term = "gas(p = 1, q = 1, time = t, scaling = 0.5)",
#'                   plot_type = plot_news_impact())
#' @export
plot_news_impact <- function(range = c(-6, 6), at = NULL, n_points = 401,
                             legend = TRUE) {
  if (!is.numeric(range) || length(range) != 2L || anyNA(range) ||
      !all(is.finite(range)) || range[1L] >= range[2L]) {
    stop("'range' must be two increasing finite numbers.", call. = FALSE)
  }
  if (!is.null(at) && (!is.numeric(at) || length(at) != 1L || is.na(at) ||
                       at < 1 || at != round(at))) {
    stop("'at' must be NULL or a row number.", call. = FALSE)
  }
  if (!is.numeric(n_points) || length(n_points) != 1L || !is.finite(n_points) ||
      n_points < 3 || n_points != round(n_points)) {
    stop("'n_points' must be a single whole number of at least 3.",
         call. = FALSE)
  }
  PlotNewsImpact(label = "plot_news_impact()", range = range, at = at,
                 n_points = n_points, legend = check_legend(legend))
}


S7::method(plot_data, list(PlotNewsImpact, modelterms7::GasTerm)) <-
  function(type, x, fit, ctx, ...) {
    spec <- fit@spec
    p <- ctx$param
    d <- spec@distrib
    n <- spec@n_obs
    if (is.matrix(spec@response)) {
      stop(paste0("the news impact needs a response that is a vector; this ",
                  "one is a matrix."), call. = FALSE)
    }
    if (length(x@submodels)) {
      stop(paste0("the loading of a term with subformulas on its parameters ",
                  "varies with the observation, and the news impact has one ",
                  "curve."), call. = FALSE)
    }
    kappa <- fit@structural[[ctx$name]]$parameter[["kappa1"]]
    scal <- filter_scaling(x)
    link <- d@link_params[[p]]
    f <- fit@fitted[[p]]
    at <- type@at
    if (is.null(at)) {
      at <- which.min(abs(f - stats::median(f)))
    } else if (at > n) {
      stop(sprintf("'at' is row %d and the fit has %d observations.", at, n),
           call. = FALSE)
    }
    th <- lapply(fit@fitted, function(v) v[at])
    z <- seq(type@range[1L], type@range[2L], length.out = type@n_points)
    if (type@range[1L] < 0 && type@range[2L] > 0) z <- sort(unique(c(z, 0)))
    y <- as.numeric(distributions7::distrib_quantile(d, stats::pnorm(z), th))
    ok <- is.finite(y)
    k <- distributions7::distrib_kernel(d, p)
    eta <- linkfunctions7::linkfun(link, th[[p]])
    s <- vapply(y[ok], function(v) as.numeric(k$score(v, th, eta)), numeric(1))
    info <- as.numeric(k$information(y[ok][1L], th, eta))
    impact <- rep(NA_real_, length(z))
    impact[ok] <- kappa * s * info^(-scal)
    out <- list(
      curve = data.frame(z = z, y = y, impact = impact)[ok, ],
      reference = data.frame(row = at, kappa = kappa, scaling = scal,
                             information = info,
                             parameter = th[[p]]))
    attr(out, "discrete") <- S7::S7_inherits(d, distributions7::discrete_distrib)
    attr(out, "labels") <- list(
      title = sprintf("%s: news impact", p),
      x = "standardized observation z",
      y = sprintf("update of the predictor of %s", p))
    out
  }


S7::method(plot_layers, PlotNewsImpact) <- function(type, data, ...) {
  m <- if (isTRUE(attr(data, "discrete"))) {
    ggplot2::geom_step(data = data$curve,
                       mapping = ggplot2::aes(x = .data$z, y = .data$impact),
                       colour = "#1f5fa8", linewidth = 0.9)
  } else {
    ggplot2::geom_line(data = data$curve,
                       mapping = ggplot2::aes(x = .data$z, y = .data$impact),
                       colour = "#1f5fa8", linewidth = 0.9)
  }
  list(ggplot2::geom_hline(yintercept = 0, colour = "grey50",
                           linetype = "dashed"),
       ggplot2::geom_vline(xintercept = 0, colour = "grey50",
                           linetype = "dashed"),
       m)
}


# ---------------------------------------------------------------------------
# plot_states() and plot_regime_density()
# ---------------------------------------------------------------------------

#' The Entry of the Evaluated Model for a Regime Term
#'
#' @description
#' Evaluates the predictors of the fit and returns the entry of the regime term
#' `name`, with the static predictor, the callbacks of the density and the
#' parameters of the term.
#'
#' @param fit A [StatmodFit()].
#' @param name The name of the term.
#'
#' @return A list as in the `regimes` element of [statmod_eta()].
#'
#' @keywords internal
regime_entry <- function(fit, name) {
  spec <- fit@spec
  ev <- statmod_eta(spec, statmod_design(spec), fit@coefficients)
  hit <- Filter(function(r) identical(r$term, name), ev$regimes)
  if (!length(hit)) {
    stop(sprintf("the fit has no regime term named '%s'.", name), call. = FALSE)
  }
  hit[[1L]]
}


#' The Parameters of a Regime Term
#'
#' @description
#' Returns the levels of the regimes, the transition matrix and the
#' stationary distribution at the fitted parameters of the term.
#'
#' @param r An entry of the `regimes` element of [statmod_eta()].
#'
#' @return A list with `mu` (the levels), `P` (the transition matrix, whose
#'   rows are the current regime) and `delta` (the stationary distribution).
#'
#' @keywords internal
regime_parameters <- function(r) {
  tm <- r$tm
  k <- tm@k
  psi <- unlist(r$psi)
  mu <- cumsum(c(psi[["level1"]], psi[paste0("gap", seq.int(2L, k))]))
  P <- as.matrix(parameters7::param_value(tm@transition,
                                          psi[tm@transition@free_names]))
  A <- diag(k) - P
  A[, k] <- 1
  delta <- as.numeric(solve(t(A), c(rep(0, k - 1L), 1)))
  list(mu = as.numeric(mu), P = P, delta = delta)
}


#' The Filtered State Probabilities of a Regime Term
#'
#' @description
#' Runs the forward recursion of the term at its fitted parameters:
#' \eqn{\alpha_t \propto f_t \circ (\alpha_{t-1} P)}, started from the
#' stationary distribution in each group, normalized at every step. The
#' normalizing constants are the contributions of the observations to the
#' log-likelihood of the term.
#'
#' @param r An entry of the `regimes` element of [statmod_eta()].
#'
#' @return A list with `prob`, a matrix with one row per observation and one
#'   column per regime, and `loglik`, the contributions of the observations.
#'
#' @keywords internal
regime_filtered <- function(r) {
  tm <- r$tm
  rp <- regime_parameters(r)
  k <- tm@k
  n <- tm@build_info$n
  idx <- seq_len(n)
  LF <- vapply(seq_len(k), function(j)
    as.numeric(r$cb$logdens(r$eta_static + rp$mu[j], idx)), numeric(n))
  LF <- matrix(LF, n, k)
  prob <- matrix(NA_real_, n, k)
  ll <- numeric(n)
  for (g in tm@build_info$order) {
    pred <- rp$delta
    for (t in g) {
      m <- max(LF[t, ])
      a <- pred * exp(LF[t, ] - m)
      ll[t] <- m + log(sum(a))
      prob[t, ] <- a / sum(a)
      pred <- as.numeric(prob[t, ] %*% rp$P)
    }
  }
  list(prob = prob, loglik = ll)
}


#' @title S7 Class of the State-Probability Graph
#'
#' @description
#' The class that [plot_states()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,smoothed See [plot_states()].
#'
#' @return A `PlotStates` object, inheriting from [term_plot].
#'
#' @name PlotStates
#' @keywords internal
PlotStates <- S7::new_class("PlotStates", parent = term_plot,
  properties = list(smoothed = S7::class_logical))


#' Graph of the State Probabilities of a Regime Term
#'
#' @description
#' Selects the graph of the probability of each regime of a
#' [modelterms7::regime()] term at each observation, drawn as stacked areas
#' over the order of the series.
#'
#' @details
#' With `smoothed = TRUE` (the default) the probabilities are the smoothed
#' ones, \eqn{P(S_t = j \mid y_1, \dots, y_n)}, as [statmod_latent()] reports
#' them. With `smoothed = FALSE` they are the filtered ones,
#' \eqn{P(S_t = j \mid y_1, \dots, y_t)}, from the forward recursion of the
#' term at its fitted parameters. The two coincide at the last observation of
#' each group. The horizontal axis is the variable named by `time` in the
#' term, or the row order when it has none, and a panel is drawn for each
#' group of a panel (at most twelve).
#'
#' @param smoothed Whether to draw the smoothed probabilities (`TRUE`) or the
#'   filtered ones (`FALSE`).
#'
#' @param legend Whether to draw the legend that names the elements of
#'   the graph.
#'
#' @return A `PlotStates` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_regime_density()],
#'   [statmod_latent()]
#'
#' @examples
#' set.seed(4)
#' n <- 80
#' s <- integer(n)
#' s[1] <- 1L
#' for (t in 2:n) s[t] <- if (runif(1) < 0.9) s[t - 1] else 3L - s[t - 1]
#' dd <- data.frame(t = seq_len(n), y = c(0, 3)[s] + rnorm(n, sd = 0.7))
#' fit <- statmod(y ~ regime(k = 2, time = t),
#'                distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, term = "regime(k = 2, time = t)")
#' ggplot2::autoplot(fit, term = "regime(k = 2, time = t)",
#'                   plot_type = plot_states(smoothed = FALSE))
#' @export
plot_states <- function(smoothed = TRUE, legend = TRUE) {
  if (!is.logical(smoothed) || length(smoothed) != 1L || is.na(smoothed)) {
    stop("'smoothed' must be TRUE or FALSE.", call. = FALSE)
  }
  PlotStates(label = "plot_states()", smoothed = smoothed,
             legend = check_legend(legend))
}


S7::method(plot_data, list(PlotStates, modelterms7::RegimeTerm)) <-
  function(type, x, fit, ctx, ...) {
    r <- regime_entry(fit, ctx$name)
    tm <- r$tm
    bp <- tm@build_info
    n <- bp$n
    k <- tm@k
    prob <- if (type@smoothed) {
      as.matrix(modelterms7::term_latent(tm, r$eta_static, fit@spec@response,
                                         r$cb$logdens, r$psi))
    } else regime_filtered(r)$prob
    ord <- bp$order
    if (length(ord) > 12L) {
      stop(sprintf(paste0("the term has %d groups and the graph draws at ",
                          "most twelve."), length(ord)), call. = FALSE)
    }
    grp <- integer(n)
    for (g in seq_along(ord)) grp[ord[[g]]] <- g
    tmv <- if (is.null(tm@time)) {
      pos <- integer(n)
      for (g in ord) pos[g] <- seq_along(g)
      pos
    } else {
      v <- tryCatch(eval(tm@time, plot_time_env(ctx, fit)), error = function(e) NULL)
      if (is.numeric(v) && length(v) == n) v else {
        pos <- integer(n)
        for (g in ord) pos[g] <- seq_along(g)
        pos
      }
    }
    lab <- sprintf("regime %d", seq_len(k))
    states <- data.frame(
      x = rep(tmv, k),
      group = factor(rep(grp, k), levels = seq_along(ord)),
      regime = factor(rep(lab, each = n), levels = lab),
      prob = as.numeric(prob))
    out <- list(states = states)
    attr(out, "labels") <- list(
      title = sprintf("%s: %s state probabilities", ctx$param,
                      if (type@smoothed) "smoothed" else "filtered"),
      x = if (is.null(tm@time)) "observation" else
        paste(deparse(tm@time), collapse = ""),
      y = "probability", fill = "regime")
    out
  }


#' The Environment in Which the Time of a Term Is Evaluated
#'
#' @description
#' Returns the data of the graph as an environment's data when it was given,
#' and an empty list otherwise, so that the `time` expression of a term is
#' read from the data where it can be.
#'
#' @param ctx The context described at [plot_data()].
#' @param fit A [StatmodFit()].
#'
#' @return A list or a data frame.
#'
#' @keywords internal
plot_time_env <- function(ctx, fit) {
  if (is.null(ctx$data)) list() else ctx$data
}


S7::method(plot_layers, PlotStates) <- function(type, data, ...) {
  s <- data$states
  out <- list(
    ggplot2::geom_area(data = s,
                       mapping = ggplot2::aes(x = .data$x, y = .data$prob,
                                              fill = .data$regime),
                       position = "stack", alpha = 0.8),
    ggplot2::scale_fill_brewer(palette = "Set1"),
    ggplot2::coord_cartesian(ylim = c(0, 1), expand = FALSE))
  if (nlevels(droplevels(s$group)) > 1L) {
    out <- c(out, list(ggplot2::facet_wrap(ggplot2::vars(.data$group),
                                           scales = "free_x")))
  }
  out
}


#' @title S7 Class of the Regime-Density Graph
#'
#' @description
#' The class that [plot_regime_density()] instantiates. Use the constructor,
#' which checks the properties.
#'
#' @param label,at,n_grid See [plot_regime_density()].
#'
#' @return A `PlotRegimeDensity` object, inheriting from [term_plot].
#'
#' @name PlotRegimeDensity
#' @keywords internal
PlotRegimeDensity <- S7::new_class("PlotRegimeDensity", parent = term_plot,
  properties = list(at = S7::class_any, n_grid = S7::class_numeric))


#' Graph of the Density of the Response under Each Regime
#'
#' @description
#' Selects the graph of the fitted density (the probability mass for a
#' discrete family) of the response under each regime of a
#' [modelterms7::regime()] term, at the covariates of a reference observation,
#' with one curve for each regime.
#'
#' @details
#' Under regime \eqn{j} the predictor of the equation of the term is the static
#' predictor plus the level \eqn{\mu_j} of the regime, and the other parameters
#' of the family are the fitted ones. The curves are the density at those
#' parameters for the reference observation `at`, which by default is the one
#' with the median static predictor. The horizontal range extends from the
#' 0.1 per cent quantile of the lowest regime to the 99.9 per cent quantile of
#' the highest.
#'
#' @param at The row of the fitting data that supplies the covariates, or
#'   `NULL` for the observation with the median static predictor.
#' @param n_grid The number of points of the grid of the response of a
#'   continuous family.
#'
#' @param legend Whether to draw the legend that names the elements of
#'   the graph.
#'
#' @return A `PlotRegimeDensity` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_states()]
#'
#' @examples
#' set.seed(4)
#' n <- 80
#' s <- integer(n)
#' s[1] <- 1L
#' for (t in 2:n) s[t] <- if (runif(1) < 0.9) s[t - 1] else 3L - s[t - 1]
#' dd <- data.frame(t = seq_len(n), y = c(0, 3)[s] + rnorm(n, sd = 0.7))
#' fit <- statmod(y ~ regime(k = 2, time = t),
#'                distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, term = "regime(k = 2, time = t)",
#'                   plot_type = plot_regime_density())
#' @export
plot_regime_density <- function(at = NULL, n_grid = 300, legend = TRUE) {
  if (!is.null(at) && (!is.numeric(at) || length(at) != 1L || is.na(at) ||
                       at < 1 || at != round(at))) {
    stop("'at' must be NULL or a row number.", call. = FALSE)
  }
  if (!is.numeric(n_grid) || length(n_grid) != 1L || !is.finite(n_grid) ||
      n_grid < 10 || n_grid != round(n_grid)) {
    stop("'n_grid' must be a single whole number of at least 10.",
         call. = FALSE)
  }
  PlotRegimeDensity(label = "plot_regime_density()", at = at, n_grid = n_grid,
                    legend = check_legend(legend))
}


S7::method(plot_data, list(PlotRegimeDensity, modelterms7::RegimeTerm)) <-
  function(type, x, fit, ctx, ...) {
    spec <- fit@spec
    d <- spec@distrib
    p <- ctx$param
    n <- spec@n_obs
    if (is.matrix(spec@response)) {
      stop(paste0("the regime densities need a response that is a vector; ",
                  "this one is a matrix."), call. = FALSE)
    }
    r <- regime_entry(fit, ctx$name)
    rp <- regime_parameters(r)
    at <- type@at
    if (is.null(at)) {
      at <- which.min(abs(r$eta_static - stats::median(r$eta_static)))
    } else if (at > n) {
      stop(sprintf("'at' is row %d and the fit has %d observations.", at, n),
           call. = FALSE)
    }
    link <- d@link_params[[p]]
    k <- length(rp$mu)
    ths <- lapply(seq_len(k), function(j) {
      th <- lapply(fit@fitted, function(v) v[at])
      th[[p]] <- linkfunctions7::linkinv(link, r$eta_static[at] + rp$mu[j])
      th
    })
    discrete <- S7::S7_inherits(d, distributions7::discrete_distrib)
    qs <- vapply(ths, function(th)
      as.numeric(distributions7::distrib_quantile(d, c(0.001, 0.999), th)),
      numeric(2))
    lo <- min(qs[1L, ], na.rm = TRUE)
    hi <- max(qs[2L, ], na.rm = TRUE)
    ys <- if (discrete) seq(floor(lo), ceiling(hi)) else
      seq(lo, hi, length.out = type@n_grid)
    lab <- sprintf("regime %d", seq_len(k))
    dens <- do.call(rbind, lapply(seq_len(k), function(j)
      data.frame(y = ys,
                 density = as.numeric(distributions7::distrib_pdf(d, ys,
                                                                  ths[[j]])),
                 regime = factor(lab[j], levels = lab))))
    out <- list(density = dens,
                reference = data.frame(row = at, regime = factor(lab, levels = lab),
                                       level = rp$mu))
    attr(out, "discrete") <- discrete
    attr(out, "labels") <- list(
      title = sprintf("%s: density under each regime", p),
      x = "response", y = if (discrete) "probability" else "density")
    attr(out, "labels")[[if (discrete) "fill" else "colour"]] <- "regime"
    out
  }


S7::method(plot_layers, PlotRegimeDensity) <- function(type, data, ...) {
  if (isTRUE(attr(data, "discrete"))) {
    list(ggplot2::geom_col(data = data$density,
                           mapping = ggplot2::aes(x = .data$y,
                                                  y = .data$density,
                                                  fill = .data$regime),
                           position = "dodge", alpha = 0.8),
         ggplot2::scale_fill_brewer(palette = "Set1"))
  } else {
    list(ggplot2::geom_line(data = data$density,
                            mapping = ggplot2::aes(x = .data$y,
                                                   y = .data$density,
                                                   colour = .data$regime),
                            linewidth = 0.9),
         ggplot2::scale_colour_brewer(palette = "Set1"))
  }
}


# ---------------------------------------------------------------------------
# Default graphs
# ---------------------------------------------------------------------------

S7::method(default_plot, modelterms7::SegTerm) <- function(term, ...) {
  plot_breakpoint()
}

S7::method(default_plot, modelterms7::GasTerm) <- function(term, ...) {
  plot_filter()
}

S7::method(default_plot, modelterms7::GasJoint) <- function(term, ...) {
  plot_filter()
}

S7::method(default_plot, modelterms7::RegimeTerm) <- function(term, ...) {
  plot_states()
}
