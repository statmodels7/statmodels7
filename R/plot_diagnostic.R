#' Checks of the Arguments of a Diagnostic Graph
#'
#' @description
#' Signal an error unless the argument has the stated form.
#'
#' @param seed A single number or `NULL`.
#' @param x The value to check.
#' @param name The name of the argument in the message.
#' @param min The smallest admissible whole number.
#' @param null Whether `NULL` is admissible.
#'
#' @return `NULL`, invisibly.
#'
#' @name diag_checks
#' @keywords internal
diag_check_seed <- function(seed) {
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                         !is.finite(seed))) {
    stop("'seed' must be NULL or a single number.", call. = FALSE)
  }
  invisible(NULL)
}

#' @rdname diag_checks
#' @keywords internal
diag_check_count <- function(x, name, min = 1, null = FALSE) {
  if (null && is.null(x)) return(invisible(NULL))
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < min ||
      x != round(x)) {
    stop(sprintf("'%s' must be %sa whole number of at least %d.", name,
                 if (null) "NULL or " else "", as.integer(min)),
         call. = FALSE)
  }
  invisible(NULL)
}

#' @rdname diag_checks
#' @keywords internal
diag_check_name <- function(x, name, null = TRUE) {
  if (null && is.null(x)) return(invisible(NULL))
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("'%s' must be %sa single string.", name,
                 if (null) "NULL or " else ""), call. = FALSE)
  }
  invisible(NULL)
}

#' @rdname diag_checks
#' @keywords internal
diag_check_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)
  }
  invisible(NULL)
}


#' Evaluate an Expression under a Seed
#'
#' @description
#' Evaluates `expr` after `set.seed(seed)` when `seed` is not `NULL`, and
#' restores the random stream of the caller afterwards.
#'
#' @param seed A single number or `NULL`.
#' @param expr The expression.
#'
#' @return The value of `expr`.
#'
#' @keywords internal
diag_with_seed <- function(seed, expr) {
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1)
    }
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  expr
}


#' Quantile Residuals of a Response
#'
#' @description
#' Returns the quantile residuals of [residuals.StatmodFit()] for the response
#' of the fit, or for another response `y` evaluated under the fitted
#' parameters (no refit).
#'
#' @param fit A [StatmodFit()].
#' @param y A numeric vector of the length of the response, or `NULL`.
#' @param seed Passed to [residuals.StatmodFit()].
#'
#' @return A numeric vector.
#'
#' @keywords internal
diag_resid <- function(fit, y = NULL, seed = NULL) {
  if (!is.null(y)) fit@spec@response <- y
  residuals(fit, residual_type = "quantile", seed = seed)
}


#' Quantile Residuals of Simulated Responses
#'
#' @description
#' Draws `n_sim` responses with [simulate.StatmodFit()] and returns the
#' quantile residual of each under the fitted parameters.
#'
#' @param fit A [StatmodFit()].
#' @param n_sim The number of replicates.
#' @param seed A single number or `NULL`.
#'
#' @return A matrix with one row per observation and `n_sim` columns.
#'
#' @keywords internal
diag_sim_resid <- function(fit, n_sim, seed) {
  diag_with_seed(seed, {
    sims <- simulate(fit, nsim = n_sim)
    vapply(seq_len(n_sim), function(k) diag_resid(fit, y = sims[[k]]),
           numeric(fit@spec@n_obs))
  })
}


#' Classes of a Covariate with Equal Counts
#'
#' @description
#' Divides the observations by a covariate into `n_inter` classes of about
#' equal count (the quantile intervals of a numeric variable), or into the
#' levels of a factor.
#'
#' @param v The covariate.
#' @param n_inter The number of classes of a numeric variable.
#'
#' @return A factor with one element per observation.
#'
#' @keywords internal
diag_classes <- function(v, n_inter) {
  if (is.matrix(v)) {
    stop("a matrix covariate cannot divide the observations into classes.",
         call. = FALSE)
  }
  if (is.numeric(v)) {
    br <- unique(stats::quantile(v, probs = seq(0, 1, length.out = n_inter + 1),
                                 na.rm = TRUE, names = FALSE))
    if (length(br) < 3L) {
      stop("the covariate has too few distinct values to form classes.",
           call. = FALSE)
    }
    return(cut(v, breaks = br, include.lowest = TRUE, dig.lab = 3))
  }
  droplevels(as.factor(v))
}


#' Pointwise Band of the Normal Order Statistics
#'
#' @description
#' For a sample of size `n` from the standard normal distribution, returns the
#' expected plotting positions \eqn{\Phi^{-1}(p_i)}, with \eqn{p_i} from
#' [stats::ppoints()], and the pointwise interval of level `level` of the
#' order statistics. The \eqn{i}-th order statistic of \eqn{n} uniform
#' variables is \eqn{\text{Beta}(i, n + 1 - i)}, and the limits are the
#' quantiles of that law mapped by \eqn{\Phi^{-1}}.
#'
#' @param n The sample size.
#' @param level The level.
#'
#' @return A data frame with the columns `x`, `lower` and `upper`.
#'
#' @keywords internal
diag_order_band <- function(n, level) {
  i <- seq_len(n)
  a <- (1 - level) / 2
  data.frame(x = stats::qnorm(stats::ppoints(n)),
             lower = stats::qnorm(stats::qbeta(a, i, n + 1 - i)),
             upper = stats::qnorm(stats::qbeta(1 - a, i, n + 1 - i)))
}


#' The Data of a Quantile-Quantile or Worm Graph
#'
#' @description
#' Sorts the quantile residuals (within the classes of `over`, when given) and
#' pairs them with the normal plotting positions and the pointwise band; for a
#' worm plot the plotting position is subtracted from the residuals and from
#' the band.
#'
#' @param type A `PlotQq` or `PlotWorm` object.
#' @param fit A [StatmodFit()].
#' @param ctx The context described at [plot_data()].
#' @param worm Whether to detrend.
#' @param over The name of the covariate that divides the observations, or
#'   `NULL`.
#' @param n_inter The number of classes.
#'
#' @return A list with the data frames `points` and `band`.
#'
#' @keywords internal
diag_qq_data <- function(type, fit, ctx, worm, over = NULL, n_inter = 4) {
  r <- diag_resid(fit, seed = type@seed)
  cl <- NULL
  if (!is.null(over)) {
    plot_needs_data(ctx, fit, type, vars = over)
    cl <- diag_classes(ctx$data[[over]], n_inter)
  }
  groups <- if (is.null(cl)) list(seq_along(r)) else
    lapply(levels(cl), function(l) which(cl == l))
  R <- if (type@band == "simulated") {
    diag_sim_resid(fit, type@n_sim, type@seed)
  } else NULL
  a <- (1 - type@level) / 2
  pts <- list()
  bds <- list()
  for (k in seq_along(groups)) {
    idx <- groups[[k]]
    n <- length(idx)
    if (!n) next
    x <- stats::qnorm(stats::ppoints(n))
    y <- sort(r[idx], na.last = TRUE)
    band <- if (is.null(R)) diag_order_band(n, type@level) else {
      S <- apply(R[idx, , drop = FALSE], 2, sort, na.last = TRUE)
      S <- matrix(S, nrow = n)
      data.frame(x = x,
                 lower = apply(S, 1, stats::quantile, probs = a, na.rm = TRUE,
                               names = FALSE),
                 upper = apply(S, 1, stats::quantile, probs = 1 - a,
                               na.rm = TRUE, names = FALSE))
    }
    if (worm) {
      y <- y - x
      band$lower <- band$lower - band$x
      band$upper <- band$upper - band$x
    }
    p <- data.frame(x = x, y = y)
    p <- p[is.finite(p$y), , drop = FALSE]
    if (!is.null(cl)) {
      p$class <- rep(levels(cl)[k], nrow(p))
      band$class <- rep(levels(cl)[k], nrow(band))
    }
    pts[[k]] <- p
    bds[[k]] <- band
  }
  pts <- do.call(rbind, pts)
  bds <- do.call(rbind, bds)
  if (!is.null(cl)) {
    pts$class <- factor(pts$class, levels = levels(cl))
    bds$class <- factor(bds$class, levels = levels(cl))
  }
  list(points = pts, band = bds)
}


#' The Layers of a Quantile-Quantile or Worm Graph
#'
#' @description
#' Draws the band, the reference line, the points and, when the data carry a
#' `class` column, one facet for each class.
#'
#' @param data The list returned by [plot_data()].
#' @param worm Whether the graph is detrended.
#'
#' @return A list of ggplot2 components.
#'
#' @keywords internal
diag_qq_layers <- function(data, worm) {
  out <- list(
    ggplot2::geom_ribbon(data = data$band,
                         mapping = ggplot2::aes(x = .data$x,
                                                ymin = .data$lower,
                                                ymax = .data$upper),
                         fill = "#1f5fa8", alpha = 0.18),
    if (worm) ggplot2::geom_hline(yintercept = 0, colour = "grey50")
    else ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey50"),
    ggplot2::geom_point(data = data$points,
                        mapping = ggplot2::aes(x = .data$x, y = .data$y),
                        colour = "grey30", alpha = 0.35, size = 0.9))
  if ("class" %in% names(data$points)) {
    out <- c(out, list(ggplot2::facet_wrap(ggplot2::vars(.data$class))))
  }
  out
}


#' The Method of a Diagnostic Graph for a Multivariate Family
#'
#' @description
#' Signals an error: quantile residuals are defined for a univariate response.
#'
#' @param type A diagnostic [plot_type] object.
#' @param x,fit,ctx,... As in [plot_data()].
#'
#' @return Nothing: the function signals an error.
#'
#' @keywords internal
diag_multivariate <- function(type, x, fit, ctx, ...) {
  stop(sprintf(paste0("the graph %s reads quantile residuals, which are not ",
                      "defined for a multivariate family."), type@label),
       call. = FALSE)
}


# ---- plot_qq ---------------------------------------------------------------

#' @title S7 Class of the Quantile-Quantile Graph
#'
#' @description
#' The class that [plot_qq()] instantiates. Use the constructor, which checks
#' the properties.
#'
#' @param label,band,level,n_sim,seed See [plot_qq()].
#'
#' @return A `PlotQq` object, inheriting from [diagnostic_plot].
#'
#' @name PlotQq
#' @keywords internal
PlotQq <- S7::new_class("PlotQq", parent = diagnostic_plot,
  properties = list(band = S7::class_character, level = S7::class_numeric,
                    n_sim = S7::class_numeric, seed = S7::class_any))


#' Quantile-Quantile Graph of the Residuals
#'
#' @description
#' Selects the graph of the quantile residuals of a fit against the quantiles
#' of the standard normal distribution, with the line \eqn{y = x} and a
#' pointwise band.
#'
#' @details
#' The quantile residual \eqn{r_i = \Phi^{-1}(F(y_i; \hat\theta_i))} is
#' standard normal under a correct model ([residuals.StatmodFit()]), and the
#' sorted residuals are drawn against \eqn{\Phi^{-1}(p_i)}, with \eqn{p_i}
#' from [stats::ppoints()]. With `band = "order"` the band at position \eqn{i}
#' is the interval of level `level` of the \eqn{i}-th order statistic of \eqn{n}
#' uniform variables, a \eqn{\text{Beta}(i, n + 1 - i)} variable, mapped by
#' \eqn{\Phi^{-1}}. With `band = "simulated"` the band is the pointwise
#' interval of the sorted quantile residuals of `n_sim` responses drawn by
#' [simulate.StatmodFit()], each evaluated under the fitted parameters
#' without refitting. For a discrete family the residuals are randomized, and
#' `seed` makes them reproducible. The graph is not defined for a
#' multivariate family.
#'
#' @param band `"order"` (the default) or `"simulated"`.
#' @param level The level of the band.
#' @param n_sim The number of simulated responses of `band = "simulated"`.
#' @param seed A single number that seeds the randomization and the
#'   simulation, or `NULL` for the ambient state.
#'
#' @return A `PlotQq` object, inheriting from [diagnostic_plot].
#'
#' @seealso [plot_worm()], [diagnostics()], [autoplot.StatmodFit()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(120))
#' dd$y <- 1 + dd$x + rnorm(120, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, type = plot_qq())
#' ggplot2::autoplot(fit, type = plot_qq(band = "simulated", n_sim = 49,
#'                                       seed = 1))
#' @export
plot_qq <- function(band = "order", level = 0.95, n_sim = 99, seed = NULL) {
  band <- match.arg(band, c("order", "simulated"))
  diag_check_level(level)
  diag_check_count(n_sim, "n_sim", min = 2)
  diag_check_seed(seed)
  PlotQq(label = "plot_qq()", band = band, level = level, n_sim = n_sim,
         seed = seed)
}

S7::method(plot_data, list(PlotQq, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) {
    out <- diag_qq_data(type, fit, ctx, worm = FALSE)
    attr(out, "labels") <- list(title = "Normal quantile-quantile",
                                x = "theoretical quantiles",
                                y = "quantile residuals")
    out
  }

S7::method(plot_layers, PlotQq) <- function(type, data, ...) {
  diag_qq_layers(data, worm = FALSE)
}


# ---- plot_worm -------------------------------------------------------------

#' @title S7 Class of the Worm Plot
#'
#' @description
#' The class that [plot_worm()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,over,n_inter,band,level,n_sim,seed See [plot_worm()].
#'
#' @return A `PlotWorm` object, inheriting from [diagnostic_plot].
#'
#' @name PlotWorm
#' @keywords internal
PlotWorm <- S7::new_class("PlotWorm", parent = diagnostic_plot,
  properties = list(over = S7::class_any, n_inter = S7::class_numeric,
                    band = S7::class_character, level = S7::class_numeric,
                    n_sim = S7::class_numeric, seed = S7::class_any))


#' Worm Plot of the Residuals
#'
#' @description
#' Selects the detrended quantile-quantile graph of the quantile residuals: the
#' difference between each sorted residual and its theoretical quantile,
#' against the theoretical quantile, with the band of [plot_qq()] detrended
#' in the same way. With `over`, the observations are divided into classes of a
#' covariate and one worm is drawn for each class.
#'
#' @details
#' A worm plot (van Buuren and Fredriks, 2001) turns the departures from the
#' line of the quantile-quantile graph into departures from the horizontal
#' axis. An intercept different from zero shows a wrong mean of the residuals,
#' a slope a wrong variance, a U shape a wrong skewness and an S shape a wrong
#' kurtosis. With `over` a numeric covariate is divided into `n_inter` classes
#' of equal count (a factor into its levels), which requires `data` in
#' [autoplot.StatmodFit()]. The band is that of [plot_qq()] computed within
#' each class.
#'
#' @param over The name of a covariate, or `NULL`.
#' @param n_inter The number of classes of a numeric covariate.
#' @param band,level,n_sim,seed See [plot_qq()].
#'
#' @return A `PlotWorm` object, inheriting from [diagnostic_plot].
#'
#' @references
#' van Buuren, S. and Fredriks, M. (2001). Worm plot: a simple diagnostic
#' device for modelling growth reference curves. *Statistics in Medicine*
#' 20(8), 1259--1277.
#'
#' @seealso [plot_qq()], [plot_calibration()], [diagnostics()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(160))
#' dd$y <- 1 + dd$x + rnorm(160, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, type = plot_worm())
#' ggplot2::autoplot(fit, type = plot_worm(over = "x", n_inter = 3), data = dd)
#' @export
plot_worm <- function(over = NULL, n_inter = 4, band = "order", level = 0.95,
                      n_sim = 99, seed = NULL) {
  diag_check_name(over, "over")
  diag_check_count(n_inter, "n_inter", min = 2)
  band <- match.arg(band, c("order", "simulated"))
  diag_check_level(level)
  diag_check_count(n_sim, "n_sim", min = 2)
  diag_check_seed(seed)
  PlotWorm(label = "plot_worm()", over = over, n_inter = n_inter, band = band,
           level = level, n_sim = n_sim, seed = seed)
}

S7::method(plot_data, list(PlotWorm, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) {
    out <- diag_qq_data(type, fit, ctx, worm = TRUE, over = type@over,
                        n_inter = type@n_inter)
    attr(out, "labels") <- list(title = "Worm plot",
                                x = "theoretical quantiles",
                                y = "deviation of the quantile residuals")
    out
  }

S7::method(plot_layers, PlotWorm) <- function(type, data, ...) {
  diag_qq_layers(data, worm = TRUE)
}


# ---- plot_residuals --------------------------------------------------------

#' @title S7 Class of the Residual Graph
#'
#' @description
#' The class that [plot_residuals()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,against,param,smooth,seed See [plot_residuals()].
#'
#' @return A `PlotResiduals` object, inheriting from [diagnostic_plot].
#'
#' @name PlotResiduals
#' @keywords internal
PlotResiduals <- S7::new_class("PlotResiduals", parent = diagnostic_plot,
  properties = list(against = S7::class_character, param = S7::class_any,
                    smooth = S7::class_logical, seed = S7::class_any))


#' Graph of the Residuals against a Quantity
#'
#' @description
#' Selects the graph of the quantile residuals against the predicted
#' conditional median (in ranks), a fitted parameter, the index of the
#' observation or a covariate, with a smoother and the horizontal lines at 0
#' and \eqn{\pm 2}.
#'
#' @details
#' With `against = "median"` the abscissa is the rank, among the
#' observations, of the conditional median
#' \eqn{F^{-1}(1/2; \hat\theta_i)}; the ranks spread the points evenly
#' whatever the scale of the response, and the median is defined for every
#' family, including those without a location parameter. With
#' `against = "fitted"` the abscissa is the fitted value of the parameter
#' `param` on the parameter scale, and with `"index"` the row number. Any
#' other string names a covariate and requires `data`. A numeric abscissa
#' carries a loess smoother of the residuals (\code{\link[stats]{loess}}) when
#' `smooth = TRUE`, which should lie near zero under a correct model; a factor
#' is drawn with jittered points. For a discrete family the residuals are
#' randomized and `seed` makes them reproducible.
#'
#' @param against `"median"` (the default), `"fitted"`, `"index"`, or the name
#'   of a covariate.
#' @param param The parameter of `against = "fitted"`, otherwise `NULL`.
#' @param smooth Whether to draw the smoother.
#' @param seed A single number that seeds the randomization, or `NULL`.
#'
#' @return A `PlotResiduals` object, inheriting from [diagnostic_plot].
#'
#' @seealso [plot_qq()], [plot_acf()], [diagnostics()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150))
#' dd$y <- 1 + dd$x + rnorm(150, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, type = plot_residuals())
#' ggplot2::autoplot(fit, type = plot_residuals(against = "x"), data = dd)
#' @export
plot_residuals <- function(against = "median", param = NULL, smooth = TRUE,
                           seed = NULL) {
  diag_check_name(against, "against", null = FALSE)
  diag_check_name(param, "param")
  if (identical(against, "fitted") && is.null(param)) {
    stop("'against = \"fitted\"' needs 'param', the parameter whose fitted ",
         "value is the abscissa.", call. = FALSE)
  }
  if (!is.null(param) && !identical(against, "fitted")) {
    stop("'param' is read only with 'against = \"fitted\"'.", call. = FALSE)
  }
  if (!is.logical(smooth) || length(smooth) != 1L || is.na(smooth)) {
    stop("'smooth' must be TRUE or FALSE.", call. = FALSE)
  }
  diag_check_seed(seed)
  PlotResiduals(label = "plot_residuals()", against = against, param = param,
                smooth = smooth, seed = seed)
}


#' A Smoother of y on x
#'
#' @description
#' Fits [stats::loess()] and predicts it on 100 points of the range of `x`;
#' where the fit fails, the means of `y` within 20 classes of `x` of equal
#' count are returned.
#'
#' @param x,y Numeric vectors.
#'
#' @return A data frame with the columns `x` and `y`.
#'
#' @keywords internal
diag_smooth <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(unique(x)) < 5L) return(data.frame(x = numeric(0), y = numeric(0)))
  g <- seq(min(x), max(x), length.out = 100)
  out <- tryCatch({
    fit <- suppressWarnings(stats::loess(y ~ x,
                                         data = data.frame(x = x, y = y)))
    p <- suppressWarnings(stats::predict(fit, newdata = data.frame(x = g)))
    if (anyNA(p)) stop("loess")
    data.frame(x = g, y = as.numeric(p))
  }, error = function(e) {
    br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = 21),
                                 names = FALSE))
    cl <- cut(x, breaks = br, include.lowest = TRUE)
    data.frame(x = as.numeric(tapply(x, cl, mean)),
               y = as.numeric(tapply(y, cl, mean)))
  })
  out[stats::complete.cases(out), , drop = FALSE]
}

S7::method(plot_data, list(PlotResiduals, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) {
    r <- diag_resid(fit, seed = type@seed)
    n <- length(r)
    th <- fit@fitted
    against <- type@against
    if (against == "median") {
      med <- as.numeric(distributions7::distrib_quantile(x, 0.5, th))
      xv <- rank(med, ties.method = "average", na.last = "keep")
      xlab <- "rank of the predicted median"
    } else if (against == "fitted") {
      if (!type@param %in% names(th)) {
        stop(sprintf("'param' must be one of %s.",
                     paste0("\"", names(th), "\"", collapse = ", ")),
             call. = FALSE)
      }
      xv <- th[[type@param]]
      xlab <- sprintf("fitted %s", type@param)
    } else if (against == "index") {
      xv <- seq_len(n)
      xlab <- "index"
    } else {
      d <- plot_needs_data(ctx, fit, type, vars = against)
      xv <- d[[against]]
      if (is.matrix(xv)) {
        stop(sprintf("'%s' is a matrix; give a vector covariate.", against),
             call. = FALSE)
      }
      xlab <- against
    }
    pts <- data.frame(x = xv, y = r)
    pts <- pts[is.finite(r) & !is.na(xv), , drop = FALSE]
    sm <- if (type@smooth && is.numeric(pts$x)) diag_smooth(pts$x, pts$y)
      else data.frame(x = numeric(0), y = numeric(0))
    out <- list(points = pts, smooth = sm)
    attr(out, "labels") <- list(title = "Quantile residuals", x = xlab,
                                y = "quantile residuals")
    out
  }

S7::method(plot_layers, PlotResiduals) <- function(type, data, ...) {
  pa <- ggplot2::aes(x = .data$x, y = .data$y)
  pts <- if (is.numeric(data$points$x)) {
    ggplot2::geom_point(data = data$points, mapping = pa, colour = "grey30",
                        alpha = 0.35, size = 0.9)
  } else {
    ggplot2::geom_jitter(data = data$points, mapping = pa, colour = "grey30",
                         alpha = 0.35, size = 0.9, width = 0.15, height = 0)
  }
  out <- list(
    ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dashed",
                        colour = "grey50"),
    ggplot2::geom_hline(yintercept = 0, colour = "grey50"), pts)
  if (nrow(data$smooth)) {
    out <- c(out, list(ggplot2::geom_line(data = data$smooth, mapping = pa,
                                          colour = "#1f5fa8", linewidth = 0.9)))
  }
  out
}


# ---- plot_pit --------------------------------------------------------------

#' @title S7 Class of the PIT Histogram
#'
#' @description
#' The class that [plot_pit()] instantiates. Use the constructor, which checks
#' the properties.
#'
#' @param label,bins,seed,level See [plot_pit()].
#'
#' @return A `PlotPit` object, inheriting from [diagnostic_plot].
#'
#' @name PlotPit
#' @keywords internal
PlotPit <- S7::new_class("PlotPit", parent = diagnostic_plot,
  properties = list(bins = S7::class_numeric, seed = S7::class_any,
                    level = S7::class_numeric))


#' Histogram of the Probability Integral Transform
#'
#' @description
#' Selects the histogram of the values \eqn{F(y_i; \hat\theta_i)} of the
#' distribution function at the observations, on the density scale, with the
#' uniform density and a pointwise binomial band for each bin.
#'
#' @details
#' Under a correct model the values are uniform on \eqn{(0, 1)}. Each of the
#' `bins` classes of equal width receives a binomial number of the \eqn{n}
#' values, with probability \eqn{1/\text{bins}}, and the band of level `level`
#' is the interval of that binomial law divided by \eqn{n} times the bin
#' width. For a discrete family the value is randomized between
#' \eqn{F(y_i^-)} and \eqn{F(y_i)} (the values are \eqn{\Phi} of the quantile
#' residuals of [residuals.StatmodFit()]), and `seed` makes it reproducible.
#' A U shape indicates a response more dispersed than the model, a hump one
#' less dispersed, and a slope a bias in the location.
#'
#' @param bins The number of bins.
#' @param seed A single number that seeds the randomization, or `NULL`.
#' @param level The level of the band.
#'
#' @return A `PlotPit` object, inheriting from [diagnostic_plot].
#'
#' @seealso [plot_calibration()], [plot_qq()], [diagnostics()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150))
#' dd$y <- rpois(150, exp(1 + dd$x))
#' fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
#' ggplot2::autoplot(fit, type = plot_pit(bins = 10, seed = 1))
#' @export
plot_pit <- function(bins = 20, seed = NULL, level = 0.95) {
  diag_check_count(bins, "bins", min = 2)
  diag_check_seed(seed)
  diag_check_level(level)
  PlotPit(label = "plot_pit()", bins = bins, seed = seed, level = level)
}

S7::method(plot_data, list(PlotPit, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) {
    u <- stats::pnorm(diag_resid(fit, seed = type@seed))
    u <- u[is.finite(u)]
    n <- length(u)
    k <- type@bins
    br <- seq(0, 1, length.out = k + 1)
    cnt <- tabulate(findInterval(u, br, rightmost.closed = TRUE), nbins = k)
    w <- 1 / k
    a <- (1 - type@level) / 2
    bars <- data.frame(xmin = br[-(k + 1)], xmax = br[-1],
                       density = cnt / (n * w))
    band <- data.frame(xmin = bars$xmin, xmax = bars$xmax,
                       lower = stats::qbinom(a, n, w) / (n * w),
                       upper = stats::qbinom(1 - a, n, w) / (n * w))
    out <- list(bars = bars, band = band, values = data.frame(u = u))
    attr(out, "labels") <- list(title = "Probability integral transform",
                                x = "F(y)", y = "density")
    out
  }

S7::method(plot_layers, PlotPit) <- function(type, data, ...) {
  list(
    ggplot2::geom_rect(data = data$band,
                       mapping = ggplot2::aes(xmin = .data$xmin,
                                              xmax = .data$xmax,
                                              ymin = .data$lower,
                                              ymax = .data$upper),
                       fill = "#1f5fa8", alpha = 0.18),
    ggplot2::geom_rect(data = data$bars,
                       mapping = ggplot2::aes(xmin = .data$xmin,
                                              xmax = .data$xmax, ymin = 0,
                                              ymax = .data$density),
                       fill = "grey70", colour = "grey30", alpha = 0.6),
    ggplot2::geom_hline(yintercept = 1, colour = "#b2182b"))
}


# ---- plot_calibration ------------------------------------------------------

#' @title S7 Class of the Calibration Graph
#'
#' @description
#' The class that [plot_calibration()] instantiates. Use the constructor,
#' which checks the properties.
#'
#' @param label,probs,over,n_inter,level,seed See [plot_calibration()].
#'
#' @return A `PlotCalibration` object, inheriting from [diagnostic_plot].
#'
#' @name PlotCalibration
#' @keywords internal
PlotCalibration <- S7::new_class("PlotCalibration", parent = diagnostic_plot,
  properties = list(probs = S7::class_numeric, over = S7::class_any,
                    n_inter = S7::class_numeric, level = S7::class_numeric,
                    seed = S7::class_any))


#' Calibration of the Predicted Quantiles
#'
#' @description
#' Selects the graph of the observed coverage of the predicted quantiles
#' against their nominal level, with the diagonal and a pointwise band. With
#' `over`, one graph is drawn for each class of a covariate.
#'
#' @details
#' For a probability \eqn{p} the predicted quantile of observation \eqn{i} is
#' \eqn{F^{-1}(p; \hat\theta_i)}, and the observed coverage is the fraction of
#' observations at or below it, \eqn{\hat c(p) = n^{-1} \sum_i 1\{u_i \le p\}},
#' where \eqn{u_i = F(y_i; \hat\theta_i)} is the probability integral
#' transform, randomized for a discrete family (the values are \eqn{\Phi} of
#' the quantile residuals). Under a correct model \eqn{\hat c(p)} has mean
#' \eqn{p} and variance \eqn{p(1 - p)/n}, and the band is
#' \eqn{p \pm z_{(1 + \text{level})/2}\sqrt{p(1 - p)/n}}, limited to
#' \eqn{[0, 1]}. With `over`, the classes are formed as in [plot_worm()] and
#' require `data`.
#'
#' @param probs The nominal levels, numbers in \eqn{(0, 1)}.
#' @param over The name of a covariate, or `NULL`.
#' @param n_inter The number of classes of a numeric covariate.
#' @param level The level of the band.
#' @param seed A single number that seeds the randomization, or `NULL`.
#'
#' @return A `PlotCalibration` object, inheriting from [diagnostic_plot].
#'
#' @seealso [plot_pit()], [plot_worm()], [diagnostics()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150))
#' dd$y <- 1 + dd$x + rnorm(150, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, type = plot_calibration())
#' ggplot2::autoplot(fit, type = plot_calibration(over = "x", n_inter = 3),
#'                   data = dd)
#' @export
plot_calibration <- function(probs = seq(0.05, 0.95, by = 0.05), over = NULL,
                             n_inter = 4, level = 0.95, seed = NULL) {
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) ||
      any(probs <= 0 | probs >= 1)) {
    stop("'probs' must be numbers between 0 and 1.", call. = FALSE)
  }
  diag_check_name(over, "over")
  diag_check_count(n_inter, "n_inter", min = 2)
  diag_check_level(level)
  diag_check_seed(seed)
  PlotCalibration(label = "plot_calibration()", probs = sort(probs),
                  over = over, n_inter = n_inter, level = level, seed = seed)
}

S7::method(plot_data, list(PlotCalibration, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) {
    u <- stats::pnorm(diag_resid(fit, seed = type@seed))
    cl <- NULL
    if (!is.null(type@over)) {
      plot_needs_data(ctx, fit, type, vars = type@over)
      cl <- diag_classes(ctx$data[[type@over]], type@n_inter)
    }
    groups <- if (is.null(cl)) list(seq_along(u)) else
      lapply(levels(cl), function(l) which(cl == l))
    z <- stats::qnorm((1 + type@level) / 2)
    p <- type@probs
    pts <- lapply(seq_along(groups), function(k) {
      uk <- u[groups[[k]]]
      uk <- uk[is.finite(uk)]
      n <- length(uk)
      se <- z * sqrt(p * (1 - p) / n)
      d <- data.frame(p = p, observed = vapply(p, function(q) mean(uk <= q),
                                               numeric(1)),
                      lower = pmax(p - se, 0), upper = pmin(p + se, 1))
      if (!is.null(cl)) d$class <- rep(levels(cl)[k], length(p))
      d
    })
    pts <- do.call(rbind, pts)
    if (!is.null(cl)) pts$class <- factor(pts$class, levels = levels(cl))
    out <- list(points = pts)
    attr(out, "labels") <- list(title = "Calibration of the predicted quantiles",
                                x = "nominal level",
                                y = "observed coverage")
    out
  }

S7::method(plot_layers, PlotCalibration) <- function(type, data, ...) {
  out <- list(
    ggplot2::geom_ribbon(data = data$points,
                         mapping = ggplot2::aes(x = .data$p,
                                                ymin = .data$lower,
                                                ymax = .data$upper),
                         fill = "#1f5fa8", alpha = 0.18),
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey50"),
    ggplot2::geom_line(data = data$points,
                       mapping = ggplot2::aes(x = .data$p, y = .data$observed),
                       colour = "#1f5fa8"),
    ggplot2::geom_point(data = data$points,
                        mapping = ggplot2::aes(x = .data$p,
                                               y = .data$observed),
                        colour = "#1f5fa8", size = 1.4))
  if ("class" %in% names(data$points)) {
    out <- c(out, list(ggplot2::facet_wrap(ggplot2::vars(.data$class))))
  }
  out
}


# ---- plot_score ------------------------------------------------------------

#' @title S7 Class of the Score Graph
#'
#' @description
#' The class that [plot_score()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,param,against See [plot_score()].
#'
#' @return A `PlotScore` object, inheriting from [diagnostic_plot].
#'
#' @name PlotScore
#' @keywords internal
PlotScore <- S7::new_class("PlotScore", parent = diagnostic_plot,
  properties = list(param = S7::class_character, against = S7::class_character))


#' Graph of the Scores of an Equation
#'
#' @description
#' Selects the graph of the contribution of each observation to the score of
#' one equation, on the scale of the link, against the index of the
#' observation or a covariate. The five largest absolute contributions are
#' highlighted and labeled with the row number.
#'
#' @details
#' The contribution of observation \eqn{i} to the score of the equation
#' `param` is \eqn{s_i = \omega_i\, \partial \ell_i / \partial \eta_i}, with
#' \eqn{\ell_i} the log-likelihood of the observation, \eqn{\eta_i} the
#' predictor of the equation and \eqn{\omega_i} the prior weight. The
#' contributions sum, against the design, to the gradient of the fit, so a
#' large \eqn{|s_i|} marks an observation that strains that equation. They
#' are the numerator of the working residuals of [plot_effect()]. A covariate
#' as abscissa requires `data`.
#'
#' @param param The equation, a parameter of the family.
#' @param against `"index"` (the default) or the name of a covariate.
#'
#' @return A `PlotScore` object, inheriting from [diagnostic_plot].
#'
#' @seealso [plot_influence()], [plot_residuals()], [diagnostics()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(120))
#' dd$y <- 1 + dd$x + rnorm(120, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, type = plot_score(param = "mu"))
#' ggplot2::autoplot(fit, type = plot_score(param = "sigma", against = "x"),
#'                   data = dd)
#' @export
plot_score <- function(param, against = "index") {
  diag_check_name(param, "param", null = FALSE)
  diag_check_name(against, "against", null = FALSE)
  PlotScore(label = "plot_score()", param = param, against = against)
}


#' Check the Equation of a Diagnostic Graph
#'
#' @description
#' Signals an error unless `param` is a parameter of the fitted family.
#'
#' @param param A single string.
#' @param fit A [StatmodFit()].
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
diag_check_param <- function(param, fit) {
  ps <- names(fit@fitted)
  if (!param %in% ps) {
    stop(sprintf("'param' must be one of %s.",
                 paste0("\"", ps, "\"", collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
}

#' The Largest Absolute Values
#'
#' @description
#' Returns the data frame of the `k` largest absolute values of `y`, with the
#' abscissa and the index as a label.
#'
#' @param x,y Numeric vectors.
#' @param k How many.
#'
#' @return A data frame with the columns `x`, `y` and `label`.
#'
#' @keywords internal
diag_top <- function(x, y, k = 5L) {
  ok <- which(is.finite(y))
  o <- ok[order(abs(y[ok]), decreasing = TRUE)]
  o <- utils::head(o, k)
  data.frame(x = x[o], y = y[o], label = as.character(o))
}

S7::method(plot_data, list(PlotScore, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) {
    diag_check_param(type@param, fit)
    spec <- fit@spec
    g <- distributions7::distrib_gradient(x, spec@response, fit@fitted,
                                          scale = "link")[[type@param]]
    s <- as.numeric(g) * spec@weights
    n <- length(s)
    if (type@against == "index") {
      xv <- seq_len(n)
      xlab <- "index"
    } else {
      d <- plot_needs_data(ctx, fit, type, vars = type@against)
      xv <- d[[type@against]]
      if (!is.numeric(xv) || is.matrix(xv)) {
        stop(sprintf("'%s' must be a numeric vector to be an abscissa.",
                     type@against), call. = FALSE)
      }
      xlab <- type@against
    }
    out <- list(points = data.frame(x = xv, y = s), top = diag_top(xv, s))
    attr(out, "labels") <- list(
      title = sprintf("Score of %s", type@param), x = xlab,
      y = sprintf("score on the predictor of %s", type@param))
    out
  }

S7::method(plot_layers, PlotScore) <- function(type, data, ...) {
  list(
    ggplot2::geom_hline(yintercept = 0, colour = "grey50"),
    ggplot2::geom_point(data = data$points,
                        mapping = ggplot2::aes(x = .data$x, y = .data$y),
                        colour = "grey30", alpha = 0.35, size = 0.9),
    ggplot2::geom_point(data = data$top,
                        mapping = ggplot2::aes(x = .data$x, y = .data$y),
                        colour = "#b2182b", size = 1.8),
    ggplot2::geom_text(data = data$top,
                       mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                              label = .data$label),
                       colour = "#b2182b", size = 3, vjust = -0.7))
}


# ---- plot_influence --------------------------------------------------------

#' @title S7 Class of the Influence Graph
#'
#' @description
#' The class that [plot_influence()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,param,measure,adjust,information See [plot_influence()].
#'
#' @return A `PlotInfluence` object, inheriting from [diagnostic_plot].
#'
#' @name PlotInfluence
#' @keywords internal
PlotInfluence <- S7::new_class("PlotInfluence", parent = diagnostic_plot,
  properties = list(param = S7::class_character, measure = S7::class_character,
                    adjust = S7::class_logical, information = S7::class_any))


#' Graph of the Influence of the Observations on an Equation
#'
#' @description
#' Selects the graph of a one-step measure of the influence of each
#' observation on the coefficients of one equation: Cook's distance or the
#' leverage. The five largest values are highlighted and labeled with the row
#' number.
#'
#' @details
#' Let \eqn{X} be the design of the equation `param`, \eqn{s_i} the
#' contribution of observation \eqn{i} to the score of the equation on the
#' scale of the link, \eqn{w_i} its information on that scale, and \eqn{V} the
#' block of the Bayesian variance ([bayesian_vcov()]) of the coefficients of
#' the equation, the inverse of the penalized information \eqn{K}. Deleting
#' the observation moves the coefficients, to first order, by
#' \deqn{\Delta_i = V x_i s_i,}
#' with \eqn{x_i} the row of the design. The leverage is
#' \eqn{h_i = w_i x_i' V x_i}, the diagonal of the hat matrix, and its sum is
#' the effective number of coefficients \eqn{p_\text{eff}} of the equation. Cook's
#' distance is
#' \deqn{D_i = \Delta_i' K \Delta_i / p_\text{eff}.}
#' With `adjust = TRUE` (the default) \eqn{\Delta_i} is divided by
#' \eqn{1 - h_i}, which is the Newton step from the full fit on the
#' objective without the observation (the information \eqn{K - w_i x_i x_i'}
#' by the Sherman-Morrison formula). For a gaussian mean with constant
#' variance this step is the exact change of the least-squares coefficients,
#' and \eqn{D_i} is the Cook distance of [stats::cooks.distance()] multiplied
#' by \eqn{\hat s^2/\hat\sigma^2}, the ratio of the residual variance of
#' `lm` to the variance of the fit. With `adjust = FALSE` the step is
#' \eqn{V x_i s_i} as stated. The coefficients of random effects and
#' smooths are part of the equation, and the penalty is part of \eqn{K}.
#'
#' @param param The equation, a parameter of the family.
#' @param measure `"cook"` (the default) or `"leverage"`.
#' @param adjust Whether to divide the step by \eqn{1 - h_i}.
#' @param information An object from [distributions7::observed_info()] (the
#'   default) or [distributions7::expected_info()], the matrix of the
#'   variance and of the weights \eqn{w_i}.
#'
#' @return A `PlotInfluence` object, inheriting from [diagnostic_plot].
#'
#' @seealso [plot_score()], [diagnostics()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(100))
#' dd$y <- 1 + dd$x + rnorm(100, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, type = plot_influence(param = "mu"))
#' ggplot2::autoplot(fit, type = plot_influence(param = "mu",
#'                                              measure = "leverage"))
#' @export
plot_influence <- function(param, measure = "cook", adjust = TRUE,
                           information = distributions7::observed_info()) {
  diag_check_name(param, "param", null = FALSE)
  measure <- match.arg(measure, c("cook", "leverage"))
  if (!is.logical(adjust) || length(adjust) != 1L || is.na(adjust)) {
    stop("'adjust' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!S7::S7_inherits(information, distributions7::ExpectedInfo) &&
      !S7::S7_inherits(information, distributions7::ObservedInfo)) {
    stop("'information' must be expected_info() or observed_info().",
         call. = FALSE)
  }
  PlotInfluence(label = "plot_influence()", param = param, measure = measure,
                adjust = adjust, information = information)
}


#' One-Step Deletion Quantities of an Equation
#'
#' @description
#' Computes, for each observation, the first-order change of the coefficients
#' of an equation on deleting it, its leverage and Cook's distance, as
#' described at [plot_influence()].
#'
#' @param fit A [StatmodFit()].
#' @param param The equation.
#' @param information An object from [distributions7::observed_info()] or
#'   [distributions7::expected_info()].
#' @param adjust Whether to divide the step by one minus the leverage.
#'
#' @return A list with `delta` (a matrix with one row per observation and one
#'   column per coefficient of the equation that the variance reports),
#'   `leverage`, `cook` and `p_eff`.
#'
#' @keywords internal
influence_deltas <- function(fit, param, information, adjust = TRUE) {
  spec <- fit@spec
  des <- statmod_design(spec)[[param]]
  X <- as.matrix(des$X)
  nms <- paste0(param, ":", des$coef_names)
  V <- vcov(fit, vcov_type = bayesian_vcov(information = information),
            readable = FALSE)
  ok <- nms %in% rownames(V)
  ok[ok] <- apply(V[nms[ok], nms[ok], drop = FALSE], 1, function(v)
    all(is.finite(v)))
  Xk <- X[, ok, drop = FALSE]
  Vk <- V[nms[ok], nms[ok], drop = FALSE]
  d <- spec@distrib
  y <- spec@response
  th <- fit@fitted
  g <- distributions7::distrib_gradient(d, y, th, scale = "link")[[param]]
  key <- paste(param, param, sep = "_")
  h <- if (S7::S7_inherits(information, distributions7::ExpectedInfo)) {
    distributions7::distrib_expected_hessian(d, y, th, scale = "link",
                                             method = information@method)
  } else {
    distributions7::distrib_hessian(d, y, th, scale = "link")
  }
  wt <- spec@weights
  s <- as.numeric(g) * wt
  w <- -as.numeric(h[[key]]) * wt
  G <- Xk * s
  VG <- G %*% Vk
  lev <- unname(w * rowSums((Xk %*% Vk) * Xk))
  p_eff <- sum(lev[is.finite(lev)])
  cook <- unname(rowSums(VG * G)) / p_eff
  delta <- VG
  if (adjust) {
    delta <- delta / (1 - lev)
    cook <- cook / (1 - lev)^2
  }
  colnames(delta) <- des$coef_names[ok]
  list(delta = delta, leverage = lev, cook = cook, p_eff = p_eff)
}

S7::method(plot_data, list(PlotInfluence, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) {
    diag_check_param(type@param, fit)
    inf <- influence_deltas(fit, type@param, type@information, type@adjust)
    v <- if (type@measure == "cook") inf$cook else inf$leverage
    idx <- seq_along(v)
    out <- list(points = data.frame(x = idx, y = v), top = diag_top(idx, v))
    attr(out, "labels") <- list(
      title = sprintf("%s of %s", if (type@measure == "cook") "Cook distance"
                      else "Leverage", type@param),
      x = "index",
      y = if (type@measure == "cook") "Cook distance" else "leverage")
    out
  }

S7::method(plot_layers, PlotInfluence) <- function(type, data, ...) {
  list(
    ggplot2::geom_segment(data = data$points,
                          mapping = ggplot2::aes(x = .data$x, xend = .data$x,
                                                 y = 0, yend = .data$y),
                          colour = "grey60"),
    ggplot2::geom_point(data = data$top,
                        mapping = ggplot2::aes(x = .data$x, y = .data$y),
                        colour = "#b2182b", size = 1.8),
    ggplot2::geom_text(data = data$top,
                       mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                              label = .data$label),
                       colour = "#b2182b", size = 3, vjust = -0.7))
}


# ---- plot_acf --------------------------------------------------------------

#' @title S7 Class of the Autocorrelation Graph
#'
#' @description
#' The class that [plot_acf()] instantiates. Use the constructor, which checks
#' the properties.
#'
#' @param label,lag_max,squared,seed See [plot_acf()].
#'
#' @return A `PlotAcf` object, inheriting from [diagnostic_plot].
#'
#' @name PlotAcf
#' @keywords internal
PlotAcf <- S7::new_class("PlotAcf", parent = diagnostic_plot,
  properties = list(lag_max = S7::class_numeric, squared = S7::class_logical,
                    seed = S7::class_any))


#' Autocorrelation of the Residuals
#'
#' @description
#' Selects the graph of the sample autocorrelation of the quantile residuals,
#' or of their squares, in the order of the observations, with the limits
#' \eqn{\pm 1.96/\sqrt{n}}.
#'
#' @details
#' The autocorrelations are those of [stats::acf()] at the lags 1 to
#' `lag_max`. Under a correct model with independent observations the quantile
#' residuals are independent standard normal variables, the autocorrelations
#' have standard deviation near \eqn{1/\sqrt{n}}, and about five percent of
#' them fall outside the limits. The autocorrelation of the squares
#' (`squared = TRUE`) detects a variance that depends on the past. The order is
#' the row order of the data, so the graph is meaningful for a series, or for
#' data sorted within a panel; for a discrete family the residuals are
#' randomized and `seed` makes them reproducible.
#'
#' @param lag_max The largest lag.
#' @param squared Whether to use the squared residuals.
#' @param seed A single number that seeds the randomization, or `NULL`.
#'
#' @return A `PlotAcf` object, inheriting from [diagnostic_plot].
#'
#' @seealso [plot_residuals()], [diagnostics()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:150)
#' dd$y <- as.numeric(arima.sim(list(ar = 0.5), n = 150))
#' fit <- statmod(y ~ 1, distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, type = plot_acf(lag_max = 20))
#' @export
plot_acf <- function(lag_max = 30, squared = FALSE, seed = NULL) {
  diag_check_count(lag_max, "lag_max", min = 1)
  if (!is.logical(squared) || length(squared) != 1L || is.na(squared)) {
    stop("'squared' must be TRUE or FALSE.", call. = FALSE)
  }
  diag_check_seed(seed)
  PlotAcf(label = "plot_acf()", lag_max = lag_max, squared = squared,
          seed = seed)
}

S7::method(plot_data, list(PlotAcf, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) {
    r <- diag_resid(fit, seed = type@seed)
    if (type@squared) r <- r^2
    n <- sum(is.finite(r))
    r[!is.finite(r)] <- NA_real_
    lag_max <- min(type@lag_max, length(r) - 1L)
    a <- stats::acf(r, lag.max = lag_max, plot = FALSE,
                    na.action = stats::na.pass)
    out <- list(acf = data.frame(lag = as.numeric(a$lag)[-1L],
                                 acf = as.numeric(a$acf)[-1L]),
                limit = data.frame(value = 1.96 / sqrt(n)))
    attr(out, "labels") <- list(
      title = if (type@squared) "Autocorrelation of the squared residuals"
        else "Autocorrelation of the residuals",
      x = "lag", y = "autocorrelation")
    out
  }

S7::method(plot_layers, PlotAcf) <- function(type, data, ...) {
  list(
    ggplot2::geom_hline(yintercept = 0, colour = "grey50"),
    ggplot2::geom_hline(yintercept = c(-1, 1) * data$limit$value,
                        linetype = "dashed", colour = "#1f5fa8"),
    ggplot2::geom_segment(data = data$acf,
                          mapping = ggplot2::aes(x = .data$lag,
                                                 xend = .data$lag, y = 0,
                                                 yend = .data$acf),
                          colour = "grey30", linewidth = 0.8))
}


# ---- plot_ppc --------------------------------------------------------------

#' @title S7 Class of the Posterior Predictive Check
#'
#' @description
#' The class that [plot_ppc()] instantiates. Use the constructor, which checks
#' the properties.
#'
#' @param label,statistic,n_sim,seed See [plot_ppc()].
#'
#' @return A `PlotPpc` object, inheriting from [diagnostic_plot].
#'
#' @name PlotPpc
#' @keywords internal
PlotPpc <- S7::new_class("PlotPpc", parent = diagnostic_plot,
  properties = list(statistic = S7::class_any, n_sim = S7::class_numeric,
                    seed = S7::class_any))


#' Comparison of the Response with Simulated Responses
#'
#' @description
#' Selects the graph that compares the observed response with responses
#' simulated from the fitted distribution. Without `statistic`, a continuous
#' family gives the observed density with the densities of the simulated
#' responses, and a discrete family the observed counts of each value with the
#' mean and the 90 percent range of the simulated counts. With a function
#' `statistic`, the graph is the histogram of its value on the simulated
#' responses with the observed value as a line.
#'
#' @details
#' The responses are drawn by [simulate.StatmodFit()] at the fitted
#' parameters, `n_sim` of them, so they carry the variation of the response
#' and not the uncertainty of the estimates. The densities are computed with
#' [stats::density()] on a common grid over the range of the observed
#' response. The statistic is a function of a numeric vector that returns one
#' number, such as `max` or `function(y) mean(y == 0)`. The graph is not
#' defined for a multivariate family.
#'
#' @param statistic `NULL`, or a function of a numeric vector that returns a
#'   number.
#' @param n_sim The number of simulated responses.
#' @param seed A single number that seeds the simulation, or `NULL`.
#'
#' @return A `PlotPpc` object, inheriting from [diagnostic_plot].
#'
#' @seealso [simulate.StatmodFit()], [plot_qq()], [diagnostics()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(120))
#' dd$y <- rpois(120, exp(1 + dd$x))
#' fit <- statmod(y ~ x, distributions7::poisson_distrib(), dd)
#' ggplot2::autoplot(fit, type = plot_ppc(n_sim = 50, seed = 1))
#' ggplot2::autoplot(fit, type = plot_ppc(statistic = max, n_sim = 50,
#'                                        seed = 1))
#' @export
plot_ppc <- function(statistic = NULL, n_sim = 200, seed = NULL) {
  if (!is.null(statistic) && !is.function(statistic)) {
    stop("'statistic' must be NULL or a function.", call. = FALSE)
  }
  diag_check_count(n_sim, "n_sim", min = 2)
  diag_check_seed(seed)
  PlotPpc(label = "plot_ppc()", statistic = statistic, n_sim = n_sim,
          seed = seed)
}

S7::method(plot_data, list(PlotPpc, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) {
    y <- as.numeric(fit@spec@response)
    sims <- simulate(fit, nsim = type@n_sim, seed = type@seed)
    if (!is.null(type@statistic)) {
      f <- type@statistic
      one <- function(v) {
        out <- f(v)
        if (!is.numeric(out) || length(out) != 1L) {
          stop("'statistic' must return a single number.", call. = FALSE)
        }
        out
      }
      out <- list(sim = data.frame(value = vapply(sims, one, numeric(1))),
                  observed = data.frame(value = one(y)))
      attr(out, "kind") <- "statistic"
      attr(out, "labels") <- list(title = "Simulated statistic",
                                  x = "statistic", y = "count")
      return(out)
    }
    if (S7::S7_inherits(x, distributions7::discrete_distrib)) {
      vals <- sort(unique(y))
      if (diff(range(vals)) <= 200) vals <- seq(min(vals), max(vals))
      cnt <- function(v) vapply(vals, function(k) sum(v == k), numeric(1))
      M <- vapply(sims, cnt, numeric(length(vals)))
      M <- matrix(M, nrow = length(vals))
      out <- list(counts = data.frame(
        value = vals, observed = cnt(y), mean = rowMeans(M),
        lower = apply(M, 1, stats::quantile, probs = 0.05, names = FALSE),
        upper = apply(M, 1, stats::quantile, probs = 0.95, names = FALSE)))
      attr(out, "kind") <- "discrete"
      attr(out, "labels") <- list(title = "Observed and simulated counts",
                                  x = "value", y = "count")
      return(out)
    }
    rg <- range(y)
    pad <- 0.1 * diff(rg)
    from <- rg[1L] - pad
    to <- rg[2L] + pad
    dens <- function(v) {
      v <- v[is.finite(v)]
      d <- stats::density(v, from = from, to = to, n = 256)
      data.frame(x = d$x, y = d$y)
    }
    sim <- lapply(seq_along(sims), function(k) {
      d <- dens(sims[[k]])
      d$sim <- k
      d
    })
    out <- list(sim = do.call(rbind, sim), observed = dens(y))
    attr(out, "kind") <- "continuous"
    attr(out, "labels") <- list(title = "Observed and simulated densities",
                                x = "response", y = "density")
    out
  }

S7::method(plot_layers, PlotPpc) <- function(type, data, ...) {
  kind <- attr(data, "kind")
  if (kind == "statistic") {
    return(list(
      ggplot2::geom_histogram(data = data$sim,
                              mapping = ggplot2::aes(x = .data$value),
                              bins = 30, fill = "grey70", colour = "grey30"),
      ggplot2::geom_vline(data = data$observed,
                          mapping = ggplot2::aes(xintercept = .data$value),
                          colour = "#b2182b", linewidth = 0.9)))
  }
  if (kind == "discrete") {
    return(list(
      ggplot2::geom_linerange(data = data$counts,
                              mapping = ggplot2::aes(x = .data$value,
                                                     ymin = .data$lower,
                                                     ymax = .data$upper),
                              colour = "#1f5fa8", linewidth = 2, alpha = 0.4),
      ggplot2::geom_point(data = data$counts,
                          mapping = ggplot2::aes(x = .data$value,
                                                 y = .data$mean),
                          colour = "#1f5fa8", shape = 95, size = 5),
      ggplot2::geom_point(data = data$counts,
                          mapping = ggplot2::aes(x = .data$value,
                                                 y = .data$observed),
                          colour = "#b2182b", size = 1.8)))
  }
  list(
    ggplot2::geom_line(data = data$sim,
                       mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                              group = .data$sim),
                       colour = "grey70", alpha = 0.5, linewidth = 0.3),
    ggplot2::geom_line(data = data$observed,
                       mapping = ggplot2::aes(x = .data$x, y = .data$y),
                       colour = "#1f5fa8", linewidth = 0.9))
}


# ---- multivariate families -------------------------------------------------

for (cls in list(PlotQq, PlotWorm, PlotResiduals, PlotPit, PlotCalibration,
                 PlotScore, PlotInfluence, PlotAcf, PlotPpc)) {
  S7::method(plot_data, list(cls, distributions7::multivariate_distrib)) <-
    diag_multivariate
}
rm(cls)
