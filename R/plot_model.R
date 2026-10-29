#' @title S7 Classes of the Graphs of the Fitted Distribution
#'
#' @description
#' The classes that [plot_param()], [plot_quantiles()], [plot_density()] and
#' [plot_moments()] instantiate. Use the constructors, which check the
#' properties.
#'
#' @param label,param,over,at,level,vcov_type,n_grid,probs,points,moments,scale
#'   See the constructors.
#'
#' @return An object inheriting from [model_plot].
#'
#' @name PlotParam
#' @aliases PlotQuantiles PlotDensity PlotMoments
#' @keywords internal
PlotParam <- S7::new_class("PlotParam", parent = model_plot,
  properties = list(
    param = S7::class_character,
    over = S7::class_character,
    at = S7::class_any,
    level = S7::class_numeric,
    vcov_type = S7::class_any,
    n_grid = S7::class_numeric
  ))

#' @rdname PlotParam
#' @keywords internal
PlotQuantiles <- S7::new_class("PlotQuantiles", parent = model_plot,
  properties = list(
    over = S7::class_character,
    probs = S7::class_numeric,
    at = S7::class_any,
    points = S7::class_logical,
    n_grid = S7::class_numeric
  ))

#' @rdname PlotParam
#' @keywords internal
PlotDensity <- S7::new_class("PlotDensity", parent = model_plot,
  properties = list(
    at = S7::class_any,
    over = S7::class_any,
    n_grid = S7::class_numeric,
    scale = S7::class_numeric
  ))

#' @rdname PlotParam
#' @keywords internal
PlotMoments <- S7::new_class("PlotMoments", parent = model_plot,
  properties = list(
    over = S7::class_character,
    moments = S7::class_character,
    at = S7::class_any,
    n_grid = S7::class_numeric
  ))


#' Check the Covariate of a Model Graph
#'
#' @description
#' Signals an error unless `over` is a single string.
#'
#' @param over The name of a covariate.
#'
#' @return The string.
#'
#' @keywords internal
check_over <- function(over) {
  if (missing(over) || !is.character(over) || length(over) != 1L ||
      is.na(over)) {
    stop("'over' must be the name of a covariate, a single string.",
         call. = FALSE)
  }
  over
}

#' Check the Values of the Other Covariates
#'
#' @description
#' Signals an error unless `at` is `NULL` or a named list.
#'
#' @param at A named list or `NULL`.
#'
#' @return `at`.
#'
#' @keywords internal
check_at_list <- function(at) {
  if (!is.null(at) && (!is.list(at) || is.null(names(at)) ||
                       any(!nzchar(names(at))))) {
    stop("'at' must be NULL or a named list.", call. = FALSE)
  }
  at
}


#' Graph of a Fitted Parameter Against a Covariate
#'
#' @description
#' Selects the graph of a parameter of the fitted distribution as a function
#' of one covariate, on the scale of the parameter, with a pointwise
#' confidence band. For a numeric covariate the parameter is drawn as a curve
#' over the range of the covariate; for a factor, as one estimate with its
#' interval per level.
#'
#' @details
#' The values are those of [predict.StatmodFit()] with `what = param` and
#' `se = TRUE` on a grid of `over`. The other covariates are held at `at`, or
#' at their median (numeric) or first level (factor), and the random effects
#' are set to zero. The band is the interval of the parameter that
#' [predict.StatmodFit()] returns for the variance selected by `vcov_type`,
#' that is, the interval of the predictor mapped through the inverse link.
#'
#' @param param The name of a parameter of the family, a single string.
#' @param over The name of the covariate on the horizontal axis.
#' @param at A named list of values for the other covariates, or `NULL`.
#' @param level The confidence level of the band.
#' @param vcov_type A [vcov_method] object, the variance of the band.
#' @param n_grid The number of points of the grid of a numeric covariate.
#'
#' @return A `PlotParam` object, inheriting from [model_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_quantiles()], [plot_moments()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150))
#' dd$y <- 1 + dd$x + rnorm(150, sd = exp(-1 + dd$x))
#' fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, plot_type = plot_param(param = "sigma", over = "x"),
#'                   data = dd)
#' @export
plot_param <- function(param, over, at = NULL, level = 0.95,
                       vcov_type = bayesian_vcov(), n_grid = 200) {
  if (missing(param) || !is.character(param) || length(param) != 1L ||
      is.na(param)) {
    stop("'param' must be the name of a parameter, a single string.",
         call. = FALSE)
  }
  over <- check_over(over)
  at <- check_at_list(at)
  check_level_grid(level, n_grid)
  vcov_type <- as_vcov_method(vcov_type)
  PlotParam(label = "plot_param()", param = param, over = over, at = at,
            level = level, vcov_type = vcov_type, n_grid = n_grid)
}


#' Graph of the Conditional Quantiles Against a Covariate
#'
#' @description
#' Selects the graph of the quantile curves of the fitted distribution as
#' functions of one covariate, with the observed points behind them. The
#' median is drawn thicker and the outermost curves dashed.
#'
#' @details
#' The curve of probability \eqn{p} is
#' \eqn{q_p(x) = F^{-1}(p; \theta(x))}, where \eqn{\theta(x)} is the vector of
#' fitted parameters at the covariate value \eqn{x}, from
#' [predict.StatmodFit()] with the random effects set to zero, and
#' \eqn{F^{-1}} is [distributions7::distrib_quantile()]. The other covariates
#' are held at `at`, or at their median (numeric) or first level (factor).
#' For a discrete family the curves are step functions, since the quantile
#' takes integer values. The graph is not available for a multivariate
#' family.
#'
#' @param over The name of the covariate on the horizontal axis.
#' @param probs The probabilities of the curves, numbers in \eqn{(0, 1)}.
#' @param at A named list of values for the other covariates, or `NULL`.
#' @param points Whether to draw the observed values of the response against
#'   `over`.
#' @param n_grid The number of points of the grid of a numeric covariate.
#'
#' @return A `PlotQuantiles` object, inheriting from [model_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_param()], [plot_density()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150))
#' dd$y <- 1 + dd$x + rnorm(150, sd = exp(-1 + dd$x))
#' fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, plot_type = plot_quantiles(over = "x"), data = dd)
#' @export
plot_quantiles <- function(over, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                           at = NULL, points = TRUE, n_grid = 200) {
  over <- check_over(over)
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) ||
      any(probs <= 0 | probs >= 1)) {
    stop("'probs' must be numbers strictly between 0 and 1.", call. = FALSE)
  }
  at <- check_at_list(at)
  if (!is.logical(points) || length(points) != 1L || is.na(points)) {
    stop("'points' must be TRUE or FALSE.", call. = FALSE)
  }
  check_level_grid(0.5, n_grid)
  PlotQuantiles(label = "plot_quantiles()", over = over,
                probs = sort(unique(probs)), at = at, points = points,
                n_grid = n_grid)
}


#' Graph of Fitted Densities at Chosen Covariate Values
#'
#' @description
#' Selects the graph of the density, or the mass function, that the fit gives
#' at each row of a data frame of covariate values. With `over`, each density
#' is drawn vertically at its position on the horizontal axis, over the
#' observed data; without it, the densities are overlaid and coloured by row.
#'
#' @details
#' Each row of `at` is completed with the other covariates of the model, held
#' at their median (numeric) or first level (factor), and the random effects
#' are set to zero. The response range of each row runs from the quantile
#' 0.001 to the quantile 0.999 of its fitted distribution. For a continuous
#' family the density is evaluated on `n_grid` points of that range. For a
#' discrete family the mass is drawn at the integers from the smallest to the
#' largest of those quantiles over the rows, as bars or, beside the data, as
#' horizontal segments. The graph is not available for a multivariate family.
#'
#' @param at A data frame of covariate values, one fitted distribution per
#'   row.
#' @param over `NULL`, or the name of a numeric covariate of `at` that gives
#'   the position of each density on the horizontal axis.
#' @param n_grid The number of points of the response grid of a continuous
#'   family.
#' @param scale With `over`, the width of the widest density as a fraction of
#'   the smallest gap between positions (of 15 per cent of the range of the
#'   covariate for a single row).
#'
#' @return A `PlotDensity` object, inheriting from [model_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_quantiles()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150))
#' dd$y <- 1 + dd$x + rnorm(150, sd = exp(-1 + dd$x))
#' fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, plot_type = plot_density(at = data.frame(x = c(0.1, 0.5,
#'                                                                   0.9)),
#'                                            over = "x"), data = dd)
#' @export
plot_density <- function(at, over = NULL, n_grid = 200, scale = 0.9) {
  if (missing(at) || !is.data.frame(at) || !nrow(at)) {
    stop("'at' must be a data frame of covariate values with at least one row.",
         call. = FALSE)
  }
  if (!is.null(over)) {
    if (!is.character(over) || length(over) != 1L || is.na(over) ||
        !over %in% names(at)) {
      stop("'over' must be NULL or the name of a column of 'at'.",
           call. = FALSE)
    }
    if (!is.numeric(at[[over]])) {
      stop("'over' must name a numeric column of 'at'.", call. = FALSE)
    }
  }
  check_level_grid(0.5, n_grid)
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) ||
      scale <= 0) {
    stop("'scale' must be a single positive number.", call. = FALSE)
  }
  PlotDensity(label = "plot_density()", at = at, over = over,
              n_grid = n_grid, scale = scale)
}


#' Graph of the Moments of the Fitted Distribution Against a Covariate
#'
#' @description
#' Selects the graph of moments of the fitted distribution as functions of
#' one covariate, one panel for each moment with its own vertical scale.
#'
#' @details
#' The moments are those of [predict.StatmodFit()] with the random effects set
#' to zero, which evaluates the closed forms of the family or, for a family
#' without them, a quadrature or a series. The available moments are
#' `"mean"`, `"variance"`, `"std_dev"`, `"skewness"` and `"kurtosis"`. A
#' moment that does not exist for some covariate value (the mean of a Cauchy,
#' the variance of a Student t with two degrees of freedom or fewer) is
#' returned as `NaN` or `Inf` by the family and is left out of the graph.
#' The other covariates are held at `at`, or at their median (numeric) or
#' first level (factor). The graph is not available for a multivariate
#' family.
#'
#' @param over The name of the covariate on the horizontal axis.
#' @param moments A character vector of moments.
#' @param at A named list of values for the other covariates, or `NULL`.
#' @param n_grid The number of points of the grid of a numeric covariate.
#'
#' @return A `PlotMoments` object, inheriting from [model_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_param()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150))
#' dd$y <- rgamma(150, shape = 3, rate = 3 / exp(1 + dd$x))
#' fit <- statmod(y ~ x, distributions7::gamma1_distrib(), dd)
#' ggplot2::autoplot(fit, plot_type = plot_moments(over = "x",
#'                   moments = c("mean", "skewness")), data = dd)
#' @export
plot_moments <- function(over, moments = c("mean", "variance"), at = NULL,
                         n_grid = 200) {
  over <- check_over(over)
  ok <- c("mean", "variance", "std_dev", "skewness", "kurtosis")
  if (!is.character(moments) || !length(moments) || !all(moments %in% ok)) {
    stop(sprintf("'moments' must be among %s.",
                 paste0("\"", ok, "\"", collapse = ", ")), call. = FALSE)
  }
  at <- check_at_list(at)
  check_level_grid(0.5, n_grid)
  PlotMoments(label = "plot_moments()", over = unique(over),
              moments = unique(moments), at = at, n_grid = n_grid)
}


#' New Data on a Grid of One Covariate
#'
#' @description
#' Builds the grid of `over` (`n_grid` points over the range of a numeric
#' covariate, the levels of a factor) and completes it with the other
#' variables of the model, held at `at` or at their median or first level.
#'
#' @param type A graph type with the property `over`.
#' @param fit A [StatmodFit()].
#' @param ctx The context described at [plot_data()].
#'
#' @return A list with `nd` (the new data) and `x` (the values of `over`).
#'
#' @keywords internal
model_grid <- function(type, fit, ctx) {
  over <- type@over
  data <- plot_needs_data(ctx, fit, type, vars = over)
  v <- data[[over]]
  if (is.matrix(v) || !(is.numeric(v) || is.factor(v) || is.character(v))) {
    stop(sprintf("'over' must name a numeric or factor covariate; '%s' is not.",
                 over), call. = FALSE)
  }
  if (over %in% names(type@at)) {
    stop(sprintf("'at' must not name '%s', which is the variable on the axis.",
                 over), call. = FALSE)
  }
  grid <- effect_grid(data, over, over, type@n_grid)
  list(nd = response_newdata(fit, data, grid, type@at), x = grid[[over]])
}

#' Fitted Parameters on New Data
#'
#' @description
#' Calls [predict.StatmodFit()] for every parameter with the random effects
#' set to zero, and recycles each parameter to the number of rows.
#'
#' @param fit A [StatmodFit()].
#' @param nd A data frame.
#'
#' @return A named list of numeric vectors.
#'
#' @keywords internal
model_theta <- function(fit, nd) {
  th <- predict(fit, what = "parameter", newdata = nd, random = "zero")
  lapply(th, rep_len, length.out = nrow(nd))
}

#' Error for a Multivariate Family
#'
#' @description
#' Signals the error of a graph that is not available for a multivariate family.
#'
#' @param type,x The graph type and the family.
#'
#' @return Nothing: it always signals an error.
#'
#' @keywords internal
check_univariate <- function(type, x) {
  stop(sprintf(paste0("the graph %s is available for univariate families; ",
                      "'%s' is multivariate. predict(what = \"mean\") and ",
                      "predict(what = \"parameter\") give its fitted ",
                      "quantities."), type@label, class(x)[1L]),
       call. = FALSE)
}

# a factor covariate is drawn with points joined by lines per group
#' Kind of an Axis Variable
#'
#' @description
#' Returns `"numeric"` for a numeric vector and `"factor"` otherwise.
#'
#' @param x A vector.
#'
#' @return A string.
#'
#' @keywords internal
kind_of <- function(x) if (is.numeric(x)) "numeric" else "factor"


# plot_param -----------------------------------------------------------------

S7::method(plot_data, list(PlotParam, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) {
    if (!type@param %in% x@params) {
      stop(sprintf("'param' must be one of %s.",
                   paste0("\"", x@params, "\"", collapse = ", ")),
           call. = FALSE)
    }
    g <- model_grid(type, fit, ctx)
    pr <- predict(fit, what = type@param, newdata = g$nd, se = TRUE,
                  level = type@level, random = "zero",
                  vcov_type = type@vcov_type)
    curve <- data.frame(x = g$x, fit = pr$fit, lower = pr$lower,
                        upper = pr$upper)
    out <- list(curve = curve)
    attr(out, "kind") <- kind_of(g$x)
    attr(out, "labels") <- list(title = sprintf("%s against %s", type@param,
                                                type@over),
                                x = type@over, y = type@param)
    out
  }

S7::method(plot_layers, PlotParam) <- function(type, data, ...) {
  cv <- data$curve
  k_fit <- sprintf("fitted %s", type@param)
  if (identical(attr(data, "kind"), "factor")) {
    k_ci <- sprintf("%s with %s", k_fit,
                    band_label(type@level, "confidence interval"))
    return(list(ggplot2::geom_pointrange(
      data = cv, mapping = ggplot2::aes(x = .data$x, y = .data$fit,
                                        ymin = .data$lower, ymax = .data$upper,
                                        colour = k_ci)),
      plot_key(colour = stats::setNames("#1f5fa8", k_ci))))
  }
  k_band <- band_label(type@level)
  list(ggplot2::geom_ribbon(
         data = cv, mapping = ggplot2::aes(x = .data$x, ymin = .data$lower,
                                           ymax = .data$upper, fill = k_band),
         alpha = 0.18),
       ggplot2::geom_line(
         data = cv, mapping = ggplot2::aes(x = .data$x, y = .data$fit,
                                           colour = k_fit),
         linewidth = 0.9),
       plot_key(colour = stats::setNames("#1f5fa8", k_fit),
                fill = stats::setNames("#1f5fa8", k_band)))
}


# plot_quantiles -------------------------------------------------------------

#' The Data of a Quantile Graph
#'
#' @description
#' Computes the quantile curves and the observed points of [plot_quantiles()].
#'
#' @param type,x,fit,ctx,discrete The graph type, the family, the fit, the context of [plot_data()] and whether the family is discrete.
#'
#' @return A list with the data frames `curves` and `points`.
#'
#' @keywords internal
quantiles_data <- function(type, x, fit, ctx, discrete) {
  g <- model_grid(type, fit, ctx)
  th <- model_theta(fit, g$nd)
  curves <- do.call(rbind, lapply(type@probs, function(p) {
    data.frame(x = g$x, prob = p,
               q = as.numeric(distributions7::distrib_quantile(x, p = p,
                                                               theta = th)))
  }))
  y <- fit@spec@response
  pts <- if (type@points && is.null(dim(y)) && is.numeric(y)) {
    keep <- fit@spec@weights > 0
    data.frame(x = ctx$data[[type@over]], y = as.numeric(y))[keep, ,
                                                             drop = FALSE]
  } else data.frame(x = numeric(0), y = numeric(0))
  out <- list(curves = curves, points = pts)
  attr(out, "kind") <- kind_of(g$x)
  attr(out, "discrete") <- discrete
  attr(out, "labels") <- list(
    title = sprintf("quantiles of the response against %s", type@over),
    x = type@over, y = "response")
  out
}

S7::method(plot_data, list(PlotQuantiles, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) quantiles_data(type, x, fit, ctx, FALSE)

S7::method(plot_data, list(PlotQuantiles, distributions7::discrete_distrib)) <-
  function(type, x, fit, ctx, ...) quantiles_data(type, x, fit, ctx, TRUE)

S7::method(plot_data, list(PlotQuantiles, distributions7::multivariate_distrib)) <-
  function(type, x, fit, ctx, ...) check_univariate(type, x)

S7::method(plot_layers, PlotQuantiles) <- function(type, data, ...) {
  cv <- data$curves
  ps <- sort(unique(cv$prob))
  outer <- if (length(ps) > 1L) c(min(ps), max(ps)) else numeric(0)
  discrete <- isTRUE(attr(data, "discrete"))
  factor_x <- identical(attr(data, "kind"), "factor")
  k_pts <- "observations"
  out <- list()
  if (nrow(data$points)) {
    pm <- ggplot2::aes(x = .data$x, y = .data$y, colour = k_pts)
    out <- c(out, list(
      if (factor_x) ggplot2::geom_jitter(
        data = data$points, mapping = pm, width = 0.15, height = 0,
        alpha = 0.35, size = 0.9)
      else ggplot2::geom_point(
        data = data$points, mapping = pm, alpha = 0.35, size = 0.9)))
  }
  # the legend names each kind of line by the probabilities that it draws
  name_of <- function(d) {
    p <- sort(unique(d$prob))
    if (identical(p, 0.5)) return("median")
    sprintf("%s %s", if (length(p) == 1L) "quantile" else "quantiles",
            if (length(p) <= 2L) paste(format(p), collapse = " and ")
            else paste(format(p), collapse = ", "))
  }
  add <- function(d, lw, lty) {
    if (!nrow(d)) return(NULL)
    k <- name_of(d)
    m <- ggplot2::aes(x = .data$x, y = .data$q, group = .data$prob,
                      colour = k)
    if (factor_x) {
      list(ggplot2::geom_line(data = d, mapping = m, linewidth = lw,
                              linetype = lty),
           ggplot2::geom_point(data = d, mapping = m))
    } else if (discrete) {
      list(ggplot2::geom_step(data = d, mapping = m, linewidth = lw,
                              linetype = lty))
    } else {
      list(ggplot2::geom_line(data = d, mapping = m, linewidth = lw,
                              linetype = lty))
    }
  }
  is_out <- cv$prob %in% outer
  is_med <- cv$prob == 0.5
  inner <- !is_out & !is_med
  parts <- list(cv[is_med, ], cv[inner, ], cv[is_out & !is_med, ])
  keys <- unlist(lapply(Filter(nrow, parts), name_of))
  c(out, add(parts[[1L]], 1.4, "solid"), add(parts[[2L]], 0.6, "solid"),
    add(parts[[3L]], 0.6, "dashed"),
    plot_key(colour = stats::setNames(c(rep("#1f5fa8", length(keys)),
                                        "grey30"), c(keys, k_pts))))
}


# plot_density ---------------------------------------------------------------

#' The Rows of a Density Graph
#'
#' @description
#' Completes the covariate rows of [plot_density()] with the other variables of the model.
#'
#' @param type,fit,ctx The graph type, the fit and the context of [plot_data()].
#'
#' @return A list with `nd`, `at` and `data`.
#'
#' @keywords internal
density_rows <- function(type, fit, ctx) {
  data <- plot_needs_data(ctx, fit, type, vars = character(0))
  at <- type@at
  for (v in intersect(names(at), names(data))) {
    if (is.factor(data[[v]]) && !is.factor(at[[v]])) {
      at[[v]] <- factor(as.character(at[[v]]), levels = levels(data[[v]]))
    }
  }
  if (!is.null(type@over) && !type@over %in% names(data)) {
    stop(sprintf("'data' has no variable '%s'.", type@over), call. = FALSE)
  }
  nd <- response_newdata(fit, data, at, NULL)
  list(nd = nd, at = at, data = data)
}

#' Labels of the Rows of a Density Graph
#'
#' @description
#' Labels each row of `at` by the covariates that vary, or by its position.
#'
#' @param at A data frame.
#'
#' @return A factor.
#'
#' @keywords internal
row_labels <- function(at) {
  vary <- Filter(function(v) length(unique(at[[v]])) > 1L, names(at))
  lab <- if (length(vary)) {
    do.call(paste, c(lapply(vary, function(v) paste0(v, "=",
                                                    format(at[[v]], digits = 3,
                                                           trim = TRUE))),
                     sep = ", "))
  } else as.character(seq_len(nrow(at)))
  factor(lab, levels = unique(lab))
}

#' The Data of a Density Graph
#'
#' @description
#' Computes the density or mass of each row of [plot_density()].
#'
#' @param type,x,fit,ctx,discrete The graph type, the family, the fit, the context of [plot_data()] and whether the family is discrete.
#'
#' @return A list with the data frames `dens` and `points`.
#'
#' @keywords internal
density_data <- function(type, x, fit, ctx, discrete) {
  rw <- density_rows(type, fit, ctx)
  th <- model_theta(fit, rw$nd)
  n <- nrow(rw$nd)
  lab <- row_labels(rw$at)
  th_i <- function(i) lapply(th, `[`, i)
  lo <- vapply(seq_len(n), function(i)
    as.numeric(distributions7::distrib_quantile(x, p = 0.001, theta = th_i(i))),
    numeric(1))
  hi <- vapply(seq_len(n), function(i)
    as.numeric(distributions7::distrib_quantile(x, p = 0.999, theta = th_i(i))),
    numeric(1))
  if (discrete) {
    ys <- seq(min(lo), max(hi))
    dens <- do.call(rbind, lapply(seq_len(n), function(i) {
      data.frame(row = lab[i], index = i, y = ys,
                 dens = as.numeric(distributions7::distrib_pdf(
                   x, y = ys, theta = th_i(i))))
    }))
  } else {
    dens <- do.call(rbind, lapply(seq_len(n), function(i) {
      ys <- seq(lo[i], hi[i], length.out = type@n_grid)
      data.frame(row = lab[i], index = i, y = ys,
                 dens = as.numeric(distributions7::distrib_pdf(
                   x, y = ys, theta = th_i(i))))
    }))
  }
  pts <- data.frame(x = numeric(0), y = numeric(0))
  if (!is.null(type@over)) {
    pos <- rw$at[[type@over]]
    dens$x0 <- pos[dens$index]
    # the widest density spans `scale` times the smallest gap between
    # positions (or 15 per cent of the range of the data for a single row)
    gaps <- diff(sort(unique(pos)))
    span <- if (length(gaps)) min(gaps) else
      0.15 * diff(range(rw$data[[type@over]], na.rm = TRUE))
    dens$x <- dens$x0 + type@scale * span * dens$dens / max(dens$dens)
    y <- fit@spec@response
    if (is.null(dim(y)) && is.numeric(y)) {
      pts <- data.frame(x = rw$data[[type@over]], y = as.numeric(y))
    }
  }
  out <- list(dens = dens, points = pts)
  attr(out, "discrete") <- discrete
  attr(out, "over") <- type@over
  attr(out, "labels") <- if (is.null(type@over)) {
    # the discrete graph fills its bars and the continuous one colours its
    # lines, and ggplot2 reports a label for an aesthetic that is not mapped
    lab <- list(title = "fitted densities", x = "response", y = "density")
    if (discrete) c(lab, list(fill = "covariates"))
    else c(lab, list(colour = "covariates"))
  } else list(title = "fitted densities", x = type@over, y = "response")
  out
}

S7::method(plot_data, list(PlotDensity, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) density_data(type, x, fit, ctx, FALSE)

S7::method(plot_data, list(PlotDensity, distributions7::discrete_distrib)) <-
  function(type, x, fit, ctx, ...) density_data(type, x, fit, ctx, TRUE)

S7::method(plot_data, list(PlotDensity, distributions7::multivariate_distrib)) <-
  function(type, x, fit, ctx, ...) check_univariate(type, x)

S7::method(plot_layers, PlotDensity) <- function(type, data, ...) {
  d <- data$dens
  discrete <- isTRUE(attr(data, "discrete"))
  out <- list()
  if (is.null(attr(data, "over"))) {
    if (discrete) {
      return(list(ggplot2::geom_col(
        data = d, mapping = ggplot2::aes(x = .data$y, y = .data$dens,
                                         fill = .data$row),
        position = "dodge", alpha = 0.8, width = 0.8)))
    }
    return(list(ggplot2::geom_line(
      data = d, mapping = ggplot2::aes(x = .data$y, y = .data$dens,
                                       colour = .data$row),
      linewidth = 0.9)))
  }
  over <- attr(data, "over")
  k_pts <- "observations"
  k_dens <- sprintf("fitted %s at a row of 'at', drawn from its value of %s",
                    if (discrete) "probabilities" else "density", over)
  if (nrow(data$points)) {
    out <- c(out, list(ggplot2::geom_point(
      data = data$points, mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                                 colour = k_pts),
      alpha = 0.35, size = 0.9)))
  }
  key <- plot_key(colour = stats::setNames(c("#1f5fa8", "grey30"),
                                           c(k_dens, k_pts)))
  if (discrete) {
    c(out, list(ggplot2::geom_segment(
      data = d[d$dens > 0, ],
      mapping = ggplot2::aes(x = .data$x0, xend = .data$x, y = .data$y,
                             yend = .data$y, group = .data$index,
                             colour = k_dens),
      linewidth = 1.2)), key)
  } else {
    c(out, list(
      ggplot2::geom_ribbon(
        data = d, mapping = ggplot2::aes(y = .data$y, xmin = .data$x0,
                                         xmax = .data$x, group = .data$index),
        orientation = "y", fill = "#1f5fa8", alpha = 0.18),
      ggplot2::geom_path(
        data = d, mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                         group = .data$index,
                                         colour = k_dens),
        linewidth = 0.7)), key)
  }
}


# plot_moments ---------------------------------------------------------------

S7::method(plot_data, list(PlotMoments, distributions7::distrib)) <-
  function(type, x, fit, ctx, ...) {
    g <- model_grid(type, fit, ctx)
    curve <- do.call(rbind, lapply(type@moments, function(m) {
      v <- predict(fit, what = m, newdata = g$nd, random = "zero")
      data.frame(x = g$x, moment = m, value = as.numeric(v))
    }))
    curve$moment <- factor(curve$moment, levels = type@moments)
    curve <- curve[is.finite(curve$value), , drop = FALSE]
    out <- list(curve = curve)
    attr(out, "kind") <- kind_of(g$x)
    attr(out, "labels") <- list(
      title = sprintf("moments of the response against %s", type@over),
      x = type@over, y = NULL)
    out
  }

S7::method(plot_data, list(PlotMoments, distributions7::multivariate_distrib)) <-
  function(type, x, fit, ctx, ...) check_univariate(type, x)

S7::method(plot_layers, PlotMoments) <- function(type, data, ...) {
  cv <- data$curve
  m <- ggplot2::aes(x = .data$x, y = .data$value)
  geom <- if (identical(attr(data, "kind"), "factor")) {
    list(ggplot2::geom_point(data = cv, mapping = m, colour = "#1f5fa8",
                             size = 2))
  } else {
    list(ggplot2::geom_line(data = cv, mapping = m, colour = "#1f5fa8",
                            linewidth = 0.9))
  }
  c(geom, list(ggplot2::facet_wrap(ggplot2::vars(.data$moment),
                                   scales = "free_y")))
}
