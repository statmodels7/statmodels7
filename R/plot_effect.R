#' @title S7 Class of the Effect Graph
#'
#' @description
#' The class that [plot_effect()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,residuals,binned,scale,at,rug,level,vcov_type,information,n_grid
#'   See [plot_effect()].
#'
#' @return A `PlotEffect` object, inheriting from [term_plot].
#'
#' @name PlotEffect
#' @keywords internal
PlotEffect <- S7::new_class("PlotEffect", parent = term_plot,
  properties = list(
    residuals = S7::class_logical,
    binned = S7::class_any,
    scale = S7::class_character,
    at = S7::class_any,
    rug = S7::class_logical,
    level = S7::class_numeric,
    vcov_type = S7::class_any,
    information = S7::class_any,
    n_grid = S7::class_numeric
  ))


#' Graph of the Effect of a Term
#'
#' @description
#' Selects the graph of the contribution of a term, or of a variable, to the
#' predictor of its equation, with a pointwise confidence band and the partial
#' residuals. For a numeric variable the contribution is drawn as a curve over
#' the range of the variable; for a factor, as one estimate with its interval
#' per level.
#'
#' @details
#' On the scale of the link the contribution of a term is \eqn{f(x) = X b},
#' evaluated on a grid of the variable, and the band is
#' \eqn{f(x) \pm z_{(1+\text{level})/2}\,\sqrt{x' V x}}, where \eqn{V} is the
#' block of [vcov.StatmodFit()] for the coefficients of the term. A smooth is
#' drawn as estimated, under the constraint that identifies it. The partial
#' residuals are \eqn{f(x_i) + (z_i - \eta_i)}, where
#' \eqn{z_i = \eta_i + s_i / w_i} is the working response of the equation,
#' \eqn{s_i} the score with respect to the predictor and \eqn{w_i} the
#' information selected by `information`. For a gaussian mean with the
#' identity link they are the ordinary partial residuals. For other equations
#' they are often skewed (for a scale parameter, the cloud lies below the
#' curve), and `binned` adds the means of the partial residuals within
#' classes of the variable.
#'
#' With `scale = "response"` the graph shows the parameter itself,
#' \eqn{g^{-1}(\eta)}, from [predict.StatmodFit()] with the other covariates
#' at `at` (or at their median and first level) and the random effects at
#' zero; on this scale the effect depends on the other covariates, and the
#' partial residuals are not drawn.
#'
#' Other variables that the term reads are held at their median or first
#' level, except for a factor with at most twelve levels, which is crossed
#' with the grid and drawn as one curve per level (`s(x, by = g)`).
#'
#' @param residuals Whether to draw the partial residuals.
#' @param binned `NULL`, or the number of classes of the variable within
#'   which the means of the partial residuals are drawn.
#' @param scale `"link"` (the default), the contribution to the predictor, or
#'   `"response"`, the parameter.
#' @param at A named list of values for the other covariates on the response
#'   scale, or `NULL`.
#' @param rug Whether to draw the observed values of the variable on the
#'   horizontal axis.
#' @param level The confidence level of the band.
#' @param vcov_type A [vcov_method] object, the variance of the band.
#' @param information An object from [distributions7::expected_info()] (the
#'   default) or [distributions7::observed_info()], the weight of the working
#'   response.
#' @param n_grid The number of points of the grid of a numeric variable.
#'
#' @param legend Whether to draw the legend that names the elements of
#'   the graph.
#'
#' @return A `PlotEffect` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_derivative()], [plot_surface()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(200), z = rnorm(200))
#' dd$y <- sin(2 * pi * dd$x) + 0.5 * dd$z + rnorm(200, sd = exp(-1 + dd$x))
#' fit <- statmod(y ~ z + s(x) | sigma ~ x, distributions7::gaussian1_distrib(),
#'                dd)
#' ggplot2::autoplot(fit, term = "s(x)", plot_type = plot_effect(), data = dd)
#' ggplot2::autoplot(fit, term = "x", param = "sigma",
#'                   plot_type = plot_effect(binned = 10), data = dd)
#' ggplot2::autoplot(fit, term = "x", param = "sigma",
#'                   plot_type = plot_effect(scale = "response"), data = dd)
#' @export
plot_effect <- function(residuals = TRUE, binned = NULL, scale = "link",
                        at = NULL, rug = TRUE, level = 0.95,
                        vcov_type = bayesian_vcov(),
                        information = distributions7::expected_info(),
                        n_grid = 200, legend = TRUE) {
  if (!is.logical(residuals) || length(residuals) != 1L || is.na(residuals)) {
    stop("'residuals' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(binned) && (!is.numeric(binned) || length(binned) != 1L ||
                           binned < 2 || binned != round(binned))) {
    stop("'binned' must be NULL or a whole number of at least 2.",
         call. = FALSE)
  }
  scale <- match.arg(scale, c("link", "response"))
  if (!is.null(at) && (!is.list(at) || is.null(names(at)))) {
    stop("'at' must be NULL or a named list.", call. = FALSE)
  }
  if (!is.logical(rug) || length(rug) != 1L || is.na(rug)) {
    stop("'rug' must be TRUE or FALSE.", call. = FALSE)
  }
  check_level_grid(level, n_grid)
  vcov_type <- as_vcov_method(vcov_type)
  if (!S7::S7_inherits(information, distributions7::ExpectedInfo) &&
      !S7::S7_inherits(information, distributions7::ObservedInfo)) {
    stop("'information' must be expected_info() or observed_info().",
         call. = FALSE)
  }
  PlotEffect(label = "plot_effect()", residuals = residuals, binned = binned,
             scale = scale, at = at, rug = rug, level = level,
             vcov_type = vcov_type, information = information,
             n_grid = n_grid, legend = check_legend(legend))
}


S7::method(plot_data, list(PlotEffect, modelterms7::additive_term)) <-
  function(type, x, fit, ctx, ...) {
    parts <- list(list(name = ctx$name, term = x, cols = ctx$cols))
    effect_data(type, parts, fit, ctx)
  }

S7::method(plot_data, list(PlotEffect, TermSum)) <-
  function(type, x, fit, ctx, ...) {
    bad <- Filter(function(p) !S7::S7_inherits(p$term,
                                               modelterms7::additive_term),
                  x@parts)
    if (length(bad)) {
      stop(sprintf(paste0("the effect of '%s' sums the terms in which it ",
                          "appears, and '%s' is not additive; draw that ",
                          "term by its name."), ctx$name, bad[[1L]]$name),
           call. = FALSE)
    }
    effect_data(type, x@parts, fit, ctx)
  }


#' The Data of an Effect Graph
#'
#' @description
#' Computes the curve (or the estimates of the levels), the band, the partial
#' residuals and their binned means for one or more parts of an equation.
#'
#' @param type A `PlotEffect` object.
#' @param parts A list of parts, each with `name`, `term` and `cols`.
#' @param fit A [StatmodFit()].
#' @param ctx The context described at [plot_data()].
#'
#' @return A list with the data frames `curve`, `points`, `bins` and `rug`,
#'   and the attributes `"labels"` and `"kind"` (`"numeric"` or `"factor"`).
#'
#' @keywords internal
effect_data <- function(type, parts, fit, ctx) {
  vars_all <- unique(unlist(lapply(parts, function(p)
    term_variables(p$term, p$name))))
  plotted <- ctx$vars
  data <- plot_needs_data(ctx, fit, type, vars = vars_all)
  # the variable on the axis: the first that is a numeric vector or a factor
  usable <- Filter(function(v) {
    vv <- data[[v]]
    !is.matrix(vv) && (is.numeric(vv) || is.factor(vv) || is.character(vv))
  }, plotted)
  if (!length(usable)) {
    stop(sprintf(paste0("the term '%s' reads no numeric or factor variable ",
                        "to put on an axis; plot_shrinkage() or ",
                        "plot_caterpillar() may describe it."), ctx$name),
         call. = FALSE)
  }
  num <- Filter(function(v) is.numeric(data[[v]]), usable)
  var <- if (length(num)) num[[1L]] else usable[[1L]]
  if (length(num) > 1L) {
    stop(sprintf(paste0("the term '%s' reads the numeric variables %s; ",
                        "plot_surface() draws a function of two."), ctx$name,
                 paste(num, collapse = " and ")), call. = FALSE)
  }
  # a factor is crossed with the grid only where it enters the drawn
  # contribution: the variables of a whole term, or those that share an
  # interaction with the variable in the columns of the linear part
  cross <- unique(unlist(lapply(parts, function(p)
    if (is.null(p$cols)) term_variables(p$term, p$name)
    else linpar_partners(p$term, var))))
  grid <- effect_grid(data, var, vars_all, type@n_grid, cross = cross)
  group <- attr(grid, "group")
  z <- stats::qnorm((1 + type@level) / 2)
  kind <- if (is.numeric(data[[var]])) "numeric" else "factor"

  if (type@scale == "response") {
    nd <- response_newdata(fit, data, grid, type@at)
    pr <- predict(fit, what = ctx$param, newdata = nd, se = TRUE,
                  level = type@level, random = "zero",
                  vcov_type = type@vcov_type)
    curve <- data.frame(x = grid[[var]], fit = pr$fit, lower = pr$lower,
                        upper = pr$upper)
    if (!is.null(group)) curve$group <- grid[[group]]
    points <- data.frame(x = numeric(0), y = numeric(0))
    bins <- points
    ylab <- sprintf("%s", ctx$param)
  } else {
    gfit <- parts_contribution(fit, ctx$param, parts, grid, type@vcov_type)
    curve <- data.frame(x = grid[[var]], fit = gfit$fit,
                        lower = gfit$fit - z * gfit$se,
                        upper = gfit$fit + z * gfit$se)
    if (!is.null(group)) curve$group <- grid[[group]]
    points <- data.frame(x = numeric(0), y = numeric(0))
    bins <- points
    if (type@residuals || !is.null(type@binned)) {
      ofit <- parts_contribution(fit, ctx$param, parts, NULL, type@vcov_type)
      wr <- working_residuals(fit, ctx$param, type@information)
      keep <- fit@spec@weights > 0 & is.finite(wr)
      pts <- data.frame(x = data[[var]], y = ofit$fit + wr)[keep, ]
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
                              x = var, y = ylab,
                              colour = group, fill = group)
  out
}


#' The Contribution of Several Parts of an Equation
#'
#' @description
#' Evaluates the sum of the contributions of `parts` on `newdata` (or on the
#' fitting rows), with the standard error that the joint variance of their
#' coefficients gives.
#'
#' @param fit A [StatmodFit()].
#' @param param The equation.
#' @param parts A list of parts, each with `name`, `term` and `cols`.
#' @param newdata A data frame, or `NULL`.
#' @param vcov_type A [vcov_method] object.
#'
#' @return A list with `fit` and `se`.
#'
#' @keywords internal
parts_contribution <- function(fit, param, parts, newdata, vcov_type) {
  des <- statmod_design(fit@spec)
  Xs <- list()
  idx <- integer(0)
  value <- 0
  for (p in parts) {
    # the block of a term that moves with its coefficients (seg(), nl()) is
    # the Jacobian of its contribution, so the contribution is term_value()
    # and the block enters only the delta method
    X <- if (is.null(newdata)) modelterms7::term_matrix(p$term)
      else modelterms7::term_predict(p$term, newdata = newdata)
    X <- as.matrix(X)
    j <- des[[param]]$blocks[[p$name]]
    if (!is.null(p$cols)) {
      X <- X[, p$cols, drop = FALSE]
      j <- j[p$cols]
      value <- value + drop(X %*% fit@coefficients[[param]][j])
    } else {
      value <- value + modelterms7::term_value(
        p$term, coef = fit@coefficients[[param]][j], newdata = newdata)
    }
    Xs[[length(Xs) + 1L]] <- X
    idx <- c(idx, j)
  }
  X <- do.call(cbind, Xs)
  nms <- paste0(param, ":", des[[param]]$coef_names[idx])
  V <- vcov(fit, vcov_type = vcov_type, readable = FALSE)[nms, nms,
                                                          drop = FALSE]
  list(fit = as.numeric(value), se = sqrt(pmax(rowSums((X %*% V) * X), 0)))
}


#' New Data for an Effect on the Response Scale
#'
#' @description
#' Completes a grid with every variable of the model, those of the terms and
#' those of the offsets: the variables in `at` take the given values, and the
#' others their median or first level in `data`.
#'
#' @param fit A [StatmodFit()].
#' @param data The data frame of the fit.
#' @param grid The grid of the effect.
#' @param at A named list, or `NULL`.
#'
#' @return A data frame with the rows of `grid`.
#'
#' @keywords internal
response_newdata <- function(fit, data, grid, at) {
  spec <- fit@spec
  vars <- unique(c(unlist(lapply(names(spec@terms), function(p)
    unlist(lapply(names(spec@terms[[p]]), function(k)
      term_variables(spec@terms[[p]][[k]], k))))),
    offset_variables(spec, data)))
  bad <- setdiff(names(at), vars)
  if (length(bad)) {
    stop(sprintf("'at' names %s, which the model does not read.",
                 paste0("'", bad, "'", collapse = ", ")), call. = FALSE)
  }
  nd <- grid
  for (v in setdiff(vars, names(grid))) {
    vo <- data[[v]]
    nd[[v]] <- if (v %in% names(at)) at[[v]]
      else if (is.matrix(vo)) matrix(apply(vo, 2, stats::median), nrow(nd),
                                     ncol(vo), byrow = TRUE)
      else if (is.numeric(vo)) stats::median(vo)
      else {
        f <- as.factor(vo)
        factor(levels(f)[1L], levels = levels(f))
      }
  }
  nd
}


#' Variables of the Offsets
#'
#' @description
#' Returns the columns of `data` that the offsets written in the formula of
#' the fit read.
#'
#' @details
#' The offsets are re-evaluated on new data when the fit predicts, so a graph
#' on the response scale needs a value for each of their variables.
#'
#' @param spec A [StatmodSpec()].
#' @param data The data frame of the fit.
#'
#' @return A character vector, possibly empty.
#'
#' @keywords internal
offset_variables <- function(spec, data) {
  eqs <- statmod_equations(spec@formula, spec@distrib@params)$equations
  v <- unlist(lapply(eqs, function(e)
    unlist(lapply(split_offsets(e)$offsets, all.vars))))
  intersect(unique(v), names(data))
}


S7::method(plot_layers, PlotEffect) <- function(type, data, ...) {
  kind <- attr(data, "kind")
  grouped <- "group" %in% names(data$curve)
  k_pts <- "partial residuals"
  k_fit <- est_label("estimated effect", type@level)
  k_bins <- sprintf("mean of the partial residuals in %d classes",
                    nrow(data$bins))
  out <- list()
  if (grouped) {
    # the colour names the levels of the crossed factor, so the elements are
    # named by aesthetics that no group uses: shape and linetype
    pa <- ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group,
                       shape = k_pts)
    key <- list(
      ggplot2::scale_shape_manual(
        name = NULL, values = stats::setNames(c(16, 23), c(k_pts, k_bins)),
        breaks = c(k_pts, k_bins)),
      ggplot2::scale_linetype_manual(
        name = NULL, values = stats::setNames("solid", k_fit)))
  } else {
    pa <- ggplot2::aes(x = .data$x, y = .data$y, colour = k_pts)
    key <- plot_key(
      colour = stats::setNames(c("#1f5fa8", "grey30"), c(k_fit, k_pts)),
      fill = if (nrow(data$bins))
        stats::setNames("#b2182b", k_bins))
  }
  if (identical(kind, "factor")) {
    if (nrow(data$points)) {
      out <- c(out, list(ggplot2::geom_jitter(data = data$points, mapping = pa,
                                              width = 0.15, height = 0,
                                              alpha = 0.35, size = 0.9)))
    }
    dodge <- ggplot2::position_dodge(width = if (grouped) 0.4 else 0)
    pr <- if (grouped) {
      ggplot2::geom_pointrange(
        data = data$curve,
        mapping = ggplot2::aes(x = .data$x, y = .data$fit, ymin = .data$lower,
                               ymax = .data$upper, colour = .data$group,
                               linetype = k_fit),
        position = dodge)
    } else {
      ggplot2::geom_pointrange(
        data = data$curve,
        mapping = ggplot2::aes(x = .data$x, y = .data$fit, ymin = .data$lower,
                               ymax = .data$upper, colour = k_fit))
    }
    # under treatment contrasts the first level is absorbed by the intercept:
    # its contribution is zero with no variance, and it is labeled as such
    ref <- data$curve[data$curve$upper - data$curve$lower == 0 &
                        data$curve$fit == 0, , drop = FALSE]
    lab <- if (nrow(ref)) {
      list(ggplot2::geom_text(
        data = ref, mapping = ggplot2::aes(x = .data$x, y = .data$fit),
        label = "reference level", vjust = 2, size = 3, colour = "grey20",
        position = dodge))
    }
    return(c(out, list(pr), lab, key))
  }
  if (nrow(data$points)) {
    out <- c(out, list(ggplot2::geom_point(data = data$points, mapping = pa,
                                           alpha = 0.35, size = 0.9)))
  }
  if (grouped) {
    out <- c(out, list(
      ggplot2::geom_ribbon(data = data$curve,
                           mapping = ggplot2::aes(x = .data$x,
                                                  ymin = .data$lower,
                                                  ymax = .data$upper,
                                                  fill = .data$group),
                           alpha = 0.18, show.legend = FALSE),
      ggplot2::geom_line(data = data$curve,
                         mapping = ggplot2::aes(x = .data$x, y = .data$fit,
                                                colour = .data$group,
                                                linetype = k_fit),
                         linewidth = 0.9, key_glyph = draw_key_band)))
  } else {
    out <- c(out, list(
      ggplot2::geom_ribbon(data = data$curve,
                           mapping = ggplot2::aes(x = .data$x,
                                                  ymin = .data$lower,
                                                  ymax = .data$upper),
                           fill = "#1f5fa8", alpha = 0.18),
      ggplot2::geom_line(data = data$curve,
                         mapping = ggplot2::aes(x = .data$x, y = .data$fit,
                                                colour = k_fit),
                         linewidth = 0.9, key_glyph = draw_key_band)))
  }
  if (nrow(data$bins)) {
    bins <- if (grouped) {
      ggplot2::geom_point(data = data$bins,
                          mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                                 shape = k_bins),
                          fill = "#b2182b", size = 2.6)
    } else {
      ggplot2::geom_point(data = data$bins,
                          mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                                 fill = k_bins),
                          shape = 23, size = 2.6)
    }
    out <- c(out, list(bins))
  }
  if (nrow(data$rug)) {
    out <- c(out, list(ggplot2::geom_rug(
      data = data$rug, mapping = ggplot2::aes(x = .data$x), sides = "b",
      alpha = 0.3)))
  }
  c(out, key)
}
