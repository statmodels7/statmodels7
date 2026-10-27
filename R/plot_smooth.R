#' S7 Class of the Derivative Graph
#'
#' @description
#' The class that [plot_derivative()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,order,level,simultaneous,n_sim,seed,vcov_type,n_grid
#'   See [plot_derivative()].
#'
#' @return A `PlotDerivative` object, inheriting from [term_plot].
#'
#' @name PlotDerivative
#' @keywords internal
PlotDerivative <- S7::new_class("PlotDerivative", parent = term_plot,
  properties = list(
    order = S7::class_numeric,
    level = S7::class_numeric,
    simultaneous = S7::class_logical,
    n_sim = S7::class_numeric,
    seed = S7::class_any,
    vcov_type = S7::class_any,
    n_grid = S7::class_numeric
  ))


#' Graph of the Derivative of a Smooth
#'
#' @description
#' Selects the graph of the first or second derivative of the effect of a
#' smooth of one numeric variable, with a confidence band. The zones of the
#' variable in which the band excludes zero, that is, where the effect
#' increases or decreases (first derivative) or is convex or concave (second
#' derivative), are shaded.
#'
#' @details
#' The effect of the smooth is \eqn{f(x) = X(x) b}, where \eqn{X(x)} is the
#' row of the block of the term at \eqn{x} and \eqn{b} its coefficients, so
#' that \eqn{f^{(m)}(x) = X^{(m)}(x) b} is linear in \eqn{b}. The rows
#' \eqn{X^{(m)}(x)} are obtained from the design of the term by one
#' finite-difference stencil of order \eqn{m} and accuracy four, with the
#' step and the weights of \pkg{numericals7} ([numericals7::fd_step()],
#' [numericals7::fd_weights()]); near the ends of the range of the variable
#' the nodes are shifted inside the range. The standard error of the
#' derivative is the exact delta-method value
#' \eqn{\sqrt{X^{(m)}(x)\, V\, X^{(m)}(x)'}}, where \eqn{V} is the block of
#' [vcov.StatmodFit()] for the coefficients of the term. The graph has the
#' derivative with respect to the variable on the axis, on the scale of the
#' link.
#'
#' With `simultaneous = FALSE` the band is the derivative
#' \eqn{\pm z_{(1+\text{level})/2}} standard errors. With
#' `simultaneous = TRUE` the factor \eqn{z} is replaced by the `level`
#' quantile of \eqn{\max_x |X^{(m)}(x)\, e| / \mathrm{se}(x)}, where \eqn{e}
#' is drawn `n_sim` times from \eqn{N(0, V)}, so that the band covers the whole
#' derivative over the grid with probability `level` (the simultaneous bands
#' of Ruppert, Wand and Carroll). A smooth with a `by` factor is drawn as one
#' curve for each level, each with its own critical value.
#'
#' @param order The order of the derivative, 1 or 2.
#' @param level The confidence level of the band.
#' @param simultaneous Whether the band is simultaneous over the grid.
#' @param n_sim The number of draws that give the simultaneous critical value.
#' @param seed `NULL`, or a single number that seeds the draws; the state of
#'   the random number generator is restored afterwards.
#' @param vcov_type A [vcov_method] object, the variance of the band.
#' @param n_grid The number of points of the grid of the variable.
#'
#' @return A `PlotDerivative` object, inheriting from [term_plot].
#'
#' @seealso [plot_effect()], [autoplot.StatmodFit()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150))
#' dd$y <- sin(2 * pi * dd$x) + rnorm(150, sd = 0.3)
#' fit <- statmod(y ~ s(x), distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, term = "s(x)", type = plot_derivative(), data = dd)
#' ggplot2::autoplot(fit, term = "s(x)",
#'                   type = plot_derivative(order = 2, simultaneous = TRUE,
#'                                          n_sim = 2000, seed = 1),
#'                   data = dd)
#' @export
plot_derivative <- function(order = 1, level = 0.95, simultaneous = FALSE,
                            n_sim = 10000, seed = NULL,
                            vcov_type = bayesian_vcov(), n_grid = 200) {
  if (!is.numeric(order) || length(order) != 1L || !order %in% c(1, 2)) {
    stop("'order' must be 1 or 2.", call. = FALSE)
  }
  check_level_grid(level, n_grid)
  if (!is.logical(simultaneous) || length(simultaneous) != 1L ||
      is.na(simultaneous)) {
    stop("'simultaneous' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1L || !is.finite(n_sim) ||
      n_sim < 100 || n_sim != round(n_sim)) {
    stop("'n_sim' must be a single whole number of at least 100.",
         call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                         !is.finite(seed))) {
    stop("'seed' must be NULL or a single number.", call. = FALSE)
  }
  vcov_type <- as_vcov_method(vcov_type)
  PlotDerivative(label = "plot_derivative()", order = order, level = level,
                 simultaneous = simultaneous, n_sim = n_sim, seed = seed,
                 vcov_type = vcov_type, n_grid = n_grid)
}


S7::method(plot_data, list(PlotDerivative, modelterms7::SmoothTerm)) <-
  function(type, x, fit, ctx, ...) {
    vars <- smooth_variables(x)
    if (length(vars) != 1L) {
      stop(sprintf(paste0("the term '%s' is a smooth of %d variables; ",
                          "plot_derivative() draws a smooth of one, and ",
                          "plot_surface() draws a function of two."),
                   ctx$name, length(vars)), call. = FALSE)
    }
    var <- vars[[1L]]
    vars_all <- unique(c(var, term_variables(x, ctx$name)))
    data <- plot_needs_data(ctx, fit, type, vars = vars_all)
    if (!is.numeric(data[[var]]) || is.matrix(data[[var]])) {
      stop(sprintf("the derivative needs a numeric variable, and '%s' is not.",
                   var), call. = FALSE)
    }
    grid <- effect_grid(data, var, vars_all, type@n_grid)
    group <- attr(grid, "group")
    lo <- min(data[[var]], na.rm = TRUE)
    hi <- max(data[[var]], na.rm = TRUE)
    D <- design_derivative(x, grid, var, type@order, lo, hi)
    cb <- term_coef_block(fit, ctx$param, ctx$name, ctx$cols, type@vcov_type)
    if (!is.null(ctx$cols)) D <- D[, ctx$cols, drop = FALSE]
    est <- drop(D %*% cb$coef)
    se <- sqrt(pmax(rowSums((D %*% cb$vcov) * D), 0))
    zp <- stats::qnorm((1 + type@level) / 2)
    groups <- if (is.null(group)) list(seq_len(nrow(grid))) else
      split(seq_len(nrow(grid)), grid[[group]])
    crit <- rep(zp, nrow(grid))
    if (type@simultaneous) {
      draws <- with_seed(type@seed, mvn_draws(cb$vcov, type@n_sim))
      for (g in groups) {
        crit[g] <- simultaneous_critical(D[g, , drop = FALSE], se[g], draws,
                                         type@level)
      }
    }
    curve <- data.frame(x = grid[[var]], fit = est, se = se,
                        lower = est - crit * se, upper = est + crit * se)
    if (!is.null(group)) curve$group <- grid[[group]]
    curve$significant <- curve$lower > 0 | curve$upper < 0
    zones <- do.call(rbind, lapply(groups, function(g)
      significant_zones(curve$x[g], curve$significant[g])))
    if (is.null(zones)) zones <- data.frame(xmin = numeric(0), xmax = numeric(0))
    critical <- data.frame(
      group = if (is.null(group)) "" else names(groups),
      pointwise = zp,
      used = vapply(groups, function(g) crit[g[1L]], numeric(1)),
      stringsAsFactors = FALSE)
    out <- list(curve = curve, zones = zones, critical = critical)
    attr(out, "labels") <- list(
      title = sprintf("%s: %s, %s derivative", ctx$param,
                      if (is.null(ctx$title)) ctx$name else ctx$title,
                      if (type@order == 1) "first" else "second"),
      x = var,
      y = sprintf("%s derivative of the effect",
                  if (type@order == 1) "first" else "second"),
      colour = group, fill = group)
    out
  }


#' The Variables of a Smooth
#'
#' @description
#' Returns the names of the variables that a smooth term smooths over, in
#' the order of the call, without its `by` variable.
#'
#' @param tm A [modelterms7::SmoothTerm].
#'
#' @return A character vector with one name for each variable; a variable
#'   written as an expression is named by its text.
#'
#' @keywords internal
smooth_variables <- function(tm) {
  vapply(tm@vars, function(v) paste(deparse(v), collapse = ""), character(1))
}


#' Rows of the Derivative of a Block
#'
#' @description
#' Computes the derivative of order `order` with respect to `var` of every
#' row of the block of a term at the rows of `grid`, by one central
#' finite-difference stencil of accuracy four on the rows that
#' [modelterms7::term_predict()] returns.
#'
#' @details
#' The step is \eqn{h = 0.05\,\mathrm{fd\_step}(1, m, 4)\,(b - a)}, where \eqn{m} is
#' the order and \eqn{[a, b]} the range of the variable. The center of the
#' stencil is moved inside \eqn{[a + rh, b - rh]}, where \eqn{r} is the reach
#' of the stencil, and the weights of the shifted nodes are computed for the
#' point itself, so that the nodes stay inside the range over which the basis
#' is defined.
#'
#' @param tm A built [modelterms7::SmoothTerm].
#' @param grid A data frame with every variable the term reads.
#' @param var The variable that is differentiated.
#' @param order The order of the derivative.
#' @param lo,hi The range of the variable.
#'
#' @return A matrix with one row for each row of `grid` and one column for
#'   each coefficient of the term.
#'
#' @keywords internal
design_derivative <- function(tm, grid, var, order, lo, hi) {
  off <- numericals7::fd_offsets(order, accuracy = 4L)
  r <- off$reach
  ks <- off$central
  h <- 0.05 * numericals7::fd_step(1, order, accuracy = 4L) * (hi - lo)
  x <- grid[[var]]
  centre <- pmin(pmax(x, lo + r * h), hi - r * h)
  shift <- round((centre - x) / h, 10)
  n <- nrow(grid)
  nodes <- grid[rep(seq_len(n), length(ks)), , drop = FALSE]
  nodes[[var]] <- rep(centre, length(ks)) + rep(ks, each = n) * h
  X <- as.matrix(modelterms7::term_predict(tm, newdata = nodes))
  p <- ncol(X)
  D <- matrix(0, n, p)
  for (s in unique(shift)) {
    rows <- which(shift == s)
    w <- numericals7::fd_weights(ks + s, order)
    for (j in seq_along(ks)) {
      if (w[j] != 0) {
        D[rows, ] <- D[rows, ] + w[j] * X[(j - 1L) * n + rows, , drop = FALSE]
      }
    }
  }
  D / h^order
}


#' Run an Expression Under a Seed
#'
#' @description
#' Evaluates `expr` after `set.seed(seed)` when `seed` is not `NULL`, and
#' restores the state of the random number generator afterwards.
#'
#' @param seed `NULL` or a single number.
#' @param expr The expression.
#'
#' @return The value of `expr`.
#'
#' @keywords internal
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  set.seed(seed)
  expr
}


#' Draws from a Multivariate Normal Distribution
#'
#' @description
#' Draws `n` vectors from \eqn{N(0, V)} through the eigendecomposition of
#' \eqn{V}, which also covers a singular \eqn{V}.
#'
#' @param V A symmetric positive semi-definite matrix.
#' @param n The number of draws.
#'
#' @return A matrix with `n` rows and `ncol(V)` columns.
#'
#' @keywords internal
mvn_draws <- function(V, n) {
  ev <- eigen((V + t(V)) / 2, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  Z <- matrix(stats::rnorm(n * ncol(V)), n, ncol(V))
  (Z %*% (sqrt(lam) * t(ev$vectors)))
}


#' The Critical Value of a Simultaneous Band
#'
#' @description
#' Returns the `level` quantile of \eqn{\max_x |D(x) e| / \mathrm{se}(x)} over
#' the draws \eqn{e}, with the maximum over the rows of `D` whose standard
#' error is positive.
#'
#' @param D The matrix of the rows of the derivative of the block.
#' @param se The standard errors of the rows.
#' @param draws The draws of the coefficients, one in each row.
#' @param level The level.
#'
#' @return A single number.
#'
#' @keywords internal
simultaneous_critical <- function(D, se, draws, level) {
  ok <- se > 1e-12 * max(se, 1e-300)
  Tm <- abs(D[ok, , drop = FALSE] %*% t(draws)) / se[ok]
  stats::quantile(apply(Tm, 2L, max), level, names = FALSE, type = 7)
}


#' The Zones in Which a Band Excludes Zero
#'
#' @description
#' Returns the maximal runs of consecutive grid points at which `flag` is
#' `TRUE`, as intervals of the variable.
#'
#' @param x The grid, increasing.
#' @param flag A logical vector as long as `x`.
#'
#' @return A data frame with `xmin` and `xmax`, one row for each run.
#'
#' @keywords internal
significant_zones <- function(x, flag) {
  r <- rle(as.vector(flag))
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  keep <- r$values
  data.frame(xmin = x[starts[keep]], xmax = x[ends[keep]])
}


S7::method(plot_layers, PlotDerivative) <- function(type, data, ...) {
  grouped <- "group" %in% names(data$curve)
  out <- list()
  if (nrow(data$zones)) {
    out <- c(out, list(ggplot2::geom_rect(
      data = data$zones,
      mapping = ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax),
      ymin = -Inf, ymax = Inf, fill = "#b2182b", alpha = 0.08,
      inherit.aes = FALSE)))
  }
  out <- c(out, list(ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                         colour = "grey40")))
  if (grouped) {
    out <- c(out, list(
      ggplot2::geom_ribbon(data = data$curve,
                           mapping = ggplot2::aes(x = .data$x,
                                                  ymin = .data$lower,
                                                  ymax = .data$upper,
                                                  fill = .data$group),
                           alpha = 0.18),
      ggplot2::geom_line(data = data$curve,
                         mapping = ggplot2::aes(x = .data$x, y = .data$fit,
                                                colour = .data$group),
                         linewidth = 0.9)))
  } else {
    out <- c(out, list(
      ggplot2::geom_ribbon(data = data$curve,
                           mapping = ggplot2::aes(x = .data$x,
                                                  ymin = .data$lower,
                                                  ymax = .data$upper),
                           fill = "#1f5fa8", alpha = 0.18),
      ggplot2::geom_line(data = data$curve,
                         mapping = ggplot2::aes(x = .data$x, y = .data$fit),
                         colour = "#1f5fa8", linewidth = 0.9)))
  }
  out
}


#' S7 Class of the Surface Graph
#'
#' @description
#' The class that [plot_surface()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,mask,slices,level,vcov_type,n_grid See [plot_surface()].
#'
#' @return A `PlotSurface` object, inheriting from [term_plot].
#'
#' @name PlotSurface
#' @keywords internal
PlotSurface <- S7::new_class("PlotSurface", parent = term_plot,
  properties = list(
    mask = S7::class_logical,
    slices = S7::class_any,
    level = S7::class_numeric,
    vcov_type = S7::class_any,
    n_grid = S7::class_numeric
  ))


#' Graph of the Surface of a Smooth of Two Variables
#'
#' @description
#' Selects the graph of the effect of a smooth of two numeric variables
#' (a tensor product [modelterms7::te()] or `s(x1, x2)`) as a raster with
#' contour lines over the ranges of the two variables. With `slices`, the
#' graph instead shows curves of the effect against the first variable at
#' fixed values of the second, each with a pointwise band.
#'
#' @details
#' The effect is \eqn{f(x_1, x_2) = X(x_1, x_2) b}, evaluated on a regular
#' grid of `n_grid` points for each variable, with a standard error
#' \eqn{\sqrt{X V X'}} from the block of [vcov.StatmodFit()] for the
#' coefficients of the term. The raster is colored on a scale that is
#' symmetric about zero. With `mask = TRUE`, a point is drawn over every cell
#' in which \eqn{|f| < z_{(1+\text{level})/2}\,\mathrm{se}}, that is, where the
#' pointwise interval of the effect contains zero. The slices are the same
#' function at the values `slices` of the second variable, and the bands are
#' \eqn{f \pm z\,\mathrm{se}}. A smooth with a `by` factor is drawn in one
#' panel for each level.
#'
#' @param mask Whether to mark the cells in which the interval of the effect
#'   contains zero.
#' @param slices `NULL`, or a numeric vector of values of the second variable
#'   at which curves of the effect against the first variable are drawn.
#' @param level The confidence level.
#' @param vcov_type A [vcov_method] object, the variance of the effect.
#' @param n_grid The number of points of the grid of each variable.
#'
#' @return A `PlotSurface` object, inheriting from [term_plot].
#'
#' @seealso [plot_effect()], [plot_derivative()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x1 = runif(200), x2 = runif(200))
#' dd$y <- sin(2 * pi * dd$x1) * dd$x2 + rnorm(200, sd = 0.2)
#' fit <- statmod(y ~ te(x1, x2), distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, term = "te(x1, x2)", type = plot_surface(n_grid = 30),
#'                   data = dd)
#' ggplot2::autoplot(fit, term = "te(x1, x2)",
#'                   type = plot_surface(slices = c(0.2, 0.8)), data = dd)
#' @export
plot_surface <- function(mask = TRUE, slices = NULL, level = 0.95,
                         vcov_type = bayesian_vcov(), n_grid = 60) {
  if (!is.logical(mask) || length(mask) != 1L || is.na(mask)) {
    stop("'mask' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(slices) && (!is.numeric(slices) || !length(slices) ||
                           any(!is.finite(slices)))) {
    stop("'slices' must be NULL or a numeric vector of finite values.",
         call. = FALSE)
  }
  check_level_grid(level, n_grid)
  vcov_type <- as_vcov_method(vcov_type)
  PlotSurface(label = "plot_surface()", mask = mask, slices = slices,
              level = level, vcov_type = vcov_type, n_grid = n_grid)
}


S7::method(plot_data, list(PlotSurface, modelterms7::SmoothTerm)) <-
  function(type, x, fit, ctx, ...) {
    vars <- smooth_variables(x)
    if (length(vars) != 2L) {
      stop(sprintf(paste0("the term '%s' is a smooth of %d variable%s; ",
                          "plot_surface() draws a function of two, and ",
                          "plot_effect() or plot_derivative() draws one of ",
                          "one."), ctx$name, length(vars),
                   if (length(vars) == 1L) "" else "s"), call. = FALSE)
    }
    vars_all <- unique(c(vars, term_variables(x, ctx$name)))
    data <- plot_needs_data(ctx, fit, type, vars = vars_all)
    for (v in vars) {
      if (!is.numeric(data[[v]]) || is.matrix(data[[v]])) {
        stop(sprintf("the surface needs numeric variables, and '%s' is not.",
                     v), call. = FALSE)
      }
    }
    z <- stats::qnorm((1 + type@level) / 2)
    rng <- lapply(vars, function(v) range(data[[v]], na.rm = TRUE))
    sliced <- !is.null(type@slices)
    g1 <- seq(rng[[1L]][1L], rng[[1L]][2L], length.out = type@n_grid)
    g2 <- if (sliced) sort(unique(type@slices)) else
      seq(rng[[2L]][1L], rng[[2L]][2L], length.out = type@n_grid)
    if (sliced && (any(g2 < rng[[2L]][1L]) || any(g2 > rng[[2L]][2L]))) {
      stop(sprintf(paste0("'slices' must lie in the range of '%s' in the ",
                          "data, [%s, %s]."), vars[[2L]],
                   format(rng[[2L]][1L]), format(rng[[2L]][2L])),
           call. = FALSE)
    }
    grid <- expand.grid(g1, g2)
    names(grid) <- vars
    others <- setdiff(vars_all, vars)
    group <- NULL
    for (o in others) {
      vo <- data[[o]]
      if (!is.numeric(vo) && is.null(group) &&
          nlevels(as.factor(vo)) <= 12L) group <- o
    }
    if (!is.null(group)) {
      f <- as.factor(data[[group]])
      lev <- factor(levels(f), levels = levels(f))
      grid <- grid[rep(seq_len(nrow(grid)), length(lev)), , drop = FALSE]
      grid[[group]] <- rep(lev, each = nrow(grid) / length(lev))
    }
    for (o in setdiff(others, group)) {
      vo <- data[[o]]
      grid[[o]] <- if (is.matrix(vo)) {
        matrix(apply(vo, 2, stats::median, na.rm = TRUE), nrow(grid), ncol(vo),
               byrow = TRUE)
      } else if (is.numeric(vo)) stats::median(vo, na.rm = TRUE) else {
        f <- as.factor(vo)
        factor(levels(f)[1L], levels = levels(f))
      }
    }
    cc <- term_contribution(fit, ctx$param, ctx$name, x, ctx$cols, grid,
                            type@vcov_type)
    cur <- data.frame(x = grid[[vars[[1L]]]], y = grid[[vars[[2L]]]],
                      fit = cc$fit, se = cc$se)
    if (!is.null(group)) cur$group <- grid[[group]]
    ttl <- sprintf("%s: %s", ctx$param,
                   if (is.null(ctx$title)) ctx$name else ctx$title)
    if (sliced) {
      cur$lower <- cur$fit - z * cur$se
      cur$upper <- cur$fit + z * cur$se
      cur$slice <- factor(cur$y, levels = g2)
      out <- list(curve = cur)
      attr(out, "kind") <- "slices"
      attr(out, "labels") <- list(title = ttl, x = vars[[1L]],
                                  y = sprintf("effect on the predictor of %s",
                                              ctx$param),
                                  colour = vars[[2L]], fill = vars[[2L]])
    } else {
      cur$insignificant <- abs(cur$fit) < z * cur$se
      out <- list(surface = cur,
                  masked = if (type@mask) cur[cur$insignificant, , drop = FALSE]
                    else cur[0, , drop = FALSE])
      attr(out, "kind") <- "surface"
      attr(out, "labels") <- list(title = ttl, x = vars[[1L]], y = vars[[2L]],
                                  fill = "effect")
    }
    attr(out, "group") <- group
    out
  }


S7::method(plot_layers, PlotSurface) <- function(type, data, ...) {
  group <- attr(data, "group")
  if (identical(attr(data, "kind"), "slices")) {
    out <- list(
      ggplot2::geom_ribbon(data = data$curve,
                           mapping = ggplot2::aes(x = .data$x,
                                                  ymin = .data$lower,
                                                  ymax = .data$upper,
                                                  fill = .data$slice,
                                                  group = .data$slice),
                           alpha = 0.18),
      ggplot2::geom_line(data = data$curve,
                         mapping = ggplot2::aes(x = .data$x, y = .data$fit,
                                                colour = .data$slice,
                                                group = .data$slice),
                         linewidth = 0.9))
  } else {
    out <- list(
      ggplot2::geom_raster(data = data$surface,
                           mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                                  fill = .data$fit)),
      ggplot2::scale_fill_gradient2(low = "#1f5fa8", mid = "white",
                                    high = "#b2182b", midpoint = 0),
      ggplot2::geom_contour(data = data$surface,
                            mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                                   z = .data$fit),
                            colour = "grey20", linewidth = 0.25))
    if (nrow(data$masked)) {
      out <- c(out, list(ggplot2::geom_point(
        data = data$masked,
        mapping = ggplot2::aes(x = .data$x, y = .data$y), colour = "grey30",
        size = 0.3, alpha = 0.5)))
    }
  }
  if (!is.null(group)) {
    out <- c(out, list(ggplot2::facet_wrap(ggplot2::vars(.data$group))))
  }
  out
}


#' S7 Class of the Criterion Graph
#'
#' @description
#' The class that [plot_criterion()] instantiates. Use the constructor, which
#' checks the properties.
#'
#' @param label,log_lambda,method,hyper See [plot_criterion()].
#'
#' @return A `PlotCriterion` object, inheriting from [term_plot].
#'
#' @name PlotCriterion
#' @keywords internal
PlotCriterion <- S7::new_class("PlotCriterion", parent = term_plot,
  properties = list(
    log_lambda = S7::class_any,
    method = S7::class_character,
    hyper = S7::class_any
  ))


#' Graph of the Outer Criterion Against a Hyperparameter
#'
#' @description
#' Selects the graph of the outer criterion of the fit (the criterion that
#' estimated the hyperparameters: [reml()], [ml()], [aic()] or [bic()]) as a
#' function of the logarithm of one hyperparameter of a term, with the other
#' hyperparameters held at their fitted values. The fitted value is marked, and
#' the effective degrees of freedom of the term are drawn on a second axis.
#'
#' @details
#' The hyperparameter is a smoothing parameter, the precision \eqn{\lambda} of
#' a quadratic penalty, or the scale of a random effect, as the term names it
#' (`hyper` selects it when the term has several). The horizontal axis is
#' \eqn{\log \theta}, with \eqn{\theta} the value of the hyperparameter, so
#' that a hyperparameter whose value is not positive has no graph. The
#' default grid has 41 points within five units of the fitted logarithm.
#'
#' With `method = "refit"` the criterion is evaluated at each grid point
#' exactly: the model is refitted with the hyperparameter at that value, the
#' coefficients are warm-started at the fit, and the criterion is read at the
#' refitted mode, as the outer search itself reads it. The effective degrees
#' of freedom of the term are read at the same refit. The criterion that
#' [reml()] and [ml()] estimate also estimates the coefficients of an equation
#' for a dispersion (the `sigma` of a gaussian); these stay at their fitted
#' values. With `method = "quadratic"` the criterion is the second-order
#' expansion
#' \eqn{c(\eta_0) + g\,\delta - \tfrac{1}{2} A\,\delta^2} in the free
#' coordinate \eqn{\eta = h(\theta)}, where \eqn{h} is the link of the
#' hyperparameter, \eqn{\delta = \eta - \eta_0}, \eqn{g} is the exact outer
#' gradient and \eqn{-A} the exact diagonal element of the outer Hessian at the
#' fit; it needs [reml()] or [ml()] and shows no degrees of freedom. The
#' criterion is shown as the fit reports it (for [reml()] and [ml()] the log
#' marginal likelihood, to be maximized; for [aic()] and [bic()] the
#' criterion, to be minimized).
#'
#' @param log_lambda `NULL`, or a numeric vector of values of the logarithm
#'   of the hyperparameter.
#' @param method `"refit"` (the default) or `"quadratic"`.
#' @param hyper `NULL`, or the name of the hyperparameter of the term.
#'
#' @return A `PlotCriterion` object, inheriting from [term_plot].
#'
#' @seealso [plot_effect()], [hyper()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(120))
#' dd$y <- sin(2 * pi * dd$x) + rnorm(120, sd = 0.3)
#' fit <- statmod(y ~ s(x), distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, term = "s(x)",
#'                   type = plot_criterion(log_lambda = seq(-3, 8, length.out = 12)))
#' @export
plot_criterion <- function(log_lambda = NULL, method = "refit", hyper = NULL) {
  if (!is.null(log_lambda) && (!is.numeric(log_lambda) ||
                               length(log_lambda) < 2L ||
                               any(!is.finite(log_lambda)))) {
    stop("'log_lambda' must be NULL or a numeric vector of at least two finite values.",
         call. = FALSE)
  }
  method <- match.arg(method, c("refit", "quadratic"))
  if (!is.null(hyper) && (!is.character(hyper) || length(hyper) != 1L ||
                          is.na(hyper))) {
    stop("'hyper' must be NULL or a single string.", call. = FALSE)
  }
  PlotCriterion(label = "plot_criterion()", log_lambda = log_lambda,
                method = method, hyper = hyper)
}


#' The Data of a Criterion Graph
#'
#' @description
#' Computes the curve of the outer criterion against the logarithm of one
#' hyperparameter of a term, the fitted point and, for a refit, the effective
#' degrees of freedom of the term.
#'
#' @param type A `PlotCriterion` object.
#' @param x The term.
#' @param fit A [StatmodFit()].
#' @param ctx The context described at [plot_data()].
#' @param ... Unused.
#'
#' @return A list with the data frames `curve` and `fitted`, and the
#'   attributes `"axis"` and `"labels"`.
#'
#' @keywords internal
criterion_plot_data <- function(type, x, fit, ctx, ...) {
  cd <- criterion_setup(fit, ctx, type@hyper)
  eta0 <- cd$eta0[[cd$row]]
  th0 <- cd$value
  u0 <- log(th0)
  u <- if (is.null(type@log_lambda)) seq(u0 - 5, u0 + 5, length.out = 41L) else
    sort(unique(type@log_lambda))
  link <- attr(cd$idx, "links")[[cd$row]]
  eta <- vapply(u, function(v) linkfunctions7::linkfun(link, exp(v)),
                numeric(1))
  c0 <- fit@criterion
  if (identical(type@method, "quadratic")) {
    q <- criterion_quadratic(fit, cd)
    d <- eta - eta0
    val <- c0 + q$g * d - 0.5 * q$A * d^2
    edf <- rep(NA_real_, length(u))
    f0 <- c0
  } else {
    ev <- criterion_evaluator(fit, cd)
    at0 <- ev(cd$eta0)
    if (!is.finite(at0$value)) {
      stop(paste0("the criterion could not be evaluated at the fitted ",
                  "hyperparameters."), call. = FALSE)
    }
    val <- edf <- rep(NA_real_, length(u))
    for (i in seq_along(u)) {
      e <- cd$eta0
      e[[cd$row]] <- eta[[i]]
      r <- ev(e)
      val[i] <- r$value
      edf[i] <- r$edf
    }
    f0 <- at0$value
  }
  ok <- is.finite(val)
  if (!any(ok)) {
    stop("the criterion could not be evaluated at any value of the grid.",
         call. = FALSE)
  }
  curve <- data.frame(x = u, value = val, edf = edf)[ok, , drop = FALSE]
  at <- data.frame(x = u0, value = f0, edf = if (identical(type@method, "refit"))
    at0$edf else NA_real_)
  axis <- NULL
  if (any(is.finite(curve$edf)) && diff(range(curve$edf, na.rm = TRUE)) > 0) {
    b <- diff(range(curve$value)) / diff(range(curve$edf, na.rm = TRUE))
    a <- min(curve$value) - b * min(curve$edf, na.rm = TRUE)
    if (is.finite(b) && b != 0) axis <- c(a = a, b = b)
  }
  out <- list(curve = curve, fitted = at)
  attr(out, "axis") <- axis
  kind <- fit@methods$outer@kind
  attr(out, "labels") <- list(
    title = sprintf("%s: %s, %s", ctx$param,
                    if (is.null(ctx$title)) ctx$name else ctx$title,
                    cd$name),
    x = sprintf("log %s", cd$name),
    y = sprintf("%s criterion", toupper(kind)))
  out
}

for (.cls in list(modelterms7::SmoothTerm, modelterms7::PenalizedTerm,
                  modelterms7::RandomTerm, modelterms7::RandomSmoothTerm)) {
  S7::method(plot_data, list(PlotCriterion, .cls)) <- criterion_plot_data
}
rm(.cls)


#' The Hyperparameter of a Term for the Criterion Graph
#'
#' @description
#' Finds the row of the outer index that belongs to the hyperparameter drawn,
#' and gathers what the evaluation of the criterion needs.
#'
#' @details
#' An error is signaled where the term has no hyperparameter that the outer
#' criterion estimates (a held value, a lasso, scad or mcp penalty, which is
#' chosen along a path), where the fit ran no outer criterion, where the
#' hyperparameter is not positive, and where the model has a structural term.
#'
#' @param fit A [StatmodFit()].
#' @param ctx The context described at [plot_data()].
#' @param hyper `NULL`, or the name of the hyperparameter.
#'
#' @return A list with the design, the blocks, the outer index, the free
#'   coordinates `eta0` of every hyperparameter, the selected `row`, its
#'   `name` and its fitted `value`.
#'
#' @keywords internal
criterion_setup <- function(fit, ctx, hyper) {
  method <- fit@methods$outer
  if (is.null(method)) {
    stop(paste0("the fit ran no outer criterion, so its hyperparameters ",
                "were not estimated and there is no criterion to draw."),
         call. = FALSE)
  }
  spec <- fit@spec
  design <- statmod_design(spec)
  if (length(attr(design, "structural"))) {
    stop(paste0("the criterion graph is not available for a model with a ",
                "structural term (gas(), regime())."), call. = FALSE)
  }
  if (length(fit@methods$pinned)) attr(design, "pinned") <- fit@methods$pinned
  blocks <- statmod_blocks(spec, design)
  idx <- outer_hyper_index(spec, blocks)
  squash <- function(s) gsub("[[:space:]]", "", s)
  mem <- if (nrow(idx)) index_members(idx) else idx
  mine <- which(mem$parameter == ctx$param & squash(mem$term) == squash(ctx$name))
  if (!length(mine)) {
    stop(sprintf(paste0("the term '%s' has no hyperparameter that the outer ",
                        "criterion estimates: its hyperparameters are held ",
                        "at the values written on the term, or its penalty ",
                        "has a kink and is chosen along a path."), ctx$name),
         call. = FALSE)
  }
  if (!is.null(hyper)) {
    mine <- mine[mem$name[mine] == hyper]
    if (!length(mine)) {
      stop(sprintf("the term '%s' has no estimated hyperparameter '%s'.",
                   ctx$name, hyper), call. = FALSE)
    }
  } else if (length(unique(mem$row[mine])) > 1L) {
    stop(sprintf(paste0("the term '%s' has several hyperparameters, %s; ",
                        "name one with 'hyper'."), ctx$name,
                 paste0("'", mem$name[mine], "'", collapse = ", ")),
         call. = FALSE)
  }
  m <- mine[[1L]]
  row <- mem$row[[m]]
  value <- hyper_value(fit@hyper, idx, row)
  if (!is.finite(value) || value <= 0) {
    stop(sprintf(paste0("the hyperparameter '%s' of '%s' is %s, and the graph ",
                        "is drawn against its logarithm."), mem$name[[m]],
                 ctx$name, format(value)), call. = FALSE)
  }
  list(design = design, blocks = blocks, idx = idx,
       eta0 = hyper_to_eta(fit@hyper, idx), row = row, name = mem$name[[m]],
       value = value, param = ctx$param, term = ctx$name)
}


#' The Specification with the Criterion's Coefficients Held
#'
#' @description
#' Returns the specification in which the coefficients that the outer
#' criterion estimates beside the hyperparameters (the dispersion of the
#' marginal criteria) are held at their fitted values, together with the
#' description of those coordinates.
#'
#' @param fit A [StatmodFit()].
#' @param design The design.
#'
#' @return A list with `spec` and `gam`.
#'
#' @keywords internal
criterion_held_spec <- function(fit, design) {
  spec <- fit@spec
  method <- fit@methods$outer
  gam <- if (outer_minimize(method)) marginal_coords(spec, design, NULL) else
    marginal_coords(spec, design, method)
  if (length(gam$where)) {
    beta <- unlist(fit@coefficients[spec@distrib@params], use.names = FALSE)
    spec <- statmod_hold(spec, gam, beta[gam$where])
  }
  list(spec = spec, gam = gam)
}


#' Evaluate the Outer Criterion at Given Hyperparameters
#'
#' @description
#' Returns a function of the free coordinates of the hyperparameters that
#' refits the model there and returns the criterion and the effective degrees
#' of freedom of the term.
#'
#' @details
#' The refit is the one of `criterion_differenced()`: the inner fit is
#' warm-started at the coefficients of `fit`, with the inner method, the
#' budget and the information of the fit, and the criterion is read at the
#' refitted mode by the function the outer search uses
#' ([statmod_marginal()], or [statmod_pe()] for [aic()] and [bic()]). The
#' value is `NA` where the refit gives no usable criterion.
#'
#' @param fit A [StatmodFit()].
#' @param cd The list returned by `criterion_setup()`.
#'
#' @return A function of the free coordinates that returns a list with
#'   `value` and `edf`.
#'
#' @keywords internal
criterion_evaluator <- function(fit, cd) {
  design <- cd$design
  hs <- criterion_held_spec(fit, design)
  spec <- hs$spec
  gam <- hs$gam
  method <- fit@methods$outer
  pe <- outer_minimize(method)
  io <- fit@methods$smooth
  budget <- method_budget(io)
  cfg <- inner_settings(io)
  obj <- statmod_objective(spec, fit@hyper, design, cfg$expected, cfg$emethod)
  warm <- obj$stack(fit@coefficients)
  basis <- if (pe) NULL else
    integrated_basis(spec, design, method@kind, gamma = length(gam$where) > 0L)
  info <- fit@methods$edf_information
  edf_exp <- S7::S7_inherits(info, distributions7::ExpectedInfo)
  edf_em <- if (edf_exp) info@method else NULL
  function(eta) {
    hy <- eta_to_hyper(eta, cd$idx, fit@hyper)
    res <- tryCatch(
      statmod_alternate(spec, design, cd$blocks, hy, io, warm, cfg$expected,
                        cfg$emethod, budget$maxit, budget$tol, verbosity(0),
                        hold_refresh = TRUE),
      error = function(e) NULL)
    bad <- list(value = NA_real_, edf = NA_real_)
    if (is.null(res)) return(bad)
    cf <- res$obj$split(res$par)
    m <- tryCatch(
      if (pe) {
        statmod_pe(spec, design, cf, hy, method, cfg$emethod,
                   statmod_active(spec, cd$blocks, res$par, hy))
      } else {
        ctx <- outer_context(spec, design, cf, hy, cfg$emethod)
        statmod_marginal(spec, design, cf, hy, method, cfg$emethod, basis, ctx)
      }, error = function(e) NULL)
    if (is.null(m) || !is.finite(m$value)) return(bad)
    ed <- tryCatch(statmod_edf(spec, cf, design, hy, edf_exp, edf_em),
                   error = function(e) NULL)
    e <- NA_real_
    if (!is.null(ed)) {
      v <- ed$edf[ed$parameter == cd$param & ed$term == cd$term]
      if (length(v)) e <- v[[1L]]
    }
    list(value = m$value, edf = e)
  }
}


#' The Quadratic Expansion of the Outer Criterion
#'
#' @description
#' Reads the exact outer gradient and the exact outer Hessian of a marginal
#' criterion at the fit and returns the entries that belong to one
#' hyperparameter.
#'
#' @param fit A [StatmodFit()].
#' @param cd The list returned by `criterion_setup()`.
#'
#' @return A list with the gradient `g` and the curvature `A` (minus the
#'   Hessian), both on the free coordinate of the hyperparameter.
#'
#' @keywords internal
criterion_quadratic <- function(fit, cd) {
  method <- fit@methods$outer
  if (!method@kind %in% c("reml", "ml")) {
    stop(paste0("method = \"quadratic\" needs the exact outer gradient and ",
                "Hessian of reml() or ml(), and the fit used ", method@kind,
                "(); use method = \"refit\"."), call. = FALSE)
  }
  design <- cd$design
  hs <- criterion_held_spec(fit, design)
  spec <- hs$spec
  gam <- hs$gam
  ng <- length(gam$where)
  basis <- integrated_basis(spec, design, method@kind, gamma = ng > 0L)
  cf <- fit@coefficients
  hy <- fit@hyper
  ctx <- outer_context(spec, design, cf, hy, NULL)
  g <- tryCatch(statmod_marginal_grad(spec, design, cf, hy, method, cd$idx,
                                      basis, ctx = ctx, gam = gam),
                error = function(e) NULL)
  cur <- outer_curvature(spec, design, cf, hy, method, cd$idx, basis, NULL, gam)
  if (is.null(g) || is.null(cur$A) || !is.finite(g[[cd$row]]) ||
      !is.finite(cur$A[cd$row, cd$row])) {
    stop(paste0("the exact outer gradient and Hessian could not be read at the ",
                "fit; use method = \"refit\"."), call. = FALSE)
  }
  list(g = g[[cd$row]], A = cur$A[cd$row, cd$row])
}


S7::method(plot_layers, PlotCriterion) <- function(type, data, ...) {
  axis <- attr(data, "axis")
  cv <- data$curve
  out <- list()
  if (!is.null(axis)) {
    cv$edf_scaled <- axis[["a"]] + axis[["b"]] * cv$edf
    a <- axis[["a"]]
    b <- axis[["b"]]
    out <- c(out, list(
      ggplot2::geom_line(data = cv[is.finite(cv$edf_scaled), , drop = FALSE],
                         mapping = ggplot2::aes(x = .data$x,
                                                y = .data$edf_scaled),
                         colour = "grey55", linetype = "dashed"),
      ggplot2::scale_y_continuous(
        sec.axis = ggplot2::sec_axis(~ (. - a) / b, name = "edf of the term"))))
  }
  c(out, list(
    ggplot2::geom_vline(xintercept = data$fitted$x, colour = "#b2182b",
                        alpha = 0.4, linetype = "dotted"),
    ggplot2::geom_line(data = cv, mapping = ggplot2::aes(x = .data$x,
                                                         y = .data$value),
                       colour = "#1f5fa8", linewidth = 0.9),
    ggplot2::geom_point(data = cv, mapping = ggplot2::aes(x = .data$x,
                                                          y = .data$value),
                        colour = "#1f5fa8", size = 1.1),
    ggplot2::geom_point(data = data$fitted,
                        mapping = ggplot2::aes(x = .data$x, y = .data$value),
                        colour = "#b2182b", size = 2.8)))
}
