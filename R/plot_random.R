#' @title S7 Classes of the Graphs of a Random Effect
#'
#' @description
#' The classes that [plot_caterpillar()], [plot_prior()], [plot_random_cov()]
#' and [plot_curves()] instantiate. Use the constructors, which check the
#' properties.
#'
#' @param label,raw,sort,level,vcov_type,information,which,n_grid,population
#'   See the constructors.
#'
#' @return An object inheriting from [term_plot].
#'
#' @name PlotCaterpillar
#' @aliases PlotPrior PlotRandomCov PlotCurves
#' @keywords internal
PlotCaterpillar <- S7::new_class("PlotCaterpillar", parent = term_plot,
  properties = list(
    raw = S7::class_logical,
    sort = S7::class_logical,
    level = S7::class_numeric,
    vcov_type = S7::class_any,
    information = S7::class_any
  ))

#' @rdname PlotCaterpillar
#' @keywords internal
PlotPrior <- S7::new_class("PlotPrior", parent = term_plot)

#' @rdname PlotCaterpillar
#' @keywords internal
PlotRandomCov <- S7::new_class("PlotRandomCov", parent = term_plot,
  properties = list(
    level = S7::class_numeric,
    which = S7::class_numeric
  ))

#' @rdname PlotCaterpillar
#' @keywords internal
PlotCurves <- S7::new_class("PlotCurves", parent = term_plot,
  properties = list(
    n_grid = S7::class_numeric,
    population = S7::class_logical
  ))


#' Graph of the Predicted Effects of a Random Effect
#'
#' @description
#' Selects the graph with one row for each level of the grouping factor of a
#' [random()] term: the predicted effect with its confidence interval and,
#' optionally, the unshrunk estimate of the level as an open circle. A term
#' with several coefficients for each level, such as a random slope, gives one
#' panel for each coefficient.
#'
#' @details
#' The predicted effect is the coefficient of the level, and the interval is
#' \eqn{\hat b \pm z_{(1+\text{level})/2}\,\sqrt{V_{jj}}}, where \eqn{V} is the
#' variance matrix of the coefficients selected by `vcov_type`. The unshrunk
#' estimate of a level is the least squares estimate of the coefficients of
#' that level from the partial residuals
#' \eqn{Z_i \hat b_j + (z_i - \eta_i)}, over the observations of the level,
#' where \eqn{Z_i} is the row of the design of the term, \eqn{\hat b_j} the
#' predicted effect of the level, and \eqn{z_i - \eta_i} the working residual
#' of the equation (see [plot_effect()]). With one coefficient for each level
#' it is the mean of the partial residuals of the level, and for a gaussian
#' random intercept it is the mean of the residuals of the level plus its
#' predicted effect. A level with too few observations to estimate all its
#' coefficients has no open circle.
#'
#' @param raw Whether to draw the unshrunk estimates.
#' @param sort Whether to order the levels by the predicted effect of the
#'   first coefficient.
#' @param level The confidence level of the intervals.
#' @param vcov_type A [vcov_method] object, the variance of the intervals.
#' @param information An object from [distributions7::expected_info()] (the
#'   default) or [distributions7::observed_info()], the weight of the working
#'   residuals.
#'
#' @return A `PlotCaterpillar` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_prior()], [plot_random_cov()],
#'   [plot_shrinkage()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(120), id = factor(sample(1:12, 120, TRUE)))
#' dd$y <- 1 + dd$x + rnorm(12)[dd$id] + rnorm(120, sd = 0.5)
#' fit <- statmod(y ~ x + random(~ 1 | id), distributions7::gaussian1_distrib(),
#'                dd)
#' ggplot2::autoplot(fit, term = "random(~1|id)", plot_type = plot_caterpillar())
#' @export
plot_caterpillar <- function(raw = TRUE, sort = TRUE, level = 0.95,
                             vcov_type = bayesian_vcov(),
                             information = distributions7::expected_info()) {
  if (!is.logical(raw) || length(raw) != 1L || is.na(raw)) {
    stop("'raw' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(sort) || length(sort) != 1L || is.na(sort)) {
    stop("'sort' must be TRUE or FALSE.", call. = FALSE)
  }
  check_level_grid(level, 2)
  vcov_type <- as_vcov_method(vcov_type)
  check_plot_information(information)
  PlotCaterpillar(label = "plot_caterpillar()", raw = raw, sort = sort,
                  level = level, vcov_type = vcov_type,
                  information = information)
}


#' Graph of the Prior of a Random Effect
#'
#' @description
#' Selects the graph that compares the predicted effects of a [random()] term
#' with the prior density that the fit estimated for them: a histogram of the
#' predicted effects on the density scale, with the density on top. A term with
#' several coefficients for each level gives one panel for each coefficient.
#'
#' @details
#' The density is the marginal density of the coefficient under the
#' distribution of the effects: the gaussian density with the fitted standard
#' deviation for the default prior, and the density of the distribution passed
#' to `random()` (a Student t, for example) with its fitted hyperparameters for
#' a univariate prior. For a multivariate gaussian prior the marginal is the
#' gaussian density with the fitted standard deviation of the coefficient. The
#' predicted effects are shrunk towards zero, so their histogram is narrower
#' than the prior.
#'
#' @return A `PlotPrior` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_caterpillar()],
#'   [plot_random_cov()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200), id = factor(sample(1:25, 200, TRUE)))
#' dd$y <- 1 + dd$x + rnorm(25)[dd$id] + rnorm(200, sd = 0.5)
#' fit <- statmod(y ~ x + random(~ 1 | id), distributions7::gaussian1_distrib(),
#'                dd)
#' ggplot2::autoplot(fit, term = "random(~1|id)", plot_type = plot_prior())
#' @export
plot_prior <- function() {
  PlotPrior(label = "plot_prior()")
}


#' Graph of the Covariance of a Random Effect
#'
#' @description
#' Selects the graph of the predicted effects of two coefficients of a
#' [random()] term, such as an intercept and a slope, one point for each level,
#' with the ellipse of the fitted covariance of the effects.
#'
#' @details
#' The ellipse is the set \eqn{\{b : b' \Sigma^{-1} b = c\}}, centred at zero,
#' where \eqn{\Sigma} is the fitted covariance matrix of the two effects, built
#' from their fitted standard deviation and correlation, and
#' \eqn{c} is the quantile `level` of the chi-squared distribution with two
#' degrees of freedom. It holds the fraction `level` of the effects when the
#' prior is gaussian, which is the only prior for which the graph is defined.
#' The predicted effects are shrunk, so the points lie inside the ellipse more
#' often than that fraction.
#'
#' @param level The probability content of the ellipse.
#' @param which The two coefficients within a level, by position.
#'
#' @return A `PlotRandomCov` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_caterpillar()], [plot_prior()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(300), id = factor(sample(1:20, 300, TRUE)))
#' dd$y <- 1 + dd$x + rnorm(20)[dd$id] + 0.5 * rnorm(20)[dd$id] * dd$x +
#'   rnorm(300, sd = 0.5)
#' fit <- statmod(y ~ x + random(~ 1 + x | id), distributions7::gaussian1_distrib(),
#'                dd)
#' ggplot2::autoplot(fit, term = "random(~1+x|id)", plot_type = plot_random_cov())
#' @export
plot_random_cov <- function(level = 0.95, which = c(1, 2)) {
  check_level_grid(level, 2)
  if (!is.numeric(which) || length(which) != 2L || any(which < 1) ||
      any(which != round(which)) || which[1L] == which[2L]) {
    stop("'which' must be two different positions of coefficients.",
         call. = FALSE)
  }
  PlotRandomCov(label = "plot_random_cov()", level = level, which = which)
}


#' Graph of the Curves of a Random Smooth
#'
#' @description
#' Selects the graph of a random smooth, `random(~ s(x) | g)`: one thin curve
#' for each level of the grouping factor, over the range of the variable
#' within that level. The graph reads the covariates and requires `data`.
#'
#' @details
#' The curve of a level is its contribution to the predictor of the equation,
#' \eqn{g_j(x) = a_{0j} + N(x) a_j + Z(x) \delta_j}, the departure of the
#' level from the curve that the other terms of the equation describe. It is
#' centred around zero and does not include the population curve. With
#' `population = TRUE` the population curve, taken from the single smooth of
#' the same variable in the same equation, is added to each curve and drawn
#' as a thick line, so that the curves are the curves of the levels.
#'
#' @param n_grid The number of points of each curve.
#' @param population Whether to add the population smooth of the variable.
#'
#' @return A `PlotCurves` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_caterpillar()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(240), id = factor(rep(1:8, each = 30)))
#' dd$y <- sin(2 * pi * dd$x) + rnorm(8, sd = 0.4)[dd$id] * dd$x +
#'   rnorm(240, sd = 0.3)
#' fit <- statmod(y ~ s(x) + random(~ s(x) | id),
#'                distributions7::gaussian1_distrib(), dd)
#' ggplot2::autoplot(fit, term = "random(~s(x)|id)", plot_type = plot_curves(),
#'                   data = dd)
#' @export
plot_curves <- function(n_grid = 100, population = FALSE) {
  if (!is.numeric(n_grid) || length(n_grid) != 1L || !is.finite(n_grid) ||
      n_grid < 2 || n_grid != round(n_grid)) {
    stop("'n_grid' must be a single whole number of at least 2.", call. = FALSE)
  }
  if (!is.logical(population) || length(population) != 1L || is.na(population)) {
    stop("'population' must be TRUE or FALSE.", call. = FALSE)
  }
  PlotCurves(label = "plot_curves()", n_grid = n_grid, population = population)
}


S7::method(default_plot, modelterms7::RandomTerm) <- function(term, ...) {
  plot_caterpillar()
}

S7::method(default_plot, modelterms7::RandomSmoothTerm) <- function(term, ...) {
  plot_curves()
}


#' The Coefficients of a Random Term by Level
#'
#' @description
#' Returns the coefficients of a random term as a matrix with one row for each
#' level of the grouping factor and one column for each coefficient within the
#' level, with their positions in the stacked coefficients of the model.
#'
#' @param fit A [StatmodFit()].
#' @param param,name The equation and the name of the term.
#' @param tm The [modelterms7::RandomTerm].
#'
#' @return A list with `levels`, `wnames`, `idx` (positions within the
#'   equation's coefficients, group by group), `coef_names` (their names in
#'   [vcov.StatmodFit()]) and `b`, the matrix.
#'
#' @keywords internal
random_block <- function(fit, param, name, tm) {
  bp <- tm@build_info
  des <- statmod_design(fit@spec)
  idx <- des[[param]]$blocks[[name]]
  glev <- bp$glevels
  m <- length(glev)
  d <- length(idx) / m
  if (is.null(glev) || d != round(d)) {
    stop(sprintf("the term '%s' does not have a coefficient set for each level.",
                 name), call. = FALSE)
  }
  wn <- bp$wnames
  if (is.null(wn) || length(wn) != d) wn <- paste0("effect", seq_len(d))
  b <- fit@coefficients[[param]][idx]
  list(levels = glev, wnames = wn, idx = idx,
       coef_names = paste0(param, ":", des[[param]]$coef_names[idx]),
       b = matrix(b, nrow = m, ncol = d, byrow = TRUE))
}


#' Reject a Random Smooth in a Graph of a Random Effect
#'
#' @param type The [plot_type] object.
#' @param x The term.
#' @param ctx The context of the graph.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
check_plain_random <- function(type, x, ctx) {
  if (S7::S7_inherits(x, modelterms7::RandomSmoothTerm)) {
    stop(sprintf(paste0("the graph %s is not available for the term '%s'; ",
                        "plot_curves() draws the curves of a random smooth."),
                 type@label, ctx$name), call. = FALSE)
  }
  invisible(NULL)
}


S7::method(plot_data, list(PlotCaterpillar, modelterms7::RandomTerm)) <-
  function(type, x, fit, ctx, ...) {
    check_plain_random(type, x, ctx)
    rb <- random_block(fit, ctx$param, ctx$name, x)
    m <- nrow(rb$b)
    d <- ncol(rb$b)
    V <- vcov(fit, vcov_type = type@vcov_type, readable = FALSE)
    se <- sqrt(pmax(diag(V[rb$coef_names, rb$coef_names, drop = FALSE]), 0))
    se <- matrix(se, nrow = m, ncol = d, byrow = TRUE)
    z <- stats::qnorm((1 + type@level) / 2)
    raw <- matrix(NA_real_, m, d)
    if (type@raw) raw <- raw_effects(fit, ctx$param, x, rb, type@information)
    est <- data.frame(
      level = rep(rb$levels, times = d),
      coef = factor(rep(rb$wnames, each = m), levels = rb$wnames),
      estimate = as.vector(rb$b), se = as.vector(se),
      raw = as.vector(raw), stringsAsFactors = FALSE)
    est$lower <- est$estimate - z * est$se
    est$upper <- est$estimate + z * est$se
    ord <- if (type@sort) rb$levels[order(rb$b[, 1L])] else rb$levels
    est$level <- factor(est$level, levels = ord)
    out <- list(effects = est)
    attr(out, "labels") <- list(
      title = sprintf("%s: %s", ctx$param, if (is.null(ctx$title)) ctx$name
                      else ctx$title),
      x = "predicted effect", y = NULL,
      caption = if (type@raw) "open circles: unshrunk estimates" else NULL)
    out
  }


#' Unshrunk Estimates of the Effects of a Random Term
#'
#' @description
#' For each level, the least squares estimate of its coefficients from the
#' partial residuals of its observations, weighted by the prior weights.
#'
#' @param fit A [StatmodFit()].
#' @param param The equation.
#' @param tm The random term.
#' @param rb The list from [random_block()].
#' @param information An object from [distributions7::expected_info()] or
#'   [distributions7::observed_info()].
#'
#' @return A matrix with the shape of `rb$b`; `NA` where a level cannot
#'   identify its coefficients.
#'
#' @keywords internal
raw_effects <- function(fit, param, tm, rb, information) {
  Z <- as.matrix(modelterms7::term_matrix(tm))
  wr <- working_residuals(fit, param, information)
  w <- fit@spec@weights
  m <- nrow(rb$b)
  d <- ncol(rb$b)
  out <- matrix(NA_real_, m, d)
  for (j in seq_len(m)) {
    cols <- (j - 1L) * d + seq_len(d)
    Zj <- Z[, cols, drop = FALSE]
    rows <- which(rowSums(abs(Zj)) > 0 & is.finite(wr) & w > 0)
    if (length(rows) < d) next
    Zr <- Zj[rows, , drop = FALSE]
    wt <- w[rows]
    A <- crossprod(Zr * sqrt(wt))
    if (!is.finite(rc <- rcond(A)) || rc < 1e-10) next
    out[j, ] <- rb$b[j, ] + drop(solve(A, crossprod(Zr, wt * wr[rows])))
  }
  out
}


S7::method(plot_layers, PlotCaterpillar) <- function(type, data, ...) {
  d <- data$effects
  lv <- levels(d$level)
  keep <- if (length(lv) > 40L) lv[round(seq(1, length(lv), length.out = 40L))]
    else lv
  out <- list(
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50"),
    ggplot2::geom_linerange(
      data = d, mapping = ggplot2::aes(y = .data$level, xmin = .data$lower,
                                       xmax = .data$upper),
      orientation = "y", colour = "#1f5fa8"),
    ggplot2::geom_point(
      data = d, mapping = ggplot2::aes(x = .data$estimate, y = .data$level),
      colour = "#1f5fa8", size = 1.8))
  if (type@raw && any(is.finite(d$raw))) {
    out <- c(out, list(ggplot2::geom_point(
      data = d[is.finite(d$raw), , drop = FALSE],
      mapping = ggplot2::aes(x = .data$raw, y = .data$level),
      shape = 1, colour = "grey30", size = 1.8)))
  }
  out <- c(out, list(ggplot2::scale_y_discrete(breaks = keep)))
  if (nlevels(d$coef) > 1L) {
    out <- c(out, list(ggplot2::facet_wrap(ggplot2::vars(.data$coef),
                                           scales = "free_x")))
  }
  out
}


#' The Fitted Distribution of the Effects of a Random Term
#'
#' @description
#' Reads the penalty of the term and the fitted hyperparameters. Returns the
#' marginal density of each coefficient within a level and, for a gaussian
#' prior, the fitted covariance matrix of the coefficients.
#'
#' @param fit A [StatmodFit()].
#' @param param,name The equation and the name of the term.
#' @param tm The random term.
#'
#' @return A list with `wnames`, `dens` (a list of functions, one for each
#'   coefficient within a level, `NULL` where the density is not available),
#'   `sd` (the standard deviation of the gaussian prior of each coefficient,
#'   `NA` otherwise) and `Sigma` (the covariance matrix, `NULL` unless the
#'   prior is gaussian).
#'
#' @keywords internal
random_prior_fit <- function(fit, param, name, tm) {
  if (!is.na(modelterms7::term_tag(tm))) {
    stop(sprintf(paste0("the term '%s' shares its covariance with other ",
                        "terms, and its prior is not drawn by itself."),
                 name), call. = FALSE)
  }
  bp <- tm@build_info
  m <- length(bp$glevels)
  d <- length(fit@coefficients[[param]][
    statmod_design(fit@spec)[[param]]$blocks[[name]]]) / m
  wn <- bp$wnames
  if (is.null(wn) || length(wn) != d) wn <- paste0("effect", seq_len(d))
  units <- statmod_penalized(fit@spec, statmod_design(fit@spec))
  units <- Filter(function(u) unit_of_term(u, param, name), units)
  dens <- vector("list", d)
  sds <- rep(NA_real_, d)
  Sigma <- matrix(0, d, d)
  gauss <- rep(FALSE, d)
  full_cov <- NULL
  for (u in units) {
    pen <- u$penalty
    th <- as.list(fit@hyper[[u$param]][[u$key]])
    cs <- sort(unique(((u$within - 1L) %% d) + 1L))
    blk <- tryCatch(as.integer(pen@block), error = function(e) NA_integer_)
    is_gauss <- grepl("gaussian", pen@parent@distrib_name, fixed = TRUE)
    if (identical(blk, 1L)) {
      for (k in cs) {
        local({
          parent <- pen@parent
          theta <- th
          dens[[k]] <<- function(x) distributions7::distrib_pdf(parent, x, theta)
        })
        if (is_gauss && !is.null(th[["sigma"]])) {
          sds[k] <- th[["sigma"]]
          Sigma[k, k] <- th[["sigma"]]^2
          gauss[k] <- TRUE
        }
      }
    } else if (is_gauss && identical(blk, length(cs))) {
      rd <- penalties7::penalty_readable(pen, th)
      sdv <- rd$value[grep("^sd_", names(rd$value))]
      R <- diag(length(cs))
      cv <- rd$value[grep("^cor_", names(rd$value))]
      for (nm in names(cv)) {
        ij <- as.integer(strsplit(sub("^cor_", "", gsub("v", "", nm)), "_")[[1L]])
        R[ij[1L], ij[2L]] <- R[ij[2L], ij[1L]] <- cv[[nm]]
      }
      S <- diag(sdv, length(sdv)) %*% R %*% diag(sdv, length(sdv))
      Sigma[cs, cs] <- S
      for (i in seq_along(cs)) {
        local({
          s <- sdv[[i]]
          dens[[cs[i]]] <<- function(x) stats::dnorm(x, 0, s)
        })
        sds[cs[i]] <- sdv[[i]]
      }
      gauss[cs] <- TRUE
    }
  }
  list(wnames = wn, dens = dens, sd = sds,
       Sigma = if (all(gauss)) Sigma else NULL)
}


S7::method(plot_data, list(PlotPrior, modelterms7::RandomTerm)) <-
  function(type, x, fit, ctx, ...) {
    check_plain_random(type, x, ctx)
    rb <- random_block(fit, ctx$param, ctx$name, x)
    pr <- random_prior_fit(fit, ctx$param, ctx$name, x)
    d <- ncol(rb$b)
    missing <- which(vapply(pr$dens, is.null, logical(1)))
    if (length(missing)) {
      stop(sprintf(paste0("the prior density of the coefficient '%s' of the ",
                          "term '%s' is not available."),
                   rb$wnames[missing[1L]], ctx$name), call. = FALSE)
    }
    eff <- data.frame(
      coef = factor(rep(rb$wnames, each = nrow(rb$b)), levels = rb$wnames),
      value = as.vector(rb$b))
    dd <- lapply(seq_len(d), function(k) {
      v <- rb$b[, k]
      lim <- range(v)
      if (is.finite(pr$sd[k])) lim <- range(lim, c(-4, 4) * pr$sd[k])
      else lim <- lim + c(-1, 1) * diff(lim)
      xs <- seq(lim[1L], lim[2L], length.out = 400L)
      data.frame(coef = factor(rb$wnames[k], levels = rb$wnames), x = xs,
                 density = pr$dens[[k]](xs))
    })
    out <- list(effects = eff, density = do.call(rbind, dd))
    attr(out, "bins") <- min(30L, max(8L, ceiling(2 * sqrt(nrow(rb$b)))))
    attr(out, "labels") <- list(
      title = sprintf("%s: %s", ctx$param, if (is.null(ctx$title)) ctx$name
                      else ctx$title),
      x = "effect", y = "density",
      caption = "bars: predicted effects; line: fitted prior")
    out
  }


S7::method(plot_layers, PlotPrior) <- function(type, data, ...) {
  out <- list(
    ggplot2::geom_histogram(
      data = data$effects,
      mapping = ggplot2::aes(x = .data$value, y = ggplot2::after_stat(density)),
      bins = attr(data, "bins"), fill = "grey80", colour = "white"),
    ggplot2::geom_line(data = data$density,
                       mapping = ggplot2::aes(x = .data$x, y = .data$density),
                       colour = "#1f5fa8", linewidth = 0.9))
  if (nlevels(data$effects$coef) > 1L) {
    out <- c(out, list(ggplot2::facet_wrap(ggplot2::vars(.data$coef),
                                           scales = "free")))
  }
  out
}


S7::method(plot_data, list(PlotRandomCov, modelterms7::RandomTerm)) <-
  function(type, x, fit, ctx, ...) {
    check_plain_random(type, x, ctx)
    rb <- random_block(fit, ctx$param, ctx$name, x)
    d <- ncol(rb$b)
    if (d < 2L) {
      stop(sprintf(paste0("the term '%s' has one coefficient for each level; ",
                          "the covariance graph needs at least two."),
                   ctx$name), call. = FALSE)
    }
    w <- type@which
    if (any(w > d)) {
      stop(sprintf("'which' selects a position beyond the %d coefficients of a level.",
                   d), call. = FALSE)
    }
    pr <- random_prior_fit(fit, ctx$param, ctx$name, x)
    if (is.null(pr$Sigma)) {
      stop(sprintf(paste0("the covariance of the effects of the term '%s' is ",
                          "drawn for a gaussian prior only."), ctx$name),
           call. = FALSE)
    }
    S <- pr$Sigma[w, w]
    r <- sqrt(stats::qchisq(type@level, 2))
    th <- seq(0, 2 * pi, length.out = 181L)
    el <- r * cbind(cos(th), sin(th)) %*% chol(S)
    out <- list(points = data.frame(level = rb$levels, x = rb$b[, w[1L]],
                                    y = rb$b[, w[2L]]),
                ellipse = data.frame(x = el[, 1L], y = el[, 2L]))
    attr(out, "Sigma") <- S
    attr(out, "labels") <- list(
      title = sprintf("%s: %s", ctx$param, if (is.null(ctx$title)) ctx$name
                      else ctx$title),
      x = rb$wnames[w[1L]], y = rb$wnames[w[2L]],
      caption = sprintf("ellipse: %g%% of the fitted prior", 100 * type@level))
    out
  }


S7::method(plot_layers, PlotRandomCov) <- function(type, data, ...) {
  list(
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50"),
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50"),
    ggplot2::geom_path(data = data$ellipse,
                       mapping = ggplot2::aes(x = .data$x, y = .data$y),
                       colour = "#b2182b"),
    ggplot2::geom_point(data = data$points,
                        mapping = ggplot2::aes(x = .data$x, y = .data$y),
                        colour = "#1f5fa8", size = 2))
}


S7::method(plot_data, list(PlotCurves, modelterms7::RandomSmoothTerm)) <-
  function(type, x, fit, ctx, ...) {
    bp <- x@build_info
    gvars <- all.vars(bp$gexpr)
    if (!is.symbol(bp$gexpr)) {
      stop(sprintf(paste0("the grouping of the term '%s' is not a variable ",
                          "of the data."), ctx$name), call. = FALSE)
    }
    xvars <- setdiff(ctx$vars, gvars)
    data <- plot_needs_data(ctx, fit, type, vars = c(xvars, gvars))
    if (length(xvars) != 1L || !is.numeric(data[[xvars]])) {
      stop(sprintf(paste0("the graph %s draws a random smooth of one numeric ",
                          "variable."), type@label), call. = FALSE)
    }
    xv <- xvars
    gv <- gvars
    des <- statmod_design(fit@spec)
    idx <- des[[ctx$param]]$blocks[[ctx$name]]
    b <- fit@coefficients[[ctx$param]][idx]
    lev <- bp$glevels
    g <- as.character(data[[gv]])
    grid <- do.call(rbind, lapply(lev, function(l) {
      xs <- data[[xv]][g == l & !is.na(g)]
      if (!length(xs)) return(NULL)
      r <- range(xs, na.rm = TRUE)
      out <- data.frame(seq(r[1L], r[2L], length.out = type@n_grid))
      names(out) <- xv
      out[[gv]] <- factor(l, levels = lev)
      out
    }))
    f <- drop(as.matrix(modelterms7::term_predict(x, newdata = grid)) %*% b)
    pop <- data.frame(x = numeric(0), fit = numeric(0))
    if (type@population) {
      terms <- fit@spec@terms[[ctx$param]]
      cand <- Filter(function(k) {
        tk <- terms[[k]]
        S7::S7_inherits(tk, modelterms7::SmoothTerm) &&
          identical(term_variables(tk, k), xv)
      }, names(terms))
      if (length(cand) != 1L) {
        stop(sprintf(paste0("'population' needs one smooth of '%s' in the ",
                            "equation '%s', and there are %d."), xv, ctx$param,
                     length(cand)), call. = FALSE)
      }
      pgrid <- grid
      pf <- parts_contribution_values(fit, ctx$param, cand, terms[[cand]], pgrid)
      f <- f + pf
      xs <- seq(min(data[[xv]]), max(data[[xv]]), length.out = type@n_grid)
      pg <- data.frame(xs)
      names(pg) <- xv
      for (v in setdiff(names(grid), xv)) pg[[v]] <- grid[[v]][1L]
      pop <- data.frame(x = xs, fit = parts_contribution_values(
        fit, ctx$param, cand, terms[[cand]], pg))
    }
    out <- list(curves = data.frame(x = grid[[xv]], fit = f,
                                    group = grid[[gv]]),
                population = pop)
    attr(out, "labels") <- list(
      title = sprintf("%s: %s", ctx$param, if (is.null(ctx$title)) ctx$name
                      else ctx$title),
      x = xv, y = if (type@population) "curve of the level"
        else "departure from the population")
    out
  }


#' The Contribution of a Smooth at Fitted Coefficients
#'
#' @param fit A [StatmodFit()].
#' @param param,name,tm The equation, the name and the smooth term.
#' @param newdata A data frame.
#'
#' @return A numeric vector with one value for each row of `newdata`.
#'
#' @keywords internal
parts_contribution_values <- function(fit, param, name, tm, newdata) {
  des <- statmod_design(fit@spec)
  b <- fit@coefficients[[param]][des[[param]]$blocks[[name]]]
  drop(as.matrix(modelterms7::term_predict(tm, newdata = newdata)) %*% b)
}


S7::method(plot_layers, PlotCurves) <- function(type, data, ...) {
  out <- list(
    ggplot2::geom_line(data = data$curves,
                       mapping = ggplot2::aes(x = .data$x, y = .data$fit,
                                              group = .data$group),
                       colour = "grey55", alpha = 0.7, linewidth = 0.4))
  if (nrow(data$population)) {
    out <- c(out, list(ggplot2::geom_line(
      data = data$population, mapping = ggplot2::aes(x = .data$x, y = .data$fit),
      colour = "#1f5fa8", linewidth = 1.2)))
  }
  out
}
