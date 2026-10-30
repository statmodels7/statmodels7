#' @title S7 Classes of the Shrinkage and Path Graphs
#'
#' @description
#' The classes that [plot_shrinkage()] and [plot_path()] instantiate. Use the
#' constructors, which check the properties.
#'
#' @param label,sort,information,log_lambda See [plot_shrinkage()] and
#'   [plot_path()].
#'
#' @return A `PlotShrinkage` or `PlotPath` object, inheriting from
#'   [term_plot].
#'
#' @name PlotShrinkage
#' @aliases PlotPath
#' @keywords internal
PlotShrinkage <- S7::new_class("PlotShrinkage", parent = term_plot,
  properties = list(
    sort = S7::class_logical,
    information = S7::class_any
  ))

#' @rdname PlotShrinkage
#' @keywords internal
PlotPath <- S7::new_class("PlotPath", parent = term_plot,
  properties = list(
    log_lambda = S7::class_any,
    information = S7::class_any
  ))


#' Graph of the Shrinkage of a Penalized Block
#'
#' @description
#' Selects the graph that compares the coefficients of a penalized block with
#' the coefficients that the same model would give without the penalty of that
#' term. Each row is one coefficient: a filled point marks the penalized
#' estimate, an open circle marks the approximate unpenalized estimate, and
#' an arrow runs from the second to the first. The graph applies to the
#' terms [ridge()], [lasso()], [enet()], [scad()] and [mcp()], to a smooth
#' and to a random effect.
#'
#' @details
#' Let \eqn{b} be the penalized estimate of all the coefficients of the model,
#' \eqn{g(b)} the gradient of the log-likelihood (without any penalty) at
#' \eqn{b}, \eqn{I(b)} the information of the log-likelihood over all the
#' coefficients, and \eqn{S_o} the Hessian of the penalties of the other terms.
#' The unpenalized estimate of the coefficients of the term is one Newton step
#' on the model in which only the penalty of that term is removed,
#' \deqn{b_u = b + (I + S_o)^{-1} r,}
#' where \eqn{r} equals \eqn{g(b)} on the coefficients of the term and is zero
#' elsewhere. No model is refitted. Where the model has no other penalized
#' term the step is \eqn{b + I^{-1} g(b)}, which for a gaussian mean with the
#' identity link is the least squares estimate exactly. The coefficients of
#' the term are the only ones shown, but the step is taken on all the
#' coefficients, so it accounts for their correlation with the rest of the
#' model.
#'
#' For a kinked penalty (a lasso, a SCAD or an MCP) the coefficients that the
#' penalty holds at zero are part of the step: their unpenalized estimates
#' are those of the model without the penalty, and the arrow shows the
#' selection. The coefficients that a kinked penalty of another term holds at
#' zero stay at zero in the step. A coefficient that the fit left out because
#' its column is aliased is not drawn.
#'
#' The step needs the information \eqn{I + S_o} to be invertible. When it is
#' not (more coefficients than observations, or a random effect beside an
#' intercept, whose group means are not identified without the penalty) the
#' graph signals an error.
#'
#' @param sort Whether to order the coefficients by their penalized estimate.
#' @param information An object from [distributions7::expected_info()] (the
#'   default) or [distributions7::observed_info()], the information of the
#'   log-likelihood in the step.
#'
#' @param legend Whether to draw the legend that names the elements of
#'   the graph.
#'
#' @return A `PlotShrinkage` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_path()], [plot_caterpillar()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
#' dd$y <- 1 + dd$x1 + 0.5 * dd$x2 + rnorm(100)
#' fit <- statmod(y ~ ridge(~ x1 + x2 + x3), distributions7::gaussian1_distrib(),
#'                dd)
#' ggplot2::autoplot(fit, term = "ridge(~x1 + x2 + x3)",
#'                   plot_type = plot_shrinkage())
#' @export
plot_shrinkage <- function(sort = TRUE,
                           information = distributions7::expected_info(),
                           legend = TRUE) {
  if (!is.logical(sort) || length(sort) != 1L || is.na(sort)) {
    stop("'sort' must be TRUE or FALSE.", call. = FALSE)
  }
  check_plot_information(information)
  PlotShrinkage(label = "plot_shrinkage()", sort = sort,
                information = information, legend = check_legend(legend))
}


#' Graph of the Coefficient Paths of a Penalized Block
#'
#' @description
#' Selects the graph of the coefficients of a penalized block against the
#' logarithm of the penalty parameter \eqn{\lambda}, with the fitted value of
#' \eqn{\lambda} marked by a dashed line. For a kinked penalty (a lasso, an
#' elastic net, a SCAD or an MCP) a second axis counts the nonzero
#' coefficients.
#'
#' @details
#' For a kinked penalty the path is computed by refitting the model at each
#' value of \eqn{\lambda}, from the largest to the smallest, each fit starting
#' from the previous one. The other hyperparameters of the penalty and all the
#' coefficients of the other equations are held at their fitted values, so
#' that the path at the fitted \eqn{\lambda} is the fit. By default the grid
#' runs from the smallest \eqn{\lambda} that empties the block, found by
#' doubling the fitted value, down to a tenth of the fitted value.
#'
#' For a penalty that is quadratic in the coefficients, \eqn{\lambda S} with
#' \eqn{S} fixed, the path is the closed form
#' \deqn{b(\lambda) = (I + S_o + \lambda S)^{-1}
#'   (I + S_o + \hat\lambda S)\,\hat b,}
#' where \eqn{\hat b} and \eqn{\hat\lambda} are the fitted coefficients and
#' penalty parameter, \eqn{I} is the information of the log-likelihood over
#' all the coefficients and \eqn{S_o} the Hessian of the penalties of the other
#' terms. It is the exact mode for a gaussian mean with the identity link, with
#' the other equations held, and the first-order approximation obtained by
#' replacing the log-likelihood with its quadratic expansion at \eqn{\hat b}
#' otherwise. By default \eqn{\log\lambda} runs over five units on each side
#' of the fitted value.
#'
#' @param log_lambda `NULL`, or an increasing numeric vector of values of
#'   \eqn{\log\lambda}.
#' @param information An object from [distributions7::expected_info()] (the
#'   default) or [distributions7::observed_info()], the information of the
#'   log-likelihood in the closed form of a quadratic penalty.
#'
#' @param legend Whether to draw the legend that names the elements of
#'   the graph.
#'
#' @return A `PlotPath` object, inheriting from [term_plot].
#'
#' @seealso [autoplot.StatmodFit()], [plot_shrinkage()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
#' dd$y <- 1 + dd$x1 + 0.5 * dd$x2 + rnorm(100)
#' fit <- statmod(y ~ lasso(~ x1 + x2 + x3), distributions7::gaussian1_distrib(),
#'                dd)
#' ggplot2::autoplot(fit, term = "lasso(~x1 + x2 + x3)", plot_type = plot_path())
#' @export
plot_path <- function(log_lambda = NULL,
                      information = distributions7::expected_info(),
                      legend = TRUE) {
  if (!is.null(log_lambda) &&
      (!is.numeric(log_lambda) || length(log_lambda) < 2L ||
       any(!is.finite(log_lambda)) || is.unsorted(log_lambda, strictly = TRUE))) {
    stop("'log_lambda' must be NULL or an increasing numeric vector of at ",
         "least two finite values.", call. = FALSE)
  }
  check_plot_information(information)
  PlotPath(label = "plot_path()", log_lambda = log_lambda,
           information = information, legend = check_legend(legend))
}


#' Check the Information Argument of a Graph
#'
#' @param information The value to check.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
check_plot_information <- function(information) {
  if (!S7::S7_inherits(information, distributions7::ExpectedInfo) &&
      !S7::S7_inherits(information, distributions7::ObservedInfo)) {
    stop("'information' must be expected_info() or observed_info().",
         call. = FALSE)
  }
  invisible(NULL)
}


S7::method(default_plot, modelterms7::PenalizedTerm) <- function(term, ...) {
  plot_shrinkage()
}


#' Does a Penalized Unit Belong to a Term
#'
#' @param u A unit from `statmod_penalized()`.
#' @param param The equation of the term.
#' @param name The name of the term.
#'
#' @return A single logical.
#'
#' @keywords internal
unit_of_term <- function(u, param, name) {
  if (!is.null(u$pieces)) {
    return(any(vapply(u$pieces, function(pc)
      identical(pc$param, param) && identical(pc$term, name), logical(1))))
  }
  identical(u$param, param) && identical(u$term, name)
}


#' The Penalized System Around One Term
#'
#' @description
#' Collects what the graphs of a penalized block need: the stacked
#' coefficients, the positions of the term's coefficients, the information of
#' the log-likelihood over all the coefficients, the Hessian of the penalties
#' of the other terms, and the coordinates that are held out of a Newton step.
#'
#' @param fit A [StatmodFit()].
#' @param param,name,cols The equation, the term and the columns of the term
#'   that are drawn (`NULL` for all).
#' @param information An object from [distributions7::expected_info()] or
#'   [distributions7::observed_info()].
#'
#' @return A list with `beta`, `idx` (positions of the drawn coefficients),
#'   `names` (their names), `grad` (the gradient of the log-likelihood), `info`
#'   (the information), `S_other`, `units` (the penalized units of the term),
#'   `free` (the coordinates of a Newton step), `spec`, `design` and `labels`.
#'
#' @keywords internal
penalty_system <- function(fit, param, name, cols, information) {
  spec <- fit@spec
  design <- statmod_design(spec)
  if (length(attr(design, "structural"))) {
    stop("this graph does not apply to a model with a structural term.",
         call. = FALSE)
  }
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  beta <- unlist(fit@coefficients[params], use.names = FALSE)
  blk <- design[[param]]$blocks[[name]]
  if (is.null(blk)) {
    stop(sprintf("the term '%s' has no coefficients in '%s'.", name, param),
         call. = FALSE)
  }
  if (!is.null(cols)) blk <- blk[cols]
  idx <- offs[match(param, params)] + blk
  labs <- rownames(coef_labels(spec, design))
  units <- statmod_penalized(spec, design)
  mine <- vapply(units, unit_of_term, logical(1), param = param, name = name)

  I <- -as.matrix(hessian(fit, information = information))
  g <- unlist(gradient(fit), use.names = FALSE)

  So <- matrix(0, total, total)
  for (u in units[!mine]) {
    if (isTRUE(u$structural) || isTRUE(u$mixed) || is.null(u$index)) next
    th <- as.list(fit@hyper[[u$param]][[u$key]])
    blkh <- as.matrix(penalties7::penalty_hessian(u$penalty, beta[u$index], th))
    blkh[!is.finite(blkh)] <- 0
    So[u$index, u$index] <- So[u$index, u$index] + blkh
  }

  out_of_step <- which(labs %in% fit@aliased)
  kinked <- setdiff(kinked_coords(spec, design), idx)
  out_of_step <- union(out_of_step, kinked[beta[kinked] == 0])
  free <- setdiff(seq_len(total), out_of_step)
  list(beta = beta, idx = idx, names = design[[param]]$coef_names[blk],
       grad = g, info = I, S_other = So, units = units[mine], free = free,
       spec = spec, design = design, labels = labs, mine = mine)
}


#' A Newton Step Without the Penalty of a Term
#'
#' @param sys The list from [penalty_system()].
#'
#' @return The unpenalized estimates of the drawn coefficients, `NA` where the
#'   coefficient is out of the step.
#'
#' @keywords internal
unpenalized_estimates <- function(sys) {
  free <- sys$free
  M <- sys$info[free, free, drop = FALSE] + sys$S_other[free, free, drop = FALSE]
  if (!all(is.finite(M)) || !is.finite(rc <- rcond(M)) || rc < 1e-12) {
    stop(paste0("the information of the log-likelihood cannot be inverted ",
                "(more coefficients than the data identify, or columns that ",
                "are aliased once the penalty of the term is removed), so ",
                "the unpenalized estimates are not defined."), call. = FALSE)
  }
  r <- numeric(length(free))
  at <- match(sys$idx, free)
  r[at[!is.na(at)]] <- sys$grad[sys$idx[!is.na(at)]]
  step <- drop(solve(M, r))
  out <- rep(NA_real_, length(sys$idx))
  out[!is.na(at)] <- sys$beta[sys$idx[!is.na(at)]] + step[at[!is.na(at)]]
  out
}


S7::method(plot_data, list(PlotShrinkage, modelterms7::PenalizedTerm)) <-
  function(type, x, fit, ctx, ...) shrinkage_data(type, x, fit, ctx)

S7::method(plot_data, list(PlotShrinkage, modelterms7::SmoothTerm)) <-
  function(type, x, fit, ctx, ...) shrinkage_data(type, x, fit, ctx)

S7::method(plot_data, list(PlotShrinkage, modelterms7::RandomTerm)) <-
  function(type, x, fit, ctx, ...) {
    if (S7::S7_inherits(x, modelterms7::RandomSmoothTerm)) {
      stop(sprintf(paste0("the graph %s is not available for the term '%s'; ",
                          "plot_curves() draws the curves of a random smooth."),
                   type@label, ctx$name), call. = FALSE)
    }
    shrinkage_data(type, x, fit, ctx)
  }


#' The Data of a Shrinkage Graph
#'
#' @param type A `PlotShrinkage` object.
#' @param x The term.
#' @param fit A [StatmodFit()].
#' @param ctx The context described at [plot_data()].
#'
#' @return A list with the data frame `estimates` and the attribute
#'   `"labels"`.
#'
#' @keywords internal
shrinkage_data <- function(type, x, fit, ctx) {
  sys <- penalty_system(fit, ctx$param, ctx$name, ctx$cols, type@information)
  bu <- unpenalized_estimates(sys)
  bp <- sys$beta[sys$idx]
  keep <- !is.na(bu)
  est <- data.frame(coef = sys$names[keep], penalized = bp[keep],
                    unpenalized = bu[keep], stringsAsFactors = FALSE)
  est$shrinkage <- est$penalized - est$unpenalized
  if (type@sort) est <- est[order(est$penalized), , drop = FALSE]
  est$coef <- factor(est$coef, levels = unique(est$coef))
  rownames(est) <- NULL
  out <- list(estimates = est)
  attr(out, "labels") <- list(
    title = sprintf("%s: %s", ctx$param, if (is.null(ctx$title)) ctx$name
                    else ctx$title),
    x = "estimate", y = NULL)
  out
}


S7::method(plot_layers, PlotShrinkage) <- function(type, data, ...) {
  d <- data$estimates
  k_pen <- "penalized estimate"
  k_unp <- "non-penalized estimate"
  k_arr <- "shrinkage"
  list(
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50"),
    ggplot2::geom_segment(
      data = d, mapping = ggplot2::aes(x = .data$unpenalized, xend = .data$penalized,
                                       y = .data$coef, yend = .data$coef,
                                       colour = k_arr),
      arrow = grid::arrow(length = grid::unit(0.08, "inches"), type = "closed")),
    ggplot2::geom_point(data = d, mapping = ggplot2::aes(x = .data$unpenalized,
                                                         y = .data$coef,
                                                         colour = k_unp),
                        shape = 1, size = 2.4),
    ggplot2::geom_point(data = d, mapping = ggplot2::aes(x = .data$penalized,
                                                         y = .data$coef,
                                                         colour = k_pen),
                        shape = 16, size = 2.4),
    plot_key(colour = stats::setNames(c("#1f5fa8", "grey30", "grey40"),
                                      c(k_pen, k_unp, k_arr)))
  )
}


S7::method(plot_data, list(PlotPath, modelterms7::PenalizedTerm)) <-
  function(type, x, fit, ctx, ...) {
    sys <- penalty_system(fit, ctx$param, ctx$name, ctx$cols, type@information)
    unit <- sys$units
    unit <- Filter(function(u) "lambda" %in% names(fit@hyper[[u$param]][[u$key]]),
                   unit)
    if (length(unit) != 1L) {
      stop(sprintf(paste0("the term '%s' has no single penalty parameter ",
                          "'lambda', which the path varies."), ctx$name),
           call. = FALSE)
    }
    u <- unit[[1L]]
    th <- fit@hyper[[u$param]][[u$key]]
    lam_hat <- unname(th[["lambda"]])
    if (!is.finite(lam_hat) || lam_hat <= 0) {
      stop("the fitted penalty parameter is not a positive number.",
           call. = FALSE)
    }
    kink <- isTRUE(penalty_has_kink(u$penalty))
    path <- if (kink) {
      kinked_path(fit, sys, u, lam_hat, type@log_lambda)
    } else {
      quadratic_path(fit, sys, u, lam_hat, type@log_lambda)
    }
    paths <- data.frame(
      log_lambda = rep(path$log_lambda, times = length(sys$idx)),
      coef = factor(rep(sys$names, each = length(path$log_lambda)),
                    levels = sys$names),
      estimate = as.vector(path$coef))
    out <- list(paths = paths,
                fitted = data.frame(log_lambda = log(lam_hat)),
                count = data.frame(log_lambda = numeric(0), y = numeric(0)))
    if (kink) {
      cnt <- rowSums(path$coef != 0, na.rm = TRUE)
      rng <- range(paths$estimate, na.rm = TRUE)
      p <- length(sys$idx)
      out$count <- data.frame(log_lambda = path$log_lambda,
                              n = cnt,
                              y = rng[1L] + cnt / p * diff(rng))
      attr(out, "count_scale") <- list(ymin = rng[1L], span = diff(rng), p = p)
    }
    attr(out, "labels") <- list(
      title = sprintf("%s: %s", ctx$param, if (is.null(ctx$title)) ctx$name
                      else ctx$title),
      x = "log lambda", y = "coefficient", colour = "coefficient")
    out
  }


#' The Path of a Kinked Penalty by Refitting
#'
#' @description
#' Fits the model at each value of \eqn{\lambda} with the other
#' hyperparameters and the coefficients of the other equations held at their
#' fitted values, from the largest value to the smallest and each fit from the
#' previous one.
#'
#' @param fit A [StatmodFit()].
#' @param sys The list from [penalty_system()].
#' @param u The penalized unit of the term.
#' @param lam_hat The fitted penalty parameter.
#' @param log_lambda `NULL` or the grid of \eqn{\log\lambda}.
#'
#' @return A list with `log_lambda` (increasing) and `coef`, a matrix with one
#'   row per value and one column per drawn coefficient.
#'
#' @keywords internal
kinked_path <- function(fit, sys, u, lam_hat, log_lambda) {
  spec0 <- fit@spec
  params <- spec0@distrib@params
  design0 <- sys$design
  # the coefficients of the other equations are held, so that the path at the
  # fitted value is the fit (a dispersion estimated on the marginal criterion
  # is not the joint mode)
  others <- setdiff(params, u$param)
  hc <- spec0@held_coef
  for (p in others) {
    if (design0[[p]]$npar == 0L) next
    hc[[p]] <- stats::setNames(as.numeric(fit@coefficients[[p]]),
                               design0[[p]]$coef_names)
  }
  spec <- S7::set_props(spec0, held_coef = hc)
  design <- statmod_design(spec)
  blocks <- statmod_blocks(spec, design)
  method <- fit@methods$smooth
  cfg <- inner_settings(method)
  beta0 <- sys$beta
  hy <- fit@hyper
  fit_at <- function(lam, beta) {
    h <- hy
    v <- h[[u$param]][[u$key]]
    v[["lambda"]] <- lam
    h[[u$param]][[u$key]] <- v
    statmod_alternate(spec, design, blocks, h, method, beta, cfg$expected,
                      cfg$emethod, cfg$maxit, cfg$tol, verbosity(0))
  }
  cols <- sys$idx
  if (is.null(log_lambda)) {
    top <- lam_hat
    beta <- beta0
    for (k in seq_len(30L)) {
      r <- fit_at(top, beta)
      beta <- r$par
      if (all(beta[cols] == 0)) break
      top <- top * 2
    }
    log_lambda <- seq(log(lam_hat) - log(10), log(top), length.out = 50L)
    if (log(top) > max(log_lambda) - 1e-12) log_lambda[length(log_lambda)] <- log(top)
    log_lambda <- sort(c(log_lambda[abs(log_lambda - log(lam_hat)) > 1e-8],
                         log(lam_hat)))
  }
  lam <- exp(log_lambda)
  ord <- order(lam, decreasing = TRUE)
  res <- matrix(NA_real_, length(lam), length(cols))
  beta <- beta0
  beta[cols] <- 0
  for (i in ord) {
    r <- fit_at(lam[i], beta)
    beta <- r$par
    res[i, ] <- beta[cols]
  }
  list(log_lambda = log_lambda, coef = res)
}


#' The Path of a Quadratic Penalty in Closed Form
#'
#' @param fit A [StatmodFit()].
#' @param sys The list from [penalty_system()].
#' @param u The penalized unit of the term.
#' @param lam_hat The fitted penalty parameter.
#' @param log_lambda `NULL` or the grid of \eqn{\log\lambda}.
#'
#' @return A list with `log_lambda` and `coef`, as [kinked_path()].
#'
#' @keywords internal
quadratic_path <- function(fit, sys, u, lam_hat, log_lambda) {
  if (is.null(log_lambda)) {
    log_lambda <- seq(log(lam_hat) - 5, log(lam_hat) + 5, length.out = 61L)
    log_lambda[31L] <- log(lam_hat)
  }
  th <- fit@hyper[[u$param]][[u$key]]
  total <- length(sys$beta)
  S_at <- function(lam) {
    h <- as.list(th)
    h[["lambda"]] <- lam
    S <- matrix(0, total, total)
    blkh <- as.matrix(penalties7::penalty_hessian(u$penalty, sys$beta[u$index],
                                                  h))
    S[u$index, u$index] <- blkh
    S
  }
  free <- sys$free
  base <- sys$info[free, free, drop = FALSE] + sys$S_other[free, free, drop = FALSE]
  b <- sys$beta[free]
  rhs <- drop((base + S_at(lam_hat)[free, free, drop = FALSE]) %*% b)
  if (!all(is.finite(rhs)) || !is.finite(rcond(base))) {
    stop("the information of the log-likelihood cannot be inverted, so the ",
         "path is not defined.", call. = FALSE)
  }
  at <- match(sys$idx, free)
  res <- matrix(NA_real_, length(log_lambda), length(sys$idx))
  for (i in seq_along(log_lambda)) {
    M <- base + S_at(exp(log_lambda[i]))[free, free, drop = FALSE]
    res[i, !is.na(at)] <- solve(M, rhs)[at[!is.na(at)]]
  }
  list(log_lambda = log_lambda, coef = res)
}


S7::method(plot_layers, PlotPath) <- function(type, data, ...) {
  k_fit <- "fitted lambda"
  k_cnt <- "number of nonzero coefficients"
  out <- list(
    ggplot2::geom_hline(yintercept = 0, colour = "grey70"),
    ggplot2::geom_line(data = data$paths,
                       mapping = ggplot2::aes(x = .data$log_lambda,
                                              y = .data$estimate,
                                              colour = .data$coef),
                       linewidth = 0.7),
    ggplot2::geom_vline(data = data$fitted,
                        mapping = ggplot2::aes(xintercept = .data$log_lambda,
                                               linetype = k_fit),
                        colour = "#b2182b"),
    # the colour names the coefficients, so the other lines are named by
    # their line type
    ggplot2::scale_linetype_manual(
      name = NULL, values = stats::setNames(c("dashed", "dotted"),
                                            c(k_fit, k_cnt)),
      breaks = c(k_fit, k_cnt)),
    # the coefficients are named in the legend only while it stays readable
    ggplot2::guides(colour = if (nlevels(factor(data$paths$coef)) <= 12L)
      ggplot2::guide_legend(order = 1) else "none")
  )
  if (nrow(data$count)) {
    sc <- attr(data, "count_scale")
    out <- c(out, list(
      ggplot2::geom_step(data = data$count,
                         mapping = ggplot2::aes(x = .data$log_lambda,
                                                y = .data$y,
                                                linetype = k_cnt),
                         colour = "grey40"),
      ggplot2::scale_y_continuous(
        sec.axis = ggplot2::sec_axis(
          function(y) (y - sc$ymin) / sc$span * sc$p,
          name = "nonzero coefficients"))))
  }
  out
}
