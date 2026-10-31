#' @include spec.R
NULL

#' Simulate a Response From a Written Model
#'
#' @description
#' Takes a formula and a distribution, draws coefficients or uses the ones
#' given, and returns the data drawn from that model together with the truth
#' behind it.
#'
#' @details
#' # What it is for
#'
#' Simulated data have a known truth: write the model, draw from it, fit it
#' back, and see whether the fit recovers what was put in. That is the shape
#' of a simulation study, and of a check on a term that has just been written.
#'
#' This is not [stats::simulate()], which draws from a model already fitted.
#' The `r` prefix is the one that R uses for a random draw, so the two names
#' stay apart.
#'
#' A covariate needs no declaration. A factor becomes its contrasts and a
#' numeric variable stays as it is, and the design is built by the interpreter
#' that a fit uses.
#'
#' # The returned truth
#'
#' A simulation study compares against the coefficients, the parameters that
#' they gave and whatever a term with state drew, so the result is a list
#' holding all of them. Pass `sim$data` where a data frame is wanted.
#'
#' The truth is kept as elements of a list and not as attributes of the data
#' frame, because an attribute is not subset with the rows: `sim[1:10, ]`
#' would keep a `theta` of the original length, while `subset()` and
#' `merge()` drop it.
#'
#' # Assembly of the predictor
#'
#' The predictor is built through [statmod_design_at()], so the simulated data
#' come from the model that was written and not from a linearization of it.
#'
#' A term whose block moves with its coefficients, [modelterms7::seg()],
#' [modelterms7::jseg()] or [modelterms7::nl()], contributes
#' `term_value()` at the coefficients supplied, not its block times them. The
#' two differ by the whole nonlinearity, and the simulation uses the first.
#'
#' # Data, or a row count
#'
#' `data` carries the covariates. A model with none, a pure time series say,
#' needs only `n`. One of the two must be given, and where both are given
#' they must agree.
#'
#' # Fixed covariates or drawn ones
#'
#' `covariates` takes one function of the observation count per column, as an
#' entry of `par` may be, and the columns are drawn afresh at every replicate.
#' A column of `data` under the same name is overwritten, and the columns of
#' `data` that are not drawn are the same in every replicate.
#'
#' The choice determines what a study measures. With `data` the covariates are
#' the same throughout, so the study measures the behavior of the estimator
#' conditional on that design; with `covariates` it measures that behavior
#' over the design as well. The two choices are equally correct, and a study
#' should state which one it ran.
#'
#' A drawn factor should have its levels fixed. A factor drawn with free
#' levels can lose a level at some replicate. The design then has another
#' shape than the one for which the coefficients were drawn, and the
#' simulation stops at that replicate with an error.
#'
#' # Several replicates
#'
#' `n_sim` draws that many data sets. The truth is drawn once and shared by
#' the replicates: a study over replicates measures the variability of an
#' estimator at a given set of parameters, so the replicates differ in what is
#' random and not in what is being estimated. Varying the truth as well
#' requires a loop over calls and reads differently.
#'
#' With `n_sim > 1` the per-replicate fields `data` and `theta` come back as
#' lists of that length, and `latent` does so when the model has a term with
#' state (otherwise it is `NULL`), while `par`, `structural` and `hyper` stay
#' single.
#'
#' # The quantities that are drawn
#'
#' If the user names no quantity, the call draws the whole truth, so writing
#' the model is enough to get data from it. Each quantity is drawn by the
#' component that defines its meaning:
#'
#' - a coefficient of a design column has no other owner and comes from
#'   `rnorm(1, 0, sd)`, which on the link scale gives predictors of order
#'   one;
#' - a coordinate that some penalty covers is drawn from that penalty read as
#'   a prior, through [penalties7::penalty_draw()]. A Gaussian random effect
#'   gives Gaussian effects, a lasso Laplace ones and a heavy-tailed prior
#'   heavy-tailed ones. The prior's own scale is drawn as well, on the chart
#'   that the penalty carries for it, and comes back in `hyper`;
#' - a structural term's own parameters are drawn by the term, through
#'   [modelterms7::term_draw()], on the chart that each one rides. A loading
#'   stays positive and a persistence stationary whatever comes out;
#' - a coefficient of a design column whose meaning only the term defines is
#'   drawn by the term, through [modelterms7::term_coef_draw()], last of all,
#'   so that the values that it writes survive. A break-point is such a
#'   coordinate: it is a position on the covariate's own axis, so a normal of
#'   width `sd` often lands outside the data, and the term then confines it,
#'   which pins it to an edge of the interval, where one of the two segments
#'   holds a small share of the rows.
#'
#' A hyperparameter that the term holds is used and not drawn, so
#' `s(x, bspline_smooth(), hyper = c(lambda = 2))` simulates at the smoothing
#' parameter that it names. If a term draws the coordinates that a prior
#' covers, the prior is reported at the width that the term used, since those
#' values are no longer the prior's own draw. They are drawn Gaussian at that
#' width, so a heavy-tailed prior over a break-point keeps its family for the
#' fitting but not for the simulation.
#'
#' A coordinate that a prior does not reach falls back to the plain draw. The
#' cases are SCAD and MCP, which are improper by construction; an anisotropic
#' tensor smooth, which is flat along its null space; a covariance class
#' spanning a filter and an equation at once, which is in neither vector on
#' its own; and a coordinate for which the draw of the penalty is missing.
#'
#' # Fixing quantities by name
#'
#' `par` is one named list over the whole model, and a key may be
#'
#' - a distribution parameter, `mu`, which is that whole equation;
#' - one of its coefficients or a group of them, `mu.(Intercept)` or
#'   `mu.random`, where a group is a name that the members extend at a dot;
#' - a structural term's own parameter, `kappa1`, or a group of those,
#'   `omega.random`.
#'
#' A value is a vector of that key's own length, a single number used for all
#' of them, or a function of the count. A function writes a structured
#' truth without a vocabulary of its own:
#' `function(k) rnorm(k, 0, 0.4)` is a random effect at a scale chosen by the
#' caller and `function(k) c(1.5, -2, rep(0, k - 2))` is a sparse truth for a
#' lasso to find. A function that returns the wrong number of values or a
#' missing value is rejected, because R would recycle a vector of the wrong
#' length into a different model. A key that addresses nothing is rejected
#' with a message that lists what the model does carry.
#'
#' A structural parameter is named on the scale of
#' [modelterms7::term_params()]: a loading is the loading and not its
#' logarithm, and a persistence is the partial autocorrelation that its chart
#' carries. A parameter that a subformula develops is different. Its
#' coordinates are the coefficients of that development, which act on the
#' unconstrained scale of the parameter's own chart, so `kappa1` is a loading
#' and `kappa1.random.3` is a group's departure on the log scale that the
#' loading rides. This keeps every group's loading positive whatever the
#' departure is.
#'
#' # A term with state
#'
#' The simulation runs through [modelterms7::term_simulate()], so the
#' recursion that generates the data is the term's own. A score-driven term
#' draws the response as it runs, and its level at one time is driven by the
#' score at the time before; a latent chain draws its path from the stationary
#' law with which the likelihood is written; a marginal break-point term draws
#' each group's positions from their prior. The values that each term drew
#' come back in `latent`, and a recovery check compares a fit with them.
#'
#' The parameters of such a term are not coefficients of any equation. They
#' are named in `par` beside the equations' coefficients, and since a formula
#' holds at most one such term, the names carry no term key. They are drawn
#' like every other coefficient.
#'
#' # The response's name
#'
#' The left-hand side of the formula must be a symbol. `log(y) ~ x` is
#' rejected because the model generates values of `log(y)`, which cannot be
#' stored in a column under either name. A censored response is rejected too,
#' for the reason that [statmod()] rejects one: there is no censored
#' likelihood.
#'
#' @param formula The model formula, as [statmod()] takes it, with the
#'   equations separated by `|`. Its left-hand side must be a symbol.
#' @param distrib A \pkg{distributions7} distribution object, which
#'   determines how many equations there are and the distribution from which
#'   the data are drawn.
#' @param data A data frame of covariates, or `NULL` where `n` is given.
#' @param n The number of observations, a positive whole number. It is
#'   required where `data` is `NULL`, and where `data` is given it must equal
#'   the number of rows.
#' @param n_sim How many data sets to draw, a positive whole number. The
#'   default is `1`.
#' @param par Optional named list of what to hold, over the whole model: a
#'   distribution parameter, one of its coefficients or a group of them, a
#'   structural term's own parameter or a group of those. Each entry is a
#'   numeric vector, a single number or a function of the count. See the
#'   details. `NULL`, the default, draws everything.
#' @param sd The width of the draws, `1` by default. A quantity that rides a
#'   chart (a structural parameter, a prior's own scale) is drawn at half of
#'   it, because a chart carries a width of one onto most of its parameter's
#'   range.
#' @param offsets Optional named list of offsets, one per parameter, as
#'   [statmod()] takes them.
#' @param covariates Optional named list of functions of the observation
#'   count, one per covariate, drawn afresh at every replicate. A drawn
#'   factor should have its levels fixed, since a replicate that loses a
#'   level stops the simulation with an error.
#'
#' @return An object of class `"StatmodSim"`, a list with the elements:
#'   \describe{
#'     \item{`data`}{the data frame, to which the response column has been
#'       added under the name of the formula's left-hand side.}
#'     \item{`par`}{the coefficients used, drawn or given, named as the
#'       design names them.}
#'     \item{`theta`}{the distribution's parameters at every observation, a
#'       named list with one vector per parameter.}
#'     \item{`latent`}{what a term with state drew (the filtered level of a
#'       score-driven term, the states of a latent chain, or the positions of
#'       a marginal break-point term), or `NULL`.}
#'     \item{`structural`}{such a term's own parameters, or `NULL`.}
#'     \item{`hyper`}{a data frame of the hyperparameters, with one row per
#'       penalty and name; each row gives the value drawn or held. Its columns
#'       are those of [hyper()] as far as they mean the same thing, and `value`
#'       takes the place of `estimate`, so a study compares the two
#'       directly.}
#'     \item{`n_sim`}{the number of data sets, as an integer.}
#'     \item{`call`}{the matched call.}
#'   }
#'   With `n_sim > 1` the fields `data` and `theta` are lists of that length,
#'   and `latent` is such a list when the model has a term with state and
#'   `NULL` otherwise; `par`, `structural` and `hyper` stay single.
#'
#' @seealso [statmod()] to fit what this draws,
#'   [simulate.StatmodFit()] to draw from a model already fitted,
#'   [print.StatmodSim()] for the printed form.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(50), g = factor(rep(c("a", "b"), 25)))
#'
#' # The coefficients are drawn.
#' sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd)
#' sim$par
#'
#' # The coefficients are given, and a fit recovers them.
#' sim2 <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
#'                  par = list(mu = c(1, 2), sigma = log(0.3)))
#' coef(statmod(y ~ x, distributions7::gaussian1_distrib(), sim2$data))
#'
#' # A sparse truth is written as a function of the coefficient count.
#' dd2 <- as.data.frame(matrix(rnorm(50 * 6), 50, 6))
#' sim3 <- rstatmod(y ~ lasso(~ V1 + V2 + V3 + V4 + V5 + V6),
#'                  distributions7::gaussian1_distrib(), dd2,
#'                  par = list(mu = function(k) c(2, -1.5, rep(0, k - 2)),
#'                             sigma = log(0.3)))
#' head(sim3$data$y, 3)
#'
#' # This model has no covariates at all.
#' sim4 <- rstatmod(y ~ 1, distributions7::gaussian1_distrib(), n = 20)
#' nrow(sim4$data)
#'
#' # A study over replicates draws the covariates afresh at each one.
#' study <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), n = 80,
#'                   n_sim = 5, par = list(mu = c(1, 2), sigma = log(0.5)),
#'                   covariates = list(x = function(n) runif(n, -2, 2)))
#' length(study$data)
#' vapply(study$data, function(d) coef(statmod(
#'   y ~ x, distributions7::gaussian1_distrib(), d),
#'   readable = FALSE)$mu[[2L]], numeric(1))
#'
#' # In a score-driven series, the term's own parameters are named beside the
#' # equations'.
#' sim5 <- rstatmod(y ~ 0 + gas(p = 1, q = 1, time = t),
#'                  distributions7::gaussian1_distrib(),
#'                  data.frame(t = 1:100),
#'                  par = list(sigma = 0, omega = 0.4, kappa1 = 0.3,
#'                             pacf1 = 0.6))
#' head(sim5$latent, 3)
#'
#' # In a panel whose level varies by group, naming nothing draws the
#' # deviations from the prior that random() declares, at a scale drawn too.
#' pan <- data.frame(id = factor(rep(1:6, each = 10)), t = rep(1:10, 6))
#' sim6 <- rstatmod(y ~ 0 + gas(p = 1, q = 1, omega ~ 1 + random(~1 | id),
#'                              by = id, time = t),
#'                  distributions7::gaussian1_distrib(), pan)
#' sim6$hyper
#' length(unique(round(sim6$theta$mu, 6))) > 6
#'
#' @export
rstatmod <- function(formula, distrib, data = NULL, n = NULL, n_sim = 1,
                     par = NULL, sd = 1, offsets = NULL,
                     covariates = NULL) {
  if (!S7::S7_inherits(distrib, distributions7::distrib)) {
    stop("'distrib' must be a distributions7 distribution object.",
         call. = FALSE)
  }
  if (!is_count(n_sim)) {
    stop("'n_sim' must be a positive whole number.", call. = FALSE)
  }
  n_sim <- as.integer(n_sim)
  covariates <- check_covariates(covariates)
  # the supplied columns, before any draw: every replicate after the first
  # draws its covariates into this frame, so the columns that are not drawn
  # are the same in every replicate
  data0 <- rstatmod_data(data, n, NULL)
  data <- rstatmod_data(data, n, covariates)
  n <- nrow(data)
  params <- distrib@params
  split <- statmod_equations(formula, params)
  nm <- rstatmod_response_name(split$response)

  build <- function(dat) {
    # THE SPECIFICATION IS THE REAL ONE, built against a placeholder
    # response, so that the design, the offsets, the interpreter and every
    # term's `build_info` are the ones a fit would produce. What the
    # placeholder is used for is narrow: a break-point term whose starting
    # position the caller did not name chooses it on a least-squares profile
    # of the response, and that choice is a starting value the coefficients
    # then overwrite.
    sim <- dat
    sim[[nm]] <- stats::rnorm(nrow(dat))
    sp <- statmod_spec(formula, distrib, sim, offsets = offsets)
    list(spec = sp, design = statmod_design(sp))
  }

  b <- build(data)
  # THE TRUTH IS DRAWN ONCE. What a study over replicates measures is the
  # variability of an estimator AT a set of parameters, so the replicates
  # differ in what is random -- the response, and the covariates where they
  # are simulated -- and not in what is being estimated. Varying the truth
  # as well is a loop over calls, and reads differently.
  truth <- rstatmod_truth(b$spec, b$design, par, sd)
  coef <- truth$coef
  names_at <- function(d) lapply(d[params], `[[`, "coef_names")
  cn <- names_at(b$design)
  # the structural term's parameters are drawn too, so a replicate whose
  # covariates changed the term's own shape would recycle another model's
  # truth exactly as a changed design would
  sn <- names(truth$psi)

  dats <- vector("list", n_sim)
  thetas <- vector("list", n_sim)
  latents <- vector("list", n_sim)
  psi <- NULL
  for (r in seq_len(n_sim)) {
    if (r > 1L && !is.null(covariates)) {
      data <- draw_covariates(covariates, n, data0)
      b <- build(data)
      # a design that changed shape cannot take coefficients drawn against
      # the first one, and recycling them would fit a different model
      su2 <- attr(b$design, "structural")
      sn2 <- if (length(su2))
        modelterms7::term_params(
          b$spec@terms[[su2[[1L]]$param]][[su2[[1L]]$term]]) else NULL
      if (!identical(names_at(b$design), cn) || !identical(sn2, sn)) {
        stop("the simulated covariates gave a design of another shape at ",
             "replicate ", r, ".\n  A factor that lost a level is the ",
             "ordinary cause; draw the covariates so that\n  every ",
             "replicate spans the same columns.", call. = FALSE)
      }
    }
    ep <- rstatmod_eta(b$spec, b$design, coef, truth$psi)
    y <- ep$y
    if (is.null(y)) y <- distributions7::distrib_rng(distrib, n, ep$theta)
    out <- data
    out[[nm]] <- y
    dats[[r]] <- out
    thetas[[r]] <- lapply(ep$theta, rep_len, n)
    latents[r] <- list(ep$latent)
    psi <- ep$structural
  }

  one <- identical(n_sim, 1L)
  structure(list(data = if (one) dats[[1L]] else dats,
                 par = Map(stats::setNames, coef, cn),
                 theta = if (one) thetas[[1L]] else thetas,
                 latent = if (one) latents[[1L]] else
                   if (all(vapply(latents, is.null, logical(1)))) NULL
                   else latents,
                 structural = psi,
                 hyper = truth$hyper,
                 n_sim = n_sim,
                 call = match.call()),
            class = "StatmodSim")
}


#' The Covariate Generators of a Simulation
#'
#' @description
#' Validates `covariates` and returns it, or `NULL`.
#'
#' @details
#' Each entry is a function of the observation count, as an entry of
#' `par` may be, so the two arguments read alike. A value that is not a
#' function is rejected, since a constant column belongs in `data` and the
#' covariates of this argument are drawn afresh at every replicate.
#'
#' @param covariates A named list, or `NULL`.
#'
#' @return The list, or `NULL`.
#'
#' @seealso [rstatmod()]
#'
#' @keywords internal
check_covariates <- function(covariates) {
  if (is.null(covariates)) return(NULL)
  if (!is.list(covariates) || is.null(names(covariates)) ||
      !all(nzchar(names(covariates)))) {
    stop("'covariates' must be a named list, one entry per covariate.",
         call. = FALSE)
  }
  ok <- vapply(covariates, is.function, logical(1))
  if (!all(ok)) {
    stop(sprintf(paste0("'covariates$%s' is not a function. Each entry is a ",
                        "function of the\n  observation count, drawn afresh ",
                        "at every replicate; a column that does not\n  ",
                        "change belongs in 'data'."),
                 names(covariates)[[which(!ok)[[1L]]]]), call. = FALSE)
  }
  covariates
}


#' Draw One Replicate's Covariates
#'
#' @description
#' Evaluates each generator at the observation count and writes the columns
#' into the data frame.
#'
#' @details
#' A generator that returns the wrong number of values is reported and not
#' recycled, for the reason given for a coefficient function: R would recycle
#' the vector without a message and the replicate would be of another model.
#' A column of `data` under the same name is overwritten, because naming the
#' column in `covariates` requests that it be drawn.
#'
#' @param covariates The generators.
#' @param n The observation count.
#' @param data The frame into which the columns are written.
#'
#' @return A data frame.
#'
#' @seealso [rstatmod()]
#'
#' @keywords internal
draw_covariates <- function(covariates, n, data) {
  for (v in names(covariates)) {
    col <- covariates[[v]](n)
    if (length(col) != n) {
      stop(sprintf(paste0("'covariates$%s' returned %d value%s where the ",
                          "simulation has %d\n  observation%s."),
                   v, length(col), if (length(col) == 1L) "" else "s", n,
                   if (n == 1L) "" else "s"), call. = FALSE)
    }
    data[[v]] <- col
  }
  data
}


#' @title Printing a Simulation
#' @name print.StatmodSim
#' @description
#' Prints the model that was written, the size of what came out, and the truth
#' behind it.
#' @param x A [rstatmod()] result.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @seealso [rstatmod()]
#' @examples
#' set.seed(1)
#' rstatmod(y ~ x, distributions7::gaussian1_distrib(),
#'          data.frame(x = runif(20)))
#' @export
print.StatmodSim <- function(x, ...) {
  cat("A simulation from a written model\n\n")
  one <- identical(x$n_sim, 1L) || is.null(x$n_sim)
  d1 <- if (one) x$data else x$data[[1L]]
  nc <- ncol(d1)
  cat(sprintf("  %d observations, %d column%s%s\n", nrow(d1), nc,
              if (identical(nc, 1L)) "" else "s",
              if (one) "" else sprintf(", %d replicates", x$n_sim)))
  for (p in names(x$par)) {
    v <- x$par[[p]]
    if (!length(v)) next
    cat(sprintf("  %-10s %s\n", p, sim_values(v)))
  }
  if (length(x$structural)) {
    cat("\n  the term's own parameters\n")
    psi <- unlist(x$structural)
    for (q in sim_groups(names(psi))) {
      cat(sprintf("  %-18s %s\n", q$label, sim_values(psi[q$idx])))
    }
  }
  if (!is.null(x$hyper) && nrow(x$hyper)) {
    cat("\n  the hyperparameters\n")
    for (i in seq_len(nrow(x$hyper))) {
      h <- x$hyper[i, ]
      cat(sprintf("  %-18s %s%s\n",
                  paste0(short_keys(h$term),
                         if (nzchar(h$term)) ":" else "", h$name),
                  format(round(h$value, 4)),
                  if (isTRUE(h$held)) "  (held)" else ""))
    }
  }
  if (!is.null(x$latent)) {
    l1 <- if (one) x$latent else x$latent[[1L]]
    cat(sprintf("\n  latent: %s\n",
                if (is.data.frame(l1))
                  sprintf("%d drawn positions", nrow(l1))
                else sprintf("%d values", length(l1))))
  }
  cat(if (one)
        "\n  the data is in $data, the truth in $par and $theta\n"
      else
        "\n  the data sets are in $data[[k]], the truth in $par\n")
  invisible(x)
}


#' A Vector of a Simulated Truth, Short Enough to Read
#'
#' @description
#' Returns the values where there are six or fewer, and otherwise the first
#' two values with the count and the standard deviation.
#'
#' @details
#' A random effect over a hundred groups is a hundred numbers, and printing
#' them hides the few numbers that a reader needs. For a block of that size
#' the printed form gives the count and the spread.
#'
#' @param v A numeric vector.
#'
#' @return A single string.
#'
#' @seealso [print.StatmodSim()]
#'
#' @keywords internal
sim_values <- function(v) {
  if (length(v) <= 6L) {
    return(paste(format(round(v, 4)), collapse = "  "))
  }
  sprintf("%s ... %d values, sd %s",
          paste(format(round(v[1:2], 4)), collapse = "  "), length(v),
          format(round(stats::sd(v), 4)))
}


#' The Groups That a Printed Truth Collapses
#'
#' @description
#' Returns one entry per scalar parameter and one per family of names that
#' share a stem, in the order in which the names arrive.
#'
#' @details
#' A structural term's parameters are `omega`, `kappa1` and the like where a
#' parameter is not developed, and `omega.(Intercept)`, `omega.random.1` and
#' so on where it is. The first are read one by one and the second are a
#' block, so the printed form is grouped at the first dot. The scalars then
#' come out as groups of one and need no special case.
#'
#' @param nm The names.
#'
#' @return A list of entries, each with `label` and `idx`.
#'
#' @seealso [print.StatmodSim()]
#'
#' @keywords internal
sim_groups <- function(nm) {
  stem <- sub("[.].*$", "", nm)
  # a development's coordinates share the parameter's own name, so the stem
  # groups them and a scalar is a group of one
  lapply(unique(stem), function(g) {
    idx <- which(stem == g)
    lab <- if (length(idx) == 1L) nm[[idx]] else
      sprintf("%s (%d)", g, length(idx))
    list(label = lab, idx = idx)
  })
}


#' A Positive Whole Number
#'
#' @description
#' Tests that a value is a single finite number, at least one, with no
#' fractional part, before it is converted with [as.integer()], which would
#' truncate it.
#'
#' @param v The value.
#'
#' @return `TRUE` or `FALSE`.
#'
#' @keywords internal
is_count <- function(v) {
  is.numeric(v) && length(v) == 1L && is.finite(v) && v >= 1 && v == round(v)
}


#' The Data Frame Against Which a Simulation Runs
#'
#' @description
#' Resolves `data` and `n` into one data frame of covariates.
#'
#' @details
#' A model may have no covariates at all (a pure time series is an example),
#' and then a row count is all that a simulation needs, so `data` and `n` are
#' alternatives and one of them must be given. Where both are given they must
#' agree. The function checks this and does not resolve a disagreement in
#' favor of one of them, since a disagreement points to a mistake of the
#' caller. `n` must be a positive whole number.
#'
#' @param data A data frame or `NULL`.
#' @param n A row count or `NULL`.
#' @param covariates The generators, or `NULL`. Where there are any,
#'   their columns are written in.
#'
#' @return A data frame.
#'
#' @seealso [rstatmod()]
#'
#' @keywords internal
rstatmod_data <- function(data, n, covariates = NULL) {
  if (is.null(data)) {
    if (is.null(n)) {
      stop("give 'data', 'n' with 'covariates', or 'n' alone where the ",
           "model has no covariates.", call. = FALSE)
    }
    if (!is_count(n)) {
      stop("'n' must be a positive whole number.", call. = FALSE)
    }
    n <- as.integer(n)
    out <- data.frame(row.names = seq_len(n))
    if (!is.null(covariates)) out <- draw_covariates(covariates, n, out)
    return(out)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!is.null(n) && (!is_count(n) || !identical(as.integer(n), nrow(data)))) {
    stop(sprintf("'n' is %s and 'data' has %d rows.",
                 format(n), nrow(data)), call. = FALSE)
  }
  if (!is.null(covariates)) {
    data <- draw_covariates(covariates, nrow(data), data)
  }
  data
}


#' The Column to Which a Simulated Response Is Written
#'
#' @description
#' Returns the name on the formula's left-hand side, which must be a symbol.
#'
#' @details
#' A transformed response is rejected. Under `log(y) ~ x` the model
#' generates values of `log(y)`, so a column called `y` would hold the wrong
#' quantity and one called `"log(y)"` would not be the name that the formula
#' reads back. A censored response is rejected for the reason that
#' [statmod()] rejects one: there is no censored likelihood with which to fit
#' it back.
#'
#' @param response The formula's left-hand side.
#'
#' @return A single string.
#'
#' @seealso [rstatmod()]
#'
#' @keywords internal
rstatmod_response_name <- function(response) {
  if (is.name(response)) return(as.character(response))
  if (is.call(response) && identical(as.character(response[[1L]]), "cens")) {
    stop("a censored response is not simulated: this package has no ",
         "censored likelihood\n  with which to fit it back, so the data would be ",
         "of a model that it cannot read. Simulate\n  the uncensored response ",
         "and censor it yourself.", call. = FALSE)
  }
  stop(sprintf(paste0("the left-hand side must be a symbol, and it is '%s'.",
                      "\n  The model generates that quantity, so there is ",
                      "no column to which it can be written:\n  write the ",
                      "model in the variable that you want simulated."),
               deparse(response)), call. = FALSE)
}


#' The Predictor and Parameters of a Simulation
#'
#' @description
#' Returns the linear predictors, the distribution's parameters and, where a
#' term has state, the response that it drew and the latent quantity behind
#' it.
#'
#' @details
#' With no structural term the result is that of [statmod_eta()]: the
#' simulated data come from the assembly that a fit reads, so a term whose
#' block moves with its coefficients contributes its value and not a
#' linearization.
#'
#' With one, the static part is assembled the same way and the term completes
#' the predictor through [modelterms7::term_simulate()]. A term that draws the
#' response as it goes returns it; one that does not returns `NULL` there and
#' the caller draws at the predictor.
#'
#' @param spec The specification.
#' @param design Its design.
#' @param coef The coefficients.
#' @param psi The structural term's own parameters, drawn once by
#'   [rstatmod_truth()] so that the replicates share them, or `NULL` to read
#'   the design's own structural state.
#'
#' @return A list with `eta`, `theta`, `y`, `latent`
#'   and `structural`.
#'
#' @seealso [rstatmod()], [statmod_eta()]
#'
#' @keywords internal
rstatmod_eta <- function(spec, design, coef, psi = NULL) {
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  su <- attr(design, "structural")
  if (!length(su)) {
    ep <- statmod_eta(spec, design, coef)
    return(list(eta = ep$eta, theta = ep$theta, y = NULL, latent = NULL,
                structural = NULL))
  }
  if (length(su) > 1L) {
    stop("a formula carries at most one structural term, the joint filter of several gas() terms counting as one.",
         call. = FALSE)
  }
  u <- su[[1L]]

  # the static part, assembled as statmod_eta assembles it
  d2 <- statmod_design_at(spec, coef, design)
  eta <- stats::setNames(vector("list", length(params)), params)
  theta <- eta
  for (p in params) {
    d <- d2[[p]]
    e <- if (d$npar == 0L) rep(0, spec@n_obs) else
      as.numeric(d$X %*% coef[[p]])
    if (!is.null(d$adj)) e <- e + d$adj
    off <- spec@offsets[[p]]
    if (!is.null(off)) e <- e + off
    eta[[p]] <- e
    theta[[p]] <- linkfunctions7::linkinv(links[[p]], e)
  }

  tm <- spec@terms[[u$param]][[u$term]]
  # drawn once by rstatmod_truth(), so the replicates share it
  if (is.null(psi)) {
    psi <- structural_psi(tm, statmod_structural_state(design)$zeta[[u$term]])
  }
  if (is_joint_filter(tm)) {
    return(rstatmod_eta_joint(spec, tm, eta, theta, psi))
  }
  p <- u$param
  lk <- links[[p]]
  n <- spec@n_obs
  theta_n <- lapply(params, function(q) rep_len(theta[[q]], n))
  names(theta_n) <- params
  at <- function(i) lapply(theta_n, function(v) v[[i]])
  # the response at one observation, drawn at the predictor the term has
  # reached there: everything but this equation's parameter is the static
  # value, and this one is the predictor's own inverse link
  draw <- function(e, i) {
    th <- at(i)
    th[[p]] <- linkfunctions7::linkinv(lk, e)
    as.numeric(distributions7::distrib_rng(spec@distrib, 1L, th))
  }
  k <- distributions7::distrib_kernel(spec@distrib, p)
  # the filter is driven by the scaled score where gas(scaling = ) is not
  # zero, so the simulated path runs the same recursion the fit did
  d_sc <- filter_scaling(tm)
  score <- if (d_sc == 0) {
    function(y, e, i) as.numeric(k$score(y, at(i), e))
  } else {
    function(y, e, i) {
      th <- at(i)
      as.numeric(k$score(y, th, e)) *
        as.numeric(k$information(y, th, e))^(-d_sc)
    }
  }
  curvature <- if (d_sc == 0) {
    function(y, e, i) as.numeric(k$curvature(y, at(i), e))
  } else {
    function(y, e, i) {
      th <- at(i)
      s <- as.numeric(k$score(y, th, e))
      info <- as.numeric(k$information(y, th, e))
      as.numeric(k$curvature(y, th, e)) * info^(-d_sc) -
        d_sc * s * info^(-d_sc - 1) * as.numeric(k$dinformation(y, th, e))
    }
  }
  res <- modelterms7::term_simulate(tm, psi, eta[[p]], draw,
                                    score = score, curvature = curvature)
  eta[[p]] <- res$eta
  theta[[p]] <- linkfunctions7::linkinv(lk, res$eta)
  list(eta = eta, theta = theta, y = res$y, latent = res$latent,
       structural = as.list(psi))
}


#' Draw the Whole Truth of a Simulation
#'
#' @description
#' Draws every coefficient, every hyperparameter and, where the formula
#' carries a structural term, that term's own parameters, once, and then
#' overwrites them with the values that `par` supplies.
#'
#' @details
#' Each quantity is drawn by the component that defines its meaning. A
#' coefficient of a design column has no other owner and is drawn from a
#' normal of width `sd`. A coordinate that some penalty covers is drawn from
#' that penalty read as a prior, so a Gaussian random effect gives Gaussian
#' effects and a lasso gives Laplace ones; the prior's own scale is drawn too,
#' on the chart that its penalty carries. A structural term's parameters are
#' drawn by the term, on the chart that each one rides.
#'
#' The truth is drawn once whatever `n_sim` is, so the replicates differ in
#' what is random and not in what is being estimated.
#'
#' A hyperparameter that the term holds is used as given and not drawn: in
#' `s(x, bspline_smooth(), hyper = c(lambda = 2))` the simulation uses the
#' smoothing parameter that the term specifies.
#'
#' @param spec The specification, built against a placeholder response.
#' @param design Its design.
#' @param par A named list, or `NULL`. See [rstatmod()] for what a key may be.
#' @param sd The width of the draws.
#'
#' @return A list with `coef`, one numeric vector per distribution parameter;
#'   `hyper`, a data frame of the hyperparameters drawn or held; and `psi`,
#'   the structural term's own parameters or `NULL`.
#'
#' @seealso [rstatmod()], [penalties7::penalty_draw()],
#'   [modelterms7::term_draw()]
#'
#' @keywords internal
rstatmod_truth <- function(spec, design, par, sd) {
  params <- spec@distrib@params
  # AN EQUATION `par` FIXES WHOLE IS NOT DRAWN, so a simulation whose truth
  # is entirely written consumes the same random numbers it always did and
  # gives the same data. Every finer key came later and has nothing to
  # reproduce.
  fixed_eq <- if (is.list(par)) intersect(names(par), params) else character(0)
  coef <- stats::setNames(lapply(params, function(p) {
    if (p %in% fixed_eq) rep(NA_real_, design[[p]]$npar)
    else stats::rnorm(design[[p]]$npar, 0, sd)
  }), params)
  su <- attr(design, "structural")
  tm <- if (length(su))
    spec@terms[[su[[1L]]$param]][[su[[1L]]$term]] else NULL
  psi_z <- if (is.null(tm)) NULL else modelterms7::term_draw(tm, sd = sd)

  units <- statmod_penalized(spec, design)
  rows <- list()
  done <- list()
  for (u in units) {
    # a Jeffreys prior is improper and is not drawn from
    if (isTRUE(u$jeffreys)) next
    # A PRIOR OVER COORDINATES THAT RIDE A CHART IS CENTRED HALF AS FAR
    # from its bound, by the rule term_draw() halves its own width by: a
    # chart carries a scale of one onto most of its parameter's range.
    # Measured on a partial-autocorrelation chart, a prior scale of 1 puts
    # 10.1 per cent of the persistences past 0.95 and one of 0.5 puts 0.9
    # per cent there, the middle ninety running -0.82 to 0.82.
    th <- rstatmod_hyper(u$penalty, sd,
                         unit = if (isTRUE(u$structural)) 0.5 else 1)
    if (length(u$fixed)) th[names(u$fixed)] <- u$fixed
    v <- penalties7::penalty_draw(u$penalty, th)
    if (is.null(v) || all(is.na(v))) next
    tg <- unit_draw_targets(u, design, params)
    if (is.null(tg) || length(tg$pos) != length(v)) next
    ok <- !is.na(v)
    if (tg$space == "zeta") {
      if (!is.null(psi_z)) psi_z[tg$pos[ok]] <- v[ok]
    } else {
      for (p in setdiff(unique(tg$param[ok]), fixed_eq)) {
        j <- ok & tg$param == p
        coef[[p]][tg$pos[j]] <- v[j]
      }
    }
    # A HYPERPARAMETER IS REPORTED ONLY WHERE THE DRAW REACHED ITS
    # COORDINATES. Where it did not -- SCAD and MCP, which are densities of
    # nothing, a covariance class spanning two vectors -- those coefficients
    # came from the plain draw and no scale of the prior describes them, so
    # a row for it would name a truth the data do not have.
    done[[length(done) + 1L]] <- list(u = u, th = th, tg = tg, ok = ok, v = v)
  }

  # THE TERM DRAWS ITS OWN COORDINATES LAST, after the plain draw and after
  # the priors, so that what it writes is what survives. A break-point is
  # the case: it is a position on the covariate's axis and nothing else in
  # the model is measured in those units.
  td <- rstatmod_term_draw(spec, design, coef, params, fixed_eq, sd)
  coef <- td$coef
  # A HYPERPARAMETER WHOSE COORDINATES `par` WROTE IS NOT REPORTED. It is
  # still drawn, so the random stream and the data do not move, but it
  # governs nothing in the truth: `par = coef(fit)` simulates from a fitted
  # model, and a row for a smoothing parameter drawn beside it would name a
  # truth the data do not have.
  written <- par_written(par, design, params,
                         if (is.null(tm)) NULL else
                           names(structural_psi(tm, psi_z)))
  for (e in done) {
    if (targets_written(e$tg, e$ok, written)) next
    th <- e$th
    s <- term_drawn_scale(e$tg, e$ok, td$owned)
    # A HELD HYPERPARAMETER WINS OVER THE TERM'S OWN WIDTH, and the draw is
    # put back with it: a caller writing `psi ~ random(~1 | id, hyper =
    # c(sigma = 0.3))` has said what the break-points vary by, which the term
    # has no way of knowing and no business overruling. What the term keeps
    # there is the placement, every coordinate no prior covers.
    if (!is.null(s) && length(e$u$fixed)) {
      j <- e$ok
      for (p in setdiff(unique(e$tg$param[j]), fixed_eq)) {
        k <- j & e$tg$param == p
        coef[[p]][e$tg$pos[k]] <- e$v[k]
      }
      s <- NULL
    }
    # WHERE THE TERM OVERWROTE THEM the prior's own draw is no longer the
    # truth of those coefficients, so what is reported is the width the term
    # used, carried onto whatever the penalty calls its hyperparameter.
    if (!is.null(s)) th <- as.list(penalty_theta_start(e$u$penalty, s))
    for (h in names(th)) {
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = if (is.null(e$u$params)) e$u$param else
          paste(unique(e$u$params), collapse = ", "),
        term = e$u$key, name = h, value = as.numeric(th[[h]]),
        held = h %in% names(e$u$fixed), stringsAsFactors = FALSE)
    }
  }

  psi <- if (is.null(tm)) NULL else structural_psi(tm, psi_z)
  out <- rstatmod_named(coef, psi, design, params, par)
  hyper <- if (length(rows)) do.call(rbind, rows) else
    data.frame(parameter = character(0), term = character(0),
               name = character(0), value = numeric(0), held = logical(0),
               stringsAsFactors = FALSE)
  list(coef = out$coef, psi = out$psi, hyper = hyper)
}


#' The Coordinates That a Simulation's `par` Writes
#'
#' @description
#' Marks, for every distribution parameter and for the structural term's own
#' parameters, the positions that some key of `par` addresses.
#'
#' @details
#' It reads the keys through [resolve_par_key()], the function that
#' [rstatmod_named()] uses to write them, so a key reaches the same positions
#' in both. A key that reaches nothing signals the same error
#' here that it would signal there.
#'
#' @param par A named list, or `NULL`.
#' @param design The design.
#' @param params The distribution parameter names.
#' @param snm The structural parameter names, or `NULL`.
#'
#' @return A list with `beta`, one logical vector per distribution parameter,
#'   and `zeta`, a logical vector over the structural parameters.
#'
#' @seealso [rstatmod_truth()], [targets_written()]
#'
#' @keywords internal
par_written <- function(par, design, params, snm) {
  beta <- stats::setNames(lapply(params, function(p)
    logical(design[[p]]$npar)), params)
  zeta <- logical(length(snm))
  if (!is.list(par) || !length(par) || is.null(names(par))) {
    return(list(beta = beta, zeta = zeta))
  }
  cn <- stats::setNames(lapply(params, function(p) design[[p]]$coef_names),
                        params)
  for (k in names(par)) {
    if (!nzchar(k)) next
    tg <- resolve_par_key(k, params, cn, snm)
    if (identical(tg$space, "zeta")) {
      zeta[tg$pos] <- TRUE
    } else {
      beta[[tg$param]][tg$pos] <- TRUE
    }
  }
  list(beta = beta, zeta = zeta)
}


#' Whether `par` Wrote Every Coordinate That a Prior Drew
#'
#' @description
#' Returns `TRUE` where every coordinate that a penalty's draw reached was then
#' written by `par`, so that the penalty's hyperparameters govern nothing in
#' the truth.
#'
#' @param tg The targets, from [unit_draw_targets()].
#' @param ok Which of them the draw reached.
#' @param written The result of [par_written()].
#'
#' @return A single logical.
#'
#' @seealso [rstatmod_truth()]
#'
#' @keywords internal
targets_written <- function(tg, ok, written) {
  if (!any(ok)) return(FALSE)
  pos <- tg$pos[ok]
  if (identical(tg$space, "zeta")) return(all(written$zeta[pos]))
  prm <- tg$param[ok]
  all(vapply(seq_along(pos), function(i)
    isTRUE(written$beta[[prm[[i]]]][[pos[[i]]]]), logical(1)))
}


#' Let Every Term Draw the Coefficients That Only It Can Draw
#'
#' @description
#' Walks the terms of every equation and hands each the slice of the drawn
#' coefficients that its block owns, taking back whatever it replaced.
#'
#' @details
#' Almost every term replaces nothing: a slope is measured in the response's
#' units against a covariate's and a normal of the caller's width is as good
#' a truth as any. A break-point is the exception, because it is a position on
#' the covariate's own axis, and [modelterms7::term_coef_draw()] is the method
#' by which a term declares such coordinates.
#'
#' It runs last, after the plain draw and after the priors, because those
#' describe a coordinate that a term owns incorrectly. The coordinates
#' replaced are returned, so that a prior covering them is reported at the
#' width that the term used and not at the one that it drew.
#'
#' An equation that `par` fixes whole is skipped: its coefficients are the
#' caller's and are written over everything afterwards anyway.
#'
#' The walk is over the design's blocks and does not descend into a term's
#' subformulas, so a break-point developing another term's own parameter
#' (`nl(a ~ 0 + seg(x))`) is drawn by the plain rule. The sub-terms in a
#' break-point's own development are `linpar()` and `random()`, which replace
#' no coordinates, so the shapes that reach one are covered by the block
#' walk.
#'
#' @param spec The specification, read for its terms.
#' @param design The design, read for each term's columns.
#' @param coef The coefficients drawn so far, a list by parameter.
#' @param params The distribution parameter names.
#' @param fixed_eq The parameters that `par` fixes whole.
#' @param sd The width of the draws.
#'
#' @return A list with `coef`, the coefficients in which each term's own have
#'   been replaced, and `owned`, a list by parameter of named numeric vectors
#'   carrying one width per replaced position, where the names are the
#'   positions.
#'
#' @seealso [rstatmod_truth()], [modelterms7::term_coef_draw()]
#'
#' @keywords internal
rstatmod_term_draw <- function(spec, design, coef, params, fixed_eq, sd) {
  owned <- stats::setNames(vector("list", length(params)), params)
  for (p in setdiff(params, fixed_eq)) {
    blocks <- design[[p]]$blocks
    if (!length(blocks)) next
    got <- numeric(0)
    for (nm in names(blocks)) {
      tm <- spec@terms[[p]][[nm]]
      cols <- blocks[[nm]]
      if (is.null(tm) || !length(cols)) next
      d <- modelterms7::term_coef_draw(tm, coef[[p]][cols], sd = sd)
      if (is.null(d) || !length(d$index)) next
      coef[[p]][cols] <- as.numeric(d$coef)
      got[as.character(cols[d$index])] <- as.numeric(d$scale)
    }
    owned[[p]] <- got
  }
  list(coef = coef, owned = owned)
}


#' The Width at Which a Term Drew a Penalty's Coordinates
#'
#' @description
#' Returns one number where every coordinate that a penalty covers was written
#' by the term that owns it, and `NULL` otherwise.
#'
#' @details
#' A penalty whose coordinates were overwritten by a term no longer describes
#' them, so the truth reported for it is the width at which the term drew.
#' The test is that every covered coordinate was owned: a penalty spanning
#' owned and unowned ones together is described by neither width, and
#' reporting either would name a truth that part of the coefficients do not
#' have.
#'
#' @param tg The unit's targets, as [unit_draw_targets()] returns them.
#' @param ok Which of them the prior's draw reached.
#' @param owned The `owned` element of [rstatmod_term_draw()].
#'
#' @return A single number, or `NULL`.
#'
#' @seealso [rstatmod_term_draw()]
#'
#' @keywords internal
term_drawn_scale <- function(tg, ok, owned) {
  if (is.null(tg) || tg$space != "beta" || !any(ok)) return(NULL)
  pp <- tg$param[ok]
  ps <- tg$pos[ok]
  sc <- numeric(0)
  for (i in seq_along(ps)) {
    o <- owned[[pp[i]]]
    k <- as.character(ps[i])
    if (is.null(o) || !(k %in% names(o))) return(NULL)
    sc <- c(sc, o[[k]])
  }
  if (!length(sc)) NULL else stats::median(sc)
}


#' Where a Penalty's Coordinates Sit in the Vectors That a Simulation Draws
#'
#' @description
#' Returns one target per coordinate that the penalty covers, in the penalty's
#' own order. Each target holds the vector to which the coordinate belongs,
#' the parameter whose coefficients contain it and its position within them.
#'
#' @details
#' A penalty over a structural term's own parameters is read among those
#' parameters, on the unconstrained scale, and a penalty over an equation's
#' coefficients among those coefficients. A covariance class spans several
#' members and its penalty reads them interleaved group by group, in the
#' order that its `index` already records, so the stacked positions are
#' turned back into a parameter and a column with the offsets that the design
#' gives.
#'
#' A class whose halves are in both vectors at once is not drawn from its
#' prior: those coordinates fall back to the plain draw, which is
#' [rstatmod()]'s rule for every coordinate that a prior does not reach.
#'
#' @param u One entry of [statmod_penalized()].
#' @param design The design.
#' @param params The distribution parameter names.
#'
#' @return A list with `space`, `param` and `pos`, or `NULL` where the entry
#'   cannot be addressed.
#'
#' @seealso [rstatmod_truth()]
#'
#' @keywords internal
unit_draw_targets <- function(u, design, params) {
  if (isTRUE(u$mixed)) return(NULL)
  if (isTRUE(u$structural)) {
    if (!length(u$cols)) return(NULL)
    return(list(space = "zeta", param = rep(NA_character_, length(u$cols)),
                pos = as.integer(u$cols)))
  }
  if (!is.null(u$cols)) {
    return(list(space = "beta", param = rep(u$param, length(u$cols)),
                pos = as.integer(u$cols)))
  }
  # a covariance class: stacked positions, turned back into (parameter,
  # column) with the same offsets statmod_penalized() built them from
  if (is.null(u$index) || !length(u$index)) return(NULL)
  npar <- vapply(design[params], function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  a <- findInterval(u$index - 1L, offs)
  list(space = "beta", param = params[a], pos = as.integer(u$index - offs[a]))
}


#' Draw One Penalty's Hyperparameters
#'
#' @description
#' Draws a value per hyperparameter on the chart that the penalty carries for
#' it, around the neutral point that [penalty_theta_start()] gives.
#'
#' @details
#' The draw is on the unconstrained scale so that whatever comes out is
#' admissible: a scale stays positive, a correlation stays inside its
#' interval, and a log-Cholesky coordinate is free on the line already. The
#' width is half of `sd` for the reason that [modelterms7::term_draw()] halves
#' it: a chart maps a width of one onto most of its parameter's range.
#'
#' `unit` moves the center away from a one-sided bound and affects nothing
#' else, so a hyperparameter that is unbounded on both sides (a coordinate of
#' a multivariate prior's own chart) is drawn around the neutral point of that
#' chart whatever the term is. This is also the point from which a fit
#' starts; on `dr_prod` it corresponds to unit standard deviations and no
#' correlation.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param sd The width of the draws.
#' @param unit How far from a one-sided bound the draw is centered, passed to
#'   [penalty_theta_start()]. It is one half for a prior over a structural
#'   term's own parameters, which ride charts.
#'
#' @return A named list, empty for a penalty with no hyperparameters.
#'
#' @seealso [rstatmod_truth()], [penalty_theta_start()]
#'
#' @keywords internal
rstatmod_hyper <- function(pen, sd, unit = 1) {
  st <- penalty_theta_start(pen, unit)
  if (!length(st)) return(list())
  nm <- names(st)
  stats::setNames(lapply(nm, function(h) {
    lk <- pen@link_params[[h]]
    z <- stats::rnorm(1, 0, sd / 2)
    if (is.null(lk)) return(as.numeric(st[[h]]) + z)
    linkfunctions7::linkinv(lk, linkfunctions7::linkfun(lk,
                                                        as.numeric(st[[h]])) + z)
  }), nm)
}


#' Hold What a Simulation's Caller Named
#'
#' @description
#' Writes the entries of `par` over the drawn truth, and each entry is
#' addressed by the name under which the result reports it.
#'
#' @details
#' A single namespace covers the whole model, so the caller names the
#' quantities to fix and everything else stays drawn. A key may be a
#' distribution parameter, which is that whole equation; one of its
#' coefficients or a group of them, written `parameter.coefficient`; a
#' structural term's own parameter; or a group of those. A group is a name
#' that its members extend at a dot, so `omega.random` is every deviation in
#' a development of `omega` and `mu.s(x)` is every coordinate of that smooth.
#'
#' A value may be a vector of the group's own length, a single number used
#' for all of them, or a function of the count. A function writes a
#' structured truth: `function(k) rnorm(k, 0, 0.4)` is a random
#' effect at a scale of its own and `function(k) c(2, -1.5, rep(0, k - 2))`
#' is a sparse truth for a lasso to find.
#'
#' A structural parameter is named on the scale of
#' [modelterms7::term_params()]: a loading is the loading and not its
#' logarithm.
#'
#' @param coef The drawn coefficients.
#' @param psi The drawn structural parameters, or `NULL`.
#' @param design The design.
#' @param params The distribution parameter names.
#' @param par A named list, or `NULL`.
#'
#' @return A list with `coef` and `psi`.
#'
#' @seealso [rstatmod()], [rstatmod_truth()]
#'
#' @keywords internal
rstatmod_named <- function(coef, psi, design, params, par) {
  if (is.null(par) || !length(par)) return(list(coef = coef, psi = psi))
  if (!is.list(par) || is.null(names(par)) || !all(nzchar(names(par)))) {
    stop("'par' must be a named list. Name a distribution parameter, one of ",
         "its\n  coefficients as 'parameter.coefficient', or a structural ",
         "term's own\n  parameter; a quantity that is not named is drawn.",
         call. = FALSE)
  }
  cn <- lapply(params, function(p) design[[p]]$coef_names)
  names(cn) <- params
  snm <- names(psi)
  # A NAME THAT REACHES TWO THINGS IS REFUSED rather than resolved by an
  # order of precedence nobody would remember. It is the exact collision --
  # a structural parameter called `mu` -- and the one a group would create,
  # a structural parameter called `mu.something` against the coefficient
  # `something` of the equation `mu`. Neither occurs among the shipped
  # families and terms; both are one line to rule out.
  # `startsWith()` rejects a NULL, which is what a model with no structural
  # term has here, so the second half is asked only where there is anything
  # to ask it of
  clash <- if (!length(snm)) character(0) else
    c(intersect(params, snm),
      unlist(lapply(params, function(p)
        snm[startsWith(snm, paste0(p, "."))])))
  if (length(clash)) {
    stop(sprintf(paste0("'%s' is reachable both as a coefficient of a ",
                        "distribution parameter\n  and as a parameter of ",
                        "the structural term, so the name is ambiguous in ",
                        "'par'."), clash[[1L]]),
         call. = FALSE)
  }
  for (k in names(par)) {
    tg <- resolve_par_key(k, params, cn, snm)
    v <- par_values(par[[k]], length(tg$pos), k)
    if (identical(tg$space, "zeta")) {
      for (i in seq_along(tg$pos)) psi[[tg$pos[[i]]]] <- v[[i]]
    } else {
      coef[[tg$param]][tg$pos] <- v
    }
  }
  list(coef = coef, psi = psi)
}


#' Resolve One Key of a Simulation's `par`
#'
#' @description
#' Turns a name into the vector and the positions that it addresses.
#'
#' @details
#' A key is tried in this order: a distribution parameter; a structural
#' parameter, then a group of structural parameters; a coefficient of an
#' equation, written `parameter.coefficient`, then a group of the coefficients
#' of that equation. A name that is both a coefficient and the stem of others
#' therefore addresses the coefficient. A group is recognized at a dot and not
#' by any prefix, since a coefficient name may begin with the letters of
#' another.
#'
#' A key that addresses nothing is reported with what the model does carry. A
#' misspelled name that was ignored would leave the quantity drawn at random,
#' and the simulation would be of another model than the one intended.
#'
#' @param key The name.
#' @param params The distribution parameter names.
#' @param cn A named list of coefficient names, one per parameter.
#' @param snm The structural parameter names, or `NULL`.
#'
#' @return A list with `space`, `param` and `pos`.
#'
#' @seealso [rstatmod_named()]
#'
#' @keywords internal
resolve_par_key <- function(key, params, cn, snm) {
  if (key %in% params) {
    return(list(space = "beta", param = key,
                pos = seq_len(length(cn[[key]]))))
  }
  if (length(snm)) {
    if (key %in% snm) {
      return(list(space = "zeta", param = NA_character_,
                  pos = match(key, snm)))
    }
    g <- which(startsWith(snm, paste0(key, ".")))
    if (length(g)) {
      return(list(space = "zeta", param = NA_character_, pos = g))
    }
  }
  for (p in params) {
    pre <- paste0(p, ".")
    if (!startsWith(key, pre)) next
    rest <- substring(key, nchar(pre) + 1L)
    nmp <- cn[[p]]
    if (rest %in% nmp) {
      return(list(space = "beta", param = p, pos = match(rest, nmp)))
    }
    g <- which(startsWith(nmp, paste0(rest, ".")))
    if (length(g)) return(list(space = "beta", param = p, pos = g))
  }
  some <- function(x) {
    if (!length(x)) return("none")
    paste(c(utils::head(x, 6L), if (length(x) > 6L) "...") , collapse = ", ")
  }
  # AN EQUATION MAY HAVE NO COLUMNS AT ALL -- `y ~ 0 + gas(...)` is the
  # ordinary case -- so the example is taken from the first that has one and
  # is left out where none does. A message that itself fails is worse than a
  # terse one.
  with_cols <- params[vapply(params, function(p) length(cn[[p]]) > 0L,
                             logical(1))]
  eg <- if (length(with_cols))
    sprintf("\n  A coefficient is named 'parameter.coefficient', as in '%s'.",
            paste0(with_cols[[1L]], ".", cn[[with_cols[[1L]]]][[1L]])) else ""
  stop(sprintf(paste0("'par' names '%s', which addresses nothing in this ",
                      "model.\n  The distribution parameters are: %s.%s%s"),
               key, paste(params, collapse = ", "), eg,
               if (length(snm))
                 sprintf("\n  The structural term's own parameters are: %s.",
                         some(snm)) else ""),
       call. = FALSE)
}


#' The Values Represented by One Entry of `par`
#'
#' @description
#' Returns a numeric vector of the length that the key addresses, from a
#' vector, a single number or a function of the count.
#'
#' @details
#' A function that returns the wrong number of values is reported and not
#' recycled: R would recycle the vector without a message and the simulation
#' would be of another model. A function that returns a missing value is
#' reported as well.
#'
#' @param e The entry.
#' @param k How many values the key addresses.
#' @param key The name, for the message.
#'
#' @return A numeric vector of length `k`.
#'
#' @seealso [rstatmod_named()]
#'
#' @keywords internal
par_values <- function(e, k, key) {
  if (is.function(e)) {
    v <- as.numeric(e(k))
    if (length(v) == k && anyNA(v)) {
      stop(sprintf("'par$`%s`' is a function and returned a missing value.",
                   key), call. = FALSE)
    }
    if (length(v) != k) {
      stop(sprintf(paste0("'par$`%s`' is a function and returned %d value%s ",
                          "where that name\n  addresses %d."),
                   key, length(v), if (length(v) == 1L) "" else "s", k),
           call. = FALSE)
    }
    return(v)
  }
  v <- as.numeric(e)
  if (length(v) == 1L && k > 1L) return(rep(v, k))
  if (length(v) != k) {
    stop(sprintf("'par$`%s`' has length %d where that name addresses %d.",
                 key, length(v), k), call. = FALSE)
  }
  v
}
