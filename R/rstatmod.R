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
#' Data whose truth is known: write the model, draw from it, fit it back, and
#' see whether the fit recovers what was put in. That is the shape of a
#' simulation study, and of a check on a term one has just written.
#'
#' This is not [stats::simulate()], which draws from a model already fitted.
#' The `r` prefix is R's own for a random draw, so the two names stay apart.
#'
#' A covariate needs no declaring. A factor becomes its contrasts and a
#' numeric stays itself, the design coming from the same interpreter a fit
#' uses.
#'
#' # The truth comes back beside the data
#'
#' A simulation study compares against the coefficients, the parameters they
#' gave and whatever a term with state drew, so the result is a list holding
#' all of them. Pass `sim$data` where a data frame is wanted.
#'
#' They were attributes of the data frame until version 0.88.0, and that was
#' worse than it looks: an attribute survives a row subset without being
#' subset itself, so `sim[1:10, ]` silently kept a `theta` of the original
#' length, while `subset()` and `merge()` dropped it altogether.
#'
#' # The predictor is assembled as a fit assembles it
#'
#' Through [statmod_design_at()], so the simulated data come from the model
#' that was written, never from a linearization of it.
#'
#' A term whose block moves with its coefficients, [modelterms7::seg()],
#' [modelterms7::jseg()] or [modelterms7::nl()], contributes
#' `term_value()` at the coefficients supplied, not its block times them. The
#' two differ by the whole nonlinearity, and an earlier version of this
#' function used the second: measured, the gap is 3 on a `seg()`, 4.05 on an
#' `nl()` and a missing value on a `jseg()`.
#'
#' # Data, or a row count
#'
#' `data` carries the covariates. A model with none, a pure time series say,
#' needs only `n`. One of the two must be given, and where both are given
#' they must agree.
#'
#' # Fixed covariates or drawn ones
#'
#' `covariates` takes one function of the observation count per column, as
#' `par` takes one per equation, and they are drawn afresh at every
#' replicate.
#'
#' The choice decides what a study measures. With `data` the covariates are
#' the same throughout, so what is measured is the estimator's behavior
#' **conditional** on that design; with `covariates` it is measured over the
#' design as well. Neither is more correct, and a study should say which it
#' ran. Measured on a simple regression at \eqn{n = 40}, the slope's standard
#' deviation is 0.1465 under the first and 0.1443 under the second, so the
#' distinction is about what a study claims, never about a large numerical
#' difference.
#'
#' A drawn factor is refused unless its levels are fixed. Drawn freely it
#' loses a level on some replicate, and the coefficients drawn against the
#' first design would then be recycled into a different model.
#'
#' # Several replicates
#'
#' `n_sim` draws that many data sets. The truth is drawn **once** and shared:
#' a study over replicates measures the variability of an estimator at a set
#' of parameters, so the replicates differ in what is random, never in what
#' is being estimated. Varying the truth as well is a loop over calls and
#' reads differently.
#'
#' With `n_sim > 1` the per-replicate fields, `data`, `theta` and `latent`,
#' come back as lists of that length, while `par` and `structural` stay
#' single.
#'
#' # The coefficients
#'
#' `par = NULL` draws every one from `rnorm(1, 0, sd)`, which on the link
#' scale gives predictors of order one. A named list fixes them instead, one
#' entry per distribution parameter, and an entry may be:
#'
#' - a numeric vector, as long as that equation has coefficients;
#' - a single number, used for every coefficient of the equation;
#' - a **function** of the coefficient count, called once and returning that
#'   many values.
#'
#' The function is how a structured truth is written without a vocabulary for
#' it. `function(k) rnorm(k, 0, 0.3)` is a random effect with its own
#' standard deviation, and `function(k) c(1.5, -2, rep(0, k - 2))` is a
#' sparse truth for a lasso to find. A function answering with the wrong
#' count is refused, R being willing to recycle it into a different model. A
#' parameter left out of the list is drawn.
#'
#' # A term with state
#'
#' Simulated through [modelterms7::term_simulate()], so the recursion that
#' generates is the term's own. A score-driven term draws the response as it
#' runs, its level at one time driven by the score at the time before; a
#' latent chain draws its path from the stationary law the likelihood is
#' written with; a marginal break-point term draws each group's positions
#' from their prior. What each drew comes back in `latent`, and that is what
#' a recovery check compares against.
#'
#' Such a term's own parameters are not coefficients of any equation, so they
#' are named through `structural`, never through `par`, on the scale
#' [modelterms7::term_params()] names: a loading is the loading, not its
#' logarithm, a persistence is the partial autocorrelation the chart carries.
#' A formula holds at most one such term, so no key is needed.
#'
#' Left unnamed they take the term's own starting values, which are
#' deliberately weak. A score-driven term starts at a loading near 0.1, and
#' the series then has almost no dynamics: measured, its level ranged over
#' 0.64 against 2.40 at named parameters. Name them, or the simulation is of
#' a model close to the one with no term at all.
#'
#' # The response's name
#'
#' The formula's left-hand side, which must be a symbol. `log(y) ~ x` is
#' refused: the model generates values of `log(y)`, and no column could
#' honestly be called either name. A censored response is refused too, for
#' the reason [statmod()] refuses one.
#'
#' @param formula The model formula, as [statmod()] takes it, with the
#'   equations separated by `|`. Its left-hand side must be a symbol.
#' @param distrib A \pkg{distributions7} distribution object, which decides
#'   how many equations there are and what is drawn from.
#' @param data A data frame of covariates, or `NULL` where `n` is given.
#' @param n The number of observations, where `data` is `NULL`. One of the
#'   two is required.
#' @param n_sim How many data sets to draw. `1` by default.
#' @param par Optional named list, one entry per distribution parameter, each
#'   a numeric vector, a single number or a function of the coefficient
#'   count. See the details. `NULL` draws every coefficient.
#' @param structural Optional named list of a structural term's own
#'   parameters, on the scale [modelterms7::term_params()] names. Strongly
#'   recommended when the formula carries such a term.
#' @param sd The standard deviation of the drawn coefficients, `1` by
#'   default. Read only for the coefficients `par` does not fix.
#' @param offsets Optional named list of offsets, one per parameter, as
#'   [statmod()] takes them.
#' @param covariates Optional named list of functions of the observation
#'   count, one per covariate, drawn afresh at every replicate. A drawn
#'   factor must have its levels fixed.
#'
#' @return An object of class `"StatmodSim"`, a list of seven:
#'   \describe{
#'     \item{`data`}{the data frame with the response column added, named
#'       after the formula's left-hand side.}
#'     \item{`par`}{the coefficients used, drawn or given, named as the
#'       design names them.}
#'     \item{`theta`}{the distribution's parameters at every observation, a
#'       named list with one vector per parameter.}
#'     \item{`latent`}{what a term with state drew, or `NULL`.}
#'     \item{`structural`}{such a term's own parameters, or `NULL`.}
#'     \item{`n_sim`}{as supplied.}
#'     \item{`call`}{the matched call.}
#'   }
#'   With `n_sim > 1` the fields `data`, `theta` and `latent` are lists of
#'   that length.
#'
#' @seealso [statmod()] to fit what this draws,
#'   [simulate.StatmodFit()] to draw from a model already fitted,
#'   [print.StatmodSim()] for the printed form.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(50), g = factor(rep(c("a", "b"), 25)))
#'
#' # coefficients drawn
#' sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd)
#' sim$par
#'
#' # or given, and recovered by a fit
#' sim2 <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
#'                  par = list(mu = c(1, 2), sigma = log(0.3)))
#' coef(statmod(y ~ x, distributions7::gaussian1_distrib(), sim2$data))
#'
#' # a sparse truth, written as a function of the coefficient count
#' dd2 <- as.data.frame(matrix(rnorm(50 * 6), 50, 6))
#' sim3 <- rstatmod(y ~ lasso(~ V1 + V2 + V3 + V4 + V5 + V6),
#'                  distributions7::gaussian1_distrib(), dd2,
#'                  par = list(mu = function(k) c(2, -1.5, rep(0, k - 2)),
#'                             sigma = log(0.3)))
#' head(sim3$data$y, 3)
#'
#' # a model with no covariates at all
#' sim4 <- rstatmod(y ~ 1, distributions7::gaussian1_distrib(), n = 20)
#' nrow(sim4$data)
#'
#' # a study over replicates, the covariates drawn afresh at each one
#' study <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), n = 80,
#'                   n_sim = 5, par = list(mu = c(1, 2), sigma = log(0.5)),
#'                   covariates = list(x = function(n) runif(n, -2, 2)))
#' length(study$data)
#' vapply(study$data, function(d) coef(statmod(
#'   y ~ x, distributions7::gaussian1_distrib(), d),
#'   readable = FALSE)$mu[[2L]], numeric(1))
#'
#' # a score-driven series, its own parameters named
#' sim5 <- rstatmod(y ~ 0 + gas(p = 1, q = 1, time = t),
#'                  distributions7::gaussian1_distrib(),
#'                  data.frame(t = 1:100), par = list(sigma = 0),
#'                  structural = list(omega = 0.4, alpha1 = 0.3,
#'                                    pacf1 = 0.6))
#' head(sim5$latent, 3)
#'
#' @export
rstatmod <- function(formula, distrib, data = NULL, n = NULL, n_sim = 1,
                     par = NULL, structural = NULL, sd = 1, offsets = NULL,
                     covariates = NULL) {
  if (!S7::S7_inherits(distrib, distributions7::distrib)) {
    stop("'distrib' must be a distributions7 distribution object.",
         call. = FALSE)
  }
  n_sim <- as.integer(n_sim)
  if (length(n_sim) != 1L || is.na(n_sim) || n_sim < 1L) {
    stop("'n_sim' must be a positive whole number.", call. = FALSE)
  }
  covariates <- check_covariates(covariates)
  data <- rstatmod_data(data, n, covariates)
  n <- nrow(data)
  params <- distrib@params
  split <- statmod_equations(formula, params)
  nm <- rstatmod_response_name(split$response)

  build <- function(dat) {
    # THE SPECIFICATION IS THE REAL ONE, built against a placeholder
    # response, so that the design, the offsets, the interpreter and every
    # term's blueprint are the ones a fit would produce. What the
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
  coef <- draw_coefficients(b$design, params, par, sd)
  names_at <- function(d) lapply(d[params], `[[`, "coef_names")
  cn <- names_at(b$design)

  dats <- vector("list", n_sim)
  thetas <- vector("list", n_sim)
  latents <- vector("list", n_sim)
  psi <- NULL
  for (r in seq_len(n_sim)) {
    if (r > 1L && !is.null(covariates)) {
      data <- draw_covariates(covariates, n, rstatmod_data(NULL, n, NULL))
      b <- build(data)
      # a design that changed shape cannot take coefficients drawn against
      # the first one, and recycling them would fit a different model
      if (!identical(names_at(b$design), cn)) {
        stop("the simulated covariates gave a design of another shape at ",
             "replicate ", r, ".\n  A factor that lost a level is the ",
             "ordinary cause; draw the covariates so that\n  every ",
             "replicate spans the same columns.", call. = FALSE)
      }
    }
    ep <- rstatmod_eta(b$spec, b$design, coef, structural)
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
#' function is rejected rather than recycled: a constant column is what
#' `data` is for, and the whole point of this argument is that the
#' covariates are drawn afresh at every replicate.
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
#' A generator that answers with the wrong length is reported rather than
#' recycled, for the reason a coefficient function is: R would recycle it
#' without a word and the replicate would be of another model. A column of
#' `data` under the same name is overwritten, as a caller
#' asking for that column to be drawn means.
#'
#' @param covariates The generators.
#' @param n The observation count.
#' @param data The frame to write into.
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
#' The model that was written, the size of what came out, and the truth
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
    cat(sprintf("  %-10s %s\n", p,
                paste(format(round(v, 4)), collapse = "  ")))
  }
  if (length(x$structural)) {
    cat("\n  the term's own parameters\n")
    for (q in names(x$structural)) {
      cat(sprintf("  %-10s %s\n", q, format(round(x$structural[[q]], 4))))
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


#' The Data Frame a Simulation Runs Against
#'
#' @description
#' Resolves `data` and `n` into one data frame of covariates.
#'
#' @details
#' A model may have no covariates at all -- a pure time series is the case --
#' and then a row count is the whole of what a simulation needs, so the two
#' arguments are alternatives rather than one being compulsory. Where both
#' are given they must agree, which is checked rather than resolved by
#' preferring one: a caller who wrote both and got them wrong wants to know.
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
    n <- as.integer(n)
    if (length(n) != 1L || is.na(n) || n < 1L) {
      stop("'n' must be a positive whole number.", call. = FALSE)
    }
    out <- data.frame(row.names = seq_len(n))
    if (!is.null(covariates)) out <- draw_covariates(covariates, n, out)
    return(out)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!is.null(n) && !identical(as.integer(n), nrow(data))) {
    stop(sprintf("'n' is %d and 'data' has %d rows.",
                 as.integer(n), nrow(data)), call. = FALSE)
  }
  if (!is.null(covariates)) {
    data <- draw_covariates(covariates, nrow(data), data)
  }
  data
}


#' The Column a Simulated Response Is Written To
#'
#' @description
#' The name on the formula's left-hand side, which must be a symbol.
#'
#' @details
#' A transformed response is rejected rather than answered. Under
#' `log(y) ~ x` the model generates values of `log(y)`, so a column
#' called `y` would hold the wrong quantity and one called
#' `"log(y)"` would not be the name the formula reads back; the earlier
#' version wrote the first of those silently. A censored response is rejected
#' for the reason [statmod()] rejects one -- there is no censored
#' likelihood to fit it back with.
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
         "censored likelihood\n  to fit it back with, so the data would be ",
         "of a model it cannot read. Simulate\n  the uncensored response ",
         "and censor it yourself.", call. = FALSE)
  }
  stop(sprintf(paste0("the left-hand side must be a symbol, and it is '%s'.",
                      "\n  The model generates that quantity, so there is ",
                      "no column to write it to:\n  write the model in the ",
                      "variable you want simulated."),
               deparse(response)), call. = FALSE)
}


#' The Predictor and Parameters of a Simulation
#'
#' @description
#' The linear predictors, the distribution's parameters, and -- where a term
#' has state -- the response it drew and the latent quantity behind it.
#'
#' @details
#' With no structural term this is [statmod_eta()] exactly, which is
#' the point: the simulated data comes from the assembly a fit reads, so a
#' term whose block moves with its coefficients contributes what it
#' contributes rather than a linearization.
#'
#' With one, the static part is assembled the same way and the term is asked
#' to finish through [modelterms7::term_simulate()]. A term that
#' draws the response as it goes returns it; one that does not returns
#' `NULL` there and the caller draws at the predictor.
#'
#' @param spec The specification.
#' @param design Its design.
#' @param coef The coefficients.
#' @param structural The structural term's own parameters, or `NULL`.
#'
#' @return A list with `eta`, `theta`, `y`, `latent`
#'   and `structural`.
#'
#' @seealso [rstatmod()], [statmod_eta()]
#'
#' @keywords internal
rstatmod_eta <- function(spec, design, coef, structural = NULL) {
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  su <- attr(design, "structural")
  if (!length(su)) {
    if (!is.null(structural)) {
      stop("'structural' names a term's own parameters and this formula ",
           "carries no\n  structural term.", call. = FALSE)
    }
    ep <- statmod_eta(spec, design, coef)
    return(list(eta = ep$eta, theta = ep$theta, y = NULL, latent = NULL,
                structural = NULL))
  }
  if (length(su) > 1L) {
    stop("a formula carries at most one structural term.", call. = FALSE)
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

  st <- statmod_structural_state(design)
  tm <- spec@terms[[u$param]][[u$term]]
  psi <- rstatmod_psi(tm, structural_psi(tm, st$zeta[[u$term]]), structural)
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
  res <- modelterms7::term_simulate(
    tm, psi, eta[[p]], draw,
    score = function(y, e, i) as.numeric(k$score(y, at(i), e)),
    curvature = function(y, e, i) as.numeric(k$curvature(y, at(i), e)))
  eta[[p]] <- res$eta
  theta[[p]] <- linkfunctions7::linkinv(lk, res$eta)
  list(eta = eta, theta = theta, y = res$y, latent = res$latent,
       structural = as.list(psi))
}


#' Draw or Validate the Coefficients of a Simulation
#'
#' @description
#' Returns one coefficient vector per distribution parameter, drawn from a
#' normal where the caller gave none.
#'
#' @details
#' An entry of `par` may be a vector of the equation's own length, a
#' single number used for all of them, or a function of the count returning
#' that many values. The function form is what expresses a structured truth
#' -- a random effect's standard deviation, a sparse vector for a lasso to
#' find -- without a vocabulary of its own, and it is checked to have
#' returned the right number of finite values, since a function that answers
#' wrongly would otherwise be recycled into a different model.
#'
#' @param design The design blocks.
#' @param params The parameter names.
#' @param par A named list, or `NULL`.
#' @param sd The standard deviation of the drawn coefficients.
#'
#' @return A named list of numeric vectors.
#'
#' @seealso [rstatmod()]
#'
#' @keywords internal
draw_coefficients <- function(design, params, par, sd) {
  if (!is.null(par)) {
    if (!is.list(par) || is.null(names(par))) {
      stop("'par' must be a named list, one entry per parameter.",
           call. = FALSE)
    }
    bad <- setdiff(names(par), params)
    if (length(bad)) {
      stop(sprintf(paste0("'par' names '%s', which is not a parameter.\n",
                          "  They are: %s."),
                   bad[1L], paste(params, collapse = ", ")), call. = FALSE)
    }
  }
  stats::setNames(lapply(params, function(p) {
    k <- design[[p]]$npar
    e <- if (is.null(par)) NULL else par[[p]]
    if (is.null(e)) return(stats::rnorm(k, 0, sd))
    if (is.function(e)) {
      v <- as.numeric(e(k))
      if (length(v) != k || anyNA(v)) {
        stop(sprintf(paste0("'par$%s' is a function and returned %d value%s ",
                            "where '%s' has %d\n  coefficient%s."),
                     p, length(v), if (length(v) == 1L) "" else "s", p, k,
                     if (k == 1L) "" else "s"), call. = FALSE)
      }
      return(v)
    }
    v <- as.numeric(e)
    if (length(v) == 1L && k > 1L) return(rep(v, k))
    if (length(v) != k) {
      stop(sprintf("'par$%s' has length %d but '%s' has %d coefficients.",
                   p, length(v), p, k), call. = FALSE)
    }
    v
  }), params)
}


#' A Structural Term's Parameters for a Simulation
#'
#' @description
#' The term's starting values with whatever the caller named written over
#' them.
#'
#' @details
#' The values are on the scale [modelterms7::term_params()] names,
#' which is the one a reader knows: a loading is the loading rather than its
#' logarithm, and a persistence is the partial autocorrelation its chart
#' carries rather than the autoregressive coefficient that chart produces. A
#' name the term does not have is reported with the ones it does, since a
#' misspelled parameter would otherwise leave the term at a starting value and
#' simulate a model with almost no dynamics.
#'
#' @param tm The built structural term.
#' @param psi Its parameters at the specification's own state.
#' @param given A named list, or `NULL`.
#'
#' @return A named list.
#'
#' @seealso [rstatmod()]
#'
#' @keywords internal
rstatmod_psi <- function(tm, psi, given) {
  if (is.null(given)) return(psi)
  if (!is.list(given) || is.null(names(given))) {
    stop("'structural' must be a named list.", call. = FALSE)
  }
  nm <- modelterms7::term_params(tm)
  bad <- setdiff(names(given), nm)
  if (length(bad)) {
    stop(sprintf(paste0("'structural' names '%s', which is not a parameter ",
                        "of this term.\n  They are: %s."),
                 bad[[1L]], paste(nm, collapse = ", ")), call. = FALSE)
  }
  for (q in names(given)) {
    v <- as.numeric(given[[q]])
    if (length(v) != 1L || !is.finite(v)) {
      stop(sprintf("'structural$%s' must be one finite number.", q),
           call. = FALSE)
    }
    psi[[q]] <- v
  }
  psi
}
