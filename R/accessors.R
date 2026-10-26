#' @importFrom stats nobs formula family weights df.residual sigma
#' @importFrom stats model.matrix simulate terms model.frame
#' @importFrom stats anova
#' @include statmod.R
NULL

#' @title The Number of Observations in a Fitted Model
#' @name nobs.StatmodFit
#' @description Returns the number of observations that entered the
#'   likelihood. This number is the sum of the prior weights, which are
#'   frequency weights.
#' @details A weight of two counts its row twice, so a fit with weights is
#'   counted as the fit of the data with every row repeated as many times as
#'   its weight. Without weights the result is the number of rows of the
#'   data, since [statmod()] drops no rows (a missing value in the response or
#'   in a variable of an equation is an error). [bic()], the cBIC of
#'   [summary.StatmodFit()], [logLik.StatmodFit()]'s `"nobs"` attribute and
#'   [df.residual.StatmodFit()] read the same number.
#' @param object A [StatmodFit()].
#' @param ... Unused.
#' @return A single number, the sum of the weights: an integer where that sum
#'   is a whole number, and a double otherwise.
#' @seealso [df.residual.StatmodFit()], which subtracts the effective degrees
#'   of freedom from it, [statmod()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' nobs(fit)
#' @keywords internal
nobs.StatmodFit <- function(object, ...) {
  check_dots("nobs", ...)
  n <- weight_total(object@spec)
  # an integer where the total is whole, as the row count was
  if (n == round(n) && n <= .Machine$integer.max) as.integer(n) else n
}
S7::method(nobs, StatmodFit) <- nobs.StatmodFit


#' The Number of Observations Counted by the Weights
#'
#' @description
#' Returns the sum of the prior weights, which are frequency weights: a weight
#' of two counts the row twice. Without weights it is the number of rows.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A single number.
#'
#' @seealso [nobs.StatmodFit()], [bic()]
#'
#' @keywords internal
weight_total <- function(spec) sum(spec@weights)


#' @title The Formula of a Fitted Model
#' @name formula.StatmodFit
#' @description
#' Returns the formula as supplied, including the equation of every
#' distribution parameter.
#' @details
#' The formula is returned whole, including its bars, because the bars are part
#' of what was written. [statmod_equations()] splits the formula into one
#' equation per distribution parameter; it takes the formula and the names of
#' the parameters.
#' @param x A [StatmodFit()].
#' @param ... Unused.
#' @return The formula object supplied to [statmod()], including its bars and
#'   its original environment.
#' @seealso [statmod_equations()] to split it into one equation per
#'   parameter, [statmod()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#' formula(fit)
#' @keywords internal
formula.StatmodFit <- function(x, ...) {
  check_dots("formula", ...)
  x@spec@formula
}
S7::method(formula, StatmodFit) <- formula.StatmodFit


#' @title The Distribution of a Fitted Model
#' @name family.StatmodFit
#' @description Returns the \pkg{distributions7} object of the model, with its
#'   links.
#' @details
#' The object returned is the family itself, so its density, its derivatives,
#' its moments and its parameters' links are available from a fit. In
#' contrast, [stats::glm()] has a `family()` that returns a description.
#' @param object A [StatmodFit()].
#' @param ... Unused.
#' @return The \pkg{distributions7} distribution object used to fit the
#'   model, including its links. A binomial family fitted to a two-column
#'   response carries the row sums as its `size`, so in that case the object
#'   differs from the one supplied. Its `@params` property names the parameters
#'   and `@link_params` holds their links.
#' @seealso [statmod()]; [distributions7::distrib_pdf()] and its siblings for
#'   the quantities that can be computed from the family
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' family(fit)@params
#' @keywords internal
family.StatmodFit <- function(object, ...) {
  check_dots("family", ...)
  object@spec@distrib
}
S7::method(family, StatmodFit) <- family.StatmodFit


#' @title The Prior Weights of a Fitted Model
#' @name weights.StatmodFit
#' @description Returns the weight with which each observation entered the
#'   likelihood.
#' @details
#' The weights multiply the log-likelihood of each observation, so they act as
#' frequency weights: a constant weight \eqn{c} on every row is \eqn{c} copies
#' of the data. The restricted likelihood then counts \eqn{n c} observations
#' against the \eqn{p} coefficients, so the REML estimate of a dispersion
#' changes under a constant rescaling of the weights, while the ML estimate
#' does not.
#' @param object A [StatmodFit()].
#' @param ... Unused.
#' @return A numeric vector with one entry per row of the data, the prior
#'   weights as supplied. The vector is all ones if no weights were given.
#'   The weights enter as given and are never normalized, so a weight of two
#'   counts an observation twice.
#' @seealso [statmod()], [nobs.StatmodFit()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' head(weights(fit))
#'
#' # The weights are frequencies, so nobs() is their sum.
#' sum(weights(fit)) == nobs(fit)
#' @keywords internal
weights.StatmodFit <- function(object, ...) {
  check_dots("weights", ...)
  object@spec@weights
}
S7::method(weights, StatmodFit) <- weights.StatmodFit


#' @title Residual Degrees of Freedom of a Fitted Model
#' @name df.residual.StatmodFit
#' @description
#' Returns the number of observations, [nobs.StatmodFit()], less the
#' effective degrees of freedom of the model.
#' @details
#' The count subtracted is the effective one, the sum over the terms of their
#' effective degrees of freedom, as [logLik.StatmodFit()] reports it. A
#' penalized block has fewer effective degrees of freedom than the number of
#' coefficients that it carries, so the result is not an integer in general.
#' @param object A [StatmodFit()].
#' @param ... Unused.
#' @return A single number, `nobs(object)` less the total effective degrees
#'   of freedom. The result is not an integer in general. Where the effective
#'   degrees of freedom of a term could not be computed, [logLik.StatmodFit()]
#'   uses the number of coefficients of that term (the number of free
#'   parameters, for a structural term) in its place, so the result is always a
#'   number.
#' @seealso [statmod_edf()] for the per-term counts that this function sums,
#'   [logLik.StatmodFit()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- sin(3 * dd$x) + rnorm(60, sd = 0.3)
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 6)), distributions7::gaussian1_distrib(), dd)
#' df.residual(fit)
#' @keywords internal
df.residual.StatmodFit <- function(object, ...) {
  check_dots("df.residual", ...)
  ll <- stats::logLik(object)
  weight_total(object@spec) - as.numeric(attr(ll, "df"))
}
S7::method(df.residual, StatmodFit) <- df.residual.StatmodFit


#' @title The Standard Deviation Implied by a Fitted Model
#' @name sigma.StatmodFit
#' @description
#' Returns the fitted standard deviation at each observation, where the family
#' has one.
#' @details
#' The result is a vector with one entry per observation, because a scale may
#' be modeled. A single residual standard deviation exists only where the
#' equation of the scale is an intercept, and the first value of the vector
#' would be a different quantity in any other model.
#'
#' The value is the standard deviation of the response under the fitted
#' distribution, computed through [distributions7::std_dev()] from the fitted
#' parameters. It is not the parameter that happens to be named `sigma`: for a
#' Gamma written by its mean and dispersion, the two are different quantities.
#' Where the second moment does not exist at the fitted parameters (the Cauchy
#' distribution, for example), the entries are `NaN`.
#' @param object A [StatmodFit()].
#' @param ... Unused.
#' @return A numeric vector of length `nobs(object)`, the standard deviation
#'   of the response under the fitted distribution at each observation. It is
#'   constant when no equation models a scale and the variance of the family
#'   does not depend on the location; for a Poisson family, for example, it
#'   varies with the fitted mean.
#' @seealso [fitted.StatmodFit()] for the fitted parameters themselves,
#'   [distributions7::std_dev()] for the quantity
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#' head(sigma(fit))
#' @keywords internal
sigma.StatmodFit <- function(object, ...) {
  check_dots("sigma", ...)
  stats::predict(object, "std_dev")
}
S7::method(sigma, StatmodFit) <- sigma.StatmodFit


#' @title The Design of One Equation
#' @name model.matrix.StatmodFit
#' @description
#' Returns the block of columns used to fit the equation of a distribution
#' parameter.
#' @details
#' A fit has one design per parameter, and `what` selects the parameter.
#' `NULL` selects the first, as [fitted.StatmodFit()] does. A term whose block
#' moves with its coefficients is returned at the fitted coefficients, which
#' is the block on which the fit ended.
#'
#' A structural term contributes no columns at all, so an equation that
#' carries only a structural term gives a matrix of no columns. The
#' contribution of such a term is a recursion, reported by
#' [predict.StatmodFit()].
#' @param object A [StatmodFit()].
#' @param what Which distribution parameter, or `NULL` for the first.
#' @param ... Unused.
#' @return An `n x p` matrix, where `p` is the width of the design of that
#'   equation, with the coefficient names as its column names. The result is a
#'   \pkg{Matrix} object where the blocks of the equation are sparse and a
#'   base matrix otherwise. The matrix has zero columns for an equation that
#'   carries only a structural term.
#' @seealso [statmod_design()] for every equation at once,
#'   [coef.StatmodFit()] for the coefficients that these columns multiply
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#' dim(model.matrix(fit))
#' colnames(model.matrix(fit, "sigma"))
#' @keywords internal
model.matrix.StatmodFit <- function(object, what = NULL, ...) {
  check_dots("model.matrix", ...)
  spec <- object@spec
  params <- spec@distrib@params
  p <- if (is.null(what)) params[[1L]] else what
  if (!is.character(p) || length(p) != 1L || !p %in% params) {
    stop(sprintf("'what' must name one of: %s.",
                 paste(params, collapse = ", ")), call. = FALSE)
  }
  d <- statmod_design_at(spec, object@coefficients,
                         statmod_design(spec))[[p]]
  if (!d$npar) return(matrix(0, spec@n_obs, 0L))
  X <- d$X
  colnames(X) <- d$coef_names
  X
}
S7::method(model.matrix, StatmodFit) <- model.matrix.StatmodFit


#' @title Draws from a Fitted Model
#' @name simulate.StatmodFit
#' @description
#' Draws responses from the fitted distribution, one column per replicate.
#' @details
#' The draws are taken at the fitted parameters, so they carry only the
#' variation of the response. The uncertainty of the estimates is not
#' in them; a parametric bootstrap adds it by refitting each column.
#'
#' For a model carrying a structural term the fitted parameters are the ones
#' that the filter reached along the observed series, so the draws are
#' conditional on that series and are not a fresh path of the process.
#' @param object A [StatmodFit()].
#' @param nsim How many replicates.
#' @param seed Passed to [base::set.seed()] if given; the stream of the caller
#'   is restored afterwards.
#' @param ... Unused.
#' @return A data frame with `nobs(object)` rows and `nsim` columns, where
#'   each column is one draw of the whole response vector.
#' @seealso [rstatmod()] to simulate from a model that was never fitted,
#'   [predict.StatmodFit()] for the fitted parameters from which the draws
#'   come
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' dim(simulate(fit, nsim = 3))
#' @keywords internal
simulate.StatmodFit <- function(object, nsim = 1, seed = NULL, ...) {
  check_dots("simulate", ...)
  if (!is_count(nsim)) {
    stop("'nsim' must be a positive whole number.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1)
    }
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  spec <- object@spec
  theta <- object@fitted
  n <- spec@n_obs
  out <- vector("list", nsim)
  for (k in seq_len(nsim)) {
    out[[k]] <- distributions7::distrib_rng(spec@distrib, n, theta)
  }
  names(out) <- paste0("sim_", seq_len(nsim))
  as.data.frame(out)
}
S7::method(simulate, StatmodFit) <- simulate.StatmodFit


#' @title Generics Not Available for a Fit
#' @name statmod_unavailable
#' @description
#' Three generics of \pkg{stats} and the generic [base::plot()] signal an
#' error on a statmod fit, and each error names the function to use instead.
#' @details
#' `terms()` has no single answer where a fit has one set of terms per
#' distribution parameter, and the formula of a fit is not a `terms` object,
#' because the bars separating the equations are not `stats` syntax.
#' [formula.StatmodFit()] gives what was written and [statmod_design()] gives
#' what it produced.
#'
#' `model.frame()` would have to return the data frame used to fit the model,
#' which a fit does not keep. A fit keeps the response, the weights, the
#' offsets and each term's blueprint, so that a basis, a set of levels or a
#' set of contrasts is reapplied to new data and not learned from it again.
#'
#' `anova()` would have to compare models by a test, and a penalized fit
#' whose hyperparameters were chosen from the same data has no null
#' distribution against which to compare. [logLik.StatmodFit()], `AIC()` and
#' `BIC()` report the fit with the sum of the effective degrees of freedom of
#' its terms, which is not corrected for the estimation of the smoothing
#' parameters; `summary(fit, correct = TRUE)` applies that correction.
#'
#' `plot()` has no single picture to draw for a model with an equation
#' per parameter, any number of terms and any number of covariates. The
#' fitted curves come from [predict.StatmodFit()] on a grid of new data, and
#' the residuals from [residuals.StatmodFit()].
#' @param object,x,formula A [StatmodFit()].
#' @param ... Unused.
#' @return Nothing. Each method signals an error that names the function to
#'   use instead.
#' @seealso [formula.StatmodFit()],
#'   [model.matrix.StatmodFit()], [logLik.StatmodFit()]
#' @keywords internal
NULL

#' @rdname statmod_unavailable
#' @keywords internal
terms.StatmodFit <- function(x, ...) {
  stop("a statmod fit has one set of terms per distribution parameter, not ",
       "one.\n  formula(fit) gives what was written and ",
       "statmod_design(fit@spec) what it produced.", call. = FALSE)
}
S7::method(terms, StatmodFit) <- terms.StatmodFit

#' @rdname statmod_unavailable
#' @keywords internal
model.frame.StatmodFit <- function(formula, ...) {
  stop("a statmod fit does not keep the fitting data; it keeps each ",
       "term's\n  blueprint, so that new data is reapplied rather than ",
       "relearned. model.matrix(fit)\n  gives an equation's design.",
       call. = FALSE)
}
S7::method(model.frame, StatmodFit) <- model.frame.StatmodFit

#' @rdname statmod_unavailable
#' @keywords internal
anova.StatmodFit <- function(object, ...) {
  stop("a penalized fit whose hyperparameters were chosen from the same ",
       "data has no\n  null distribution against which to test. logLik(), AIC() ",
       "and BIC() report the fit with the\n  uncorrected effective degrees ",
       "of freedom; summary(fit, correct = TRUE) corrects them.",
       call. = FALSE)
}
S7::method(anova, StatmodFit) <- anova.StatmodFit

#' @rdname statmod_unavailable
#' @keywords internal
plot.StatmodFit <- function(x, ...) {
  stop("a statmod fit has no plot() method: a model with an equation per ",
       "parameter has
  no single picture. predict(fit, newdata = grid) ",
       "gives the fitted curves and
  residuals(fit) the quantile residuals ",
       "to draw.", call. = FALSE)
}
S7::method(plot, StatmodFit) <- plot.StatmodFit


#' Reject an Argument That a Method Does Not Take
#'
#' @description
#' Signals an error when a method that has `...` for compatibility with its
#' generic receives an argument in it. The message names each argument, as in
#' `vcov() takes no argument 'x'.`, and shows `<unnamed>` for an argument
#' without a name. The check keeps a misspelled or nonexistent argument from
#' being accepted and ignored.
#'
#' @param fn The name of the method's generic, for the message.
#' @param ... The dots of the method.
#'
#' @return `NULL`, invisibly, when the dots are empty.
#'
#' @keywords internal
check_dots <- function(fn, ...) {
  if (!...length()) return(invisible(NULL))
  dn <- ...names()
  bad <- if (is.null(dn)) rep("", ...length()) else dn
  stop(sprintf("%s() takes no argument %s.", fn,
               paste0("'", ifelse(nzchar(bad), bad, "<unnamed>"), "'",
                      collapse = ", ")), call. = FALSE)
}
