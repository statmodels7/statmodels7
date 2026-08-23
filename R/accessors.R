#' @importFrom stats nobs formula family weights df.residual sigma
#' @importFrom stats model.matrix simulate terms model.frame
#' @importFrom stats anova
#' @include statmod.R
NULL

#' @title The Number of Observations a Model Was Fitted To
#' @name nobs.StatmodFit
#' @description The row count of the fitting data.
#' @param object A [StatmodFit()].
#' @param ... Unused.
#' @return An integer.
#' @seealso [logLik.StatmodFit()], [statmod()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' nobs(fit)
#' @keywords internal
nobs.StatmodFit <- function(object, ...) object@spec@n_obs
S7::method(nobs, StatmodFit) <- nobs.StatmodFit


#' @title The Formula a Model Was Written With
#' @name formula.StatmodFit
#' @description
#' The formula as supplied, every distribution parameter's equation included.
#' @details
#' It is returned whole rather than split into one formula per parameter: the
#' bars are part of what was written, and a caller wanting the equations
#' separately gets them from the fit's specification.
#' @param x A [StatmodFit()].
#' @param ... Unused.
#' @return A formula.
#' @seealso [statmod()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#' formula(fit)
#' @keywords internal
formula.StatmodFit <- function(x, ...) x@spec@formula
S7::method(formula, StatmodFit) <- formula.StatmodFit


#' @title The Distribution a Model Was Fitted With
#' @name family.StatmodFit
#' @description The \pkg{distributions7} object, with its links.
#' @details
#' It is the family itself rather than a description of one, so everything
#' the family can do is available from a fit: its density, its derivatives,
#' its moments and its parameters' links.
#' @param object A [StatmodFit()].
#' @param ... Unused.
#' @return A \pkg{distributions7} distribution.
#' @seealso [statmod()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' family(fit)@params
#' @keywords internal
family.StatmodFit <- function(object, ...) object@spec@distrib
S7::method(family, StatmodFit) <- family.StatmodFit


#' @title The Prior Weights of a Fitted Model
#' @name weights.StatmodFit
#' @description The weights each observation entered the likelihood with.
#' @param object A [StatmodFit()].
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [statmod()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' head(weights(fit))
#' @keywords internal
weights.StatmodFit <- function(object, ...) object@spec@weights
S7::method(weights, StatmodFit) <- weights.StatmodFit


#' @title What a Fit Leaves Unspent
#' @name df.residual.StatmodFit
#' @description
#' The observation count less the effective degrees of freedom the model
#' spent.
#' @details
#' The count subtracted is the EFFECTIVE one, the trace of the model's
#' smoother, and not the number of coefficients: a penalized block spends
#' less than it carries, which is the whole reason a smoothing parameter is
#' estimated. It is therefore not an integer, and for a model whose degrees
#' of freedom could not be counted it is `NA`.
#' @param object A [StatmodFit()].
#' @param ... Unused.
#' @return A number.
#' @seealso [logLik.StatmodFit()], [statmod_edf()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- sin(3 * dd$x) + rnorm(60, sd = 0.3)
#' fit <- statmod(y ~ s(x, k = 6), distributions7::gaussian1_distrib(), dd)
#' df.residual(fit)
#' @keywords internal
df.residual.StatmodFit <- function(object, ...) {
  ll <- stats::logLik(object)
  as.numeric(object@spec@n_obs) - as.numeric(attr(ll, "df"))
}
S7::method(df.residual, StatmodFit) <- df.residual.StatmodFit


#' @title The Standard Deviation a Fitted Model Implies
#' @name sigma.StatmodFit
#' @description
#' The fitted standard deviation at each observation, where the family has
#' one.
#' @details
#' A VECTOR rather than a number, because the whole point of the framework is
#' that a scale may be modelled: a single residual standard deviation exists
#' only where the scale's equation is an intercept, and returning its first
#' value would silently answer a different question everywhere else.
#'
#' What is returned is the standard deviation of the response under the
#' fitted distribution, through [distributions7::std_dev()], and not
#' whichever parameter happens to be spelled `sigma`: for a Gamma
#' written by its mean and dispersion the two are different quantities. A
#' family with no second moment signals an error rather than reporting one.
#' @param object A [StatmodFit()].
#' @param ... Unused.
#' @return A numeric vector.
#' @seealso [predict.StatmodFit()], [fitted.StatmodFit()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#' head(sigma(fit))
#' @keywords internal
sigma.StatmodFit <- function(object, ...) {
  stats::predict(object, "std_dev")
}
S7::method(sigma, StatmodFit) <- sigma.StatmodFit


#' @title The Design of One Equation
#' @name model.matrix.StatmodFit
#' @description
#' The block of columns a distribution parameter's equation was fitted with.
#' @details
#' A fit has one design per parameter, so which one is asked for is an
#' argument rather than something to be guessed; `NULL` gives the first,
#' as [fitted.StatmodFit()] does. A term whose block moves with its
#' coefficients is returned AT the fitted ones, which is the block the fit
#' ended on.
#'
#' A structural term contributes no columns at all, so an equation carrying
#' one alone gives a matrix of no columns; what such a term contributes is a
#' recursion, reported by [predict.StatmodFit()].
#' @param object A [StatmodFit()].
#' @param what Which distribution parameter, or `NULL` for the first.
#' @param ... Unused.
#' @return A matrix, sparse where the equation's blocks are.
#' @seealso [statmod_design()], [coef.StatmodFit()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#' dim(model.matrix(fit))
#' colnames(model.matrix(fit, "sigma"))
#' @keywords internal
model.matrix.StatmodFit <- function(object, what = NULL, ...) {
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
#' Responses drawn from the fitted distribution, one column per replicate.
#' @details
#' The draws are taken at the FITTED parameters, so what they carry is the
#' variation of the response and not the uncertainty of the estimates; a
#' parametric bootstrap adds the second by refitting each column.
#'
#' For a model carrying a structural term the fitted parameters are the ones
#' the filter reached ALONG the observed series, so the draws are conditional
#' on that series rather than a fresh path of the process.
#' @param object A [StatmodFit()].
#' @param nsim How many replicates.
#' @param seed Passed to [base::set.seed()] if given, the caller's
#'   stream being restored afterwards.
#' @param ... Unused.
#' @return A data frame of `nsim` columns.
#' @seealso [rstatmod()], [predict.StatmodFit()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' dim(simulate(fit, nsim = 3))
#' @keywords internal
simulate.StatmodFit <- function(object, nsim = 1, seed = NULL, ...) {
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L) {
    stop("'nsim' must be a positive whole number.", call. = FALSE)
  }
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


#' @title What a Fit Does Not Answer
#' @name statmod_refusals
#' @description
#' Three generics of \pkg{stats} signal an error on a statmod fit, each
#' naming what to ask instead.
#' @details
#' `terms()` would have to report one set of terms where a fit has one
#' per distribution parameter, and the formula it was written with is not a
#' `terms` object -- the bars separating the equations are not
#' `stats`' syntax. `formula()` gives what was written and
#' [statmod_design()] gives what it produced.
#'
#' `model.frame()` would have to return the fitting data, which a fit
#' does not keep: what it keeps is each term's blueprint, so that new data is
#' reapplied rather than relearned.
#'
#' `anova()` would have to compare models by a test, and a penalized fit
#' whose hyperparameters were chosen from the same data has no null
#' distribution to compare against. [logLik.StatmodFit()],
#' `AIC` and `BIC` are what this package reports, with the
#' effective degrees of freedom corrected for the smoothing parameters
#' having been estimated.
#' @param object,x,formula A [StatmodFit()].
#' @param ... Unused.
#' @return Nothing; each method signals an error.
#' @seealso [formula.StatmodFit()],
#'   [model.matrix.StatmodFit()], [logLik.StatmodFit()]
#' @keywords internal
NULL

#' @rdname statmod_refusals
#' @keywords internal
terms.StatmodFit <- function(x, ...) {
  stop("a statmod fit has one set of terms per distribution parameter, not ",
       "one.\n  formula(fit) gives what was written and ",
       "statmod_design(fit@spec) what it produced.", call. = FALSE)
}
S7::method(terms, StatmodFit) <- terms.StatmodFit

#' @rdname statmod_refusals
#' @keywords internal
model.frame.StatmodFit <- function(formula, ...) {
  stop("a statmod fit does not keep the fitting data: what it keeps is each ",
       "term's\n  blueprint, so that new data is reapplied rather than ",
       "relearned. model.matrix(fit)\n  gives an equation's design.",
       call. = FALSE)
}
S7::method(model.frame, StatmodFit) <- model.frame.StatmodFit

#' @rdname statmod_refusals
#' @keywords internal
anova.StatmodFit <- function(object, ...) {
  stop("a penalized fit whose hyperparameters were chosen from the same ",
       "data has no\n  null distribution to test against. logLik(), AIC() ",
       "and BIC() report the fit,\n  with the degrees of freedom corrected ",
       "for the smoothing parameters having been\n  estimated.",
       call. = FALSE)
}
S7::method(anova, StatmodFit) <- anova.StatmodFit
