#' @include report.R
#' @importFrom stats predict fitted
NULL

#' The Quantities a Fit Can Predict
#'
#' @description
#' The table [predict.StatmodFit()] resolves a moment name against: the five
#' names it understands, each mapped to the \pkg{distributions7} generic that
#' computes it.
#'
#' @details
#' Written once as a table so that the recognized names and the functions
#' they call cannot disagree, and so that [unknown_what()] can list them in
#' its message.
#'
#' @return A named list of five functions, keyed `"mean"`, `"variance"`,
#'   `"std_dev"`, `"skewness"` and `"kurtosis"`. Each takes a distribution
#'   and a parameter list and returns one value per observation.
#'
#' @seealso [predict.StatmodFit()], the only caller.
#'
#' @keywords internal
predict_moments <- function() {
  list(
    mean = function(d, th) mean(d, th),
    variance = function(d, th) distributions7::variance(d, th),
    std_dev = function(d, th) distributions7::std_dev(d, th),
    skewness = function(d, th) distributions7::skewness(d, th),
    kurtosis = function(d, th) distributions7::kurtosis(d, th)
  )
}


#' @title Predict From a Fitted Model
#' @name predict.StatmodFit
#' @description
#' Predicts from a fitted model: any one of the distribution's parameters,
#' any of its moments, or every parameter or linear predictor at once, at the
#' fitting data or at new data, with standard errors and intervals on
#' request.
#' @details
#' # What can be asked for
#'
#' `what` takes
#' \describe{
#'   \item{a parameter's name}{`"mu"`, `"sigma"`, `"alpha"` --
#'     whatever the family calls them. Always available, whatever the family:
#'     a parameter is what the model fits, and it exists even where a moment
#'     does not.}
#'   \item{a moment's name}{`"mean"`, `"variance"`, `"std_dev"`,
#'     `"skewness"`, `"kurtosis"`. Available where the family has one, and
#'     answering `NaN` or `NA` where it does not exist. A Cauchy's mean is
#'     `NaN`, which is the correct answer for it.}
#'   \item{`"parameter"`}{every parameter at once, as a named list. The
#'     default.}
#'   \item{`"link"`}{every linear predictor at once, before the inverse
#'     link.}
#'   \item{`"response"`}{the response itself, with `interval =
#'     "prediction"`: its median and the ends of an interval that a new
#'     observation falls in. See the section on intervals.}
#' }
#' A parameter's name may be prefixed by `"link:"` to ask for its
#' predictor instead of its value, as `"link:sigma"`.
#'
#' # The argument order departs from [stats::predict()]
#'
#' There the second argument is `newdata`; here it is `what`. A statmod fit
#' has several parameters and several moments, so choosing among them is the
#' ordinary variation and predicting on new data is the occasional one.
#'
#' Passing a data frame second is caught and named, instead of failing
#' somewhere inside.
#'
#' # New data
#'
#' Goes through each term's blueprint, so a factor keeps the levels and
#' contrasts it was fitted with, a spline its knots, and a basis its
#' reparametrization. Nothing is rebuilt from whatever the new frame happens
#' to contain.
#'
#' # A score-driven term is predicted past the series
#'
#' Such a term's contribution at one row is the state a recursion has
#' reached, so new rows without the response continue the series instead of
#' being read on their own. Each row is placed by its own time within its own
#' group, and must come after every observed time of that group; a row
#' falling inside the observed series is rejected, since there the response
#' is known and the filter must be run, never continued.
#'
#' New rows that carry the response are a series of their own: the
#' structural term is rebuilt on them with the fitted parameters, and they
#' are read as the fitting rows are, standard errors included. With the
#' fitting data as `newdata` the result is `predict(fit)`'s; this holds for a
#' [modelterms7::regime()] term as well, whose prediction is the
#' posterior-weighted predictor.
#'
#' Beyond the data the score sits at its conditional mean of zero, which the
#' model's own definition guarantees, so the continuation is the
#' deterministic recursion and involves no simulation.
#'
#' # A new group, and the typical one
#'
#' A random-effect term predicts a group the fit saw with that group's
#' estimated effect. A group it never saw has none, and `random` says what to
#' put in its place. `"zero"` sets every effect of the term to zero, at every
#' row, which is the prediction for the typical group (lme4's and
#' glmmTMB's `re.form = NA`). `"marginal"` averages over the prior instead,
#' \eqn{E_b[h^{-1}(\eta_0 + z^\top b)]}, which is the population average.
#' With an identity link the two agree for the mean; with a log link the
#' marginal mean is \eqn{\exp(\eta_0 + \sigma_b^2/2)}, and with a logit it is
#' pulled towards one half. A moment is averaged by the law of total
#' expectation and the variance adds the variance of the conditional mean,
#' so a random effect on the scale enters the marginal variance of the
#' response.
#'
#' The average is a Gauss-Hermite product grid where every prior is Gaussian
#' and there are at most three coordinates in all; adaptive quadrature
#' against the prior's density ([numericals7::quad_vec()]) where the only
#' effect set aside is one coordinate with a prior that is not Gaussian; and
#' 10000 draws from the priors otherwise, on a seed of their own. A prior that
#' is not Gaussian is integrated only where the inverse link is bounded: a
#' Student t has no moment generating function, so under a log link the
#' average does not exist, and such a request signals an error. Under a link
#' whose predictor has a restricted domain (the square root and inverse
#' links) the average does not exist either, since a Gaussian effect leaves
#' the domain with positive probability; it is reported where that
#' probability is below 1e-8 and is `NA`, with a warning, elsewhere. A term that shares a covariance
#' with others through a label has no prior of its own and is read at
#' `"zero"` only. The standard error of a marginal parameter is the delta
#' method on the fixed part, conditional on the prior's scale.
#'
#' A random effect written inside a term's subformula, such as
#' `nl(~ Asym / (1 + exp((xmid - age) / scal)), Asym ~ 1 + random(~ 1 | g))`
#' or `seg(t, psi ~ random(~ 1 | id))`, is set aside by the same modes. Its
#' key is the outer term's key, the parameter and the sub-term joined by
#' `"::"`, and an unnamed mode applies to it as well. The effect enters the
#' term's parameter, and the term is not linear in it, so the term is
#' evaluated at the effect's value rather than having its columns removed:
#' `"zero"` evaluates it with the effect at zero at every row, and the
#' delta method reads its Jacobian there. `"marginal"` averages the term's
#' contribution, carried through the inverse link, over the nodes of the
#' product grid of a Gaussian prior; a prior of another family signals an
#' error. A predictor (`"link"`) is the average of the predictors at the
#' nodes. The standard error of a marginal parameter is the delta method of
#' that average, whose gradient sums the design at each node, and its
#' interval is the delta method's on the link scale carried through the
#' link. `interval = "group"` and `"prediction"` are not available for such
#' an effect, and neither is a standard error where effects of both kinds
#' are averaged over.
#'
#' # Three intervals
#'
#' `interval` says what an interval is for.
#' \describe{
#'   \item{`"confidence"`}{the default, with `se = TRUE`: the uncertainty of
#'     the estimates alone. For a term set aside by `random` it is the
#'     interval of the typical group's value, which a new group's value
#'     falls outside of far more often than the level says.}
#'   \item{`"group"`}{the interval of a new group's parameter, with `random =
#'     "zero"` or `"marginal"`. A new group's parameter is
#'     \eqn{\theta^* = h^{-1}(\eta^*)}, with \eqn{\eta^*} the typical
#'     group's predictor plus its estimation error plus \eqn{z^\top b}, the
#'     effects of a new group drawn from the prior, in the parameter's own
#'     equation, so a random effect on `sigma` widens the interval for
#'     `sigma`. Under Gaussian priors the variance of \eqn{\eta^*} is
#'     \eqn{se^2 + z^\top \Sigma_b z}, the standard error glmmTMB reports at
#'     a new level on the link scale. The interval's ends are quantiles of
#'     \eqn{\eta^*} carried through the inverse link, which is monotone, so
#'     they are the quantiles of \eqn{\theta^*} under any link; under a prior
#'     that is not Gaussian they are the quantiles of the mixture of
#'     [predictive_mixture()]. The standard error is the standard deviation
#'     of \eqn{\theta^*}, not the delta method, and is `NA` where it is not
#'     known to exist (see [group_interval()]). The fit is the typical
#'     group's value under `random = "zero"` and the population average
#'     under `"marginal"`. Implies `se = TRUE`.}
#'   \item{`"prediction"`}{an interval for a new observation of the response,
#'     with `what = "response"`. The family is averaged over the predictors of
#'     every parameter, taken jointly Gaussian -- correlated across equations
#'     where a label ties them -- and the interval's ends are quantiles of that
#'     average and the fit is its median, which every family has, so no
#'     location parameter is needed; `se` is the standard deviation of the
#'     average, `NA` where the family's variance does not exist. For a
#'     discrete family the ends are values of the support, so the coverage is
#'     at least the level rather than equal to it. Effects read
#'     `"conditional"` enter with their own estimates, so this is also the
#'     interval of a new observation of a group the fit saw. What the average
#'     is over is `predictive`, below.}
#' }
#'
#' # The predictive distribution
#'
#' A new observation varies for three reasons: the family itself, the effects
#' of a new group where `random` sets them aside, and the error in the
#' estimates. The first two are random in the model and always enter. The
#' third is what `predictive` decides, for every parameter alike, whether the
#' family has one parameter (a Poisson) or several, modelled or not.
#' \describe{
#'   \item{`"averaged"`}{the default. The family is averaged over the
#'     predictors at the fit's estimates with the covariance of the estimates
#'     from [vcov.StatmodFit()] plus the one the effects add. This is the
#'     normal approximation to the parametric bootstrap predictive
#'     distribution of Harris (1989): the estimates are averaged over their
#'     approximate sampling distribution, centred at the fit. For a Gaussian
#'     with a constant \eqn{\sigma = e^{\gamma_0}}, \eqn{\hat\gamma_0} of
#'     variance \eqn{v}, the variance of a new observation is
#'     \eqn{\mathrm{se}_0^2 + z^\top\Sigma_b z + \hat\sigma^2 e^{2v}}, the
#'     last term being the average of \eqn{\hat\sigma^{*2}} over that
#'     distribution.}
#'   \item{`"plugin"`}{the family at the estimates, averaged over the effects
#'     alone: the estimative interval, which ignores the error in the
#'     estimates and so covers less than the level in a small sample.}
#'   \item{`"bootstrap"`}{the parametric bootstrap predictive distribution
#'     itself: `n_boot` responses simulated from the fit, with the random
#'     effects drawn afresh, each refitted, and the family at each refit
#'     averaged over the effects. `boot_refit = "coefficients"` refits at the
#'     fit's hyperparameters; `"full"` chooses them again on every replica, so
#'     their uncertainty enters too, at the cost of a whole fit each. See
#'     [predictive_bootstrap()].}
#' }
#' Measured at 95 per cent over 400 samples, the coverage of a new
#' observation is 0.887 (plug-in) against 0.925 (averaged) for a Gaussian
#' regression at \eqn{n = 12} (over a further 200 samples, 0.900, 0.950 and
#' 0.945 with the bootstrap), and 0.935 against 0.948 for a Gamma with its
#' dispersion modelled at \eqn{n = 25}.
#' Both intervals for a new group are conditional on the covariance the fit
#' estimated: its uncertainty is not propagated. Measured at 95 per cent over
#' 200 fits of a random intercept at eight observations a group, a new
#' group's predictor is covered 0.898 of the time over 10 groups and 0.943
#' over 40 by `"group"`, and a new observation 0.939 and 0.951 by
#' `"prediction"`, where the confidence interval of the typical group covers
#' 0.45 and 0.26; on a Poisson response the same read 0.909 and 0.937 and,
#' the support being discrete, 0.975 and 0.976. Under a prior that is not
#' Gaussian both `interval = "group"` and `interval = "prediction"` read
#' its quantiles, a Student t prior as a scale mixture of Gaussians and any
#' other prior by Monte Carlo (see [predictive_mixture()]). A term that shares a label with an
#' effect read with its own estimate, or with one written inside a
#' subformula, is rejected.
#'
#' A forecast's standard error, with `se = TRUE`, is the uncertainty of the
#' estimated parameters alone, carried by the delta method through the
#' continued recursion, whose derivative is exact. A forecast also carries
#' the uncertainty of the future scores, which the delta method does not
#' reach, so the standard error and the interval leave it out, and a warning
#' says so.
#' @param object A [StatmodFit()].
#' @param what What to predict: a parameter's name, optionally prefixed
#'   `"link:"`; a moment's name; `"parameter"` (the default) or `"link"`. An
#'   unrecognized name signals an error listing what is available.
#' @param newdata A data frame, or `NULL` for the fitting data. Needs the
#'   covariates the model names but not the response.
#' @param se `TRUE` to report the standard error and an interval as well.
#'   `FALSE` by default.
#' @param level The interval's level, `0.95` by default. Read only where `se`
#'   is `TRUE`.
#' @param random How a [modelterms7::random()] term is read:
#'   `"conditional"` (the default) with each group's own estimated effect, so
#'   a level the fit never saw signals an error; `"zero"` with every effect at
#'   zero, the typical group; `"marginal"` averaged over the prior the fit
#'   estimated, the population average \eqn{E_b[h^{-1}(\eta + z^\top b)]}.
#'   A single string applies to every such term; a character vector named by
#'   the terms' keys chooses term by term, the rest staying conditional. See
#'   the section on new groups.
#' @param interval `"confidence"` (the default), `"group"` or
#'   `"prediction"`. See the section on intervals.
#' @param predictive With `interval = "prediction"`, what the family is
#'   averaged over: `"averaged"` (the default), `"plugin"` or `"bootstrap"`.
#'   See the section on the predictive distribution.
#' @param n_boot The number of bootstrap replicas, `200` by default. Read only
#'   where `predictive = "bootstrap"`.
#' @param boot_refit `"coefficients"` (the default), refitting each replica
#'   at the fit's hyperparameters, or `"full"`, choosing them again. Read only
#'   where `predictive = "bootstrap"`.
#' @param ... Passed to [vcov.StatmodFit()] where `se` is `TRUE`. That is
#'   where `type` chooses between the Bayesian variance, the frequentist one
#'   and the unconditional one. A band around a penalized term is where the
#'   last of the three differs most from the others, since it is the fitted
#'   values that a smoothing parameter moves: measured on a univariate smooth
#'   at \eqn{n = 200}, `type = "unconditional"` widens this interval by 1.1
#'   per cent on average and 7.6 per cent at its widest point.
#' @return With `se = FALSE`, a numeric vector of `nrow(newdata)` values when
#'   `what` names one quantity, and a named list of such vectors for
#'   `"parameter"` and `"link"`.
#'
#'   With `se = TRUE`, a data frame with columns `fit`, `se`, `lower` and
#'   `upper` in place of each vector. `se` is `NA` for an observation whose
#'   predictor reads a coefficient that has no variance, which is the truth
#'   about such a fit, and no gap in the arithmetic.
#'
#'   With `interval = "prediction"`, a data frame with columns `fit` (the
#'   predictive median), `se` (the predictive standard deviation), `lower`
#'   and `upper`.
#' @seealso [fitted.StatmodFit()] for one parameter's fitted values,
#'   [residuals.StatmodFit()] for the matched diagnostic,
#'   [vcov.StatmodFit()] for the variance the standard errors come from.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#'
#' # One parameter, and one of the family's moments.
#' head(predict(fit, "mu"))
#' head(predict(fit, "variance"))
#'
#' # A parameter's predictor instead of its value.
#' head(predict(fit, "link:sigma"))
#'
#' # For a Gaussian the mean is mu and the variance is sigma squared, which
#' # is what the moments come to.
#' all.equal(predict(fit, "mean"), predict(fit, "mu"))
#' all.equal(predict(fit, "variance"), predict(fit, "sigma")^2)
#'
#' # With an interval.
#' head(predict(fit, "mu", se = TRUE))
#'
#' # Every parameter at once, on either scale.
#' str(predict(fit, "parameter"))
#'
#' # A new group: the typical one, and the population average.
#' gg <- data.frame(g = factor(rep(letters[1:12], each = 10)),
#'                  x = rnorm(120))
#' gg$y <- rpois(120, exp(0.5 + 0.3 * gg$x + rnorm(12, sd = 0.6)[gg$g]))
#' fg <- statmod(y ~ x + random(~ 1 | g), distributions7::poisson_distrib(),
#'               gg)
#' new <- data.frame(g = "new", x = 0)
#' predict(fg, "mu", new, random = "zero")
#' predict(fg, "mu", new, random = "marginal")
#'
#' # The interval of a new group's mean, and of a new count from it.
#' predict(fg, "mu", new, random = "zero", interval = "group")
#' predict(fg, "response", new, random = "zero", interval = "prediction",
#'         predictive = "plugin")
#' predict(fg, "response", new, random = "zero", interval = "prediction")
#'
#' # A name the family does not have is refused, and the message says what
#' # is available.
#' try(predict(fit, "median"))
#' @references Harris, I. R. (1989). Predictive fit for natural exponential
#'   families. *Biometrika*, 76, 675--684.
#' @keywords internal
predict.StatmodFit <- function(object, what = "parameter", newdata = NULL,
                               se = FALSE, level = 0.95,
                               random = "conditional",
                               interval = c("confidence", "group",
                                            "prediction"),
                               predictive = c("averaged", "plugin",
                                              "bootstrap"),
                               n_boot = 200L,
                               boot_refit = c("coefficients", "full"), ...) {
  if (is.data.frame(what)) {
    stop(paste0("The second argument of statmod's predict() is 'what', not\n",
                "  'newdata': a fit has several parameters and several\n",
                "  moments, so choosing among them is the ordinary variation.\n",
                "  Write predict(fit, \"mu\", newdata) or\n",
                "  predict(fit, newdata = your_data)."), call. = FALSE)
  }
  if (!is.character(what) || length(what) != 1L) {
    stop("'what' must be a single string.", call. = FALSE)
  }
  interval <- match.arg(interval)
  predictive <- match.arg(predictive)
  boot_refit <- match.arg(boot_refit)
  if (!identical(predictive, "averaged") && !identical(interval, "prediction")) {
    stop(paste0("'predictive' says what a prediction interval averages over, ",
                "and is read only
  with interval = \"prediction\"."),
         call. = FALSE)
  }
  if (identical(predictive, "bootstrap") &&
      (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) ||
       n_boot < 1 || n_boot != round(n_boot))) {
    stop("'n_boot' must be a single positive whole number.", call. = FALSE)
  }
  rm <- random_modes(object@spec, random)
  aside <- rm[rm$mode != "conditional", , drop = FALSE]
  nest <- attr(rm, "nested")
  nest <- nest[nest$mode != "conditional", , drop = FALSE]
  if (nrow(nest) && !identical(interval, "confidence")) {
    stop(sprintf(paste0("interval = \"%s\" is not available for a random ",
                        "effect written inside a\n  term's subformula (%s). ",
                        "random = \"zero\" and \"marginal\" are, with the ",
                        "default\n  interval."), interval, nest$key[1L]),
         call. = FALSE)
  }
  if (identical(what, "response") && !identical(interval, "prediction")) {
    stop(paste0("The response is predicted with interval = \"prediction\": ",
                "what it has is a\n  distribution, not a value to report ",
                "alone."), call. = FALSE)
  }
  if (identical(interval, "prediction") && !identical(what, "response")) {
    stop(paste0("interval = \"prediction\" is an interval for a new ",
                "observation of the\n  response: ask for what = ",
                "\"response\". A parameter's interval for a new group\n",
                "  is interval = \"group\"."), call. = FALSE)
  }
  if (identical(interval, "group") && !nrow(aside)) {
    stop(paste0("interval = \"group\" is the interval of a new group, which ",
                "needs a random\n  effect set aside: random = \"zero\" or ",
                "\"marginal\"."), call. = FALSE)
  }
  if (!identical(interval, "confidence")) {
    if (S7::S7_inherits(object@spec@distrib,
                        distributions7::multivariate_distrib)) {
      stop("interval = \"", interval, "\" is not available for a ",
           "multivariate family.", call. = FALSE)
    }
    if (length(statmod_structural(object@spec))) {
      stop("interval = \"", interval, "\" is not available beside a ",
           "structural term.", call. = FALSE)
    }
  }
  if (identical(interval, "group")) se <- TRUE
  if ((nrow(aside) || nrow(nest)) &&
      length(statmod_structural(object@spec))) {
    stop("random = \"zero\" or \"marginal\" is not available beside a ",
         "structural term: its\n  level is a recursion read at the ",
         "predictor, so setting a random effect aside\n  would move the ",
         "recursion as well.", call. = FALSE)
  }
  spec <- spec_at(object, newdata, need_response = FALSE)
  unseen <- if (nrow(aside)) split(aside$key, aside$param) else NULL
  coef_use <- object@coefficients
  pr <- NULL
  design <- withCallingHandlers(
    if (nrow(nest)) {
      # a random effect inside a subformula enters a term that is not linear
      # in it: the term is evaluated with its effects read in one group's
      # columns, at zero here and at each node of the prior below
      pr <- nested_prepare(object, spec, nest, unseen)
      spec <- pr$spec
      coef_use <- pr$coef
      pr$design
    } else statmod_design(spec, unseen),
    error = function(e) {
      if (grepl("was not present at build time", conditionMessage(e),
                fixed = TRUE)) {
        stop(paste0(conditionMessage(e), "\n  A group the fit never saw ",
                    "has no effect of its own: predict(random = \"zero\")\n",
                    "  reads it as the typical group and random = ",
                    "\"marginal\" averages over the prior."),
             call. = FALSE)
      }
    })
  # a random effect set aside contributes nothing to the predictor at ANY
  # row, a group the fit saw included: its columns are zero, so its
  # coefficients enter neither the predictor nor the delta method
  for (i in seq_len(nrow(aside))) {
    idx <- design[[aside$param[i]]]$blocks[[aside$key[i]]]
    if (length(idx)) design[[aside$param[i]]]$X[, idx] <- 0
  }
  mt <- aside[aside$mode == "marginal", , drop = FALSE]
  nm <- which(nest$mode == "marginal")
  nodes <- NULL
  eta_node <- NULL
  if (nrow(mt) || length(nm)) {
    priors <- c(lapply(seq_len(nrow(mt)), function(i)
      random_prior(spec, design, object, mt$param[i], mt$key[i])),
      lapply(nm, function(j) nested_prior(spec, design, object, nest[j, ])))
    nodes <- random_nodes(priors)
  }
  # the coefficients and the predictors at node k of the nested effects
  # averaged over: every row reads its effect in the first level's columns
  if (length(nm)) {
    coef_node <- function(k) {
      cf <- coef_use
      for (j in seq_along(nm)) {
        cl <- pr$cols[[nm[j]]]
        cf[[cl$param]][cl$first] <- nodes$b[[nrow(mt) + j]][k, ]
      }
      cf
    }
    eta_node <- function(k) statmod_eta(spec, design, coef_node(k))$eta
  }
  # A STRUCTURAL TERM HAS NO BLOCK TO REAPPLY: its contribution is the state
  # a recursion has reached, so new rows CONTINUE the series rather than
  # being read on their own. Running the ordinary assembly there returned
  # the fitting data's values whatever `newdata` held, which is why the two
  # paths are separated rather than merged.
  # rows carrying the response were rebuilt as a series of their own by
  # statmod_respec(), and are read as the fitting rows are
  cont <- !is.null(newdata) && length(attr(design, "structural")) &&
    !isTRUE(statmod_response_known(spec@response))
  ep <- if (cont) {
    statmod_eta_continued(object, spec, design, deriv = isTRUE(se))
  } else statmod_eta(spec, design, coef_use)
  # A FORECAST'S STANDARD ERROR is the uncertainty of the parameters alone,
  # carried by the delta method through the continued recursion. A forecast
  # also carries the uncertainty of the future scores, which no delta method
  # reaches, so the number is reported with a warning saying what it leaves
  # out (Giovanni, 2026-10-03: a warning, not the error it used to be).
  if (cont && isTRUE(se)) {
    warning("The standard error of a prediction past the series is that of ",
            "the estimated\n  parameters alone. A forecast also carries the ",
            "uncertainty of the future\n  scores, which this standard error ",
            "and its interval leave out.", call. = FALSE)
  }
  params <- spec@distrib@params
  if (identical(interval, "prediction")) {
    check_trials_known(spec@distrib, "A prediction interval")
    comps <- switch(predictive,
      averaged = predictive_mixture(
        object, spec, design, aside, ep$eta,
        predictive_cov(object, spec, design, list(), ...)$fixed),
      plugin = predictive_mixture(object, spec, design, aside, ep$eta, NULL),
      bootstrap = predictive_bootstrap(object, spec, design, aside,
                                       as.integer(n_boot), boot_refit))
    return(predictive_response(spec, comps, level))
  }
  if (isTRUE(se)) {
    su <- predict_se(object, spec, design, ep, level, ...)
    # A NEW GROUP'S PARAMETER varies around the typical group's by what its
    # effects add, z' Sigma_b z in its own equation, on top of the estimates'
    # own uncertainty
    # A NEW GROUP'S PARAMETER varies around the typical group's by what its
    # effects add, on top of the estimates' own uncertainty: its interval and
    # standard deviation are those of h^-1(eta*), see group_interval()
    if (identical(interval, "group")) {
      su <- group_interval(object, spec, design, aside, ep, su, level, ...)
    }
    # a marginal parameter is averaged over the effects: the fit and the ends
    # of the interval are the fixed part's carried through that average,
    # which is monotone in the predictor, and the standard error multiplies
    # the predictor's by the averaged derivative of the inverse link
    for (p in unique(mt$param)) {
      # for a new group the interval is its parameter's, already set, and only
      # the fit is the population average
      if (identical(interval, "group")) {
        su[[p]]$fit <- random_marginal(spec, design, ep, mt, nodes, p)
        next
      }
      eta_at <- function(col) {
        e <- ep$eta
        e[[p]] <- su[[p]][[col]]
        e
      }
      su[[p]]$fit <- random_marginal(spec, design, ep, mt, nodes, p)
      su[[p]]$se <- abs(random_marginal(spec, design, ep, mt, nodes, p,
                                        deriv = TRUE)) * su[[p]]$se_eta
      ends <- cbind(
        random_marginal(spec, design, ep, mt, nodes, p, eta_at("eta_lower")),
        random_marginal(spec, design, ep, mt, nodes, p, eta_at("eta_upper")))
      su[[p]]$lower <- pmin(ends[, 1L], ends[, 2L])
      su[[p]]$upper <- pmax(ends[, 1L], ends[, 2L])
    }
    # with nested effects averaged over, the predictor is not a shift of the
    # fixed part, so the delta method reads the design at every node
    if (length(nm)) {
      if (nrow(mt)) {
        stop(paste0("se = TRUE is not available where random effects ",
                    "written in an equation and\n  inside a subformula are ",
                    "both averaged over. Each kind alone is."),
             call. = FALSE)
      }
      for (p in unique(nest$param[nm])) {
        fitv <- random_marginal(spec, design, ep, mt, nodes, p,
                                eta_node = eta_node)
        ms <- nested_marginal_se(object, pr, coef_node, eta_node, nodes$w, p,
                                 fitv, level, ...)
        su[[p]][c("fit", "se", "lower", "upper")] <- ms
      }
    }
    return(se_answer(su, what, params, spec))
  }
  if ((nrow(mt) || length(nm)) && !identical(what, "link") &&
      !startsWith(what, "link:")) {
    return(random_marginal(spec, design, ep, mt, nodes, what,
                           eta_node = eta_node))
  }
  # a predictor averaged over nested effects is the average of the
  # predictors at the nodes, the term not being linear in them
  if (length(nm)) {
    n <- spec@n_obs
    eb <- lapply(ep$eta, function(e) numeric(n))
    for (k in seq_along(nodes$w)) {
      ek <- eta_node(k)
      for (p in names(eb)) {
        eb[[p]] <- eb[[p]] + nodes$w[k] * rep_len(as.numeric(ek[[p]]), n)
      }
    }
    ep$eta <- eb
  }

  if (identical(what, "link")) return(ep$eta)
  if (identical(what, "parameter")) return(ep$theta)

  if (startsWith(what, "link:")) {
    p <- substring(what, 6L)
    if (!p %in% params) stop(unknown_what(p, params), call. = FALSE)
    return(ep$eta[[p]])
  }
  if (what %in% params) return(rep_len(ep$theta[[what]], spec@n_obs))

  mom <- predict_moments()
  if (what %in% names(mom)) {
    check_trials_known(spec@distrib, sprintf("The %s of the response", what))
    # only a missing method becomes the friendly message: a catch-all here
    # would report any failure as "this family has no such moment", which is
    # the shape of error the toolkit records as worse than none
    v <- tryCatch(mom[[what]](spec@distrib, ep$theta),
                  error = function(e) {
                    if (grepl("method", conditionMessage(e), fixed = TRUE)) {
                      stop(sprintf(paste0(
                        "'%s' does not implement %s(). Its parameters are\n",
                        "  available: ask for one of %s."),
                        spec@distrib@distrib_name, what,
                        paste(params, collapse = ", ")), call. = FALSE)
                    }
                    stop(e)
                  })
    return(rep_len(v, spec@n_obs))
  }
  stop(unknown_what(what, params), call. = FALSE)
}
S7::method(predict, StatmodFit) <- predict.StatmodFit


#' The Message for an Unrecognized Prediction Target
#'
#' @description
#' Builds the error [predict.StatmodFit()] signals when `what` names
#' nothing it can compute, listing this family's own parameter names, the
#' five moments and the two collective targets.
#'
#' @details
#' The family's parameters are listed by name rather than described, since
#' they differ from family to family and are the commonest thing a caller
#' means. A data frame passed where `what` belongs, which is the mistake
#' [stats::predict()]'s argument order invites, is recognized and named
#' separately.
#'
#' @param what What was asked for, for the message.
#' @param params The family's parameter names, in the family's order.
#'
#' @return A single string, ready for [stop()].
#'
#' @seealso [predict.StatmodFit()], the caller,
#'   [predict_moments()] for the moment names listed.
#'
#' @keywords internal
unknown_what <- function(what, params) {
  sprintf(paste0("'%s' is neither a parameter of this distribution nor a\n",
                 "  moment. The parameters are: %s.\n",
                 "  The moments are: %s.\n",
                 "  \"parameter\" and \"link\" give all of them at once."),
          what, paste(params, collapse = ", "),
          paste(names(predict_moments()), collapse = ", "))
}


#' @title The Fitted Values of a Model
#' @name fitted.StatmodFit
#' @description
#' One distribution parameter's fitted values, as a vector of the data's
#' length.
#' @details
#' The result is one vector, never the whole set, so that this and
#' [residuals.StatmodFit()] are a matched pair and a diagnostic drawn from
#' them needs no unpacking.
#'
#' The default is the **first** parameter, never the mean. A family may
#' have
#' no mean, a Cauchy being one, and a default that fails on a legitimate
#' family is worse than one that always answers.
#'
#' The whole set at once is `predict(fit, "parameter")`, and the mean, where
#' it exists, is `predict(fit, "mean")`.
#' @param object A [StatmodFit()].
#' @param what Which distribution parameter, a string naming one of the
#'   family's, or `NULL` for the first.
#' @param ... Unused.
#' @return A numeric vector of length `nobs(object)`, that parameter's fitted
#'   values on its own scale.
#' @seealso [predict.StatmodFit()] for the other quantities and for new data,
#'   [residuals.StatmodFit()] for the matched diagnostic.
#' @keywords internal
fitted.StatmodFit <- function(object, what = NULL, ...) {
  th <- object@fitted
  ps <- object@spec@distrib@params
  if (is.null(what)) what <- ps[[1L]]
  if (!is.character(what) || length(what) != 1L || !what %in% ps) {
    stop(sprintf(paste0("'what' must be one of the distribution's",
                        " parameters (%s).\n  The whole set at once is",
                        " predict(fit, \"parameter\")."),
                 paste(ps, collapse = ", ")), call. = FALSE)
  }
  rep_len(as.numeric(th[[what]]), object@spec@n_obs)
}
S7::method(fitted, StatmodFit) <- fitted.StatmodFit



#' The Uncertainty of a Predicted Predictor
#'
#' @description
#' The standard error of each equation's linear predictor, and of the
#' parameter it gives, with an interval.
#'
#' @details
#' # The delta method, equation by equation
#'
#' An equation's predictor is \eqn{\eta_{ip} = x_{ip}'\beta_p}, so its
#' variance is \eqn{x_{ip}' V_{pp} x_{ip}}, with \eqn{V} the variance of the
#' coefficients **as estimated**: the coordinates, since the design is
#' written in them, never the quantities [coef.StatmodFit()] reports by
#' default.
#'
#' The equations do not mix. One equation's predictor reads that equation's
#' coefficients alone, whatever the covariance between the blocks.
#'
#' # A term whose block moves needs no special case
#'
#' Its block is the Jacobian \eqn{\partial\eta/\partial\beta} by
#' construction, that being why a linear fit on it is a Gauss-Newton step,
#' so the row already is the derivative and the delta method is exact to
#' first order. That covers [modelterms7::seg()], [modelterms7::jseg()] and
#' [modelterms7::nl()], including the parameters of a nonlinear term
#' developed over covariates.
#'
#' Measured against a numerical derivative of the predictor in the estimated
#' coefficients, which shares no arithmetic with the design row: 1.7e-12 on a
#' parametric block, 8.9e-11 on a smooth with a random effect, 1.3e-10 on a
#' `seg()` and 1.1e-11 on an `nl()` with a ridge.
#'
#' # The interval
#'
#' Built on the scale the equation is written on and mapped back through the
#' link, as every interval in this toolkit is. A scale rides a logarithm, so
#' its lower end cannot come out negative.
#'
#' # A coefficient with no variance carries none forward
#'
#' A discontinuous break-point term's block is a working linearization and is
#' held out of [vcov.StatmodFit()], so every observation whose predictor
#' reads it reports `NA` for its standard error. That is the truth about such
#' a fit, no gap in the arithmetic.
#'
#' @param object A fitted model.
#' @param spec The specification the prediction is made under.
#' @param design Its design.
#' @param ep The predictors and parameters, as [statmod_eta()]
#'   returns them.
#' @param level The interval's level.
#' @param ... Passed to [vcov()].
#'
#' @return A named list, one entry per distribution parameter, each a data
#'   frame with `fit`, `se`, `lower` and `upper` on the
#'   link scale and the same four on the parameter scale.
#'
#' @seealso [predict.StatmodFit()], [vcov.StatmodFit()]
#'
#' @keywords internal
predict_se <- function(object, spec, design, ep, level = 0.95, ...) {
  V <- vcov(object, readable = FALSE, ...)
  z <- stats::qnorm((1 + level) / 2)
  links <- spec@distrib@link_params
  out <- stats::setNames(vector("list", length(spec@distrib@params)),
                         spec@distrib@params)
  for (p in spec@distrib@params) {
    d <- design[[p]]
    n <- spec@n_obs
    eta <- rep_len(as.numeric(ep$eta[[p]]), n)
    v <- rep(NA_real_, n)
    key <- if (d$npar) paste(p, d$coef_names, sep = ":") else character(0)
    X <- if (d$npar) as.matrix(d$X) else matrix(0, n, 0L)
    # A SCORE-DRIVEN TERM ADDS ITS LEVEL to this equation's predictor, and
    # the level is no column of the design: its derivative is the forward
    # Jacobian of the recursion, in the coefficients of EVERY equation (the
    # scores are read at all the predictors) and the term's free parameters,
    # and it replaces the design row
    fl <- structural_se_columns(spec, design, ep, p, object@coefficients)
    if (!is.null(fl) && nrow(fl$X) == n) {
      X <- fl$X
      key <- fl$key
    }
    # past the series the rows come from the continued recursion instead
    cc <- ep$cont_cols[[p]]
    if (!is.null(cc) && nrow(cc$X) == n) {
      X <- cc$X
      key <- cc$key
    }
    if (length(key) && all(key %in% rownames(V))) {
      Vp <- as.matrix(V[key, key, drop = FALSE])
      if (nrow(X) == n) {
        v <- row_quad(X, Vp, X)
        v[!is.na(v) & v < 0] <- 0
      }
    }
    se_eta <- sqrt(v)
    lo <- eta - z * se_eta
    hi <- eta + z * se_eta
    g <- links[[p]]
    th <- linkfunctions7::linkinv(g, eta)
    # the interval is the predictor's, carried through the link, so it keeps
    # the parameter inside its own set; the standard error beside it is the
    # delta method, which the interval does not use
    ends <- cbind(linkfunctions7::linkinv(g, lo),
                  linkfunctions7::linkinv(g, hi))
    out[[p]] <- data.frame(
      eta = eta, se_eta = se_eta, eta_lower = lo, eta_upper = hi,
      fit = as.numeric(th),
      se = abs(linkfunctions7::dlinkinv(g, eta)) * se_eta,
      lower = pmin(ends[, 1L], ends[, 2L]),
      upper = pmax(ends[, 1L], ends[, 2L]))
  }
  out
}


#' The Shape a Prediction With Its Uncertainty Comes Back In
#'
#' @description
#' The rows [predict_se()] computed, reduced to what was asked for.
#'
#' @details
#' The vocabulary is the one a prediction without uncertainty answers, minus
#' the moments: a moment's delta method needs its derivative in every
#' parameter, which \pkg{distributions7} does not offer as a generic, and an
#' interval assembled from what is at hand instead would be a number nobody
#' could check.
#'
#' @param su The per-parameter tables.
#' @param what What was asked for.
#' @param params The distribution's parameters.
#' @param spec The specification, for the message.
#'
#' @return A data frame, or a named list of them.
#'
#' @keywords internal
se_answer <- function(su, what, params, spec) {
  cols <- c("fit", "se", "lower", "upper")
  lcols <- c("eta", "se_eta", "eta_lower", "eta_upper")
  ren <- function(d, cc) stats::setNames(d[, cc, drop = FALSE], cols)
  if (identical(what, "parameter")) {
    return(lapply(su, ren, cc = cols))
  }
  if (identical(what, "link")) {
    return(lapply(su, ren, cc = lcols))
  }
  if (what %in% params) return(ren(su[[what]], cols))
  if (startsWith(what, "link:")) {
    p <- substring(what, 6L)
    if (!p %in% params) stop(unknown_what(p, params), call. = FALSE)
    return(ren(su[[p]], lcols))
  }
  stop(sprintf(paste0("'%s' has no standard error here. A moment's would ",
                      "need its\n  derivative in every parameter, which the ",
                      "distribution does not\n  offer; ask for a parameter ",
                      "(%s), for \"parameter\", or for \"link\"."),
               what, paste(params, collapse = ", ")), call. = FALSE)
}


#' The Derivative Row of an Equation Carrying a Filter
#'
#' @description
#' The derivative of the predictor of the equation carrying a score-driven
#' term, one row per observation, in every coordinate of the variance matrix
#' it moves with: the coefficients of every equation and the term's free
#' parameters on their unconstrained scale.
#'
#' @details
#' A filter's level is a recursion, not a column, so it has no row of a
#' design. Its derivative is the forward Jacobian of the recursion, which
#' [filter_joint_jacobian()] returns. Every equation's coefficients enter it,
#' not only those of the filter's own equation: the scores that drive the
#' recursion are read at the predictors of every equation, so a coefficient
#' of the scale moves the level of a filter in the mean. Leaving those
#' columns out, on the gaussian score-driven model of the Nile flow, gave a
#' standard error of the filtered mean up to 12 per cent too large.
#'
#' A parameter that an intercept in the same equation holds is not estimated
#' and is not in that matrix, so it is not here either.
#'
#' @param spec The specification.
#' @param design Its design.
#' @param ep The predictors, as [statmod_eta()] returns them.
#' @param p The distribution parameter whose equation is being read.
#' @param coef The coefficients the predictors were evaluated at.
#'
#' @return A list with `X` (the derivative rows) and `key` (the names of
#'   their columns in the variance matrix), or `NULL` where the equation
#'   carries no filter.
#'
#' @seealso [predict_se()], [filter_joint_jacobian()]
#'
#' @keywords internal
structural_se_columns <- function(spec, design, ep, p, coef) {
  fs <- ep$filters
  if (!length(fs) ||
      !any(vapply(fs, function(f) identical(f$param, p), logical(1)))) {
    return(NULL)
  }
  fj <- filter_joint_jacobian(spec, design, coef)
  if (is.null(fj)) return(NULL)
  list(X = fj$J, key = fj$key)
}


#' The Forward Jacobian of a Filter in Every Coordinate
#'
#' @description
#' The derivative of the predictor that a score-driven term produces, in the
#' coefficients of every equation followed by the term's free parameters on
#' their unconstrained scale, at each observation, together with the static
#' rows of every equation in the same columns.
#'
#' @details
#' The rows are the ones the joint information is assembled from
#' ([statmod_full_information()]): each equation's design placed in its own
#' columns, and for the filter's equation the Jacobian that
#' [modelterms7::term_curvature()] propagates beside the state. A model
#' carries at most one structural term, so every other equation is static
#' and its design is its derivative.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#' @param coef The coefficients.
#'
#' @return `NULL` where the model carries no filter, otherwise a list with
#'   `J` (the filter equation's rows), `V` (every equation's rows, the
#'   filter's being `J`), `H` (the family's second derivatives on the link
#'   scale), `ap` (the index of the filter's parameter), `nb` (the number of
#'   coefficient columns), `key` (the names of the columns in the variance
#'   matrix) and `free` (the term's free parameters), all columns restricted
#'   to the free ones.
#'
#' @seealso [structural_se_columns()], [continued_deriv_inputs()]
#'
#' @keywords internal
filter_joint_jacobian <- function(spec, design, coef) {
  jd <- joint_design_rows(spec, design, coef)
  if (is.null(jd)) return(NULL)
  d <- spec@distrib
  th <- jd$ev$theta
  gl <- distributions7::distrib_gradient(d, spec@response, th, scale = "link",
                                         threads = spec@threads)
  H <- distributions7::distrib_hessian(d, spec@response, th, scale = "link",
                                       threads = spec@threads)
  D3 <- distributions7::distrib_deriv3(d, spec@response, th, scale = "link",
                                       threads = spec@threads)
  cv <- filter_curvature(spec, design, jd$f, jd$ap, jd$V, gl, H, D3, th)
  V <- jd$V
  V[[jd$ap]] <- cv$jacobian
  V <- lapply(V, function(x) x[, jd$keep, drop = FALSE])
  # the derivative of what drives the filter, for a continuation's rows
  Hd <- filter_driving(spec, th, jd$ap, filter_scaling(jd$f$tm), gl, H)$H
  params <- jd$params
  ckey <- unlist(lapply(params, function(q) {
    if (design[[q]]$npar) paste(q, design[[q]]$coef_names, sep = ":") else
      character(0)
  }), use.names = FALSE)
  lb <- tryCatch(jd$f$tm@label, error = function(e) "")
  tkey <- paste(jd$f$param, if (length(lb) == 1L && nzchar(lb))
    paste(lb, jd$zn, sep = ".") else jd$zn, sep = ":")
  list(J = V[[jd$ap]], V = V, H = H, Hd = Hd, ap = jd$ap, nb = jd$nb,
       key = c(ckey, tkey)[jd$keep],
       free = jd$zn[jd$keep[jd$keep > jd$nb] - jd$nb])
}


#' Predictors at Rows That Continue the Series
#'
#' @description
#' The linear predictors and the parameters at new rows for a model carrying
#' a structural term, with the term's own recursion carried forward from
#' where the fitting data left it.
#'
#' @details
#' The static part is the ordinary assembly: each term's block reapplied at
#' the new rows, the offsets re-evaluated there. What cannot be reapplied is
#' the structural term, whose contribution at one row is the state a
#' recursion has reached over every row before it. The fit is therefore run
#' once at the observed rows to recover that state -- the level and the
#' driving quantity at each of them, the level being the difference between
#' the filtered predictor and the static one -- and the term is asked to
#' continue from there through
#' [modelterms7::term_continue()].
#'
#' Only rows without the response are continued, and they must come after
#' the observed series. Rows that carry it are a series of their own, rebuilt
#' by [statmod_respec()] and read as the fitting rows are, which is why
#' `predict(fit, newdata = <the fitting data>)` returns the fitted values.
#'
#' A term whose contribution is a likelihood mixed over latent states is
#' rejected: what such a term reports at an observed row is a posterior over
#' states, which past the data is a predictive distribution, no single
#' value.
#'
#' With `deriv = TRUE` the continuation is differentiated as well, in the
#' coefficients of every equation and the term's free parameters
#' ([continued_deriv_inputs()]), which is what a standard error of the
#' forecast reads.
#'
#' @param fit The fitted model.
#' @param spec The specification at the new data.
#' @param design Its design.
#' @param deriv Whether to return the forecast's derivative rows.
#'
#' @return A list shaped as [statmod_eta()]'s, without the
#'   memoized filter objects, plus `cont_cols`: for each equation whose
#'   filter was continued with `deriv = TRUE`, a list with `X` (the
#'   derivative of the predictor at the new rows in the coefficients of
#'   every equation and the term's free parameters) and `key` (the names of
#'   its columns in the variance matrix).
#'
#' @seealso [predict.StatmodFit()]
#'
#' @keywords internal
statmod_eta_continued <- function(fit, spec, design, deriv = FALSE) {
  cont_cols <- list()
  nd <- spec@newdata
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  coef <- fit@coefficients
  su <- attr(design, "structural")
  for (u in su) {
    if (!identical(u$kind, "filter")) {
      stop(sprintf(paste0("'%s' reports a likelihood mixed over latent ",
                          "states, not a predictor, so\n  it has no value ",
                          "to continue past the series. predict(fit) at the ",
                          "observed\n  rows gives its posterior-weighted ",
                          "predictor."), u$term), call. = FALSE)
    }
  }

  # the static part at the new rows
  eta <- stats::setNames(vector("list", length(params)), params)
  theta <- eta
  d2 <- statmod_design_at(spec, coef, design)
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

  # the state the observed part ended at
  ospec <- fit@spec
  odesign <- statmod_design(ospec)
  oep <- statmod_eta(ospec, odesign, coef)
  ost <- statmod_structural_state(odesign)
  for (f in oep$filters) {
    tm <- ospec@terms[[f$param]][[f$term]]
    psi <- structural_psi(tm, ost$zeta[[f$term]])
    p <- f$param
    f_past <- as.numeric(f$eta) - as.numeric(f$eta_static)
    # the driving quantity, read at the predictor the recursion produced --
    # the same callback the filter itself was handed
    s_past <- vapply(seq_along(f_past),
                     function(i) f$cb$score(f$eta[[i]], i), numeric(1))
    dv <- if (isTRUE(deriv)) {
      continued_deriv_inputs(ospec, odesign, coef, f, ost, tm)
    } else NULL
    cont <- modelterms7::term_continue(tm, psi, f_past, s_past, nd,
                                       deriv = dv$deriv)
    eta[[p]] <- eta[[p]] + as.numeric(cont)
    if (!is.null(dv)) {
      # the forecast's derivative row: the static part in the equation's
      # own columns plus what the continued level adds in every column
      jac <- attr(cont, "jacobian")
      if (d2[[p]]$npar) {
        jac[, dv$col] <- jac[, dv$col] + as.matrix(d2[[p]]$X)
      }
      cont_cols[[p]] <- list(X = jac, key = dv$key)
    }
    theta[[p]] <- linkfunctions7::linkinv(links[[p]], eta[[p]])
  }
  list(eta = eta, theta = theta, filters = list(), regimes = list(),
       eta_static = eta, cont_cols = cont_cols)
}


#' What the Derivative of a Continued Filter Starts From
#'
#' @description
#' The derivatives of the level and of the score at the observed rows, and of
#' the term's parameters, in the coordinates the variance matrix of the fit
#' is written in: the coefficients of every equation followed by the term's
#' free parameters on their unconstrained scale.
#'
#' @details
#' The coordinates are the ones [structural_se_columns()] uses for a
#' prediction at the observed rows, so a forecast and a fitted value are
#' read against the same rows and columns of the variance. The derivative of
#' the filtered predictor at an observed row is the forward Jacobian of
#' [filter_joint_jacobian()]; the level's derivative is that minus the
#' equation's design row. The score drives the recursion and depends on the
#' predictors of every equation, so its derivative is
#' \eqn{\sum_b \ell_{pb,t} V_{b,t}}, with \eqn{V_{b,t}} the derivative row
#' of equation \eqn{b} and \eqn{\ell_{pb,t}} the family's second derivative
#' on the link scale. A parameter's derivative in its own unconstrained
#' coordinate is its link's, and one for a coefficient of a development.
#'
#' @param ospec,odesign The fit's specification and design.
#' @param coef The coefficients.
#' @param f The filter object at the observed rows.
#' @param ost The structural state.
#' @param tm The term.
#'
#' @return A list with `deriv` (as [modelterms7::term_continue()] takes
#'   it), `col` (the columns of the filter equation's coefficients) and
#'   `key` (the names of the columns in the variance matrix), or `NULL`
#'   where the fit carries no filter.
#'
#' @seealso [statmod_eta_continued()], [structural_se_columns()]
#'
#' @keywords internal
continued_deriv_inputs <- function(ospec, odesign, coef, f, ost, tm) {
  fj <- filter_joint_jacobian(ospec, odesign, coef)
  if (is.null(fj)) return(NULL)
  params <- ospec@distrib@params
  ap <- fj$ap
  p <- params[ap]
  n <- ospec@n_obs
  npar <- vapply(odesign, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  col <- offs[ap] + seq_len(npar[ap])
  # the equation's static rows in the joint columns
  S <- matrix(0, n, ncol(fj$J))
  if (npar[ap]) S[, col] <- as_dense(odesign[[p]]$X)
  ds <- matrix(0, n, ncol(fj$J))
  for (b in seq_along(params)) {
    ds <- ds + rep_len(fj$Hd[[hess_key(params, ap, b)]], n) * fj$V[[b]]
  }
  z <- ost$zeta[[f$term]]
  nm <- names(z)
  links <- modelterms7::term_links(tm)
  dpsi <- matrix(0, length(nm), ncol(fj$J))
  for (i in seq_along(fj$free)) {
    j <- fj$free[i]
    dpsi[match(j, nm), fj$nb + i] <-
      linkfunctions7::dlinkinv(links[[j]], z[[j]])
  }
  list(deriv = list(df_past = fj$J - S, ds_past = ds, dpsi = dpsi),
       col = col, key = fj$key)
}


#' Whether New Rows Carry the Response
#'
#' @description
#' `TRUE` where every row has one, `FALSE` where none does, and an
#' error where some do.
#'
#' @details
#' It is what separates a re-reading of a model on another series from a
#' continuation of the one it was fitted to, and the separation has to be
#' all-or-nothing: a filter's recursion at one row reads the rows before it,
#' so a frame carrying the response on some rows only describes neither
#' operation, and answering it would mean choosing a reading the caller did
#' not ask for.
#'
#' @param y The response as the specification carries it.
#'
#' @return `TRUE` or `FALSE`.
#'
#' @seealso [statmod_eta_continued()]
#'
#' @keywords internal
statmod_response_known <- function(y) {
  if (is.null(y)) return(FALSE)
  ok <- if (is.matrix(y)) stats::complete.cases(y) else !is.na(y)
  if (all(ok)) return(TRUE)
  if (!any(ok)) return(FALSE)
  stop("the new data carries the response on some rows and not others, and ",
       "the two\n  readings differ: rows with it are read by running the ",
       "filter over them, rows\n  without it by continuing the fitted ",
       "series. Ask for one or the other.", call. = FALSE)
}
