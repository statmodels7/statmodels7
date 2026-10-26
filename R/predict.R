#' @include report.R predict_objects.R
#' @importFrom stats predict fitted
NULL

#' The Quantities That a Fit Can Predict
#'
#' @description
#' Returns the table against which [predict.StatmodFit()] resolves a moment
#' name. The table maps each recognized name to the \pkg{distributions7}
#' generic that computes it.
#'
#' @details
#' The table is written once so that the recognized names and the functions
#' that they call cannot disagree, and so that [unknown_what()] can list the
#' names in its message.
#'
#' @return A named list of functions, keyed `"mean"`, `"variance"`,
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
#' any of its moments, every parameter or linear predictor at once, or a new
#' observation of the response, at the fitting data or at new data, with
#' standard errors and intervals on request.
#' @details
#' # Quantities that can be predicted
#'
#' `what` takes
#' \describe{
#'   \item{a parameter's name}{`"mu"`, `"sigma"` or `"alpha"`, whatever the
#'     family calls them. A parameter is available for every family, because
#'     the model fits it, and it exists even where a moment does not.}
#'   \item{a moment's name}{`"mean"`, `"variance"`, `"std_dev"`,
#'     `"skewness"`, `"kurtosis"`. A moment is available where the family has
#'     one, and the result is `NaN` or `NA` where the moment does not exist.
#'     A Cauchy's mean is
#'     `NaN`, which is the correct answer for it. A moment has no standard
#'     error. For a binomial fitted with `cbind()`, a moment needs the
#'     number of trials of each row, so `newdata` has to carry the variables
#'     with which the response is written.}
#'   \item{`"parameter"`}{every parameter at once, as a named list. This is
#'     the default.}
#'   \item{`"link"`}{every linear predictor at once, before the inverse
#'     link.}
#'   \item{`"response"`}{the response itself, with `interval =
#'     prediction_interval()`, which it requires: its median and the ends of an
#'     interval in which a new observation falls. See the section on
#'     intervals.}
#' }
#' A parameter's name may be prefixed by `"link:"` to return its
#' predictor instead of its value, as `"link:sigma"`.
#'
#' # The argument order departs from [stats::predict()]
#'
#' There the second argument is `newdata`; here it is `what`. A statmod fit
#' has several parameters and several moments, so choosing among them is the
#' ordinary variation and predicting on new data is the occasional one.
#'
#' A data frame passed second is detected, and the error names the argument
#' order.
#'
#' # New data
#'
#' New data is passed through each term's blueprint, so a factor keeps the
#' levels and contrasts with which it was fitted, a spline its knots, and a
#' basis its reparametrization. The new frame supplies the values at which the
#' terms are evaluated and does not change how they are built. Offsets are
#' evaluated again at the new rows.
#'
#' The new frame needs the variables of the equations that are read. Where
#' `what` names one parameter and `interval` is [confidence_interval()], the other
#' equations are not read, and the variables that only they use may be absent
#' from `newdata`. The response is not needed, except that a binomial fitted
#' with `cbind()` reads the number of trials from it for a moment or a
#' prediction interval.
#'
#' # A score-driven term is predicted past the series
#'
#' The contribution of such a term at one row is the state that a recursion
#' has reached, so new rows without the response continue the series instead
#' of being read on their own. Each row is placed by its own time within its
#' own group, and must come after every observed time of that group. A row
#' that falls inside the observed series is rejected, since there the response
#' is known and the filter has to be run, not continued. A group that is
#' absent from the fitting data is rejected as well, and a term built without
#' `time` cannot be continued, because the place of a row in the series is
#' then its position in the data frame.
#'
#' New rows that carry the response are a series of their own: the
#' structural term is rebuilt on them with the fitted parameters, and they
#' are read as the fitting rows are, including the standard errors. If the
#' fitting data are passed as `newdata`, the result is that of `predict(fit)`;
#' this holds for a [modelterms7::regime()] term as well, whose prediction is
#' the posterior-weighted predictor. A frame that carries the response on some
#' rows only is rejected. A [modelterms7::regime()] term cannot be continued
#' past the series, since what it reports at an observed row is a posterior
#' over states.
#'
#' Beyond the data the score sits at its conditional mean of zero, which the
#' model's own definition guarantees, so the continuation is the
#' deterministic recursion and involves no simulation.
#'
#' # A new group, and the typical one
#'
#' A random-effect term predicts a group that is present in the fitting data
#' with that group's estimated effect. A group that is absent from the
#' fitting data has none, and what is put in
#' its place depends on `random`. `"zero"` sets every effect of the
#' term to zero, at every row, which is the prediction for the typical group
#' (lme4's and glmmTMB's `re.form = NA`). `"marginal"` averages over the prior
#' instead,
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
#' probability is below 1e-8 and is `NA`, with a warning, elsewhere. A term
#' that shares a covariance with others through a label has no prior of its
#' own, so `"marginal"` cannot average over it; it can be read
#' `"conditional"` or set aside at `"zero"`. The standard error of a marginal
#' parameter is the delta method on the fixed part, conditional on the
#' prior's scale.
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
#' contribution, carried through the inverse link, over the nodes that the
#' rule above builds for the priors set aside; a prior of another family
#' signals an error. A predictor (`"link"`) is the average of the predictors
#' at the nodes. The standard error of a marginal parameter is the delta method of
#' that average, whose gradient sums the design at each node, and its
#' interval is the delta method's on the link scale carried through the
#' link. [group_interval()] and [prediction_interval()] are not available for such
#' an effect, and a standard error is not available where an average is taken
#' over effects of both kinds.
#'
#' # Intervals and predictive distributions
#'
#' The object given to `interval` selects the uncertainty that the interval
#' describes: [confidence_interval()] the estimates alone, [group_interval()]
#' the estimates and the effects of a new group, and [prediction_interval()] a
#' new observation of the response. The topic [predict_intervals] compares
#' the three and gives the formulas. A prediction interval holds a predictive
#' object that sets how the estimation uncertainty enters: [pred_averaged()]
#' (the default), [pred_plugin()] or [pred_bootstrap()]. The topic
#' [predictive_distributions] compares them and states what each costs.
#'
#' [group_interval()] and [prediction_interval()] are rejected for a
#' multivariate family and for a model with a structural term. The structural
#' term also rejects `random = "zero"` and `random = "marginal"`: its level is
#' a recursion read at the predictor, so setting a random effect aside would
#' move the recursion as well.
#'
#' A forecast's standard error, with `se = TRUE`, is only the uncertainty of
#' the estimated parameters, carried by the delta method through the
#' continued recursion, whose derivative is exact. A forecast also carries
#' the uncertainty of the future scores, which the delta method does not
#' reach, so the standard error and the interval leave it out, and a warning
#' is issued for the equations that the continued filter reaches.
#' @param object A [StatmodFit()].
#' @param what What to predict: a parameter's name, optionally prefixed
#'   `"link:"`; a moment's name; `"parameter"` (the default) or `"link"`; or
#'   `"response"`, which requires `interval = prediction_interval()` and is the
#'   only value that this interval accepts. An unrecognized name signals an error
#'   listing what is available.
#' @param newdata A data frame, or `NULL` for the fitting data. It needs the
#'   variables of the equations that are read: where `what` names one
#'   parameter and `interval` is [confidence_interval()], the variables that only the
#'   other equations use may be absent. It does not need the response, except
#'   that a moment or a prediction interval of a binomial fitted with
#'   `cbind()` reads the number of trials from it. For a model with a
#'   structural term, rows that all carry the response are read as a series of
#'   their own, and rows without it continue the fitted series. A data frame
#'   with no rows signals an error.
#' @param se `TRUE` to report the standard error and an interval as well.
#'   The default is `FALSE`. A moment has no standard error, and `se = TRUE`
#'   with a moment signals an error. [group_interval()] sets it to `TRUE`, and
#'   [prediction_interval()] does not read it, because its result carries
#'   its own `se` column.
#' @param level The interval's level, `0.95` by default, a single number
#'   strictly between 0 and 1. It is read wherever an interval is returned:
#'   with `se = TRUE`, with [group_interval()] and with
#'   [prediction_interval()], which gives its interval without `se`.
#' @param random How a [modelterms7::random()] term is read:
#'   `"conditional"` (the default) with each group's own estimated effect, so
#'   a level that is absent from the fitting data signals an error; `"zero"`
#'   with every effect at zero, the typical group; `"marginal"` averaged over
#'   the prior that the fit estimated, the population average
#'   \eqn{E_b[h^{-1}(\eta + z^\top b)]}. A single string applies to every such
#'   term; a character vector named by the terms' keys sets the mode term by
#'   term, and the terms that it does not name stay conditional. See the section
#'   on new groups.
#' @param interval An interval object that selects the uncertainty that the
#'   interval describes. [confidence_interval()] (the default) gives the
#'   interval of the parameter, the predictor or the moment from the
#'   uncertainty of the estimates; [group_interval()] gives the interval of the
#'   parameter of a new group; [prediction_interval()] gives the interval of a
#'   new observation of the response, and takes the predictive distribution as
#'   its argument (`predictive = pred_averaged()`, [pred_plugin()] or
#'   [pred_bootstrap()]). A character string is rejected with an error. See
#'   [predict_intervals] and [predictive_distributions].
#' @param ... Only `type`, the variance type of [vcov.StatmodFit()]: the
#'   Bayesian variance, the frequentist one or the unconditional one, which
#'   also accounts for the estimated smoothing parameters. It is read where
#'   `se` is `TRUE`, with [group_interval()], and with
#'   [prediction_interval()] holding [pred_averaged()], and it signals an
#'   error elsewhere. Any other argument signals an error that names it. The
#'   unconditional variance differs most from the others in a band around a
#'   penalized term, since a smoothing parameter moves the fitted values.
#' @return With `se = FALSE`, a numeric vector with one value per row of
#'   `newdata` (of the fitting data where `newdata` is `NULL`) when `what`
#'   names one quantity, and a named list of such vectors for `"parameter"`
#'   and `"link"`.
#'
#'   With `se = TRUE`, a data frame with columns `fit`, `se`, `lower` and
#'   `upper` in place of each vector, or, for `"parameter"` and `"link"`, a
#'   named list of such data frames, one per parameter. `se` is `NA` for an
#'   observation whose predictor reads a coefficient that has no variance.
#'
#'   With [prediction_interval()], a data frame with columns `fit` (the
#'   predictive median), `se` (the predictive standard deviation, `NA` where
#'   it does not exist), `lower` and `upper`, all `NA` for a row whose
#'   predictor reads a coefficient that has no variance.
#' @seealso [fitted.StatmodFit()] for one parameter's fitted values,
#'   [residuals.StatmodFit()] for the matched diagnostic,
#'   [vcov.StatmodFit()] for the variance from which the standard errors are
#'   computed.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#'
#' # One parameter and one of the family's moments are predicted.
#' head(predict(fit, "mu"))
#' head(predict(fit, "variance"))
#'
#' # The predictor of a parameter is returned instead of its value.
#' head(predict(fit, "link:sigma"))
#'
#' # For a Gaussian the mean is mu and the variance is sigma squared, and the
#' # moments of the fit return these values.
#' all.equal(predict(fit, "mean"), predict(fit, "mu"))
#' all.equal(predict(fit, "variance"), predict(fit, "sigma")^2)
#'
#' # The result carries an interval as well.
#' head(predict(fit, "mu", se = TRUE))
#'
#' # Every parameter is returned at once, on either scale.
#' str(predict(fit, "parameter"))
#'
#' # For a new group, the prediction is the typical one or the population
#' # average.
#' gg <- data.frame(g = factor(rep(letters[1:12], each = 10)),
#'                  x = rnorm(120))
#' gg$y <- rpois(120, exp(0.5 + 0.3 * gg$x + rnorm(12, sd = 0.6)[gg$g]))
#' fg <- statmod(y ~ x + random(~ 1 | g), distributions7::poisson_distrib(),
#'               gg)
#' new <- data.frame(g = "new", x = 0)
#' predict(fg, "mu", new, random = "zero")
#' predict(fg, "mu", new, random = "marginal")
#'
#' # The calls return the interval of a new group's mean and of a new count
#' # from it.
#' predict(fg, what = "mu", newdata = new, random = "zero",
#'         interval = group_interval())
#' predict(fg, what = "response", newdata = new, random = "zero",
#'         interval = prediction_interval(predictive = pred_plugin()))
#' predict(fg, what = "response", newdata = new, random = "zero",
#'         interval = prediction_interval())
#'
#' # A name that the family does not have is rejected, and the message lists
#' # what is available.
#' try(predict(fit, "median"))
#' @references Harris, I. R. (1989). Predictive fit for natural exponential
#'   families. *Biometrika*, 76, 675--684.
#' @keywords internal
predict.StatmodFit <- function(object, what = "parameter", newdata = NULL,
                               se = FALSE, level = 0.95,
                               random = "conditional",
                               interval = confidence_interval(), ...) {
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
  ps <- object@spec@distrib@params
  known <- what %in% c(ps, "parameter", "link", "response",
                       names(predict_moments())) ||
    (startsWith(what, "link:") && substring(what, 6L) %in% ps)
  if (!known) {
    stop(unknown_what(if (startsWith(what, "link:")) substring(what, 6L) else
      what, ps), call. = FALSE)
  }
  ispec <- predict_interval_spec(interval)
  interval <- ispec$interval
  predictive <- ispec$predictive
  n_boot <- ispec$n_boot
  boot_refit <- ispec$boot_refit
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  # the dots carry the variance type to vcov(), and nothing else: an
  # argument accepted here and read by nobody would change nothing
  dots <- list(...)
  dn <- names(dots)
  if (length(dots) && (is.null(dn) || any(!nzchar(dn)) ||
                       any(dn != "type"))) {
    bad <- if (is.null(dn)) rep("", length(dots)) else dn
    bad <- bad[bad != "type"]
    stop(sprintf(paste0("predict() takes no argument %s. The only further ",
                        "argument is 'type',\n  the variance type of the ",
                        "standard error."),
                 paste0("'", ifelse(nzchar(bad), bad, "<unnamed>"), "'",
                        collapse = ", ")), call. = FALSE)
  }
  if (length(dots) &&
      ((!isTRUE(se) && identical(interval, "confidence")) ||
       (identical(interval, "prediction") &&
        !identical(predictive, "averaged")))) {
    stop(paste0("'type' is the variance type of the standard error, and is ",
                "read only with\n  se = TRUE, with group_interval(), or ",
                "with prediction_interval(predictive = pred_averaged())."),
         call. = FALSE)
  }
  rm <- random_modes(object@spec, random)
  aside <- rm[rm$mode != "conditional", , drop = FALSE]
  nest <- attr(rm, "nested")
  nest <- nest[nest$mode != "conditional", , drop = FALSE]
  if (nrow(nest) && !identical(interval, "confidence")) {
    stop(sprintf(paste0("%s is not available for a random ",
                        "effect written inside a\n  term's subformula (%s). ",
                        "random = \"zero\" and \"marginal\" are, with the ",
                        "default\n  confidence_interval()."),
                 interval_call(interval), nest$key[1L]),
         call. = FALSE)
  }
  if (identical(what, "response") && !identical(interval, "prediction")) {
    stop(paste0("The response is predicted with interval = prediction_interval(): ",
                "it has a\n  distribution, not a single value ",
                "to report."), call. = FALSE)
  }
  if (identical(interval, "prediction") && !identical(what, "response")) {
    stop(paste0("prediction_interval() is an interval for a new ",
                "observation of the\n  response: use what = ",
                "\"response\". A parameter's interval for a new group\n",
                "  is group_interval()."), call. = FALSE)
  }
  if (identical(interval, "group") && !nrow(aside)) {
    stop(paste0("group_interval() is the interval of a new group, which ",
                "needs a random\n  effect set aside: random = \"zero\" or ",
                "\"marginal\"."), call. = FALSE)
  }
  if (!identical(interval, "confidence")) {
    if (S7::S7_inherits(object@spec@distrib,
                        distributions7::multivariate_distrib)) {
      stop(interval_call(interval), " is not available for a ",
           "multivariate family.", call. = FALSE)
    }
    if (length(statmod_structural(object@spec))) {
      stop(interval_call(interval), " is not available beside a ",
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
  if (is.data.frame(newdata) && !nrow(newdata)) {
    stop("'newdata' has no rows.", call. = FALSE)
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
    } else statmod_design(spec, unseen,
                          predict_unread(object, what, newdata, interval)),
    error = function(e) {
      if (grepl("was not present at build time", conditionMessage(e),
                fixed = TRUE)) {
        stop(paste0(conditionMessage(e), "\n  A group that was not in the ",
                    "fitting data has no effect of its own. Where the model has no\n",
                    "  structural term, predict(random = \"zero\") reads it as ",
                    "the typical group and\n  random = \"marginal\" averages ",
                    "over the prior."),
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
  reached <- if (cont) unlist(lapply(attr(design, "structural"),
                                      function(u) u$param)) else NULL
  returned <- if (what %in% ps) what else if (startsWith(what, "link:"))
    substring(what, 6L) else ps
  if (cont && isTRUE(se) && any(reached %in% returned)) {
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
    # standard deviation are those of h^-1(eta*), see group_interval_ends()
    if (identical(interval, "group")) {
      su <- group_interval_ends(object, spec, design, aside, ep, su, level, ...)
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
                    "both averaged over. Each kind is available on its own."),
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
                        "  available: use one of %s."),
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

#' The Call That Names an Interval Kind
#'
#' @param kind `"confidence"`, `"group"` or `"prediction"`.
#' @return The constructor call as a string.
#' @keywords internal
interval_call <- function(kind) {
  switch(kind, confidence = "confidence_interval()", group = "group_interval()",
         prediction = "prediction_interval()")
}

#' The Equations That a Prediction Does Not Read
#'
#' @description
#' Returns the distribution parameters whose equations
#' [predict.StatmodFit()] can leave unbuilt. These are the parameters other
#' than `what` whose variables `newdata` does not carry, where `what` names
#' one parameter.
#'
#' @details
#' A prediction of one parameter reads only that parameter's equation, so
#' the covariates of the other equations need not be in `newdata`: on
#' `accel ~ s(times) | sigma ~ z`, `predict(fit, what = "mu", newdata)` does
#' not need `z`. Every equation is built where `what` is not the name of a
#' parameter (`"parameter"`, `"link"`, `"link:<parameter>"` or a moment),
#' where `interval` is not [confidence_interval()], or where the model carries a
#' structural term, whose recursion reads the predictors of every equation.
#'
#' @param object A [StatmodFit()].
#' @param what,newdata,interval As in [predict.StatmodFit()].
#'
#' @return A character vector of parameter names, possibly empty.
#'
#' @keywords internal
predict_unread <- function(object, what, newdata, interval = "confidence") {
  params <- object@spec@distrib@params
  if (is.null(newdata) || !(what %in% params) ||
      !identical(interval, "confidence") ||
      length(statmod_structural(object@spec))) {
    return(character(0))
  }
  eq <- statmod_equations(object@spec@formula, params)$equations
  # a variable of an unread equation that is not a column of newdata is not
  # looked up anywhere else: the equation is not read, so leaving it unbuilt
  # costs nothing, where reading the formula's environment could find a
  # vector of the fitting rows
  missing_vars <- function(q) {
    v <- all.vars(q[[length(q)]])
    v[!(v %in% names(newdata))]
  }
  setdiff(params[vapply(params, function(p)
    !is.null(eq[[p]]) && length(missing_vars(eq[[p]])) > 0L, logical(1))], what)
}



#' The Message for an Unrecognized Prediction Target
#'
#' @description
#' Builds the error that [predict.StatmodFit()] signals when `what` names
#' nothing that it can compute, listing this family's own parameter names,
#' the moments and the collective targets.
#'
#' @details
#' The family's parameters are listed by name rather than described, since
#' they differ from family to family. A data frame passed where `what`
#' belongs, which is the mistake that the argument order of
#' [stats::predict()] invites, is recognized and named separately.
#'
#' @param what The name that was passed, for the message.
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
#' Returns one distribution parameter's fitted values, as a vector with one
#' value per row of the fitting data.
#' @details
#' The result is one vector, never the whole set, so that this and
#' [residuals.StatmodFit()] are a matched pair and a diagnostic drawn from
#' them needs no unpacking.
#'
#' The default is the first parameter, not the mean. A family may have no
#' mean (a Cauchy has none), whereas every family has a first parameter.
#'
#' The whole set at once is `predict(fit, "parameter")`, and the mean, where
#' it exists, is `predict(fit, "mean")`.
#' @param object A [StatmodFit()].
#' @param what Which distribution parameter, a string naming one of the
#'   family's, or `NULL` for the first.
#' @param ... Not used. Any further argument signals an error, since the
#'   values at new data are given by `predict(fit, what, newdata = new_data)`.
#' @return A numeric vector with one value per row of the fitting data, that
#'   parameter's fitted values on its own scale.
#' @seealso [predict.StatmodFit()] for the other quantities and for new data,
#'   [residuals.StatmodFit()] for the matched diagnostic.
#' @keywords internal
fitted.StatmodFit <- function(object, what = NULL, ...) {
  if (...length()) {
    stop(paste0("fitted() takes 'what' alone. The values at new data are ",
                "predict(fit, what,\n  newdata = new_data)."), call. = FALSE)
  }
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
#' Computes the standard error of each equation's linear predictor, and of
#' the parameter that the predictor gives, with an interval.
#'
#' @details
#' # The delta method, equation by equation
#'
#' An equation's predictor is \eqn{\eta_{ip} = x_{ip}'\beta_p}, so its
#' variance is \eqn{x_{ip}' V_{pp} x_{ip}}, with \eqn{V} the variance of the
#' coefficients as estimated: the coordinates, since the design is written in
#' them, not the quantities that [coef.StatmodFit()] reports by default.
#'
#' The equations do not mix. One equation's predictor reads only the
#' coefficients of that equation, whatever the covariance between the
#' blocks. A
#' predictor that reads no coefficient (an equation written as `~ 0`, or an
#' equation with only an offset) is a known constant, with a standard error
#' of zero.
#'
#' # A term whose block moves needs no special case
#'
#' Its block is the Jacobian \eqn{\partial\eta/\partial\beta} by
#' construction, which is why a linear fit on it is a Gauss-Newton step.
#' The row is therefore already the derivative, and the delta method is exact
#' to first order. That covers [modelterms7::seg()], [modelterms7::nl()] and a
#' smoothed [modelterms7::jump()] or [modelterms7::jseg()], including the
#' parameters of a nonlinear term developed over covariates. A sharp
#' [modelterms7::jump()] or [modelterms7::jseg()] is fitted as an ordinary
#' block at its settled break-point, so the standard error of a prediction
#' that reads it is finite and is conditional on that position.
#'
#' # The interval
#'
#' The interval is built on the scale on which the equation is written and
#' mapped back through the link, as every interval in this toolkit is. A
#' scale parameter under a log link, for example, cannot have a negative
#' lower end.
#'
#' # A coefficient with no variance carries none forward
#'
#' An entry of the variance that is not finite blanks the standard error of
#' the rows whose design reaches it, and does not blank that of any other
#' row. An aliased
#' coefficient is such a case. The position of a sharp break-point term is
#' held and has no variance either, but its column of the design is zero, so
#' it blanks no row.
#'
#' A coefficient that a lasso, an elastic net, a SCAD or an MCP holds at
#' exactly zero is treated as fixed: the model is selected without it, its
#' row of the variance is `NA`, and its column of the design is left out of
#' the quadratic form, so it contributes nothing. Otherwise a prediction at
#' new data would read every such coefficient against a non-zero covariate
#' and return `NA` on every row.
#'
#' @param object A fitted model.
#' @param spec The specification the prediction is made under.
#' @param design Its design.
#' @param ep The predictors and parameters, as [statmod_eta()]
#'   returns them.
#' @param level The interval's level.
#' @param ... Passed to [vcov()], which also receives `readable = FALSE` from
#'   this function.
#'
#' @return A named list, one entry per distribution parameter, each a data
#'   frame with the columns `eta`, `se_eta`, `eta_lower` and `eta_upper` on
#'   the link scale and `fit`, `se`, `lower` and `upper` on the parameter
#'   scale.
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
  # the coefficients a kinked penalty holds at zero, by key: fixed at zero,
  # they carry no variance into the prediction (see the details)
  stacked <- unlist(lapply(spec@distrib@params, function(p) {
    if (design[[p]]$npar) paste(p, design[[p]]$coef_names, sep = ":") else
      character(0)
  }))
  zk <- tryCatch(zero_kinked(spec, design, object@coefficients),
                 error = function(e) integer(0))
  held0 <- stacked[zk]
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
    # a predictor that reads no coefficient is a known constant (an equation
    # written as `~ 0`, an offset alone): its variance is zero
    if (!length(key)) v <- rep(0, n)
    if (length(key) && all(key %in% rownames(V))) {
      keep <- !(key %in% held0)
      Vp <- as.matrix(V[key[keep], key[keep], drop = FALSE])
      if (nrow(X) == n) {
        Xk <- X[, keep, drop = FALSE]
        v <- row_quad(Xk, Vp, Xk)
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


#' The Shape in Which a Prediction With Its Uncertainty Is Returned
#'
#' @description
#' Reduces the rows that [predict_se()] computed to the quantity requested.
#'
#' @details
#' The accepted values of `what` are those of a prediction without
#' uncertainty, minus the moments: a moment's delta method needs its
#' derivative in every parameter, which \pkg{distributions7} does not provide
#' as a generic, so a moment signals an error.
#'
#' @param su The per-parameter tables.
#' @param what The quantity requested, as in [predict.StatmodFit()].
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
  stop(sprintf(paste0("'%s' has no standard error here. A moment's standard ",
                      "error would need its\n  derivative in every parameter, ",
                      "which the distribution does not offer;\n  use a ",
                      "parameter (%s), \"parameter\" or \"link\"."),
               what, paste(params, collapse = ", ")), call. = FALSE)
}


#' The Derivative Row of an Equation Carrying a Structural Term
#'
#' @description
#' Returns the derivative of the predictor of the equation that carries a
#' structural term, one row per observation, in every coordinate of the
#' variance matrix on which the predictor depends: the coefficients of every
#' equation and the term's free parameters on their unconstrained scale. A
#' score-driven term is read through [filter_joint_jacobian()] and a term of
#' the likelihood shape through [loglik_shift_jacobian()].
#'
#' @details
#' A filter's level is a recursion, not a column, so it has no row of a
#' design. Its derivative is the forward Jacobian of the recursion, which
#' [filter_joint_jacobian()] returns. Every equation's coefficients enter it,
#' not only those of the filter's own equation: the scores that drive the
#' recursion are read at the predictors of every equation, so a coefficient
#' of the scale moves the level of a filter in the mean. Leaving those
#' columns out gives a wrong standard error for the filtered mean.
#'
#' A parameter that an intercept in the same equation holds is not estimated
#' and is not in that matrix, so it is not here either.
#'
#' @param spec The specification.
#' @param design Its design.
#' @param ep The predictors, as [statmod_eta()] returns them.
#' @param p The distribution parameter whose equation is being read.
#' @param coef The coefficients at which the predictors were evaluated.
#'
#' @return A list with `X` (the derivative rows) and `key` (the names of
#'   their columns in the variance matrix), or `NULL` where the equation
#'   carries no structural term.
#'
#' @seealso [predict_se()], [filter_joint_jacobian()], [loglik_shift_jacobian()]
#'
#' @keywords internal
structural_se_columns <- function(spec, design, ep, p, coef) {
  rg <- ep$regimes
  if (length(rg) &&
      any(vapply(rg, function(r) identical(r$param, p), logical(1)))) {
    return(loglik_shift_jacobian(spec, design, coef))
  }
  fs <- ep$filters
  in_joint <- any(vapply(fs, function(f) isTRUE(f$joint) && p %in% f$params,
                         logical(1)))
  if (in_joint) {
    fj <- joint_filter_jacobian(spec, design, coef)
    return(list(X = fj$V[[p]], key = fj$key))
  }
  if (!length(fs) ||
      !any(vapply(fs, function(f) identical(f$param, p), logical(1)))) {
    return(NULL)
  }
  fj <- filter_joint_jacobian(spec, design, coef)
  if (is.null(fj)) return(NULL)
  list(X = fj$J, key = fj$key)
}


#' The Derivative Row of an Equation Carrying a Likelihood Mixed Over States
#'
#' @description
#' Returns the derivative of the posterior-weighted predictor of the equation
#' that a structural term of the likelihood shape (a regime term or a
#' marginal break-point term) shifts, one row per observation, in the
#' coefficients of every equation and the term's free parameters on their
#' unconstrained scale.
#'
#' @details
#' The predictor is the static part plus the posterior mean of the term's
#' shift. Its derivative is the equation's design row plus the derivative of
#' that mean, which [modelterms7::term_shift_jacobian()] returns: the
#' posterior moves with every coordinate, through the emissions of every
#' equation and through the term's own parameters, and the shifts move with
#' the term's own parameters. The family's score at each component's
#' predictor is supplied here, as for [statmod_regime_information()].
#'
#' @param spec The specification.
#' @param design Its design.
#' @param coef The coefficients.
#'
#' @return A list with `X` (the derivative rows) and `key` (the names of
#'   their columns in the variance matrix), restricted to the free columns,
#'   or `NULL` where the model carries no such term.
#'
#' @seealso [structural_se_columns()], [statmod_regime_information()]
#'
#' @keywords internal
loglik_shift_jacobian <- function(spec, design, coef) {
  ev <- statmod_eta(spec, design, coef)
  if (!length(ev$regimes)) return(NULL)
  r <- ev$regimes[[1L]]
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  ap <- match(r$param, params)
  npar <- vapply(design[params], function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  nb <- sum(npar)
  n <- spec@n_obs
  zn <- modelterms7::term_params(r$tm)
  np <- length(zn)
  m <- nb + np
  Vs <- lapply(seq_along(params), function(a) {
    M <- matrix(0, n, m)
    if (npar[a] > 0L) {
      M[, offs[a] + seq_len(npar[a])] <- as_dense(design[[params[a]]]$X)
    }
    M
  })
  theta_at <- function(e, i) {
    stats::setNames(lapply(seq_along(params), function(q) {
      x <- if (q == ap) e else rep_len(ev$eta_static[[params[q]]], n)[i]
      linkfunctions7::linkinv(links[[params[q]]], x)
    }), params)
  }
  yof <- function(i) {
    if (is.matrix(spec@response)) spec@response[i, , drop = FALSE]
    else spec@response[i]
  }
  ld <- function(e, i) {
    distributions7::distrib_pdf(spec@distrib, yof(i), theta_at(e, i),
                                log = TRUE, threads = spec@threads)
  }
  gr <- function(e, i) {
    g <- distributions7::distrib_gradient(spec@distrib, yof(i),
                                          theta_at(e, i), scale = "link",
                                          threads = spec@threads)
    vapply(params, function(q) rep_len(g[[q]], length(i)), numeric(length(i)))
  }
  sj <- modelterms7::term_shift_jacobian(r$tm, r$eta_static, spec@response,
                                         logdens = ld, grad = gr, psi = r$psi,
                                         seed = Vs, cols = nb + seq_len(np),
                                         level = ap)
  J <- Vs[[ap]] + sj$jacobian
  ckey <- unlist(lapply(params, function(q) {
    if (design[[q]]$npar) paste(q, design[[q]]$coef_names, sep = ":") else
      character(0)
  }), use.names = FALSE)
  lb <- tryCatch(r$tm@label, error = function(e) "")
  tkey <- paste(r$param, if (length(lb) == 1L && nzchar(lb))
    paste(lb, zn, sep = ".") else zn, sep = ":")
  held <- statmod_structural_state(design)$held[[r$term]]
  keep <- if (length(held)) setdiff(seq_len(m), nb + match(held, zn)) else
    seq_len(m)
  list(X = J[, keep, drop = FALSE], key = c(ckey, tkey)[keep])
}


#' The Forward Jacobian of a Filter in Every Coordinate
#'
#' @description
#' Returns the derivative of the predictor that a score-driven term produces,
#' in the coefficients of every equation followed by the term's free
#' parameters on their unconstrained scale, at each observation, together
#' with the static rows of every equation in the same columns.
#'
#' @details
#' The rows are those from which the joint information is assembled
#' ([statmod_full_information()]): each equation's design placed in its own
#' columns, and for the filter's equation the Jacobian that
#' [modelterms7::term_curvature()] propagates beside the state. A model
#' carries at most one such filter, so every other equation is static and
#' its design is its derivative; a joint filter is read by
#' [joint_filter_jacobian()].
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#' @param coef The coefficients.
#'
#' @return `NULL` where the model carries no filter, otherwise a list with
#'   `J` (the filter equation's rows), `V` (every equation's rows, of which
#'   the filter's is `J`), `H` (the family's second derivatives on the link
#'   scale), `Hd` (the derivatives of the filter's driving quantity, which
#'   are `H` where the filter is not scaled), `ap` (the index of the filter's
#'   parameter), `nb` (the number of coefficient columns), `key` (the names of
#'   the columns in the variance matrix) and `free` (the term's free
#'   parameters). All columns are restricted to the free ones.
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
#' Computes the linear predictors and the parameters at new rows for a model
#' that carries a structural term, whose own recursion continues from where
#' the fitting data left it.
#'
#' @details
#' The static part is the ordinary assembly: each term's block is reapplied
#' at the new rows, and the offsets are re-evaluated there. The structural
#' term cannot be reapplied, because its contribution at one row is the state
#' that a recursion has reached over every row before it. The fit is
#' therefore run once at the observed rows to recover that state, which
#' consists of the level and the driving quantity at each of them (the level
#' is the difference between the filtered predictor and the static one), and
#' the term continues from there through [modelterms7::term_continue()].
#'
#' Only rows without the response are continued, and they must come after
#' the observed series. Rows that carry it are a series of their own, rebuilt
#' by [statmod_respec()] and read as the fitting rows are, which is why
#' `predict(fit, newdata = <the fitting data>)` returns the fitted values.
#'
#' A term whose contribution is a likelihood mixed over latent states is
#' rejected, because such a term reports a posterior over states at an
#' observed row, and past the data this posterior is a predictive
#' distribution and not a single value.
#'
#' With `deriv = TRUE` the continuation is differentiated as well, in the
#' coefficients of every equation and the term's free parameters
#' ([continued_deriv_inputs()]), and a standard error of the forecast reads
#' this derivative.
#'
#' @param fit The fitted model.
#' @param spec The specification at the new data.
#' @param design Its design.
#' @param deriv Whether to return the forecast's derivative rows.
#'
#' @return A list shaped as [statmod_eta()]'s, without the memoized filter
#'   objects and with `eta_static`, the static predictor at the new rows
#'   before the continued levels are added, plus `cont_cols`: for each
#'   equation whose filter was continued with `deriv = TRUE`, a list with `X`
#'   (the derivative of the predictor at the new rows in the coefficients of
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

  # the static part at the new rows, kept apart as `eta_static` before the
  # continued levels are added
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

  eta_static <- eta
  # the state the observed part ended at
  ospec <- fit@spec
  odesign <- statmod_design(ospec)
  oep <- statmod_eta(ospec, odesign, coef)
  ost <- statmod_structural_state(odesign)
  for (f in oep$filters) {
    tm <- ospec@terms[[f$param]][[f$term]]
    psi <- structural_psi(tm, ost$zeta[[f$term]])
    if (isTRUE(f$joint)) {
      # the levels and the driving vector at the observed rows, from which
      # the joint recursion continues
      f_past <- f$eta - f$eta_static
      s_past <- t(vapply(seq_len(nrow(f$eta)), function(i)
        as.numeric(f$cb$score(as.numeric(f$eta[i, ]), i)),
        numeric(length(f$params))))
      s_past <- matrix(s_past, nrow(f$eta), length(f$params))
      dv <- if (isTRUE(deriv)) {
        joint_continued_deriv_inputs(ospec, odesign, coef, ost)
      } else NULL
      cont <- modelterms7::term_continue(tm, psi, f_past, s_past, nd,
                                         deriv = dv$deriv)
      for (k in seq_along(f$params)) {
        pk <- f$params[k]
        eta[[pk]] <- eta[[pk]] + as.numeric(cont[, k])
        theta[[pk]] <- linkfunctions7::linkinv(links[[pk]], eta[[pk]])
        if (!is.null(dv)) {
          jac <- attr(cont, "jacobian")[[k]]
          if (d2[[pk]]$npar) {
            jac[, dv$cols[[k]]] <- jac[, dv$cols[[k]]] + as.matrix(d2[[pk]]$X)
          }
          cont_cols[[pk]] <- list(X = jac, key = dv$key)
        }
      }
      next
    }
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
       eta_static = eta_static, cont_cols = cont_cols)
}


#' The Inputs to the Derivative of a Continued Filter
#'
#' @description
#' Returns the derivatives of the level and of the score at the observed
#' rows, and of the term's parameters, in the coordinates in which the
#' variance matrix of the fit is written: the coefficients of every equation
#' followed by the term's free parameters on their unconstrained scale.
#'
#' @details
#' The coordinates are those that [structural_se_columns()] uses for a
#' prediction at the observed rows, so a forecast and a fitted value are
#' read against the same rows and columns of the variance. The derivative of
#' the filtered predictor at an observed row is the forward Jacobian of
#' [filter_joint_jacobian()]; the level's derivative is that minus the
#' equation's design row. The quantity that drives the recursion depends on
#' the predictors of every equation, so its derivative is
#' \eqn{\sum_b \ell_{pb,t} V_{b,t}}, with \eqn{V_{b,t}} the derivative row
#' of equation \eqn{b} and \eqn{\ell_{pb,t}} the derivative of the driving
#' quantity with respect to the predictor of equation \eqn{b}. Where the
#' filter is not scaled, the driving quantity is the score and
#' \eqn{\ell_{pb,t}} is the family's second derivative on the link scale;
#' where it is scaled, \eqn{\ell_{pb,t}} is the derivative that
#' `filter_driving()` returns. The derivative of a parameter in its own
#' unconstrained coordinate is that of its link, and for a coefficient of a
#' development it is one.
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
#' Returns `TRUE` where every row has the response, `FALSE` where every row
#' lacks it, and signals an error where some rows have it and others do not.
#'
#' @details
#' The distinction separates a re-reading of a model on another series from a
#' continuation of the series to which the model was fitted, and it has to be
#' all-or-nothing. The recursion of a filter at one row reads the rows before
#' it, so a frame that carries the response on some rows only describes
#' neither operation, and any answer would require choosing one of the two
#' readings for the caller.
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
       "series. Pass rows with the response or rows without it.",
       call. = FALSE)
}
