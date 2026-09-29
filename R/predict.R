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
#' reached, so new rows continue the series instead of being read on their
#' own. Each row is placed by its own time within its own group, and must
#' come after every observed time of that group; a row falling inside the
#' observed series is refused, since there the response is known and the
#' filter must be run, never continued.
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
#' and there are at most three coordinates in all, and 10000 draws from the
#' prior otherwise, on a seed of their own. A prior that is not Gaussian is
#' integrated only where the inverse link is bounded: a Student t has no
#' moment generating function, so under a log link the average does not
#' exist, and such a request signals an error. A term that shares a covariance
#' with others through a label has no prior of its own and is read at
#' `"zero"` only. The standard error of a marginal parameter is the delta
#' method on the fixed part, conditional on the prior's scale.
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
#'     "zero"` or `"marginal"`. The variance of its predictor is the
#'     estimates' plus the one a new group's effects add,
#'     \eqn{se^2 + z^\top \Sigma_b z}, in the parameter's own equation, so
#'     a random effect on `sigma` widens the interval for `sigma`. On the link
#'     scale this is the standard error glmmTMB reports at a new level. The
#'     interval is the predictor's, carried through the link; the fit stays
#'     what `random` asks for, the typical group's value or the population
#'     average. Implies `se = TRUE`.}
#'   \item{`"prediction"`}{an interval for a new observation of the response,
#'     with `what = "response"`. Every predictor is taken jointly Gaussian,
#'     with the covariance of the estimates plus the one the effects set aside
#'     add -- correlated across equations where a label ties them -- and
#'     \eqn{Y} is the family averaged over it. The interval's ends are
#'     quantiles of that average and the fit is its median, which every family
#'     has, so no location parameter is needed; `se` is the standard deviation
#'     of the average, `NA` where the family's variance does not exist. For a
#'     discrete family the ends are values of the support, so the coverage is
#'     at least the level rather than equal to it. Effects read
#'     `"conditional"` enter with their own estimates, so this is also the
#'     interval of a new observation of a group the fit saw.}
#' }
#' Both intervals for a new group are conditional on the covariance the fit
#' estimated: its uncertainty is not propagated. Measured at 95 per cent over
#' 200 fits of a random intercept at eight observations a group, a new
#' group's predictor is covered 0.898 of the time over 10 groups and 0.943
#' over 40 by `"group"`, and a new observation 0.939 and 0.951 by
#' `"prediction"`, where the confidence interval of the typical group covers
#' 0.45 and 0.26; on a Poisson response the same read 0.909 and 0.937 and,
#' the support being discrete, 0.975 and 0.976. A prior that is not Gaussian
#' has no covariance to read, and
#' a term that shares a label with an effect read with its own estimate, or
#' with one written inside a subformula, is rejected.
#'
#' A forecast reports **no standard error**. `se = TRUE` gives the
#' uncertainty of the parameters, while a forecast carries the uncertainty of
#' the future scores as well, which is the larger part and is no delta
#' method. Reporting the smaller half under the name of the whole would
#' mislead.
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
#' predict(fg, "response", new, random = "zero", interval = "prediction")
#'
#' # A name the family does not have is refused, and the message says what
#' # is available.
#' try(predict(fit, "median"))
#' @keywords internal
predict.StatmodFit <- function(object, what = "parameter", newdata = NULL,
                               se = FALSE, level = 0.95,
                               random = "conditional",
                               interval = c("confidence", "group",
                                            "prediction"), ...) {
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
  rm <- random_modes(object@spec, random)
  aside <- rm[rm$mode != "conditional", , drop = FALSE]
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
  if (nrow(aside) && length(statmod_structural(object@spec))) {
    stop("random = \"zero\" or \"marginal\" is not available beside a ",
         "structural term: its\n  level is a recursion read at the ",
         "predictor, so setting a random effect aside\n  would move the ",
         "recursion as well.", call. = FALSE)
  }
  spec <- spec_at(object, newdata, need_response = FALSE)
  design <- withCallingHandlers(
    statmod_design(spec, if (nrow(aside))
      split(aside$key, aside$param) else NULL),
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
  nodes <- NULL
  if (nrow(mt)) {
    priors <- lapply(seq_len(nrow(mt)), function(i)
      random_prior(spec, design, object, mt$param[i], mt$key[i]))
    nodes <- random_nodes(priors)
  }
  # A STRUCTURAL TERM HAS NO BLOCK TO REAPPLY: its contribution is the state
  # a recursion has reached, so new rows CONTINUE the series rather than
  # being read on their own. Running the ordinary assembly there returned
  # the fitting data's values whatever `newdata` held, which is why the two
  # paths are separated rather than merged.
  cont <- !is.null(newdata) && length(attr(design, "structural"))
  ep <- if (cont) statmod_eta_continued(object, spec, design) else
    statmod_eta(spec, design, object@coefficients)
  if (cont && isTRUE(se)) {
    stop("The standard error of a prediction past the series is not ",
         "reported. What\n  predict(se = TRUE) gives is the uncertainty of ",
         "the parameters; a forecast\n  carries the uncertainty of the ",
         "future scores as well, which is the larger part\n  and is no ",
         "delta method. predict(fit, se = TRUE) at the observed rows is ",
         "exact.", call. = FALSE)
  }
  params <- spec@distrib@params
  if (identical(interval, "prediction")) {
    blocks <- random_blocks(spec, design, object, aside)
    pc <- predictive_cov(object, spec, design, blocks, ...)
    return(predictive_response(spec, ep$eta, pc$fixed + pc$random, level))
  }
  if (isTRUE(se)) {
    su <- predict_se(object, spec, design, ep, level, ...)
    # A NEW GROUP'S PARAMETER varies around the typical group's by what its
    # effects add, z' Sigma_b z in its own equation, on top of the estimates'
    # own uncertainty
    if (identical(interval, "group")) {
      blocks <- random_blocks(spec, design, object, aside)
      add <- predictive_cov(object, spec, design, blocks, fixed = FALSE)$random
      zq <- stats::qnorm((1 + level) / 2)
      for (p in params) {
        a <- add[p, p, ]
        if (!any(a != 0)) next
        s <- sqrt(su[[p]]$se_eta^2 + a)
        g <- spec@distrib@link_params[[p]]
        su[[p]]$se_eta <- s
        su[[p]]$eta_lower <- su[[p]]$eta - zq * s
        su[[p]]$eta_upper <- su[[p]]$eta + zq * s
        ends <- cbind(linkfunctions7::linkinv(g, su[[p]]$eta_lower),
                      linkfunctions7::linkinv(g, su[[p]]$eta_upper))
        su[[p]]$lower <- pmin(ends[, 1L], ends[, 2L])
        su[[p]]$upper <- pmax(ends[, 1L], ends[, 2L])
        su[[p]]$se <- abs(linkfunctions7::dlinkinv(g, su[[p]]$eta)) * s
      }
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
    return(se_answer(su, what, params, spec))
  }
  if (nrow(mt) && !identical(what, "link") && !startsWith(what, "link:")) {
    return(random_marginal(spec, design, ep, mt, nodes, what))
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
    # the level is no column of the design: the filter returns its exact
    # derivative in the term's own parameters at every observation, which is
    # the only place it can be computed. Carried onto the free scale by each
    # parameter's own link -- the same chain the exact gradient uses -- those
    # become further columns of the derivative row, read against the tail of
    # the joint variance.
    fl <- structural_se_columns(spec, design, ep, p, X)
    if (!is.null(fl) && nrow(fl$J) == n) {
      X <- cbind(fl$X, fl$J)
      key <- c(key, fl$key)
    }
    if (length(key) && all(key %in% rownames(V))) {
      Vp <- as.matrix(V[key, key, drop = FALSE])
      if (nrow(X) == n && all(is.finite(Vp))) {
        v <- rowSums((X %*% Vp) * X)
        v[v < 0] <- 0
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


#' The Columns a Structural Term Adds to a Derivative Row
#'
#' @description
#' The derivative of one equation's predictor in the free parameters of the
#' score-driven term sitting in it, one row per observation, on the scale the
#' variance matrix is indexed by.
#'
#' @details
#' A filter's level is a recursion, not a column, so it has no row of a
#' design. What it has is the derivative the recursion propagates beside the
#' state, which [modelterms7::term_filter()] returns on the
#' parameter scale; multiplying by each parameter's own
#' \eqn{h'(\zeta_j)} carries it to the unconstrained scale the joint
#' matrix is written in, which is the chain the exact gradient already uses.
#'
#' A parameter an intercept in the same equation holds is not estimated and
#' is not in that matrix, so it is not here either.
#'
#' The design's own rows are corrected at the same time, through
#' [modelterms7::term_static_deriv()]: a coefficient of this
#' equation moves the level as well as the static part, because the scores
#' driving the recursion are read at the predictor the recursion produces.
#' Measured on a score-driven mean with one covariate beside it, leaving
#' that out understates the standard error by about a quarter.
#'
#' @param spec The specification.
#' @param design Its design.
#' @param ep The predictors, as [statmod_eta()] returns them.
#' @param p The distribution parameter whose equation is being read.
#' @param X The equation's design rows.
#'
#' @return A list with `X`, `J` and `key`, or `NULL`
#'   where the equation carries no filter.
#'
#' @seealso [predict_se()], [statmod_filter_at()]
#'
#' @keywords internal
structural_se_columns <- function(spec, design, ep, p, X) {
  fs <- ep$filters
  if (!length(fs)) return(NULL)
  sst <- statmod_structural_state(design)
  cols <- list()
  keys <- character(0)
  seen <- FALSE
  for (f in fs) {
    if (!identical(f$param, p)) next
    if (ncol(X)) {
      D <- modelterms7::term_static_deriv(f$tm, f$curv, X, f$psi)
      if (!is.null(D) && identical(dim(D), dim(X))) X <- D
    }
    seen <- TRUE
    z <- sst$zeta[[f$term]]
    nm <- names(z)
    hl <- sst$held[[f$term]]
    free <- setdiff(nm, hl)
    if (!length(free)) next
    links <- modelterms7::term_links(f$tm)
    ch <- vapply(free, function(j)
      linkfunctions7::dlinkinv(links[[j]], z[[j]]), numeric(1))
    J <- as.matrix(f$jacobian)[, match(free, nm), drop = FALSE]
    cols[[length(cols) + 1L]] <- J * rep(ch, each = nrow(J))
    lb <- tryCatch(f$tm@label, error = function(e) "")
    nmf <- if (length(lb) == 1L && nzchar(lb)) paste(lb, free, sep = ".") else
      free
    keys <- c(keys, paste(p, nmf, sep = ":"))
  }
  if (!seen) return(NULL)
  J <- if (length(cols)) do.call(cbind, cols) else matrix(0, nrow(X), 0L)
  list(X = X, J = J, key = keys)
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
#' Which of the two a call asks for is decided by the response, not by the
#' times. New rows carrying the response are a re-reading: the filter is run
#' over them from the term's own seed, and that is what a caller means by
#' predicting a model on another series, and is why
#' `predict(fit, newdata = <the fitting data>)` returns the fitted
#' values. New rows without it are a continuation, and must come after the
#' observed series. A frame carrying the response on some rows only is
#' rejected: the two readings differ, and picking one would answer a question
#' that was not asked.
#'
#' A term whose contribution is a likelihood mixed over latent states is
#' rejected: what such a term reports at an observed row is a posterior over
#' states, which past the data is a predictive distribution, no single
#' value.
#'
#' @param fit The fitted model.
#' @param spec The specification at the new data.
#' @param design Its design.
#'
#' @return A list shaped as [statmod_eta()]'s, without the
#'   memoized filter objects.
#'
#' @seealso [predict.StatmodFit()]
#'
#' @keywords internal
statmod_eta_continued <- function(fit, spec, design) {
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
  known <- statmod_response_known(spec@response)
  for (f in oep$filters) {
    tm <- ospec@terms[[f$param]][[f$term]]
    psi <- structural_psi(tm, ost$zeta[[f$term]])
    p <- f$param
    if (identical(known, TRUE)) {
      # A RE-READING: the response is there, so the filter is run over these
      # rows from the term's own seed rather than continued. The term is
      # rebuilt on them, which is what gives its recursion their own times
      # and groups; every other block still goes through its blueprint.
      tmn <- modelterms7::term_build(tm, nd)
      cb <- structural_callbacks(spec, theta, p)
      out <- modelterms7::term_filter(tmn, eta[[p]], spec@response,
                                      cb$score, cb$curvature, psi,
                                      fast = cb$fast, threads = spec@threads)
      eta[[p]] <- out$eta
    } else {
      f_past <- as.numeric(f$eta) - as.numeric(f$eta_static)
      # the driving quantity, read at the predictor the recursion produced --
      # the same callback the filter itself was handed
      s_past <- vapply(seq_along(f_past),
                       function(i) f$cb$score(f$eta[[i]], i), numeric(1))
      eta[[p]] <- eta[[p]] +
        modelterms7::term_continue(tm, psi, f_past, s_past, nd)
    }
    theta[[p]] <- linkfunctions7::linkinv(links[[p]], eta[[p]])
  }
  list(eta = eta, theta = theta, filters = list(), regimes = list(),
       eta_static = eta)
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
