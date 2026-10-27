#' @include statmod.R
NULL

#' @title S7 Class for the Interval of a Prediction
#'
#' @description
#' This class is the abstract parent of the interval objects that the
#' `interval` argument of [predict.StatmodFit()] takes: [confidence_interval()],
#' [group_interval()] and [prediction_interval()]. An interval object selects
#' the uncertainty that the interval describes.
#'
#' @details
#' The class is abstract and cannot be instantiated. Calling
#' `interval_type()` signals an error. The class makes a membership test
#' possible with `S7::S7_inherits()`. The three shipped intervals are compared
#' in [predict_intervals].
#'
#' @param label A short name, printed by [print.interval_type()]. The
#'   constructor of each subclass sets it, since the class itself cannot be
#'   built.
#'
#' @return Nothing: the class is abstract and the constructor signals an
#'   error. Used as a class object it has one property, `label`, which every
#'   subclass inherits.
#'
#' @seealso [predict_intervals] for the three intervals side by side,
#'   [predict.StatmodFit()] for the consumer.
#'
#' @examples
#' S7::S7_inherits(confidence_interval(), interval_type)
#' S7::S7_inherits(prediction_interval(), interval_type)
#' try(interval_type("my own"))
#'
#' @export
interval_type <- S7::new_class(
  "interval_type",
  properties = list(label = S7::class_character),
  abstract = TRUE
)

#' @title S7 Class for the Predictive Distribution of a New Response
#'
#' @description
#' This class is the abstract parent of the objects that [prediction_interval()]
#' takes in its `predictive` argument: [pred_averaged()], [pred_plugin()] and
#' [pred_bootstrap()]. A predictive object selects how the estimation
#' uncertainty enters the distribution of a new response.
#'
#' @details
#' The class is abstract and cannot be instantiated. Calling
#' `predictive_type()` signals an error. The three shipped distributions are
#' compared in [predictive_distributions].
#'
#' @param label A short name, printed by [print.predictive_type()]. The
#'   constructor of each subclass sets it, since the class itself cannot be
#'   built.
#'
#' @return Nothing: the class is abstract and the constructor signals an
#'   error. Used as a class object it has one property, `label`, which every
#'   subclass inherits.
#'
#' @seealso [predictive_distributions] for the three side by side,
#'   [prediction_interval()] for the consumer.
#'
#' @examples
#' S7::S7_inherits(pred_averaged(), predictive_type)
#' S7::S7_inherits(pred_bootstrap(n = 50), predictive_type)
#' try(predictive_type("my own"))
#'
#' @export
predictive_type <- S7::new_class(
  "predictive_type",
  properties = list(label = S7::class_character),
  abstract = TRUE
)

#' @title S7 Classes for the Shipped Intervals and Predictive Distributions
#' @description
#' These are the classes that [confidence_interval()], [group_interval()],
#' [prediction_interval()], [pred_averaged()], [pred_plugin()] and
#' [pred_bootstrap()] instantiate. `ConfidenceInterval`, `GroupInterval` and
#' `PredictionInterval` inherit from [interval_type()]; `PredAveraged`,
#' `PredPlugin` and `PredBootstrap` inherit from [predictive_type()].
#' `PredictionInterval` adds the property `predictive`, and `PredBootstrap`
#' adds `n` and `refit`. Use the constructors, which validate; these raw ones
#' do not.
#' @param label A single string, the name that the object prints.
#' @param predictive A [predictive_type()]. Only `PredictionInterval` carries
#'   this property.
#' @param n,refit The number of replicas and the refit mode. Only
#'   `PredBootstrap` carries these properties.
#' @return An S7 object inheriting from [interval_type()] or from
#'   [predictive_type()], with the properties listed above.
#' @seealso [predict_intervals], [predictive_distributions] for the
#'   constructors side by side.
#' @name ConfidenceInterval-class
#' @aliases ConfidenceInterval GroupInterval PredictionInterval PredAveraged
#'   PredPlugin PredBootstrap
#' @keywords internal
ConfidenceInterval <- S7::new_class("ConfidenceInterval", parent = interval_type)

#' @rdname ConfidenceInterval-class
#' @keywords internal
GroupInterval <- S7::new_class("GroupInterval", parent = interval_type)

#' @rdname ConfidenceInterval-class
#' @keywords internal
PredictionInterval <- S7::new_class("PredictionInterval", parent = interval_type,
  properties = list(predictive = predictive_type))

#' @rdname ConfidenceInterval-class
#' @keywords internal
PredAveraged <- S7::new_class("PredAveraged", parent = predictive_type)

#' @rdname ConfidenceInterval-class
#' @keywords internal
PredPlugin <- S7::new_class("PredPlugin", parent = predictive_type)

#' @rdname ConfidenceInterval-class
#' @keywords internal
PredBootstrap <- S7::new_class("PredBootstrap", parent = predictive_type,
  properties = list(n = S7::class_integer, refit = S7::class_character))


#' @title Intervals of a Prediction
#' @name predict_intervals
#' @aliases predict_intervals
#'
#' @description
#' The `interval` argument of [predict.StatmodFit()] takes an interval object
#' that selects the uncertainty that the interval describes. There are three,
#' and they differ in what the interval is about.
#'
#' \describe{
#'   \item{[confidence_interval()]}{the default, with `se = TRUE`: the
#'     interval of a parameter, a predictor or a moment, which carries only the
#'     uncertainty of the estimates. It costs one evaluation of the delta
#'     method.}
#'   \item{[group_interval()]}{the interval of the parameter of a new group,
#'     with `random = "zero"` or `"marginal"`. It adds the variation of the
#'     effects of a new group to the uncertainty of the estimates. It costs
#'     the delta method plus a quantile of the distribution of the
#'     predictor, by closed form under Gaussian priors and by quadrature or
#'     Monte Carlo otherwise.}
#'   \item{[prediction_interval()]}{the interval of a new observation of the
#'     response, with `what = "response"`. It adds the variation of the
#'     family itself. The way that the estimation uncertainty enters is set
#'     by its argument `predictive`, described in [predictive_distributions].
#'     The cost is that of the predictive distribution chosen.}
#' }
#'
#' @details
#' # Confidence interval
#'
#' With the default, an effect set aside by `random` gives the interval of the
#' typical group's value, and the value of a new group falls outside that
#' interval far more often than the level implies. The interval is built on
#' the scale of the equation and mapped back through the link.
#'
#' # Group interval
#'
#' A new group's parameter is \eqn{\theta^* = h^{-1}(\eta^*)}, with
#' \eqn{\eta^*} the typical group's predictor plus its estimation error plus
#' \eqn{z^\top b}, the effects of a new group drawn from the prior, in the
#' parameter's own equation, so a random effect on `sigma` widens the interval
#' for `sigma`. Under Gaussian priors the variance of \eqn{\eta^*} is
#' \eqn{se^2 + z^\top \Sigma_b z}, the standard error that glmmTMB reports at
#' a new level on the link scale. The ends of the interval are quantiles of
#' \eqn{\eta^*} carried through the inverse link, which is monotone, so they
#' are the quantiles of \eqn{\theta^*} under any link; under a prior that is
#' not Gaussian they are the quantiles of the mixture of
#' [predictive_mixture()]. Where the interval of \eqn{\eta^*} leaves the
#' domain of the link (a square root link and a negative predictor), the ends
#' are `NA` at that row, with a warning. The standard error is the standard
#' deviation of \eqn{\theta^*}, not the delta method, and is `NA` where it is
#' not known to exist. The fit is the typical group's value under
#' `random = "zero"` and the population average under `"marginal"`. The
#' interval implies `se = TRUE`.
#'
#' # Prediction interval
#'
#' The family is averaged over the predictors of every parameter, which are
#' taken to be jointly Gaussian (correlated across equations where a label
#' ties them). The ends of the interval are quantiles of that average, and the
#' fit is its median, which every family has, so the family does not need a
#' location parameter. `se` is the standard deviation of the average, and is
#' `NA` where the family's variance does not exist. For a discrete family the
#' ends are values of the support, so the coverage is at least the level
#' rather than equal to it. Effects read `"conditional"` enter with their own
#' estimates, so this is also the interval of a new observation of a group
#' that is present in the fitting data.
#'
#' # Restrictions
#'
#' [group_interval()] and [prediction_interval()] are rejected for a
#' multivariate family and for a model with a structural term. The structural
#' term also rejects `random = "zero"` and `random = "marginal"`: its level is
#' a recursion read at the predictor, so setting a random effect aside would
#' move the recursion as well.
#'
#' Under a prior that is not Gaussian, both [group_interval()] and
#' [prediction_interval()] read the quantiles of the prior: a Student t prior
#' as a scale mixture of Gaussians; a single univariate prior of another
#' family, where it is the only prior that is not Gaussian, on quantile nodes,
#' which do not depend on the seed; and any other combination of priors by
#' Monte Carlo, which does (see [predictive_mixture()]). A term that shares a
#' label with an effect read with its own estimate, or with one written inside
#' a subformula, is rejected.
#'
#' The level of every interval is the argument `level` of
#' [predict.StatmodFit()].
#'
#' @return Not applicable: this topic describes the constructors. Each returns
#'   an S7 object that inherits from [interval_type()].
#'
#' @seealso [confidence_interval()], [group_interval()],
#'   [prediction_interval()], [predictive_distributions],
#'   [predict.StatmodFit()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' nd <- data.frame(x = c(0.2, 0.8))
#' predict(fit, what = "mu", newdata = nd, se = TRUE,
#'         interval = confidence_interval())
#' predict(fit, what = "response", newdata = nd,
#'         interval = prediction_interval())
NULL

#' @title Predictive Distributions of a New Response
#' @name predictive_distributions
#' @aliases predictive_distributions
#'
#' @description
#' The argument `predictive` of [prediction_interval()] takes an object that
#' selects how the estimation uncertainty enters the distribution of a new
#' response. There are three.
#'
#' \describe{
#'   \item{[pred_averaged()]}{the default: the family averaged over the
#'     approximate sampling distribution of the estimates. It costs one
#'     evaluation of the covariance of the estimates.}
#'   \item{[pred_plugin()]}{the family at the estimates, ignoring their
#'     error. It costs one evaluation of the predictors.}
#'   \item{[pred_bootstrap()]}{the parametric bootstrap predictive
#'     distribution, by simulation and refit of `n` responses. It costs `n`
#'     refits, each as expensive as a fit with the hyperparameters held
#'     (`refit = "coefficients"`) or as a whole fit (`refit = "full"`).}
#' }
#'
#' @details
#' A new observation varies for three reasons: the family itself, the effects
#' of a new group where `random` sets them aside, and the error in the
#' estimates. The first two are random in the model and always enter. The
#' third enters according to the predictive object, for every parameter alike,
#' whether the family has one parameter (a Poisson) or several, and whether or
#' not a parameter has covariates.
#'
#' # Averaged
#'
#' The family is averaged over the predictors at the fit's estimates with the
#' covariance of the estimates from [vcov.StatmodFit()] plus the covariance
#' that the effects add. This is the normal approximation to the parametric
#' bootstrap predictive distribution of Harris (1989): the estimates are
#' averaged over their approximate sampling distribution, centered at the
#' fit. For a Gaussian with a constant \eqn{\sigma = e^{\gamma_0}},
#' \eqn{\hat\gamma_0} of variance \eqn{v}, the variance of a new observation
#' is \eqn{\mathrm{se}_0^2 + z^\top\Sigma_b z + \hat\sigma^2 e^{2v}}, where
#' the last term is the average of \eqn{\hat\sigma^{*2}} over that
#' distribution. The argument `type` of [predict.StatmodFit()], the variance
#' type of [vcov.StatmodFit()], is read only with this predictive.
#'
#' # Plug-in
#'
#' The family at the estimates, averaged only over the effects: the
#' estimative interval, which ignores the error in the estimates and so covers
#' less than the level in a small sample.
#'
#' # Bootstrap
#'
#' The parametric bootstrap predictive distribution itself: `n` responses are
#' simulated from the fit with the random effects drawn afresh, each response
#' is refitted, and the family at each refit is averaged over the effects.
#' `refit = "coefficients"` refits at the fit's hyperparameters; `"full"`
#' chooses them again on every replica, so their uncertainty enters too, at
#' the cost of a whole fit each. The bootstrap is rejected for a model with a
#' block that moves with its coefficients ([modelterms7::seg()],
#' [modelterms7::jump()], [modelterms7::jseg()], [modelterms7::nl()]), for a
#' response that is not one number per observation, and, with
#' `refit = "full"`, for a cross-validated criterion. A replica whose refit
#' fails is left out, with a warning that reports how many were. The seed is
#' the caller's: `set.seed()` before the call fixes the replicas. See
#' [predictive_bootstrap()].
#'
#' # Hyperparameters
#'
#' The interval for a new group ([group_interval()]) is conditional on the
#' covariance that the fit estimated, and so is the interval for a new
#' observation, so the uncertainty of the hyperparameters is not propagated.
#' The exception is [pred_bootstrap()] with `refit = "full"`, where the
#' hyperparameters are chosen again on every replica.
#'
#' @return Not applicable: this topic describes the constructors. Each returns
#'   an S7 object that inherits from [predictive_type()].
#'
#' @references Harris, I. R. (1989). Predictive fit for natural exponential
#'   families. *Biometrika*, 76, 675--684.
#'
#' @seealso [pred_averaged()], [pred_plugin()], [pred_bootstrap()],
#'   [prediction_interval()], [predict_intervals].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' nd <- data.frame(x = 0.5)
#' for (p in list(pred_averaged(), pred_plugin())) {
#'   print(predict(fit, what = "response", newdata = nd,
#'                 interval = prediction_interval(predictive = p)))
#' }
NULL


#' @title The Confidence Interval of a Prediction
#'
#' @description
#' Selects the interval of a parameter, a predictor or a moment that carries
#' only the uncertainty of the estimates. It is the default `interval` of
#' [predict.StatmodFit()] and is reported with `se = TRUE`. Choose it for the
#' typical group's value; for the value of a new group use [group_interval()],
#' and for a new response use [prediction_interval()].
#'
#' @details
#' The standard error is the delta method on the scale of the equation, and the
#' interval is carried through the link. For an effect set aside by `random`
#' the interval is that of the typical group's value. See [predict_intervals].
#'
#' @return An S7 object of class `ConfidenceInterval`, inheriting from
#'   [interval_type()], with the one property `label`.
#'
#' @seealso [group_interval()], [prediction_interval()], [predict_intervals],
#'   [predict.StatmodFit()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' predict(fit, what = "mu", newdata = data.frame(x = 0.5), se = TRUE,
#'         interval = confidence_interval())
#'
#' @export
confidence_interval <- function() {
  ConfidenceInterval(label = "confidence")
}

#' @title The Interval of a New Group's Parameter
#'
#' @description
#' Selects the interval of the parameter of a new group, that is, the
#' parameter at a level of a [modelterms7::random()] term that is not in the
#' fitting data. It adds the variation of the effects of a new group to the
#' uncertainty of the estimates. It requires `random = "zero"` or
#' `"marginal"` in [predict.StatmodFit()] and implies `se = TRUE`.
#'
#' @details
#' The ends are quantiles of the new group's predictor carried through the
#' inverse link; the standard error is the standard deviation of the new
#' group's parameter. See [predict_intervals] for the formulas and the
#' restrictions.
#'
#' @return An S7 object of class `GroupInterval`, inheriting from
#'   [interval_type()], with the one property `label`.
#'
#' @seealso [confidence_interval()], [prediction_interval()],
#'   [predict_intervals], [predict.StatmodFit()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(g = factor(rep(1:8, each = 5)), x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(8, sd = 0.5)[dd$g] + rnorm(40, sd = 0.3)
#' fit <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
#'                dd)
#' predict(fit, what = "mu", newdata = data.frame(x = 0.5, g = "new"),
#'         random = "zero", interval = group_interval())
#'
#' @export
group_interval <- function() {
  GroupInterval(label = "group")
}

#' @title The Interval of a New Observation
#'
#' @description
#' Selects the interval in which a new observation of the response falls. It
#' requires `what = "response"` in [predict.StatmodFit()] and does not read
#' `se`, because its result carries its own `se` column. It adds the variation
#' of the family to that of the estimates and of the effects set aside; the
#' argument `predictive` sets how the estimation uncertainty enters.
#'
#' @details
#' The family is averaged over the predictors of every parameter. See
#' [predict_intervals] for the interval and [predictive_distributions] for
#' the three values of `predictive`.
#'
#' @param predictive A [predictive_type()] object:
#'   [pred_averaged()] (the default) averages the family over the approximate
#'   sampling distribution of the estimates, [pred_plugin()] uses the
#'   estimates alone, and [pred_bootstrap()] simulates and refits. See
#'   [predictive_distributions].
#'
#' @return An S7 object of class `PredictionInterval`, inheriting from
#'   [interval_type()], with the properties `label` and `predictive`.
#'
#' @seealso [confidence_interval()], [group_interval()],
#'   [predict_intervals], [predictive_distributions],
#'   [predict.StatmodFit()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' nd <- data.frame(x = 0.5)
#' predict(fit, what = "response", newdata = nd,
#'         interval = prediction_interval())
#' predict(fit, what = "response", newdata = nd,
#'         interval = prediction_interval(predictive = pred_plugin()))
#'
#' @export
prediction_interval <- function(predictive = pred_averaged()) {
  if (!S7::S7_inherits(predictive, predictive_type)) {
    stop("'predictive' must be a predictive object: pred_averaged(), ",
         "pred_plugin() or pred_bootstrap().", call. = FALSE)
  }
  PredictionInterval(label = "prediction", predictive = predictive)
}

#' @title Predictive Distribution Averaged Over the Estimates
#'
#' @description
#' The default `predictive` of [prediction_interval()]: the family averaged
#' over the predictors at the fit's estimates, with the covariance of the
#' estimates plus the covariance that the set-aside effects add. It is the
#' normal approximation to the parametric bootstrap predictive distribution
#' (Harris, 1989), and it covers a new observation at close to the level in a
#' small sample, where [pred_plugin()] covers less. It costs one evaluation of
#' the covariance of the estimates, so it is the choice when [pred_bootstrap()]
#' is too expensive.
#'
#' @details
#' It is the only predictive that reads the argument `type` of
#' [predict.StatmodFit()]. See [predictive_distributions].
#'
#' @return An S7 object of class `PredAveraged`, inheriting from
#'   [predictive_type()], with the one property `label`.
#'
#' @references Harris, I. R. (1989). Predictive fit for natural exponential
#'   families. *Biometrika*, 76, 675--684.
#'
#' @seealso [pred_plugin()], [pred_bootstrap()], [prediction_interval()],
#'   [predictive_distributions].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' predict(fit, what = "response", newdata = data.frame(x = 0.5),
#'         interval = prediction_interval(predictive = pred_averaged()))
#'
#' @export
pred_averaged <- function() {
  PredAveraged(label = "averaged")
}

#' @title Plug-in Predictive Distribution
#'
#' @description
#' The family at the estimates, averaged only over the effects that
#' `random` sets aside. The interval is the estimative one, which ignores the
#' error in the estimates and so covers less than the level in a small sample.
#' It costs one evaluation of the predictors and is the cheapest of the three
#' predictive objects.
#'
#' @details
#' See [predictive_distributions].
#'
#' @return An S7 object of class `PredPlugin`, inheriting from
#'   [predictive_type()], with the one property `label`.
#'
#' @seealso [pred_averaged()], [pred_bootstrap()], [prediction_interval()],
#'   [predictive_distributions].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' predict(fit, what = "response", newdata = data.frame(x = 0.5),
#'         interval = prediction_interval(predictive = pred_plugin()))
#'
#' @export
pred_plugin <- function() {
  PredPlugin(label = "plugin")
}

#' @title Parametric Bootstrap Predictive Distribution
#'
#' @description
#' The parametric bootstrap predictive distribution of Harris (1989): `n`
#' responses are simulated from the fit with the random effects drawn afresh,
#' each response is refitted, and the family at each refit is averaged over the
#' effects. It does not rely on the normal approximation of [pred_averaged()],
#' and it costs `n` refits.
#'
#' @details
#' With `refit = "coefficients"` each replica is refitted at the fit's
#' hyperparameters. With `refit = "full"` the hyperparameters are chosen again
#' on every replica, so their uncertainty enters the interval as well, at the
#' cost of a whole fit per replica. The bootstrap is rejected for a model with
#' a block that moves with its coefficients, for a response that is not one
#' number per observation, and, with `refit = "full"`, for a cross-validated
#' criterion. The replicas use the random number generator of the session, so
#' `set.seed()` before the call reproduces them. If the fit holds a
#' stored parametric bootstrap ([add_bootstrap()]) with the same `refit`, its
#' first `n` replicates are used and nothing is refitted; an entry with fewer
#' than `n` replicates is extended and stored again. With
#' [modelterms7::random()] terms the entry is used when `random` sets every
#' one of them aside and the entry holds the refitted random effects. See
#' [predictive_distributions] and [predictive_bootstrap()].
#'
#' @param n The number of replicas, a single positive whole number. The default
#'   is `200`.
#' @param refit `"coefficients"` (the default), refitting each replica at the
#'   fit's hyperparameters, or `"full"`, choosing them again.
#'
#' @return An S7 object of class `PredBootstrap`, inheriting from
#'   [predictive_type()], with the properties `label`, `n` (an integer) and
#'   `refit`.
#'
#' @references Harris, I. R. (1989). Predictive fit for natural exponential
#'   families. *Biometrika*, 76, 675--684.
#'
#' @seealso [pred_averaged()], [pred_plugin()], [prediction_interval()],
#'   [predictive_distributions].
#'
#' @examples
#' set.seed(2)
#' dd <- data.frame(x = runif(30))
#' dd$y <- 1 + dd$x + rnorm(30, sd = 0.5)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' predict(fit, what = "response", newdata = data.frame(x = 0.5),
#'         interval = prediction_interval(
#'           predictive = pred_bootstrap(n = 20, refit = "coefficients")))
#'
#' @export
pred_bootstrap <- function(n = 200L, refit = c("coefficients", "full")) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
      n != round(n)) {
    stop("'n' must be a single positive whole number.", call. = FALSE)
  }
  refit <- match.arg(refit)
  PredBootstrap(label = paste0("bootstrap, ", as.integer(n), " replicas, ",
                               "refit ", refit),
                n = as.integer(n), refit = refit)
}


#' @title Print an Interval or a Predictive Distribution
#' @name print.interval_type
#' @aliases print.predictive_type
#' @description
#' Prints a one-line description of an [interval_type()] or a
#' [predictive_type()] object: its label, and for [prediction_interval()] the
#' predictive distribution that it holds.
#' @param x An [interval_type()] or [predictive_type()] object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' print(confidence_interval())
#' print(prediction_interval(predictive = pred_bootstrap(n = 50)))
#' @keywords internal
S7::method(print, interval_type) <- function(x, ...) {
  cat("<interval> ", x@label, "\n", sep = "")
  invisible(x)
}

S7::method(print, PredictionInterval) <- function(x, ...) {
  cat("<interval> ", x@label, " (predictive: ", x@predictive@label, ")\n",
      sep = "")
  invisible(x)
}

S7::method(print, predictive_type) <- function(x, ...) {
  cat("<predictive> ", x@label, "\n", sep = "")
  invisible(x)
}


#' @title Check the Interval Argument of predict()
#'
#' @description
#' Validates the `interval` argument of [predict.StatmodFit()] and returns the
#' kind of interval and, for a prediction interval, the kind and the settings
#' of the predictive distribution, as the strings that the body of the
#' prediction reads.
#'
#' @param interval The argument as passed.
#'
#' @return A list with `interval` (`"confidence"`, `"group"` or
#'   `"prediction"`), `predictive` (`"averaged"`, `"plugin"` or
#'   `"bootstrap"`), `n_boot` and `boot_refit`.
#'
#' @keywords internal
predict_interval_spec <- function(interval) {
  if (!S7::S7_inherits(interval, interval_type)) {
    stop("'interval' must be an interval object: confidence_interval(), ",
         "group_interval() or prediction_interval().", call. = FALSE)
  }
  out <- list(interval = "confidence", predictive = "averaged",
              n_boot = 200L, boot_refit = "coefficients")
  if (S7::S7_inherits(interval, GroupInterval)) {
    out$interval <- "group"
  } else if (S7::S7_inherits(interval, PredictionInterval)) {
    out$interval <- "prediction"
    p <- interval@predictive
    if (S7::S7_inherits(p, PredPlugin)) {
      out$predictive <- "plugin"
    } else if (S7::S7_inherits(p, PredBootstrap)) {
      out$predictive <- "bootstrap"
      out$n_boot <- p@n
      out$boot_refit <- p@refit
    }
  } else if (!S7::S7_inherits(interval, ConfidenceInterval)) {
    stop("'interval' is an interval of an unknown kind; use ",
         "confidence_interval(), group_interval() or ",
         "prediction_interval().", call. = FALSE)
  }
  out
}
