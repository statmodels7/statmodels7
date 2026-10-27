#' @include inference.R
NULL

# The hyperparameters, estimated by a marginal criterion.
#
# The penalized objective is F(b; t) = -l(b) + sum_j rho_j(b; t), and every
# penalties7 penalty keeps its normalizing constant, so rho IS minus a log
# prior density and nothing has to be added by hand. The marginal likelihood
# is then
#
#   L(t) = integral exp(-F(b; t)) db
#
# and its Laplace approximation at the penalized mode b(t) is
#
#   log L(t) = l(b) - rho(b; t) + (q/2) log(2 pi) - (1/2) log|A'(H + S)A|
#
# with H the information of the likelihood, S the penalty's Hessian, and A an
# orthonormal basis of whichever subspace is integrated over. Writing it out
# with the constants in place reproduces Wood's (2011) REML criterion term for
# term, the -(M_p/2) log(2 pi) of that paper appearing here as the difference
# between the rank of the penalty and the dimension of the coefficient vector.

#' How the Hyperparameters Are Estimated
#'
#' @description
#' Every outer criterion is an instance of this class: [reml()], [ml()],
#' [aic()], [bic()] and [cv()] all return one. It carries which subspace of
#' the coefficients is integrated over, which information matrix enters the
#' determinant, and the settings that a path or a cross-validation needs.
#'
#' One class serves every criterion, so a property that a given criterion does
#' not read is present and unused. The class validates the values of its
#' properties. The constructors [reml()], [ml()], [aic()], [bic()] and [cv()]
#' also check `information`, match `rule` against its choices and check `k`
#' and `marginal`.
#'
#' @param kind `"reml"`, `"ml"`, `"aic"`, `"bic"` or `"cv"`. Determines which
#'   criterion is evaluated and, through [outer_minimize()], whether it is
#'   made small or large.
#' @param hessian An object that inherits from [Information]: [observed_info()],
#'   [expected_info()] or [auto_info()]. Which information enters the Laplace
#'   determinant. [statmod()] settles [auto_info()] against the family through
#'   [outer_resolve()]. The constructors set it from their argument
#'   `information`.
#' @param k The price of one effective degree of freedom, for a
#'   prediction-error criterion. It is `NA_real_` for [bic()], which resolves
#'   it to \eqn{\log n} at fit time, and for the criteria that do not charge
#'   per degree of freedom.
#' @param nfolds How many folds [cv()] uses. Only `"cv"` reads it.
#' @param rule `"min"` or `"1se"`. Only `"cv"` reads it.
#' @param folds A fold number per observation, or `integer(0)` for folds
#'   drawn at fit time. Only `"cv"` reads it.
#' @param marginal Which distribution parameters have their unpenalized
#'   coefficients estimated on the marginal criterion: `NULL` for every
#'   parameter except the position, `"none"`, `"all"`, or a vector of
#'   names. Only `"reml"` and `"ml"` read it; see [reml()].
#'
#' @return An object of class `OuterMethod` with one property per argument
#'   above.
#'
#' @seealso [reml()] and [ml()] for the marginal criteria, [aic()] and
#'   [bic()] for the prediction-error ones, [cv()] for cross-validation,
#'   [statmod()] for where one is passed.
#'
#' @examples
#' reml()
#' ml(information = observed_info())
#' aic()
#' cv(nfolds = 5)
#'
#' @name OuterMethod-class
#' @aliases OuterMethod
#' @keywords internal
#' @export
OuterMethod <- S7::new_class("OuterMethod",
  properties = list(
    kind = S7::class_character,
    hessian = Information,
    k = S7::class_numeric,
    nfolds = S7::class_numeric,
    rule = S7::class_character,
    folds = S7::class_numeric,
    marginal = S7::new_property(S7::new_union(NULL, S7::class_character),
                                default = NULL)
  ),
  validator = function(self) {
    if (!identical(length(self@kind), 1L) ||
        !self@kind %in% c("reml", "ml", "aic", "bic", "cv")) {
      return(paste0("Property 'kind' must be \"reml\", \"ml\", \"aic\", ",
                    "\"bic\" or \"cv\"."))
    }
    if (!identical(length(self@rule), 1L) ||
        !self@rule %in% c("min", "1se")) {
      return("Property 'rule' must be \"min\" or \"1se\".")
    }
    if (length(self@nfolds) != 1L || self@nfolds < 2) {
      return("Property 'nfolds' must be a single number, at least 2.")
    }
    if (!is.null(self@marginal) &&
        (!length(self@marginal) || anyNA(self@marginal) ||
         any(!nzchar(self@marginal)))) {
      return(paste0("Property 'marginal' must be NULL or a non-empty ",
                    "character vector of parameter names, \"none\" or ",
                    "\"all\"."))
    }
    NULL
  }
)


#' The Properties That Every Criterion Carries
#'
#' @description
#' Returns the defaults for `nfolds`, `rule` and `folds`, which only [cv()]
#' reads, so that one class can serve every criterion and each constructor
#' names only what it needs.
#'
#' @details
#' The settings of a path are not among them. How many values a path
#' visits, how far down it reaches, and whether a term's own hyperparameters
#' are combined or swept one at a time all belong to the term. [aic()] and
#' [bic()] are also used for a model's smooth hyperparameters, which are read
#' at the mode and not swept, so a path setting on the criterion would not
#' apply to them. [cv()] is accepted only as a `sparse_criterion`.
#'
#' @return A named list with `nfolds`, `rule` and `folds`, to be spliced into
#'   an [OuterMethod()] call.
#'
#' @seealso [path_fallbacks()] for the settings to which a term falls back
#'   when it declares a kinked penalty without offering a grid argument,
#'   [modelterms7::term_search()] for where a path's shape is set.
#'
#' @keywords internal
outer_path_defaults <- function() {
  list(nfolds = 10, rule = "min", folds = numeric(0))
}


#' Estimate the Hyperparameters by a Marginal Likelihood
#'
#' @description
#' Build the criterion object that [statmod()] uses to choose a model's smooth
#' hyperparameters by a marginal likelihood: the coefficients are integrated
#' out by a Laplace approximation at the penalized mode, and the remainder is
#' maximized in the hyperparameters. Pass the result as `outer_criterion`.
#'
#' `reml()` integrates every coefficient that `marginal` does not name. `ml()`
#' integrates only the penalized directions and profiles the rest. `reml()` is
#' the default of [statmod()].
#'
#' @details
#' # The criterion
#'
#' At the penalized mode \eqn{\hat\beta(\theta)},
#' \deqn{\log L(\theta) = \ell(\hat\beta) - \rho(\hat\beta;\theta)
#'   + \frac{q}{2}\log 2\pi - \frac12\log|A'(H+S)A|,}
#' with \eqn{H} the information of the log-likelihood, \eqn{S} the penalty's
#' second derivative in the coefficients, and \eqn{A} an orthonormal basis of
#' the subspace integrated over. The term \eqn{\rho} needs no correction: a
#' \pkg{penalties7} penalty keeps its normalizing constant, so it is exactly
#' minus a log prior density, and for a quadratic penalty that constant
#' carries the \eqn{-\frac{r}{2}\log\lambda} and the log pseudo-determinant
#' that a marginal criterion needs. Written out, the expression reproduces
#' Wood's (2011) REML criterion term for term.
#'
#' # The effective degrees of freedom of the fit
#'
#' The effective degrees of freedom stored in the fitted model (`fit@edf`),
#' the share of each term in \eqn{\mathrm{tr}[(H+S)^{-1}H]}, are read on the
#' `information` of the criterion: [distributions7::observed_info()] by
#' default, or the [distributions7::expected_info()] object with its `method`.
#' The count is not part of the criterion here, but the standard errors of the
#' same fit read the same default matrix, so the count and the standard errors
#' of a summary agree. [summary.StatmodFit()] recomputes the count when its
#' `vcov_type` names the other matrix. For a positive semidefinite \eqn{H} the
#' total lies in \eqn{[0, p]}, with \eqn{p} the number of coefficients; an
#' indefinite observed information can give a direction of negative curvature
#' a negative share, and the summary then adds a note.
#'
#' # What each one integrates
#'
#' Write \eqn{\gamma} for the unpenalized coefficients of the parameters that
#' `marginal` names, and \eqn{u} for every other coefficient. The criterion
#' is maximized in \eqn{\gamma} together with the hyperparameters. At each
#' value of \eqn{\gamma}, the coefficients \eqn{u} are set to the penalized
#' mode.
#'
#' `reml()` integrates every coefficient in \eqn{u}. An unpenalized one is
#' integrated under a flat prior; the absence of a penalty amounts to a flat
#' prior. `ml()` integrates only the directions in the range space of the
#' penalty, so an unpenalized coefficient in \eqn{u} is read at the mode. An
#' ordinary covariate is unpenalized. The linear component of a
#' Demmler-Reinsch smooth is unpenalized too, because the penalty of that
#' smooth does not act on it.
#'
#' With `marginal = "all"`, \eqn{\gamma} holds every unpenalized
#' coefficient. The two criteria then integrate the same directions, and they
#' are the same function.
#'
#' Where the criterion estimates coefficients, a coefficient that the model
#' does not identify at the first fit is left out of the determinant for the
#' whole search and stays free in the inner fit. This corresponds to the REML
#' log-likelihood of [stats::lm()] (`logLik(fit, REML = TRUE)`), which uses the
#' rank of the design in place of its number of columns.
#'
#' # Where a dispersion and a shape are estimated
#'
#' By default `marginal` names every parameter except the position, which is
#' the family's first parameter: the location, or the scale for a family with
#' no location, such as [distributions7::weibull1_distrib()]. For a gaussian
#' mixed model this is the REML of `lme` and `glmmTMB`. The coefficients of
#' the mean are integrated, and \eqn{\sigma} is maximized on the criterion.
#' Without any penalty, `reml()` returns the \eqn{\sigma} of [stats::lm()],
#' \eqn{\sqrt{\mathrm{rss}/(n-p)}}, and its criterion equals
#' `logLik(lm(...), REML = TRUE)`. `ml()` returns \eqn{\sqrt{\mathrm{rss}/n}}.
#' If `agec` is defined as `age - 11`, the model
#' `distance ~ agec + random(~ agec | Subject)` on `nlme::Orthodont` gives the
#' estimates and the criterion of `lme`.
#'
#' `marginal = "none"` gives the convention of `gamlss` and of mgcv's
#' `gaulss`: every parameter of the distribution is read at the joint mode of
#' the penalized likelihood. That estimate ignores the degrees of freedom of
#' the coefficients estimated beside it, so a dispersion read there is biased
#' downwards. Where a dispersion carries a random effect of its own, the bias
#' at the joint mode is larger than the bias on the criterion.
#'
#' Where some coefficients of the location are integrated, the estimate
#' depends on the parametrization of the family. The dispersion of
#' [distributions7::gamma1_distrib()] and the variance of
#' [distributions7::gamma2_distrib()], \eqn{\sigma^2 = \phi\mu^2}, differ
#' under `reml()`, because the second contains the integrated mean. With
#' `"all"` the criterion does not depend on the parametrization.
#'
#' Some configurations are not covered, and [marginal_coords()] lists them: a
#' structural term anywhere in the model, a block that is a working
#' linearization (a sharp [modelterms7::jump()] or [modelterms7::jseg()])
#' anywhere in the model, a block that moves with its coefficients in the
#' parameter's own equation, and a quadratic penalty in that equation whose
#' null space is not spanned by coordinates. In each case a parameter that the
#' default names keeps the convention of `"none"`, and a parameter named
#' explicitly is rejected. Beside a penalty with a kink the default names no
#' parameter, and a parameter named explicitly is estimated on the criterion,
#' and the coordinates that the kinked penalty covers are excluded. For a
#' multivariate family the default names no parameter, and an explicit
#' `marginal` other than `"none"` is an error.
#'
#' # Which hyperparameters
#'
#' The criterion estimates the hyperparameters of the terms fitted in one
#' system, meaning the terms whose penalty is twice differentiable. A lasso, a
#' SCAD or an MCP has a kink, its coefficients are estimated by a method of
#' their own, and a Laplace approximation is undefined at a point where the
#' second derivative does not exist. Those hyperparameters are not on the outer
#' search. [statmod()] chooses them by a path under `sparse_criterion`, which is
#' [bic()] by default and can also be [aic()] or [cv()].
#'
#' # The exact gradient
#'
#' The criterion has an exact gradient on the observed and on the expected
#' information, for every penalty whose derivatives in the coefficients and the
#' hyperparameters are available, as those of [modelterms7::s()],
#' [modelterms7::te()] and any [penalties7::quadratic_penalty()] are; a model
#' with a [modelterms7::regime()] term has none. The gradient, and the Hessian
#' where it is exact, are supplied to the search. The default optimizer is
#' [optimizers7::newton()] with a bounded step where the Hessian is exact and
#' steers the search, [optimizers7::lbfgs()] where only the gradient is exact,
#' [optimizers7::brent()] for a single coordinate without a gradient, and
#' [optimizers7::nelder_mead()] otherwise. A simplex needs one vertex per
#' dimension and a quasi-Newton method does not, so an exact gradient saves
#' criterion evaluations once there are two or more hyperparameters.
#'
#' # ML needs a null basis
#'
#' `ml()` reads the null space of every penalty over which it integrates,
#' because the null space determines which directions are profiled.
#' [penalties7::is_proper()] identifies a penalty with no null space, and
#' [penalties7::penalty_null_basis()] returns the null space of a quadratic,
#' structured, additive or separable penalty. A penalty that provides neither
#' is rejected by name, and the criterion does not assume a subspace for it.
#'
#' @param information The information that enters the determinant, an object
#'   that inherits from [Information]:
#'   - [observed_info()] (the default): the curvature of the log-likelihood at
#'     the data, which makes the criterion the Laplace approximation;
#'   - [expected_info()]: the Fisher information, a function of the parameters
#'     alone. Its `method` selects how the expectation is obtained: `NULL`
#'     (the closed form of the family where it has one, and
#'     [distributions7::expect_opg()] otherwise), [distributions7::expect_opg()],
#'     [distributions7::expect_integrate()] or [distributions7::expect_mc()]
#'     (see [distributions7::expected_derivative_methods]);
#'   - [auto_info()]: settled against the family by [statmod()] through
#'     [outer_resolve()], as the expected information for a Student t, a skew t
#'     or a Cauchy response, whose observed information is not positive
#'     definite at an outlier, and as the observed information otherwise.
#'
#'   The expected information is used where the caller names it, as for a
#'   Student t whose degrees of freedom are modeled. The fitted model records
#'   the information used.
#'
#'   On the expected information, the `method` object defines the matrix and
#'   every derivative of it that the exact outer gradient and Hessian read, so
#'   that the derivatives are those of the criterion that is evaluated. With the
#'   closed form of the family, the derivative is the closed form of the
#'   derivative (every shipped family; the Hessian also needs the second
#'   derivative of the information, which every shipped family supplies except
#'   [distributions7::enet_distrib()]). With [distributions7::expect_integrate()]
#'   or [distributions7::expect_mc()], it is the derivative of the expectation,
#'   which moves the density as well as the integrand. With
#'   [distributions7::expect_opg()], the matrix is the outer product of the
#'   scores at the data and the derivative is its exact derivative at fixed
#'   data. A family without a closed form reads [distributions7::expect_opg()]
#'   for the matrix and for the derivatives when `method` is `NULL`. The outer
#'   product depends on the response, so the criterion built on it has exact
#'   derivatives and is not the Laplace approximation of an expectation.
#'   With `method = NULL` on a family without a closed form, the matrix of the
#'   criterion is therefore an estimate of the expected information, noisier
#'   than a computed expectation and about 25 times cheaper: on a density-only
#'   family (n = 300, ten basis functions) a fit took 0.8 s against 25.3 s for
#'   `method = expect_integrate(identity = "score")`, with smoothing
#'   parameters 12.57 and 13.00, against 13.02 for the observed information.
#'   The computed expectation is selected with
#'   `method = expect_integrate(identity = "score")`.
#'   [distributions7::expect_mc()] draws a new sample at each call of the
#'   matrix and of each derivative, so that its derivative is that of the
#'   expectation and not of the sampled matrix, and the outer search is noisy at
#'   the level of the Monte Carlo error. [information_objects] compares the three
#'   objects.
#'
#' @param marginal Which distribution parameters have their unpenalized
#'   coefficients estimated by maximizing this criterion, instead of being
#'   read at the joint mode (`ml()`) or integrated (`reml()`). `NULL`, the
#'   default, names every parameter except the position, which is the
#'   family's first parameter. `"none"` names no parameter, `"all"` names
#'   every one, and a character vector names the parameters listed. With
#'   `"all"`, `reml()` and `ml()` are the same criterion. For a multivariate
#'   family the default names no parameter, and naming one explicitly is an
#'   error.
#'   Where the model carries a penalty with a kink (a lasso, an elastic
#'   net, a SCAD or an MCP, in any equation), `NULL` names no parameter:
#'   the path that chooses the kinked penalty scores the model at the joint
#'   mode, and the fit returned is the model that it scored. A parameter
#'   named explicitly is still estimated on the criterion, at the value that
#'   the path chose. Where the family lacks the two derivatives that the
#'   criterion reads in some parameter (a Laplace response) and the model
#'   carries no smooth hyperparameter, the default criterion is not run at
#'   all, and a `reml()` or `ml()` passed by name is rejected.
#'
#' @return An [OuterMethod()] object of kind `"reml"` or `"ml"`, whose
#'   `hessian` holds the object given as `information` and whose `marginal` is
#'   as supplied; its path settings are unused.
#'
#' @references
#' Wood, S. N. (2011). Fast stable restricted maximum likelihood and marginal
#' likelihood estimation of semiparametric generalized linear models.
#' *Journal of the Royal Statistical Society, Series B*, 73(1), 3--36.
#'
#' @seealso [aic()], [bic()] and [cv()] for the prediction-error criteria,
#'   [information_objects] for the argument `information`, [statmod()] for
#'   where these are passed, [hyper()] for reading back what they chose.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(200, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(200, sd = 0.3)
#'
#' fit <- statmod(y ~ s(x, bspline_smooth(k = 10)), distributions7::gaussian1_distrib(), dd,
#'                outer_criterion = reml())
#'
#' # The smoothing parameter was estimated, and hyper() reports it.
#' hyper(fit)
#'
#' # ML profiles the unpenalized directions instead of integrating them, so
#' # it shrinks a little less. The gap is small here because only two of the
#' # ten coefficients are unpenalized; it widens with the fixed effects.
#' fml <- statmod(y ~ s(x, bspline_smooth(k = 10)), distributions7::gaussian1_distrib(), dd,
#'                outer_criterion = ml())
#' c(reml = hyper(fit)$estimate, ml = hyper(fml)$estimate)
#' c(reml = sum(fit@edf$edf), ml = sum(fml@edf$edf))
#'
#' # The marginal log-likelihood, which these criteria maximize, is only
#' # available where one of them ran.
#' logLik(fit, loglik_type = "marginal")
#'
#' @export
reml <- function(information = observed_info(), marginal = NULL) {
  check_information(information, arg = "information", fn = "reml")
  do.call(OuterMethod, c(list(kind = "reml", hessian = information,
                             k = NA_real_,
                             marginal = check_marginal_arg(marginal)),
                         outer_path_defaults()))
}

#' @rdname reml
#' @export
ml <- function(information = observed_info(), marginal = NULL) {
  check_information(information, arg = "information", fn = "ml")
  do.call(OuterMethod, c(list(kind = "ml", hessian = information,
                             k = NA_real_,
                             marginal = check_marginal_arg(marginal)),
                         outer_path_defaults()))
}


#' Settle a Marginal Criterion's Information Against the Family
#'
#' @description
#' Replaces [auto_info()] by the information that the criterion uses for
#' `distrib`: [expected_info()] for a Student t, a skew t or a Cauchy response,
#' wrapped or not, and [observed_info()] for every other family. The
#' approximation that the [auto_info()] object carries passes to the
#' [expected_info()] object. Any other information is returned unchanged.
#'
#' @details
#' The log-density of these families is not concave in the location, so
#' their observed information is not positive definite at an outlier, and
#' the Laplace determinant built on it can fail there. The expected
#' information is therefore used for these families, and the observed
#' information for the others.
#'
#' @param method An [OuterMethod()].
#' @param distrib The response family.
#'
#' @return `method`, in which `hessian` is settled.
#'
#' @keywords internal
outer_resolve <- function(method, distrib) {
  if (!info_auto(method@hessian)) return(method)
  method@hessian <- if (observed_indefinite(distrib))
    expected_info(method = info_emethod(method@hessian)) else observed_info()
  method
}


#' Whether a Family's Observed Information Can Be Indefinite
#'
#' @description
#' Returns `TRUE` for a Student t, a skew t or a Cauchy response, including a
#' family that has one of them as its parent. The log-density of these
#' families is not concave in the location.
#'
#' @param distrib A family.
#'
#' @return A single logical.
#'
#' @keywords internal
observed_indefinite <- function(distrib) {
  t_like <- paste0("distributions7::", c("StudentT1Distrib", "StudentT2Distrib",
                                         "SkewTDistrib", "CauchyDistrib"))
  if (any(class(distrib) %in% t_like)) return(TRUE)
  if (S7::prop_exists(distrib, "parent_distrib")) {
    return(observed_indefinite(distrib@parent_distrib))
  }
  FALSE
}


#' Check the `marginal` Argument of a Marginal Criterion
#'
#' @description
#' Checks the shape of `marginal` before any family is available: the value
#' must be `NULL` or a character vector with no missing and no empty entry.
#' `"none"` and `"all"` must stand alone. Whether the names belong to the
#' family is checked by [marginal_params()], at fit time.
#'
#' @param marginal The argument as given.
#'
#' @return `marginal`, in which repeated names are removed, or an error that
#'   names the argument.
#'
#' @keywords internal
check_marginal_arg <- function(marginal) {
  if (is.null(marginal)) return(NULL)
  if (!is.character(marginal) || !length(marginal) || anyNA(marginal) ||
      any(!nzchar(marginal))) {
    stop(paste0("'marginal' must be NULL, \"none\", \"all\" or a character",
                " vector of\n  distribution parameter names."),
         call. = FALSE)
  }
  if (length(marginal) > 1L && any(marginal %in% c("none", "all"))) {
    stop("'marginal' takes \"none\" or \"all\" alone, not beside names.",
         call. = FALSE)
  }
  unique(marginal)
}


#' @title Print an Outer Method
#' @name print.OuterMethod
#' @description
#' Prints one line naming the criterion and the information on which it is
#' built, as `<REML>  observed information`, and for [cv()] the number of
#' folds and the selection rule.
#' @param x An [OuterMethod()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @seealso [reml()], [aic()] and [cv()] for the objects printed.
#' @keywords internal
print.OuterMethod <- function(x, ...) {
  if (identical(x@kind, "cv")) {
    cat(sprintf("<CV>  %s folds, rule \"%s\"\n",
                if (length(x@folds)) "given" else format(x@nfolds), x@rule))
    return(invisible(x))
  }
  cat(sprintf("<%s>  %s information\n", toupper(x@kind),
              info_label(x@hessian)))
  invisible(x)
}
S7::method(print, OuterMethod) <- print.OuterMethod


#' The Hyperparameters That an Outer Method Estimates
#'
#' @description
#' Returns one row per hyperparameter of a penalized term that is fitted in the
#' joint system, with the link that carries it onto the whole line.
#'
#' @details
#' A term whose penalty has a kink is left out. Its coefficients are
#' estimated by a coordinate descent while everything else is held fixed, and
#' the criterion is a Laplace approximation, which needs a second derivative
#' that does not exist at a kink. Those hyperparameters are chosen by a path
#' instead.
#'
#' A hyperparameter that its own term holds is left out too, whatever the kind
#' of penalty to which it belongs, and so is every member of a shared label
#' when one member holds it.
#'
#' The result is the index of the hyperparameters on the outer search. The
#' search also moves the coefficients that the criterion estimates (see
#' [marginal_coords()]), so its dimension is the number of rows plus the number
#' of those coefficients. If there are no rows and no such coefficients, the
#' criterion does not run.
#'
#' # Shared hyperparameters
#'
#' Hyperparameters carrying the same label through
#' [modelterms7::term_ids()] are estimated at one value, so they occupy one
#' row: the search has one coordinate for the group. The row is identified by
#' its first member, so that every lookup that reads `parameter`, `term` and
#' `name` keeps working, and the `"members"` attribute lists the
#' hyperparameters that the row represents.
#'
#' Sharing does not merge the penalties, which remain two objects estimated at
#' one value, so the derivative of the row is the sum of the derivatives of its
#' members and the value is written back under the key of each member. If no
#' hyperparameter is shared, the member table is the index itself, one line per
#' row.
#'
#' Two members must carry the same link, since one free value has to lie in
#' both of their domains, and a label may not span a smooth penalty and a
#' kinked one, which are estimated by two different procedures and would be
#' given two different values. Both cases are rejected here, where both sides
#' are visible.
#'
#' @param spec A [StatmodSpec()].
#' @param blocks The block split, as [statmod_blocks()] returns it.
#'
#' @return A data frame with one row per estimated hyperparameter, or per
#'   group of shared ones, and three character columns, `parameter`, `term`
#'   and `name`. The links that carry each onto the whole line are in the
#'   `"links"` attribute, a list with one element per row of the frame (a list
#'   column would not survive the subsetting that the search does). The
#'   `"members"` attribute is a data frame of the same three columns plus
#'   `row`, which lists the hyperparameters that each row represents.
#'
#' @seealso [hyper_to_eta()] for the map that those links define,
#'   [statmod_held_ids()] for the shared labels that a term holds,
#'   [index_members()] for the member table.
#'
#' @keywords internal
outer_hyper_index <- function(spec, blocks) {
  kinked <- vapply(blocks$sparse, function(b)
    paste(b$param, b$term, sep = "\r"), character(1))
  rows <- list()
  links <- list()
  labs <- character(0)
  kink_labs <- character(0)
  units <- statmod_penalty_keys(spec)
  # a label one member holds is held for the group -- see statmod_held_ids()
  held_ids <- names(statmod_held_ids(units))
  for (u in units) {
    lab_of <- function(h) if (h %in% names(u$ids)) u$ids[[h]] else NA_character_
    if (paste(u$param, u$key, sep = "\r") %in% kinked) {
      # its label is recorded but its hyperparameters are not, so that a
      # label spanning the two machines can be caught below
      for (h in u$penalty@params) {
        if (h %in% names(u$fixed)) next
        l <- lab_of(h)
        if (!is.na(l)) kink_labs <- c(kink_labs, l)
      }
      next
    }
    for (h in u$penalty@params) {
      # HELD by the term, so there is nothing here to estimate. The term is
      # where the penalty is named and where the caller says so.
      if (h %in% names(u$fixed)) next
      if (lab_of(h) %in% held_ids) next
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = u$param, term = u$key, name = h,
        stringsAsFactors = FALSE)
      links[[length(links) + 1L]] <- u$penalty@link_params[[h]]
      labs <- c(labs, lab_of(h))
    }
  }
  if (!length(rows)) return(outer_hyper_index_none())
  mem <- do.call(rbind, rows)
  index_group(mem, links, labs, kink_labs)
}


#' The Outer Index With No Hyperparameter
#'
#' @description
#' Returns the empty index that [outer_hyper_index()] returns where no
#' hyperparameter is estimated, with the same columns and attributes as a full
#' one.
#'
#' @return A data frame of zero rows with columns `parameter`, `term` and
#'   `name`, and the attributes `links` and `members`.
#'
#' @examples
#' nrow(statmodels7:::outer_hyper_index_none())
#'
#' @keywords internal
outer_hyper_index_none <- function() {
  structure(data.frame(parameter = character(0), term = character(0),
                       name = character(0), stringsAsFactors = FALSE),
            links = list(),
            members = data.frame(row = integer(0), parameter = character(0),
                                 term = character(0), name = character(0),
                                 stringsAsFactors = FALSE))
}


#' Collapse the Shared Hyperparameters of an Index
#'
#' @description
#' Turns the member table into the index on which the search runs. The index
#' has one row per group, identified by the group's first member, and carries
#' the member table and the links alongside.
#'
#' @details
#' A label ties its members into one coordinate. Two checks run here because
#' this is where both sides are visible: every member of a group must carry
#' the same link (compared as objects, so two links of one class with different
#' parameters differ), because one free value has to lie in the domain of each
#' member, and a label may not span a smooth penalty and a kinked one, which
#' are estimated by two different procedures and would be given two values
#' under one name.
#'
#' If no hyperparameter is shared, each row is its own group and the member
#' table is the index.
#'
#' @param mem A data frame of `parameter`, `term` and `name`, one line per
#'   estimated hyperparameter.
#' @param links The links, one per line of `mem`.
#' @param labs The sharing labels, one per line of `mem`, `NA` where the
#'   hyperparameter is not shared.
#' @param kink_labs The labels carried by hyperparameters that a path
#'   estimates.
#'
#' @return The index, with its `"links"` and `"members"` attributes.
#'
#' @keywords internal
index_group <- function(mem, links, labs, kink_labs = character(0)) {
  bad <- intersect(stats::na.omit(labs), kink_labs)
  if (length(bad)) {
    stop(sprintf(paste0("the label '%s' is written on a smooth penalty and",
                        " on a kinked one.\n  The two are estimated by",
                        " different machines -- a marginal criterion at the",
                        "\n  mode, and a path over a grid -- so one value",
                        " cannot serve both."), bad[1L]), call. = FALSE)
  }
  # an unshared hyperparameter is its own group, and is given a key nothing
  # else can collide with rather than a label of its own
  key <- ifelse(is.na(labs), paste0("\r", seq_along(labs)), labs)
  first <- !duplicated(key)
  row <- match(key, key[first])
  for (g in which(first)) {
    same <- which(row == row[g])
    if (length(same) < 2L) next
    # the links themselves and not their classes: two power links of
    # different exponents are one class and two domains
    cls <- vapply(links[same], function(l) {
      pr <- tryCatch(unlist(l@link_params), error = function(e) NULL)
      if (length(pr)) {
        paste0(class(l)[[1L]], "(",
               paste(names(pr), format(pr), sep = " = ", collapse = ", "), ")")
      } else class(l)[[1L]]
    }, "")
    same_link <- vapply(links[same], identical, logical(1),
                        y = links[[same[[1L]]]])
    if (!all(same_link)) {
      k <- which(!same_link)[[1L]]
      stop(sprintf(paste0("the label '%s' ties '%s' of '%s' to '%s' of",
                          " '%s', whose\n  links differ (%s against %s). One",
                          " value cannot lie in both domains."),
                   labs[same[[1L]]], mem$name[same[[1L]]],
                   mem$term[same[[1L]]], mem$name[same[[k]]],
                   mem$term[same[[k]]], cls[[1L]], cls[[k]]),
           call. = FALSE)
    }
  }
  idx <- mem[first, , drop = FALSE]
  rownames(idx) <- NULL
  members <- cbind(row = row, mem, stringsAsFactors = FALSE)
  rownames(members) <- NULL
  structure(idx, links = links[first], members = members)
}


#' The Hyperparameters Behind Each Row of the Outer Index
#'
#' @description
#' Returns the member table: one line per estimated hyperparameter, with the
#' row of the index that stands for it.
#'
#' @details
#' If no hyperparameter is shared, the member table is the index itself with a
#' `row` column counting up. Where a label ties several hyperparameters, they
#' carry the same `row`, and a quantity accumulated over members lands on the
#' coordinate of the group.
#'
#' @param idx The index, as [outer_hyper_index()] returns it.
#'
#' @return A data frame with columns `row`, `parameter`, `term` and `name`.
#'
#' @keywords internal
index_members <- function(idx) {
  m <- attr(idx, "members")
  if (!is.null(m)) return(m)
  cbind(row = seq_len(nrow(idx)), idx, stringsAsFactors = FALSE)
}


#' Move Between the Hyperparameters and the Free Vector
#'
#' @description
#' `hyper_to_eta()` carries the estimated hyperparameters onto the whole
#' line through their links; `eta_to_hyper()` puts a free vector back.
#'
#' @details
#' The outer search runs on the free scale for the reason that every other
#' search in the toolkit does: a smoothing parameter is positive, and an
#' optimizer that has no information about the domain can step outside it. The
#' link of each hyperparameter is read from the `"links"` attribute of the
#' index that [outer_hyper_index()] returns.
#'
#' The two are inverses over the estimated hyperparameters:
#' `eta_to_hyper(hyper_to_eta(h, idx), idx, h)` returns `h`. The
#' hyperparameters that `idx` does not name are untouched by both.
#'
#' @param hyper The hyperparameter structure, keyed by distribution parameter
#'   and then by term.
#' @param idx The index, as [outer_hyper_index()] returns it, whose rows
#'   identify the hyperparameters on the search and whose `"links"` attribute
#'   gives the link that carries each.
#' @param eta A free vector, as long as `idx` has rows.
#'
#' @return `hyper_to_eta()` gives an unnamed numeric vector with one entry
#'   per row of `idx`. `eta_to_hyper()` gives the hyperparameter structure,
#'   in which the estimated entries are replaced and every other entry is as
#'   it was.
#'
#' @seealso [outer_hyper_index()] for the index and its links,
#'   [outer_fit()] for the search that runs on this scale.
#'
#' @keywords internal
hyper_to_eta <- function(hyper, idx) {
  links <- attr(idx, "links")
  vapply(seq_len(nrow(idx)), function(i) {
    v <- hyper[[idx$parameter[i]]][[idx$term[i]]][[idx$name[i]]]
    linkfunctions7::linkfun(links[[i]], v)
  }, numeric(1))
}

#' @rdname hyper_to_eta
#' @keywords internal
eta_to_hyper <- function(eta, idx, hyper) {
  links <- attr(idx, "links")
  # ONE value per row, written under EVERY member's own key: sharing does not
  # merge the penalties, so each goes on being a penalty with a hyperparameter
  # of its own and every reader of the store finds what it expects there. With
  # nothing shared there is one member per row and this is the loop over rows.
  mem <- index_members(idx)
  for (i in seq_len(nrow(mem))) {
    r <- mem$row[i]
    hyper[[mem$parameter[i]]][[mem$term[i]]][[mem$name[i]]] <-
      linkfunctions7::linkinv(links[[r]], eta[[r]])
  }
  hyper
}


#' The Subspace Over Which a Marginal Criterion Integrates
#'
#' @description
#' Returns `NULL` for REML, which integrates everything, and an orthonormal
#' basis of the penalty's range space for ML.
#'
#' @details
#' The total penalty is block diagonal over the terms, and each term owns its
#' own columns, so the null space of the total penalty is the direct sum of the
#' unpenalized columns and each penalty's own null space. The basis is
#' assembled that way.
#'
#' It is not read off the eigenvalues of the assembled matrix. A count of
#' eigenvalues above a tolerance depends on the scale of the smoothing
#' parameters: if two smoothing parameters are many orders of magnitude apart,
#' the small contributions fall below the tolerance, while the null space
#' itself does not change. A penalty with a kink contributes no column, and a
#' penalty over the own parameters of a structural term is appended by
#' [statmod_marginal_full()].
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param kind `"reml"` or `"ml"`.
#' @param gamma `TRUE` where the criterion also estimates coefficients
#'   ([marginal_coords()]); the function then returns an empty basis, instead
#'   of an error, for a model with no penalized direction.
#'
#' @return `NULL` for `"reml"`, which integrates every direction and needs no
#'   basis. For `"ml"`, it is a `p x r` matrix with orthonormal columns
#'   spanning the range space of the penalties, where `p` is the total
#'   coefficient count and `r` is the total rank. The matrix has no columns
#'   (\eqn{r = 0}) if no direction of the design is penalized and either the
#'   criterion estimates coefficients or a structural term is penalized.
#'   Otherwise an error is signaled.
#'
#' @seealso [penalty_range_basis()] for one penalty's contribution,
#'   [reml()] and [ml()] for the distinction that this implements.
#'
#' @keywords internal
integrated_basis <- function(spec, design, kind, gamma = FALSE) {
  if (identical(kind, "reml")) return(NULL)
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  cols <- list()
  for (u in statmod_penalized(spec, design)) {
    # a penalty over a structural term's own parameters owns no column of the
    # design; its directions are appended by statmod_marginal_full()
    if (isTRUE(u$structural)) next
    pen <- u$penalty
    # a kinked penalty's coordinates stay at the joint mode and are not
    # integrated, which is laplace_pinned()'s rule for reml()
    if (!isTRUE(u$mixed) &&
        isTRUE(tryCatch(penalty_has_kink(pen), error = function(e) TRUE))) {
      next
    }
    # A MIXED class owns columns for HALF its coordinates, and the split is
    # exact only because a class prior is proper: penalty_range_basis() is
    # then the identity, so each coordinate is its own direction and the two
    # halves go to the two readers unchanged. A class prior that was not
    # proper would have a range mixing the two spaces, which neither reader
    # can express, and ml() says so rather than projecting onto half of it.
    idx_u <- u$index
    if (isTRUE(u$mixed)) {
      if (!isTRUE(penalties7::is_proper(pen))) {
        stop(sprintf(paste0(
          "ml() cannot read the null space of the covariance block '%s',",
          " which\n  spans both the coefficients and a filter's own",
          " parameters. Use reml(),\n  which integrates every coordinate",
          " and needs no such basis."), u$key), call. = FALSE)
      }
      idx_u <- u$beta_index
    }
    k <- length(idx_u)
    R <- penalty_range_basis(pen, k, u$param, u$key)
    if (!ncol(R)) next
    M <- matrix(0, total, ncol(R))
    M[idx_u, ] <- R
    cols[[length(cols) + 1L]] <- M
  }
  # nothing penalized and coefficients the criterion estimates: ml() is then
  # the profile likelihood, and the subspace it integrates over is empty
  if (!length(cols) && (structural_penalized(spec, design) || isTRUE(gamma))) {
    return(matrix(0, total, 0))
  }
  if (!length(cols)) {
    stop(paste0("ml() has no penalized direction over which to integrate.\n",
                "  Every penalty here is entirely null, which is not a model",
                " that ml() can\n  read; reml() integrates the whole coefficient",
                " vector instead."), call. = FALSE)
  }
  do.call(cbind, cols)
}


#' An Orthonormal Basis of a Penalty's Range Space
#'
#' @description
#' Returns an orthonormal basis of the directions that a penalty constrains.
#' [ml()] integrates over these directions, and [reml()] integrates them
#' together with the others.
#'
#' @details
#' A proper penalty constrains every direction, so the basis is the identity
#' and [penalties7::is_proper()] is the only function that is called. A
#' penalty with a null space supplies it through
#' [penalties7::penalty_null_basis()], and the range basis is the orthogonal
#' complement.
#'
#' A penalty that supplies neither is rejected by name, because integrating
#' over a guessed subspace would not give the marginal likelihood.
#'
#' @param pen A \pkg{penalties7} penalty object.
#' @param k The number of coefficients in the term's block.
#' @param p The distribution parameter, named in the error message.
#' @param nm The term's name, named in the error message.
#'
#' @return A `k x r` matrix with orthonormal columns, where `r` is the
#'   penalty's rank. It is the `k x k` identity for a proper penalty.
#'
#' @seealso [integrated_basis()], which assembles these blockwise.
#'
#' @keywords internal
penalty_range_basis <- function(pen, k, p, nm) {
  if (isTRUE(penalties7::is_proper(pen))) return(diag(1, k, k))
  N <- tryCatch(penalties7::penalty_null_basis(pen), error = function(e) NULL)
  if (is.null(N)) {
    stop(sprintf(paste0("ml() cannot read the null space of the penalty on",
                        " '%s' in '%s'.\n  Its class does not expose one, so",
                        " the profiled directions cannot be\n  guessed.",
                        " Use reml(), which integrates every",
                        "\n  coefficient and needs no such basis."), nm, p),
         call. = FALSE)
  }
  N <- as.matrix(N)
  if (!ncol(N)) return(diag(1, k, k))
  # the range is the orthogonal complement of the null space inside the block
  Q <- qr.Q(qr(N), complete = TRUE)
  Q[, seq.int(ncol(N) + 1L, k), drop = FALSE]
}


#' Does a Structural Term Carry a Penalty of Its Own?
#'
#' @description
#' Reports whether any penalty of the model covers the own parameters of a
#' structural term, as opposed to a block of design columns, as with the
#' deviations of a score-driven filter over a panel. A covariance class that
#' spans both the coefficients and those parameters counts as well.
#'
#' @details
#' The answer determines the matrix over which the Laplace determinant is
#' taken. `TRUE` sends the criterion to [statmod_marginal_full()], which spans
#' the coefficients and the term's free parameters; `FALSE` leaves it over
#' the coefficients only.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, whose structural state holds the terms.
#'
#' @return A single logical. It is `TRUE` where a penalty covers a structural
#'   term's own parameters, or a covariance class spans both the coefficients
#'   and those parameters, and `FALSE` for a model with no structural term and
#'   for one whose structural term has no such penalty. `mixed_penalized()`
#'   returns `TRUE` for the second case only.
#'
#' @seealso [statmod_marginal_full()] for the matrix that it selects,
#'   [statmod_structural_penalty()] for the penalties that it seeks.
#'
#' @keywords internal
structural_penalized <- function(spec, design) {
  for (u in statmod_penalized(spec, design)) {
    # a MIXED class shrinks a filter's own parameters too, half of its
    # coordinates being among them, so the determinant has to span them for
    # the same reason. Answering FALSE here sent the criterion to the route
    # over the coefficients alone, where a mixed class contributes nothing at
    # all -- statmod_penalty_at() leaves it to joint_penalty_at() -- so the
    # coefficient block came back unpenalized and singular.
    if (isTRUE(u$structural) || isTRUE(u$mixed)) return(TRUE)
  }
  FALSE
}

#' @rdname structural_penalized
#' @keywords internal
mixed_penalized <- function(spec, design) {
  for (u in statmod_penalized(spec, design)) {
    if (isTRUE(u$mixed)) return(TRUE)
  }
  FALSE
}


#' The Penalized Curvature Over the Coefficients and a Filter's Parameters
#'
#' @description
#' Returns the matrix whose determinant the Laplace approximation reads, where
#' a penalty covers a structural term's own parameters instead of a block of
#' design columns.
#'
#' @details
#' # Why the determinant has to span them
#'
#' A marginal criterion integrates the quantities that its penalty shrinks.
#' Where the penalty is a structural term's own, as with the deviations of a
#' panel, those quantities are not coefficients.
#'
#' Taken over the coefficients only, the determinant then does not depend on
#' the hyperparameter at all, and the criterion reduces to the penalized
#' likelihood, whose maximum in a shrinkage parameter is at no shrinkage. The
#' criterion then differs from the marginal likelihood.
#'
#' # The matrix
#'
#' [statmod_full_information()] spans the coefficients followed by the term's
#' free parameters, which is the order that the joint fit uses, and
#' [statmod_structural_penalty()] gives the Hessian of the penalty in those same
#' parameters. A covariance class split between the two halves is placed by
#' [joint_penalty_at()].
#'
#' The information here is the observed one whatever the `hessian` of the
#' method is. The expected information over the own parameters of a filter is
#' not a quantity that the toolkit carries, since the state of the recursion
#' enters the expectation.
#'
#' # Under ML
#'
#' The tail is integrated only over the penalized coordinates, through the
#' penalty's own range basis, so a parameter that the penalty does not cover,
#' such as a population value, is profiled instead of integrated. An
#' unpenalized coefficient is profiled too.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients, at the penalized mode.
#' @param hyper The hyperparameters.
#' @param basis The integrated subspace over the coefficients, as
#'   [integrated_basis()] returns it, or `NULL` for [reml()].
#'
#' @return A square symmetric matrix spanning the integrated coefficients
#'   followed by the integrated structural parameters. The value is `NULL`
#'   where the structural term could not be run at the given point, which a
#'   caller reads as an unusable hyperparameter.
#'
#' @seealso [statmod_marginal()], the caller,
#'   [structural_penalized()] for the predicate that selects this route.
#'
#' @keywords internal
statmod_marginal_full <- function(spec, design, coef, hyper, basis = NULL) {
  sst <- statmod_structural_state(design)
  su <- Filter(function(u) identical(u$kind, "filter"),
               attr(design, "structural"))
  if (is.null(sst) || !length(su)) return(NULL)
  key <- su[[1L]]$term
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  nb <- sum(vapply(design, function(d) d$npar, integer(1)))
  ix <- nb + seq_along(free)

  K <- tryCatch(statmod_full_information(spec, coef, design),
                error = function(e) NULL)
  if (is.null(K) || nrow(K) != nb + length(free)) return(NULL)
  S <- matrix(0, nrow(K), ncol(K))
  # as_dense() for the same reason the joint objective's own Hessian gives:
  # this matrix spans the coefficients AND the structural term's parameters
  # and is dense by construction, so a sparse block cannot be written into a
  # slice of it. This site is the twin of the one in statmod_fit_joint(), and
  # only that one failed when the penalty's accumulator followed the design --
  # a formula carrying both a filter and a random effect is what reaches it.
  S[seq_len(nb), seq_len(nb)] <-
    as_dense(statmod_penalty_at(spec, coef, hyper, design, "hessian"))
  ps <- structural_penalty_block(spec, design, hyper, length(free))
  if (!is.null(ps)) S[ix, ix] <- ps
  # a class split between the two halves is placed whole, cross block and all,
  # by the one function the joint step and the variance also read
  S <- S + joint_penalty_at(spec, design, coef, hyper, "hessian", nrow(S))
  S <- zap_nonfinite(S)
  M <- K + S
  if (is.null(basis)) return(M)
  # the same subspace the gradient projects onto, composed once: two callers
  # building it separately would agree only by accident
  A <- structural_joint_basis(spec, design, key, free, nb, basis)
  crossprod(A, M %*% A)
}


#' Which of a Structural Term's Free Parameters a Penalty Covers
#'
#' @description
#' Returns the positions, among a structural term's free parameters, that
#' some penalty of the model shrinks. Those are the directions over which
#' [ml()] integrates; the rest are profiled.
#'
#' @details
#' A term declares its penalties through `modelterms7::term_penalties()`,
#' and each entry names a subset of the term's own parameters. This collects
#' those subsets and maps them onto positions in `free`.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param key The term's key in the specification, its call as written.
#' @param free The term's free parameter names, in the order in which the term
#'   holds them.
#'
#' @return An integer vector of positions into `free`, sorted and without
#'   repeats. The vector is empty if the term is not covered by any penalty.
#'
#' @seealso [statmod_marginal_full()], the caller.
#'
#' @keywords internal
structural_range_cols <- function(spec, design, key, free) {
  sst <- statmod_structural_state(design)
  nm <- names(sst$zeta[[key]])
  out <- integer(0)
  for (u in statmod_penalized(spec, design)) {
    # a MIXED class shrinks free parameters too -- its structural half -- and
    # those directions are integrated over exactly as a structural unit's are
    if (isTRUE(u$mixed) && identical(u$class$sterm, key)) {
      zc <- unlist(lapply(Filter(function(z) isTRUE(z$structural), u$pieces),
                          function(z) z$zcols), use.names = FALSE)
      out <- c(out, match(nm[zc], free))
      next
    }
    if (!isTRUE(u$structural) || !identical(u$term, key)) next
    out <- c(out, match(nm[u$cols], free))
  }
  sort(unique(out[!is.na(out)]))
}


#' The Marginal Criterion at Given Coefficients and Hyperparameters
#'
#' @description
#' Computes the Laplace approximation to the log marginal likelihood,
#' evaluated at the penalized mode: the log-likelihood minus the penalty, plus
#' \eqn{(q/2)\log 2\pi}, minus half the log determinant of the penalized
#' information over the integrated subspace. A coordinate with a non-finite
#' curvature is pinned to the identity (see [pin_boundary()]), which lowers
#' \eqn{q} by one. The determinant is read through [pd_logdet()].
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients, at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#' @param basis The integrated subspace, from [integrated_basis()].
#' @param ctx An evaluation context from [outer_context()] built at the same
#'   coefficients and hyperparameters, or `NULL` to assemble the information
#'   and the penalty afresh.
#'
#' @return A list with `value`, `loglik`, `penalty`,
#'   `logdet` and `q`, or `NULL` where the determinant does not
#'   exist. Where the integrated subspace is empty, `value` is the profile
#'   log-likelihood minus the penalty, and `logdet` and `q` are zero.
#'
#' @keywords internal
statmod_marginal <- function(spec, design, coef, hyper, method,
                             emethod = NULL, basis = NULL, ctx = NULL) {
  expected <- info_expected(method@hessian)
  ll <- statmod_loglik_at(spec, coef, design)
  rho <- statmod_penalty_at(spec, coef, hyper, design, "value")
  if (structural_penalized(spec, design)) {
    M <- statmod_marginal_full(spec, design, coef, hyper, basis)
    if (is.null(M)) return(NULL)
  } else {
    # ⚠️ Sharing the factorization with ctx_penalized() was MEASURED AND NOT
    # TAKEN. Where the criterion's matrix is the one that context holds -- the
    # observed information, nothing projected away -- the same matrix is
    # factorized here and there, which at p = 503 is 12.4 ms spent twice, and
    # reading the determinant from the context removes one of them. End to
    # end it is worth nothing: 1.01x at p = 503 and 0.92x to 1.04x over the
    # four shapes, because the criterion is evaluated at many points the
    # gradient never reaches (every trial point of a line search) and there
    # is nothing to share at those. The cheap route below is sparse where the
    # matrix is, which serves the criterion-only points as well.
    H <- ctx_information(ctx, spec, design, coef, hyper, expected, emethod)
    S <- ctx_penalty(ctx, spec, design, coef, hyper)
    M <- H + S
    if (!is.null(basis)) M <- crossprod(basis, M %*% basis)
  }
  # a determinant that does not exist is a hyperparameter the search should
  # step away from, not an error: at a far-out value the penalized information
  # can lose definiteness while the fit itself is sound.
  #
  # It must be a statement about the MATRIX, though, and reading whether
  # chol() raised is not one. Measured on a hierarchical score-driven panel,
  # K+S reaches a condition number of 8.0e15 while still being positive
  # definite, and there whether the factorization succeeds is decided by
  # rounding: the outer search saw a dozen consecutive points as unavailable
  # while backtracking towards one that had been available a moment before,
  # which is a hole in the domain that is not there. pd_logdet() keeps the
  # cheap route where it is safe and pays for the eigendecomposition only
  # where the condition estimate says the cheap one cannot be trusted.
  # A COORDINATE AT A BOUNDARY IS NOT ONE THE LAPLACE APPROXIMATION INTEGRATES
  # OVER. Its curvature there is NaN, which denies the whole matrix a
  # factorization, and it is pinned to the identity instead: the determinant
  # gains log(1) = 0 and the dimension q drops by one, so the coordinate
  # contributes neither a determinant term nor a 2*pi. Measured on a Student
  # t whose nu had reached double.xmax, the criterion went from unavailable
  # -- which stopped the whole fit at its first evaluation -- to a finite
  # value at a converged mode. See pin_boundary().
  M <- pin_boundary(M, if (is.null(basis)) laplace_pinned(spec, design) else
                      integer(0))
  held <- attr(M, "held")
  # nothing left to integrate: ml() with every unpenalized direction
  # estimated by the criterion and no penalty, the profile likelihood
  if (!nrow(M)) {
    return(list(value = ll - rho, loglik = ll, penalty = rho, logdet = 0,
                q = 0L))
  }
  ld <- pd_logdet(M)
  if (!isTRUE(ld$ok)) return(NULL)
  q <- nrow(M) - (if (is.null(held)) 0L else held)
  logdet <- ld$logdet
  list(value = ll - rho + q / 2 * log(2 * pi) - logdet / 2,
       loglik = ll, penalty = rho, logdet = logdet, q = q)
}


#' Estimate the Hyperparameters
#'
#' @description
#' Runs `outer_optimizer` on the marginal criterion, refitting the
#' coefficients at every hyperparameter that it tries.
#'
#' @details
#' Each evaluation is a whole inner fit, warm-started from the incumbent, the
#' best usable point found so far. Near the optimum the coefficients move
#' little, so the inner fit converges in few iterations.
#'
#' The optimizer is chosen by [outer_default_optimizer()] from the derivatives
#' that the criterion supplies: [optimizers7::newton()] where the exact Hessian
#' is used, [optimizers7::lbfgs()] where the gradient is exact (or a covariance
#' class spans a structural term's parameters), [optimizers7::brent()] for a
#' single coordinate without a gradient, and [optimizers7::nelder_mead()]
#' otherwise. An optimizer given explicitly is used as given, and an optimizer
#' that needs a gradient that is not available signals its own error.
#'
#' @param spec The specification.
#' @param design The design.
#' @param blocks The block split.
#' @param hyper The starting hyperparameters.
#' @param inner_optimizer How the smooth block is fitted.
#' @param method An [OuterMethod()].
#' @param optimizer An \pkg{optimizers7} optimizer, or `NULL` to choose one with
#'   [outer_default_optimizer()] from the derivatives that the criterion
#'   supplies. With `NULL`, the stopping rule and the line search are also set
#'   by the function.
#' @param beta The starting coefficients, stacked.
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#' @param maxit,tol The alternation's budget and tolerance.
#' @param vb The resolved verbosity.
#' @param hold_hyper `TRUE` to hold the hyperparameters at `hyper` and
#'   search over the criterion's own coefficients only.
#'
#' @return A list with `par`, `hyper`, `value`, `criterion`, `converged`,
#'   `settled`, `obj`, `hist_blocks`, `hist_inner`, `aliased`, `spec`,
#'   `hist_outer`, `iterations`, `evaluations`, `exact_gradient`,
#'   `exact_hessian`, `optimizer`, `held_coef`, `marginal` and `pinned`.
#'   `hist_outer` holds the trace of the search.
#'
#' @seealso [reml()], [statmod()]
#'
#' @keywords internal
outer_fit <- function(spec, design, blocks, hyper, inner_optimizer, method,
                      optimizer, beta, emethod, maxit, tol, vb,
                      hold_hyper = FALSE) {
  # with the hyperparameters held the search runs over the criterion's own
  # coefficients alone, which is what a bootstrap replica refitted at the
  # fit's hyperparameters needs
  idx <- outer_hyper_index(spec, blocks)
  if (isTRUE(hold_hyper)) idx <- outer_hyper_index_none()
  # the coefficients the criterion estimates beside the hyperparameters: the
  # outer vector is (eta, gamma), eta first, so every consumer that reads the
  # first nrow(idx) entries as hyperparameters goes on doing so
  gam <- if (outer_minimize(method)) marginal_coords(spec, design, NULL) else
    marginal_coords(spec, design, method)
  nh <- nrow(idx)
  ng <- length(gam$where)
  nv <- nh + ng
  if (!nh && !ng) {
    # statmod() does not reach this: a criterion applies to the smooth
    # penalties and comes into play only where there is one, so a model with
    # none leaves it unused rather than failing. What this guards is the
    # internal contract, for a caller reaching the outer loop directly.
    stop(paste0("there is no hyperparameter for a marginal criterion to\n",
                "  estimate: the terms here carry no penalty that is twice\n",
                "  differentiable. A lasso, a scad or an mcp is chosen by a",
                "\n  path over its own values instead -- bic(), aic() or cv() --",
                " because its penalty\n  has a kink."),
         call. = FALSE)
  }
  expected <- info_expected(method@hessian)
  basis <- integrated_basis(spec, design, method@kind, gamma = ng > 0L)
  labels <- c(paste(idx$parameter, idx$term, idx$name, sep = "/"),
              paste(gam$param, gam$name, sep = "/"))
  # what the TRACE prints. The full label names the column of hist_outer,
  # where it has to stay unique and reconstructible; a reader of a running
  # fit needs the term and not its whole specification.
  shown <- c(paste(idx$parameter, short_keys(idx$term), idx$name, sep = "/"),
             paste(gam$param, gam$name, sep = "/"))
  # the value of the j-th coordinate of v on the scale a reader reads: a
  # hyperparameter on its own scale, a coefficient as it is
  v_value <- function(hy, v, j) {
    if (j <= nh) hyper_value(hy, idx, j) else v[[j]]
  }
  # the specification with the estimated coefficients held at v's values
  spec_at <- function(v) statmod_hold(spec, gam, v[nh + seq_len(ng)])
  # the inner objective's gradient with the held coordinates left out: at the
  # held mode their score is not zero, and a mode error or a resolution is a
  # statement about the free coordinates alone
  free_score <- function(g) {
    if (!is.null(g) && ng) g[gam$where] <- 0
    g
  }
  pe <- outer_minimize(method)
  basis <- if (pe) NULL else basis
  exact <- outer_gradient_ok(spec, design, idx, method, 1L, gamma = ng > 0L)
  exact2 <- exact && outer_gradient_ok(spec, design, idx, method, 2L,
                                       gamma = ng > 0L)
  # whether the STOPPING RULE is this package's business too. An optimizer
  # given by name comes with its own rule and keeps it; one chosen here is
  # chosen whole, the rule included, because only this package knows what the
  # criterion it is being pointed at can resolve.
  chose_optimizer <- is.null(optimizer)
  if (chose_optimizer) {
    optimizer <- outer_default_optimizer(
      exact, outer_newton_ok(spec, design, exact2, method),
      mixed_penalized(spec, design), nv)
  }
  # Whether the TRACE prints a gradient. It reports what the search is using,
  # so where the search uses none there is none to report, and asking would
  # compute the quantity the laziness below exists to avoid.
  trace_grad <- exact &&
    "gradient" %in% optimizers7::optimizer_provides(optimizer)

  # the search minimizes; a marginal likelihood is maximized and goes in with
  # its sign turned, a prediction-error criterion as it is. The method says
  # which it is rather than the search assuming.
  sgn <- if (pe) 1 else -1

  state <- new.env(parent = emptyenv())
  state$beta <- beta
  state$inner <- NULL
  # the INCUMBENT: the best usable point so far, which is the iterate a
  # monotone line search stands on. See the warm start in evaluate().
  state$best_v <- Inf
  state$best_beta <- NULL
  state$best_zeta <- NULL
  state$rows <- list()
  state$evals <- 0L
  state$key <- NULL
  # the worst value seen ON THE SEARCH'S OWN SCALE, which is what an
  # unavailable point is made worse than
  state$worst <- NA_real_
  # What the criterion can resolve here: every usable reading so far, and the
  # summary of them resolution_summary() takes. See there for why a summary
  # and not the smallest.
  state$res_seen <- numeric(0)
  state$resolution <- NA_real_

  # WHICH criterion this search is running, written once. A prediction-error
  # method and a marginal one are read by different functions and only one of
  # them takes the integrated basis, so a second reader that picked for itself
  # would answer for the wrong quantity -- which is what
  # criterion_resolution() did on an aic() fit before this existed.
  criterion_at <- function(cf, hy, par, ctx = NULL, sp = state$sp) {
    if (is.null(sp)) sp <- spec
    if (pe) {
      statmod_pe(sp, design, cf, hy, method, emethod,
                 statmod_active(sp, blocks, par, hy))
    } else {
      statmod_marginal(sp, design, cf, hy, method, emethod, basis, ctx)
    }
  }

  # one inner fit serves the value and the gradient: an optimizer asks for
  # both at the same point, and refitting for the second would double the cost
  # of every step
  evaluate <- function(v) {
    key <- paste(format(v, digits = 17), collapse = ",")
    if (identical(state$key, key)) return(state$last)
    eta <- v[seq_len(nh)]
    hy <- eta_to_hyper(eta, idx, hyper)
    sp <- spec_at(v)
    state$sp <- sp
    # a hyperparameter far enough out takes the model somewhere its own
    # machinery cannot go -- a variance beyond what a density can represent,
    # a curvature that is no longer definite. That is a step the search should
    # take back, not an error. At the STARTING point it is neither: nothing
    # has moved yet, so a failure there is the caller's and is raised.
    m <- NULL
    res <- NULL
    cf <- NULL
    ctx <- NULL
    err <- NULL
    # A step the search takes back must leave NO TRACE, and until now only
    # half of it did. The coefficients are protected -- `state$beta` is
    # written below, after the point is known to be usable -- but a
    # structural term's own parameters live in the design's structural state,
    # which is an ENVIRONMENT the inner fit writes into as it goes, and
    # `statmod_fit_structural()` stores its optimizer's last point whether
    # that optimizer converged or not. So an unavailable point moved the
    # filter's parameters and the next evaluation started from wherever the
    # failure had left them, which ratchets: measured on a panel of twenty
    # groups, the search took one long first step, and then every one of the
    # thirty backtracked points was unavailable, including points beside the
    # start that had converged at the first evaluation.
    sst <- statmod_structural_state(design)
    # ⚠️ EVERY TRIAL WARM-STARTS FROM THE INCUMBENT, not from the last point
    # tried. The search's accept/reject happens inside optimizers7 and is not
    # visible here, so writing the warm start at every USABLE point let a
    # rejected trial move it; where the inner problem is multimodal the chain
    # then left the incumbent's basin and could not come back. Measured on a
    # gas panel with a developed loading: the first evaluation, at sigma = 1,
    # read -1367.360; every later one -- the backtracks converging back onto
    # sigma = 1 included -- read -1370.656 from a rejected trial's filter
    # state, so the line search compared against a value it could no longer
    # reproduce, backtracked to exhaustion on all three routes and reported
    # the start at the WORSE of the two values.
    warm <- state$beta
    if (!is.null(state$best_beta)) {
      warm <- state$best_beta
      if (!is.null(sst) && !is.null(state$best_zeta)) {
        sst$zeta <- state$best_zeta
        sst$key <- NULL
        sst$value <- NULL
      }
    }
    z_save <- if (is.null(sst)) NULL else sst$zeta
    if (ng) warm[gam$where] <- v[nh + seq_len(ng)]
    body <- function() {
      res <<- statmod_alternate(sp, design, blocks, hy, inner_optimizer,
                                warm, expected, emethod, maxit, tol,
                                vb_inner(vb), hold_refresh = TRUE)
      cf <<- res$obj$split(res$par)
      # Both criteria are read AT THE MODE -- a Laplace approximation there,
      # or a log-likelihood and an edf there -- so where the inner fit did not
      # reach one there is no criterion to compare, and the point is
      # unavailable rather than merely worse. This used to be enforced by
      # accident: a wild hyperparameter made the inner fit RAISE, the catch
      # below turned that into an unavailable point, and the search stepped
      # back. Once the inner fit stopped raising and started returning its
      # last usable point, the criterion at that point looked ordinary and the
      # search walked to a smoothing parameter of 1.8e308.
      # one context per point: the criterion, the gradient and the Hessian
      # read the same information, the same penalty and the same factorization
      # instead of each assembling its own. It is built BEFORE the point is
      # judged, because what judges it reads that same factorization.
      ctx <<- outer_context(sp, design, cf, hy, emethod)
      # ⚠️ AND THE FLAG IS NOT WHAT DECIDES, because it answers another
      # question. The flag says whether the inner stopping rule fired;
      # availability asks whether the criterion, a Laplace expansion AT THE
      # MODE, is valid here. Measured on `y ~ s(x) | sigma ~ s(z)`, of 38
      # inner fits in one search 38 are at their mode by the second reading
      # -- 1e-09 to 3e-09 against a limit of 1e-03 -- and FOUR by the first.
      # The other 34 stopped on the objective-stall guard with the objective
      # already fixed to twelve significant digits and the score oscillating
      # between 2.5e-06 and 3.3e-06, just above its absolute tolerance of
      # 1e-06. Read as unavailable, they made this line search backtrack
      # eleven times an iteration and accept a step of 0.0026 where the
      # Newton step is 1.4, so the search moved 0.005 in eta over 38
      # evaluations and stopped 4.0 criterion units short of the optimum its
      # own gradient was correctly pointing at.
      #
      # The rule ADDS points and never removes one -- a converged run stays
      # usable whatever the mode error reads -- so no model that fitted
      # before can stop fitting. Making the flag STRICTER is the other
      # direction, and piano_stabilita.txt 13d measured it and withdrew it:
      # it cost a false negative on a good fit.
      if (!isTRUE(res$converged)) {
        q <- inner_mode_error(ctx, sp, design, cf, hy,
                              free_score(tryCatch(res$obj$gr(res$par),
                                                  error = function(e) NULL)),
                              expected)
        if (!is.finite(q) || q > mode_error_limit()) {
          m <<- NULL
          return(invisible(NULL))
        }
      }
      m <<- criterion_at(cf, hy, res$par, ctx)
      if (is.null(m)) return(invisible(NULL))
      # The derivatives are NOT computed here. Which of them the search reads
      # is the OPTIMIZER's business and it varies: nelder_mead reads neither,
      # lbfgs and bfgs read the gradient and never the Hessian, newton reads
      # both. Computing them because the criterion COULD supply them paid for
      # the dearest quantity in the fit -- the Hessian is 49 to 78 per cent of
      # it -- and threw it away. They are produced on demand instead, which
      # needs no predicate about the optimizer and follows what it actually
      # does rather than what its class declares.
      invisible(NULL)
    }
    if (state$evals == 0L) {
      body()
      # A COORDINATE THE MODEL DOES NOT IDENTIFY IS LEFT OUT OF THE
      # DETERMINANT, as lm() computes its REML criterion with the rank rather
      # than with the number of columns. Where the criterion estimates a
      # coefficient, it runs on a model with no smooth hyperparameter as
      # well, and there a flat direction of K -- which the fit only warned
      # about before -- makes log|K| undefined at the start. The set is read
      # ONCE, at this first fit, and held for the whole search, so the
      # criterion does not change dimension from one point to the next.
      if (is.null(m) && ng && !is.null(res)) {
        Kf <- tryCatch(statmod_penalized_at(sp, cf, design, hy, expected,
                                            emethod),
                       error = function(e) NULL)
        fl <- if (is.null(Kf)) integer(0) else
          sort(unique(c(deficient_coords(Kf),
                        vanished_coords(sp, cf, design, Kf))))
        fl <- setdiff(fl, gam$where)
        if (length(fl)) {
          # left free in the inner fit, where the objective still depends on
          # them, and out of the determinant alone: the design carries them,
          # and pinned_coords() reads them wherever K is pinned
          attr(design, "pinned") <<- fl
          m <<- NULL
          body()
        }
      }
    } else {
      err <- tryCatch({
        body()
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(err)) m <- NULL
    }
    state$evals <- state$evals + 1L
    if (is.null(m)) {
      # the half that was missing: put the structural parameters back where
      # the last USABLE point left them, as the coefficients already are
      if (!is.null(sst) && !is.null(z_save)) {
        sst$zeta <- z_save
        sst$key <- NULL
        sst$value <- NULL
      }
      if (vb$outer) {
        # the PARAMETER scale, as the accepted line uses: printing eta here
        # made a prior scale look negative
        hv <- tryCatch(paste(signif(vapply(seq_along(labels), function(j)
          hyper_value(hy, idx, j), numeric(1)), 5), collapse = ", "),
          error = function(e) paste(signif(eta, 5), collapse = ", "))
        cat(sprintf("  %6d %16s %12s   %s%s\n", state$evals, "unavailable",
                    "--", hv,
                    if (is.null(err)) "" else paste0("  [", err, "]")))
      }
      # The intent stated above, completed for the branch that does not
      # raise: at the STARTING point an unavailable criterion is the
      # caller's problem and is reported with its cause. Handing the
      # optimizer a non-finite value instead has it stop with "the
      # objective is not finite at the starting value", which names the
      # point and not the reason -- measured on a gamma model whose inner
      # fit CONVERGED to a degenerate parameter, so the penalized
      # information had no Cholesky factor at any hyperparameter and the
      # message pointed nowhere.
      if (state$evals == 1L) {
        cause <- if (!is.null(err)) {
          paste0(": ", err)
        } else if (!is.null(res) && !isTRUE(res$converged)) {
          ": the inner fit did not converge there"
        } else {
          paste0(": the inner fit converged to a point where the ",
                 "penalized information has no Cholesky factor, as happens ",
                 "when a parameter is degenerate (a coordinate run to a ",
                 "boundary, a direction that the data do not identify",
                 ")")
        }
        stop(paste0("the ", toupper(method@kind), " criterion is ",
                    "unavailable at the starting hyperparameters", cause,
                    ". Fit with outer_criterion = NULL and inspect the ",
                    "coefficients."), call. = FALSE)
      }
      nh <- nrow(idx)
      # A value the search moves away from, and it has to be FINITE. An
      # infinite one is worse than useless to a method that differences its
      # own gradient: the probe lands in the unavailable region, the
      # difference comes back non-finite, the direction is meaningless and the
      # line search stops -- measured on a panel whose outer bfgs, having no
      # analytic gradient over a structural penalty, died at its 73rd
      # evaluation with the search still improving. A barrier strictly worse
      # than everything seen keeps the objective's own scale, so a difference
      # taken across it points back into the feasible region and the step is
      # simply rejected.
      w <- if (is.finite(state$worst)) state$worst else 1
      bar <- w + abs(w) + 1
      out <- list(value = sgn * bar, grad = rep(0, nv),
                  hess = if (exact2) diag(if (pe) 1 else -1, nv, nv) else NULL,
                  ok = FALSE)
      state$key <- key
      state$last <- out
      return(out)
    }
    state$beta <- res$par
    state$inner <- res
    state$hyper <- hy
    vfn <- sgn * m$value
    if (is.finite(vfn) && vfn < state$best_v) {
      state$best_v <- vfn
      state$best_beta <- res$par
      state$best_zeta <- if (is.null(sst)) NULL else sst$zeta
    }
    if (is.finite(vfn)) {
      state$worst <- if (is.finite(state$worst)) max(state$worst, vfn) else vfn
    }
    row <- data.frame(evaluation = state$evals, criterion = m$value,
                      loglik = m$loglik, penalty = m$penalty)
    for (j in seq_along(labels)) {
      row[[labels[j]]] <- v_value(hy, v, j)
    }
    state$rows[[length(state$rows) + 1L]] <- row
    if (vb$outer) {
      # the same columns the inner optimizer prints, and the hyperparameters
      # on ONE scale: the PARAMETER scale, which is what a reader of a fit
      # wants and what the returned object carries. The unavailable line used
      # to print the free scale instead, so a prior scale read as negative --
      # it was log(sigma) under a heading that looked like sigma.
      # the trace reports the gradient the search is using, so it asks for one
      # only where the search has one to use: printing it under a
      # derivative-free method would compute the very quantity nobody wants.
      gv <- if (exact && trace_grad) derivs(1L) else NULL
      gmax <- if (!is.null(gv) && all(is.finite(gv)))
        sprintf("%12.4g", max(abs(gv))) else sprintf("%12s", "--")
      cat(sprintf("  %6d %16.6f %s   %s\n", state$evals, m$value, gmax,
                  paste(signif(vapply(seq_along(labels), function(j)
                    v_value(hy, v, j), numeric(1)), 5), collapse = ", ")))
    }
    # `par`, `split` and `score` are carried for criterion_resolution(), which
    # needs the mode the criterion was read at and the score the inner fit
    # stopped short by. The score costs one gradient of the inner objective
    # against a whole inner fit, so it is not worth deferring.
    out <- list(value = m$value, ok = TRUE, ctx = ctx, cf = cf, hy = hy,
                par = res$par, split = res$obj$split,
                score = free_score(res$obj$gr(res$par)), sp = sp)
    state$key <- key
    state$last <- out
    # RECOMPUTED AT EVERY USABLE POINT rather than once at the start. The mode
    # is least well located at the first evaluation, which is a cold start, so
    # a reading taken there is the worst of the run and governs every step
    # after it. One criterion assembly per evaluation against a whole inner
    # fit is what it costs.
    if (chose_optimizer) {
      r <- criterion_resolution(out, sp, design, method, criterion_at)
      if (is.finite(r) && r > 0) {
        state$res_seen <- c(state$res_seen, r)
        state$resolution <- resolution_summary(state$res_seen)
      } else if (!is.null(attr(r, "mode_error"))) {
        # REFUSED because the inner fit is not at a mode, and that is worth
        # saying: the search then runs with no resolution at all, which is
        # what it did before one existed, and the reason is a property of the
        # fit rather than of the search. Said ONCE -- it is the same fit
        # every evaluation and a line per evaluation would be noise.
        if (vb$outer && !isTRUE(state$said_mode_error)) {
          state$said_mode_error <- TRUE
          vb_say(paste0("no resolution: the inner fit sits %.3g above its own",
                        " minimum, so the displacement that it implies is not a",
                        " correction"), attr(r, "mode_error"), indent = 5L)
        }
      }
    }
    out
  }

  # The derivatives, produced when the search asks and cached beside the point
  # they belong to. `order` is 1 for the gradient and 2 for the Hessian.
  derivs <- function(order) {
    st <- state$last
    slot <- if (order == 1L) "grad" else "hess"
    if (!is.null(st[[slot]])) return(st[[slot]])
    if (!isTRUE(st$ok)) {
      # an unavailable point keeps the barrier's own answers
      return(if (order == 1L) rep(0, nv) else
        if (exact2) diag(if (pe) 1 else -1, nv, nv) else NULL)
    }
    fallback <- function(o) if (o == 1L) rep(0, nv) else
      diag(if (pe) 1 else -1, nv, nv)
    if (pe) {
      # the prediction-error route computes the gradient on the way to the
      # Hessian, so both are kept rather than the second order being asked for
      # twice
      d <- statmod_pe_derivs(spec, design, st$cf, st$hy, method, idx, order,
                             emethod)
      if (is.null(st$grad)) {
        st$grad <- if (is.null(d$grad)) fallback(1L) else d$grad
      }
      if (order >= 2L) st$hess <- if (is.null(d$hess)) fallback(2L) else d$hess
    } else if (order == 1L) {
      v <- statmod_marginal_grad(st$sp, design, st$cf, st$hy, method, idx,
                                 basis, ctx = st$ctx, gam = gam)
      st$grad <- if (is.null(v)) fallback(1L) else v
    } else {
      v <- statmod_marginal_hess(st$sp, design, st$cf, st$hy, method, idx,
                                 basis, ctx = st$ctx, gam = gam)
      st$hess <- if (is.null(v)) fallback(2L) else v
    }
    state$last <- st
    st[[slot]]
  }

  fn <- function(eta) sgn * evaluate(eta)$value
  gr <- function(eta) { evaluate(eta); sgn * derivs(1L) }
  he <- function(eta) { evaluate(eta); sgn * derivs(2L) }

  # optimizers7 checks a caller-supplied gradient against a directional
  # difference of the objective, and on a criterion this curved the step it
  # uses is too long: measured on a variance component, it read a rate of
  # 0.996 where the gradient is 0.4786 and numDeriv agrees with the gradient
  # to 1.3e-4. The guard is for a gradient that may not match the objective;
  # here both come from this package and the agreement is what the tests
  # establish, so it is turned off for this search and for nothing else.
  if (exact) {
    old_check <- getOption("optimizers7.check_gradient")
    options(optimizers7.check_gradient = FALSE)
    on.exit(options(optimizers7.check_gradient = old_check), add = TRUE)
  }

  # The legend, printed ONCE. The names are what identify a hyperparameter and
  # they are long -- a term's key is the call that produced it -- so repeating
  # them on every evaluation is what made a trace of 130 evaluations
  # unreadable. Here they are said once and the lines below carry the values.
  if (vb$outer) {
    vb_rule(sprintf("outer search: %s over %d hyperparameter%s",
                    toupper(method@kind), length(shown),
                    if (length(shown) == 1L) "" else "s"),
            vb_name(optimizer), indent = 2L)
    for (j in seq_along(shown)) {
      vb_say("h%-3d %s", j, shown[j], indent = 5L)
    }
    # the column header the evaluation lines below fill, in the shape the
    # inner optimizer prints its own. `criterion` is what the search
    # maximizes and `|grad|max` its exact gradient where there is one; the
    # hyperparameters are on the PARAMETER scale throughout.
    cat(sprintf("  %6s %16s %12s   %s\n", "eval", toupper(method@kind),
                "|grad|max",
                if (length(shown) == 1L) "h" else
                  paste0("h1..h", length(shown))))
  }

  eta0 <- c(hyper_to_eta(hyper, idx), beta[gam$where])

  # A CORRELATION'S ANGLE THAT RUNS PAST 9 IS SEARCHED AGAIN INSIDE A BOX OF
  # HALF-WIDTH 9, where the package chose the optimizer. The criterion's slope
  # in an angle decays like exp(-2|z|), so a correlation estimated at one runs
  # out along a direction that is flat, and where it stops is platform
  # arithmetic. Past about 10 the covariance is singular to working
  # precision: measured on a shared (mu, sigma) intercept block whose truth is
  # a correlation of one, a search started anywhere between -8.5 and -11 ran
  # to -14.55, ten coefficients of sigma became unidentified, the certificate
  # read `unknown` and the summary lost the standard errors of both standard
  # deviations; Windows stopped at -9.14 and Linux at -9.37 (both `boundary`)
  # or -10.24 (`unknown`). The box sits one unit past the certificate's edge
  # of 8, whose test is strict, so a search that runs out is still read as at
  # the edge, at a correlation within 3e-8 of one.
  #
  # ⚠️ THE BOX IS ENTERED ONLY WHERE THE FREE SEARCH LEAVES IT. A box is a
  # reparametrization, and it changes the search's path everywhere and not
  # only near the edge: applied from the start, the battery's
  # gas-mixed-class, whose free search stops at an angle of -0.74, went to
  # -7.15 and a criterion 1.06 higher but with an indefinite curvature
  # (state `unknown`). Run only after the free search has left (-9, 9), every
  # fit that stays inside is unchanged.
  #
  # An angle is a free value parameters7 names z<i>.<j>, or a hyperparameter
  # on a rhobit link.
  angle <- grepl("(^|_)z[0-9]+\\.[0-9]+$", idx$name) |
    vapply(attr(idx, "links"), function(l) any(grepl("Rhobit", class(l))),
           logical(1))
  box_lo <- rep(-Inf, length(eta0))
  box_up <- rep(Inf, length(eta0))
  box_lo[which(angle)] <- -9
  box_up[which(angle)] <- 9
  has_box <- chose_optimizer && any(angle)
  use_box <- FALSE
  mini <- function(opt, start, ...) {
    args <- list(opt, fn, start, ...)
    if (use_box && isTRUE(optimizers7::optimizer_bounded(opt))) {
      args$lower <- box_lo
      args$upper <- box_up
      out <- which(angle & abs(start[seq_along(angle)]) >= 9)
      args[[3L]][out] <- sign(start[out]) * 8.5
    }
    do.call(optimizers7::minimize, args)
  }
  ran_out <- function(r) {
    has_box && !use_box && any(abs(r@par[which(angle)]) > 9)
  }

  # THE CRITERION HAS A RESOLUTION AND THE STOPPING RULE HAS TO KNOW IT.
  # Every evaluation refits the coefficients from the RUNNING warm start, so
  # the value at one hyperparameter depends on the path taken to it, and a
  # line search cannot verify a decrease smaller than that. Both halves of an
  # optimizer's default rule are then unreachable: `crit_grad()` asks 1e-6 of
  # a gradient that bottoms out where the decrease stops being verifiable, and
  # `crit_rel_obj()` asks 1e-12 of a value of order 1e4, which is an absolute
  # 1e-8. The run reports failure at a point it does not leave -- measured on
  # a t-prior random effect, three consecutive evaluations identical to seven
  # digits with the flag FALSE, while the inner fit at those hyperparameters
  # converges and the log-likelihood agrees to nine digits.
  #
  # The resolution is COMPUTED at this fit rather than declared from the inner
  # tolerance, and the difference is not one of accuracy. Declared as the
  # criterion's scale times the tolerance it cannot be both safe and useful:
  # measured, the resolution divided by that scale ranges over 140 times across
  # shapes at ONE tolerance, so a bound covering a random intercept over 500
  # levels is four orders too large for a gaussian smooth, where it stops the
  # search with the outer gradient still 1.8e-2. `criterion_resolution()`
  # displaces the mode by the error the inner fit's own score implies and reads
  # how far the criterion moves, which costs one assembly and no refit; see
  # there for what it was measured against.
  #
  # It is ADDED to the optimizer's own rule rather than replacing it, so the
  # run can only stop earlier and never later, and it goes to the LINE SEARCH
  # as well, the two covering different exits. The rule is read between
  # iterations; a search that exhausts its backtracking budget never gets
  # there, the loop asking the rule with no previous iterate, where anything
  # reading a CHANGE in the objective returns FALSE by construction. Told what
  # the criterion can resolve, the search returns at the first step whose
  # predicted improvement is below it instead of paying for thirty trials,
  # each of them a whole inner fit.
  # the optimizer's own rule, before the resolution is added to it: the
  # restart below reads its own unit and not this one
  own_criterion <- optimizer@criterion
  if (chose_optimizer) {
    fn(eta0)
    resolution <- state$resolution
    # ⚠️ NOT FOR A DERIVATIVE-FREE SEARCH. A simplex or a compass search
    # reports the best vertex, which stays where it is over every iteration
    # that only reshapes the simplex, so a rule reading a CHANGE in the
    # objective fires at the first such iteration whatever the resolution.
    # Measured on a multivariate t prior over the random slopes of
    # `sleepstudy`, the default nelder_mead() stopped after 10 evaluations at a
    # criterion of -898.17, two iterations from its start, where the same
    # search without the rule runs 773 evaluations. Such a method declares
    # "stationarity" rather than "gradient", which is what is read here.
    derivative_free <- identical(optimizers7::optimizer_provides(optimizer),
                                 "stationarity")
    if (is.finite(resolution) && resolution > 0 && !derivative_free) {
      # the RULE keeps the reading from the starting point, being a property of
      # the optimizer object and settled before the run; only the line search,
      # which is asked again at every iteration, follows the running summary
      optimizer <- S7::set_props(
        optimizer,
        criterion = optimizers7::crit_any(
          optimizer@criterion, optimizers7::crit_abs_obj(resolution)))
    }
    # ⚠️ THE LINE SEARCH IS GIVEN A TENTH OF IT, and the two reasons both point
    # the same way. What is computed is the spread between evaluations reached
    # from DIFFERENT warm starts, while a line search compares trials taken one
    # after another from nearly the same state, where the criterion is far more
    # reproducible -- repeated at one hyperparameter from a running warm start
    # its spread is exactly zero. And the reading is the more aggressive
    # consumer, ending a whole iteration rather than one comparison. Measured
    # at the full number on a gaussian smooth: the search stops after 3
    # evaluations against 29 and gives up 3.2e-06 of criterion against a
    # resolution of 9.7e-07, so it stops just before it has to; at a tenth it
    # reaches the same optimum as a search that was told nothing.
    #
    # A derivative-free search has no line search, so the property is asked for
    # rather than assumed.
    #
    # ⚠️ AND IT IS GIVEN ITS CLOSURE WHETHER OR NOT THE FIRST READING WAS
    # USABLE, which is the whole difference between the two consumers.
    # The criterion's rule is a number settled before the run, so it can only
    # be built from the reading at the starting point -- a cold start, the
    # worst-located mode of the whole fit. The line search reads a CLOSURE over
    # resolution_summary() of every reading so far, so it needs no reading now
    # and picks up whatever the search refines later.
    #
    # Nesting it inside the test above threw that away: measured on
    # `seg(x, psi ~ random(~1|id))`, whose first reading is refused because the
    # cold-started mode sits 0.046 above its own minimum, the fit went from 5
    # evaluations to 18 and from converged to not -- at an identical answer,
    # cor 0.9932 and rmse 0.0674 either way. Outside the test the closure
    # serves that fit from the good readings that follow.
    #
    # Where no reading is ever usable the closure returns NA, which
    # optimizers7 reads as no resolution at all -- which is the answer for a
    # fit whose mode is never located.
    #
    if ("line_search" %in% S7::prop_names(optimizer)) {
      optimizer <- S7::set_props(
        optimizer,
        line_search = S7::set_props(
          optimizer@line_search,
          resolution = function() state$resolution / 10))
    }
    # THE BACKTRACKING BUDGET, and it is a COST decision rather than a
    # numerical one, which is why it sits outside the resolution's block and
    # is set whether or not one could be computed.
    #
    # optimizers7 defaults to 30 backtracks, which is right for an objective
    # costing microseconds. Here a trial is a penalized refit, and a line
    # search that is going to fail spends the whole budget discovering it --
    # the traces that started this work show two blocks of 31 identical
    # evaluations, which is that budget burned twice.
    #
    # ⚠️ WHAT IT BUYS IS SMALLER THAN THE EVALUATION COUNT SUGGESTS, and the
    # reason is in evaluate() above: every trial of a line search warm-starts
    # from the INCUMBENT, the best usable point so far, and as the step
    # shrinks the trial begins at very nearly its own answer. Measured, removing 22 of 38 evaluations
    # removed 2.8 s of 30.8 -- 0.13 s each against an average evaluation's
    # 0.81. An evaluation count is not a cost here.
    #
    # ⚠️ AND IT MUST BE MEASURED ON THIS PATH AND NOT ANOTHER. Swept with the
    # optimizer NAMED -- which turns the resolution above off -- a short budget
    # flips the convergence flag from TRUE to FALSE on healthy shapes, and that
    # would have been a bad trade. It does not happen here: measured at 30, 12
    # and 8 backtracks on a smooth, two smooths with a random effect, and a
    # random intercept, all three are unchanged in evaluations, criterion to
    # six decimals, effective degrees of freedom and flag, the resolution rule
    # stopping the search before the budget is ever reached. What moves is the
    # expensive shape: a hierarchical break-point model goes from 31
    # evaluations and 25.6 s to 13 and 20.3, with the criterion 1.3e-04 BETTER
    # and the edf identical.
    #
    # 12 rather than 8 because 12 takes 5.3 s of the 6.6 available and 8 takes
    # the rest at a criterion the named-optimizer sweep shows starting to
    # degrade: against a budget of 30 the criterion given up there is 4.8e-05
    # at 12, 7.8e-04 at 8 and 3.3e-03 at 6, and the last is past the 1.6e-03
    # that criterion can resolve.
    if ("line_search" %in% S7::prop_names(optimizer)) {
      optimizer <- S7::set_props(
        optimizer,
        line_search = S7::set_props(optimizer@line_search,
                                    max_step = outer_backtracks()))
    }
    # A TYPICAL SIZE PER COORDINATE where the criterion also estimates
    # coefficients. newton() reads its lengths in the parameters' own units --
    # the bound of 5, and the length-one scaling of a direction it had to
    # repair -- and a hyperparameter's free value carries no units while a
    # coefficient carries those of its parameter and of its covariate.
    # Measured on `dist ~ speed` with gaussian2 on the identity link, started
    # at a variance of 651 where the restricted likelihood is not concave in
    # it, the repaired steps moved the variance by one unit per evaluation,
    # 185 evaluations in all; a standard deviation of 1538 on the identity
    # link took 201 and stopped 0.9 per cent short. A coefficient's typical
    # size is its conditional standard error at the start, where that exceeds
    # one, which is the same size in any units; the Hessian is the one the
    # first iteration reads anyway.
    if (ng && exact2 && "typical" %in% S7::prop_names(optimizer) &&
        is.null(optimizer@typical)) {
      H0 <- tryCatch(he(eta0), error = function(e) NULL)
      if (!is.null(H0)) {
        se0 <- 1 / sqrt(abs(diag(as.matrix(H0)))[nh + seq_len(ng)])
        se0[!is.finite(se0)] <- 1
        if (any(se0 > 1)) {
          optimizer <- S7::set_props(
            optimizer, typical = c(rep(1, nh), pmax(1, se0)))
        }
      }
    }
  }

  res <- if (exact2) {
    mini(optimizer, eta0, gr = gr, he = he)
  } else if (exact) {
    mini(optimizer, eta0, gr = gr)
  } else {
    mini(optimizer, eta0)
  }

  if (ran_out(res)) {
    use_box <- TRUE
    res_box <- tryCatch(
      if (exact2) {
        mini(optimizer, res@par, gr = gr, he = he)
      } else if (exact) {
        mini(optimizer, res@par, gr = gr)
      } else {
        mini(optimizer, res@par)
      }, error = function(e) NULL)
    if (!is.null(res_box) && is.finite(res_box@value)) res <- res_box
  }

  # A LINE SEARCH THAT FAILED AT THE ROUNDING OF THE CRITERION IS A
  # RESOLUTION STOP. Where no reading of the criterion's resolution exists (a
  # fit with no penalty, whose hyperparameters are a dispersion's
  # coefficients) or the reading is smaller than the arithmetic, a Newton
  # step at the optimum predicts a decrease the criterion cannot show, every
  # trial is rejected on the last bits, and the run reports no acceptable
  # step at the optimum. Measured on y ~ x | sigma ~ z on arm64 macOS: twelve
  # backtracks and the flag FALSE where x86_64 accepted the first trial by a
  # few ulps. The test is made AFTER the search has failed, so a run that
  # succeeds is untouched: a floor read before the trials stopped ordinary
  # gaussian fits at a relative 4e-7 in sigma. The floor is 1e3 eps |f|, and
  # a wrong direction, which predicts a large decrease, stays a failure.
  if (chose_optimizer && exact2 && !isTRUE(res@converged) &&
      identical(res@message, "the line search found no acceptable step")) {
    floor_res <- 1e3 * .Machine$double.eps *
      max(1, if (is.finite(res@value)) abs(res@value) else 1)
    at_floor <- isTRUE(tryCatch({
      evaluate(res@par)
      g <- derivs(1L)
      A <- sgn * derivs(2L)
      dec <- joint_decrement(g, (A + t(A)) / 2)
      is.finite(dec) && dec <= floor_res
    }, error = function(e) FALSE))
    if (at_floor) {
      res <- S7::set_props(
        res, converged = TRUE,
        criterion_met = "no decrease above the objective's resolution",
        message = "")
    }
  }

  # WHETHER THE REPORTED POINT IS THE OPTIMUM, which the search's flag does
  # not say: the flag says whether its stopping rule fired. Beside a kinked
  # block the criterion carries the coordinate descent's own rounding, of
  # order 1e-8, which criterion_resolution() does not see, and a line search
  # asked for a smaller decrease runs out of backtracks at the optimum. The
  # Newton decrement at the point answers the question in the criterion's
  # units, as statmod_certificate() does, and a path reads it to decide
  # whether a point is scored -- the rule inner_mode_error() applies to an
  # inner fit, one level up. Measured on y ~ s(x) + lasso(...) over six
  # samples, the path lost 3 to 11 of its 25 points to that message.
  settled <- isTRUE(res@converged)
  if (!settled && exact2) {
    settled <- isTRUE(tryCatch({
      evaluate(res@par)
      g <- derivs(1L)
      A <- sgn * derivs(2L)
      dec <- joint_decrement(g, (A + t(A)) / 2)
      is.finite(dec) && dec <= eval(formals(statmod_certificate)$tol)
    }, error = function(e) FALSE))
  }

  # A NEWTON SEARCH THIS PACKAGE CHOSE AND THAT DID NOT SETTLE IS RUN AGAIN BY
  # LBFGS, from the same start, and the better of the two is kept.
  # Where the criterion's curvature is indefinite the Newton direction is a
  # repaired one and the search can crawl: on GAGurine with a Student t and
  # nu ~ s(Age) it stopped at a REML criterion of -845.09 with the curvature
  # not negative definite, where lbfgs() and nelder_mead() reach -816.59 --
  # the smoothing parameter of nu at its edge, so nu constant in age.
  if (chose_optimizer && !settled && exact &&
      inherits(optimizer, "optimizers7::Newton")) {
    # with its own stopping rule and line search: the measured search that
    # reaches the maximum is lbfgs() as a caller names it, and the budget and
    # resolution this function gives its own choice stopped it at -818.63
    alt <- optimizers7::lbfgs()
    res2 <- tryCatch(mini(alt, eta0, gr = gr),
                     error = function(e) NULL)
    if (!is.null(res2) && is.finite(res2@value) &&
        res2@value < res@value - 1e-8 * max(1, abs(res@value))) {
      res <- res2
      optimizer <- alt
      settled <- isTRUE(res@converged)
      if (!settled && exact2) {
        settled <- isTRUE(tryCatch({
          evaluate(res@par)
          g <- derivs(1L)
          A <- sgn * derivs(2L)
          dec <- joint_decrement(g, (A + t(A)) / 2)
          is.finite(dec) && dec <= eval(formals(statmod_certificate)$tol)
        }, error = function(e) FALSE))
      }
    }
  }

  # A SEARCH THIS PACKAGE CHOSE THAT LEFT ON THE RESOLUTION RULE WITH THE
  # CRITERION STILL RISING IS RUN AGAIN FROM ITS POINT, with an empty
  # quasi-Newton memory. The rule reads the decrease the search's own model
  # predicts, and a model fed one bad curvature pair predicts nothing: on the
  # battery's gas-panel-omega an unconverged inner point read -654.9 entered
  # lbfgs()'s memory and the search left after three iterations with the
  # gradient at 1.2 and the criterion at -609.528, where it rises to
  # -609.490. The Newton decrement at the point is what decides, in the
  # criterion's units, against mode_error_limit(); the certificate's 0.01
  # left that run 0.0028 short. At most three restarts.
  #
  # ⚠️ THE RESTART STOPS IN THE DECREMENT'S UNIT AND NOT IN THE RESOLUTION'S.
  # It is run only where the exact decrement exceeds mode_error_limit(), so
  # its rule is a gain below that limit and its line search is given a tenth
  # of it. Given the search's own rule instead, it stopped wherever the first
  # step's predicted gain fell against a tenth of the resolution read at the
  # start, which on gas-panel-omega is about 2 REML units and moves by 45 per
  # cent with a perturbation of one ulp in the rhobit link's derivatives: at
  # one reading the step went through (0.217 against 0.199) and the criterion
  # reached -609.4905, at the other it did not (0.217 against 0.287) and the
  # fit stopped at -609.4932 with a decrement of 2.6e-3. With its own unit the
  # restart reaches -609.490415 from either, with the decrement at 1e-6.
  restart_opt <- optimizer
  if (chose_optimizer && "line_search" %in% S7::prop_names(optimizer)) {
    restart_opt <- S7::set_props(
      optimizer,
      criterion = optimizers7::crit_any(
        own_criterion, optimizers7::crit_abs_obj(mode_error_limit())),
      line_search = S7::set_props(
        optimizer@line_search,
        resolution = function() mode_error_limit() / 10))
  }
  for (attempt in seq_len(3L)) {
    if (!(chose_optimizer && exact2 && isTRUE(res@converged) &&
          (attempt > 1L ||
           identical(res@criterion_met,
                     "no decrease above the objective's resolution")))) break
    rising <- isTRUE(tryCatch({
      evaluate(res@par)
      g <- derivs(1L)
      A <- sgn * derivs(2L)
      dec <- joint_decrement(g, (A + t(A)) / 2)
      is.finite(dec) && dec > mode_error_limit()
    }, error = function(e) FALSE))
    if (!rising) break
    res2 <- tryCatch(mini(restart_opt, res@par, gr = gr, he = he),
                     error = function(e) NULL)
    if (is.null(res2) || !is.finite(res2@value) ||
        !(res2@value < res@value - mode_error_limit())) {
      tryCatch(evaluate(res@par), error = function(e) NULL)
      break
    }
    res <- res2
    settled <- isTRUE(res@converged)
  }

  # A HYPERPARAMETER AT THE EDGE OF ITS CHART IS PROBED FROM INSIDE. On the
  # free scale the criterion's slope carries the chart's own derivative, which
  # vanishes at the edge, so a search can stop there with a vanishing gradient
  # at a point that is not a maximum. Measured on nlme::Machines with
  # random(~1 | Worker) + random(~1 | Worker:Machine), the default search
  # stopped with the Worker sd at 0.00025 and a REML criterion of -110.63,
  # where lme() reaches -107.84; with the Worker sd held, the criterion rises
  # steadily from 0 to 4.78. The probe moves each coordinate past `edge` back
  # towards its starting value, at four fractions of the distance, with the
  # rest held, and restarts the search from the best probe where it gains more
  # than mode_error_limit(). Four fractions and not one: an AR(1) correlation
  # that ran to rho = 1 gains nothing at rho = 0.985 and loses 23 at its start
  # of zero, and gains 0.21 at rho = 0.84. See edge_probe_hyper().
  if (chose_optimizer && nh) {
    for (attempt in seq_len(3L)) {
      # a probe that gains less than the threshold is not taken, and must
      # leave no trace: evaluate() makes any better point the incumbent the
      # final refit warm-starts from, and on a model with a multimodal inner
      # problem that moved the convergence flag of a fit whose point and
      # criterion were unchanged
      keep <- mget(c("best_v", "best_beta", "best_zeta", "beta"), envir = state)
      restore <- function() {
        for (k in names(keep)) assign(k, keep[[k]], envir = state)
        # and the point the search reported is evaluated once more, so the
        # design's own state (the blocks refreshed at the last point read)
        # is that point's and the history ends on it, as it did before
        tryCatch(evaluate(res@par), error = function(e) NULL)
      }
      pr <- edge_probe_hyper(res@par, eta0, nh, res@value, fn)
      if (is.null(pr)) {
        restore()
        break
      }
      res2 <- tryCatch(
        if (exact2) {
          mini(optimizer, pr, gr = gr, he = he)
        } else if (exact) {
          mini(optimizer, pr, gr = gr)
        } else {
          mini(optimizer, pr)
        }, error = function(e) NULL)
      if (is.null(res2) || !is.finite(res2@value) ||
          !(res2@value < res@value - mode_error_limit())) {
        restore()
        break
      }
      res <- res2
      settled <- isTRUE(res@converged)
      if (!settled && exact2) {
        settled <- isTRUE(tryCatch({
          evaluate(res@par)
          g <- derivs(1L)
          A <- sgn * derivs(2L)
          dec <- joint_decrement(g, (A + t(A)) / 2)
          is.finite(dec) && dec <= eval(formals(statmod_certificate)$tol)
        }, error = function(e) FALSE))
      }
    }
  }

  # the last evaluation is not necessarily the optimum, so the fit is taken at
  # the reported point rather than at whatever was tried last
  hy <- eta_to_hyper(res@par[seq_len(nh)], idx, hyper)
  sp <- spec_at(res@par)
  # from the incumbent, as every trial was: the reported point is the one the
  # search stands on, and refitting it from a rejected trial's state is what
  # reported a worse mode than the search had found
  warm <- state$beta
  if (!is.null(state$best_beta)) {
    warm <- state$best_beta
    sst <- statmod_structural_state(design)
    if (!is.null(sst) && !is.null(state$best_zeta)) {
      sst$zeta <- state$best_zeta
      sst$key <- NULL
      sst$value <- NULL
    }
  }
  if (ng) warm[gam$where] <- res@par[nh + seq_len(ng)]
  inner <- statmod_alternate(sp, design, blocks, hy, inner_optimizer,
                             warm, expected, emethod, maxit, tol, vb,
                             hold_refresh = TRUE)
  cff <- inner$obj$split(inner$par)
  # a scaled SCAD or MCP is read at the curvature the refit settled at
  if (!is.null(inner$spec)) sp <- inner$spec
  m <- if (pe) statmod_pe(sp, design, cff, hy, method, emethod,
                          statmod_active(sp, blocks, inner$par, hy)) else
    statmod_marginal(sp, design, cff, hy, method, emethod, basis)
  list(par = inner$par, hyper = hy, value = inner$value,
       criterion = if (is.null(m)) NA_real_ else m$value,
       converged = res@converged && inner$converged,
       settled = settled && inner$converged,
       obj = inner$obj, hist_blocks = inner$hist_blocks,
       hist_inner = inner$hist_inner,
       # the columns the final refit's pivot left out, which statmod() reports
       # as not identified; without them a fit that went through the search
       # was read only on the information at the mode, where a column the
       # pivot drops at its tolerance can still be inverted
       aliased = inner$aliased, spec = inner$spec,
       hist_outer = if (length(state$rows)) do.call(rbind, state$rows) else
         NULL,
       iterations = res@iterations, evaluations = state$evals,
       exact_gradient = exact, exact_hessian = exact2,
       # the optimizer that actually searched, which is the caller's where
       # one was given and otherwise the one chosen from what the criterion
       # can supply: a fit says what fitted it rather than leaving a reader
       # to reconstruct the default
       optimizer = optimizer,
       # the coefficients the criterion estimated, which a refit after the
       # search has to go on holding: see statmod()
       held_coef = if (ng) sp@held_coef else list(),
       marginal = gam,
       pinned = as.integer(attr(design, "pinned")))
}


#' Probe the Hyperparameters at the Edge of Their Chart From Inside
#'
#' @description
#' Moves each hyperparameter whose free value has run past `edge` back
#' towards its starting value while the other coordinates are held, and returns
#' the free vector of the best probe where it improves the criterion by more
#' than `tol`.
#'
#' @details
#' A hyperparameter lives on a chart such as \eqn{\sigma = e^\eta} or
#' \eqn{\rho = \tanh(\eta)}, and the slope of the criterion in \eqn{\eta}
#' is its slope in the bounded quantity times \eqn{d\sigma/d\eta} or
#' \eqn{d\rho/d\eta}, which tend to zero at the edge. A search can therefore
#' stop at the edge with a vanishing gradient where the criterion still rises
#' towards the inside. A finite step reads that rise where a derivative
#' cannot. Each coordinate past the edge is probed at the points
#' \eqn{s + f(\hat\eta - s)}, with \eqn{s} its starting value (or zero where
#' the start itself is past the edge), \eqn{\hat\eta} the reported value and
#' \eqn{f} equal to three quarters, one half, one quarter and zero, so that the
#' last probe is the starting value itself. A coordinate that sits at an edge,
#' such as the smoothing parameter of a smooth of noise, loses criterion at
#' every probe and is left where it is. The cost is four evaluations of the
#' criterion per coordinate at the edge, one inner fit each.
#'
#' @param par The point that the search reported: the free values of the
#'   hyperparameters first, then any coefficients that the criterion
#'   estimates.
#' @param start The starting point of the search, of the same length.
#' @param nh The number of hyperparameters at the head of `par`.
#' @param value The value of the search's objective at `par`, which it
#'   minimizes.
#' @param fn The search's objective.
#' @param edge The free value past which a coordinate is at the edge, the
#'   default of [statmod_certificate()].
#' @param tol The gain, in the criterion's units, that counts:
#'   [mode_error_limit()].
#'
#' @return The free vector of the best probe, or `NULL` if all coordinates lie
#'   inside the edge or every probe gains at most `tol`.
#'
#' @keywords internal
edge_probe_hyper <- function(par, start, nh, value, fn, edge = 8,
                             tol = mode_error_limit()) {
  at_edge <- which(abs(par[seq_len(nh)]) > edge)
  if (!length(at_edge) || !is.finite(value)) return(NULL)
  best <- NULL
  best_v <- value - tol
  for (j in at_edge) {
    s <- if (abs(start[[j]]) <= edge) start[[j]] else 0
    for (f in c(0.75, 0.5, 0.25, 0)) {
      v <- par
      v[[j]] <- s + f * (par[[j]] - s)
      fv <- tryCatch(fn(v), error = function(e) NA_real_)
      if (is.finite(fv) && fv < best_v) {
        best_v <- fv
        best <- v
      }
    }
  }
  best
}


#' One Hyperparameter of an Index
#'
#' @description
#' Returns the value to which the `j`-th row of an outer index points.
#'
#' @param hyper The hyperparameter structure.
#' @param idx The index.
#' @param j A row.
#'
#' @return A single number.
#'
#' @keywords internal
hyper_value <- function(hyper, idx, j) {
  as.numeric(hyper[[idx$parameter[j]]][[idx$term[j]]][[idx$name[j]]])[1L]
}


#' The Verbosity of an Inner Fit Inside the Outer Search
#'
#' @description
#' Returns the same switches, but keeps the block trace only where the inner
#' switch is also on, since one line per pass per hyperparameter tried is too
#' long a trace for the outer search.
#'
#' @param vb The resolved verbosity.
#'
#' @return A list of switches.
#'
#' @keywords internal
vb_inner <- function(vb) {
  vb$blocks <- vb$blocks && vb$inner
  vb
}


#' The Outer Line Search's Backtracking Budget
#'
#' @description
#' Returns the number of backtracks that a line search inside [outer_fit()]
#' may take before it stops, where the optimizer is the one that
#' [outer_default_optimizer()] selects.
#'
#' @details
#' \pkg{optimizers7} defaults to 30, which suits an objective that costs
#' microseconds. Every trial here is a penalized refit, and a line search that
#' fails spends its whole budget before it stops.
#'
#' Every trial warm-starts from the incumbent, the best usable point found so
#' far, so as the step shrinks the refit begins near its own answer and a
#' rejected backtrack costs less than an average evaluation. The budget applies
#' only to the line search of an optimizer selected by
#' [outer_default_optimizer()]; the resolution that [criterion_resolution()]
#' supplies ends the search before the budget on most shapes.
#'
#' An optimizer that the caller named keeps its own budget, as it keeps its own
#' stopping rule.
#'
#' @return A single number, 12.
#'
#' @seealso [outer_fit()], [criterion_resolution()]
#'
#' @keywords internal
outer_backtracks <- function() 12


#' Which Optimizer the Outer Search Uses When the Caller Names None
#'
#' @description
#' The choice is made from what the criterion can supply: its exact Hessian,
#' its exact gradient, or neither.
#'
#' A covariance class split between the coefficients and a filter's own
#' parameters takes [optimizers7::lbfgs()] even where the exact gradient is
#' unavailable, because such a class carries at least three hyperparameters
#' (two scales and a correlation) and a simplex search is unreliable at that
#' dimension.
#'
#' [optimizers7::newton()] is given `max_length = 5`, which bounds every
#' component of a step on the free scale, as `maxNstep` does in mgcv. Where the
#' criterion flattens as a smoothing parameter grows, the curvature at the
#' start can be small enough for the Newton step to jump past the maximum onto
#' that plateau, where the gradient is close to zero and the search stops.
#'
#' A search over one coordinate with no exact gradient takes
#' [optimizers7::brent()], which brackets the maximum and reduces the bracket
#' by parabolic steps, stopping when its half-width on the free scale is
#' below \eqn{10^{-3}}, a change of one part in a thousand in the
#' hyperparameter. On a single coordinate it needs fewer evaluations than
#' [optimizers7::nelder_mead()].
#'
#' @param exact Whether the criterion has an exact gradient.
#' @param use_hess Whether the search should steer by the exact Hessian, which
#'   is not the same question as whether one exists: see [outer_newton_ok()].
#' @param mixed Whether the model carries a covariance class spanning the
#'   coefficients and a structural term's own parameters.
#' @param n The number of coordinates that the search moves.
#'
#' @return An \pkg{optimizers7} optimizer.
#'
#' @seealso [outer_fit()], [outer_gradient_ok()]
#'
#' @keywords internal
outer_default_optimizer <- function(exact, use_hess, mixed = FALSE,
                                    n = NA_integer_) {
  # the step is bounded at 5 on the free scale, mgcv's maxNstep: a Newton
  # step set by a small curvature can otherwise jump past the maximum onto
  # the plateau where a smoothing parameter runs to infinity. See the
  # details above.
  if (use_hess) return(optimizers7::newton(max_length = 5))
  if (exact || mixed) return(optimizers7::lbfgs())
  # one coordinate and no gradient: a bracket reduced by parabolic steps.
  # See the details above.
  if (isTRUE(n == 1L)) return(optimizers7::brent(criterion =
                                optimizers7::crit_stationary(1e-3)))
  optimizers7::nelder_mead()
}


#' Whether the Default Search Steers by the Exact Hessian
#'
#' @description
#' Returns whether a search that [outer_fit()] selects by default should take
#' [outer_gradient_ok()]'s order-2 answer as its direction, which is a narrower
#' question than whether the Hessian exists.
#'
#' @details
#' The two are separated because a criterion may have an exact second
#' derivative that is a poor direction for the search. Over a covariance class
#' the Hessian of the criterion is strongly indefinite wherever the angle of
#' the chart approaches its boundary, so that [optimizers7::newton()] repairs
#' an indefinite Hessian and the repaired step can end elsewhere than a
#' quasi-Newton search. The default therefore uses [optimizers7::lbfgs()] for a
#' covariance class, whose Hessian is exact but not a reliable direction.
#'
#' The function returns `FALSE` for the expected information, for a model with
#' a structural term, and for a penalty with a covariance block wider than one.
#' In the first two cases the cost of the Hessian at each evaluation outweighs
#' the saving in evaluations.
#'
#' The availability of the Hessian still governs the rest: it is supplied to an
#' optimizer that the caller names, `statmod_hyper_vcov()` and
#' `vcov(vcov_type = unconditional_vcov())` read it, and it is recorded as
#' `exact_hessian`. Only the default direction is held back.
#'
#' The block width is read from the `block` property of the penalty where it
#' has one; a penalty without that property is treated as univariate.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, from [statmod_design()].
#' @param exact2 Whether the criterion has an exact Hessian at all.
#' @param method An [OuterMethod()], or `NULL`. Where its `hessian` is
#'   `"expected"` the function returns `FALSE`.
#'
#' @return `TRUE` or `FALSE`.
#'
#' @seealso [outer_default_optimizer()], [outer_gradient_ok()]
#'
#' @keywords internal
outer_newton_ok <- function(spec, design, exact2, method = NULL) {
  if (!isTRUE(exact2)) return(FALSE)
  # THE EXPECTED ROUTE keeps lbfgs(), measured before its order 2 was opened
  # and decided with that measurement (piano_derivate_esatte.txt, D5). With the
  # exact Hessian emulated by the stencil and the search chosen by this
  # package, newton() reaches the criterion lbfgs() reaches to 3e-6 on six
  # models and is slower or level in wall time on all six, each of its
  # evaluations paying a Hessian: 0.39 against 0.34 s on a gaussian smooth,
  # 0.41 against 0.25 on a Poisson one, 1.11 against 0.99 on a gamma, 0.80
  # against 0.69 on a negative binomial, 1.26 against 1.25 on a beta, 4.50
  # against 3.44 on three smooths beside a random effect.
  if (!is.null(method) && info_expected(method@hessian)) return(FALSE)
  # A MODEL CARRYING A STRUCTURAL TERM keeps lbfgs(), and that is a
  # measurement rather than a limitation. Its Hessian exists now --
  # statmod_structural_hess() assembles it on the joint vector -- and each
  # evaluation of it costs one term_fourth() per PAIR of hyperparameters,
  # which is the dominant term. Measured over five panels, every route
  # reaches the same criterion to six decimals and newton() takes fewer
  # evaluations on three of them, yet is slower in wall time on ALL five:
  # 8 evaluations in 1.54 s against 5 in 2.05, 5 in 1.11 against 5 in 2.92,
  # 45 in 5.88 against 11 in 10.62, 7 in 2.38 against 9 in 11.10, and 69 in
  # 3.02 against 45 in 4.06. A Hessian at two hyperparameters costs 0.805 s
  # where a criterion evaluation costs 0.34, which accounts for the gap
  # exactly.
  #
  # What the Hessian's availability still governs is untouched, exactly as
  # for the covariance class below: it is supplied to an optimizer the
  # caller NAMES, it is what statmod_hyper_vcov() and
  # vcov(vcov_type = unconditional_vcov()) read, and it is recorded as `exact_hessian`.
  if (length(attr(design, "structural"))) return(FALSE)
  units <- tryCatch(statmod_penalized(spec, design),
                    error = function(e) list())
  for (u in units) {
    pen <- u$penalty
    if (is.null(pen)) next
    b <- if ("block" %in% S7::prop_names(pen)) pen@block else 1L
    if (isTRUE(as.integer(b) > 1L)) return(FALSE)
  }
  TRUE
}
