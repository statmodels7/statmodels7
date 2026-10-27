#' @include predict.R inference_objects.R
#' @importFrom stats vcov confint
NULL

# What a fitted model reports about its own uncertainty.
#
# The coefficients live on the LINK scale by construction -- they are the
# coefficients of a linear predictor, free on the whole line -- so a symmetric
# Wald interval is built where the quantity actually ranges and needs no
# mapping back. That is the one place this differs from
# distributions7::fit_distrib(), which estimates a bounded parameter and has to
# carry its interval through the link.

#' Where Each Stacked Coefficient Comes From
#'
#' @description
#' Returns one row per stacked coefficient. Each row names the distribution
#' parameter, the term and the coefficient, and records whether the term carries
#' a penalty and whether that penalty has a kink.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#'
#' @return A data frame with as many rows as there are coefficients.
#'
#' @keywords internal
coef_labels <- function(spec, design) {
  params <- spec@distrib@params
  units <- statmod_penalized(spec, design)
  rows <- list()
  for (p in params) {
    d <- design[[p]]
    if (d$npar == 0L) next
    term <- rep(NA_character_, d$npar)
    pen <- rep(FALSE, d$npar)
    kink <- rep(FALSE, d$npar)
    for (nm in names(d$blocks)) term[d$blocks[[nm]]] <- nm
    for (u in units) {
      # a class is considered under EVERY equation one of its members sits
      # in; `param` on it is a convention for the hyperparameter store
      here <- if (is.null(u$pieces)) identical(u$param, p) else
        any(vapply(u$pieces, function(z) identical(z$param, p), TRUE))
      if (!here) next
      # a structural unit's penalty sits on the term's OWN parameters,
      # which contribute no column: its cols index the term's parameter
      # vector, and writing them here grew pen past the design and
      # recycled the labels into duplicate rows
      if (isTRUE(u$structural)) next
      # a covariance class has no columns of its own: its members do, each in
      # its own equation, so the marks go on the pieces that live in this one.
      # A member INSIDE a structural term has none either, and carries `cols`
      # NULL for that reason, so a class split between the two marks its
      # coefficient half here and nothing for the other -- writing the
      # structural half's parameter positions would mark whichever design
      # columns they happen to number.
      cs <- if (is.null(u$pieces)) u$cols else
        unlist(lapply(Filter(function(z) identical(z$param, p), u$pieces),
                      function(z) z$cols), use.names = FALSE)
      pen[cs] <- TRUE
      kink[cs] <- penalty_has_kink(u$penalty)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = p, term = term, coefficient = d$coef_names,
      penalized = pen, kinked = kink, stringsAsFactors = FALSE)
  }
  if (!length(rows)) {
    return(data.frame(parameter = character(0), term = character(0),
                      coefficient = character(0), penalized = logical(0),
                      kinked = logical(0), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- paste(out$parameter, out$coefficient, sep = ":")
  out
}


#' Name the Aliased Coefficients of a Fit
#'
#' @description
#' Converts the coordinates that the pivot of the scoring step dropped into
#' coefficient labels.
#'
#' @details
#' A design of less than full rank has no estimate for the columns that
#' repeat information already carried: only combinations are estimable, and
#' the pivot determines which column is left out, as in `lm()` and `glm()`.
#' Those columns come back from the solve as coordinates of the coefficient
#' vector; this function converts them into the labels that [coef_labels()]
#' supplies, so that a later step can identify a coefficient without rebuilding
#' a design.
#'
#' A coordinate that a penalty covers is never among them, because the pivot
#' runs on the augmented system (design and penalty factor together), in which
#' a column that the design alone does not identify is identified. Under
#' `ridge(~ 0 + x + v)` with two identical columns, the effect is split evenly
#' between the columns and the variance matrix is complete.
#'
#' @param spec A `statmod_spec`.
#' @param design The design, as [statmod_design()] builds it.
#' @param idx Integer coordinates of the coefficient vector, as the
#'   alternation reports them. `NULL` or empty gives `character(0)`.
#'
#' @return A character vector of coefficient labels, possibly empty.
#'
#' @seealso [coef_labels()], which supplies the names, and [vcov()], which
#'   reports these coordinates as missing and inverts the rest of the matrix.
#' @keywords internal
aliased_labels <- function(spec, design, idx) {
  if (is.null(idx) || !length(idx)) return(character(0))
  nm <- rownames(coef_labels(spec, design))
  idx <- idx[idx >= 1L & idx <= length(nm)]
  if (!length(idx)) return(character(0))
  nm[sort(unique(as.integer(idx)))]
}


#' @title The Variance Matrix of a Fit
#' @name vcov.StatmodFit
#' @description
#' Returns the variance of the estimated coefficients, over every distribution
#' parameter's block at once.
#' @details
#' Writing \eqn{H} for the information of the log-likelihood and \eqn{S} for
#' the second derivative of the penalty,
#' \deqn{V_b = (H + S)^{-1}, \qquad V_f = (H+S)^{-1} H (H+S)^{-1}.}
#' \eqn{V_b} is the posterior variance under the prior whose negative
#' logarithm is the penalty, and it is the matrix from which an interval around
#' a penalized term should be built: it carries the smoothing bias as variance,
#' and for that reason such intervals cover at about their nominal rate.
#' \eqn{V_f} is the sampling variance of the penalized estimator at a fixed
#' penalty, which is smaller and covers less. With no penalty \eqn{S = 0} and
#' both are \eqn{H^{-1}}, for a coefficient that the outer criterion does not
#' estimate. Under [reml()] or [ml()] with the default `marginal = NULL`, in a
#' model without a structural term, the unpenalized coefficients of every
#' distribution parameter other than the first (for example the log standard
#' deviation of a Gaussian) are estimated on the criterion, and [bayesian_vcov()]
#' and [unconditional_vcov()] take the curvature of that criterion for them.
#' These two conventions then differ from [frequentist_vcov()] for those
#' coefficients even without a penalty.
#'
#' For a prior that is not Gaussian, \eqn{V_b} is computed in the same way for
#' every penalty that is twice differentiable at the mode. \eqn{H + S} is then
#' the curvature of the log-posterior at the mode, and \eqn{V_b} is the Laplace
#' approximation of the posterior variance. The statement about coverage above
#' holds for a Gaussian prior, whose \eqn{S} does not depend on the
#' coefficients. For another prior \eqn{S} is evaluated at the estimate and can
#' be negative. Under a Student t prior with scale \eqn{s} and \eqn{\nu}
#' degrees of freedom it is negative for every effect larger than
#' \eqn{\sqrt{\nu}\, s} in absolute value, so such an effect is shrunk little
#' and its variance comes mostly from the likelihood. The approximation is
#' local, and it can be poor where the posterior is skewed or has more than one
#' mode. Where \eqn{H + S} is not positive definite, a coordinate whose
#' diagonal entry is not positive is held as described below, and otherwise an
#' error is signaled.
#'
#' The first two matrices are conditional on the hyperparameters, which
#' are evaluated at the value where the outer search stopped, as though that
#' value were known. It was estimated from the same data, and the third matrix
#' adds the variance that this estimation contributes:
#' \deqn{V' = V_b + J V_\theta J', \qquad
#'   J = -(H + S)^{-1} \frac{\partial^2 \rho}{\partial\beta \partial\theta},}
#' with \eqn{V_\theta} the variance of the estimated hyperparameters on their
#' own free scale. It is the delta method applied to the map from the
#' hyperparameter to the penalized mode (Wood, Pya and Safken, 2016), and the
#' matrix that mgcv returns as `unconditional = TRUE`. The correction is
#' positive semi-definite, so \eqn{V'} is never narrower than \eqn{V_b}. It
#' matters where a smoothing parameter is poorly determined, and it affects the
#' penalized coordinates much more than the unpenalized ones. It requires the
#' curvature of the outer criterion, which the two conditional matrices do not
#' read, so [unconditional_vcov()] is more expensive than the other two
#' conventions.
#'
#' Where no hyperparameter was estimated by a differentiable criterion there
#' is nothing to propagate, the correction is exactly zero and
#' [unconditional_vcov()] returns \eqn{V_b} itself. Where a hyperparameter was estimated and the
#' curvature of the criterion cannot be read at the fitted point, the
#' conditional matrix is returned with a warning, so that the narrower matrix
#' is not mistaken for the wider one. Where the curvature can be read for only
#' some of the hyperparameters, the correction is a lower bound and a warning
#' states it. See [hyper_correction()].
#'
#' A coefficient that a kinked penalty has set to zero has `NA` in its row and
#' column. At zero the penalty is not twice differentiable, so \eqn{S} does
#' not exist there and the curvature cannot be read. The coefficients that a
#' lasso or an MCP leaves non-zero have a variance, which is conditional on that
#' selection; [summary.StatmodFit()] reports this in a note. Away from zero the
#' lasso contributes nothing to \eqn{S}, and SCAD and MCP contribute a negative
#' curvature on their concave part, which makes the variance of a coefficient
#' there larger than that of the likelihood alone. For these penalties \eqn{V_b}
#' is therefore not the posterior variance under the corresponding prior,
#' because a curvature read at the mode does not include the mass that the prior
#' puts near zero.
#'
#' The coefficients that the marginal criterion of [reml()] or [ml()] estimates
#' take the curvature of that criterion in the Bayesian and the unconditional
#' matrices, and the criterion reads the information of its own `information`
#' argument. The `information` of `vcov_type` selects the matrix used as
#' \eqn{H} everywhere else, in the sandwich of the frequentist matrix and in
#' the penalized block of the other two.
#'
#' The `information` of `vcov_type` selects the matrix used as \eqn{H}. Its
#' default is `NULL`, which is resolved for each fit by the family of the fit
#' with the rule of [distributions7::vcov.distrib_fit()]: the matrix is
#' [distributions7::observed_info()], the curvature of the likelihood at the
#' data in hand, when the log-likelihood is twice differentiable in every
#' parameter of the family (every [distributions7::params_order()] is at least
#' 2), and [distributions7::expected_info()] with `method = NULL`, for the
#' whole matrix, otherwise. The choice does not depend on the information with
#' which the fit was computed (`inner_information` of [statmod()]). An object
#' in `information` overrides the rule. [distributions7::expected_info()]
#' selects the Fisher information, the information averaged over the model. The
#' two agree asymptotically and not in a sample. A model that carries a
#' structural term is always inverted on the observed information, because a
#' filter and a mixture over states do not have an expected information.
#' [summary.StatmodFit()] states the matrix when it is not the observed one.
#'
#' The observed information is singular where the second derivative of the
#' log-likelihood in a coefficient is zero at the estimates, as for the
#' location of a Laplace response. The family has a kink in that parameter, so
#' the default `NULL` reads the expected information, which is finite there,
#' and the variance is available. With `information = observed_info()` the rows
#' and columns of the coordinates that it does not identify are missing, as
#' described above for a flat direction.
#'
#' The `method` of [distributions7::expected_info()] selects how the expected
#' information is obtained. With `method = NULL` the closed form of the family
#' is used where the family has one, and [distributions7::expect_opg()], the
#' outer product of the observed scores, which costs one gradient, otherwise.
#' [distributions7::expect_integrate()] evaluates the expectation itself, a sum
#' over the support for a discrete family and a quadrature for a continuous
#' one, and is orders of magnitude more expensive;
#' [distributions7::expect_mc()] simulates it.
#' [distributions7::expected_derivative_methods] compares the three methods.
#'
#' @param object A [StatmodFit()].
#' @param vcov_type A [vcov_method] object, which selects the convention of the
#'   matrix and, inside the object, the information \eqn{H}. [bayesian_vcov()]
#'   (the default) returns \eqn{V_b = (H+S)^{-1}}, the posterior variance
#'   conditional on the hyperparameters. [frequentist_vcov()] returns the
#'   sandwich \eqn{V_f = (H+S)^{-1}H(H+S)^{-1}}, conditional on the
#'   hyperparameters as well. [unconditional_vcov()] returns
#'   \eqn{V' = V_b + J V_\theta J'}, which adds the uncertainty of the
#'   hyperparameters. Each takes `information`, which is `NULL` by default (the
#'   observed information, or the expected information for a family with a kink
#'   in a parameter, resolved for each fit) and may be [distributions7::observed_info()] or
#'   [distributions7::expected_info()]. The conventions are compared at
#'   [vcov_method].
#' @param readable Whether the matrix is over the quantities that the model
#'   reports (`TRUE`, the default) or over the coordinates on which it was
#'   estimated (`FALSE`). A structural term declares its quantities, for
#'   example a persistence, and their variance is obtained by the delta method.
#' @param param The distribution parameters whose rows and columns are
#'   kept, a character vector of parameter names. `NULL` keeps all of them.
#' @param ... Unused. Any argument signals an error.
#' @return A square matrix with dimnames `parameter:name`. By default it is over
#'   the quantities that the model reports, which include the own parameters of
#'   a structural term. With `readable = FALSE` it is over the coefficients of
#'   every equation, the free parameters of a structural term following those
#'   of its own equation. `param` restricts it to the rows and columns of
#'   the named equations. A coordinate that is held or not identified has `NA`
#'   in its row and column.
#' @references
#' Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and model
#' selection for general smooth models. *Journal of the American
#' Statistical Association*, 111(516), 1548--1563.
#' @seealso [vcov_method] for the three conventions side by side,
#'   [confint.StatmodFit()], [summary.StatmodFit()],
#'   [hyper_correction()], which builds the third matrix's correction
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' sqrt(diag(vcov(fit)))
#'
#' # With a penalized term the three differ for a penalized coefficient, and
#' # the widest is the one that does not read the smoothing parameter as known.
#' ds <- data.frame(x = runif(200))
#' ds$y <- sin(2 * pi * ds$x) + rnorm(200, sd = 0.3)
#' fs <- statmod(y ~ s(x, bspline_smooth(k = 10)), distributions7::gaussian1_distrib(), ds)
#' vapply(list(frequentist_vcov(), bayesian_vcov(), unconditional_vcov()),
#'        function(ty) sqrt(diag(vcov(fs, vcov_type = ty)))[["mu:s(x).z2"]], 0)
#'
#' # the information is chosen inside the convention
#' sqrt(diag(vcov(fit, vcov_type = bayesian_vcov(
#'   information = distributions7::expected_info()))))
#' @keywords internal
vcov.StatmodFit <- function(object,
                            vcov_type = bayesian_vcov(),
                            readable = TRUE,
                            param = NULL, ...) {
  check_dots("vcov", ...)
  vcov_type <- as_vcov_method(vcov_type)
  convention <- vcov_type@name
  info <- vcov_information(vcov_type, object@spec@distrib)
  expected <- info_expected(info)
  emethod <- info_emethod(info)
  spec <- object@spec
  design <- statmod_design(spec)
  if (convention == "frequentist" && !is.null(vcov_type@sandwich)) {
    sandwich_check(object, design)
  }
  coef <- object@coefficients
  lab <- coef_labels(spec, design)
  total <- nrow(lab)
  nm <- rownames(lab)

  # A model carrying a structural term is inverted over the coefficients AND
  # the term's own parameters together, and the coefficient block of that
  # inverse is taken. Inverting the coefficient block alone would report the
  # variance that would hold with the term's parameters known, which is not
  # what was estimated. The information there is the observed one: neither a
  # filter nor a mixture over states has an expected information to offer.
  fil <- length(attr(design, "structural")) > 0L
  H <- if (fil) statmod_full_information(spec, coef, design) else
    statmod_information_at(spec, coef, design, expected, emethod)
  S <- statmod_penalty_at(spec, coef, object@hyper, design, "hessian")
  nz <- nrow(H) - total
  if (nz > 0L) {
    # A structural term's own information is assembled dense, so the padded
    # penalty is dense too rather than a sum of two kinds whose result would
    # depend on which side of the bind it arrived on.
    S <- as_dense(S)
    S <- rbind(cbind(S, matrix(0, total, nz)),
               matrix(0, nz, total + nz))
    # A penalty over the term's OWN parameters belongs in that tail. Leaving
    # it out is not a conservative choice: the deviations of a panel are not
    # identified without it -- a constant added to a population value and
    # taken off every deviation leaves the filter unchanged -- so the matrix
    # is singular along exactly that direction and nothing can be reported.
    ps <- structural_penalty_block(spec, design, object@hyper, nz)
    if (!is.null(ps)) S[total + seq_len(nz), total + seq_len(nz)] <- ps
    # and a class split between the coefficients and those parameters is
    # placed whole by the one function the fit and the criterion also read:
    # its cross block is what carries the correlation between the two, and
    # leaving it out would report each half's variance as though the other
    # were not there
    S <- S + joint_penalty_at(spec, design, coef, object@hyper, "hessian",
                              nrow(S))
  }

  keep <- rep(TRUE, total)
  beta <- unlist(coef[spec@distrib@params], use.names = FALSE)
  keep[lab$kinked & beta == 0] <- FALSE
  frz <- frozen_block(spec, lab)
  keep[frz] <- FALSE
  # the slot of a held break-point is a column of zeros by design and is
  # held rather than reported as a direction the information misses
  keep[term_held_stack(spec, design)] <- FALSE
  # An ALIASED coordinate carries no information of its own and, unlike the
  # flat directions handled below, which one it is has already been settled
  # by the pivot that fitted the model. Leaving it in makes the matrix
  # singular along a direction that is a COMBINATION of two columns, which
  # nothing can hold: the whole matrix was then refused, losing every
  # standard error for one that does not exist. Dropping it here reports the
  # rest and leaves its row and column NA, which is what `vcov()` on an
  # aliased `lm()` returns.
  keep[nm %in% object@aliased] <- FALSE
  # a kinked penalty contributes no curvature away from its kink either, and
  # any non-finite entry would be the kink itself reached by a hair
  S <- zap_nonfinite(S)

  out <- matrix(NA_real_, total, total, dimnames = list(nm, nm))
  if (!any(keep)) return(out)
  keep_full <- c(keep, rep(TRUE, nz))
  A <- (H + S)[keep_full, keep_full, drop = FALSE]
  lb <- c(nm[keep], rep("", nz))
  # what a SMALL eigenvalue means is decided inside solve_pd(), on the
  # equilibrated matrix: a smoothing parameter at 1e15 and a break-point
  # term's annealed working columns both separate the scales without
  # flattening any direction, and per-direction scaling forgives both
  Vb <- tryCatch(solve_pd(A, "the penalized information", lb),
                 error = function(e) e)
  # A COORDINATE THE INFORMATION CARRIES NOTHING ABOUT IS HELD, and the rest
  # is reported. A combination c'beta has a finite variance exactly where c
  # is orthogonal to the null space, so a coefficient the null space does
  # not touch keeps its variance whatever happens to the others: refusing
  # the whole matrix loses those, which on a count model whose dispersion
  # ran to its Poisson limit in one province was 56 standard errors thrown
  # away for one that does not exist. It runs ONLY where solve_pd() has
  # already refused, so a fit whose matrix is invertible is untouched, and
  # where nothing can be held the original message stands.
  if (inherits(Vb, "error")) {
    flat <- uninformative_coords(A)
    if (!length(flat)) stop(conditionMessage(Vb), call. = FALSE)
    warning(held_condition(sprintf(paste0(
      "The penalized information carries nothing about %s, so",
      " %s no\n  variance and %s row%s missing here. The rest of the matrix",
      " stands:\n  a coefficient is estimable exactly where it is orthogonal",
      " to the flat\n  direction, and these coordinates carry the whole of",
      " it."),
      paste(ifelse(nzchar(lb[flat]), lb[flat], "a structural parameter"),
            collapse = ", "),
      if (length(flat) == 1L) "it has" else "they have",
      if (length(flat) == 1L) "its" else "their",
      if (length(flat) == 1L) " is" else "s are")))
    keep_full[which(keep_full)[flat]] <- FALSE
    A <- (H + S)[keep_full, keep_full, drop = FALSE]
    Vb <- solve_pd(A, "the penalized information", lb[-flat])
  }
  # the coefficients the marginal criterion estimated carry its curvature in
  # place of their Schur complement; see marginal_vcov()
  Vm <- if (convention %in% c("bayesian", "unconditional")) {
    tryCatch(marginal_vcov(object, design, A, keep_full, convention),
             error = function(e) NULL)
  } else NULL
  V <- if (!is.null(Vm)) Vm else switch(convention,
    bayesian = Vb,
    # as_dense() because the information follows the DESIGN's storage: an
    # equation carrying a random effect has a sparse one, the sandwich then
    # comes back an S4 Matrix, and writing that into a slice of the base
    # matrix assembled below is a length error rather than a conversion. The
    # product is dense whatever H is, Vb being dense, so nothing is given up.
    frequentist = if (is.null(vcov_type@sandwich)) {
      as_dense(Vb %*% H[keep_full, keep_full, drop = FALSE] %*% Vb)
    } else {
      sandwich_vcov(object, design, as_dense(Vb), which(keep_full),
                    vcov_type@sandwich)
    },
    unconditional = {
      # THE OTHER TWO ARE CONDITIONAL ON THE HYPERPARAMETERS, both read at
      # the value the search stopped at as though it had been known. This
      # one adds the movement of the mode under them, and it is Vb that is
      # handed over rather than recomputed so that the two halves of the sum
      # describe one model -- the same information, the same held
      # coordinates, the same aliasing.
      # keep_full WHOLE, not its coefficient head plus a count: the tail is
      # part of the vector the mode moves in, and the flat-direction branch
      # above may have dropped a coordinate from it, which a count cannot say
      cc <- tryCatch(hyper_correction(spec, design, coef, object@hyper,
                                      object@methods$outer, Vb,
                                      keep_full, nz),
                     error = function(e) list(C = NULL, n_hyper = NA_integer_,
                                              complete = FALSE))
      if (is.null(cc$C)) {
        # zero and unavailable are different answers and the caller is told
        # which: with nothing estimated by a differentiable criterion there
        # is nothing to propagate and Vb IS the unconditional variance.
        if (!identical(cc$n_hyper, 0L)) {
          warning(conditional_condition(paste0(
            "The hyperparameters' own uncertainty could not be read here, ",
            "so the\n  conditional variance is returned. The ",
            "outer criterion's\n  Hessian could not be read at this point, ",
            "or its inverse is not a variance.")))
        }
        Vb
      } else {
        if (!isTRUE(cc$complete)) {
          warning(conditional_condition(paste0(
            "Some of the hyperparameters' uncertainty could not be read, so ",
            "the\n  correction is a lower bound: a coordinate whose own ",
            "curvature carries no\n  variance contributes nothing, and so ",
            "does a penalty over a structural\n  term's own parameters.")))
        }
        Vb + cc$C
      }
    })
  # THE JOINT INVERSE IS KEPT WHOLE. Its coefficient block is not the
  # inverse of the coefficient block wherever the two are correlated, which
  # is why the matrix is built jointly; and the structural block is the
  # variance of the term's own parameters, which used to be computed here
  # and dropped, leaving a model whose predictor is a filter with no
  # variance to report at all.
  if (any(frz)) {
    warning(frozen_condition(sprintf(paste0(
      "The block of %s is a working linearization with",
                           " a frozen weight,\n  not a Jacobian, so its",
                           " curvature is not the model's: every one of its",
                           "\n  coefficients would be reported with a",
                           " standard error of zero. They are\n  missing",
                           " here instead. The rest of the matrix stands and",
                           " is conditional\n  on the break-points where",
                           " they are; resample for their own",
      " uncertainty."),
      paste(unique(lab$term[frz]), collapse = ", "))))
  }
  out <- matrix(NA_real_, total + nz, total + nz)
  out[keep_full, keep_full] <- V
  rn <- c(nm, if (nz) rep("", nz) else character(0))
  who <- c(lab$parameter, if (nz) rep(NA_character_, nz) else character(0))
  if (nz) {
    tl <- structural_tail_names(spec, design)
    if (length(tl) == nz) {
      rn[total + seq_len(nz)] <- unname(tl)
      who[total + seq_len(nz)] <- names(tl)
    }
  }
  dimnames(out) <- list(rn, rn)
  if (isTRUE(readable)) {
    out <- readable_vcov(spec, design, object, out)
  } else if (nz) {
    # BY EQUATION, as the readable route and `coef()` both are. The joint
    # order puts every design block first and the structural tail last, so a
    # model whose predictor is a filter listed the scale's intercept before
    # the mean's own parameters and did not line up with `coef()`.
    ord <- order(match(who, spec@distrib@params))
    out <- out[ord, ord, drop = FALSE]
  }
  if (!is.null(param)) out <- select_parameter(out, param, spec)
  out
}

#' The Names of the Structural Tail of the Joint Information
#'
#' @description
#' Returns the free parameters of every structural term, under the names that a
#' coefficient of the same term would carry.
#'
#' @details
#' A level that an intercept in the same equation carries is held and is not in
#' the joint information, so it is not here either. The tail consists of the
#' free parameters in the order of the term, which is the order in which the
#' information was assembled.
#'
#' @param spec The fitted specification.
#' @param design The design.
#'
#' @return A character vector, empty where there is no structural term.
#'
#' @keywords internal
structural_tail_names <- function(spec, design) {
  su <- attr(design, "structural")
  if (!length(su)) return(character(0))
  sst <- statmod_structural_state(design)
  out <- character(0)
  for (u in su) {
    tm <- spec@terms[[u$param]][[u$term]]
    free <- setdiff(modelterms7::term_params(tm), sst$held[[u$term]])
    ds <- structural_display(tm, u$param, free)
    out <- c(out, stats::setNames(paste(ds$param, ds$name, sep = ":"),
                                  ds$param))
  }
  out
}

#' The Variance of the Quantities That a Fit Reports
#'
#' @description
#' Applies the delta method to the joint variance, with the Jacobian that
#' [readable_joint()] supplies.
#'
#' @details
#' A quantity that reads a coordinate whose variance is missing has no
#' variance either, and its row and column are left missing. Such a
#' coordinate is one that a kinked penalty has set to zero, or a parameter that
#' an intercept holds. Computing the entry from the coordinates that do have a
#' variance would report a number for something that the fit did not estimate.
#'
#' @param spec The fitted specification.
#' @param design The design.
#' @param fit The fit.
#' @param V The joint variance over the coordinates.
#'
#' @return The variance over the quantities.
#'
#' @keywords internal
readable_vcov <- function(spec, design, fit, V) {
  rj <- tryCatch(readable_joint(spec, design, fit), error = function(e) NULL)
  if (is.null(rj) || nrow(rj$jacobian) == 0L ||
      ncol(rj$jacobian) != nrow(V)) {
    return(V)
  }
  J <- rj$jacobian
  ok <- !is.na(diag(V))
  bad <- rj$held | apply(J[, !ok, drop = FALSE] != 0, 1L, any)
  Vz <- V
  Vz[is.na(Vz)] <- 0
  out <- J %*% Vz %*% t(J)
  # SYMMETRIZED, and not symmetric by construction. An entry and its
  # transpose collect the same products in a different ORDER, and
  # floating-point addition is not associative, so the two can differ
  # in the last bit -- measured, they did on Linux and not on Windows
  # or macOS, which is what a bit-for-bit assertion turns into a red
  # job on one platform of five. A variance matrix is symmetric, so
  # the halving is the contract rather than a repair.
  out <- (out + t(out)) / 2
  out[bad, ] <- NA_real_
  out[, bad] <- NA_real_
  nm <- paste(rj$parameter, rj$name, sep = ":")
  dimnames(out) <- list(nm, nm)
  out
}

#' One Equation's Submatrix of a Variance Matrix
#'
#' @description
#' Returns the rows and columns of the named distribution parameters.
#'
#' @param V A variance matrix whose names are `parameter:name`.
#' @param param The distribution parameters to keep.
#' @param spec The fitted specification, for the message.
#'
#' @return The submatrix.
#'
#' @keywords internal
select_parameter <- function(V, param, spec) {
  who <- sub(":.*$", "", rownames(V))
  keep <- who %in% param
  if (!any(keep)) {
    stop(sprintf(paste0("'param' matched nothing. The model's ",
                        "distribution parameters are %s."),
                 paste(spec@distrib@params, collapse = ", ")), call. = FALSE)
  }
  V[keep, keep, drop = FALSE]
}
S7::method(vcov, StatmodFit) <- vcov.StatmodFit

#' The Quantities That a Fit Reports, and the Map From Its Coordinates
#'
#' @description
#' Returns one row per quantity that the model reports, over the coefficients
#' of every equation and the own parameters of a structural term, together
#' with the Jacobian from the coordinates on which those were estimated.
#'
#' @details
#' This function builds the readable view that [vcov.StatmodFit()] and
#' [confint.StatmodFit()] share. [coef.StatmodFit()] builds its own readable
#' vector from the same [modelterms7::term_readable()] declarations.
#'
#' A term declares its quantities through
#' [modelterms7::term_readable()], which returns the quantities and
#' the Jacobian. The coordinates that the Jacobian touches are replaced by the
#' quantities; a coordinate that the quantities do not read stands where it
#' is, with a unit row of its own. A developed parameter is therefore left
#' intact: its development is a vector of coefficients over covariates with no
#' single value, so the term declares no quantity for it and the function
#' replaces no coordinate.
#'
#' The joint coordinate vector is the design coefficients of every equation
#' in order, then the free parameters of each structural term. A level that
#' an intercept in the same equation carries is held and is absent from the
#' information from which the variance is computed. A quantity that reads a held
#' parameter is marked: its value stands, and its variance would be only that of
#' the remaining coordinates, so it is not reported. A quantity that reads a
#' coordinate past the edge of its chart is marked in the same way
#' ([structural_edge_rows()]).
#'
#' @param spec The fitted specification.
#' @param design The design.
#' @param fit The fit.
#'
#' @return A list with `name`, `value`, `parameter`,
#'   `term`, `scale` (one link per quantity), `held` (whether
#'   the quantity reads a parameter that is not estimated), `jacobian`
#'   (quantities by joint coordinates) and `n_design`.
#'
#' @seealso [modelterms7::term_readable()],
#'   [statmod_structural_table()]
#'
#' @keywords internal
readable_joint <- function(spec, design, fit) {
  params <- spec@distrib@params
  npar <- vapply(params, function(p) design[[p]]$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  su <- attr(design, "structural")
  sst <- if (length(su)) statmod_structural_state(design) else NULL

  tails <- list()
  at <- total
  for (u in su) {
    tm <- spec@terms[[u$param]][[u$term]]
    nmp <- modelterms7::term_params(tm)
    hl <- if (is.null(sst)) character(0) else sst$held[[u$term]]
    free <- setdiff(nmp, hl)
    tails[[length(tails) + 1L]] <- list(
      u = u, tm = tm, nm = nmp, held = hl, free = free, at = at)
    at <- at + length(free)
  }
  njoint <- at
  ident <- linkfunctions7::identity_link()

  rows <- list()
  add <- function(name, value, parameter, term, scale, held, j) {
    rows[[length(rows) + 1L]] <<- list(
      name = name, value = value, parameter = parameter, term = term,
      scale = scale, held = held, j = j)
  }

  for (pi in seq_along(params)) {
    p <- params[[pi]]
    d <- design[[p]]
    if (!d$npar) next
    cf <- fit@coefficients[[p]]
    term_of <- rep(NA_character_, d$npar)
    for (nm in names(d$blocks)) term_of[d$blocks[[nm]]] <- nm
    # which coordinates a term's quantities cover, and what they are
    swap <- vector("list", d$npar)
    drop <- rep(FALSE, d$npar)
    for (nm in names(d$blocks)) {
      idx <- d$blocks[[nm]]
      if (!length(idx)) next
      term <- spec@terms[[p]][[nm]]
      rd <- tryCatch(modelterms7::term_readable(term, cf[idx]),
                     error = function(e) NULL)
      if (is.null(rd) || !length(rd$name)) next
      sup <- which(apply(rd$jacobian != 0, 2L, any))
      if (!length(sup)) next
      lb <- tryCatch(term@label, error = function(e) "")
      drop[idx[sup]] <- TRUE
      swap[[idx[sup][1L]]] <- list(rd = rd, idx = idx, term = nm, lb = lb)
    }
    for (i in seq_len(d$npar)) {
      s <- swap[[i]]
      if (!is.null(s)) {
        for (k in seq_along(s$rd$name)) {
          j <- numeric(njoint)
          j[offs[[pi]] + s$idx] <- s$rd$jacobian[k, ]
          nmk <- if (nzchar(s$lb)) paste(s$lb, s$rd$name[[k]], sep = ".") else
            s$rd$name[[k]]
          sc <- if (is.null(s$rd$scale)) ident else s$rd$scale[[k]]
          add(nmk, s$rd$value[[k]], p, s$term, sc, FALSE, j)
        }
      }
      if (drop[[i]]) next
      j <- numeric(njoint)
      j[offs[[pi]] + i] <- 1
      add(d$coef_names[[i]], cf[[i]], p, term_of[[i]], ident, FALSE, j)
    }
  }

  for (tl in tails) {
    z <- sst$zeta[[tl$u$term]]
    rd <- tryCatch(modelterms7::term_readable(tl$tm, z),
                   error = function(e) NULL)
    if (is.null(rd) || !length(rd$name)) {
      vals <- unlist(z)[tl$free]
      ds <- structural_display(tl$tm, tl$u$param, tl$free)
      for (k in seq_along(tl$free)) {
        j <- numeric(njoint)
        j[tl$at + k] <- 1
        add(ds$name[[k]], vals[[k]], ds$param[[k]], tl$u$term, ident, FALSE,
            j)
      }
      next
    }
    # the Jacobian is over the term's WHOLE parameter vector and the joint
    # matrix carries the free ones alone, so the held columns are left out
    # and a quantity that reads one is marked rather than given the variance
    # of the rest
    fi <- match(tl$free, tl$nm)
    # a quantity that reads a coordinate past the edge of its chart has no
    # variance to report either: see statmod_structural_table()
    edge_q <- tryCatch(structural_edge_rows(tl$tm, z, tl$nm),
                       error = function(e) rep(FALSE, length(rd$name)))
    ds <- structural_display(tl$tm, tl$u$param, rd$name)
    for (k in seq_along(rd$name)) {
      j <- numeric(njoint)
      j[tl$at + seq_along(fi)] <- rd$jacobian[k, fi]
      touches <- (length(tl$held) > 0L &&
        any(rd$jacobian[k, match(tl$held, tl$nm)] != 0)) || isTRUE(edge_q[[k]])
      sc <- if (is.null(rd$scale)) ident else rd$scale[[k]]
      add(ds$name[[k]], rd$value[[k]], ds$param[[k]], tl$u$term, sc, touches,
          j)
    }
  }

  # BY EQUATION, and stably. The joint coordinate order puts every design
  # block first and the structural tail last, so a model whose predictor is a
  # filter listed the scale's intercept before the mean's own parameters. The
  # Jacobian's COLUMNS keep the joint order, which is the order the variance
  # matrix is in; only the rows move.
  ord <- order(match(vapply(rows, function(z) z$parameter, character(1)),
                     params))
  rows <- rows[ord]
  list(name = vapply(rows, function(z) z$name, character(1)),
       value = vapply(rows, function(z) as.numeric(z$value), numeric(1)),
       parameter = vapply(rows, function(z) z$parameter, character(1)),
       term = vapply(rows, function(z) as.character(z$term)[[1L]],
                     character(1)),
       scale = lapply(rows, function(z) z$scale),
       held = vapply(rows, function(z) z$held, logical(1)),
       jacobian = do.call(rbind, lapply(rows, function(z) z$j)),
       n_design = total)
}


#' Invert a Matrix That Ought to Be Positive Definite
#'
#' @description
#' Returns the inverse through the Cholesky factor, and signals an error naming
#' the matrix when a direction is flat.
#'
#' @details
#' A failure here is a statement about the fit and not about the arithmetic:
#' at a maximum the penalized information is positive definite, so a matrix
#' that is not positive definite indicates that the run stopped away from a
#' maximum.
#'
#' The matrix is first scaled to a unit diagonal,
#' \eqn{A_e = D^{-1/2} A D^{-1/2}} with \eqn{D} the diagonal of \eqn{A}. The
#' test is `lmin > 1e-12 * ||A_e||_1` on the smallest eigenvalue of
#' \eqn{A_e}, and never whether `chol()` raised. On an exactly singular matrix
#' the pivot that should be zero comes out positive or negative according to
#' rounding, so the second answer differs between platforms for one matrix.
#'
#' The scaling separates a flat direction from the scale separation that a
#' large smoothing parameter legitimately produces, because a direction is
#' flat only against its own scale. Returning a pseudo-inverse instead would
#' give a standard error for a direction that the data do not identify.
#'
#' The smallest eigenvalue is estimated and not computed, from LAPACK's
#' condition estimator (`dpocon`) read on the Cholesky factor that the
#' inverse needs anyway: `rcond` is
#' \eqn{1/(\lVert A\rVert_1\lVert A^{-1}\rVert_1)}, so
#' `rcond * ||A||_1` is \eqn{1/\lVert A^{-1}\rVert_1}. For a symmetric
#' matrix the exact value of this product lies between
#' \eqn{\lambda_{\min}/\sqrt{p}} and \eqn{\lambda_{\min}}. The estimator
#' of \eqn{\lVert A^{-1}\rVert_1} is a lower bound, so the value read can
#' exceed the exact product, and can exceed \eqn{\lambda_{\min}} by a small
#' factor. The two cases to be separated differ by many
#' orders of magnitude, so the error of the estimator does not move either case
#' into the range of the other. The estimate replaces a full
#' eigendecomposition, which gives the exact answer to the same question and
#' costs \eqn{O(p^3)} with a larger constant than the Cholesky factorization.
#'
#' The message names the coefficients that load on the eigenvector of the
#' smallest eigenvalue, and no cause. A flat direction can arise from a run
#' that has not reached a maximum, from two columns of the design carrying the
#' same information, or from a parameter that has drifted to where its
#' information vanishes.
#'
#' @param A A square matrix.
#' @param what What the matrix is, for the message.
#' @param labels The names of the coefficients that index `A`.
#'
#' @return The inverse.
#'
#' @keywords internal
solve_pd <- function(A, what, labels = NULL) {
  # The verdict comes from the smallest eigenvalue and not from whether
  # chol() raised. On a matrix with an exactly zero eigenvalue -- two columns
  # of the design carrying the same information is the ordinary way to get
  # one -- the pivot that should be zero comes out positive or negative
  # according to rounding, so the same fit was refused on Windows and
  # accepted on Linux and macOS. Reading the condition estimate off the
  # factor is a statement about the matrix instead: where the factorization
  # succeeds on a singular matrix by luck, the estimate is at the rounding
  # scale and the answer is the same on every platform.
  #
  # The test runs on the JACOBI-EQUILIBRATED matrix, D^-1/2 A D^-1/2 with
  # D its diagonal, which is what tells a flat direction from scale
  # separation whatever produced the separation. A reference scale passed
  # by the caller (the unpenalized information's largest diagonal) used to
  # do that, and it covered one source only: a smoothing parameter at 1e15
  # separates the scales and the reference forgave it, but the committed
  # working block of a break-point term carries auxiliary columns near
  # 1/(2 c d) -- 1e8 at the annealed floor, 1e16 on the information's
  # diagonal -- so the DESIGN itself set the reference and ordinary
  # curvature at 242 in the ordinary-scale directions read as flat
  # (measured, on the flagship jseg's own summary). Equilibration makes
  # the test per-direction: unit diagonal, so the smallest eigenvalue is
  # small only where a direction is flat AGAINST ITS OWN SCALE, and an
  # exact collinearity stays exactly singular. The inverse is recovered
  # through the same scaling, D^-1/2 (D^-1/2 A D^-1/2)^-1 D^-1/2.
  #
  # A non-positive diagonal entry means that coordinate carries no
  # information at all, which is the flat case stated directly.
  # A non-finite entry is what a parameter run out of its range leaves behind,
  # and it is not a matrix to decompose: the factorization raises its own
  # error there, which would replace the message below with one naming
  # neither the fit nor the direction. flat_directions() reports those rows
  # instead.
  # The test and the message's flat direction both need a dense matrix, and
  # there is no cheap sparse counterpart to replace them with. Densifying
  # HERE is deliberate rather than a lapse: this runs once, at vcov(), on a
  # p by p matrix, where the sparsity that matters is in the n by p design
  # and the per-iteration products taken against it. The same judgement is
  # recorded for the observed Hessian of a regime mixture, which is also
  # computed once and left in R.
  A <- as_dense(A)
  if (ncol(A) > 0L && all(is.finite(A)) && all(diag(A) > 0)) {
    s <- 1 / sqrt(diag(A))
    Ae <- A * tcrossprod(s)
    ch <- tryCatch(chol(Ae), error = function(e) NULL)
    if (!is.null(ch)) {
      # the 1-norm, which is what the condition estimator is expressed in
      anorm <- max(colSums(abs(Ae)))
      rc <- chol_rcond_cpp(ch, anorm)
      lmin <- if (is.na(rc)) 0 else rc * anorm
      if (is.finite(lmin) && lmin > 1e-12 * anorm) {
        # Inverted through the factor already in hand. A condition number of
        # 1e15 born of scale separation costs the well-determined directions
        # nothing here and the shrunk ones simply report variances near zero.
        out <- chol2inv(ch) * tcrossprod(s)
        return((out + t(out)) / 2)
      }
    }
  }
  stop(sprintf(paste0("%s is not positive definite, so there is no variance",
                      "\n  matrix at this point.%s\n  A fit can reach such a",
                      " point without failing: a parameter that runs\n  to",
                      " the edge of its space leaves no information behind,",
                      " and the score\n  is then small because the surface is",
                      " flat, not because it is a maximum.\n  Compare the",
                      " fitted parameters against the data before reading",
                      " anything\n  else."),
               what, flat_directions(A, labels)), call. = FALSE)
}


#' The Condition That a Frozen Block Raises
#'
#' @description
#' Returns a warning of its own class, so that a caller reporting the same
#' thing in its own channel can muffle it in place of repeating it.
#'
#' @details
#' [summary.StatmodFit()] calls [vcov()] more than
#' once and would raise the warning once per call; the summary records it
#' once, as a note.
#'
#' @param msg The message.
#'
#' @return A condition.
#'
#' @keywords internal
frozen_condition <- function(msg) {
  structure(class = c("statmod_frozen_block", "warning", "condition"),
            list(message = msg, call = NULL))
}


#' The Condition Through Which a Held Coordinate Warns
#'
#' @description
#' Returns a classed warning, so that a caller can silence or catch the one
#' raised when
#' [vcov.StatmodFit()] holds a coordinate about which the information carries
#' nothing, without silencing every other warning of a fit.
#'
#' @param msg The message.
#'
#' @return A condition of class `statmod_held_coord`.
#'
#' @seealso [uninformative_coords()], which finds them, and
#'   [vcov.StatmodFit()], which raises this.
#'
#' @keywords internal
held_condition <- function(msg) {
  structure(class = c("statmod_held_coord", "warning", "condition"),
            list(message = msg, call = NULL))
}


#' The Condition Through Which an Unavailable Correction Warns
#'
#' @description
#' Returns a classed warning, so a caller can catch the one raised when
#' `vcov(vcov_type = unconditional_vcov())` cannot read the hyperparameters' own
#' uncertainty and returns the conditional variance instead.
#'
#' @details
#' It carries a class of its own for the same reason as the condition of a
#' frozen block ([frozen_condition()]): [summary.StatmodFit()] calls
#' [vcov.StatmodFit()] more than once, so without a class the same message
#' would reach the reader several times, and muffling it by position would
#' suppress whatever else was raised on the way.
#'
#' @param msg The message.
#'
#' @return A condition of class `statmod_conditional_variance`.
#'
#' @seealso [hyper_correction()], which determines whether it fires, and
#'   [vcov.StatmodFit()], which raises it.
#'
#' @keywords internal
conditional_condition <- function(msg) {
  structure(class = c("statmod_conditional_variance", "warning", "condition"),
            list(message = msg, call = NULL))
}


#' Which Coefficients Belong to a Block That Is Not a Jacobian
#'
#' @description
#' Returns the positions of every coefficient of a term whose design block is
#' a working linearization with a frozen weight.
#'
#' @details
#' A block whose term yields `FALSE` from
#' [modelterms7::term_jacobian_block()] is a working linearization with a
#' weight held at the previous iterate, so the curvature that the block carries
#' is that of the working model and not that of the model. The standard errors
#' computed from it would be zero, which would be reported as a measurement, and
#' those coefficients are left missing.
#'
#' The check reads the term through [modelterms7::term_jacobian_block()] and
#' never through its class, so a construction whose block is a Jacobian keeps
#' its inference: a continuous [modelterms7::seg()], and a discontinuous one
#' smoothed by an [numericals7::abs_smoother()], both return `TRUE`.
#'
#' @param spec The fitted specification.
#' @param lab The coefficient labels.
#'
#' @return A logical vector over the coefficients.
#'
#' @keywords internal
frozen_block <- function(spec, lab) {
  out <- rep(FALSE, nrow(lab))
  for (p in names(spec@terms)) {
    for (nm in names(spec@terms[[p]])) {
      ok <- tryCatch(modelterms7::term_jacobian_block(spec@terms[[p]][[nm]]),
                     error = function(e) TRUE)
      if (isTRUE(ok)) next
      out[lab$parameter == p & !is.na(lab$term) & lab$term == nm] <- TRUE
    }
  }
  out
}


#' When a Matrix Is Factorized Sparsely
#'
#' @description
#' Returns whether a symmetric matrix is large enough and sparse enough that a
#' sparse Cholesky factorization is used in place of a dense one. The decision
#' depends only on the matrix, and not on any term, any model or any
#' preference of the caller.
#'
#' @details
#' A matrix is factorized sparsely when it is stored as an S4 sparse matrix,
#' its order is at least `min_dim` and its fraction of nonzero entries is at
#' most `max_density`. Below about a hundred coefficients the sparse route is
#' slower, because its fixed cost is the coercion and the S4 dispatch around it,
#' which does not shrink with the matrix. The size condition prevents this on
#' the fully dense penalized information of a single smooth.
#'
#' Both quantities are read from the matrix, and the first is its storage.
#' A matrix held as a base matrix is never factorized sparsely, whatever its
#' zeros. This is not a test of the term that produced the matrix.
#' [statmod_information_at()] accumulates into the design's own kind, so the
#' penalized matrix is stored sparsely exactly when the design is. \pkg{modelterms7}
#' settles the storage of a block at build when `sparse = NULL`, the default,
#' and builds it sparse where the product of the number of rows and the number
#' of indicator columns exceeds \eqn{10^5}; `sparse = TRUE` and `sparse = FALSE`
#' force the sparse and the dense storage. Where the storage is dense the
#' factorization is a small share of the time of a fit, which is spent in the
#' \eqn{O(np^2)} products against the design, and coercing a dense
#' \eqn{p \times p} matrix here to save a share of that size would cost more
#' than it returns.
#'
#' @param M A symmetric matrix.
#' @param min_dim The smallest order for which a sparse factorization is used.
#' @param max_density The largest fraction of nonzero entries for which a
#'   sparse factorization is used.
#'
#' @return A single logical.
#'
#' @seealso [pd_factor()]
#'
#' @keywords internal
worth_sparse <- function(M, min_dim = 100L, max_density = 0.10) {
  if (!isS4(M)) return(FALSE)
  p <- ncol(M)
  if (!length(p) || p < min_dim) return(FALSE)
  nz <- tryCatch(length(M@x), error = function(e) NA_integer_)
  if (is.na(nz)) return(FALSE)
  nz / (as.numeric(p) * p) <= max_density
}


#' The Smallest Eigenvalue of a Sparse Factor's Matrix, Estimated
#'
#' @description
#' Returns \eqn{1/\lVert A^{-1}\rVert_1} from a sparse Cholesky factor, which is
#' the same quantity that LAPACK's `dpocon` produces from a dense one.
#'
#' @details
#' The sparse route needs a condition estimate of its own, because
#' `chol_rcond_cpp` reads a dense triangular factor. `Matrix::rcond` is not
#' used, since it costs more than the dense factorization that the sparse route
#' is meant to replace. The estimate is Higham's one-norm estimator
#' ([Matrix::onenormest()]) applied to the solves with the factor, which costs
#' a handful of triangular solves.
#'
#' For a symmetric matrix \eqn{\lVert A^{-1}\rVert_1 \ge \lVert
#' A^{-1}\rVert_2 = 1/\lambda_{\min}}, so the quantity returned is at or below
#' the smallest eigenvalue, and the error of the estimator is a further
#' underestimate of the norm in the other direction. It is used as the dense
#' estimate is: to separate a matrix comfortably positive definite from one
#' that is not, two situations whose smallest eigenvalues differ by many orders
#' of magnitude. An error of a factor of two in the estimate does not change
#' that verdict.
#'
#' @param L A `CHMfactor`.
#' @param p The order of the matrix.
#'
#' @return A single number, or `NA_real_` where the estimate failed.
#'
#' @seealso [pd_factor()], [pd_logdet()]
#'
#' @keywords internal
sparse_lmin <- function(L, p) {
  ax <- function(x) as.matrix(Matrix::solve(L, as.matrix(x)))
  est <- tryCatch(
    Matrix::onenormest(A.x = ax, At.x = ax, n = p, t = 1, silent = TRUE)$est,
    error = function(e) NA_real_)
  if (!is.finite(est) || est <= 0) return(NA_real_)
  1 / est
}


#' Factorize a Penalized Information Once
#'
#' @description
#' Returns the Cholesky factor of a matrix that a Laplace approximation
#' requires to be positive definite, together with its log-determinant, in the
#' storage of the matrix itself.
#'
#' @details
#' This function is the single place where the penalized matrix is factorized.
#' The criterion needs its log-determinant, the gradient needs the movement of
#' the mode and the Hessian needs both and the inverse, so one factorization
#' serves all three.
#'
#' Whether the matrix is accepted never depends on whether a factorization
#' raised an error: where the cheap test is inconclusive, the
#' eigendecomposition of the matrix settles the question. The sparse route
#' carries its own condition estimate ([sparse_lmin()]) in place of the dense
#' one, and falls back to the dense route where that estimate cannot be formed,
#' so the matrix is judged by the same criterion on either storage.
#'
#' @param M A symmetric matrix, sparse or dense.
#' @param scale A reference magnitude, as [pd_logdet()] takes.
#'
#' @return A list with `logdet`, `ok`, `factor` and
#'   `sparse`. The factor is `NULL` where the answer came from the
#'   eigendecomposition.
#'
#' @seealso [pd_logdet()], [ctx_penalized()]
#'
#' @keywords internal
pd_factor <- function(M, scale = NULL) {
  no <- function() list(logdet = NA_real_, ok = FALSE, factor = NULL,
                        sparse = FALSE)
  if (!length(ncol(M)) || !ncol(M)) return(no())

  if (worth_sparse(M)) {
    p <- ncol(M)
    if (!all(is.finite(M@x))) return(no())
    anorm <- max(Matrix::colSums(abs(M)))
    ref <- if (is.null(scale) || !is.finite(scale) || scale <= 0) anorm
           else min(scale, anorm)
    Ms <- tryCatch(Matrix::forceSymmetric(M), error = function(e) NULL)
    # CHOLMOD reports a matrix it cannot factorize as a WARNING and not as an
    # error, where the dense route raises and is caught. A hyperparameter the
    # search should step away from is not something to warn a caller about --
    # the criterion returns NULL there and the search backtracks -- so the
    # warning is suppressed and the verdict is read off the return value, as
    # it is on the dense route.
    L <- if (is.null(Ms)) NULL else
      suppressWarnings(tryCatch(Matrix::Cholesky(Ms, LDL = FALSE, super = NA),
                                error = function(e) NULL))
    if (!is.null(L)) {
      lmin <- sparse_lmin(L, p)
      if (is.finite(lmin) && lmin > sqrt(.Machine$double.eps) * ref) {
        ld <- tryCatch(2 * as.numeric(Matrix::determinant(
          L, logarithm = TRUE, sqrt = TRUE)$modulus),
          error = function(e) NA_real_)
        if (is.finite(ld)) {
          return(list(logdet = ld, ok = TRUE, factor = L, sparse = TRUE))
        }
      }
    }
    # inconclusive: the dense route asks the matrix, which is what settles it
  }

  d <- pd_logdet_dense(as_dense(M), scale)
  c(d, list(sparse = FALSE))
}


#' The Log-Determinant of a Penalized Information, Robustly
#'
#' @description
#' Returns \eqn{\log|M|} for a matrix that a Laplace approximation requires to
#' be positive definite, by the cheap route where that is safe and by a more
#' expensive one where it is not.
#'
#' @details
#' At a condition number near the rounding floor, the arithmetic and not the
#' matrix determines whether the factorization raises an error. Reading the
#' determinant off `chol(M)` alone would
#' therefore report the criterion as nonexistent at points where it exists,
#' and the outer search would backtrack from them. The same consideration
#' applies to [solve_pd()] and to the rank tests of basis7.
#'
#' The factorization is tried first, because it costs O(p^3/3) and is the
#' common case. Where it succeeds, LAPACK's condition estimator reads
#' the smallest eigenvalue off the factor already in hand for O(p^2), and a
#' matrix comfortably away from the floor is accepted with the determinant that
#' the factor gives. The eigendecomposition is computed only where that test is
#' inconclusive, or where the factorization raised an error. It costs more and
#' settles the question from the matrix itself: a factorization that failed by
#' rounding on a matrix that is in fact positive definite is recovered there,
#' and a matrix that is genuinely rank deficient is rejected deterministically
#' and not according to the platform.
#'
#' @param M A symmetric matrix.
#' @param scale A reference magnitude: the unpenalized information's own scale,
#'   so that a hyperparameter legitimately sent to 1e15 can be distinguished
#'   from a flat direction. `NULL` uses the 1-norm of `M`, and a larger value is
#'   capped at that norm.
#'
#' @return A list with `logdet` and `ok`, or `ok = FALSE` where
#'   the matrix is not positive definite.
#'
#' @seealso [solve_pd()], [statmod_marginal()]
#'
#' @keywords internal
pd_logdet <- function(M, scale = NULL) {
  r <- pd_factor(M, scale)
  r$factor <- NULL
  r$sparse <- NULL
  r
}


#' The Dense Route of pd_logdet
#'
#' @description
#' Applies the three routes described at [pd_logdet()] to a dense matrix.
#'
#' @details
#' This function is split out so that [pd_factor()] can reach it as the
#' fallback of the sparse route without restating the verdict: it is the one
#' place that determines whether a matrix is positive definite, and it holds
#' the one set of thresholds.
#'
#' @param M A dense symmetric matrix.
#' @param scale A reference magnitude.
#'
#' @return A list with `logdet` and `ok`; `factor`, the Cholesky factor, when
#'   the Cholesky route accepted the matrix; and `min_ev` and `max_ev` when the
#'   eigendecomposition rejected it.
#'
#' @seealso [pd_logdet()]
#'
#' @keywords internal
pd_logdet_dense <- function(M, scale = NULL) {
  if (!ncol(M) || !all(is.finite(M))) return(list(logdet = NA_real_, ok = FALSE))
  anorm <- max(colSums(abs(M)))
  ref <- if (is.null(scale) || !is.finite(scale) || scale <= 0) anorm
         else min(scale, anorm)
  # ⚠️ The two thresholds are NOT solve_pd()'s, and the difference is the
  # quantity. Inverting a matrix at a condition number of 1e14 loses most of
  # the answer, so solve_pd() refuses there; a log-DETERMINANT is a sum of
  # logarithms and survives it -- measured, chol and eigen agree to six
  # significant figures at 1e14 (-644.725631 against -644.725128). What
  # matters here is only that every eigenvalue is genuinely POSITIVE, and an
  # eigenvalue below p eps times the largest is not distinguishable from zero
  # in double precision whatever the factorization reports.
  ch <- tryCatch(chol(M), error = function(e) NULL)
  if (!is.null(ch)) {
    rc <- tryCatch(chol_rcond_cpp(ch, anorm), error = function(e) NA_real_)
    lmin <- if (is.na(rc)) 0 else rc * anorm
    # comfortably clear of the floor: chol's answer is its own evidence, and
    # this is the common case, so the eigendecomposition is never computed
    if (is.finite(lmin) && lmin > sqrt(.Machine$double.eps) * ref) {
      # the factor travels with the answer so that a caller wanting the
      # inverse as well does not compute a second one of the same matrix
      return(list(logdet = 2 * sum(log(diag(ch))), ok = TRUE, factor = ch))
    }
  }
  # inconclusive or refused: ask the MATRIX rather than the arithmetic
  ev <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev) || !all(is.finite(ev))) {
    return(list(logdet = NA_real_, ok = FALSE))
  }
  # THE FLOOR SCALES WITH THE DIMENSION, as the rank tests in augmented.R and
  # outer_criteria.R already read it. An exactly singular K + S comes back
  # from eigen() with a smallest eigenvalue of either sign whose size is set
  # by rounding: measured on 216 penalized informations with a duplicated,
  # rescaled or summed column (gaussian and Poisson, n 500 to 20000, p 33 to
  # 305), |min|/max reached 9 eps and never 0.27 p eps, and a floor of eps
  # alone accepted 52 of them. Over the 579 matrices the suite accepts on
  # this route the smallest legitimate ratio is 1.9e4 eps.
  if (min(ev) > ncol(M) * .Machine$double.eps * max(ev)) {
    return(list(logdet = sum(log(ev)), ok = TRUE))
  }
  list(logdet = NA_real_, ok = FALSE, min_ev = min(ev), max_ev = max(ev))
}




#' Which Coefficients a Singular Curvature Is Flat In
#'
#' @description
#' Returns a string that names the coefficients that load on the eigenvector
#' of the smallest eigenvalue.
#'
#' @param A A square matrix.
#' @param labels Its coefficient names.
#'
#' @return A single string, empty when no direction can be named.
#'
#' @keywords internal
flat_directions <- function(A, labels) {
  if (is.null(labels) || !length(labels)) return("")
  # a parameter far enough out makes its own row non-finite rather than merely
  # small, and then there is no eigenvector to read: the rows themselves are
  # the answer, and they are the more direct one
  bad <- labels[apply(!is.finite(A), 1L, any)]
  if (length(bad)) {
    return(sprintf(paste0("\n  Its entries are not finite in the rows of: %s",
                          "\n  Such rows appear when a parameter has run out of",
                          " its range."),
                   paste(bad, collapse = ", ")))
  }
  e <- tryCatch(eigen(A, symmetric = TRUE), error = function(e) NULL)
  if (is.null(e)) return("")
  k <- which.min(e$values)
  v <- abs(e$vectors[, k])
  hit <- labels[v > 0.2 * max(v)]
  # A COORDINATE WITH NO NAME IS A STRUCTURAL TERM'S OWN PARAMETER. The
  # matrix spans those beside the coefficients and only the coefficients are
  # named, so where the flat direction lies in the tail the list came out
  # empty and the message read "carried by: , ". Naming what there is to name
  # and saying where the rest is beats a row of commas.
  span <- sum(v > 0.2 * max(v))
  hit <- hit[!is.na(hit) & nzchar(hit)]
  ev <- sprintf("\n  (its smallest eigenvalue is %s, against %s at the largest).",
                format(signif(e$values[k], 3)), format(signif(max(e$values), 3)))
  if (!length(hit)) {
    if (!span) return("")
    return(paste0("\n  It is flat along a direction that lies outside the ",
                  "coefficients, and that\n  involves only a structural",
                  " term's own parameters: they are inverted together with",
                  " the coefficients.", ev))
  }
  paste0("\n  It is flat along a direction carried by: ",
         paste(hit, collapse = ", "),
         if (span > length(hit))
           " (and a structural term's own parameters)" else "", ev)
}


#' Which Coordinates a Singular Information Carries Nothing About
#'
#' @description
#' Returns the positions whose row and column of a penalized information are
#' empty, so that holding them and inverting the rest reports the variance of
#' every other coordinate exactly.
#'
#' @details
#' A combination \eqn{c'\hat\beta} has a finite asymptotic variance exactly
#' where \eqn{c} is orthogonal to the null space \eqn{\mathcal{N}} of
#' \eqn{K}, and then \eqn{\mathrm{Var}(c'\hat\beta) = c'K^{+}c} with
#' \eqn{K^{+}} the Moore-Penrose inverse. Reading `diag(K^+)` is therefore
#' not enough on its own: for a \eqn{c} that is not orthogonal it returns a
#' number where the truth is infinite. The rule is in two parts, invert away
#' from \eqn{\mathcal{N}} and report nothing on it, and this function is the
#' second part.
#'
#' A candidate is a coordinate with no curvature of its own. Its diagonal is
#' not finite, as a parameter that has run out of its range leaves it;
#' or it is not positive, which is the reading of [solve_pd()] for that entry
#' and is a statement about the coordinate regardless of its scale; or, among
#' the coordinates whose diagonal is positive, it carries the whole of a null
#' direction of the Jacobi-equilibrated matrix. The equilibration
#' distinguishes a flat direction from scale separation, and requiring one
#' coordinate to carry the whole direction excludes the collinear case, where
#' the null vector is \eqn{(e_1 - e_2)/\sqrt{2}} and each of the two carries
#' half of it.
#'
#' A candidate is held only where holding it disturbs nothing, and the
#' amount that it disturbs is the Schur correction that it removes from the
#' others, \eqn{K_{Aj}K_{jj}^{-1}K_{jA}}. That test has two branches, because it
#' has two scales.
#'
#' Where \eqn{K_{jj}} is positive the coordinate has a scale of its own and
#' the correction, read against \eqn{K_{AA}} in equilibrated units where that
#' block has a unit diagonal, is
#' \deqn{\max_k \Bigl(\frac{|K_{jk}|}{\sqrt{K_{kk}}}\Bigr)^{2}
#'       \Big/ K_{jj} \;\le\; \texttt{schur}.}
#' It is dimensionless and needs no reference: a coordinate measured a
#' thousand times smaller than its neighbors has a small row and a small
#' diagonal, and the ratio is unchanged by the scaling. Scaling one coordinate
#' down by \eqn{10^{-14}} leaves the ratio unchanged, and the function holds no
#' coordinate.
#'
#' Where \eqn{K_{jj}} is not positive the coordinate has no scale of its own
#' (the information in that direction is zero and not merely small), so that
#' ratio is \eqn{0/0} and carries no information. The only scale left is that
#' of the matrix, and the question becomes whether the row sits at the level
#' to which its entries were accumulated:
#' \eqn{\max_k |K_{jk}| \le \texttt{row\_tol}\max|K|}. For a count model whose
#' dispersion has left its range, \eqn{K_{jj}} is exactly zero and the entries
#' of the row of that coordinate are at the rounding level of the matrix.
#'
#' Where both branches fail, dropping the coordinate would leave
#' \eqn{K_{AA}^{-1}}, the variance conditional on that coefficient being
#' known, which is smaller than the marginal one; the function holds no
#' coordinate and the error signaled by the caller stands.
#'
#' @param A A penalized information, over the coefficients and any structural
#'   tail.
#' @param tol The relative eigenvalue below which a direction is flat.
#' @param share How much of the null space one coordinate must carry to be
#'   held, as a share of one.
#' @param schur The largest Schur correction, in equilibrated units, that
#'   counts as no disturbance to the coordinates that are kept, where the
#'   candidate has a positive diagonal.
#' @param row_tol The share of the matrix's largest entry below which a row
#'   is read as accumulated rounding, where the candidate does not.
#'
#' @return An integer vector of positions in `A`, possibly empty.
#'
#' @seealso [solve_pd()], which signals an error for the matrix from which
#'   these are taken, and [vcov.StatmodFit()], which holds them.
#'
#' @keywords internal
uninformative_coords <- function(A, tol = 1e-10, share = 1 - 1e-6,
                                 schur = 1e-8, row_tol = 1e-12) {
  A <- as_dense(A)
  p <- ncol(A)
  if (!p) return(integer(0))
  # BY THE DIAGONAL, which is boundary_coords()' rule: a coordinate at a
  # bound makes its whole row non-finite, cross terms included, so a row
  # test marks its neighbours as well. Whatever is left non-finite among
  # the survivors goes too, there being no matrix there either.
  bad <- boundary_coords(A)
  rest <- setdiff(seq_len(p), bad)
  if (length(rest)) {
    extra <- rest[apply(!is.finite(A[rest, rest, drop = FALSE]), 1L, any)]
    bad <- c(bad, extra)
    rest <- setdiff(rest, extra)
  }
  if (!length(rest)) return(sort(bad))
  d <- diag(A)
  # A NON-POSITIVE DIAGONAL is the flat case stated directly: that
  # coordinate has no curvature of its own, whatever scale it is on, and
  # the sign of an entry at the rounding level is a coin toss. A positive
  # one equilibrates and is read as a direction of the correlation matrix.
  cand <- rest[d[rest] <= 0]
  rest <- setdiff(rest, cand)
  if (length(rest) >= 2L) {
    B <- A[rest, rest, drop = FALSE]
    s <- 1 / sqrt(diag(B))
    Be <- B * tcrossprod(s)
    e <- tryCatch(eigen(Be, symmetric = TRUE), error = function(e) NULL)
    if (!is.null(e)) {
      ev <- abs(e$values)
      flat <- ev <= tol * max(ev)
      if (any(flat)) {
        # each coordinate's share of the null space: one where the
        # direction IS that coordinate, a half where two of them share it
        load <- rowSums(e$vectors[, flat, drop = FALSE]^2)
        cand <- c(cand, rest[load > share])
        rest <- setdiff(rest, rest[load > share])
      }
    }
  }
  if (!length(cand)) return(sort(bad))
  # HOLDING IS EXACT ONLY WHERE IT DISTURBS NOTHING, and the amount it
  # disturbs is the Schur correction A_Aj A_jj^-1 A_jA it removes from the
  # others. Read in equilibrated units, where the kept block has a unit
  # diagonal, that is max_k (|A_jk|/sqrt(A_kk))^2 / |A_jj|, and it needs no
  # reference scale: a coordinate measured a thousand times smaller than
  # its neighbours has a small row AND a small diagonal, and the ratio sees
  # through both. Where the correction is not negligible, dropping the
  # coordinate would report the variance CONDITIONAL on it as if it were
  # the marginal one, so nothing is held and the caller's refusal stands.
  if (length(rest)) {
    sk <- sqrt(d[rest])
    scale <- max(abs(A[rest, rest]))
    ok <- vapply(cand, function(j) {
      row <- abs(A[j, rest])
      if (!all(is.finite(row))) return(FALSE)
      if (d[j] > 0) {
        # a scale of its own: the correction is r^2/A_jj, dimensionless
        r <- max(row / sk)
        r * r <= schur * d[j]
      } else {
        # NO scale of its own, so r^2/A_jj is 0/0 and says nothing. The
        # only scale left is the matrix's, and the question is whether the
        # row sits at the level its entries were accumulated to.
        max(row) <= row_tol * scale
      }
    }, logical(1))
    cand <- cand[ok]
  }
  sort(c(bad, cand))
}


#' The Coordinates That a Fit Does Not Identify, Found After the Fact
#'
#' @description
#' Returns the coordinates that a pivoted decomposition of the penalized
#' information at the mode leaves out, for a fit whose own solve named none.
#'
#' @details
#' Only [iwls()] on a pivoting decomposition reports which column it
#' dropped, which is a by-product of the solve that fitted the model. An
#' `optimizers7` method does not solve by a pivot, and the `chol` and `svd`
#' routes do not either, so a design of less than full
#' rank left every coordinate unnamed and [vcov.StatmodFit()] signaled an
#' error for the whole matrix. The question that a pivot settles during the
#' fit is settled here on the matrix that `vcov()` inverts, so that the alias
#' and the variance refer to the same model.
#'
#' The matrix is \eqn{K = H + S}, which makes a penalized coordinate safe
#' without a clause of its own: the curvature of the penalty is in \eqn{S}, so
#' a column that the design alone does not identify is identified in \eqn{K},
#' as it is identified in the augmented system that the pivoted route
#' factorizes. Under `ridge(~ 0 + x1 + x3)` with two identical columns, the
#' function names no coordinate.
#'
#' \eqn{K} is equilibrated to a unit diagonal before the pivot runs, as in
#' [solve_pd()] and the sparse rank test: a smoothing parameter sent to 1e15
#' and the annealed columns of a break-point term both separate the scales
#' without flattening a direction, and per-direction scaling accommodates both.
#' A coordinate whose curvature is small only because of its scale is not
#' named.
#'
#' A coordinate is a candidate only where its own curvature is finite and
#' positive, which is the rule of [boundary_coords()]: a parameter at the clamp
#' of its link makes its whole row non-finite, and a coordinate about which the
#' information carries nothing has an empty row. Both cases differ from
#' aliasing. In them the estimate stands and only the variance does not, a
#' situation that [uninformative_coords()] handles, whereas aliasing reports
#' the coefficient itself as missing.
#'
#' [solve_pd()] determines whether there is anything to name, and a tolerance
#' of this function does not, so a fit whose information inverts is untouched
#' and the alias cannot disagree with the variance reported beside it. The gate
#' is necessary because \eqn{K} is \eqn{X'X} up to the weights, so it squares
#' the conditioning of the design, and a pivoted QR read at the tolerance of
#' `dqrdc2` is twice as strict here as it is on the augmented system. For two
#' columns collinear to within 1e-4 (an ordinary pair of correlated
#' covariates), the bare pivot would name one of them while `solve_pd()`
#' inverts the matrix without difficulty. With the gate, the function leaves
#' such a case untouched and names a coordinate only where the whole matrix
#' would otherwise have been rejected.
#'
#' Which coordinate is then named is determined by the pivot, as in
#' [augmented_solve()].
#'
#' A coordinate that a kinked penalty holds at zero is left out of the test,
#' as [vcov.StatmodFit()] leaves it out of the matrix that it inverts. The
#' kink identifies such a coordinate although \eqn{K} carries no curvature for
#' it: under `1 + lasso(~ g)` with every level of `g` coded, the direction
#' that raises the intercept and lowers every level is flat in \eqn{K} and
#' not in the objective, whose minimum over it is unique. Reading the pivot on
#' the whole matrix would name the last level, so that coordinate is excluded.
#'
#' @param K The penalized information at the mode, \eqn{H + S}.
#' @param held The coordinates that a kinked penalty holds at zero, from
#'   [zero_kinked()], which are not candidates.
#'
#' @return An integer vector of coordinates of the coefficient vector,
#'   possibly empty.
#'
#' @seealso [uninformative_coords()], which holds a coordinate about which the
#'   information carries nothing instead of aliasing it, and
#'   [augmented_solve()], whose pivot settles the same question during the
#'   fit.
#'
#' @keywords internal
deficient_coords <- function(K, held = integer(0)) {
  K <- tryCatch(as_dense(K), error = function(e) NULL)
  if (is.null(K) || !is.matrix(K) || ncol(K) < 2L) return(integer(0))
  d <- diag(K)
  # by the DIAGONAL and not by the row: one non-finite column would
  # otherwise make every row non-finite and leave nothing to test, which
  # is what a first version did on a fit carrying both a clamp and a
  # duplicated column
  rest <- setdiff(which(is.finite(d) & d > 0), held)
  if (length(rest) < 2L) return(integer(0))
  A <- K[rest, rest, drop = FALSE]
  if (!all(is.finite(A))) return(integer(0))
  # nothing is named where the matrix inverts: the pivot on K is twice as
  # strict as the one on the design, K squaring the conditioning, so read
  # alone it would alias an ordinary pair of correlated covariates
  inverts <- tryCatch({ solve_pd(A, "the penalized information"); TRUE },
                      error = function(e) FALSE)
  if (inverts) return(integer(0))
  s <- 1 / sqrt(d[rest])
  q <- tryCatch(qr(A * tcrossprod(s)), error = function(e) NULL)
  if (is.null(q) || q$rank >= length(rest)) return(integer(0))
  sort(rest[setdiff(seq_along(rest), q$pivot[seq_len(q$rank)])])
}

#' Coordinates Whose Design Column Has Vanished
#'
#' @description
#' Returns the coordinates of the stacked coefficient vector whose column of
#' the design, read at the fitted coefficients, has vanished: its norm is at
#' most `eps` times the largest column norm of its own equation, and the
#' penalized information carries nothing on its diagonal either.
#'
#' @details
#' This function covers the one case that [deficient_coords()] leaves unnamed
#' by design. That function reads only \eqn{K}, and an empty row of \eqn{K} can
#' arise in two different ways, which the matrix does not distinguish. At a
#' parameter's boundary the design column is alive and the working weight
#' vanishes: the estimate stands and only its variance does not, a case that
#' [uninformative_coords()] handles.
#' Where a coefficient's own column vanishes, the data do not move with it
#' and it is not identified at all: a break-point run out of the data has
#' change and position columns that are identically zero there. A sharp step
#' or a truncated line written into [modelterms7::nl()] can park a break-point
#' beyond the data in this way, and without this rule the function would name
#' one of the two coordinates, or would name none.
#'
#' Both conditions are required. If only the design column were used, the
#' function would name a penalized coordinate whose column is empty (a level of
#' a random effect with no observations, which its prior identifies), and if
#' only the diagonal were used, it would name a coordinate at a boundary. The
#' column is compared with its own
#' equation's largest, since two equations' designs need not share a scale,
#' and the diagonal with \eqn{\epsilon^2} times the largest, the square of the
#' same resolution.
#'
#' @param spec The fitted specification.
#' @param coef The coefficients, one vector per distribution parameter.
#' @param design The design that [statmod_design()] built.
#' @param K The penalized information at `coef`, as [deficient_coords()]
#'   reads it.
#'
#' @return An integer vector of coordinates of the stacked coefficient
#'   vector, possibly empty.
#'
#' @seealso [deficient_coords()] for a column that the others span,
#'   [uninformative_coords()] for a coordinate at a boundary.
#'
#' @keywords internal
vanished_coords <- function(spec, coef, design, K) {
  K <- tryCatch(as_dense(K), error = function(e) NULL)
  if (is.null(K) || !is.matrix(K)) return(integer(0))
  d <- diag(K)
  fin <- is.finite(d)
  top <- max(c(0, d[fin & d > 0]))
  flat_k <- fin & d >= 0 & d <= .Machine$double.eps^2 * top
  if (!any(flat_k)) return(integer(0))
  da <- statmod_design_at(spec, coef, design)
  out <- integer(0)
  off <- 0L
  for (p in spec@distrib@params) {
    m <- design[[p]]$npar
    if (m == 0L) next
    X <- da[[p]]$X
    nr <- sqrt(as.numeric(if (inherits(X, "Matrix")) Matrix::colSums(X^2)
                          else colSums(X^2)))
    top_x <- max(c(0, nr[is.finite(nr)]))
    gone <- which(is.finite(nr) & nr <= .Machine$double.eps * top_x)
    out <- c(out, off + gone)
    off <- off + m
  }
  out <- out[out <= length(d)]
  sort(out[flat_k[out]])
}


#' Report a Model Whose Columns Are Not All Identified
#'
#' @description
#' Raises a warning naming the coordinates that the fit could not estimate, so
#' that a reader who does not print the summary learns that the model as
#' written does not identify every column of its design.
#'
#' @details
#' A coefficient that the model does not identify is one point of a flat ridge:
#' the fit converges and reports numbers, and the default print does not
#' indicate which of them are meaningful. `lm()` and `glm()` put `NA` in
#' `coef()`, and [statmod()] follows that convention. The warning also reports
#' the coordinates when `summary()` is not printed, where an `NA` among the
#' rows of a table is easy to overlook.
#'
#' The coordinates are those of [deficient_coords()], read at the fitted
#' coefficients on \eqn{K = H + S} and gated by [solve_pd()], which makes the
#' warning safe to raise. A rank test on the raw design addresses a different
#' question and has false positives that would make the warning noisy: it
#' reports a deficiency on `random(~ 1 | g)`, on `random(~ x | g)` and on
#' `s(x) + random(~ 1 | g)`, whose columns are identified by their own penalty
#' and whose fits are ordinary. Read on \eqn{K}, the curvature of the penalty
#' is in the matrix and the function names none of the three. The warning does
#' fire on `nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g))`, where the intercept of a
#' parameter sits beside a complete set of indicators that sum to it.
#'
#' Reading the matrix at the fitted coefficients and not at the start is part
#' of the rule: a Jacobian block is degenerate at a starting value of zero for
#' a reason that goes away as soon as the parameter moves, and
#' `nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g))`, which is identified, is
#' deficient there and is not named here.
#'
#' The message names no cause, because two different shapes reach it. In the
#' first the design is deficient because a column is the sum of others. In the
#' second the design has full rank and \eqn{K} does not, as for
#' `jump(x, smoothed = smooth_quintic(h = g))` at one median spacing \eqn{g},
#' and for `jseg()` under the same smoother: the block has rank 3 of 3 while
#' `mu:jump.psi1` is named, because the quintic is exact outside its own width,
#' so that where every observation lies farther than \eqn{h} from the
#' break-point its derivative column carries no curvature. The quintic
#' resolves its width at \eqn{5/(2\log 2)} median spacings, which avoids this
#' on a uniform covariate;
#' a width supplied by hand can still produce it. Both cases are genuine
#' (in each of them the coefficient has no standard error), and a message
#' asserting the first would be wrong about the second, so the message reports
#' the fact and leaves the diagnosis to the reader. The two smoothers whose
#' derivative has unbounded support, `smooth_probit()` and
#' `smooth_hyperbolic()`, are not named on the same data.
#'
#' @param alias The aliased coordinates by name, as
#'   [aliased_labels()] gives them.
#'
#' @return `NULL`, invisibly. The function is called for the warning.
#'
#' @seealso [deficient_coords()], which determines which coordinates those are,
#'   and [vcov.StatmodFit()], which reports them as missing.
#'
#' @keywords internal
warn_aliased <- function(alias) {
  if (!length(alias)) return(invisible(NULL))
  nm <- as.character(alias)
  k <- length(nm)
  shown <- if (k > 4L) c(nm[seq_len(4L)], sprintf("and %d more", k - 4L)) else nm
  warning(sprintf(paste0(
    "%d coefficient%s not identified at the fitted point: the penalized ",
    "information\n  is flat in %s direction, so the values reported there are ",
    "one point of a ridge and\n  the estimate, the standard error and the ",
    "interval are missing.\n  %s"),
    k, if (k == 1L) " is" else "s are", if (k == 1L) "its" else "their",
    paste(shown, collapse = ", ")), call. = FALSE)
  invisible(NULL)
}


#' @title Confidence Intervals for a Fit
#' @name confint.StatmodFit
#' @description
#' Returns confidence intervals for the coefficients of every distribution
#' parameter, by default the Wald ones and otherwise by inverting a likelihood
#' test.
#' @details
#' With `readable = FALSE`, the interval for a coefficient of a linear
#' predictor is symmetric about the estimate and needs no mapping back. Such a
#' coefficient is unbounded whatever the distribution parameter to which it
#' belongs, because the link has already carried that parameter onto the whole
#' line, so the scale on which the interval is built is the scale on which the
#' quantity lives. The interval does not respect a bound on the distribution
#' parameter itself; for that, map an interval for the predictor through the
#' inverse link at the covariate values of interest. With the default `readable
#' = TRUE`, each quantity that a structural term declares (a loading or a
#' persistence, for example) has an interval built on the link scale of that
#' quantity and mapped back, so that interval is not symmetric.
#'
#' The variance comes from [vcov.StatmodFit()], so the same three
#' conventions apply, and a coefficient that a kinked penalty has set to zero
#' has `NA` in place of an interval. For an interval around a penalized term,
#' `vcov_type = unconditional_vcov()` does not treat the smoothing parameter as
#' known.
#'
#' # Inverting a test instead
#'
#' The argument `test` takes a [test_statistic] object (see there for the four
#' statistics) and selects the statistic whose acceptance region is the
#' interval. The Wald interval is the set of \eqn{b} that the Wald statistic
#' does not reject, and it is symmetric by construction: that statistic is a
#' parabola in \eqn{b}, since the curvature is read once at \eqn{\hat\beta}. The
#' other three read the likelihood at each \eqn{b} in turn, so their
#' intervals follow its shape and are not symmetric where it is not. Where the
#' likelihood is nearly quadratic the four intervals agree; with a small sample
#' and a skewed likelihood the likelihood-ratio and gradient intervals shift
#' away from the Wald interval.
#'
#' Each of the three is found by [statmod_invert()], which brackets the
#' crossing outward from the Wald interval and locates it with a root finder,
#' so one interval requires a few dozen restricted refits. The `parm` argument
#' therefore selects the coefficients before the search and not after it, so
#' that only the requested intervals are computed.
#'
#' A row for which the restricted fit is not defined keeps `NA` limits.
#' [testable_coords()] identifies those rows: a coefficient under a kinked
#' penalty, an aliased one, and every coefficient of a model carrying a
#' structural term. The estimate and the standard error in such a row are those
#' of the fit and do not depend on `test`.
#'
#' A coefficient under a penalty that is twice differentiable (a smooth's own
#' coordinates, a ridge, a random effect) is inverted like any other, and the
#' interval is the credible one that [vcov.StatmodFit()] gives for it under
#' `vcov_type = bayesian_vcov()`, conditional on the hyperparameters that the fit
#' reached. See [statmod_stat_at()].
#'
#' The three restricted tests ([lr_test()], [score_test()] and
#' [gradient_test()], the ones for which [test_restricted()] is `TRUE`) need
#' `readable = FALSE`, since a test is about one coefficient and a readable
#' quantity is a function of several at once.
#' @param object A [StatmodFit()].
#' @param parm Which coefficients: a distribution parameter's name, a term's
#'   name, a vector of `parameter:coefficient` labels, or `NULL` for all of
#'   them. With `readable = TRUE`, the name of a quantity is accepted as well.
#' @param level The confidence level.
#' @param vcov_type A [vcov_method] object, passed to [vcov.StatmodFit()]:
#'   [bayesian_vcov()] (the default), [frequentist_vcov()] or
#'   [unconditional_vcov()]. The conventions are compared at [vcov_method].
#' @param readable Whether to report the quantities in which the model is written
#'   (a term's own parameters under the names of the literature) and not the
#'   raw coordinates. The value is passed to [vcov.StatmodFit()] and must be
#'   `FALSE` for a restricted test ([test_restricted()]).
#' @param test A [test_statistic] object, the test whose acceptance region is
#'   the interval. [wald_test()] (the default) gives the interval that is
#'   symmetric about the estimate, from the variance matrix and with no refit.
#'   [lr_test()], [score_test()] and [gradient_test()] are inverted by
#'   [statmod_invert()] and cost a few dozen restricted refits per interval.
#'   The four are compared at [test_statistic].
#' @param ... Passed to [vcov.StatmodFit()], which signals an error for any
#'   argument that it does not take.
#' @return A data frame with the parameter, the term, the coefficient, the
#'   estimate, its standard error and the two limits.
#' @seealso [vcov.StatmodFit()], [summary.StatmodFit()], [test_statistic],
#'   [vcov_method], [statmod_test()], which tests one coefficient against a
#'   value of its own, and [statmod_invert()] and [statmod_stat_at()], the two
#'   internal
#'   functions from which an inverted interval is built
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' confint(fit)
#' confint(fit, "sigma")
#' confint(fit, parm = "mu:x", test = lr_test(), readable = FALSE)
#' @keywords internal
confint.StatmodFit <- function(object, parm = NULL, level = 0.95,
                               vcov_type = bayesian_vcov(),
                               readable = TRUE,
                               test = wald_test(),
                               ...) {
  vcov_type <- as_vcov_method(vcov_type)
  test <- as_test_statistic(test)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  # A TEST IS ABOUT ONE COEFFICIENT and a readable quantity is a function of
  # several at once, so there is no single coordinate to hold: the request is
  # refused rather than answered on the coordinates under the other names.
  if (test_restricted(test) && isTRUE(readable)) {
    stop(sprintf(paste0("test = %s_test() inverts a test on one coefficient, ",
                        "so it needs\n  readable = FALSE. The readable ",
                        "quantities are functions of several\n  ",
                        "coefficients at once, so one coordinate cannot be ",
                        "held in their place."), test@name), call. = FALSE)
  }
  spec <- object@spec
  design <- statmod_design(spec)
  z <- stats::qnorm((1 + level) / 2)
  V <- vcov(object, readable = readable, vcov_type = vcov_type, ...)
  se <- sqrt(diag(V))
  if (isTRUE(readable)) {
    rj <- readable_joint(spec, design, object)
    est <- rj$value
    # EACH INTERVAL ON THE SCALE THAT KEEPS ITS QUANTITY IN ITS OWN SET, and
    # mapped back, as every interval in the toolkit is: a loading is
    # positive and rides a log, an autoregressive coefficient is confined to
    # no interval a scalar link expresses and rides the identity. Read on
    # the raw scale a positive quantity routinely gets a negative lower end.
    lo <- hi <- rep(NA_real_, length(est))
    for (i in seq_along(est)) {
      if (!is.finite(se[[i]])) next
      g <- rj$scale[[i]]
      t0 <- linkfunctions7::linkfun(g, est[[i]])
      st <- se[[i]] * abs(linkfunctions7::dlinkfun(g, est[[i]]))
      ends <- sort(c(linkfunctions7::linkinv(g, t0 - z * st),
                     linkfunctions7::linkinv(g, t0 + z * st)))
      lo[[i]] <- ends[[1L]]
      hi[[i]] <- ends[[2L]]
    }
    out <- data.frame(parameter = rj$parameter, term = rj$term,
                      coefficient = rj$name, estimate = est, se = se,
                      lower = lo, upper = hi, stringsAsFactors = FALSE)
    rownames(out) <- rownames(V)
    if (length(object@aliased)) {
      al <- rownames(out) %in% object@aliased
      if (any(al)) out$estimate[al] <- NA_real_
    }
    if (is.null(parm)) return(out)
    keep <- rownames(out) %in% parm | out$parameter %in% parm |
      out$term %in% parm | out$coefficient %in% parm
    if (!any(keep)) {
      stop(sprintf(paste0("'parm' matched nothing. It takes a distribution",
                          " parameter\n  (%s), a term's name, a quantity's",
                          " name, or a label of the\n  form",
                          " 'parameter:name'."),
                   paste(spec@distrib@params, collapse = ", ")),
           call. = FALSE)
    }
    return(out[keep, , drop = FALSE])
  }
  lab <- coef_labels(spec, design)
  est <- unlist(object@coefficients[spec@distrib@params], use.names = FALSE)
  # BY NAME, not by position. The matrix carries a structural term's free
  # parameters right after the coefficients of the equation it sits in, so
  # the first nrow(lab) entries of its diagonal are the coefficients only
  # where that equation is the last one. Read by position, a filter in the
  # mean gave the scale's intercept the standard error of the filter's
  # log-loading (Nile, gas(1, 1): 0.3295 printed against 0.0707).
  se <- unname(se[rownames(lab)])
  out <- data.frame(lab[, c("parameter", "term", "coefficient")],
                    estimate = est, se = se,
                    lower = est - z * se, upper = est + z * se,
                    stringsAsFactors = FALSE)
  rownames(out) <- rownames(lab)
  # an aliased coordinate has no estimate to report, and the zero the solve
  # left there would read as one: its standard error and its interval are
  # already missing, the variance matrix having held it
  if (length(object@aliased)) {
    al <- rownames(out) %in% object@aliased
    if (any(al)) out$estimate[al] <- NA_real_
  }
  # THE SUBSET IS TAKEN FIRST where a test is being inverted, because each
  # row then costs a few dozen restricted refits: filtering afterwards, as
  # the Wald route can afford to, would pay for every coefficient of the
  # model to report one of them.
  if (!is.null(parm)) {
    keep <- rownames(out) %in% parm | out$parameter %in% parm |
      out$term %in% parm
    if (!any(keep)) {
      stop(sprintf(paste0("'parm' matched nothing. It takes a distribution",
                          " parameter\n  (%s), a term's name, or a label of",
                          " the form 'parameter:coefficient'."),
                   paste(spec@distrib@params, collapse = ", ")), call. = FALSE)
    }
    out <- out[keep, , drop = FALSE]
  }
  if (test_restricted(test)) {
    out <- invert_rows(object, out, level, test, vcov_type)
  }
  out
}
S7::method(confint, StatmodFit) <- confint.StatmodFit


#' Replace a Table's Limits by an Inverted Test
#'
#' @description
#' Runs [statmod_invert()] on each row for which a restricted fit is defined
#' and writes the two limits back, leaving the rest of the table as it was.
#'
#' @details
#' A row that the restricted fit cannot hold keeps `NA` limits rather than the
#' Wald limits that it carried on arrival, so that one column never carries two
#' different intervals. [testable_coords()] determines which rows those are; if
#' it names none of the rows of the table, the function signals an error, since
#' a table that consists only of missing values is less informative than the
#' reason.
#'
#' @param fit A [StatmodFit()].
#' @param out The table that [confint.StatmodFit()] has built, already subset.
#' @param level The confidence level.
#' @param test The [test_statistic] to invert.
#' @param vcov_type The [vcov_method] that sets the starting bracket.
#'
#' @return `out`, in which the `lower` and `upper` columns are replaced.
#'
#' @seealso [statmod_invert()], which does the search.
#'
#' @keywords internal
invert_rows <- function(fit, out, level, test, vcov_type) {
  ok <- names(testable_coords(fit))
  hit <- rownames(out) %in% ok
  if (!any(hit)) {
    stop(sprintf(paste0("The requested coefficients cannot be held at a",
                        " value, so test = %s_test() has\n  nothing to invert. A ",
                        "coefficient under a kinked penalty, an aliased one",
                        "\n  and every coefficient of a model carrying a ",
                        "structural term are all outside it."), test@name),
         call. = FALSE)
  }
  out$lower <- rep(NA_real_, nrow(out))
  out$upper <- rep(NA_real_, nrow(out))
  for (i in which(hit)) {
    ci <- tryCatch(statmod_invert(fit, out$parameter[[i]],
                                  out$coefficient[[i]], level, test, vcov_type),
                   error = function(e) c(lower = NA_real_, upper = NA_real_))
    out$lower[[i]] <- ci[["lower"]]
    out$upper[[i]] <- ci[["upper"]]
  }
  out
}


#' What Kind of Block a Term Reports As
#'
#' @description
#' Returns the reading that a term receives in a summary: a structural term's
#' own parameters, the quantities of a break-point term, the coefficients of an
#' unpenalized term, a smooth's linear part and smoothing parameter, a random
#' effect's variance parameters, a selection's survivors, or the coefficients
#' of another penalized term.
#'
#' @details
#' The classification is by the term's class and by its penalties, not by its
#' label, so a term given a name of its own is read the same way. A structural
#' term (a filter or a regime chain) and a break-point term are classified
#' first and report as `"structural"` and `"breakpoint"` whatever penalties
#' they carry. A random effect that carries a label is read as `"random"`,
#' because the covariance class that collects it carries its penalty. For any
#' other term the penalties are the ones that the term declares through
#' [modelterms7::term_penalties()]: a term with none, and with no labeled
#' random effect in its subformulas, is `"parametric"`; a term penalized over
#' part of its parameters is read as penalized, never as parametric, and as a
#' selection when any of its penalties has a kink.
#'
#' @param term A built term.
#'
#' @return One of `"structural"`, `"breakpoint"`,
#'   `"parametric"`, `"smooth"`, `"random"`,
#'   `"selection"`, `"penalized"`.
#'
#' @keywords internal
term_block_kind <- function(term) {
  # a structural term is asked about FIRST OF ALL: it contributes no design
  # columns, so what it is has nothing to do with whether it carries a
  # penalty. Asked later it came back "parametric" whenever it carried none,
  # which put a term with no columns among the terms whose columns are one
  # unpenalized block.
  if (S7::S7_inherits(term, modelterms7::structural_term)) {
    return("structural")
  }
  # a break-point term is asked about next, before its penalties are
  # looked at: its block is a working linearization whose coefficients are
  # not the quantities of the model, so it wants a section of its own
  # whether or not a development of its coefficients carries a penalty
  if (S7::S7_inherits(term, modelterms7::SegTerm)) return("breakpoint")
  # a LABELLED random effect is asked about before its penalties too, and
  # for the same reason: it declares none, its coefficients being covered
  # by the covariance class's, and read after the penalties it came back
  # "parametric" -- forty grouping indicators printed one per line under
  # a heading that says they are an unpenalized block
  if (S7::S7_inherits(term, modelterms7::RandomTerm) &&
      !is.na(modelterms7::term_tag(term))) {
    return("random")
  }
  ent <- modelterms7::term_penalties(term)
  # AND THE SAME ONE LEVEL DOWN. A term developing one of its own parameters
  # over a labelled random effect -- `nl(~ a * exp(-r * x), a ~ 1 +
  # random(~ 1 | u | id))` -- declares no penalty either, the class carrying
  # it, and read after the penalties it came back parametric: twelve
  # predictions printed one per line with a z and a p beside them, under a
  # heading saying they are an unpenalized block, and the term's own
  # compartments gone. The question is asked of everything under the term,
  # so a term written later is covered without an edit here.
  if (!length(ent) && !length(term_tags_deep(term))) return("parametric")
  if (S7::S7_inherits(term, modelterms7::RandomTerm)) return("random")
  if (S7::S7_inherits(term, modelterms7::SmoothTerm)) return("smooth")
  if (any(vapply(ent, function(e) penalty_has_kink(e$penalty), logical(1)))) {
    return("selection")
  }
  "penalized"
}


#' Which Coefficients of a Smooth Are the Linear Part
#'
#' @description
#' Returns `TRUE` for the columns of a smooth's block that are not touched by
#' any penalty of the term, which carry its unpenalized (linear or polynomial)
#' effect and are the columns that the summary prints.
#'
#' @details
#' The remaining coefficients of the block are coordinates of an orthonormal
#' basis of the wiggly part. A single one of them has no interpretation, and
#' together they determine the effective degrees of freedom, which the block
#' header reports instead.
#'
#' The columns are identified from the construction and not from a suffix in a
#' coefficient's name, which is only a label. The decisive quantity is the
#' penalty: a column that the roughness matrix leaves alone belongs to the
#' null space that the smoother kept, and that null space carries the
#' unpenalized effect. The answer therefore follows the smoother: an order-2
#' penalty leaves one column free, an order-3 penalty two, and a periodic
#' basis none, because its null space contains only the constant and the
#' constant belongs to the model's intercept. Under a factor `by` the block is
#' one copy per level and the free columns of every level are marked, whether
#' the levels share a smoothing parameter or not.
#'
#' @param term A built smooth term.
#' @param k The number of columns in its block.
#'
#' @return A logical vector of length `k`.
#'
#' @keywords internal
smooth_linear_cols <- function(term, k) {
  out <- rep(FALSE, k)
  if (k < 1L) return(out)
  # ASKED OF THE ENTRIES AND NOT OF ONE PENALTY. A smooth declares one
  # penalty over its whole block by default, and there `term_penalty()`
  # answers the same thing; but it declares one per level of a factor `by`
  # under by_hyper = "level", and one over the penalized coordinates ALONE
  # where the smoother carries a penalty factory, and in neither case is
  # there a single penalty to fetch. Read through the singular question a
  # factory smooth reported its smoothing parameter and nothing else.
  ent <- tryCatch(modelterms7::term_penalties(term), error = function(e) list())
  if (!length(ent)) return(out)
  covered <- rep(FALSE, k)
  # `en` and not `e`: the handler of the tryCatch below binds `e` to the
  # condition, which would shadow the entry for anything written inside it
  for (en in ent) {
    idx <- en$index[en$index >= 1L & en$index <= k]
    if (!length(idx)) next
    pen <- en$penalty
    p <- if (!is.null(pen) && "P" %in% S7::prop_names(pen)) {
      tryCatch(as.matrix(S7::prop(pen, "P")), error = function(e) NULL)
    } else if (!is.null(pen) && "mats" %in% S7::prop_names(pen)) {
      # AN ADDITIVE PENALTY CARRIES ITS MATRICES UNDER ANOTHER NAME, one per
      # smoothing parameter, and a coordinate is free where NO component
      # touches it. Summing the absolute values answers that exactly and
      # without an argument about definiteness, which summing the components
      # themselves would need. Read through the "P" question alone an
      # adaptive smooth fell into the branch below, written for a penalty
      # that carries no matrix at all, and reported every coordinate as
      # shrunk -- losing the linear row, which is the same defect the
      # singular question produced before this one.
      tryCatch(as.matrix(Reduce(`+`, lapply(S7::prop(pen, "mats"), abs))),
               error = function(e) NULL)
    } else {
      NULL
    }
    if (is.null(p) || !nrow(p)) {
      # a penalty carrying no matrix -- a lasso, a heavy-tailed prior --
      # shrinks every coordinate it indexes, and indexes only its own
      covered[idx] <- TRUE
      next
    }
    hit <- colSums(abs(p)) != 0
    # a quadratic penalty may be ONE BLOCK repeated over the levels, its
    # matrix narrower than the coordinates it covers
    covered[idx] <- if (length(idx) %% length(hit) == 0L) {
      rep(hit, length.out = length(idx))
    } else {
      rep(TRUE, length(idx))
    }
  }
  # EVERY column no penalty touches, not only the leading run. Under a factor
  # `by` the block is one copy per level, level-major, so the free straight
  # line of level j sits at the start of its copy and not at the start of the
  # block; reading the leading run alone printed the line of the first level
  # and dropped the others, which on four diets left three of the four
  # slopes a reader of the summary wants off the page.
  !covered
}


#' A Summary of a Fitted Model
#'
#' @description
#' Holds what [summary.StatmodFit()] returns: the blocks of each
#' distribution parameter, the degrees of freedom, the information criteria
#' and the notes that qualify the numbers.
#'
#' @param call The fit's call.
#' @param distrib_name The distribution's name.
#' @param n_obs The number of rows of the data.
#' @param weight_total The sum of the prior weights, which [nobs()] reports.
#' @param tables A named list, one entry per distribution parameter, each a
#'   list of block records.
#' @param classes The covariance blocks that two or more pieces share (terms,
#'   or developments of one term), one block record each. They belong to no
#'   equation and are printed ahead of the equations.
#' @param links The link of each distribution parameter, a character vector
#'   named by parameter.
#' @param edf The per-term degrees of freedom.
#' @param structural The table of the structural term's parameters, or `NULL`.
#' @param loglik The conditional log-likelihood at the fit.
#' @param df The effective degrees of freedom in total, including the
#'   correction for estimated hyperparameters where it was requested.
#' @param aic,bic The information criteria.
#' @param converged Whether every loop stopped on its own rule.
#' @param elapsed The elapsed time, in seconds.
#' @param level The confidence level at which the intervals were built.
#' @param vcov_type Which variance convention was used.
#' @param test Which statistic the coefficient tables report: `"wald"`,
#'   `"lr"`, `"score"` or `"gradient"`.
#' @param max_coef The number of coefficient rows that each block prints, or
#'   `NULL` where the option `statmodels7.summary_max_coef` determines it.
#' @param notes A character vector of the notes that qualify the numbers.
#' @param certificate The list returned by [statmod_certificate()], or `NULL`
#'   where it could not be computed.
#'
#' @return An object of class `StatmodSummary`.
#'
#' @seealso [summary.StatmodFit()]
#'
#' @examples
#' dd <- data.frame(y = rnorm(30), x = runif(30))
#' S7::S7_inherits(summary(statmod(y ~ x,
#'                                 distributions7::gaussian1_distrib(), dd)),
#'                 StatmodSummary)
#'
#' @name StatmodSummary-class
#' @aliases StatmodSummary
#' @keywords internal
#' @export
StatmodSummary <- S7::new_class("StatmodSummary",
  properties = list(
    call = S7::class_any,
    distrib_name = S7::class_character,
    n_obs = S7::class_numeric,
    weight_total = S7::class_numeric,
    tables = S7::class_list,
    # THE BLOCKS THAT BELONG TO NO EQUATION: a covariance a label shares
    # between terms, and where the label crosses a bar, between distribution
    # parameters. Its coordinates are effects of several equations at once,
    # so it is not a block of any of them and is carried on its own.
    classes = S7::new_property(S7::class_list, default = list()),
    # the link each equation is written on, named by parameter. Every
    # coefficient in a block is a coefficient of the LINEAR PREDICTOR, so
    # what it means for the parameter depends on the link, and a summary
    # that does not say which one leaves the reader to remember the
    # family's defaults.
    links = S7::new_property(S7::class_character, default = character(0)),
    edf = S7::class_any,
    structural = S7::class_any,
    loglik = S7::class_numeric,
    df = S7::class_numeric,
    aic = S7::class_numeric,
    bic = S7::class_numeric,
    converged = S7::class_logical,
    elapsed = S7::class_numeric,
    level = S7::class_numeric,
    # the variance convention of the standard errors, as the object that was
    # passed; print() reads its name
    vcov_type = vcov_method,
    # which test the statistic column reports. The default is Wald's, which
    # is the only one read off the unrestricted fit; the other three cost a
    # refit per row and are reported as a signed root, so the column keeps
    # one meaning and the header says which. Stored as the object, which
    # carries the label that print() shows.
    test = S7::new_property(test_statistic, default = quote(wald_test())),
    # how many coefficient rows a block prints, or NULL where the caller
    # said nothing and the option decides. It is carried on the object
    # rather than asked for at print time because an R session prints a
    # summary by evaluating it, so summary() is the only place a reader
    # ever passes an argument.
    max_coef = S7::class_any,
    notes = S7::class_character,
    # the reading of statmod_certificate(), or NULL where it could not be
    # taken. A list rather than columns: its four readings are of different
    # kinds and one of them is a set of names.
    certificate = S7::new_property(S7::class_any, default = NULL),
    # the information object that the standard errors invert, after the
    # variance object's NULL has been resolved against the family; print()
    # names it when it is not the observed one
    information = S7::new_property(S7::class_any, default = NULL)
  )
)


#' @title Summarize a Fitted Model
#' @name summary.StatmodFit
#' @description
#' Returns one coefficient table per distribution parameter, carrying the
#' estimate, its standard error, the statistic of the test named by `test`
#' (Wald's by default), its p-value and an interval, together with the degrees
#' of freedom, the information criteria and the notes that qualify the numbers.
#' @details
#' Each distribution parameter is read as blocks, not as one list of
#' coefficients, because most coefficients of a fitted model are not
#' quantities of interest in themselves. The blocks are
#' \describe{
#'   \item{the parametric terms}{every unpenalized term together, one row per
#'     coefficient, which is the ordinary table.}
#'   \item{one block per smooth}{the linear component's coefficient where the
#'     construction carries one, the smoothing parameter, and the effective
#'     degrees of freedom. The coefficients of the wiggly part are not shown:
#'     they are coordinates in an orthonormal basis, a single one of them has
#'     no interpretation, and together they determine the edf.}
#'   \item{one block per random effect}{the parameters of the distribution of
#'     the effects, usually called the variance components, the edf and a
#'     summary of the predicted effects. The effects themselves, of which there
#'     is one per level, are not listed.}
#'   \item{one block per selection}{a lasso, a SCAD or an MCP: its
#'     hyperparameters, how many coefficients survived, and those coefficients.
#'     The ones set exactly to zero are counted, not listed.}
#'   \item{one block per other penalized term}{its coefficients, which stay
#'     interpretable under a ridge, together with its hyperparameters.}
#'   \item{one block per break-point term}{the quantities that the term reports
#'     (slopes, changes and positions), with standard errors from the delta
#'     method applied to the coefficients through which the term is fitted. A
#'     position carries no test.}
#'   \item{one block per structural term}{the term's own parameters (a
#'     filter's level, loadings and persistences, or the transition
#'     probabilities of a regime chain), with a level that an intercept holds
#'     marked as held.}
#' }
#'
#' A hyperparameter is the first row of its block, since it governs
#' every coefficient under it, and its name carries a mark in square brackets
#' that gives the source of the value. A hyperparameter estimated by [reml()]
#' or [ml()] maximizes a twice differentiable criterion, so it carries a
#' standard error and an interval, both read on the free scale that its link
#' defines and mapped back ([statmod_hyper_vcov()]). One chosen by [aic()],
#' [bic()] or [cv()] along a path over a kinked penalty is the argument of a
#' minimum over a grid, so the row names the criterion and leaves the
#' remaining columns empty: there is no curvature at such a point from which
#' to read a standard error. A smooth penalty under [aic()] or [bic()] is
#' minimized over its free scale, and its row also has no standard error,
#' because the criterion is not a likelihood. A hyperparameter that the term
#' holds is marked fixed, and one that stayed at the starting value of its
#' term because the outer criterion was not run is marked start.
#'
#' The statistic column reports the test named by `test`. By default it is
#' Wald's, the estimate over its standard error, which is the only one of the
#' four that needs no refit and is read only from the unrestricted fit. The
#' other three -- the likelihood ratio, Rao's score and Terrell's gradient --
#' hold the coefficient at the null value and refit, one refit per row, and
#' are reported as the signed root of their \eqn{\chi^2_1} statistic with the
#' sign of \eqn{\hat\beta_j - b}, so that the column carries one kind of
#' quantity whichever test was requested and Wald's own \eqn{z} is the case
#' where that root is exact. The p-value is the \eqn{\chi^2_1} one in every
#' case, and for Wald the two readings agree identically.
#'
#' Only the statistic and its p-value change with `test`: the interval stays
#' the Wald one, since inverting a test requires many refits for each
#' coefficient. [confint.StatmodFit()] inverts the same test through its
#' `test` argument, one coefficient at a time. A row for which the restricted
#' fit is not defined -- a coefficient under a kinked penalty, an aliased one,
#' any coefficient of a model carrying a structural term -- reports `NA`
#' instead of the Wald value, so the column never carries two tests at once.
#' Every other row is tested, including a smooth's own coordinates and a random
#' effect's, on the penalized objective and with the reading described in
#' [statmod_stat_at()]. A quantity of a break-point term that is one
#' coefficient carries the restricted statistic of that coefficient; a
#' quantity that combines several coefficients, and the position of a
#' break-point, carry none.
#'
#' The meaning of a Wald p-value depends on the row. For an unpenalized
#' coefficient it is the ordinary p-value. For a coefficient in a penalized
#' block it is conditional on the smoothing parameter, which was not estimated
#' jointly with it, and it does not account for the shrinkage of the estimate
#' toward zero. For a coefficient that a kinked penalty selected, the row
#' exists only because that coefficient survived the selection, and a naive
#' interval there under-covers.
#'
#' The degrees of freedom are the effective ones, summed over the terms,
#' so a penalized term counts the degrees of freedom that it uses instead of
#' the number of its columns. The information criteria are built on that
#' count: `aic` is \eqn{-2\ell + 2\,\mathrm{df}} and `bic` replaces the factor
#' 2 by the logarithm of the sum of the prior weights, which are frequency
#' weights.
#'
#' # Reading the foot
#'
#' Below the criteria the summary describes the point at which the fit stopped
#' with a few numbers from [statmod_certificate()]. The summary prints no
#' verdict. Each number is dimensionless, so it does not change with the units
#' of a coefficient. \describe{
#'   \item{`inner   max |grad|/se`}{the largest
#'     \eqn{|g_j|/\sqrt{A_{jj}}} over the coefficients, where \eqn{g} is the
#'     gradient of the penalized objective and \eqn{A} its curvature, the
#'     penalized information. \eqn{1/\sqrt{A_{jj}}} is the standard error of
#'     coefficient \eqn{j} conditional on the others, which is smaller than
#'     the marginal one that the tables print, and the ratio equals
#'     \eqn{\sqrt{2\Delta_j}} with \eqn{\Delta_j = g_j^2/(2A_{jj})} the gain
#'     in log-likelihood of a Newton step in that coordinate alone. At a mode
#'     it is small. For a model carrying a score-driven filter it covers the
#'     filter's own parameters as well. It leaves out a coefficient that a
#'     kinked penalty has set to zero, a frozen working block, and an aliased
#'     column, and it is read on the directions that the kinks of a
#'     break-point sitting on an observation leave free.}
#'   \item{`inner   min eigen`}{the smallest eigenvalue of
#'     \eqn{D^{-1/2}AD^{-1/2}} with \eqn{D = \mathrm{diag}(A)}, the curvature
#'     rescaled to a unit diagonal. It is at most one. Near one the
#'     coordinates are well separated; near zero some combination of them is
#'     barely identified, as two nearly collinear covariates or an intercept
#'     beside a random effect make it; below zero the point is not a mode.}
#'   \item{`zeros   max |score|/kink`}{printed only where a kinked penalty
#'     (a lasso, an elastic net, a SCAD or an MCP) set coefficients exactly to
#'     zero, together with how many. At a zero the penalty has no derivative
#'     and the coefficient no standard error, so the condition for an optimum
#'     is \eqn{|s_j| \le \kappa_j}: the derivative of the log-likelihood in
#'     that coefficient within the size of the kink, \eqn{\lambda} for a
#'     lasso. The number is the largest \eqn{|s_j|/\kappa_j}. At or below one
#'     every zero is an optimum; above one a coefficient would improve the
#'     objective by leaving zero. A reading \eqn{r} below one means that the
#'     coefficient with the largest ratio is held by the kink with a margin of
#'     \eqn{1 - r}; if the other coefficients were fixed it would leave zero
#'     at a \eqn{\lambda} smaller by the factor \eqn{r}. Along a path the others
#'     move too, and the threshold moves with them. The number is not a test
#'     of significance, only a check on whether the zeros agree with the
#'     \eqn{\lambda} of the fit, and values close to one are ordinary where
#'     \eqn{\lambda} was chosen on a grid, because the next coefficient to
#'     enter is close to its threshold by construction.}
#'   \item{`outer   max |grad|/se` and `min eigen`}{the same two readings for
#'     the criterion that estimated the hyperparameters, on the free scale
#'     that their links define, printed only where [reml()] or [ml()] ran. The
#'     curvature is the negated Hessian of the criterion. A hyperparameter
#'     that ran to the edge of its chart is left out, because its curvature
#'     collapses with its gradient; where every one did the line reads
#'     `not available: at a boundary`. A hyperparameter chosen along a path
#'     has no gradient, and there is no outer line.}
#' }
#' @param object A [StatmodFit()].
#' @param level The confidence level.
#' @param test A [test_statistic] object, the statistic that the coefficient
#'   tables report. [wald_test()] (the default) reads the variance matrix and
#'   refits nothing. [lr_test()], [score_test()] and [gradient_test()] cost one
#'   restricted refit per row and are reported as a signed root. The four are
#'   compared at [test_statistic]; see also [statmod_stat_at()]. A character
#'   string is rejected with an error.
#' @param vcov_type A [vcov_method] object, the variance matrix, passed to
#'   [vcov.StatmodFit()]: [bayesian_vcov()] (the default),
#'   [frequentist_vcov()] or [unconditional_vcov()]. The information from which
#'   the standard errors are computed is an argument of each of them
#'   ([distributions7::observed_info()] by default). The last carries the
#'   uncertainty of the hyperparameters into every standard error under it,
#'   where `correct` carries the same uncertainty into the degrees of freedom;
#'   the two are the same quantity read on two surfaces, and both can be used
#'   in one call. The conventions are compared at [vcov_method].
#' @param correct Whether the degrees of freedom include the cost of
#'   estimating the hyperparameters. The ordinary count treats them as known
#'   although they were estimated from the same data, so a criterion built on
#'   it is too generous. See [statmod_edf_correction()]. The default is
#'   `FALSE`, which leaves the count unchanged; the correction is zero where no
#'   hyperparameter was estimated by a marginal criterion.
#' @param max_coef How many coefficient rows each block of the printed
#'   summary shows, `Inf` or `NA` for all of them. `NULL`, the default,
#'   reads the option `statmodels7.summary_max_coef`, and 10 where that is
#'   unset. A hyperparameter row is never hidden, whatever this is, and a
#'   block of twelve rows or fewer is never abridged. The value is carried
#'   on the result, so `summary(fit, max_coef = Inf)` prints in full at the
#'   console; [print.StatmodSummary()] takes the same argument for an object
#'   already in hand.
#' @param ... Passed to [vcov.StatmodFit()].
#' @return A [StatmodSummary()].
#' @seealso [vcov.StatmodFit()], [confint.StatmodFit()], [test_statistic],
#'   [vcov_method], [statmod_test()], which computes the same four statistics for one
#'   coefficient against a value that need not be zero, and
#'   [statmod_stat_at()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(120))
#' dd$y <- 1 + 2 * dd$x + rnorm(120, sd = 0.4)
#' summary(statmod(y ~ x | sigma ~ x,
#'                 distributions7::gaussian1_distrib(), dd))
#' @keywords internal
summary.StatmodFit <- function(object, level = 0.95,
                               vcov_type = bayesian_vcov(),
                               correct = FALSE,
                               test = wald_test(),
                               max_coef = NULL,
                               ...) {
  vcov_type <- as_vcov_method(vcov_type)
  test <- as_test_statistic(test)
  # A TERM WHOSE BLOCK IS A WORKING LINEARIZATION says so ONCE, as a note,
  # rather than once per call to vcov() -- of which this function makes more
  # than one. The condition carries a class of its own, so muffling it cannot
  # swallow another warning raised on the way. A variance that could not be
  # made unconditional is caught the same way and for the same reason.
  frozen_msg <- character(0)
  catch_frozen <- function(expr) {
    withCallingHandlers(expr,
      statmod_frozen_block = function(cnd) {
        frozen_msg <<- unique(c(frozen_msg, conditionMessage(cnd)))
        invokeRestart("muffleWarning")
      },
      statmod_conditional_variance = function(cnd) {
        frozen_msg <<- unique(c(frozen_msg, conditionMessage(cnd)))
        invokeRestart("muffleWarning")
      })
  }
  # ON THE RAW COORDINATES, both of them. This function does its own reading
  # of what a term is about -- readable_rows() for a design block, the
  # structural table for a term with none -- and it INDEXES the variance
  # matrix by the design's own coefficient names. Given the readable matrix
  # those keys are not in it, and every delta-method standard error a
  # break-point term reports came back missing.
  ci <- catch_frozen(confint(object, level = level, vcov_type = vcov_type,
                             readable = FALSE, ...))
  spec <- object@spec
  design <- statmod_design(spec)
  # A HELD BREAK-POINT says once why its position has no standard error
  held_terms <- unlist(lapply(names(spec@terms), function(p) {
    Filter(function(nm) length(tryCatch(
      modelterms7::term_held(spec@terms[[p]][[nm]]),
      error = function(e) integer(0))) > 0L, names(spec@terms[[p]]))
  }))
  if (length(held_terms)) {
    frozen_msg <- c(frozen_msg, sprintf(paste0(
      "The break-points of %s are held at the minimum of their profile,",
      " which is\n  constant between consecutive observations and has no",
      " curvature, so they carry\n  no standard error and the other",
      " coefficients are conditional on them. A\n  smoothed term",
      " (smoothed = numericals7::smooth_probit()) gives them one."),
      paste(unique(held_terms), collapse = ", ")))
  }
  ci$statistic <- ci$estimate / ci$se
  ci$p_value <- 2 * stats::pnorm(-abs(ci$statistic))
  test_msg <- character(0)
  if (test_restricted(test)) {
    rs <- restricted_stat_rows(object, ci, test, vcov_type, spec, design)
    ci$statistic <- rs$statistic
    ci$p_value <- rs$p_value
    test_msg <- rs$note
  }
  lab <- coef_labels(spec, design)

  # READ BEFORE THE BLOCKS ARE BUILT, and not only because a note below asks
  # whether the point is a maximum: a hyperparameter the certificate finds at
  # a boundary carries no standard error, so the blocks have to be given the
  # verdict rather than reached after it. It costs one outer gradient and one
  # solve, which is nothing beside a summary that already inverts the
  # penalized information, and a fold of cv() or a path point that never
  # prints itself pays nothing at all.
  cert <- tryCatch(statmod_certificate(object), error = function(e) NULL)
  at_edge <- if (is.null(cert)) character(0) else cert$boundary_key

  # one variance matrix for every block that needs one: a term reported by
  # the quantities it is about carries them across by the delta method
  V <- catch_frozen(tryCatch(vcov(object, readable = FALSE, vcov_type = vcov_type,
                                  ...),
                             error = function(e) NULL))
  # A RIDGE OF THE OUTER CRITERION costs every coefficient that moves along it
  # its standard error, its interval and its test: the criterion does not
  # change along the direction, so the variance there describes rounding. The
  # coefficients the criterion estimated are named by the certificate; the
  # ones held in the inner fit are those the variance matrix moves with the
  # ridge's combination (measured on a Burr XII whose shape runs to its
  # Weibull limit in one level of a factor: the three coefficients of log k
  # and, with them, the three of log mu, by 1/c of each).
  ridge_rows <- if (is.null(cert) || !length(cert$ridge)) integer(0) else
    tryCatch(ridge_coef_rows(cert$ridge, V, rownames(ci)),
             error = function(e) integer(0))
  ridge_msg <- character(0)
  if (length(ridge_rows)) {
    ci[ridge_rows, c("se", "statistic", "p_value", "lower", "upper")] <-
      NA_real_
    ridge_msg <- sprintf(paste0(
      "The outer criterion is flat along %s (%s), so the coefficients that ",
      "move\n  along it carry no standard error, interval or test: %s. ",
      "Their estimates are\n  one point of the ridge. A ridge of this kind ",
      "appears when a parameter runs to\n  a limit at which the family ",
      "becomes another one."),
      if (length(cert$ridge) == 1L) "a direction" else "several directions",
      paste(vapply(cert$ridge, function(x) paste0(
        ifelse(x < 0, "- ", "+ "), format(abs(x), digits = 3), " ",
        names(x), collapse = " "), ""), collapse = "; "),
      paste(rownames(ci)[ridge_rows], collapse = ", "))
  }
  # A COEFFICIENT THAT THE CRITERION ESTIMATED AND LEFT AT A BOUNDARY has no
  # curvature there, so its standard error describes rounding: measured on a
  # Student t fitted to gaussian data, the intercept of log nu stopped at
  # 20.8 with a standard error of 1e4 and an interval from -2e4 to 2e4.
  edge_rows <- if (is.null(cert) || !length(cert$boundary_coef)) integer(0)
    else stats::na.omit(match(cert$boundary_coef, rownames(ci)))
  edge_msg <- character(0)
  if (length(edge_rows)) {
    ci[edge_rows, c("se", "statistic", "p_value", "lower", "upper")] <-
      NA_real_
    edge_msg <- sprintf(paste0(
      "%s %s estimated on the outer criterion and %s at a boundary of it: ",
      "the
  criterion has no curvature there, so %s no standard error, ",
      "interval or test.
  A boundary of this kind appears when a parameter ",
      "runs to a limit at which the
  family becomes another one."),
      paste(rownames(ci)[edge_rows], collapse = ", "),
      if (length(edge_rows) == 1L) "was" else "were",
      if (length(edge_rows) == 1L) "sits" else "sit",
      if (length(edge_rows) == 1L) "it carries" else "they carry")
  }
  strc <- tryCatch(statmod_structural_table(object, level),
                   error = function(e) NULL)
  tables <- lapply(spec@distrib@params, function(p)
    summary_blocks(object, spec, design, p, ci, level, V, strc, at_edge,
                   test))
  names(tables) <- spec@distrib@params
  # a covariance shared between terms is the property of none of them, and
  # summary_blocks() held its rows apart for that reason
  cls <- tryCatch(summary_class_blocks(spec, design, tables, object@edf,
                                       object@coefficients, object@hyper),
                  error = function(e) list())
  # A LINE SAYING "REPORTED ABOVE" IS FALSE WHERE THE SECTION IS NOT THERE.
  # It cannot happen as the two are written -- every member of a shared class
  # reaches the block builder and exactly one of them carries the rows -- and
  # a member that has said where its numbers are and does not show them is
  # the worst of the three states, so the two views are made to agree here
  # rather than left to agree by construction. A class block's own note is
  # its legend and carries no key.
  shown <- vapply(cls, function(z) z$term, "")
  tables <- lapply(tables, function(bl) lapply(bl, function(b) {
    if (!length(b$note) || is.null(names(b$note))) return(b)
    b$note <- unname(b$note[names(b$note) %in% shown])
    b
  }))

  ll <- logLik.StatmodFit(object)
  df <- attr(ll, "df")
  aliased_msg <- if (length(object@aliased)) {
    sprintf(paste0("%d coefficient%s not estimated: the design carries the",
                   " same information in\n  other columns, so the pivot that",
                   " fitted the model left %s out. %s estimate,\n  standard",
                   " error and interval are reported as missing (%s)."),
            length(object@aliased),
            if (length(object@aliased) == 1L) " is" else "s are",
            if (length(object@aliased) == 1L) "it" else "them",
            if (length(object@aliased) == 1L) "Its" else "Their",
            paste(object@aliased, collapse = ", "))
  } else character(0)
  notes <- c(character(0), frozen_msg, aliased_msg, ridge_msg, edge_msg,
             test_msg,
             tryCatch(class_notes(spec, design), error = function(e) character(0)))
  # WHERE THE PARAMETERS ENDED UP, for a fit that did not converge. It
  # qualifies the fit rather than describing it, so it is a note and not a
  # line of print(): the measured case is a scale that ran to 1e-15 while
  # 380 of 400 coefficients survived, which is read once when something
  # looks wrong and never otherwise.
  # ⚠️ ASKED OF THE CERTIFICATE, as the note on the maximum below is, and for
  # the same reason: a regime chain whose transition coordinate ran to the
  # edge (geyser, alr1.1 at -27) leaves the search's flag FALSE at a point
  # the certificate reads as converged, and the note said "The fit did not
  # converge" under that verdict.
  cert_state <- if (is.null(cert) || is.null(cert$state)) NA_character_ else
    cert$state
  cert_bad <- is.na(cert_state) || !cert_state %in% c("converged", "boundary")
  # a certificate that is UNKNOWN by construction (a criterion that is not a
  # marginal likelihood, hyperparameters held, a path, no criterion) says
  # nothing against the fit's own flag, so the ranges are printed only where
  # the certificate read the point as not converged, or where there is no
  # certificate and the flag is FALSE
  not_conv <- identical(cert_state, "not converged") ||
    ((is.na(cert_state) || identical(cert_state, "unknown")) &&
       !isTRUE(object@converged))
  if (not_conv) {
    r <- tryCatch(fitted_ranges(object), error = function(e) "")
    if (nzchar(r)) notes <- c(notes, r)
  }

  # The count above reads the hyperparameters as though they were known,
  # and they were estimated from the same data. Adding what that costs is
  # off by default because it changes a number a reader may be comparing
  # with an earlier fit.
  # the matrix the standard errors invert: the one that the variance object
  # names, or the family's default where it names none, and the observed one
  # where a structural term has no expected information
  info_used <- if (length(attr(design, "structural"))) {
    distributions7::observed_info()
  } else vcov_information(vcov_type, spec@distrib)
  corr <- 0
  if (isTRUE(correct)) {
    cc <- tryCatch(statmod_edf_correction(spec, object@coefficients,
                                          object@hyper, design,
                                          object@methods$outer,
                                          expected = info_expected(info_used),
                                          emethod = info_emethod(info_used),
                                          pinned = object@methods$pinned),
                   error = function(e)
                     list(total = 0, per = numeric(0), n_hyper = NA_integer_))
    corr <- cc$total
    df <- df + corr
    nh <- if (is.null(cc$n_hyper)) 0L else cc$n_hyper
    # an error is a correction that could not be computed, not one that is
    # zero for want of a hyperparameter
    if (is.na(nh)) nh <- 1L
    notes <- c(notes, if (corr > 0) sprintf(paste0(
      "The degrees of freedom carry %.3g for the hyperparameters, which ",
      "were\n  estimated rather than given. Without it the criteria are ",
      "too generous."), corr) else if (nh > 0L) paste0(
      "The correction for the estimated hyperparameters could not be ",
      "computed\n  here, so the degrees of freedom read them as though they ",
      "were given.") else paste0(
      "A marginal criterion estimated no hyperparameter here, so ",
      "there is\n  nothing for the correction to propagate and it is zero."))
  }
  if (any(lab$penalized)) {
    # WHICH criteria were in force, so a reader knows whether the number at
    # the head of a penalized block was chosen or given. The two kinds carry
    # different guarantees and are named separately.
    # ⚠️ A SHARED BLOCK'S ROWS ARE HELD APART, so a walk over `tables` alone
    # does not reach them -- and a fit whose only penalized rows are a
    # covariance a label collects printed NONE of the three notes below, not
    # even the one saying its hyperparameters carry a standard error at all.
    # Measured before this line existed: such a summary carried one note, the
    # one naming the shared block. The two are gathered once here and every
    # question below is asked of both.
    tabs <- c(unlist(lapply(tables, function(bl) lapply(bl, `[[`, "table")),
                     recursive = FALSE),
              lapply(cls, `[[`, "table"))
    src <- unique(unlist(lapply(tabs, function(tb)
      tb$source[tb$role %in% c("fixed", "estimated")]), use.names = FALSE))
    marg <- setdiff(intersect(src, c("reml", "ml")), NA)
    # a PATH chooses the hyperparameters of a kinked penalty, and the sparse
    # criterion records which it reached; a smooth penalty under aic() or
    # bic() is minimized by a search over its free scale instead, which is
    # not a grid
    spc <- object@methods$sparse_criterion
    path <- if (length(object@methods$sparse_hyper) && !is.null(spc))
      spc@kind else character(0)
    om <- object@methods$outer
    pu <- tryCatch(statmod_penalized(spec, design), error = function(e) list())
    smooth_est <- any(vapply(pu, function(u)
      !isTRUE(tryCatch(penalty_has_kink(u$penalty, u$key),
                       error = function(e) FALSE)) &&
        length(setdiff(u$penalty@params, names(u$fixed))) > 0L, logical(1)))
    search <- if (!is.null(om) && om@kind %in% c("aic", "bic") && smooth_est)
      om@kind else character(0)
    # WHETHER THE STANDARD ERROR IS THERE, and not merely whether the
    # criterion was marginal. A hyperparameter the search left before
    # reaching a maximum has a curvature of the wrong sign in its own
    # direction and no variance follows from it, so the row carries the
    # estimate and nothing else -- and a note promising an interval beside
    # it would be false of the fit in front of the reader.
    mse <- unlist(lapply(tabs, function(tb) {
      r <- tb$role %in% "estimated" & tb$source %in% c("reml", "ml")
      if (!any(r)) return(logical(0))
      is.finite(tb$se[r])
    }), use.names = FALSE)
    if (length(marg) && any(mse)) {
      notes <- c(notes, sprintf(paste0(
        "A hyperparameter marked %s was estimated by that criterion, and its",
        " standard\n  error and interval are read on the free scale that its link ",
        "defines and mapped\n  back. Every coefficient beside it is still ",
        "conditional on the value reached."),
        paste(toupper(marg), collapse = " or ")))
    }
    # ⚠️ THE CURVATURE IS NOT THE ONLY REASON A STANDARD ERROR IS MISSING, and
    # this note used to be the only one offered. A coordinate at a boundary has
    # a curvature that is a proper maximum, so printing the sentence below over
    # its blank cell would state the opposite of what happened -- the shape
    # this file already records for a note that contradicted the table beside
    # it. Where a boundary coordinate exists the boundary note answers instead.
    # ⚠️ EXCEPT WHERE NOTHING AT ALL CARRIES A STANDARD ERROR, which a
    # boundary cannot explain: it holds the coordinates it names and leaves
    # the rest of the matrix computed, so a table in which every row is blank
    # was refused for the other reason and the sentence below is the one that
    # answers for it. Measured on a class at the boundary whose outer
    # curvature `hyper_variance()` could not read at all: three blank rows,
    # of which the boundary note explains one. Both notes are printed there.
    if (length(marg) && sum(!mse) > length(at_edge)) {
      notes <- c(notes, sprintf(paste0(
        "A hyperparameter marked %s with no standard error was estimated at ",
        "a point\n  where the criterion's curvature in its own direction is ",
        "not negative: the\n  search did not reach a maximum there, so a ",
        "variance cannot be derived from it. Its\n  estimate stands and every ",
        "coefficient beside it is conditional on it."),
        paste(toupper(marg), collapse = " or ")))
    }
    if (length(marg) && any(!mse) && length(at_edge)) {
      notes <- c(notes, sprintf(paste0(
        "A hyperparameter marked %s that has run to an edge of its range ",
        "reports no\n  standard error, and the quantities that ",
        "depend on it do not report one either: the\n  criterion has stopped moving in that ",
        "direction, so the estimate is pinned\n  against the edge rather ",
        "than located by the data, and an interval around\n  it would ",
        "describe a curvature that is not the uncertainty. The rest of the\n",
        "  matrix is computed and reported as usual. The certificate names ",
        "the\n  coordinates: %s."),
        paste(toupper(marg), collapse = " or "),
        paste(cert$boundary, collapse = ", ")))
    }
    if (length(path)) {
      notes <- c(notes, sprintf(paste0(
        "A hyperparameter of a kinked penalty marked %s was chosen by a path ",
        "over its\n  own values. It is the argument of a minimum over a grid ",
        "rather than the root\n  of a derivative, so a standard error ",
        "cannot be derived from it; its uncertainty is a\n  resampling question."),
        toupper(path)))
    }
    if (length(search)) {
      notes <- c(notes, sprintf(paste0(
        "A hyperparameter of a smooth penalty marked %s was chosen by ",
        "minimizing that\n  criterion over its free scale. The criterion is ",
        "not a likelihood, so a\n  standard error cannot be derived from its ",
        "curvature, and every coefficient beside it\n  is conditional on ",
        "the value reached."), toupper(search)))
    }
    if ("fixed" %in% src) {
      notes <- c(notes, paste0(
        "A hyperparameter marked fixed is held at the value that it was given, not",
        " estimated,\n  so it has no standard error and every interval beside",
        " it is conditional on it."))
    }
    if ("start" %in% src) {
      notes <- c(notes, paste0(
        "A hyperparameter marked start stayed at the starting value of its ",
        "term, because\n  the fit ran no criterion to estimate it, so it has no ",
        "standard error and every\n  interval beside it is conditional on it."))
    }
  }
  if (any(lab$kinked)) {
    nz <- sum(lab$kinked &
                unlist(object@coefficients[spec@distrib@params],
                       use.names = FALSE) != 0)
    notes <- c(notes, sprintf(paste0(
      "%d of %d coefficients under a kinked penalty survived. Their rows ",
      "are\n  conditional on that selection, and the ones at zero have no ",
      "variance at\n  all: at the kink there is no curvature to read."),
      nz, sum(lab$kinked)))
  }
  # ⚠️ WHETHER THE POINT IS A MAXIMUM is the CERTIFICATE's answer and not the
  # search's flag, and this note used to read the flag: a fit whose search
  # stopped on its own rule at a point located to 2e-11 printed "the fit did
  # not converge, so everything above is read at a point that is not a
  # maximum" directly under a certificate saying CONVERGED. The two are
  # different questions, so the note is emitted only where the question it
  # answers -- is this a maximum -- has actually been answered no.
  if (cert_bad && !identical(cert_state, "unknown")) {
    notes <- c(notes, paste0(
      "The point reported is not certified as a maximum, so every estimate ",
      "and\n  interval above is read where the surface is still moving."))
  } else if (identical(cert_state, "unknown")) {
    why <- if (length(cert$reason)) paste(cert$reason, collapse = "; ") else
      "the certificate could not read it"
    notes <- c(notes, sprintf(paste0(
      "The point reported is not certified as a maximum of an outer ",
      "criterion:\n  %s."), why))
  }

  # A smoothed break-point term: the smoother and its width are part of the
  # model -- the transition's width, the bent-cable reading -- so they are
  # reported; and where a break-point carries a random development, a
  # smoother declaring a scale correction (the probit's convolution
  # identity) gets the corrected scale printed beside the apparent one.
  notes <- c(notes, tryCatch(smoothed_notes(spec, object),
                             error = function(e) character(0)))

  # A structural term contributes no columns, so its block is built from
  # what it reports rather than from a design; `strc` is computed above,
  # where summary_blocks() can be given it.
  if (!is.null(strc) && any(strc$held)) {
    notes <- c(notes, paste0(
      "A level marked held is carried by an intercept in the same equation ",
      "and is\n  not estimated: the two are exactly confounded, so only one ",
      "of them can be."))
  }
  if (!is.null(strc) && isTRUE(any(strc$edge))) {
    notes <- c(notes, sprintf(paste0(
      "%s depend%s on a coordinate that has run to the edge of its chart, ",
      "so %s reported\n  without a standard error and an interval: the ",
      "chart's derivative tends to\n  zero there and the delta method would ",
      "report a standard error that\n  collapses where the quantity stops ",
      "being identified."),
      paste(strc$name[strc$edge], collapse = ", "),
      if (sum(strc$edge) == 1L) "s" else "",
      if (sum(strc$edge) == 1L) "it is" else "they are"))
  }

  StatmodSummary(
    call = object@call, distrib_name = spec@distrib@distrib_name,
    n_obs = spec@n_obs, weight_total = weight_total(spec),
    tables = tables, classes = cls, edf = object@edf,
    links = vapply(spec@distrib@link_params,
                   function(g) g@link_name, character(1)),
    structural = strc,
    loglik = object@loglik, df = df,
    aic = -2 * object@loglik + 2 * df,
    bic = -2 * object@loglik + log(weight_total(spec)) * df,
    converged = object@converged, elapsed = object@elapsed,
    level = level, vcov_type = vcov_type, test = test, max_coef = max_coef,
    notes = notes, certificate = cert, information = info_used)
}
S7::method(summary, StatmodFit) <- summary.StatmodFit


#' A Table's Statistics From a Restricted Fit
#'
#' @description
#' Computes the likelihood-ratio, score or gradient statistic for each
#' coefficient that a restricted fit can hold, against the null that it is
#' zero, and reports it as the signed root of the \eqn{\chi^2_1} value
#' together with the p-value of that value.
#'
#' @details
#' The root carries the sign of the estimate, so the column holds one kind of
#' quantity whichever statistic produced it, and Wald's own \eqn{z} is the
#' case where the root is exact. The p-value is the \eqn{\chi^2_1} one and is
#' not recomputed from the root.
#'
#' A row that [testable_coords()] does not name keeps `NA` in both columns
#' instead of the Wald value with which it arrived, so that each column
#' carries a single test.
#'
#' @param fit A [StatmodFit()].
#' @param ci The flat interval table, with the Wald statistic already in it.
#' @param test A restricted [test_statistic]: [lr_test()], [score_test()] or
#'   [gradient_test()].
#' @param vcov_type The [vcov_method] that the Wald bracket reads, passed on.
#' @param spec,design The specification and its design.
#'
#' @return A list with `statistic`, `p_value` and `note`, the last a
#'   character vector of the notes that qualify the column.
#'
#' @seealso [statmod_stat_at()], which computes one of them.
#'
#' @keywords internal
restricted_stat_rows <- function(fit, ci, test, vcov_type, spec, design) {
  ok <- names(testable_coords(fit, spec, design))
  hit <- rownames(ci) %in% ok
  st <- rep(NA_real_, nrow(ci))
  pv <- rep(NA_real_, nrow(ci))
  bad <- 0L
  o0 <- if (identical(test@name, "lr") && any(hit)) joint_objective(fit) else NULL
  for (i in which(hit)) {
    s <- tryCatch(statmod_stat_at(fit, ci$parameter[[i]], ci$coefficient[[i]],
                               0, test, vcov_type, objective0 = o0),
                  error = function(e) NULL)
    if (is.null(s) || !isTRUE(is.finite(s$statistic))) next
    if (isFALSE(s$converged)) bad <- bad + 1L
    st[[i]] <- sign(ci$estimate[[i]]) * sqrt(max(s$statistic, 0))
    pv[[i]] <- s$p.value
  }
  nm <- c(lr = "likelihood-ratio", score = "score",
          gradient = "gradient")[[test@name]]
  note <- sprintf(paste0("The statistic is the %s one, the signed root of ",
                         "its chi-squared value\n  with the sign of the ",
                         "estimate, and its p-value is that value's. Each ",
                         "row costs\n  one refit in which the coefficient is held ",
                         "at zero, on the penalized objective and\n  at the ",
                         "hyperparameters at which the fit stopped, so a row ",
                         "under a penalty carries\n  the same conditional ",
                         "reading as its bayesian standard error. The\n  ",
                         "intervals beside it are still Wald's: ",
                         "confint(test =) inverts this test\n  instead."),
                 nm)
  n_no <- sum(!hit)
  if (n_no) {
    note <- c(note, sprintf(paste0("%d row%s no statistic: a coefficient ",
                                   "under a kinked penalty, an aliased one",
                                   "\n  and every coefficient of a model ",
                                   "carrying a structural term cannot be ",
                                   "held at a\n  value, so the %s statistic ",
                                   "is not defined for %s."),
                            n_no, if (n_no == 1L) " has" else "s have", nm,
                            if (n_no == 1L) "it" else "them"))
  }
  if (bad) {
    note <- c(note, sprintf(paste0("%d restricted refit%s not converge, so ",
                                   "the statistic read off %s is\n  the ",
                                   "value at the point where that refit ",
                                   "stopped, not at a maximum."),
                            bad, if (bad == 1L) " did" else "s did",
                            if (bad == 1L) "it" else "them"))
  }
  list(statistic = st, p_value = pv, note = note)
}


#' Which Readings a Chart's Coordinate Enters
#'
#' @description
#' Reads a readable block's Jacobian as a statement about dependence: `TRUE`
#' where a quantity depends on a coordinate, `FALSE` where the entry is the
#' floating-point image of a structural zero.
#'
#' @details
#' The reading replaces the test `== 0`, an exact comparison on a quantity
#' that the chart produces by arithmetic. On a covariance that a class
#' carries, whose chart is `parameters7::dr_prod(2)`, a standard deviation is
#' a function of its own log coordinate only, so its entry in the correlation
#' column is structurally zero. The entry that is computed can be a small
#' nonzero number from cancellation, and an exact comparison then reads it as
#' a dependence, which removes the standard error and the interval of the
#' standard deviation for a dependence that does not exist.
#'
#' A Jacobian entry is judged against the largest entry of its own row,
#' because a row is a gradient, and the tolerance is
#' `sqrt(.Machine$double.eps)`, the point beyond which a differentiated
#' quantity cannot be distinguished from the rounding of the quantities from
#' which it was formed. The two possible errors
#' are asymmetric, so the tolerance errs on the generous side: a genuine
#' dependence of relative size below `sqrt(eps)` that is read as zero changes
#' a variance contribution by a relative amount of order `1e-16`, whereas a
#' cancellation artifact read as a dependence removes the standard error of a
#' quantity.
#'
#' @param J The Jacobian of a readable block, quantities by coordinates.
#'
#' @return A logical matrix of the same shape, `TRUE` where the quantity in
#'   that row depends on the coordinate in that column. A row that is zero
#'   throughout depends on nothing and comes back all `FALSE`.
#'
#' @seealso [readable_hyper_rows()], [summary_blocks()]
#'
#' @keywords internal
jacobian_depends <- function(J) {
  if (is.null(J) || !length(J)) return(matrix(FALSE, 0L, 0L))
  J <- as.matrix(J)
  tol <- sqrt(.Machine$double.eps)
  scale <- apply(abs(J), 1L, max)
  scale[!is.finite(scale)] <- 0
  out <- abs(J) > tol * scale
  out[!is.finite(J)] <- TRUE
  dimnames(out) <- dimnames(J)
  out
}


#' The Quantities That a Penalty's Hyperparameters Describe
#'
#' @description
#' Replaces the coordinate rows of a penalty whose hyperparameters are a chart
#' with the quantities that the penalty declares through
#' [penalties7::penalty_readable()]: the standard deviations and
#' correlations of a correlated random effect, in place of the logarithms of
#' a Cholesky diagonal and the entries below it.
#'
#' @details
#' The standard error is the delta method, and it composes two Jacobians: the
#' penalty's, which is in the parameter scale of its hyperparameters, and the
#' link's, because the variance matrix is on the free scale on which the outer
#' criterion was maximized. Each interval is built on the scale that the
#' quantity declares, log for a standard deviation and Fisher's z for a
#' correlation, and mapped back, so a standard deviation cannot be given a
#' negative lower end and a correlation cannot be given an interval that
#' leaves \eqn{(-1, 1)}. The coefficients of a linear predictor keep the
#' symmetric interval of the estimate plus or minus \eqn{z} standard errors.
#'
#' The summary does not print a test here, for the same reason that it does
#' not print one for the coordinate rows: a \eqn{z} statistic tests the null
#' hypothesis that the quantity is zero, which for a standard deviation is the
#' edge of its range.
#'
#' @param rd The result of [penalties7::penalty_readable()].
#' @param th The penalty's hyperparameters, as fitted.
#' @param Vh The hyperparameter variance matrix, or `NULL`.
#' @param p The distribution parameter whose equation contains the term.
#' @param key The penalty's key.
#' @param level The confidence level.
#' @param role,src The role and the source reported for the coordinate rows.
#' @param cols The column names of a summary block.
#' @param labels One label per coordinate of the matrix, or nothing. A
#'   correlation between the second and the third coordinate of a covariance
#'   is a number about two named effects, and the family that carries the
#'   chart does not identify which; where the labels are given, they replace
#'   the positions in the printed name.
#' @param at_edge The keys of the hyperparameter coordinates that
#'   [statmod_certificate()] found at a boundary. A quantity whose
#'   Jacobian entry on one of them is not structurally zero reports no
#'   standard error and no interval; every other quantity keeps its own.
#'
#' @return A data frame of rows, in the shape of a summary block.
#'
#' @keywords internal
readable_hyper_rows <- function(rd, th, Vh, p, key, level, role, src, cols,
                                labels = character(0),
                                at_edge = character(0)) {
  nm <- names(th)
  k <- length(rd$value)
  se <- rep(NA_real_, k)
  kk <- paste(p, key, nm, sep = "\r")
  if (!is.null(Vh)) {
    j <- match(kk, rownames(Vh))
    if (!anyNA(j)) {
      vb <- as.matrix(Vh[j, j, drop = FALSE])
      lk <- attr(attr(Vh, "idx"), "links")[j]
      d <- vapply(seq_along(nm), function(i) {
        eta <- linkfunctions7::linkfun(lk[[i]], th[[i]])
        linkfunctions7::dlinkinv(lk[[i]], eta)
      }, 0)
      jac <- rd$jacobian * rep(d, each = k)
      # a coordinate at a boundary has no variance, statmod_hyper_vcov()
      # having held it; it contributes nothing here, and every quantity that
      # depends on it loses its standard error below
      edge <- kk %in% at_edge
      if (any(edge)) {
        vb[edge, ] <- 0
        vb[, edge] <- 0
        jac[, edge] <- 0
      }
      if (all(is.finite(vb)) && all(is.finite(jac))) {
        se <- sqrt(pmax(diag(jac %*% vb %*% t(jac)), 0))
      }
    }
  }
  # WHICH QUANTITIES A COORDINATE AT A BOUNDARY COSTS, and it is asked of the
  # Jacobian rather than of the chart: a quantity whose column there is
  # structurally zero does not depend on that coordinate and keeps its
  # standard error, which is the same question `keep` asks one caller up. On
  # the two-dimensional covariance a class carries, the correlation is the
  # only reading the angle enters, so the two standard deviations keep their
  # numbers and the correlation loses its own -- and a chart of another shape
  # is covered without an edit.
  # ⚠️ STRUCTURALLY zero, judged against the row's own scale and not by an
  # exact comparison: `jacobian_depends()` records what an exact test cost.
  held <- which(kk %in% at_edge)
  if (length(held) && !is.null(rd$jacobian)) {
    dep <- apply(jacobian_depends(rd$jacobian)[, held, drop = FALSE], 1L, any)
    se[dep] <- NA_real_
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  v <- as.numeric(rd$value)
  tr <- as.character(rd$transform)
  lo <- hi <- rep(NA_real_, k)
  for (i in seq_len(k)) {
    if (!is.finite(se[i])) next
    ends <- switch(tr[i],
      log = exp(log(v[i]) + c(-1, 1) * z * se[i] / v[i]),
      atanh = tanh(atanh(v[i]) + c(-1, 1) * z * se[i] / (1 - v[i]^2)),
      v[i] + c(-1, 1) * z * se[i])
    ends <- sort(ends)
    lo[i] <- ends[[1L]]
    hi[i] <- ends[[2L]]
  }
  out <- data.frame(name = readable_coord_names(names(rd$value), labels),
                    estimate = v, se = se,
                    statistic = NA_real_, p_value = NA_real_,
                    lower = lo, upper = hi,
                    role = rep(role[[1L]], k), source = rep(src[[1L]], k),
                    stringsAsFactors = FALSE)
  stats::setNames(out, c(cols, "source"))
}

#' The Blocks of One Distribution Parameter
#'
#' @description
#' Groups the terms of a distribution parameter into the blocks that a
#' summary prints: the parametric terms together, and one block for each other
#' term (smooth, random effect, selection, other penalized term, break-point
#' term or structural term).
#'
#' @param fit A [StatmodFit()].
#' @param spec The specification.
#' @param design The design.
#' @param p The distribution parameter.
#' @param ci The flat interval table that [confint.StatmodFit()] returns,
#'   to which the statistic and the p-value have been added.
#' @param level The confidence level at which the intervals are built.
#' @param V The variance matrix over the stacked coefficients, or
#'   `NULL`. It is needed only by a term reported through
#'   [modelterms7::term_readable()], whose quantities are
#'   functions of several coefficients at once and whose standard errors
#'   therefore come from the delta method and not from a diagonal entry.
#' @param test Which test the statistic column reports. Under a test other
#'   than `"wald"`, a quantity of a break-point term that is one coefficient
#'   carries the restricted statistic of that coefficient, and any other
#'   quantity carries none.
#'
#' @return A list of block records, each with `kind`, `label`, `term`,
#'   `n_coef`, `edf`, `n_zero`, `table` and `components`, and, where they
#'   apply, `classes` and `note` (the covariance blocks and the note that
#'   locates them), `n_na` (the aliased coefficients of a selection) and
#'   `spread` (the summary of the predicted effects of a random effect). A
#'   term written in parameters of its own that develops one of them over
#'   covariates reports that parameter as a compartment of its own, carrying
#'   its hyperparameter and the rows of its sub-terms, and `table` keeps only
#'   what is left. Such a term with a kinked sub-term also carries
#'   `n_kinked`, the number of coordinates under the kink, which `n_zero` and
#'   the count of selected coefficients read. The element `head` is always
#'   `NULL`.
#'
#' @param st The structural table, or `NULL`. A structural term has no
#'   design columns, so its block is built from the parameters that the term
#'   reports and not from a block of the design, and its hyperparameter is
#'   reported there instead of in a block of its own that would carry nothing
#'   else.
#' @param at_edge The `boundary_key` of [statmod_certificate()]:
#'   the hyperparameter coordinates that have run to an edge of their range
#'   and whose own gradient is already met. Such a coordinate reports no
#'   standard error and no interval, like every quantity that depends on it;
#'   the variance matrix around it is computed and reported as it always was.
#'
#' @keywords internal
summary_blocks <- function(fit, spec, design, p, ci, level = 0.95,
                           V = NULL, st = NULL, at_edge = character(0),
                           test = wald_test()) {
  rows <- ci[ci$parameter == p, , drop = FALSE]
  cols <- c("name", "estimate", "se", "statistic", "p_value", "lower",
            "upper", "role")
  # `source` says what put the number there -- a criterion by name, or the
  # caller -- which `role` alone cannot, every criterion answering "estimated"
  all_cols <- c(cols, "source")
  empty <- stats::setNames(
    data.frame(character(0), numeric(0), numeric(0), numeric(0), numeric(0),
               numeric(0), numeric(0), character(0), character(0),
               stringsAsFactors = FALSE), all_cols)

  coef_rows <- function(nm) {
    r <- rows[rows$term == nm, , drop = FALSE]
    if (!nrow(r)) return(empty)
    out <- data.frame(name = r$coefficient, estimate = r$estimate, se = r$se,
                      statistic = r$statistic, p_value = r$p_value,
                      lower = r$lower, upper = r$upper, role = "coefficient",
                      source = "", stringsAsFactors = FALSE)
    stats::setNames(out, all_cols)
  }
  # WHICH hyperparameter was estimated and by what. A marginal criterion --
  # reml(), ml() -- maximizes a twice differentiable function of it, so it
  # carries a standard error and an interval, both read on the free scale its
  # link defines and mapped back (statmod_hyper_vcov). A path -- aic(), bic(),
  # cv() over a kinked penalty -- chooses the argument of a minimum over a
  # grid, so there is no curvature to read and the row reports the value and
  # the criterion that chose it. Only a value the caller set is held.
  outer_ran <- !is.null(fit@methods$outer)
  Vh <- if (outer_ran) tryCatch(
    statmod_hyper_vcov(spec, design, fit@coefficients, fit@hyper,
                       fit@methods$outer, inner = fit@methods$smooth,
                       pinned = fit@methods$pinned, at_edge = at_edge),
    error = function(e) NULL) else NULL
  spc <- fit@methods$sparse_criterion
  spc_keys <- fit@methods$sparse_hyper
  if (is.null(spc_keys)) spc_keys <- character(0)
  # WHAT WAS HELD is the terms' answer and nobody else's: a hyperparameter
  # a term fixed is fixed whatever criterion ran beside it.
  held <- tryCatch(statmod_held(spec, design), error = function(e) character(0))
  # WHERE A TERM'S NUMBERS ARE when they are not under it. A member of a
  # shared covariance declares no penalty of its own, so its block would
  # otherwise report nothing at all and say so.
  class_note <- function(nm) {
    cl <- class_unit_of(spec, statmod_design(spec), p, nm)
    if (is.null(cl) || length(cl$pieces) < 2L) return(character(0))
    # a label written in a SUBFORMULA puts part of a term's block into the
    # class and leaves the rest where it was, so the two cases are told
    # apart rather than described by one sentence that fits neither
    mine <- Filter(function(z) identical(z$param, p) && identical(z$term, nm),
                   cl$pieces)
    whole <- length(mine) == 1L && is.null(mine[[1L]]$within)
    stats::setNames(
      sprintf("%s in the covariance block '%s', reported above",
              if (whole) "its effects are" else
                "some of its coefficients are", cl$key),
      cl$key)
  }
  # A term may carry more than one penalty, each filed under a key of its
  # own, so the rows of a term are those of every key belonging to it. Where
  # there are several the hyperparameter is named for the penalty as well:
  # two lambdas in one block, one on the slope changes and one on the jumps,
  # are not the same number and cannot appear under the same name.
  hyper_parts <- function(nm, at = p) {
    # `at` is the equation in which the term is kept: a joint filter keeps its
    # hyperparameters there, whichever equation's block is being built
    ent <- modelterms7::term_penalties(term_of_equation(spec, p, nm))
    des <- statmod_design(spec)
    # A LABELLED term declares no penalty of its own: the covariance class
    # carries it. The class's hyperparameters are reported ONCE, under its
    # first member, and a note names the equations the block spans -- a
    # class is one penalty and printing it under each member would say
    # there are several.
    cl <- class_unit_of(spec, des, p, nm)
    if (!is.null(cl)) {
      first <- cl$pieces[[1L]]
      # under the FIRST member only, so one penalty prints once
      if (identical(first$param, p) && identical(first$term, nm)) {
        # ADDED to the term's own entries, never in place of them: a term may
        # carry a penalty of its own beside a labelled sub-term, and replacing
        # them would drop it from the page
        ent <- c(ent, list(list(name = "", penalty = cl$penalty,
                                index = unit_positions(cl),
                                class_key = cl$key)))
      }
    }
    if (!length(ent)) {
      return(list(rows = list(), index = list(), penalty = list(),
                  classes = list()))
    }
    own <- length(ent) - as.integer(!is.null(cl) &&
                                    identical(cl$pieces[[1L]]$param, p) &&
                                    identical(cl$pieces[[1L]]$term, nm))
    out <- lapply(ent, function(e) {
      key <- if (is.null(e$class_key)) {
        statmod_entry_key(nm, utils::head(ent, own), e)
      } else e$class_key
      th <- fit@hyper[[at]][[key]]
      if (is.null(th) || !length(th)) return(empty)
      u <- statmod_unit(spec, des, at, key)
      if (is.null(u)) return(empty)
      # the entry's name is prefixed where the term carries several, so two
      # sigmas in one table stay apart; the prefix is recorded so that a
      # compartment, whose header already names the coefficient, can take it
      # off again
      pre <- if (length(ent) > 1L && nzchar(e$name)) paste0(e$name, ".") else ""
      lab <- paste0(pre, names(th))
      marginal <- outer_ran && !penalty_has_kink(u$penalty)
      hk <- paste(at, key, names(th), sep = "\r")
      # PER HYPERPARAMETER, and the test has to be the vector one: `ifelse`
      # returns a result the length of its TEST, so a scalar `marginal`
      # returned one answer and recycled it over a penalty carrying two
      # hyperparameters -- the elastic net's alpha, which no path varies,
      # was reported as chosen by the criterion that had chosen its lambda
      chosen <- (rep(marginal, length(hk)) | hk %in% spc_keys) &
        !(hk %in% held)
      role <- ifelse(chosen, "estimated", "fixed")
      by <- if (marginal) fit@methods$outer@kind else
        if (is.null(spc)) "estimated" else spc@kind
      # not chosen and not held: with no criterion the value is the term's
      # starting value, which the caller did not give
      src <- ifelse(chosen, by,
                    ifelse(hk %in% held | names(th) %in% names(u$fixed) |
                             outer_ran, "fixed", "start"))
      r <- data.frame(name = lab, estimate = as.numeric(th),
                      se = NA_real_, statistic = NA_real_, p_value = NA_real_,
                      lower = NA_real_, upper = NA_real_, role = role,
                      source = src, stringsAsFactors = FALSE)
      # the interval is built where the criterion was maximized and mapped
      # back, as every other interval in the toolkit is, so a positive
      # hyperparameter cannot be given a negative lower end; the standard
      # error printed beside it is the delta method onto its own scale. No
      # test accompanies it: the null a z of value/se reports on is that the
      # hyperparameter is zero, which for a smoothing parameter is the edge
      # of its range and not an interior hypothesis.
      if (!is.null(Vh) && marginal) {
        z <- stats::qnorm(1 - (1 - level) / 2)
        lk <- attr(Vh, "idx")
        for (i in seq_along(th)) {
          k <- paste(at, key, names(th)[[i]], sep = "\r")
          j <- match(k, rownames(Vh))
          # a shared hyperparameter has one row, under its first member, and
          # the member table says which row every other member reads
          if (is.na(j)) {
            mem <- attr(lk, "members")
            jm <- if (is.null(mem)) NA_integer_ else
              match(k, paste(mem$parameter, mem$term, mem$name, sep = "\r"))
            if (!is.na(jm)) j <- mem$row[[jm]]
          }
          if (is.na(j)) next
          link <- attr(lk, "links")[[j]]
          # a coordinate the criterion's curvature cannot produce a
          # variance for is NA here and keeps the row's empty se and
          # interval; sort() drops an NA, so reaching the ends below with
          # one would ask for an element that is not there
          if (!is.finite(Vh[j, j])) next
          # AND ONE AT A BOUNDARY REPORTS NONE EITHER, for a different reason
          # and with the variance matrix around it untouched. Its curvature is
          # a proper maximum, so a number does come out; what has gone is the
          # quantity's identifiability, and the delta method reads the
          # collapse the wrong way round -- the derivative onto the reported
          # scale tends to zero exactly where the estimate stops being
          # pinned down, so the standard error tends to zero with it. Measured
          # on a covariance class at the boundary, the correlation printed
          # 1.0000 with a standard error of 2.6e-05 beside an interval of
          # (-1, 1): the interval was right and the standard error claimed
          # five decimals of a quantity the data does not identify.
          if (k %in% at_edge || rownames(Vh)[[j]] %in% at_edge) next
          se_eta <- sqrt(Vh[j, j])
          eta <- linkfunctions7::linkfun(link, r$estimate[[i]])
          r$se[[i]] <- abs(linkfunctions7::dlinkinv(link, eta)) * se_eta
          ends <- sort(c(linkfunctions7::linkinv(link, eta - z * se_eta),
                         linkfunctions7::linkinv(link, eta + z * se_eta)))
          r$lower[[i]] <- ends[[1L]]
          r$upper[[i]] <- ends[[2L]]
        }
      }
      # WHERE THE HYPERPARAMETERS ARE A CHART and not the quantities, the
      # quantities are reported instead: a correlated random effect is about
      # the standard deviations and correlations of its effects, and the
      # logarithms of a Cholesky diagonal are the coordinates that produce
      # them. The penalty declares which case it is; every other branch
      # answers NULL and its rows stand.
      rd <- tryCatch(penalties7::penalty_readable(u$penalty, th),
                     error = function(e) NULL)
      if (!is.null(rd) && length(rd$value)) {
        # WHICH EFFECT EACH COORDINATE IS. The chart's coordinates are
        # numbered by the multivariate family, which is a law on R^d and has
        # no model to name them after; here they are the within-group columns
        # of one term, or of several terms across several equations where a
        # label collects them. The labels are passed only where there are as
        # many of them as the prior has dimensions, so a penalty of another
        # shape keeps the names it gives itself.
        own <- if (is.null(e$class_key))
          entry_owner(spec@terms[[p]][[nm]], e$index)
        cvl <- if (!is.null(e$class_key)) class_coords(cl)$label else
          term_coord_labels(own, if (identical(own, spec@terms[[p]][[nm]]))
            e$index)
        blk <- tryCatch(as.integer(u$penalty@block),
                        error = function(z) NA_integer_)
        if (length(blk) != 1L || is.na(blk) || blk != length(cvl)) {
          cvl <- character(0)
        }
        rr <- readable_hyper_rows(rd, th, if (marginal) Vh else NULL, p, key,
                                  level, role, src, cols, cvl, at_edge)
        # A hyperparameter the readable block does not DESCRIBE keeps its own
        # row: a multivariate Student t is about the standard deviations and
        # the correlations of its scale matrix, and its degrees of freedom are
        # none of those. The question is asked of the Jacobian -- a column that
        # is zero throughout is a coordinate no quantity depends on -- so a
        # family that declares more later is covered without an edit.
        # ⚠️ Read against each row's own scale rather than exactly: the same
        # cancellation `jacobian_depends()` records would otherwise take a
        # coordinate's own row away for an entry of relative size 1e-23.
        keep <- !apply(jacobian_depends(rd$jacobian), 2L, any)
        if (any(keep)) {
          rr <- rbind(rr, stats::setNames(r[keep, , drop = FALSE],
                                          c(cols, "source")))
        }
        attr(rr, "prefix") <- pre
        return(rr)
      }
      r <- stats::setNames(r, c(cols, "source"))
      attr(r, "prefix") <- pre
      r
    })
    # PER ENTRY and not stacked, because an entry belongs where its
    # coefficients do: a penalty covering a developed parameter's columns is
    # that parameter's hyperparameter and is reported inside its
    # compartment, while one covering the term's own columns stays with the
    # term. Which is which is read off the columns rather than off the
    # entry's name, so a term that names its entries differently is covered
    # without an edit.
    # A SHARED BLOCK IS NOT THIS TERM'S, and it is held apart here rather
    # than printed under whichever member the walk reached first. Its
    # coordinates belong to several terms and, where a label crosses a bar,
    # to several equations: printed inside one member's block a reader is
    # told the standard deviation of an effect on `sigma` under the term of
    # `mu`, and the other member reports nothing at all. The summary collects
    # these and prints them once, ahead of the equations, where the
    # coordinates can be named.
    shared <- vapply(ent, function(e) !is.null(e$class_key) &&
                       length(cl$pieces) > 1L, TRUE)
    cls <- out[shared]
    names(cls) <- vapply(ent[shared], function(e) e$class_key, "")
    # THE PENALTIES TRAVEL WITH THEIR ROWS rather than being looked up again
    # from the term: the entries here are the term's own plus, where a label
    # collects part of it, the class's, which term_penalties() does not report
    # -- a labelled sub-term declares none, the block being the class's. A
    # second walk would come back one entry shorter and pair each row with the
    # wrong penalty, or with none.
    list(rows = out[!shared], index = lapply(ent[!shared], `[[`, "index"),
         penalty = lapply(ent[!shared], `[[`, "penalty"),
         classes = cls)
  }
  hp_rows <- function(hp) {
    if (!length(hp$rows)) return(empty)
    do.call(rbind, hp$rows)
  }
  # A term whose block is a working linearization is reported by what it is
  # ABOUT and not by the coefficients it is fitted through: a break-point
  # term's auxiliary pair carries the position as -g/delta, which is a
  # number no reader wants to compute. term_readable() gives the quantities
  # and the Jacobian from the coefficients, so the variance comes across by
  # the delta method -- the same route segmented reports a break-point's
  # standard error by. A term that answers NULL, which is what a developed
  # one does, falls back to its coefficients.
  readable_rows <- function(nm) {
    term <- spec@terms[[p]][[nm]]
    idx <- design[[p]]$blocks[[nm]]
    rd <- tryCatch(modelterms7::term_readable(term, fit@coefficients[[p]][idx]),
                   error = function(e) NULL)
    if (is.null(rd) || !length(rd$name)) return(NULL)
    # V is indexed by NAME rather than by position: coef_labels() skips a
    # parameter with no coefficients, so a stacked offset computed here
    # would not be the one the matrix was built with
    key <- paste(p, design[[p]]$coef_names[idx], sep = ":")
    se <- rep(NA_real_, length(rd$name))
    if (!is.null(V) && all(key %in% rownames(V))) {
      Vb <- as.matrix(V[key, key, drop = FALSE])
      # a quantity is read only off the coordinates its row of the Jacobian
      # touches, so a coordinate with no variance -- the slot of a held
      # break-point, a column of zeros -- leaves every other quantity its
      # standard error and takes away only those that depend on it, as
      # readable_vcov() does
      ok <- is.finite(diag(Vb))
      bad <- apply(rd$jacobian[, !ok, drop = FALSE] != 0, 1L, any)
      Vz <- Vb
      Vz[!is.finite(Vz)] <- 0
      se <- sqrt(pmax(diag(rd$jacobian %*% Vz %*% t(rd$jacobian)), 0))
      se[bad] <- NA_real_
    }
    z <- stats::qnorm(1 - (1 - level) / 2)
    st <- rd$value / se
    # A break-point gets an estimate and an interval and NO test. The null
    # a z of value/se would report on is that the position is zero, which
    # is not a hypothesis anyone holds; the one a reader wants is that
    # there is no break-point at all, and under it the position is a
    # nuisance parameter that vanishes, so the classical p-value is wrong
    # by a factor of three to five (Davies' problem). segmented prints the
    # estimate and the standard error alone for the same reason.
    pv <- 2 * stats::pnorm(-abs(st))
    # under a restricted test the statistic of a quantity that IS one of the
    # term's coefficients is that coefficient's restricted statistic; a
    # quantity that mixes coefficients has no restricted statistic, and the
    # Wald value beside the header of another test would be a third thing
    if (test_restricted(test)) {
      cr <- coef_rows(nm)
      for (i in seq_along(st)) {
        nz <- which(rd$jacobian[i, ] != 0)
        one <- length(nz) == 1L && isTRUE(rd$jacobian[i, nz] == 1) &&
          nrow(cr) == ncol(rd$jacobian)
        st[[i]] <- if (one) cr$statistic[[nz]] else NA_real_
        pv[[i]] <- if (one) cr$p_value[[nz]] else NA_real_
      }
    }
    pos <- grepl("^psi[0-9]*([.]|$)", rd$name)
    st[pos] <- NA_real_
    pv[pos] <- NA_real_
    out <- data.frame(name = rd$name, estimate = rd$value, se = se,
                      statistic = st, p_value = pv,
                      lower = rd$value - z * se, upper = rd$value + z * se,
                      role = "coefficient", source = "",
                      stringsAsFactors = FALSE)
    stats::setNames(out, all_cols)
  }
  term_edf <- function(nm) {
    if (is.null(fit@edf)) return(NA_real_)
    e <- fit@edf
    v <- e$edf[e$parameter == p & e$term == nm]
    if (length(v)) v[1L] else NA_real_
  }

  # WHAT A SUB-TERM IS CALLED where a compartment header names it. The
  # term's own label is used where it has one; a parametric development is
  # named for what it carries, an intercept alone being the population
  # value of the parameter and anything else a set of covariates.
  sub_label <- function(s) {
    lb <- tryCatch(s@label, error = function(e) "")
    if (length(lb) == 1L && nzchar(lb)) return(lb)
    nms <- tryCatch(modelterms7::term_coef_names(s),
                    error = function(e) character(0))
    if (length(nms) == 1L && endsWith(nms, "(Intercept)")) "intercept" else
      "covariates"
  }
  # WHAT A HYPERPARAMETER IS, rather than which coordinate carries it: a
  # quadratic penalty's `lambda` is a precision. A random effect's
  # hyperparameters keep the names they carry in an ordinary random() block
  # (sigma, sd[...], cor[...], nu), so a development reads as the same object
  # wherever it is written. Any other name is left as it stands.
  hyper_label <- function(nm, pen) {
    if (identical(nm, "lambda") &&
        isTRUE(tryCatch(penalties7::is_quadratic(pen),
                        error = function(e) FALSE))) {
      return("precision")
    }
    nm
  }
  # WHAT A SUB-TERM CONTRIBUTES to its compartment: its rows filtered the way
  # a block of that kind is filtered at the top level, so a smooth shows its
  # unpenalized part, a kinked penalty shows what it kept, and a random
  # effect shows no coefficients at all: its predictions are summarized by
  # effect_spread(), as in an ordinary random() block, and are still in
  # `coef()`.
  #
  # `rows_at` is how the rows of a set of the term's own coefficients are
  # obtained. An additive term reads its design block's rows by column; a
  # structural one has no design and reads what it REPORTS by position in
  # its parameter vector. Everything below is common to the two.
  sub_block <- function(s, ix, rows_at) {
    sk <- term_block_kind(s)
    rr <- rows_at(ix)
    nms <- tryCatch(modelterms7::term_coef_names(s),
                    error = function(e) character(0))
    # the sub-term's own names, keeping the mark of a coordinate the fit
    # holds (structural_rows() appends it)
    if (length(nms) == nrow(rr)) {
      rr$name <- ifelse(endsWith(rr$name, " (held)"), paste0(nms, " (held)"),
                        nms)
    }
    # AN ALIASED COEFFICIENT IS NOT A SELECTED ONE. Its estimate is `NA`,
    # so `!= 0` is `NA`, and a logical `NA` indexes a row of nothing but
    # `NA`s: the block's label column then carried one, `max(nchar(labs))`
    # was `NA`, and the printer died inside `formatC()` on "missing value
    # where TRUE/FALSE needed" with the header line already mangled.
    # Reachable from any fit whose selected columns are linearly dependent,
    # which is what `nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g))` is.
    keep <- switch(sk,
      smooth = smooth_linear_cols(s, nrow(rr)),
      selection = !is.na(rr$estimate) & rr$estimate != 0,
      rep(TRUE, nrow(rr)))
    if (identical(sk, "random")) keep <- rep(FALSE, nrow(rr))
    # the predictions of a random effect are summarized as an ordinary
    # random() block summarizes them, see effect_spread()
    sp <- if (identical(sk, "random") && nrow(rr)) effect_spread(s, rr$estimate)
    list(kind = sk, table = rr[keep, , drop = FALSE], lines = character(0),
         spread = sp, n_coef = nrow(rr))
  }
  # ONE COMPARTMENT PER DEVELOPED PARAMETER, carrying its own hyperparameter
  # first and then each sub-term's rows in the order the block binds them.
  compartments <- function(term, rows_at, dev, hp) {
    lapply(dev, function(cp) {
      mine <- vapply(hp$index, function(ii) {
        length(ii) > 0L && all(ii %in% cp$index)
      }, logical(1))
      hr <- hp$rows[mine]
      pens <- hp$penalty[mine]
      # INSIDE A COMPARTMENT THE ENTRY'S PREFIX IS TAKEN OFF: the header
      # names the coefficient, so `psi1::random(~1 | id).sigma` under
      # `psi1 ~ intercept + random` reads as `effect sd`, as it does when the
      # term carries one random effect. Where two entries of one compartment
      # would then share a name the prefixes stay, since they are what tells
      # the rows apart.
      bare <- lapply(hr, function(h) {
        pre <- attr(h, "prefix")
        if (!nrow(h) || is.null(pre) || !nzchar(pre)) return(h$name)
        ifelse(startsWith(h$name, pre), substring(h$name, nchar(pre) + 1L),
               h$name)
      })
      if (!anyDuplicated(unlist(bare))) {
        for (i in seq_along(hr)) if (nrow(hr[[i]])) hr[[i]]$name <- bare[[i]]
      }
      for (i in seq_along(hr)) {
        if (!nrow(hr[[i]])) next
        hr[[i]]$name <- vapply(hr[[i]]$name, hyper_label, character(1),
                               pen = pens[[i]], USE.NAMES = FALSE)
      }
      parts <- lapply(seq_along(cp$subs), function(i) {
        sub_block(cp$subs[[i]], cp$sub_index[[i]], rows_at)
      })
      tb <- do.call(rbind, c(hr, lapply(parts, function(z) z$table)))
      ln <- unlist(lapply(parts, function(z) z$lines))
      # PRINTED AS A MIXED MODEL IS: the fixed effects of the development
      # first, then each penalized sub-term under its own name with its
      # hyperparameters, its kept rows and its summary line. A hyperparameter
      # is filed with the sub-term whose coordinates it covers.
      owner <- vapply(seq_along(hr), function(h) {
        ii <- hp$index[mine][[h]]
        w <- which(vapply(cp$sub_index, function(si) {
          length(ii) > 0L && all(ii %in% si)
        }, logical(1)))
        if (length(w)) w[[1L]] else NA_integer_
      }, integer(1))
      pen <- vapply(seq_along(cp$subs), function(i) {
        !identical(parts[[i]]$kind, "parametric") || any(owner %in% i)
      }, logical(1))
      fixed <- do.call(rbind, lapply(parts[!pen], function(z) z$table))
      loose <- do.call(rbind, hr[is.na(owner)])
      sections <- list()
      if (!is.null(fixed) && nrow(fixed)) {
        sections[[1L]] <- list(header = if (any(pen)) "fixed effects",
                               table = fixed, lines = character(0))
      }
      for (i in which(pen)) {
        nm_i <- names(cp$subs)[i]
        if (is.null(nm_i) || !nzchar(nm_i)) nm_i <- sub_label(cp$subs[[i]])
        t_i <- do.call(rbind, c(hr[which(owner == i)], list(parts[[i]]$table)))
        # a random sub-term is headed as an ordinary random() block is, by
        # its own call and the number of its coefficients
        sections[[length(sections) + 1L]] <- list(
          header = if (identical(parts[[i]]$kind, "random"))
            sprintf("%s   [%d coefficients]", nm_i, parts[[i]]$n_coef)
            else nm_i,
          table = if (is.null(t_i)) empty else t_i,
          lines = parts[[i]]$lines, spread = parts[[i]]$spread)
      }
      if (!is.null(loose) && nrow(loose)) {
        sections <- c(list(list(header = NULL, table = loose,
                                lines = character(0))), sections)
      }
      list(name = cp$name,
           header = sprintf("%s ~ %s", cp$name, dev_formula(term, cp)),
           table = if (is.null(tb)) empty else tb,
           lines = if (is.null(ln)) character(0) else ln,
           sections = sections,
           spread = Filter(Negate(is.null),
                           lapply(parts, function(z) z$spread)),
           n_coef = length(cp$index))
    })
  }
  # THE FORMULA THAT DEVELOPS A PARAMETER, as it was written: the term
  # keeps its subformulas by parameter, so `psi ~ 0 + group` is reported as
  # `psi1 ~ 0 + group`. A term that keeps none is described by its
  # sub-terms.
  dev_formula <- function(term, cp) {
    # nl() and the break-point terms keep them as `subformulas`, gas() as
    # `submodels`
    sf <- tryCatch(term@subformulas, error = function(e) NULL)
    if (!is.list(sf) || is.null(sf[[cp$name]])) {
      sf <- tryCatch(term@submodels, error = function(e) NULL)
    }
    f <- if (is.list(sf)) sf[[cp$name]] else NULL
    if (inherits(f, "formula")) {
      return(paste(deparse(f[[length(f)]], width.cutoff = 500L),
                   collapse = " "))
    }
    paste(vapply(cp$subs, sub_label, character(1)), collapse = " + ")
  }
  # THE DEVELOPED PARAMETERS READ AT A GLANCE, one line each: the population
  # value of the development, and what develops it. A parameter that is a
  # number of its own is one row of the table below and needs no line above
  # it; a developed one is spread over a compartment where its population
  # value is labeled by the development's intercept, so this is the only
  # place the parameter's own name appears beside a number.
  head_rows <- function(rows_at, dev, term) {
    out <- lapply(dev, function(cp) {
      r <- NULL
      note <- paste("~", dev_formula(term, cp))
      for (i in seq_along(cp$subs)) {
        j <- which(tryCatch(modelterms7::term_coef_names(cp$subs[[i]]),
                            error = function(e) character(0)) ==
                   "(Intercept)")
        if (length(j) == 1L) {
          r <- rows_at(cp$sub_index[[i]][[j]])
          break
        }
      }
      if (!is.null(r) && !nrow(r)) r <- NULL
      data.frame(name = cp$name,
                 estimate = if (is.null(r)) NA_real_ else r$estimate[[1L]],
                 lower = if (is.null(r)) NA_real_ else r$lower[[1L]],
                 upper = if (is.null(r)) NA_real_ else r$upper[[1L]],
                 note = note, stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
  }
  # WHICH COLUMNS ARE THE TERM'S OWN, and which hyperparameters go with
  # them: everything a developed parameter does not claim.
  own_hyper <- function(hp, dev) {
    hp$rows[!vapply(hp$index, function(ii) {
      any(vapply(dev, function(cp) length(ii) > 0L && all(ii %in% cp$index),
                 logical(1)))
    }, logical(1))]
  }
  developed <- function(term, n_own) {
    cp <- tryCatch(modelterms7::term_components(term),
                   error = function(e) list())
    if (!length(cp)) return(list())
    idx <- unlist(lapply(cp, function(z) z$index), use.names = FALSE)
    if (!length(idx) || max(idx) > n_own) return(list())
    cp
  }
  # A STRUCTURAL TERM contributes no design columns, so its block is built
  # from what it REPORTS -- the quantities term_readable() names, with the
  # variance of the joint information behind them -- and its hyperparameter
  # sits with the coordinates it shrinks rather than in a block of its own
  # carrying nothing else.
  structural_rows <- function(r) {
    # a component with no reported quantity at these positions gives no
    # rows, which the callers read as "nothing to show"
    if (!nrow(r)) return(empty)
    stats::setNames(data.frame(
      name = ifelse(r$held, paste0(r$name, " (held)"), r$name),
      estimate = r$estimate, se = r$se, statistic = NA_real_,
      p_value = NA_real_, lower = r$lower, upper = r$upper,
      role = "coefficient", source = "", stringsAsFactors = FALSE), all_cols)
  }
  structural_block <- function(nm, term, prefix = NULL, at = p) {
    r <- st[st$parameter == p & st$term == nm, , drop = FALSE]
    hp <- hyper_parts(nm, at = at)
    cp_all <- developed(term, length(modelterms7::term_params(term)))
    dev <- cp_all[vapply(cp_all, function(z) length(z$subs) > 0L, logical(1))]
    # a joint filter is reported under each of its equations: this block
    # shows the components of this equation's term, and the hyperparameters
    # of the other equations' components are shown under those
    other <- list()
    shown <- term
    if (!is.null(prefix)) {
      mine <- startsWith(names(dev), prefix)
      other <- dev[!mine]
      dev <- dev[mine]
      # inside the equation's block a component is named and described as
      # the term of that equation names it: kappa1, with its formula
      if (is_joint_filter(term)) {
        shown <- term@terms[[substring(prefix, 1L, nchar(prefix) - 1L)]]
        dev <- lapply(dev, function(cp) {
          cp$name <- substring(cp$name, nchar(prefix) + 1L)
          cp
        })
      }
    }
    rows_at <- function(ix) {
      structural_rows(r[!is.na(r$position) & r$position %in% ix, ,
                        drop = FALSE])
    }
    own <- structural_rows(r[!(r$component %in% c(names(dev), names(other))),
                             , drop = FALSE])
    tb <- do.call(rbind, c(own_hyper(hp, c(dev, other)), list(own)))
    list(kind = "structural", label = block_label("structural"), term = nm,
         n_coef = nrow(r), edf = term_edf(nm), n_zero = 0L,
         table = if (is.null(tb)) empty else tb,
         head = NULL,
         components = compartments(shown, rows_at, dev, hp),
         classes = hp$classes, note = class_note(nm))
  }

  blocks <- list()
  para <- character(0)
  for (nm in names(spec@terms[[p]])) {
    if (identical(term_block_kind(spec@terms[[p]][[nm]]), "parametric")) {
      para <- c(para, nm)
    }
  }
  if (length(para)) {
    tb <- do.call(rbind, lapply(para, coef_rows))
    blocks[[length(blocks) + 1L]] <- list(
      kind = "parametric", label = "Parametric terms", term = NA_character_,
      n_coef = nrow(tb), edf = sum(vapply(para, term_edf, numeric(1))),
      n_zero = 0L, table = tb, head = NULL, components = list())
  }

  for (nm in names(spec@terms[[p]])) {
    term <- spec@terms[[p]][[nm]]
    kind <- term_block_kind(term)
    if (identical(kind, "parametric")) next
    # a joint filter is reported under each of its equations below
    if (is_joint_filter(term)) next
    k <- length(design[[p]]$blocks[[nm]])
    if (!k) {
      # no design columns: a structural term, reported by what it says it
      # reports. Where the structural table could not be built there is
      # nothing to show and the term is passed over rather than given an
      # empty block.
      if (is.null(st) || !any(st$parameter == p & st$term == nm)) next
      blocks[[length(blocks) + 1L]] <- structural_block(nm, term)
      next
    }
    cr <- coef_rows(nm)
    hp <- hyper_parts(nm)
    rows_at <- function(ix) cr[ix, , drop = FALSE]
    # a component whose parameter is a number of its own needs no
    # compartment: it is one coefficient and belongs in the term's own table
    cp_all <- if (nrow(cr) == k) developed(term, k) else list()
    dev <- cp_all[vapply(cp_all, function(z) length(z$subs) > 0L, logical(1))]
    if (length(dev)) {
      own <- setdiff(seq_len(k),
                     unlist(lapply(dev, function(z) z$index),
                            use.names = FALSE))
      # the own table reports the same quantities a term without a
      # development reports: the position of a break-point and not the pair
      # of working coefficients it is read off
      rr <- if (identical(kind, "breakpoint")) readable_rows(nm) else NULL
      # a held break-point reports its developed position in that table, one
      # value per column of the sub-design, so the compartment of the slot
      # that once carried it is not shown
      if (!is.null(rr) &&
          isTRUE(tryCatch(term@blueprint$held, error = function(e) FALSE))) {
        dev <- dev[!grepl("^psi[0-9]+$",
                          vapply(dev, function(z) z$name, character(1)))]
      }
      tb <- do.call(rbind, c(own_hyper(hp, dev), list(
        if (is.null(rr)) cr[own, , drop = FALSE] else rr)))
      # the count of a selection reads the coordinates under a kinked
      # sub-term alone, as a lasso written in an equation does: the term's
      # own parameters and an unpenalized intercept are not selected
      kz <- unlist(lapply(dev, function(z) unlist(lapply(seq_along(z$subs),
        function(s) {
          pen <- tryCatch(z$subs[[s]]@penalty, error = function(e) NULL)
          if (!is.null(pen) && isTRUE(penalty_has_kink(pen))) z$sub_index[[s]]
          else integer(0)
        }))))
      blocks[[length(blocks) + 1L]] <- list(
        kind = kind, label = block_label(kind), term = nm, n_coef = k,
        n_kinked = if (length(kz)) length(kz) else NULL,
        edf = term_edf(nm),
        n_zero = if (identical(kind, "selection") && length(kz))
          sum(cr$estimate[kz] == 0, na.rm = TRUE) else 0L,
        n_na = if (identical(kind, "selection") && length(kz))
          sum(is.na(cr$estimate[kz])) else 0L,
        table = if (is.null(tb)) empty else tb,
        head = NULL,
        components = compartments(term, rows_at, dev, hp),
        classes = hp$classes, note = class_note(nm))
      next
    }
    keep <- switch(kind,
      smooth = smooth_linear_cols(term, nrow(cr)),
      # a coefficient a kinked penalty set to zero is counted, not listed
      # an aliased one is listed, its estimate missing, as base R lists it
      selection = is.na(cr$estimate) | cr$estimate != 0,
      rep(TRUE, nrow(cr)))
    if (identical(kind, "random")) keep <- rep(FALSE, nrow(cr))
    body <- if (identical(kind, "breakpoint")) {
      rr <- readable_rows(nm)
      if (is.null(rr)) cr else rr
    } else {
      cr[keep, , drop = FALSE]
    }
    # the hyperparameters come FIRST in every penalized block: they govern
    # everything below them, and a table that opens with a hundred selected
    # coefficients buries the one number that produced that selection
    tb <- rbind(hp_rows(hp), body)
    blocks[[length(blocks) + 1L]] <- list(
      kind = kind, label = block_label(kind), term = nm, n_coef = k,
      edf = term_edf(nm),
      n_zero = if (identical(kind, "selection"))
        sum(cr$estimate == 0, na.rm = TRUE) else 0L,
      n_na = if (identical(kind, "selection")) sum(is.na(cr$estimate)) else 0L,
      table = tb, head = NULL, components = list(), classes = hp$classes,
      note = class_note(nm),
      spread = if (identical(kind, "random")) effect_spread(term, cr$estimate))
  }
  # a joint filter spans several equations and is reported under each, by
  # the rows of the structural table that belong to this one, under the
  # name its term had in the formula and with the parameters it spends here
  for (jf in joint_filters_of(spec, p)) {
    if (is.null(st) || !any(st$parameter == p & st$term == jf$key)) next
    b <- structural_block(jf$key, jf$term, prefix = paste0(p, "."),
                          at = jf$param)
    r <- st[st$parameter == p & st$term == jf$key, , drop = FALSE]
    if (length(jf$term@keys) && !is.null(jf$term@keys[p]) &&
        !is.na(jf$term@keys[p])) {
      b$term <- unname(jf$term@keys[p])
    }
    b$edf <- sum(!r$held)
    blocks[[length(blocks) + 1L]] <- b
  }
  blocks
}


#' The Covariance Class to Which One Term Belongs
#'
#' @description
#' Returns the penalized unit of the covariance class that covers a term's
#' coefficients, or `NULL` where the term has no label.
#'
#' @details
#' A term with a label declares no penalty of its own, because the class
#' carries it, so the hyperparameters of such a term are read from the class.
#' A class's pieces record the pair of the equation and the term's name, and
#' the lookup uses that pair.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param param The distribution parameter in which the term sits.
#' @param nm The term's name in the formula.
#'
#' @return One unit, as [statmod_penalized()] returns it, or `NULL`.
#'
#' @seealso [statmod_classes()] for the classes, [summary_blocks()] for the
#'   reader that needs this.
#'
#' @keywords internal
class_unit_of <- function(spec, design, param, nm) {
  for (u in statmod_penalized(spec, design)) {
    if (is.null(u$pieces)) next
    hit <- vapply(u$pieces, function(z)
      identical(z$param, param) && identical(z$term, nm), TRUE)
    if (any(hit)) return(u)
  }
  NULL
}


#' The Coordinates of a Covariance Block
#'
#' @description
#' Returns one row per coordinate of a class's covariance, naming the equation
#' to which it belongs, the term that carries it, the within-group column that
#' it is, and the label that a summary prints for it.
#'
#' @details
#' The prior of a covariance class describes the effects of one group over
#' every column that the label collects, so its dimension is the sum of the
#' members' within-group widths. Its coordinates run over the members in
#' order, and each member contributes its own columns in order, which is
#' exactly how [class_index()] interleaves them.
#'
#' The multivariate family that carries the chart gives those coordinates the
#' names `v1`, `v2`, and so on, since it is a law on \eqn{\mathbb{R}^d} without
#' reference to the model. The equation and the term to which a coordinate
#' belongs are determined by this layer; without them a printed correlation
#' between `v1` and `v3` cannot be attributed to any effect.
#'
#' The label is the shortest one that separates the coordinates: only the
#' column where every coordinate comes from one term, the equation and the
#' column where more than one term is involved, and the term as well where two
#' terms of one equation write the same column name. If two labels are still
#' equal, a running number `#i` is appended. A coordinate reached
#' through a subformula carries the parameter that it develops in front of its
#' column, so that an effect on a break-point is not read as an effect on the
#' equation in which the break-point sits.
#'
#' @param u One penalized unit carrying `pieces`, as
#'   [statmod_penalized()] returns for a covariance class.
#'
#' @return A data frame with `param`, `term`, `column` and `label`, one row
#'   per coordinate.
#'
#' @seealso [class_notes()], [summary_class_blocks()], which print it;
#'   [modelterms7::term_group()], which supplies the column names.
#'
#' @keywords internal
class_coords <- function(u) {
  cd <- do.call(rbind, lapply(u$pieces, function(pc) {
    wn <- pc$group$names
    if (length(wn) != pc$dim) wn <- paste0("column", seq_len(pc$dim))
    # WHAT THE EFFECT IS AN EFFECT ON. A label written in a subformula puts
    # the random intercept of a break-point into the block, and a coordinate
    # called `mu:(Intercept)` there would be read as the mean's own.
    if (length(pc$path)) {
      wn <- paste0(paste(pc$path, collapse = ":"), ":", wn)
    }
    data.frame(param = rep(pc$param, pc$dim), term = rep(pc$term, pc$dim),
               column = as.character(wn), stringsAsFactors = FALSE)
  }))
  cd$label <- coord_labels_of(cd)
  cd
}

#' @rdname class_coords
#' @param cd The first three columns of the result.
#' @keywords internal
coord_labels_of <- function(cd) {
  one <- length(unique(paste(cd$param, cd$term, sep = "\r"))) == 1L
  lab <- if (one) cd$column else paste(cd$param, cd$column, sep = ":")
  if (anyDuplicated(lab)) {
    lab <- paste(cd$param, cd$term, cd$column, sep = ":")
  }
  # two columns of one term cannot share a name, so the third form separates
  # everything a design can produce; the numbering is the answer of last
  # resort and exists so that a label is never ambiguous
  if (anyDuplicated(lab)) lab <- paste0(lab, "#", seq_along(lab))
  lab
}

#' The Sub-Term to Which a Penalty Entry Belongs
#'
#' @description
#' Returns the term whose own columns an entry of
#' [modelterms7::term_penalties()] covers: the term itself where the
#' entry is over its own block, and the sub-term developing one of its
#' parameters where it is not.
#'
#' @details
#' A term that develops a parameter over another term reports that term's
#' penalties as its own, and each entry names the columns that it covers.
#' Those columns are the sub-term's block, so a question about the penalty
#' (above all, what its coordinates are) is a question about the sub-term, and
#' reading it from the parent returns nothing: `nl()` has no grouping, and a
#' correlated random effect inside it does.
#'
#' The walk recurses, translating the entry's positions into the sub-term's
#' own on the way down, so a development two levels deep is reached.
#'
#' @param term One built term.
#' @param ii The entry's columns, in the term's own block.
#'
#' @return One built term.
#'
#' @seealso [term_coord_labels()], its caller.
#'
#' @keywords internal
entry_owner <- function(term, ii) {
  if (is.null(ii) || !length(ii)) return(term)
  cp <- tryCatch(modelterms7::term_components(term),
                 error = function(e) list())
  for (c1 in cp) {
    for (k in seq_along(c1$subs)) {
      si <- c1$sub_index[[k]]
      if (length(si) && all(ii %in% si)) {
        return(entry_owner(c1$subs[[k]], match(ii, si)))
      }
    }
  }
  term
}

#' The Coordinates of a Term's Own Covariance
#'
#' @description
#' Returns the within-group column names of a random-effect term, which are
#' the coordinates of its own multivariate prior, or an empty vector for a
#' term that has no grouping.
#'
#' @details
#' Where the prior covers some of the within-group columns and not all of
#' them, as the null-space prior of a smooth random effect does, `ii` gives
#' the entry's positions in the block and the labels are those of the columns
#' that it covers, in the order of the block.
#'
#' @param term One built term.
#' @param ii The entry's positions in the term's block, or `NULL` for every
#'   within-group column.
#'
#' @return A character vector, possibly empty.
#'
#' @seealso [class_coords()] for the same question about a shared block,
#'   [entry_owner()] for the term to which this applies.
#'
#' @keywords internal
term_coord_labels <- function(term, ii = NULL) {
  gr <- tryCatch(modelterms7::term_group(term), error = function(e) NULL)
  if (is.null(gr) || length(gr$names) != gr$dim) return(character(0))
  lab <- as.character(gr$names)
  if (length(ii)) {
    lab <- lab[unique((as.integer(ii) - 1L) %% gr$dim + 1L)]
  }
  lab
}

#' A Coordinate's Name in Place of Its Number
#'
#' @description
#' Rewrites the names that a multivariate family gives to the quantities of a
#' matrix parameter (`sd_v1`, `cor_v1_v2`) into names for the same quantities
#' in terms of the coordinates in which this model is written:
#' `sd[mu:(Intercept)]`, `cor[mu:(Intercept), sigma:(Intercept)]`.
#'
#' @details
#' A multivariate family names its coordinates by position, and each of its
#' readings -- a standard deviation, a correlation, a partial correlation, a
#' conditional scale -- is named `prefix_vi` or `prefix_vi_vj`. The rewrite
#' applies only to that shape: a name that does not match, or whose index has
#' no label, is left unchanged, so a family declaring a reading of another
#' kind keeps its own name.
#'
#' @param nm The names as the family gives them.
#' @param labels One label per coordinate, in the matrix's own order.
#'
#' @return A character vector as long as `nm`.
#'
#' @seealso [class_coords()] and [term_coord_labels()] for the labels,
#'   [distributions7::mv_derived()] for the names that are rewritten.
#'
#' @keywords internal
readable_coord_names <- function(nm, labels) {
  if (!length(labels) || !length(nm)) return(nm)
  ok <- function(i) !anyNA(i) && all(i >= 1L & i <= length(labels))
  out <- nm
  for (k in seq_along(nm)) {
    z <- regmatches(nm[[k]],
                    regexec("^(.+?)_v([0-9]+)_v([0-9]+)$", nm[[k]]))[[1L]]
    if (length(z) == 4L) {
      i <- suppressWarnings(as.integer(z[3:4]))
      if (ok(i)) {
        out[[k]] <- sprintf("%s[%s, %s]", z[[2L]], labels[[i[[1L]]]],
                            labels[[i[[2L]]]])
        next
      }
    }
    z <- regmatches(nm[[k]], regexec("^(.+?)_v([0-9]+)$", nm[[k]]))[[1L]]
    if (length(z) == 3L) {
      i <- suppressWarnings(as.integer(z[[3L]]))
      if (ok(i)) out[[k]] <- sprintf("%s[%s]", z[[2L]], labels[[i]])
    }
  }
  out
}

#' The Legend of a Covariance Block
#'
#' @description
#' Returns one line per coordinate, giving its label and the term from which
#' it comes, for the head of the block that a summary prints.
#'
#' @param cd The coordinates, as [class_coords()] returns them.
#'
#' @return A character vector.
#'
#' @keywords internal
coord_legend <- function(cd) {
  w <- max(nchar(cd$label))
  c("coordinates:",
    sprintf("  %-*s  from %s", w, cd$label, cd$term))
}


#' The Notes on Covariance Classes
#'
#' @description
#' Returns one note per class with two or more pieces (terms, or developments
#' of one term), naming the label, the grouping and the terms whose
#' coefficients share the block.
#'
#' @details
#' A class's hyperparameters are printed once, at the head of the summary,
#' where each coordinate is named by the equation and the column to which it
#' belongs. The note gives the same information in one sentence, for a reader
#' who has the summary object and not the printed page.
#'
#' A class of one member has no note: there is nothing shared to report, and
#' its block is the random effect that it would have been without a label.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#'
#' @return A character vector, possibly empty.
#'
#' @seealso [summary.StatmodFit()], which collects it;
#'   [summary_class_blocks()], which prints the numbers.
#'
#' @keywords internal
class_notes <- function(spec, design) {
  out <- character(0)
  for (u in tryCatch(statmod_penalized(spec, design),
                     error = function(e) list())) {
    if (is.null(u$pieces) || length(u$pieces) < 2L) next
    who <- vapply(u$pieces, function(z) sprintf("'%s' in '%s'", z$term, z$param),
                  "")
    out <- c(out, sprintf(paste0(
      "The covariance block '%s' is shared by %s: one %d-variate prior over ",
      "the effects of one level of %s. It belongs to no single one of them ",
      "and is reported once, ahead of the equations, with each coordinate ",
      "named by its equation and its column."),
      u$key, paste(who, collapse = " and "), u$class$dim, u$class$group))
  }
  out
}


#' The Covariance Blocks That a Summary Prints Ahead of the Equations
#'
#' @description
#' Returns one block record per covariance class shared by more than one
#' member: the standard deviations and correlations of the shared prior, in
#' which each coordinate is named for the equation, the term and the column to
#' which it belongs.
#'
#' @details
#' A shared block belongs to no single member. Printed inside the block of one
#' member, it would report the standard deviation of an effect on `sigma`
#' under a term of `mu`, and the block of the other member would state that
#' there is nothing to report on its own. The block is therefore printed once,
#' before the equations, where its coordinates can be given names.
#'
#' The rows are those that [summary_blocks()] held apart, so the numbers are
#' produced in one place; this function adds the heading of the class and the
#' legend that gives the term of each coordinate.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param tables The per-parameter block lists.
#' @param edf The per-term degrees of freedom, or `NULL`.
#' @param coef The coefficients, or `NULL`. With `hyper`, they let the block
#'   report what the class itself spends rather than what its members' terms
#'   do.
#' @param hyper The hyperparameters, or `NULL`.
#'
#' @return A list of block records in the shape that [print_block()] reads.
#'   Each record carries its coordinates in `coords`. The list is empty where
#'   every class has fewer than two members.
#'
#' @seealso [class_coords()] for the naming, [class_notes()] for the note that
#'   gives the same information in prose.
#'
#' @keywords internal
summary_class_blocks <- function(spec, design, tables, edf = NULL,
                                 coef = NULL, hyper = NULL) {
  # the model's smoother, built at most once and only where a class asks for
  # it: a summary of a model carrying none must not pay for a factorization
  smoother <- NULL
  built <- FALSE
  get_smoother <- function() {
    if (built) return(smoother)
    built <<- TRUE
    if (is.null(coef) || is.null(hyper)) return(NULL)
    smoother <<- tryCatch({
      js <- joint_smoother_diag(spec, coef, design, hyper)
      if (!is.null(js) && !isTRUE(js$failed)) c(js$beta, js$zeta) else {
        H <- statmod_information_at(spec, coef, design, TRUE,
                                  NULL)
        S <- zap_nonfinite(statmod_penalty_at(spec, coef, hyper, design,
                                              "hessian"))
        diag(solve_pd(as_dense(H + S), "the penalized information") %*%
               as_dense(H))
      }
    }, error = function(e) NULL)
    smoother
  }
  rows <- list()
  for (bl in tables) {
    for (b in bl) {
      if (is.null(b$classes)) next
      for (k in names(b$classes)) rows[[k]] <- b$classes[[k]]
    }
  }
  if (!length(rows)) return(list())
  out <- list()
  for (u in statmod_penalized(spec, design)) {
    if (is.null(u$pieces) || length(u$pieces) < 2L) next
    tb <- rows[[u$key]]
    if (is.null(tb) || !nrow(tb)) next
    cd <- class_coords(u)
    # the members' counts add up to the class's: each is the trace of the
    # model's smoother over that member's own columns, and the class is
    # their union.
    #
    # WHAT THE CLASS ITSELF SPENDS, which is the trace of the model's own
    # smoother over the coordinates the class collects and not over whatever
    # else its members' terms carry. Where every member is an ordinary term
    # whose block the class covers entirely the two are the same number --
    # measured, 23.07079 + 17.37986 against 40.45065 -- and where a member is a
    # filter they are not: the term's row counts its level and its persistence
    # as well, which the class does not collect.
    #
    # ONE COUNT PER TERM in the fallback, not one per piece. Two members of a
    # class may be two developments of ONE term -- the loading and the level of
    # one filter, and measured, two labelled subformulas of one nl() would be
    # the same shape -- and the count is filed by (parameter, term), so a sum
    # over the pieces reads that term's row twice: on a filter carrying both it
    # reported 44.00 for a class inside a model whose whole effective count is
    # 24.00.
    e <- edf
    who <- unique(vapply(u$pieces, function(z)
      paste(z$param, z$term, sep = "\r"), ""))
    ed <- sum(vapply(strsplit(who, "\r", fixed = TRUE), function(z) {
      v <- if (is.null(e)) numeric(0) else
        e$edf[e$parameter == z[[1L]] & e$term == z[[2L]]]
      if (length(v)) v[[1L]] else NA_real_
    }, 0))
    sm <- get_smoother()
    if (!is.null(sm)) {
      pos <- unit_joint_positions(u, spec, design)
      if (length(pos) && all(pos <= length(sm))) ed <- sum(sm[pos])
    }
    out[[length(out) + 1L]] <- list(
      kind = "class", label = u$key, term = u$key,
      n_coef = length(unit_positions(u)), edf = NA_real_, n_zero = 0L,
      bits = c(sprintf("%d levels of %s", u$class$m, u$class$group),
               sprintf("%d coordinates per level", u$class$dim),
               if (is.finite(ed)) sprintf("edf %.2f", ed)),
      note = coord_legend(cd), coords = cd,
      table = tb, head = NULL, components = list(), classes = list())
  }
  out
}


#' The Heading Under Which a Block of Each Kind Is Printed
#'
#' @description
#' Returns the name that a summary gives a block, from the kind that
#' [term_block_kind()] returns.
#'
#' @param kind The kind of the block.
#'
#' @return A single string.
#'
#' @keywords internal
block_label <- function(kind) {
  switch(kind, smooth = "Smooth", random = "Random effect",
         selection = "Selection", breakpoint = "Break-points",
         parametric = "Parametric terms", structural = "Structural",
         "Penalized")
}


#' The Notes That a Smoothed Break-Point Term Adds to a Summary
#'
#' @description
#' Returns a note for each smoothed break-point term that names the smoother
#' and the width that was resolved at build, and a second note for each random
#' development of a break-point under a Gaussian prior, when the smoother
#' declares a scale correction, that gives the corrected scale beside the
#' apparent one.
#'
#' @details
#' The correction is declared by the smoother. The probit satisfies the exact
#' convolution identity \eqn{\tau^2_{\mathrm{apparent}} = \tau^2 + h^2}, so
#' the corrected scale is \eqn{\sqrt{\tau^2 - h^2}}. A smoother that declares
#' no correction (the hyperbolic and the quintic) receives only the first
#' note. The apparent scale is the standard deviation parameter of a Gaussian
#' prior written by its scale, or \eqn{1/\sqrt{\lambda}} for a quadratic
#' penalty with the single hyperparameter \eqn{\lambda}. A development with
#' any other prior receives no scale note.
#'
#' @param spec The fitted specification, whose terms are those that the fit
#'   left.
#' @param object The fit.
#'
#' @return A character vector, possibly empty.
#'
#' @seealso [numericals7::abs_smoother()]
#'
#' @keywords internal
smoothed_notes <- function(spec, object) {
  out <- character(0)
  for (p in names(spec@terms)) {
    for (nm in names(spec@terms[[p]])) {
      tm <- spec@terms[[p]][[nm]]
      if (!S7::S7_inherits(tm, modelterms7::SegTerm)) next
      bp <- tryCatch(tm@blueprint, error = function(e) NULL)
      smx <- if (is.list(bp)) bp$smooth else NULL
      if (is.null(smx)) next
      sm <- smx$sm
      out <- c(out, sprintf(paste0(
        "'%s' in '%s' is smoothed (%s, %s = %s%s): its break-points are ",
        "ordinary\n  parameters, the width of the transition is part of the ",
        "model, and its rows are those of the\n  smoothed model."),
        nm, p, sm@smoother_name, sm@width_name,
        format(smx$width, digits = 3),
        if (!is.null(smx$w_group)) ", per group" else ""))
      if (is.null(sm@tau_correction)) next
      ent <- tryCatch(modelterms7::term_penalties(tm),
                      error = function(e) list())
      for (e in ent) {
        if (!grepl("^psi", e$name)) next
        pen <- e$penalty
        key <- statmod_entry_key(nm, ent, e)
        th <- tryCatch(object@hyper[[p]][[key]], error = function(err) NULL)
        if (is.null(th)) next
        # The apparent scale of the break-point deviations. A gaussian
        # prior written by its scale carries it as sigma, which is what
        # random() declares; the quadratic ridge carries the precision, so
        # the scale is 1/sqrt(lambda). Any other prior is left without the
        # note rather than given a number of the wrong meaning: the
        # convolution identity composes GAUSSIAN variances.
        tau <- NULL
        if ("sigma" %in% pen@params &&
            grepl("gaussian", pen@penalty_name, fixed = TRUE)) {
          tau <- suppressWarnings(as.numeric(th[["sigma"]]))
        } else if (identical(pen@params, "lambda") &&
                   isTRUE(tryCatch(penalties7::is_quadratic(pen),
                                   error = function(err) FALSE))) {
          lam <- suppressWarnings(as.numeric(th[["lambda"]]))
          if (length(lam) == 1L && is.finite(lam) && lam > 0) {
            tau <- 1 / sqrt(lam)
          }
        }
        if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) ||
            tau <= 0) next
        tauc <- sm@tau_correction(tau, smx$width)
        out <- c(out, sprintf(paste0(
          "The scale of the random break-points of '%s' composes with the ",
          "smoothing\n  width: apparent tau %s, corrected ",
          "sqrt(tau^2 - %s^2) = %s."),
          nm, format(tau, digits = 4), sm@width_name,
          format(tauc, digits = 4)))
      }
    }
  }
  out
}

#' @title Print a Model Summary
#' @name print.StatmodSummary
#' @description
#' Prints the call, then the covariance blocks shared between equations, then
#' each distribution parameter's blocks, then the degrees of freedom, the
#' criteria and the notes.
#'
#' Below the criteria the point is described by dimensionless numbers from
#' [statmod_certificate()]: for the inner fit, and for the outer criterion
#' where there is one, the largest gradient in conditional standard errors
#' and the smallest eigenvalue of the equilibrated curvature
#' ([certificate_readings()]), and for the coefficients that a kinked penalty
#' has set to zero the largest ratio of the pull of the likelihood to the kink
#' ([zero_readings()]). The foot does not include a verdict. How each is read
#' is on the page of [summary.StatmodFit()], under "Reading the foot". The
#' line of the observations also gives the sum of the prior weights where they
#' are not all 1.
#' @param x A [StatmodSummary()].
#' @param digits Significant digits in the tables.
#' @param notes Whether to print the notes that qualify the numbers. With the
#'   default `FALSE` the foot gives only the number of notes. Most of them
#'   state conventions of the model, and some name the coefficients or the
#'   values that the fit reached, so they can differ between models. The
#'   notes are on the summary's `notes` property in either case.
#' @param max_coef How many coefficient rows a block shows, `Inf` or `NA`
#'   for all of them. `NULL`, the default, uses the value that was given to
#'   [summary.StatmodFit()]; if there is none, it uses the option
#'   `statmodels7.summary_max_coef`, and otherwise 10. A hyperparameter row is
#'   shown whatever `max_coef` is: it governs the coefficients under it, and
#'   every one of them is conditional on the value that it reached. A block
#'   short enough to fit in twelve rows is never abridged, so raising `max_coef`
#'   changes nothing for a parametric block of ordinary size.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @seealso [summary.StatmodFit()]
#' @keywords internal
print.StatmodSummary <- function(x, digits = 4L, notes = FALSE,
                                 max_coef = NULL, ...) {
  # `n` was this argument's name in 0.96.0 and reads as a sample size
  # everywhere else, so it was renamed. It does NOT reach the dots: `n` is
  # a prefix of `notes`, so R matches it there, and the old spelling
  # silently asked for the notes while leaving the row count at its
  # default -- measured, `print(s, n = 4)` sets `notes` to 4. Reading the
  # dots would therefore have caught nothing; the type is what catches it.
  if (!is.logical(notes)) {
    stop("'notes' must be TRUE or FALSE. If you meant how many\n  ",
         "coefficients to print, that argument is 'max_coef', and it is\n  ",
         "on summary() as well: summary(fit, max_coef = Inf).",
         call. = FALSE)
  }
  # what the caller says here wins over what the summary was built with,
  # which in turn wins over the option: printing an object already in hand
  # is the more specific request of the two.
  if (is.null(max_coef)) max_coef <- x@max_coef
  cat("A statmod fit\n\n")
  cat("Call:  ", paste(deparse(x@call), collapse = "\n        "), "\n\n",
      sep = "")
  cat("Distribution: ", x@distrib_name, "     Observations: ", x@n_obs,
      sep = "")
  # the weights are frequencies, so where they are not all 1 the number of
  # observations that nobs() and every criterion read is their sum
  sw <- x@weight_total
  if (length(sw) && !isTRUE(all.equal(sw, x@n_obs))) {
    cat(sprintf("  (prior weights summing to %s, not %d)",
                format(signif(sw, 6)), as.integer(x@n_obs)))
  }
  cat("\n")

  # WHAT THE STATISTIC COLUMN HOLDS is named in the header rather than left
  # to the footer: `z` is the estimate over its standard error and `r` the
  # signed root of a restricted test, and the two are different quantities
  # a reader comparing two summaries would otherwise read as one.
  tst <- x@test@name
  stat <- if (identical(tst, "wald")) "z" else "r"

  # AHEAD OF THE EQUATIONS, because a covariance a label shares belongs to
  # none of them: its coordinates are effects of two equations at once, and
  # under either one of them half of what it says is about the other.
  if (length(x@classes)) {
    cat("\n", strrep("=", 3L), " shared covariance blocks\n", sep = "")
    for (b in x@classes) print_block(b, digits, max_coef, stat)
  }

  for (p in names(x@tables)) {
    lk <- if (p %in% names(x@links)) x@links[[p]] else ""
    cat("\n", strrep("=", 3L), " ", p,
        if (nzchar(lk)) paste0("   [", lk, " link]") else "", "\n", sep = "")
    blocks <- x@tables[[p]]
    if (!length(blocks)) {
      cat("  (no coefficients)\n")
      next
    }
    for (b in blocks) print_block(b, digits, max_coef, stat)
  }

  cat(sprintf("\n%.0f%% intervals, %s variance%s%s\n", 100 * x@level,
              x@vcov_type@name,
              if (is.null(x@information) ||
                  !info_expected(x@information)) "" else
                paste0(" on the ", tolower(info_title(x@information)),
                       info_method_note(x@information)),
              if (identical(tst, "wald")) "" else
                sprintf(";  %s test, signed root", tst)))
  # WHICH log-likelihood and WHICH degrees of freedom, because the pairing
  # is what makes the criterion mean anything and the two conventions in
  # common use are not comparable. This one is conditional: the likelihood
  # is read at the fitted coefficients, a penalized coefficient among them,
  # and the count is the effective degrees of freedom. A mixed-model package
  # reporting a MARGINAL likelihood integrates its random effects out and
  # counts variance components instead, and its AIC is a different number
  # answering a different question (Vaida and Blanchard, 2005).
  cat(sprintf(paste0("conditional log-likelihood %.6f    effective df %.2f",
                     "\ncAIC %.3f    cBIC %.3f\n"),
              x@loglik, x@df, x@aic, x@bic))
  # FOUR NUMBERS AND NO VERDICT (Giovanni, 2026-09-16). The block this
  # replaces printed the search's own flag and the certificate's state with a
  # paragraph of reasons, and both are readings a reader has to interpret --
  # the flag was wrong four times in five where it fired. What is printed now
  # is the point itself, twice over where there is a criterion: how far the
  # gradient is from zero in standard errors, and whether the curvature is
  # that of an optimum, both dimensionless. The certificate still decides
  # what summary() suppresses at a boundary; it is no longer printed.
  cat(sprintf("fitted in %s\n", format_duration(x@elapsed)))
  if (!is.null(x@certificate)) {
    ct <- x@certificate
    line <- function(label, r) {
      if (is.null(r) || !any(is.finite(r))) return(invisible())
      cat(sprintf("%-7s max |grad|/se %s   min eigen %s\n", label,
                  if (is.finite(r[["gradient"]]))
                    format(signif(r[["gradient"]], 2)) else "NA",
                  if (is.finite(r[["eigen"]]))
                    format(signif(r[["eigen"]], 2)) else "NA"))
    }
    line("inner", ct$inner)
    # the coefficients a kinked penalty set to zero have no standard error and
    # no gradient to vanish, so they get a check of their own: at an optimum
    # the likelihood's pull on each is at most the kink, a ratio at most 1
    if (length(ct$zeros) && is.finite(ct$zeros[["ratio"]])) {
      cat(sprintf("zeros   max |score|/kink %s   (%d coefficients at zero)\n",
                  format(signif(ct$zeros[["ratio"]], 2)),
                  as.integer(ct$zeros[["n"]])))
    }
    line("outer", ct$outer)
    # a criterion whose every coordinate ran to an edge, or a coefficient at
    # its link's clamp, leaves nothing to read; a silent absence would read
    # the same as a model with no criterion at all
    if (!any(is.finite(ct$outer)) && length(ct$boundary)) {
      cat("outer   not available: at a boundary\n")
    }
  }
  # THE NOTES ARE COUNTED AND NOT PRINTED. Almost every one of them states a
  # CONVENTION rather than a fact of this fit -- why a hyperparameter's test
  # columns are empty, what a coefficient beside an estimated one is
  # conditional on -- so it reads the same under every model and fills the
  # foot of every summary with a paragraph nobody reads twice. Nothing is
  # dropped in silence: the count says they are there and how to see them,
  # and `summary(fit)@notes` has carried them all along.
  if (length(x@notes)) {
    if (isTRUE(notes)) {
      cat("\n")
      for (n in x@notes) {
        cat(paste0("  ", strwrap(n, width = 74L)), sep = "\n")
        cat("\n")
      }
    } else {
      cat(sprintf("%d note%s: print(summary(fit), notes = TRUE)\n",
                  length(x@notes), if (length(x@notes) > 1L) "s" else ""))
    }
  }
  invisible(x)
}
S7::method(print, StatmodSummary) <- print.StatmodSummary


#' The Rows of a Block Formatted for Printing
#'
#' @description
#' Renders the six numeric columns of a summary table as strings. The cells
#' that have no number for a hyperparameter row are left empty, and the name
#' of each such row carries the origin of its value.
#'
#' @details
#' A hyperparameter row prints numbers where there are any: one estimated by
#' a marginal criterion carries a standard error and an interval. Where there
#' is none the columns are blank. The name of every hyperparameter row carries,
#' in square brackets, the origin of its value (the role `"fixed"` or the
#' criterion that chose the value).
#'
#' The whole block is formatted in one call, its own rows and every
#' compartment's together, so that the widths to which a column is padded are
#' the same throughout and the compartments line up with the table above them.
#'
#' @param tb A summary table.
#' @param digits Significant digits.
#' @param stat The label at the head of the statistic column.
#'
#' @return A list with `cells`, a character matrix of six columns, and
#'   `name`, the row labels.
#'
#' @keywords internal
format_block_cells <- function(tb, digits = 4L, stat = "z") {
  hyp <- tb$role %in% c("fixed", "estimated")
  fixed <- hyp & !is.finite(tb$se)
  # each number on its own: a block mixes rows of very different sizes (a
  # slope of 0.002 beside a position of 450), and one format for the column
  # printed the position as 449.900000 and its standard error as 7.791e+00
  one <- function(v) {
    vapply(v, function(x) {
      if (is.na(x)) "NA" else format(signif(x, digits))
    }, character(1), USE.NAMES = FALSE)
  }
  num <- function(v) ifelse(is.na(v), "", one(v))
  out <- cbind(
    estimate = one(tb$estimate),
    se = num(tb$se),
    z = num(tb$statistic),
    p = ifelse(is.na(tb$p_value), "",
               format.pval(tb$p_value, digits = digits, eps = 1e-16)),
    lower = num(tb$lower),
    upper = num(tb$upper))
  # the header says which statistic the column holds: `z` is the estimate
  # over its standard error and `r` the signed root of a restricted test,
  # and a column headed `z` carrying the second would be a third thing
  colnames(out)[[3L]] <- stat
  out[fixed, c("se", stat, "p", "lower", "upper")] <- ""
  # an estimated one that does carry an interval still has no test: the null
  # a z would report on is that the hyperparameter is zero, the edge of its
  # range rather than an interior hypothesis
  out[hyp & !fixed, c(stat, "p")] <- ""
  src <- if (is.null(tb$source)) tb$role else
    ifelse(nzchar(tb$source), tb$source, tb$role)
  nm <- tb$name
  nm[hyp] <- sprintf("%s [%s]", nm[hyp], src[hyp])
  list(cells = out, name = nm)
}

#' Which Rows of a Table a Summary Prints
#'
#' @description
#' Selects every row of a short table and, of a long one, the
#' hyperparameters together with the first few coefficients. The summary
#' reports the remaining rows as a count.
#'
#' @details
#' A block of a few coefficients is printed whole. A long block is abridged to
#' its first `max_coef` coefficient rows, and the rows that are dropped remain
#' in [coef()]. The hyperparameters are never dropped, whatever the length,
#' because they govern every coefficient under them.
#'
#' @param tb A summary table.
#' @param max_coef How many coefficient rows a long table keeps, `Inf` or `NA`
#'   for all of them. `NULL` reads the option `statmodels7.summary_max_coef`,
#'   and then `show`.
#' @param cap The number of rows up to which a table is printed whole. A table
#'   is cut only above `max(cap, max_coef + number of hyperparameter rows)`.
#' @param show The number of coefficient rows kept when `max_coef` is `NULL`
#'   and the option is not set.
#'
#' @return An integer vector of row positions.
#'
#' @keywords internal
block_rows_shown <- function(tb, max_coef = NULL, cap = 12L, show = 10L) {
  if (is.null(max_coef)) {
    max_coef <- getOption("statmodels7.summary_max_coef", show)
  }
  if (is.na(max_coef) || !is.finite(max_coef)) return(seq_len(nrow(tb)))
  max_coef <- max(0L, as.integer(max_coef))
  # A HYPERPARAMETER ROW IS NEVER HIDDEN. It governs every coefficient under
  # it, and a block that showed ten coefficients and dropped the smoothing
  # parameter that produced them would bury the one number the rest are
  # conditional on.
  hyp <- which(tb$role %in% c("fixed", "estimated"))
  rest <- setdiff(seq_len(nrow(tb)), hyp)
  if (nrow(tb) <= max(cap, max_coef + length(hyp))) {
    return(seq_len(nrow(tb)))
  }
  sort(c(hyp, utils::head(rest, max_coef)))
}

#' Print One Block of a Model Summary
#'
#' @description
#' Prints the heading of a block, the term read at a glance where it is
#' written in parameters of its own, its own coefficients, and one indented
#' compartment per parameter developed over covariates.
#'
#' @details
#' A term that develops one of its own parameters carries columns that mean
#' different things: a break-point's population value and its per-group
#' deviations are not comparable quantities, and a table that stacks them
#' reads as a list of numbers and no longer as a model. Each developed
#' parameter is therefore printed as a compartment of its own, headed by what
#' develops it, with its fixed effects first and then each penalized sub-term,
#' rendered the way a block of that kind is rendered at the top level. A
#' random sub-term is headed by its call and its number of coefficients and
#' reports the hyperparameters of its prior under the names of an ordinary
#' random() block, followed by the summary of its predicted effects
#' ([effect_spread()]); the predictions themselves are in [coef()].
#'
#' @param b A block record from [summary_blocks()].
#' @param digits Significant digits.
#' @param max_coef How many coefficient rows a long block keeps, as
#'   [block_rows_shown()] takes it.
#' @param stat The label at the head of the statistic column: `"z"` for the
#'   estimate over its standard error, `"r"` for the signed root of a
#'   restricted test.
#'
#' @return `NULL`, invisibly. The function is called for the printing.
#'
#' @keywords internal
print_block <- function(b, digits = 4L, max_coef = NULL, stat = "z") {
  head <- if (is.na(b$term)) b$label else b$term
  bits <- character(0)
  if (!is.null(b$bits)) {
    # a block that is not a term of an equation says what it is on its own
    # terms: a covariance class is a number of coordinates over a number of
    # levels, and has neither a coefficient count nor degrees of freedom of
    # its own, those belonging to the members
    bits <- b$bits
  } else if (!identical(b$kind, "parametric")) {
    bits <- c(bits, sprintf("%d %s", b$n_coef,
                            if (identical(b$kind, "structural")) "parameters"
                            else "coefficients"))
    if (is.finite(b$edf)) bits <- c(bits, sprintf("edf %.2f", b$edf))
  }
  if (identical(b$kind, "selection")) {
    nk <- if (is.null(b$n_kinked)) b$n_coef else b$n_kinked
    nn <- if (is.null(b$n_na)) 0L else b$n_na
    bits <- c(bits, sprintf("%d selected, %d at zero", nk - b$n_zero - nn,
                            b$n_zero))
    if (nn > 0L) bits <- c(bits, sprintf("%d not identified", nn))
  }
  cat("\n", head, sep = "")
  if (length(bits)) cat("   [", paste(bits, collapse = ", "), "]", sep = "")
  cat("\n")
  # WHAT THE BLOCK IS PART OF, before its numbers rather than after them: a
  # term whose coefficients are covered by a shared covariance says so, and a
  # shared covariance says which coordinate is which.
  for (z in b$note) cat("  ", z, "\n", sep = "")

  if (!is.null(b$head) && nrow(b$head)) {
    print_block_head(b$head, digits)
    cat("\n")
  }

  comp <- if (is.null(b$components)) list() else b$components
  # the term's own table first, then one section per compartment, each
  # carrying the rows it keeps and, for a random effect, the summary of its
  # predictions
  secs <- list(list(header = NULL, indent = 2L, tb = b$table,
                    lines = character(0), spread = b$spread))
  for (cp in comp) {
    if (is.null(cp$sections)) {
      secs[[length(secs) + 1L]] <- list(header = cp$header, indent = 4L,
                                        tb = cp$table, lines = cp$lines)
      next
    }
    # the formula heads the compartment, and its sections follow without a
    # blank line, indented under their own headers where they have one
    secs[[length(secs) + 1L]] <- list(header = cp$header, indent = 4L,
                                      tb = cp$table[0L, , drop = FALSE],
                                      lines = character(0))
    for (z in cp$sections) {
      secs[[length(secs) + 1L]] <- list(
        header = z$header, indent = if (is.null(z$header)) 4L else 6L,
        tb = z$table, lines = z$lines, gap = FALSE, spread = z$spread)
    }
  }
  for (i in seq_along(secs)) {
    keep <- block_rows_shown(secs[[i]]$tb, max_coef)
    secs[[i]]$hidden <- nrow(secs[[i]]$tb) - length(keep)
    secs[[i]]$tb <- secs[[i]]$tb[keep, , drop = FALSE]
  }
  ns <- vapply(secs, function(z) nrow(z$tb), integer(1))
  if (!sum(ns)) {
    # a member of a shared covariance has said where its numbers are, and
    # "nothing to report" under it would contradict the line above
    if (!length(b$note)) cat("  (nothing to report on its own)\n")
    return(invisible(NULL))
  }
  # THE WHOLE BLOCK IS FORMATTED AT ONCE, so a compartment's numbers are
  # padded to the same widths as the table above it and the columns line up
  # down the block rather than restarting at every section
  fm <- format_block_cells(do.call(rbind, lapply(secs, function(z) z$tb)),
                           digits, stat)
  used <- apply(fm$cells, 2L, function(z) any(nzchar(z)))
  used[[1L]] <- TRUE
  fm$cells <- fm$cells[, used, drop = FALSE]
  # a name inside a compartment is the sub-term's own, and the term's prefix
  # is dropped from its own rows: the heading already says which term this is
  nm <- fm$name
  # the COEFFICIENT rows alone: a hyperparameter's name is not built from
  # the term's, so including it in the question leaves the names no prefix in
  # common and the stripping fires nowhere
  at <- 0L
  for (i in seq_along(secs)) {
    if (ns[[i]] && (i > 1L || !is.na(b$term))) {
      ic <- at + which(secs[[i]]$tb$role == "coefficient")
      if (length(ic)) nm[ic] <- drop_common_prefix(nm[ic])
    }
    at <- at + ns[[i]]
  }
  labs <- character(0)
  ind <- integer(0)
  at <- 0L
  rows <- list()
  for (i in seq_along(secs)) {
    ii <- if (ns[[i]]) at + seq_len(ns[[i]]) else integer(0)
    at <- at + ns[[i]]
    rows[[i]] <- ii
    labs <- c(labs, nm[ii], names(secs[[i]]$lines))
    ind <- c(ind, rep(secs[[i]]$indent, ns[[i]] + length(secs[[i]]$lines)))
  }
  w <- max(nchar(labs) + ind)
  cw <- pmax(nchar(colnames(fm$cells)),
             apply(nchar(fm$cells), 2L, max))
  pad <- function(x, k) {
    mapply(function(v, j) formatC(v, width = j, flag = " "), x, k,
           USE.NAMES = FALSE)
  }
  cat(sub("[ ]+$", "", paste0(strrep(" ", w + 2L),
      paste(pad(colnames(fm$cells), cw), collapse = " "))), "\n", sep = "")
  for (i in seq_along(secs)) {
    s <- secs[[i]]
    if (!is.null(s$header)) cat(if (!isFALSE(s$gap)) "\n",
                                strrep(" ", s$indent - 2L), s$header,
                                "\n", sep = "")
    for (r in rows[[i]]) {
      cat(sub("[ ]+$", "", paste0(
        strrep(" ", s$indent),
        formatC(nm[[r]], width = -(w - s$indent + 2L), flag = " "),
        paste(pad(fm$cells[r, ], cw), collapse = " "))), "\n", sep = "")
    }
    for (k in seq_along(s$lines)) {
      cat(sub("[ ]+$", "", paste0(
        strrep(" ", s$indent),
        formatC(names(s$lines)[[k]], width = -(w - s$indent + 2L),
                flag = " "), s$lines[[k]])), "\n", sep = "")
    }
    if (s$hidden > 0L) {
      cat(strrep(" ", s$indent),
          sprintf("... %d more, in coef()\n", s$hidden), sep = "")
    }
    print_effect_spread(s$spread, indent = s$indent, digits = digits)
  }
  invisible(NULL)
}

#' The Predicted Random Effects of a Term, Summarized
#'
#' @description
#' Computes the six-number summary of [base::summary()] (minimum, quartiles,
#' mean and maximum) of a random-effect term's predicted effects, with one row
#' per coordinate of the effect: the intercept, and each covariate of a random
#' slope.
#'
#' @details
#' The coefficients of one group are adjacent in the block, so the
#' predictions are read as a matrix with one row per level and one column per
#' coordinate, and the coordinates are named as [modelterms7::term_group()]
#' names them. A term with one level and several coordinates is a single draw of
#' a vector, and its table has one row, the summary across the coordinates. The
#' statistics are taken over the finite effects, because an aliased effect is
#' `NA`.
#'
#' @param term A built random-effect term.
#' @param est Its predicted effects, in the order of its coefficients.
#'
#' @return A list with `table`, a data frame with one row per coordinate and
#'   the six summary columns, `levels`, the number of groups, and `group`, the
#'   grouping expression; `NULL` where the term carries no grouping.
#'
#' @keywords internal
effect_spread <- function(term, est) {
  gr <- tryCatch(modelterms7::term_group(term), error = function(e) NULL)
  est <- as.numeric(est)
  if (is.null(gr) || !length(est)) return(NULL)
  m <- length(gr$levels)
  d <- gr$dim
  if (m * d != length(est)) return(NULL)
  M <- matrix(est, nrow = m, ncol = d, byrow = TRUE)
  # an aliased effect is NA, and summary() would add an "NA's" column to its
  # row only; the six statistics are taken over the finite effects
  six <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  row6 <- function(z) {
    z <- z[is.finite(z)]
    v <- if (length(z)) as.numeric(summary(z)) else rep(NA_real_, 6L)
    as.data.frame(as.list(v), col.names = six, check.names = FALSE)
  }
  if (m == 1L && d > 1L) {
    # one group is one draw of a vector (a latent process over the levels of
    # a within-group factor), and the spread is read across its coordinates
    tb <- row6(M[1L, ])
    rownames(tb) <- sprintf("%d coordinates", d)
  } else {
    tb <- do.call(rbind, lapply(seq_len(d), function(j) row6(M[, j])))
    rownames(tb) <- gr$names
  }
  list(table = tb, levels = m,
       group = paste(deparse(gr$expr, width.cutoff = 500L), collapse = " "))
}


#' Print the Summary of a Term's Predicted Random Effects
#'
#' @description
#' Prints what [effect_spread()] returns under the heading
#' `predicted effects (m levels of g)`, one row per coordinate.
#'
#' @param sp What [effect_spread()] returns.
#' @param indent The indentation of the heading.
#' @param digits Significant digits.
#'
#' @return `NULL`, invisibly. The function is called for the printing.
#'
#' @keywords internal
print_effect_spread <- function(sp, indent = 2L, digits = 4L) {
  if (is.null(sp)) return(invisible(NULL))
  cat(strrep(" ", indent),
      sprintf("predicted effects (%d %s of %s)\n", sp$levels,
              if (sp$levels == 1L) "level" else "levels", sp$group),
      sep = "")
  tb <- sp$table
  # each number on its own, a mean that is zero up to rounding printed as 0
  cells <- t(apply(as.matrix(tb), 1L, function(z) {
    z <- zapsmall(z, digits = digits + 3L)
    vapply(z, function(v) format(signif(v, digits)), character(1))
  }))
  cells <- matrix(cells, nrow = nrow(tb))
  cw <- pmax(nchar(colnames(tb)), apply(nchar(cells), 2L, max))
  w <- max(nchar(rownames(tb)))
  pre <- strrep(" ", indent + 2L)
  pad <- function(x) {
    mapply(function(v, j) formatC(v, width = j, flag = " "), x, cw,
           USE.NAMES = FALSE)
  }
  cat(pre, strrep(" ", w), "  ", paste(pad(colnames(tb)), collapse = " "),
      "\n", sep = "")
  for (i in seq_len(nrow(tb))) {
    cat(pre, formatC(rownames(tb)[i], width = -w, flag = " "), "  ",
        paste(pad(cells[i, ]), collapse = " "), "\n", sep = "")
  }
  invisible(NULL)
}


#' The Term Read at a Glance
#'
#' @description
#' Prints one line per parameter in which a term is written, giving its value
#' and interval where it is a number and, where it is developed over
#' covariates, the population value of that development beside what develops
#' it.
#'
#' @param hd The head record of a block.
#' @param digits Significant digits.
#'
#' @return `NULL`, invisibly. The function is called for the printing.
#'
#' @keywords internal
print_block_head <- function(hd, digits = 4L) {
  val <- ifelse(is.na(hd$estimate), "", format(signif(hd$estimate, digits)))
  ci <- ifelse(is.na(hd$lower) | is.na(hd$upper), "",
               sprintf("[%s, %s]", format(signif(hd$lower, digits)),
                       format(signif(hd$upper, digits))))
  w <- max(nchar(hd$name))
  for (i in seq_len(nrow(hd))) {
    cat(sub("[ ]+$", "", paste0(
      "  ", formatC(hd$name[[i]], width = -(w + 2L), flag = " "),
      formatC(val[[i]], width = max(nchar(val)), flag = " "), "  ",
      formatC(ci[[i]], width = -max(nchar(ci)), flag = " "), "  ",
      hd$note[[i]])), "\n", sep = "")
  }
  invisible(NULL)
}

#' Drop the Prefix That a Set of Coefficient Names Share
#'
#' @description
#' Removes the first dotted piece where every name carries the same one and
#' removing it leaves every name non-empty.
#'
#' @details
#' A term composes its coefficients' names from its own and its parameters',
#' so inside the block of one term the leading piece repeats on every row
#' and repeats what the heading has already stated. It is dropped only for
#' printing; [coef()] and the summary's own tables keep the
#' names with which the fit was built, which are the names that another call
#' can use as an index.
#'
#' Only one piece is dropped, and never every piece that the names share: the
#' term's own name is one, and a set of coefficients that agree further along,
#' such as `r.1` and `r.2` of a matrix column, would otherwise be left as bare
#' numbers.
#'
#' @param nms The names.
#'
#' @return The names, shortened where they share a prefix.
#'
#' @keywords internal
drop_common_prefix <- function(nms) {
  if (length(nms) < 1L) return(nms)
  sp <- strsplit(nms, ".", fixed = TRUE)
  if (min(lengths(sp)) < 2L) return(nms)
  if (length(unique(vapply(sp, function(z) z[[1L]], character(1)))) != 1L) {
    return(nms)
  }
  vapply(sp, function(z) paste(z[-1L], collapse = "."), character(1),
         USE.NAMES = FALSE)
}


#' The Certificate of the Reported Point
#'
#' @description
#' Reads the outer criterion at the reported point, independently of the path
#' of the search, and returns a state (`"converged"`, `"boundary"`,
#' `"not converged"` or `"unknown"`) that is determined by the Newton
#' decrement of the criterion. The result also holds the mode error of the
#' coefficients, the hyperparameters and flat directions at a boundary, and the
#' dimensionless readings that [summary.StatmodFit()] prints.
#'
#' @details
#' The flag of an optimizer records whether a search stopped on its own rule,
#' which describes the search and does not order fits by quality. The
#' certificate describes the reported point and does not depend on the path
#' that led to it.
#'
#' The state comes from the rise that the criterion would still buy, and the
#' mode error is reported beside it. That rise is the Newton decrement
#' \eqn{\tfrac{1}{2}g^\top A^{-1}g} of [joint_decrement()], in the criterion's
#' own units, and `tol` is a threshold on it. The mode error is the number of
#' log-likelihood units by which the coefficients lie above the penalized
#' mode. Where an outer criterion exists, the mode error does not enter the
#' verdict; for a model with no penalty and no coefficient estimated on a
#' criterion there is no other reading, and the mode error against
#' [mode_error_limit()] determines the state (see below).
#'
#' The verdict is a threshold on the Newton decrement and not on the gradient,
#' because the gradient of a criterion summed over \eqn{n} observations and
#' \eqn{p} penalized coefficients depends on the model, so that a single
#' threshold cannot serve every model. See [joint_decrement()]. The default
#' `tol` is 1e-2. A user can see and check a report of `"not converged"` at a
#' good point, whereas a user cannot check the certification of a point that
#' is not a maximum, so the threshold is set low.
#'
#' The cost of the certificate is one outer gradient, one curvature and one
#' solve, at a point that the fit already holds. The criterion is reconstructed
#' from `fit@spec`, so the reading refers to the fitted model. Where the form
#' carries an analytic outer Hessian, the certificate does not refit. Where it
#' does not (a criterion on the expected information for a family that does
#' not supply the second derivative of its expected information, a model with
#' a refreshable block or a structural term, or a penalty whose Hessian moves
#' with the coefficients beside a structural term), the curvature is obtained
#' by central differences of the exact gradient over the hyperparameters and
#' the coefficients that the criterion estimates, which takes
#' \eqn{4(n_h + n_\gamma)} refits. [outer_curvature()] records which route was
#' taken, and the result reports it in `curvature`.
#'
#' Where there is no outer gradient to read, the state depends on the case. A
#' model with no penalty and no coefficient that the criterion estimates (a
#' family of one parameter, or `marginal` empty), as with `linpar`, `nl`,
#' `seg`, `jump` and `jseg`, has no coordinate for which an outer gradient
#' could be read, so the only question left is whether the inner fit reached
#' its mode. The state is `"converged"` when the mode error is at most
#' [mode_error_limit()] and `"not converged"` otherwise. In every other case
#' without an outer gradient the state is `"unknown"` and `reason` names the
#' case: a model whose only hyperparameters are kinked (`lasso`, `scad`,
#' `mcp`), which are chosen along a path because a Laplace approximation at a
#' mode that sits on the kink has no meaning, so that the hyperparameter is
#' the argument of a minimum over a grid; hyperparameters that are all held;
#' hyperparameters chosen by [aic()], [bic()] or [cv()]; and a model fitted
#' with no outer criterion. `"unknown"` states that the certificate has no
#' outer gradient to read and is not a report of failure to converge. With a
#' kinked penalty the mode error is reported over the coordinates that the
#' kink leaves free ([free_of_kinks()]): at a coefficient that the penalty has
#' set to zero the score does not vanish but lies in the subdifferential, so
#' that reading the score there would give a false mode error. A model that
#' also carries a smooth penalty estimated by [reml()] or [ml()] has an outer
#' gradient, which the certificate reads as usual.
#'
#' A form whose criterion has no exact gradient ([outer_gradient_ok()] at
#' order 1) has the criterion itself differenced by refits when it has at most
#' two hyperparameters and no coefficient estimated on the criterion (three
#' refits for one coordinate, nine for two), and `curvature` is then
#' `"differenced criterion"`. Otherwise the state is `"unknown"`, and `reason`
#' gives the cause. This is consistent with [outer_curvature()], which
#' differences a Hessian where the analytic one is missing, because that
#' function differences the exact gradient.
#'
#' A hyperparameter may run to an edge and belong there: on a covariate with
#' no effect the smoothing parameter grows without bound, the criterion is
#' flat in it, and the fit has not failed to converge. A coordinate is
#' reported as sitting at a boundary where the criterion no longer moves with
#' it, that is where its own curvature \eqn{\lvert A_{jj}\rvert} is at most
#' `flat`. Such a hyperparameter has no meaningful interval, because
#' \eqn{1/\sqrt{A_{jj}}} is the conditional standard error of its free value,
#' and `boundary_key` gives [summary.StatmodFit()] the keys of these
#' hyperparameters, so that it leaves their standard errors out.
#'
#' The label is read from the curvature and not from the free value, because
#' the size of a smoothing parameter depends on how its penalty is normalized
#' and on the units of the response. The default `flat` of 2e-3 corresponds to
#' a conditional standard error of the free value above 22.
#'
#' The set of coordinates excluded from the verdict is narrower: a coordinate
#' belongs to it exactly when it is reported at a boundary and its
#' [coord_decrement()], the rise that the coordinate would still give if it
#' were the only one to move, is at most `tol`. The second condition prevents
#' `flat` from determining the state: restricting a decrement to a subset of the
#' coordinates is that same maximum under a constraint and so no larger, hence a
#' coordinate moved out of the interior set cannot raise the quantity that
#' determines the state. A coordinate at an edge that has not settled stays
#' under test and is also named. Both `"converged"` and `"boundary"` are
#' certified.
#'
#' The set of reported coordinates and the set excluded from the verdict can
#' differ: at a boundary the curvature collapses along with the gradient, so
#' the one-coordinate reading \eqn{g_j^2/(2A_{jj})} is a ratio of two
#' quantities going to zero and can be large. Such a coordinate is named and
#' stays under test.
#'
#' A coefficient that the criterion estimates beside the hyperparameters
#' (the argument `marginal` of [reml()] and [ml()]) is read the same way, but
#' its curvature is multiplied by \eqn{\max(1, \gamma_j^2)}: it is reported at a
#' boundary where its conditional standard error exceeds \eqn{22} times
#' \eqn{\max(1, \lvert\gamma_j\rvert)}. A coefficient carries the units of
#' its parameter and of its covariate, which a hyperparameter's free value
#' does not, so the bare curvature depends on those units. The ratio of the
#' coefficient to its conditional standard error does not. The degrees of
#' freedom of a Student \eqn{t} at their gaussian limit are an example of such
#' a coordinate.
#'
#' Where no curvature can be read at all (both routes fail to produce one, and
#' the state is then `"unknown"`) the label falls back on the free value,
#' `edge`.
#'
#' A ridge is a direction along which several coordinates move together
#' while the criterion does not change, although the criterion is curved in
#' each of those coordinates separately. It is
#' read on the equilibrated curvature \eqn{A/(ss')},
#' \eqn{s_j = \sqrt{|A_{jj}|}}, whose eigenvalues are dimensionless. An
#' eigenvalue of at most `flat` in absolute value names its eigenvector,
#' written on the original coordinates in `boundary` and in `ridge` (see
#' [certificate_ridges()]). A ridge whose own decrement is within `tol` leaves
#' the verdict and the state is `"boundary"`. A Burr XII shape that runs to its
#' Weibull limit in one level of a factor is an example: the eigenvector moves
#' only that level.
#'
#' A fit that [statmod()] has recorded with a coordinate at the edge of its
#' chart whose return raises the log-likelihood has the state `"not converged"`
#' whatever the readings above give, and `reason` names the coordinate.
#'
#' @param fit A [StatmodFit()].
#' @param tol The largest rise, in the criterion's own units, that a certified
#'   point may still have available; see [joint_decrement()].
#' @param flat The largest curvature \eqn{\lvert A_{jj}\rvert} of the outer
#'   criterion in a hyperparameter's free value at which that hyperparameter
#'   is reported as sitting at a boundary; for a coefficient the criterion
#'   estimates, the curvature times \eqn{\max(1, \gamma_j^2)}, and for a free
#'   value of a matrix prior that scales with the matrix, the curvature in
#'   units of the matrix's scale ([hyper_unit_scale()]). It is also the
#'   largest eigenvalue, in absolute value, of a flat direction of the
#'   equilibrated curvature. It determines only the label and never the
#'   verdict; see the details.
#' @param edge The free value beyond which a hyperparameter is reported as
#'   sitting at a boundary where no curvature can be read, and only there.
#'
#' @return A list with `state` (`"converged"`, `"boundary"`,
#'   `"not converged"` or `"unknown"`), `decrement`, `gradient`,
#'   `mode_error`, `curvature`, `boundary`, `boundary_key`, `ridge` and
#'   `reason`. `ridge` is a list with one named vector per flat direction of
#'   the outer curvature, over the coordinates that move along it.
#'   The verdict uses `decrement`, and `gradient` is reported beside it.
#'   `curvature` is `"analytic"` when the decrement was read against an analytic
#'   Hessian, `"differenced"` against one differenced from the exact gradient,
#'   `"differenced criterion"` when the criterion itself was differenced by
#'   refits, and `NA` when no curvature was read. `boundary_key` names the
#'   hyperparameters among the coordinates of `boundary`, in the key by which
#'   [statmod_hyper_vcov()] labels its rows; [summary.StatmodFit()] reads it to
#'   leave a standard error off a hyperparameter pinned at an edge.
#'   `boundary_coef` names the coefficients estimated on the outer criterion
#'   among the coordinates of `boundary`, as `parameter:name`, and
#'   [summary.StatmodFit()] reads it to leave their standard errors, tests and
#'   intervals off. The state is `"not converged"` also when the fit recorded a
#'   coordinate at the edge of its chart whose return raises the log-likelihood.
#'
#'   Three further fields are those that [summary.StatmodFit()] prints, and its
#'   page describes how they are read: `inner` and `outer`, each a named vector
#'   with `gradient` and `eigen` from [certificate_readings()] (`NA` where
#'   there is nothing to read), and `zeros`, which
#'   reads `NULL` where no coefficient sits at the kink of a penalty and is
#'   otherwise a named vector with `ratio` and `n` from [zero_readings()].
#'
#' @seealso [statmod()], [joint_decrement()], [outer_curvature()],
#'   [mode_error_limit()], [criterion_resolution()]
#'
#' @examples
#' dd <- data.frame(x = runif(120))
#' dd$y <- sin(4 * dd$x) + rnorm(120, 0, 0.3)
#' statmod_certificate(statmod(y ~ s(x, bspline_smooth(k = 8)),
#'                             distributions7::gaussian1_distrib(), dd))$state
#'
#' @export
statmod_certificate <- function(fit, tol = 1e-2, flat = 2e-3, edge = 8) {
  out <- certificate_core(fit, tol, flat, edge)
  # A COORDINATE AT THE EDGE OF ITS CHART THAT POINTS INWARD makes the point
  # something other than a maximum, whatever the readings above say: they
  # are taken on the free scale, where the chart's own derivative vanishes
  # at the edge and takes the score with it. statmod() records what
  # edge_violations() found after its restarts, and a violation that
  # survived them is reported here.
  kk <- fit@methods$kkt
  if (is.data.frame(kk) && nrow(kk)) {
    out$state <- "not converged"
    out$reason <- c(sprintf(paste0(
      "%s sits at the edge of its chart (free value %.4g) and moving it back",
      " to %g raises the log-likelihood by %.4g, so the point is not a",
      " maximum"), ifelse(nzchar(kk$term),
                          paste0(kk$param, ":", kk$term, ":", kk$name),
                          paste0(kk$param, ":", kk$name)),
      kk$eta, kk$target, kk$gain), out$reason)
  }
  out
}

#' The Certificate's Readings Before the Edge Check
#'
#' @description
#' Computes everything that [statmod_certificate()] computes except the check
#' at the edge of a chart, which that function applies to the result.
#'
#' @inheritParams statmod_certificate
#' @return The list that [statmod_certificate()] returns.
#' @keywords internal
certificate_core <- function(fit, tol, flat, edge) {
  out <- list(state = "unknown", decrement = NA_real_, gradient = NA_real_,
              mode_error = NA_real_, curvature = NA_character_,
              boundary = character(0), boundary_key = character(0),
              boundary_coef = character(0),
              ridge = list(), reason = character(0),
              inner = c(gradient = NA_real_, eigen = NA_real_),
              outer = c(gradient = NA_real_, eigen = NA_real_),
              zeros = NULL)
  method <- fit@methods$outer
  spec <- fit@spec
  design <- tryCatch(statmod_design(spec), error = function(e) NULL)
  if (is.null(design)) {
    out$reason <- "the design could not be rebuilt from the fit"
    return(out)
  }
  cf <- fit@coefficients
  hy <- fit@hyper
  # THE COEFFICIENTS THE CRITERION ESTIMATED are outer coordinates: they are
  # held in the inner fit, so their score is not a mode error, and they are
  # read beside the hyperparameters in the outer gradient and curvature. The
  # coordinates the model does not identify stay out of the determinant, as
  # they did during the search.
  gam <- tryCatch(marginal_coords(spec, design, method),
                  error = function(e) marginal_coords(spec, design, NULL))
  ng <- length(gam$where)
  if (ng) {
    beta_all <- unlist(cf[spec@distrib@params], use.names = FALSE)
    spec <- statmod_hold(spec, gam, beta_all[gam$where])
  }
  if (length(fit@methods$pinned)) attr(design, "pinned") <- fit@methods$pinned
  ctx <- tryCatch(outer_context(spec, design, cf, hy,
                       NULL),
                  error = function(e) NULL)

  # HOW FAR ABOVE THE MODE, in log-likelihood units: the decrease the mode's
  # own Newton correction predicts. Reported whether or not a criterion ran.
  if (!is.null(ctx)) {
    pen <- tryCatch(ctx_penalized(ctx, spec, design, cf, hy, FALSE),
                    error = function(e) NULL)
    if (!is.null(pen)) {
      obj <- tryCatch(statmod_objective(spec, hy, design, FALSE,
                                          NULL),
                      error = function(e) NULL)
      sc <- if (is.null(obj)) NULL else
        tryCatch(obj$gr(obj$stack(cf)), error = function(e) NULL)
      if (!is.null(sc) && ng) sc[gam$where] <- 0
      if (!is.null(sc)) sc <- free_of_kinks(sc, spec, design, cf)
      # A BREAK-POINT ON AN OBSERVATION is a linear constraint on the
      # coefficients (kink_constraints()), and the inner fit solves on its
      # null space: the score along a constrained direction is a subgradient
      # and not a residual, so the mode error and the readings are taken over
      # the directions that the constraints leave free, as mode_solver() does
      kc <- if (is.null(sc) || is.null(obj)) NULL else tryCatch({
        b0 <- obj$stack(cf)
        kink_constraints(spec, design, cf, obj$split(seq_along(b0)),
                         length(b0))
      }, error = function(e) NULL)
      if (!is.null(kc) && !nrow(kc)) kc <- NULL
      if (!is.null(sc) && all(is.finite(sc))) {
        db <- tryCatch(if (is.null(kc)) {
          as.numeric(as.matrix(pen$inv) %*% sc)
        } else {
          Nk <- null_basis(kc, length(sc))
          Kd <- as.matrix(pen$K)
          as.numeric(Nk %*% solve(crossprod(Nk, Kd %*% Nk),
                                  crossprod(Nk, sc)))
        }, error = function(e) NULL)
        if (!is.null(db) && all(is.finite(db))) {
          out$mode_error <- 0.5 * sum(sc * db)
        }
      }
      # A STRUCTURAL TERM'S OWN PARAMETERS ARE PART OF THE INNER FIT, so the
      # mode error is read over the joint vector, as alternation_settled()
      # reads it for the convergence flag. Read on the coefficients alone,
      # with the term held, it missed the direction in which the coefficients
      # and the term move together: on MASS::geyser, `waiting ~ regime(k = 2)`
      # stopped 1.3e-3 below the maximum with a joint decrement of 1.275e-3,
      # and the certificate read 2.4e-4 and `converged` while the flag read
      # FALSE. The coefficients' part of the gradient is the one read above,
      # with the held and kinked coordinates set to zero.
      su <- attr(design, "structural")
      if (length(su) && !is.null(obj) && !is.null(sc) &&
          all(is.finite(sc)) && is.finite(out$mode_error)) {
        mj <- tryCatch({
          jp <- statmod_joint_pieces(spec, design, obj, hy,
                                     kinds = c("filter", "loglik"))
          u <- c(obj$stack(cf), jp$zeta())
          if (!is.finite(jp$raw(u))) stop("not evaluable")
          z <- jp$zeta()
          gj <- c(sc, jp$gr(u)[jp$ix])
          Kj <- as_dense(jp$he(u))
          Kj <- (Kj + t(Kj)) / 2
          # a parameter of the term that a kinked penalty holds at zero has a
          # score in the subdifferential, as a kinked coefficient does, so it
          # is left out (measured on gas(by = g, omega ~ 0 + lasso(~ g)): read
          # with those parameters, the mode error was 8.6e-3 at a KKT point
          # where the coefficients read 2.1e-9)
          nm <- names(statmod_structural_state(design)$zeta[[jp$key]])
          zk <- integer(0)
          for (bl in statmod_blocks(spec, design)$sparse) {
            if (!isTRUE(bl$structural) || !identical(bl$zterm, jp$key)) next
            pos <- match(nm[bl$cols], jp$free)
            pos <- pos[!is.na(pos)]
            zk <- c(zk, pos[z[pos] == 0])
          }
          sz <- setdiff(seq_along(jp$ix), zk)
          if (is.null(kc)) {
            keep <- c(seq_along(sc), length(sc) + sz)
            dep <- deficient_coords(Kj[keep, keep, drop = FALSE])
            if (length(dep)) keep <- keep[-dep]
            0.5 * sum(gj[keep] * solve(Kj[keep, keep, drop = FALSE],
                                       gj[keep]))
          } else {
            nz <- length(jp$ix)
            Nk <- null_basis(kc, length(sc))
            Iz <- diag(nz)[, sz, drop = FALSE]
            N <- rbind(cbind(Nk, matrix(0, nrow(Nk), ncol(Iz))),
                       cbind(matrix(0, nz, ncol(Nk)), Iz))
            gn <- crossprod(N, gj)
            0.5 * sum(gn * solve(crossprod(N, Kj %*% N), gn))
          }
        }, error = function(e) NA_real_)
        if (length(mj) == 1L && is.finite(mj)) out$mode_error <- abs(mj)
      }
      # THE TWO NUMBERS summary() PRINTS for the inner fit, over the
      # coordinates vcov() reports: a coefficient a kinked penalty set to zero
      # has a score that lies in the subdifferential rather than vanishing, a
      # frozen working block is not a Jacobian, and an aliased column has no
      # curvature of its own, so each would read as a failure it is not.
      if (!is.null(sc)) {
        lab <- tryCatch(coef_labels(spec, design), error = function(e) NULL)
        keep <- seq_along(sc)
        if (!is.null(lab) && nrow(lab) == length(sc)) {
          beta <- unlist(cf[spec@distrib@params], use.names = FALSE)
          drop <- (lab$kinked & beta == 0) | frozen_block(spec, lab) |
            rownames(lab) %in% fit@aliased
          keep <- which(!drop)
        }
        out$inner <- tryCatch(kink_readings(sc, pen$K, keep, kc),
                              error = function(e) out$inner)
        if (!is.null(lab) && nrow(lab) == length(sc)) {
          out$zeros <- tryCatch(zero_readings(spec, design, obj, cf, hy, lab),
                                error = function(e) NULL)
        }
        # A STRUCTURAL TERM'S OWN PARAMETERS (a filter's or a regime's) ARE
        # ESTIMATED BESIDE THE COEFFICIENTS, in one system, so the inner
        # reading is over that joint vector: read on
        # the coefficients alone, a model whose equation is one intercept and
        # a gas() term reports an equilibrated eigenvalue of exactly 1 whatever
        # the filter's curvature is. The closures are the fit's own.
        if (length(su) && !is.null(obj) &&
            all(vapply(su, function(u) u$kind %in% c("filter", "loglik"),
                       TRUE))) {
          out$inner <- tryCatch({
            jp <- statmod_joint_pieces(spec, design, obj, hy,
                                       kinds = c("filter", "loglik"))
            u <- c(obj$stack(cf), jp$zeta())
            if (!is.finite(jp$raw(u))) stop("not evaluable")
            kink_readings(jp$gr(u), jp$he(u), c(keep, jp$ix), kc)
          }, error = function(e) c(gradient = NA_real_, eigen = NA_real_))
        }
      }
    }
  }

  blocks <- tryCatch(statmod_blocks(spec, design), error = function(e) NULL)
  idx <- if (is.null(blocks)) NULL else outer_hyper_index(spec, blocks)
  no_outer <- is.null(method) || !method@kind %in% c("ml", "reml") ||
    is.null(idx) || (!nrow(idx) && !ng)
  if (no_outer) {
    # ⚠️ TWO QUITE DIFFERENT WAYS TO HAVE NO OUTER GRADIENT, and they call for
    # different answers. Measured over the reference battery, seven cases of
    # twenty-nine land here:
    #
    #   NO PENALTY AT ALL -- linpar, nl, seg, jump, jseg sharp and smoothed.
    #     There is no hyperparameter to estimate, so there is nothing an outer
    #     gradient could say, and the only question left is whether the inner
    #     fit reached its mode. The mode error answers exactly that, and it
    #     answers well: 5.2e-11, 1.5e-10, 1.9e-10, 2.0e-06 and 7.9e-05 on fits
    #     that are right, against 1.215 on a `jump` fitted to data carrying a
    #     slope and a slope change it has no term for -- misspecified, so its
    #     break-point iteration never settles, its annealing runs to the floor
    #     and its block reaches 3.0e+13 on the information's diagonal. On data
    #     where the jump IS the truth the same term reads 4.2e-04.
    #
    #   A KINKED HYPERPARAMETER -- lasso, scad, mcp -- chosen by a PATH
    #     because a Laplace approximation at a mode sitting on the kink is
    #     arithmetic without a meaning. There is a hyperparameter, but it is
    #     the argmin over a grid rather than the root of a derivative. AND THE
    #     MODE ERROR IS NOT A READING HERE EITHER: at a coefficient the
    #     penalty has set to zero the score does not vanish, it lies in the
    #     subdifferential. Measured on a lasso, the mode error of 4.7e-03 is
    #     carried by a coordinate whose coefficient is exactly 0 and whose
    #     score is -0.715. So this stays unknown, and says why.
    pen_units <- tryCatch(statmod_penalized(spec, design),
                          error = function(e) list())
    if (length(pen_units)) {
      # ⚠️ AND CASES THE TWO ABOVE DO NOT COVER, whose reason used to be the
      # kinked one whatever the model carried: a prediction-error criterion on
      # a SMOOTH penalty, and hyperparameters held at their values. Measured on
      # a gamma smooth under aic(), the reason read "the only hyperparameters
      # here belong to a penalty with a kink" of a model carrying none. The
      # state stays unknown -- this certificate reads the outer gradient of
      # reml() and ml() alone -- and the reason says which case it is.
      kinked <- vapply(pen_units, function(u)
        isTRUE(tryCatch(penalty_has_kink(u$penalty, u$key),
                        error = function(e) FALSE)), logical(1))
      estimated <- vapply(pen_units, function(u)
        length(setdiff(u$penalty@params, names(u$fixed))) > 0L, logical(1))
      first <- if (!any(estimated)) {
        paste0("every hyperparameter here was held at the value written on its ",
               "term, so there is no outer gradient to read")
      } else if (all(kinked[estimated])) {
        paste0("the only hyperparameters here belong to a penalty with a kink, ",
               "chosen along a path rather than by a criterion with a derivative")
      } else if (is.null(method)) {
        paste0("the fit ran no outer criterion, so the hyperparameters of the smooth ",
               "penalties stayed at their starting values and there is no outer ",
               "gradient to read")
      } else if (method@kind %in% c("aic", "bic")) {
        sprintf(paste0("the hyperparameters here were chosen by %s(), %s, and ",
                       "this certificate reads the outer gradient of reml() ",
                       "and ml() alone"), method@kind,
                if (identical(method@kind, "aic")) "a prediction-error criterion"
                else "an information criterion")
      } else if (identical(method@kind, "cv")) {
        paste0("the hyperparameters here were chosen by cv(), which refits on ",
               "folds and has no gradient to read")
      } else {
        paste0("the criterion ", method@kind, "() estimates no hyperparameter ",
               "here, so there is no outer gradient to read")
      }
      second <- if (any(kinked) && is.finite(out$mode_error)) {
        sprintf(paste0("the inner fit is %.4g log-likelihood units above its ",
                       "mode over the coordinates that the kink leaves free"),
                out$mode_error)
      } else if (is.finite(out$mode_error)) {
        sprintf("the inner fit is %.4g log-likelihood units above its mode",
                out$mode_error)
      }
      out$reason <- paste(c(first, second), collapse = "; ")
      return(out)
    }
    if (!is.finite(out$mode_error)) {
      out$reason <- "the model carries no penalty, and the mode error could not be read"
      return(out)
    }
    out$state <- if (out$mode_error <= mode_error_limit()) "converged" else
      "not converged"
    out$reason <- sprintf(paste0(
      "the model carries no penalty, so there is no outer gradient; the ",
      "reading is the inner fit's own, %.4g log-likelihood units above its ",
      "mode against %g"), out$mode_error, mode_error_limit())
    return(out)
  }
  # WHERE THE FORM HAS NO EXACT OUTER GRADIENT the criterion itself is
  # differenced, one refit per point, as long as the coordinates are few
  # enough for that to cost less than the fit: three refits for one
  # coordinate and nine for two. Measured on a marginal seg over a 20 x 15
  # panel with one hyperparameter, the certificate read `unknown` beside a
  # fit of 159 s; differenced, it costs three refits warm-started at the fit.
  diffd <- NULL
  if (!outer_gradient_ok(spec, design, idx, method, 1L, gamma = ng > 0L)) {
    nv <- nrow(idx) + ng
    if (ng) {
      out$reason <- sprintf(paste0(
        "this form has no exact outer gradient, and the criterion is ",
        "differenced over\n  hyperparameters alone, while %d coefficient%s ",
        "estimated on it"), ng, if (ng == 1L) " is" else "s are")
      return(out)
    }
    if (nv > 2L) {
      out$reason <- sprintf(paste0(
        "this form has no exact outer gradient, and differencing it over %d ",
        "coordinates would take %d refits"), nv, 1L + 2L * nv^2)
      return(out)
    }
    diffd <- tryCatch(criterion_differenced(fit, spec, design, blocks, method,
                                            idx, cf, hy),
                      error = function(e) NULL)
    if (is.null(diffd)) {
      out$reason <- paste0("this form has no exact outer gradient, and the ",
                           "refits that difference the criterion did not all ",
                           "give a finite value")
      return(out)
    }
  }
  basis <- integrated_basis(spec, design, method@kind, gamma = ng > 0L)
  g <- if (!is.null(diffd)) diffd$g else
    tryCatch(statmod_marginal_grad(spec, design, cf, hy, method, idx, basis,
                                   ctx = ctx, gam = gam),
             error = function(e) NULL)
  if (is.null(g) || !all(is.finite(g))) {
    # ⚠️ A NON-FINITE OUTER GRADIENT AT A BOUNDARY IS NOT THE SAME COMPLAINT,
    # and saying which it is costs one read. Where a coefficient sits at the
    # clamp its link keeps it strictly inside, the family's THIRD and FOURTH
    # derivatives there are not finite -- measured on the Student t, all ten
    # components of the third go from nu = 1e150, which is sqrt(double.xmax)
    # and so the signature of a product of two quantities of order nu -- and
    # the outer gradient reads exactly those. The fit itself may be perfectly
    # well located: on the reference battery's `fam-studentt` the mode error
    # is 1.3e-10 against a limit of 1e-3 and the criterion is the best of the
    # three routes, so what is missing is the reading and not the answer.
    #
    # The state stays `not converged` rather than `boundary`, deliberately:
    # certifying a point whose hyperparameters were never verified would
    # claim more than has been checked. What it can honestly do is name the
    # coordinate, so the reader is not left with "not finite".
    H0 <- tryCatch(ctx_information(ctx, spec, design, cf, hy, FALSE,
                                   NULL),
                   error = function(e) NULL)
    frozen <- if (is.null(H0)) integer(0) else boundary_coords(H0)
    if (length(frozen)) {
      npar <- vapply(design, function(d) d$npar, integer(1))
      ends <- cumsum(npar)
      who <- unique(names(npar)[vapply(frozen, function(j)
        which(j <= ends)[1], integer(1))])
      out$boundary <- who
      out$reason <- sprintf(paste0(
        "the outer gradient is not finite at the reported point, because a ",
        "coefficient of %s sits at the bound inside which its link keeps it ",
        "strictly, and the family's third and fourth derivatives are not finite ",
        "there. The fit itself is %.4g log-likelihood units above its mode"),
        paste(who, collapse = ", "), out$mode_error)
    } else {
      out$reason <- "the outer gradient is not finite at the reported point"
    }
    out$state <- "not converged"
    return(out)
  }
  # A COORDINATE AT A BOUNDARY is one the criterion no longer moves with: its
  # own curvature |A_jj| is at most `flat`. The threshold decides the label
  # and never the verdict, for the reason the conjunction below states: a
  # coordinate is EXCLUDED from the test only if coord_decrement() says it has
  # ALREADY met `tol` on its own, so removing it from `interior` cannot raise
  # what decides the state -- an exact inequality since 0.127.0, see the block
  # below. Delete that second conjunct and the threshold starts excusing
  # coordinates from the test.
  #
  # ⚠️ IT READ THE FREE VALUE UNTIL 0.153.0, |eta| > 8, and a smoothing
  # parameter's size is the penalty's normalization and the response's units
  # as much as the data: the same smooth sits at eta -0.20, 9.01 and -13.96 as
  # the response is multiplied by 1, 0.01 and 1000, its curvature at 2.50,
  # 2.50 and 2.71. What `flat` separates, measured over the reference battery,
  # the lotto 0 and D1 nets and MASS::mcycle: 4.2e-07 to 2.6e-04 at a real
  # edge -- smooths on pure noise, six prior scales of hyper-8, a random
  # slope's scale at eta -6.27 that the value never named -- against 1.5e-02
  # to 168 in the interior, the smallest being a random break-point's prior
  # scale at eta -5.32. 2e-3 is the geometric middle of that gap.
  eta <- hyper_to_eta(hy, idx)
  cv <- if (!is.null(diffd)) {
    list(A = diffd$A, source = "differenced criterion")
  } else {
    outer_curvature(spec, design, cf, hy, method, idx, basis,
                    fit@methods$smooth, gam = gam)
  }
  if (is.null(cv$A)) {
    # ⚠️ NO VERDICT, BUT STILL THE BOUNDARY COORDINATES, and leaving them out
    # was a regression CI found and this machine hid. `boundary_key` is what
    # summary() reads to leave a standard error off a coordinate pinned at an
    # edge, so returning early with none takes that away from a reader on
    # exactly the fits least able to spare it -- a mixed covariance class
    # inside a filter, where the order-2 route is refused and the stencil
    # refuses too.
    #
    # With no curvature there is no `flat` to read, so the label falls back
    # on the free value, `edge` -- a stopgap, and the one place the value
    # still decides anything. The state stays "unknown" and the reason says
    # why, so nothing is certified on a value alone.
    out$gradient <- max(abs(g))
    at_edge <- which(abs(eta) > edge)
    if (length(at_edge)) {
      out$boundary <- vapply(at_edge, function(k)
        paste(idx$parameter[k], idx$term[k], idx$name[k], sep = "/"), character(1))
      out$boundary_key <- vapply(at_edge, function(k)
        paste(idx$parameter[k], idx$term[k], idx$name[k], sep = "\r"), character(1))
    }
    out$reason <- cv$why
    return(out)
  }
  out$curvature <- cv$source
  # ⚠️ TWO SETS AND NOT ONE, and conflating them is what made a coordinate at
  # a chart's edge stop being NAMED as one. What is REPORTED -- and so what
  # summary() reads to leave a standard error off a coordinate pinned there --
  # is the flat curvature alone: a hyperparameter the criterion no longer
  # moves with has no meaningful interval whatever its gradient says.
  #
  # What may be EXCLUDED from the verdict is the narrower set, and it keeps
  # the second conjunct for the reason it has always had: removing a
  # coordinate from the interior can only lower the decrement, so a set
  # chosen by the curvature alone would let `flat` decide the state.
  #
  # Measured, the two genuinely differ, and exactly where the label matters
  # most: at a boundary the curvature collapses with the gradient, so the
  # one-coordinate reading g_j^2/(2 A_jj) is a ratio of two quantities going
  # to zero and can come back large. Such a coordinate is named and stays
  # under test.
  # An estimated coefficient has an edge too: a Student t's nu estimated on
  # the criterion runs to its gaussian limit, where the curvature in it is
  # 3e-13 and the certificate otherwise read nothing at all. It is named by
  # its equation and its coefficient, and only a hyperparameter carries a
  # key, which is what summary() reads for a hyperparameter's row.
  #
  # ⚠️ A COEFFICIENT CARRIES UNITS, which a hyperparameter's free value does
  # not, so its bare curvature is not a statement about an edge: a constant
  # sigma on the identity link with the response in centimeters read 1538
  # with a standard error of 154 and was named a boundary, while the same fit
  # in meters was not. Its curvature is multiplied by max(1, gamma^2), which
  # reads the standard error against the size of the coefficient and is the
  # same in any units.
  #
  # ⚠️ SO DOES A FREE VALUE OF A MATRIX PRIOR THAT SCALES WITH THE MATRIX.
  # The off-diagonal free values of log_cholesky() are entries of a Cholesky
  # factor, in the units of the effects: a correlated random intercept and
  # slope with the response multiplied by 1000 read 1.77e-5 where the same
  # fit read 17.71, and was named a boundary at the same optimum. Such a
  # curvature is read in units of the matrix's own scale (hyper_unit_scale());
  # a log-scale or a correlation coordinate keeps its bare curvature, so the
  # edges measured above are read as before.
  rel_scale <- rep(1, length(g))
  rel_scale[seq_len(nrow(idx))] <- hyper_unit_scale(spec, design, hy, idx)
  if (ng) rel_scale[nrow(idx) + seq_len(ng)] <- pmax(1, beta_all[gam$where]^2)
  reported <- which(abs(diag(cv$A)) * rel_scale <= flat)
  rep_h <- reported[reported <= nrow(idx)]
  rep_g <- reported[reported > nrow(idx)] - nrow(idx)
  if (length(reported)) {
    out$boundary <- c(
      vapply(rep_h, function(k)
        paste(idx$parameter[k], idx$term[k], idx$name[k], sep = "/"),
        character(1)),
      vapply(rep_g, function(k) paste(gam$param[k], gam$name[k], sep = "/"),
             character(1)))
  }
  # the same coefficients by the row names of vcov() and of the summary
  # table, which summary() reads to leave their standard errors off
  if (length(rep_g)) {
    out$boundary_coef <- paste(gam$param[rep_g], gam$name[rep_g], sep = ":")
  }
  if (length(rep_h)) {
    # THE SAME COORDINATES IN THE KEY summary() reads its variance matrix by.
    # It is returned rather than recomposed there because a term's name is a
    # deparsed call and may carry a slash of its own, so parsing `boundary`
    # back into three pieces is not an inverse; and because two callers who
    # each build a key agree by accident, which this file records as a defect
    # of its own.
    out$boundary_key <- vapply(rep_h, function(k)
      paste(idx$parameter[k], idx$term[k], idx$name[k], sep = "\r"),
      character(1))
  }
  at_edge <- intersect(reported, which(coord_decrement(g, cv$A) <= tol))
  interior <- setdiff(seq_along(g), at_edge)
  # A RIDGE IS A DIRECTION, and the reading above takes one coordinate at a
  # time, so it cannot see one along which every coordinate moves. Measured on
  # a Burr XII whose shape k runs to its Weibull limit in the base level of a
  # factor: the flat direction is the intercept of log k with the opposite
  # change in the other two levels, each coordinate alone has a curvature
  # above 24, and the state was `unknown` or `converged` on the sign of an
  # eigenvalue of 2e-7. The equilibrated curvature A / (s s'), s the root of
  # its diagonal, reads such a direction directly: its smallest eigenvalue
  # was +-2e-7 there, with the eigenvector the ridge exactly, against 0.375
  # to 1 over the reference battery, where a coordinate at an edge leaves it
  # near one. A direction whose eigenvalue is at most `flat` in absolute value
  # is named; it leaves the verdict, as an edge does, only if what it would
  # buy alone is already within `tol`.
  rg <- certificate_ridges(g, cv$A, setdiff(seq_along(g), reported), flat, tol)
  if (length(rg$dirs)) {
    lab_all <- c(paste(idx$parameter, idx$term, idx$name, sep = "/"),
                 if (ng) paste(gam$param, gam$name, sep = "/"))
    key_all <- c(paste(idx$parameter, idx$term, idx$name, sep = "\r"),
                 rep(NA_character_, ng))
    where_all <- c(rep(NA_integer_, nrow(idx)), gam$where)
    out$ridge <- lapply(rg$dirs, function(x) {
      on <- which(x != 0)
      structure(stats::setNames(x[on], lab_all[on]), where = where_all[on])
    })
    out$boundary <- c(out$boundary, vapply(out$ridge, function(x)
      paste0("ridge: ", paste0(ifelse(x < 0, "- ", "+ "), format(abs(x),
                                                                 digits = 3),
                               " ", names(x), collapse = " ")), ""))
    hk <- unique(unlist(lapply(rg$dirs, function(x)
      key_all[which(x != 0)])))
    out$boundary_key <- unique(c(out$boundary_key, hk[!is.na(hk)]))
  }
  # ⚠️ WITH NO INTERIOR COORDINATE THERE IS NOTHING TO MEASURE, and 0 was the
  # wrong way to say so: printed as `outer gradient 0` it reads as the
  # gradient VANISHING -- a stationary point, the best news a reader could
  # have -- where it means that every coordinate is at a boundary and none
  # was tested. It is the empty-result shape this file records for a query
  # that returns no rows and for a check that never ran: an absence reported
  # as a pass. NA instead, and summary() then prints no gradient at all, the
  # `at a boundary` line beneath it already naming every coordinate.
  #
  # Measured on a smooth over a response carrying no signal, where REML sends
  # the one smoothing parameter to 7.6e+07 and it is the whole hyperparameter
  # vector: the line read `outer gradient 0` beside `BOUNDARY` and now reads
  # nothing, with the state and the mode error unchanged.
  out$gradient <- if (length(interior)) max(abs(g[interior])) else NA_real_
  # the outer pair summary() prints, over the same coordinates the verdict
  # reads: one that ran to an edge has a curvature collapsing with its
  # gradient, and its equilibrated reading would say nothing about the rest
  out$outer <- certificate_readings(g, cv$A, interior)
  # THE VERDICT IS THE NEWTON DECREMENT over the coordinates still under test,
  # which is how much the criterion would still rise if the point took the
  # step its own curvature calls for. `tol` is therefore in the criterion's
  # own units and means the same thing at every sample size and on every
  # shape, which a threshold on the gradient does not.
  #
  # ⚠️ REMOVING A COORDINATE FROM THE INTERIOR SET CANNOT RAISE IT, so the
  # paragraph that made `edge` safe survives the change of reading, and by an
  # exact inequality rather than by a measurement. The decrement is
  # max_x (2 g'x - x'Ax) over the whole vector; restricting to the interior
  # coordinates is that same maximum under x_E = 0, a constrained maximum of
  # the same function, hence no larger. The per-coordinate reading `at_edge`
  # uses is the same maximum restricted to one coordinate, so it too is
  # bounded by the joint one and a coordinate called an edge had already
  # passed the test the verdict applies.
  # A ridge that has met `tol` on its own is taken out the same way, as the
  # maximum over the directions orthogonal to it, which is again a
  # constrained maximum of the same function and so no larger.
  if (length(interior)) {
    dec <- if (length(rg$settled)) {
      decrement_off(g, cv$A, interior, rg$settled)
    } else {
      joint_decrement(g[interior], cv$A[interior, interior, drop = FALSE])
    }
  } else {
    dec <- 0
  }
  out$decrement <- if (length(interior)) dec else NA_real_
  if (!is.finite(dec)) {
    # the curvature is not negative definite in every direction still under
    # test, so the point is not a maximum of the criterion in those
    # coordinates and no decrement follows from it. Saying that is worth more
    # than a verdict read off a gradient in the wrong units.
    out$reason <- paste0(
      "the outer criterion's curvature at the reported point is not negative ",
      "definite in every coordinate still under test, so the rise that a ",
      "Newton step would give cannot be read")
    return(out)
  }
  if (dec <= tol) {
    out$state <- if (length(at_edge) || length(rg$settled)) "boundary" else
      "converged"
  } else {
    out$state <- "not converged"
    out$reason <- sprintf(paste0(
      "the outer criterion would still rise by %.4g at the reported point ",
      "after the Newton step of its own curvature, which exceeds the limit ",
      "%g; the gradient there is %.4g"), dec, tol, out$gradient)
  }
  if (is.finite(out$mode_error) && out$mode_error > mode_error_limit()) {
    out$reason <- c(out$reason, sprintf(
      paste0("the coefficients sit %.4g log-likelihood units above the",
             " penalized mode, so the criterion, its gradient and vcov()",
             " are all read at that point"),
      out$mode_error))
  }
  out
}



#' The Outer Criterion Differenced by Refitting
#'
#' @description
#' Computes the gradient and the curvature of the outer criterion in its free
#' coordinates by central differences of the criterion itself: the model is
#' refitted at each displaced point, warm-started at the fitted coefficients,
#' and the criterion is read at the refitted mode.
#'
#' @details
#' With one coordinate the readings are
#' \eqn{g = (c_+ - c_-)/(2h)} and \eqn{A = -(c_+ - 2c_0 + c_-)/h^2}; with two
#' the four axial points give the diagonal and the four diagonal points the
#' off-diagonal entry
#' \eqn{A_{12} = -(c_{++} - c_{+-} - c_{-+} + c_{--})/(4h^2)}. The center
#' \eqn{c_0} is refitted too, so every reading comes from the same inner
#' iteration. The step is \eqn{h = 10^{-2}} on the free scale, where the
#' logarithm of a hyperparameter moves by one percent; the truncation error is
#' of order \eqn{h^2} and the effect of an inner tolerance \eqn{\epsilon}
#' on the curvature of order \eqn{\epsilon/h^2}.
#'
#' @param fit The [StatmodFit()].
#' @param spec,design,blocks The specification, its design and its blocks.
#' @param method The outer criterion.
#' @param idx The hyperparameter index, [outer_hyper_index()].
#' @param cf,hy The fitted coefficients and hyperparameters.
#' @param h The step on the free scale.
#'
#' @return A list with the gradient `g` and the curvature `A` (minus the
#'   Hessian). An error is signaled where a refit gives no finite
#'   criterion; the caller reads it as `NULL`.
#'
#' @keywords internal
criterion_differenced <- function(fit, spec, design, blocks, method, idx, cf,
                                  hy, h = 1e-2) {
  io <- fit@methods$smooth
  budget <- method_budget(io)
  cfg <- inner_settings(io)
  emethod <- cfg$emethod
  obj <- statmod_objective(spec, hy, design, cfg$expected, emethod)
  warm <- obj$stack(cf)
  basis <- integrated_basis(spec, design, method@kind)
  eta0 <- hyper_to_eta(hy, idx)
  nv <- length(eta0)
  at <- function(d) {
    hy2 <- eta_to_hyper(eta0 + d, idx, hy)
    res <- statmod_alternate(spec, design, blocks, hy2, io, warm, cfg$expected,
                             emethod, budget$maxit, budget$tol, verbosity(0),
                             hold_refresh = TRUE)
    cf2 <- res$obj$split(res$par)
    ctx2 <- outer_context(spec, design, cf2, hy2, emethod)
    m <- statmod_marginal(spec, design, cf2, hy2, method, emethod, basis, ctx2)
    if (is.null(m) || !is.finite(m$value)) stop("no criterion")
    m$value
  }
  e <- function(j) {
    v <- numeric(nv)
    v[j] <- 1
    v
  }
  c0 <- at(numeric(nv))
  g <- numeric(nv)
  A <- matrix(0, nv, nv)
  for (j in seq_len(nv)) {
    cp <- at(h * e(j))
    cm <- at(-h * e(j))
    g[j] <- (cp - cm) / (2 * h)
    A[j, j] <- -(cp - 2 * c0 + cm) / h^2
  }
  if (nv == 2L) {
    cpp <- at(h * (e(1) + e(2)))
    cpm <- at(h * (e(1) - e(2)))
    cmp <- at(h * (-e(1) + e(2)))
    cmm <- at(-h * (e(1) + e(2)))
    A[1, 2] <- A[2, 1] <- -(cpp - cpm - cmp + cmm) / (4 * h^2)
  }
  list(g = g, A = A)
}


#' The Outer Criterion's Curvature at a Reported Point
#'
#' @description
#' Computes \eqn{A = -(H + H^\top)/2}, the matrix against which a Newton
#' decrement is read, from the analytic Hessian where the form has one and
#' from one central difference of the exact gradient where it does not.
#'
#' @details
#' [statmod_certificate()] reads the rise that the criterion would still buy,
#' which needs a curvature and not only a gradient. The analytic route is
#' [statmod_marginal_hess()], gated by [outer_gradient_ok()] at order 2. The
#' gate is required for correctness, since the assembly returns a number
#' wherever it is called and a wrong one where the order is not covered.
#'
#' Where the order is not covered, the exact gradient is differenced instead
#' by [statmod_hess_stencil()], which is also the route that
#' [statmod_marginal_hess()] takes for a model carrying a filter. The forms
#' that reach it are a criterion on the expected information for a family
#' that does not supply the second derivative of its expected information (or
#' for a model with a refreshable block or a structural term), and a penalty
#' whose Hessian moves with the coefficients, such as a heavy-tailed prior on
#' a random effect, placed beside a structural term, where
#' [statmod_structural_hess()] reads no movement of that Hessian. The same
#' prior without a structural term is analytic, through
#' [statmod_penalty_second()]. Both kinds of form carry an exact gradient,
#' and the difference is taken of that gradient.
#'
#' The two routes agree where both exist. The difference costs
#' \eqn{4(n_h + n_\gamma)} refits, where \eqn{n_\gamma} is the number of
#' coefficients that the criterion estimates; the analytic route refits
#' nothing.
#'
#' `source` is informative only. For a model carrying a filter
#' [statmod_marginal_hess()] falls back to [statmod_hess_stencil()] where the
#' analytic assembly fails, and `source` then still reads `"analytic"`,
#' because this function has checked only whether the order was covered.
#'
#' @param spec,design,coef,hyper,method,idx,basis As [statmod_marginal_hess()]
#'   takes them.
#' @param inner The inner optimizer with which the differenced route refits.
#' @param gam The coefficients that the criterion estimates, from
#'   [marginal_coords()], or `NULL`.
#'
#' @return A list with `A`, the symmetric negated Hessian or `NULL`;
#'   `source`, `"analytic"` or `"differenced"`, or `NA` where there is no `A`;
#'   and `why`, the reason where there is no `A`.
#'
#' @seealso [statmod_certificate()], [joint_decrement()]
#'
#' @keywords internal
outer_curvature <- function(spec, design, coef, hyper, method, idx,
                            basis = NULL, inner = NULL, gam = NULL) {
  ng <- if (is.null(gam)) 0L else length(gam$where)
  sym <- function(H) {
    if (is.null(H)) return(NULL)
    H <- as.matrix(H)
    if (!all(is.finite(H))) return(NULL)
    -(H + t(H)) / 2
  }
  if (outer_gradient_ok(spec, design, idx, method, 2L, gamma = ng > 0L)) {
    A <- sym(tryCatch(statmod_marginal_hess(spec, design, coef, hyper, method,
                                            idx, basis, inner = inner,
                                            gam = gam),
                      error = function(e) NULL))
    if (!is.null(A)) return(list(A = A, source = "analytic", why = character(0)))
  }
  # the differenced route perturbs the hyperparameters and the coefficients
  # that the criterion estimates
  A <- sym(tryCatch(statmod_hess_stencil(spec, design, coef, hyper, method, idx,
                                         basis, inner, gam = gam),
                    error = function(e) NULL))
  if (!is.null(A)) return(list(A = A, source = "differenced", why = character(0)))
  list(A = NULL, source = NA_character_, why = paste0(
    "the outer criterion's curvature could not be read at the reported ",
    "point, by the analytic route or by differencing the exact gradient"))
}



#' The Newton Decrement of a Criterion at a Point
#'
#' @description
#' Computes the Newton decrement \eqn{\tfrac{1}{2}g^\top A^{-1} g}, in the
#' criterion's own units.
#'
#' @details
#' Under a quadratic model of the criterion at the reported point, taking the
#' step \eqn{A^{-1}g} raises it by exactly this amount. The decrement is
#' therefore the part of the criterion that remains to be gained at the point,
#' in the units in which criteria are compared, and it approximates the gap
#' between the criterion at the point and at the optimum. A gradient does not
#' have that scale: its size depends on the number of penalized coefficients
#' and on the shape of the model, so that a single threshold on it cannot
#' serve every model.
#'
#' The result is `NA` where \eqn{A} is not positive definite, because the
#' point is then not a maximum in those directions and no rise follows from
#' it.
#'
#' @param g The outer criterion's gradient, over the coordinates under test.
#' @param A The negated symmetric curvature over the same coordinates.
#'
#' @return A single number, or `NA_real_`.
#'
#' @seealso [statmod_certificate()], [outer_curvature()], [coord_decrement()]
#'
#' @keywords internal
joint_decrement <- function(g, A) {
  if (!length(g)) return(0)
  if (any(!is.finite(g)) || any(!is.finite(A))) return(NA_real_)
  ev <- tryCatch(eigen(A, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev) || !all(ev > 0)) return(NA_real_)
  v <- tryCatch(solve(A, g), error = function(e) NULL)
  if (is.null(v) || any(!is.finite(v))) return(NA_real_)
  0.5 * sum(g * v)
}



#' The Flat Directions of the Outer Curvature
#'
#' @description
#' Computes the directions along which the equilibrated outer curvature
#' \eqn{A/(ss')}, \eqn{s_j = \sqrt{|A_{jj}|}}, has an eigenvalue of at most
#' `flat` in absolute value, over the coordinates `keep`.
#'
#' @details
#' A coordinate that the curvature does not resolve shows on the diagonal,
#' which [statmod_certificate()] reads one coordinate at a time. A direction
#' along which several coordinates move together can be flat while the
#' criterion is curved in every coordinate separately, and the equilibrated
#' matrix reads such a direction: its diagonal is one, so its smallest
#' eigenvalue measures how far the coordinates are from a linear dependence in
#' the curvature, whatever their scales. A direction is returned on the
#' original coordinates, scaled so that its largest entry is one in absolute
#' value, and entries below 1e-3 are set to zero. Its own decrement is
#' \eqn{(q'\tilde g)^2/(2|\lambda|)}, with \eqn{q} the equilibrated
#' eigenvector and \eqn{\tilde g = g/s}; the absolute value is taken because
#' the sign of an eigenvalue as small as this is the sign of rounding. The
#' function returns no direction when any eigenvalue is below `-flat`: the
#' point is then not a maximum, and a flat direction there is not a ridge of
#' one.
#'
#' @param g The outer gradient.
#' @param A The outer curvature, positive at a maximum.
#' @param keep The coordinates to read, those not already reported at an
#'   edge.
#' @param flat The largest eigenvalue, in absolute value, of a flat
#'   direction.
#' @param tol The largest decrement of a direction that leaves the verdict.
#'
#' @return A list with `dirs`, the flat directions as vectors over every
#'   coordinate, and `settled`, those whose own decrement is at most `tol`.
#'
#' @seealso [statmod_certificate()], [decrement_off()]
#'
#' @keywords internal
certificate_ridges <- function(g, A, keep, flat, tol) {
  out <- list(dirs = list(), settled = list())
  keep <- keep[is.finite(g[keep])]
  if (length(keep) < 2L) return(out)
  A <- as.matrix(A)
  Ak <- A[keep, keep, drop = FALSE]
  Ak <- (Ak + t(Ak)) / 2
  if (any(!is.finite(Ak))) return(out)
  s <- sqrt(abs(diag(Ak)))
  if (any(s == 0)) return(out)
  e <- tryCatch(eigen(Ak / tcrossprod(s), symmetric = TRUE),
                error = function(e) NULL)
  if (is.null(e)) return(out)
  # a ridge is a flat direction of a MAXIMUM: where the curvature is clearly
  # of the wrong sign in some direction the point is not one, and the
  # directions orthogonal to that one say nothing. Measured on a
  # zero-inflated negative binomial with a random effect on each equation,
  # stopped where the outer curvature is a rank-one matrix of order 1e14 and
  # the equilibrated eigenvalues are -5 and four of order 1e-8: every
  # direction but one read as a ridge.
  if (min(e$values) < -flat) return(out)
  for (k in which(abs(e$values) <= flat)) {
    q <- e$vectors[, k]
    xk <- q / s
    xk <- xk / max(abs(xk))
    xk[abs(xk) < 1e-3] <- 0
    x <- numeric(length(g))
    x[keep] <- xk
    out$dirs[[length(out$dirs) + 1L]] <- x
    gq <- sum(q * g[keep] / s)
    lam <- max(abs(e$values[k]), .Machine$double.xmin)
    if (gq^2 / (2 * lam) <= tol) {
      out$settled[[length(out$settled) + 1L]] <- x
    }
  }
  out
}


#' The Coefficients That Move Along a Ridge
#'
#' @description
#' Returns the rows of the coefficient table whose standard errors are
#' dropped because of the ridges that [statmod_certificate()] names.
#'
#' @details
#' A ridge is written over the coordinates of the outer criterion; the
#' coefficients among them are taken as they are. The coefficients held in
#' the inner fit move with the ridge through the mode, and the variance
#' matrix carries that: the change of coefficient \eqn{j} per unit along the
#' combination \eqn{x'\gamma} is \eqn{b_j = (Vx)_j / (x'Vx)}, and a
#' coefficient with \eqn{|b_j| \ge 10^{-2}} is taken too; the direction is
#' scaled so that its largest entry is one. Where the variance matrix lacks a
#' row for one of the ridge's coefficients or has a non-finite entry among
#' them (the mode reading already held one of them), only the coefficients
#' that the ridge names are taken.
#'
#' @param ridge The certificate's `ridge`, in which each element carries the
#'   attribute `where`, the stacked positions of its coefficients (`NA` for a
#'   hyperparameter).
#' @param V The variance matrix over the coefficients, or `NULL`.
#' @param rows The row names of the coefficient table, in stacked order.
#'
#' @return The row indices, sorted.
#'
#' @seealso [certificate_ridges()], [summary.StatmodFit()]
#'
#' @keywords internal
ridge_coef_rows <- function(ridge, V, rows) {
  out <- integer(0)
  Vd <- if (is.null(V)) NULL else as.matrix(V)
  for (x in ridge) {
    w <- attr(x, "where")
    ok <- !is.na(w)
    if (!any(ok)) next
    w <- w[ok]
    xv <- as.numeric(x)[ok]
    out <- c(out, w)
    if (is.null(Vd)) next
    nm <- rows[w]
    if (!all(nm %in% rownames(Vd))) next
    Vw <- Vd[, nm, drop = FALSE]
    if (any(!is.finite(Vd[nm, nm]))) next
    vr <- drop(crossprod(xv, Vd[nm, nm, drop = FALSE] %*% xv))
    if (!is.finite(vr) || vr <= 0) next
    b <- drop(Vw %*% xv) / vr
    hit <- names(b)[is.finite(b) & abs(b) >= 1e-2]
    out <- c(out, stats::na.omit(match(hit, rows)))
  }
  sort(unique(as.integer(out)))
}


#' The Decrement Over the Directions Orthogonal to Some Others
#'
#' @description
#' Computes the joint decrement over the coordinates `interior`, maximized
#' over the steps that have no component along the directions `dirs` in the
#' equilibrated coordinates.
#'
#' @details
#' A constraint on the step can only lower the maximum of
#' \eqn{2g'x - x'Ax}, so the result is at most [joint_decrement()] over the
#' same coordinates, as removing a coordinate is.
#'
#' @param g,A The outer gradient and curvature.
#' @param interior The coordinates under test.
#' @param dirs A list of directions, vectors over every coordinate.
#'
#' @return A single number, `NA` where the reduced curvature is not
#'   positive definite.
#'
#' @seealso [certificate_ridges()], [joint_decrement()]
#'
#' @keywords internal
decrement_off <- function(g, A, interior, dirs) {
  A <- as.matrix(A)[interior, interior, drop = FALSE]
  A <- (A + t(A)) / 2
  gi <- g[interior]
  s2 <- abs(diag(A))
  C <- do.call(cbind, lapply(dirs, function(x) s2 * x[interior]))
  r <- qr(C)$rank
  if (r >= length(interior)) return(0)
  N <- qr.Q(qr(C), complete = TRUE)[, -seq_len(r), drop = FALSE]
  joint_decrement(as.numeric(crossprod(N, gi)), crossprod(N, A %*% N))
}


#' The Decrement of One Coordinate
#'
#' @description
#' Computes \eqn{g_j^2/(2A_{jj})}, the decrement of the one-coordinate
#' problem.
#'
#' @details
#' [statmod_certificate()] reads it for the boundary label, because a verdict
#' over the whole vector does not show which coordinate has stopped moving. It
#' is the same maximum that [joint_decrement()] takes, restricted to one
#' coordinate, so it is bounded by the joint reading, and a coordinate that it
#' calls settled has already passed the test that the verdict applies. This
#' keeps `flat` a label and not a verdict. It is a weaker reading than the
#' joint decrement and is not used as a verdict.
#'
#' The result is `Inf` where the diagonal or the gradient entry is not finite
#' or the diagonal is not positive, so that such a coordinate is never called
#' settled.
#'
#' @param g The outer criterion's gradient.
#' @param A The negated symmetric curvature.
#'
#' @return A numeric vector as long as `g`.
#'
#' @seealso [joint_decrement()], [statmod_certificate()]
#'
#' @keywords internal
coord_decrement <- function(g, A) {
  d <- as.numeric(diag(as.matrix(A)))
  out <- ifelse(is.finite(d) & d > 0 & is.finite(g), g^2 / (2 * d), Inf)
  as.numeric(out)
}



#' The Two Readings That a Summary Prints About a Point
#'
#' @description
#' Computes the largest gradient in standard-error units, \eqn{\max_j
#' |g_j|/\sqrt{|A_{jj}|}}, and the smallest eigenvalue of the equilibrated
#' curvature \eqn{D^{-1/2} A D^{-1/2}} with \eqn{D = |\mathrm{diag}(A)|}.
#'
#' @details
#' Both are dimensionless, so both are independent of the units of a
#' coefficient and of the sample size. \eqn{|g_j|/\sqrt{A_{jj}}} is the length
#' of the one-coordinate Newton step measured in that coordinate's standard
#' error conditional on the others, \eqn{1/\sqrt{A_{jj}}}, which is smaller
#' than the marginal one; equivalently it is \eqn{\sqrt{2\Delta_j}}, with
#' \eqn{\Delta_j = g_j^2/(2A_{jj})} the amount that this one-coordinate step
#' would still gain in the objective's own units. The equilibrated matrix has
#' a unit diagonal, so its smallest eigenvalue is at most one: near zero it
#' reports a direction that the curvature barely identifies, and below zero a
#' point that is not an optimum.
#'
#' `A` is the curvature with the sign that makes an optimum positive
#' definite: the penalized information for the inner fit, the negated outer
#' Hessian for the criterion. A coordinate with zero curvature reads an
#' infinite gradient where its gradient is not zero, and keeps a unit scale
#' in the equilibration.
#'
#' @param g The gradient.
#' @param A The curvature, a matrix of the same order.
#' @param keep The coordinates to read.
#'
#' @return A named numeric vector, `gradient` and `eigen`, `NA` where there
#'   is nothing to read.
#'
#' @seealso [statmod_certificate()], [joint_decrement()]
#'
#' @keywords internal
certificate_readings <- function(g, A, keep = seq_along(g)) {
  out <- c(gradient = NA_real_, eigen = NA_real_)
  if (is.null(A)) return(out)
  keep <- keep[is.finite(g[keep])]
  if (!length(keep)) return(out)
  A <- as.matrix(A)[keep, keep, drop = FALSE]
  g <- g[keep]
  if (any(!is.finite(A))) return(out)
  A <- (A + t(A)) / 2
  s <- sqrt(abs(diag(A)))
  ratio <- ifelse(g == 0, 0, abs(g) / s)
  out[["gradient"]] <- max(ratio)
  s[s == 0] <- 1
  ev <- tryCatch(eigen(A / tcrossprod(s), symmetric = TRUE,
                       only.values = TRUE)$values,
                 error = function(e) NULL)
  if (!is.null(ev) && all(is.finite(ev))) out[["eigen"]] <- min(ev)
  out
}


#' The Inner Readings on the Directions That the Kinks Leave Free
#'
#' @description
#' Computes [certificate_readings()] over the null space of the linear
#' constraints that a break-point sitting on an observation imposes
#' ([kink_constraints()]), restricted to the coordinates in `keep`.
#'
#' @details
#' With no constraint this is [certificate_readings()] itself. A constraint
#' that holds one coordinate removes it from the reading. A constraint over
#' a developed position, a row \eqn{z_j^\top} of its design, leaves a
#' hyperplane, and the readings are those of the reduced system
#' \eqn{(N^\top g, N^\top A N)} with \eqn{N} an orthonormal basis of the
#' null space ([null_basis()]).
#'
#' @param g The score.
#' @param A The penalized information.
#' @param keep The coordinates read.
#' @param kc The constraint matrix, with as many columns as the coefficients
#'   (padded with zeros where `g` is longer), or `NULL`.
#'
#' @return The two readings, as [certificate_readings()] returns them.
#'
#' @keywords internal
kink_readings <- function(g, A, keep = seq_along(g), kc = NULL) {
  if (is.null(kc) || !nrow(kc)) return(certificate_readings(g, A, keep))
  B <- matrix(0, nrow(kc), length(g))
  B[, seq_len(ncol(kc))] <- kc
  keep <- keep[is.finite(g[keep])]
  B <- B[, keep, drop = FALSE]
  B <- B[rowSums(B != 0) > 0, , drop = FALSE]
  if (!nrow(B)) return(certificate_readings(g, A, keep))
  N <- null_basis(B, length(keep))
  if (!ncol(N)) return(c(gradient = NA_real_, eigen = NA_real_))
  Ak <- as.matrix(A)[keep, keep, drop = FALSE]
  certificate_readings(as.numeric(crossprod(N, g[keep])),
                       crossprod(N, Ak %*% N))
}



#' The Optimality Check of the Coefficients That a Kinked Penalty Sets to Zero
#'
#' @description
#' Computes the largest ratio \eqn{|s_j|/\kappa_j} over the coefficients that
#' a penalty with a kink sets exactly to zero, where \eqn{s_j} is the
#' derivative of the log-likelihood in that coefficient and \eqn{\kappa_j} is
#' the size of the kink there.
#'
#' @details
#' At a zero the penalty has no derivative, and the condition for an optimum
#' is not a vanishing gradient but \eqn{|s_j| \le \kappa_j}: the pull of the
#' likelihood must not exceed what the kink holds back. A ratio at or below one
#' therefore means that every zero is where it belongs, and a ratio above one
#' identifies a coefficient that would improve the objective by leaving zero.
#' The standard error does not enter, because a coefficient at zero has none.
#'
#' With the other coefficients fixed and a kink proportional to \eqn{\lambda},
#' the coefficient with the largest ratio \eqn{r < 1} leaves zero at a
#' \eqn{\lambda} that is the fraction \eqn{r} of the fitted one, while along a
#' path the other coefficients move and the threshold moves with them. The
#' ratio is not a test of significance, and values close to one are ordinary
#' where \eqn{\lambda} was chosen on a grid, because the next coefficient to
#' enter is close to its threshold by construction.
#'
#' The kink is read from the penalty's own gradient just to the right of
#' zero, so it is \eqn{\lambda} for a lasso, \eqn{\lambda\alpha} for an
#' elastic net and \eqn{\lambda} scaled by the spread under `standardize`,
#' with no formula per family. The penalties are separable, so every zero is
#' moved off in one evaluation.
#'
#' @param spec,design,obj,coef,hyper The fit's specification, design,
#'   objective, coefficients and hyperparameters.
#' @param lab [coef_labels()] of the design.
#'
#' @return `NULL` where no coefficient at zero has a finite ratio, otherwise a
#'   named vector with `ratio` and `n`, the number of coefficients at zero.
#'
#' @seealso [statmod_certificate()], [certificate_readings()]
#'
#' @keywords internal
zero_readings <- function(spec, design, obj, coef, hyper, lab) {
  b <- obj$stack(coef)
  z <- which(lab$kinked & b == 0)
  if (!length(z)) return(NULL)
  flat <- function(cf) unlist(statmod_penalty_at(spec, cf, hyper, design,
                                                 "gradient"),
                              use.names = FALSE)
  pull <- obj$gr(b)[z] - flat(coef)[z]
  bp <- b
  bp[z] <- sqrt(.Machine$double.eps) * .Machine$double.eps
  kink <- flat(obj$split(bp))[z]
  ok <- is.finite(pull) & is.finite(kink) & kink > 0
  if (!any(ok)) return(NULL)
  c(ratio = max(abs(pull[ok]) / kink[ok]), n = length(z))
}
