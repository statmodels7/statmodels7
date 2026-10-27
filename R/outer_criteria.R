#' @include outer_hessian.R
NULL

# The prediction-error criteria, and their exact derivatives.
#
# A marginal criterion integrates the coefficients out; a prediction-error one
# estimates how well the fit would predict, and every member is a function of
# two quantities:
#
#   the log-likelihood at the penalized mode,   l(b(t))
#   the effective degrees of freedom,           tau(t) = tr(J^-1 H)
#
# with J = H + S the penalized information and H the likelihood's own. An
# unpenalized model has J = H and tau = p, and every degree of freedom a
# penalty takes away shows up there.
#
# Both are differentiated by the implicit function theorem, exactly as the
# marginal criterion is, and out of the same pieces: b_m = -J^-1 c_m, the
# contractions T and U of the third and fourth derivatives of the objective in
# the coefficients, and the penalty's own theta-derivatives. What differs is
# ONE thing, and it is the thing that makes these criteria different objects:
# the envelope theorem does not apply. l alone is not stationary at the mode --
# the penalized objective is -- so dl/dt carries db/dt from the start, and at
# the mode dl/db equals drho/db, which is what makes it computable.

#' Prediction-Error Criteria for the Hyperparameters
#'
#' @description
#' Build the object that specifies, for [statmod()], that a model's
#' hyperparameters are chosen by an estimate of prediction error. Pass the
#' result as `outer_criterion` for the smooth hyperparameters, as
#' `sparse_criterion` for a lasso's or a SCAD's, or as both.
#'
#' The criterion is
#' \deqn{C(\theta) = -2\ell(\hat\beta(\theta)) + \kappa\,\tau(\theta),}
#' with \eqn{\tau = \mathrm{tr}[(H+S)^{-1}H]} the effective degrees of freedom
#' and \eqn{\kappa} the price of one of them. `aic()` charges 2 by default and
#' `bic()` charges \eqn{\log n}, with \eqn{n} the number of observations (the
#' total of the weights, which are frequencies). Both are minimized.
#'
#' @details
#' # Marginal and prediction-error criteria
#'
#' [reml()] and [ml()] integrate the coefficients out and maximize what is
#' left, which needs a penalty that is twice differentiable at the mode.
#' A kinked penalty puts the mode on the kink for every coefficient that it sets
#' to zero, where the penalized information does not exist. A marginal
#' criterion therefore leaves the coordinates that such a penalty covers out
#' of its determinant and does not select the hyperparameter of the kinked
#' penalty. These two criteria select it, and `bic()` is `statmod()`'s default
#' `sparse_criterion`.
#'
#' For a smooth penalty both families apply. A marginal criterion selects the
#' hyperparameter that makes the data most likely after the coefficients are
#' integrated out, and these select the one that predicts best.
#'
#' # Exact derivatives
#'
#' The derivatives follow from the implicit function theorem, using the same
#' pieces as the marginal criteria. The envelope theorem does not apply here:
#' \eqn{\ell} alone is not stationary at the penalized mode, so its derivative
#' carries \eqn{d\hat\beta/d\theta} from the first order. This is computable
#' because of the mode's own condition,
#' \eqn{\partial\ell/\partial\beta = \partial\rho/\partial\beta}.
#'
#' The gradient is exact with either information. The Hessian is exact with
#' `information = observed_info()`, which is the default, and is not available
#' with [expected_info()].
#'
#' # Which effective degrees of freedom
#'
#' The trace runs over the whole coefficient vector, so a term's contribution
#' reads the full penalized information, including the term's couplings with
#' the other blocks. [statmod_edf()] reads each term's share from the diagonal
#' block of the same full smoother matrix \eqn{(H+S)^{-1}H}, so the per-term
#' numbers add up to \eqn{\tau}. The count stored in the fitted model
#' (`fit@edf`) is read on the `information` of this criterion, so that it
#' reproduces the \eqn{\tau} that the criterion charged, and
#' [summary.StatmodFit()] keeps it whatever its `vcov_type` is. For a
#' positive semidefinite \eqn{H} the total lies in \eqn{[0, p]}, with \eqn{p}
#' the number of coefficients; an indefinite observed information can give a
#' direction of negative curvature a negative share.
#'
#' Where a penalty has a kink the trace is restricted to the active
#' coordinates, so for a lasso \eqn{\tau} is the number of surviving
#' coefficients, the unbiased count of Zou, Hastie and Tibshirani (2007).
#'
#' # GCV is not offered
#'
#' Classical GCV divides a residual sum of squares by \eqn{(n-\tau)^2}, which
#' estimates an unknown scale. Here every distribution parameter carries its
#' own equation and the fit profiles no parameter, so there is no unknown
#' scale for that ratio to estimate and the criterion degenerates to `aic()`.
#' That is the substitution that Wood (2008) makes in the other direction when
#' the scale is known.
#'
#' A GCV on the squared error of the fitted mean is a different and
#' well-defined object. It would need the derivative of that mean in the
#' parameters, which is not one of \pkg{distributions7}'s generics.
#'
#' # Kinked penalties
#'
#' The penalized mode is only piecewise smooth in such a hyperparameter,
#' turning a corner whenever a coefficient joins the active set or leaves it,
#' so the criterion inherits the corners and is swept over a grid. How the
#' path covers a term carrying more than one kinked
#' hyperparameter belongs to the term, through
#' [modelterms7::term_search()]. The same criterion object also selects the
#' model's smooth hyperparameters, which are not swept.
#'
#' @param k The price of one degree of freedom, a single non-negative finite
#'   number. Defaults to 2, which is Akaike's. `k = log(n)` by hand is `bic()`.
#'   A value that is not a single non-negative finite number is an error.
#'   `bic()` does not take it: its \eqn{\kappa} is `NA` on the object and is
#'   resolved to \eqn{\log n} when the model is fitted, because the sample
#'   size is unknown until then.
#' @param information The information on which the criterion is built, an
#'   object that inherits from [Information]: [observed_info()] (the default)
#'   or [expected_info()], whose `method` is `NULL` (the closed form of the
#'   family where it has one, and [distributions7::expect_opg()] otherwise) or
#'   one of [distributions7::expect_opg()], [distributions7::expect_integrate()]
#'   and [distributions7::expect_mc()] (see
#'   [distributions7::expected_derivative_methods]). [auto_info()] is not
#'   admissible, since the effective degrees of freedom are the trace on one
#'   named matrix, and signals an error. The exact gradient is available on
#'   either matrix; the Hessian is exact on the observed one and is not
#'   available on the expected one, where the search uses an optimizer that
#'   needs only the gradient. On the expected information, the `method` object
#'   defines the matrix and the derivative of it that the gradient reads, as
#'   described for the `information` of [reml()]. With `method = NULL` on a
#'   family without a closed form, the matrix is the
#'   [distributions7::expect_opg()] estimate of the expected information,
#'   noisier than a computed expectation and about 25 times cheaper (0.8 s
#'   against 25.3 s for `method = expect_integrate(identity = "score")` on a
#'   density-only family with n = 300); the computed expectation is selected
#'   with that method. [information_objects] compares the objects.
#'
#' @return An [OuterMethod()] object with properties `kind` (`"aic"` or
#'   `"bic"`), `hessian` (the object given as `information`), `k` (the number
#'   given, or `NA_real_`
#'   for `bic()`), `marginal` (`NULL`), and the path settings `nfolds`, `rule`
#'   and `folds`, which only [cv()] reads.
#'
#' @references
#' Zou, H., Hastie, T. and Tibshirani, R. (2007). On the degrees of freedom of
#' the lasso. *The Annals of Statistics* **35**(5), 2173--2192.
#'
#' Wood, S. N. (2008). Fast stable direct fitting and smoothness selection for
#' generalized additive models. *Journal of the Royal Statistical Society,
#' Series B* **70**(3), 495--518.
#'
#' @seealso [reml()] and [ml()] for the marginal criteria, [cv()] for
#'   cross-validation, [statmod()] for where these are passed.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(200, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(200, sd = 0.3)
#'
#' fa <- statmod(y ~ s(x, bspline_smooth(k = 10)), distributions7::gaussian1_distrib(), dd,
#'               outer_criterion = aic())
#' fb <- statmod(y ~ s(x, bspline_smooth(k = 10)), distributions7::gaussian1_distrib(), dd,
#'               outer_criterion = bic())
#'
#' # BIC charges log(n) = 5.3 per degree of freedom against AIC's 2, so it
#' # buys a smoother fit: a larger smoothing parameter and fewer edf.
#' c(aic = unlist(fa@hyper), bic = unlist(fb@hyper))
#' c(aic = sum(fa@edf$edf), bic = sum(fb@edf$edf))
#'
#' # The object is a specification and carries no data.
#' aic()
#' aic(k = 4)@k
#'
#' @export
aic <- function(k = 2, information = observed_info()) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 0) {
    stop("'k' must be a single non-negative number.", call. = FALSE)
  }
  check_information(information, arg = "information", fn = "aic",
                    auto = FALSE)
  do.call(OuterMethod, c(list(kind = "aic", hessian = information,
                             k = as.numeric(k)), outer_path_defaults()))
}

#' @rdname aic
#' @export
bic <- function(information = observed_info()) {
  check_information(information, arg = "information", fn = "bic",
                    auto = FALSE)
  do.call(OuterMethod, c(list(kind = "bic", hessian = information,
                             k = NA_real_), outer_path_defaults()))
}


#' Is a Criterion Minimized?
#'
#' @description
#' Reports whether a criterion is to be made small: `TRUE` for [aic()],
#' [bic()] and [cv()], which estimate prediction error, and `FALSE` for
#' [reml()] and [ml()], which are log-likelihoods and are maximized.
#'
#' @details
#' The outer search always minimizes, so the value returned here determines
#' whether the criterion enters as it stands or after a change of sign. Reading
#' the orientation from the method, and not from a list of kinds written into
#' the search, lets a criterion added later declare its own orientation in one
#' place.
#'
#' @param method An [OuterMethod()], as [reml()], [ml()], [aic()], [bic()] or
#'   [cv()] returns it. Only its `kind` is read.
#'
#' @return A single logical.
#'
#' @seealso [outer_k()] for the other property of a criterion that the search
#'   reads, [aic()] for the criteria for which this returns `TRUE`.
#'
#' @keywords internal
outer_minimize <- function(method) {
  method@kind %in% c("aic", "bic", "cv")
}


#' The Price of One Degree of Freedom
#'
#' @description
#' Returns \eqn{\kappa}, the price that a prediction-error criterion charges for
#' one effective degree of freedom, resolved against the sample size where
#' the method left it open. [bic()] stores `NA`, which this function replaces
#' by \eqn{\log n}; every other criterion returns the `k` with which it was
#' constructed.
#'
#' @details
#' The resolution happens at fit time because [bic()] takes no data, so
#' \eqn{n} is not available at construction. A criterion object is a
#' specification, reusable across models of different sizes.
#'
#' @param method An [OuterMethod()]. Its `kind` and `k` are read.
#' @param n The number of observations, a single positive number. The fit
#'   passes the total of the weights, which are frequencies.
#'
#' @return A single number: \eqn{\log n} for [bic()], and `method@k` for
#'   anything else, which is `NA_real_` for [reml()], [ml()] and [cv()],
#'   and these three criteria charge nothing per degree of freedom.
#'
#' @seealso [aic()] for where `k` is set, [outer_minimize()] for the other
#'   property read at the same point.
#'
#' @keywords internal
outer_k <- function(method, n) {
  # n is weight_total(): the weights are frequencies, so a weight of two is
  # two observations
  if (identical(method@kind, "bic")) return(log(n))
  method@k
}


#' The Effective Degrees of Freedom of a Whole Fit
#'
#' @description
#' Computes \eqn{\tau = \mathrm{tr}[(H+S)^{-1}H]}, the trace of the influence
#' matrix, which is the number of effective degrees of freedom charged by a
#' prediction-error criterion. An unpenalized fit has \eqn{S = 0} and
#' \eqn{\tau = p}; shrinkage lowers it, and \eqn{\tau} falls to the dimension
#' of the penalty's null space as the hyperparameter grows.
#'
#' @details
#' # Restricting to the active set
#'
#' A penalty with a kink is not twice differentiable everywhere, and the mode
#' of such a fit sits on the kink for every coefficient that it sets to zero. It
#' is twice differentiable away from the kink, so the trace is taken over the
#' active coordinates only, which form the submodel that the fit has selected.
#'
#' For a penalty that is linear there, the lasso, \eqn{S_{AA}} vanishes and
#' \eqn{\tau} is exactly the number of coefficients that survived, the
#' unbiased count of Zou, Hastie and Tibshirani (2007). The elastic net keeps
#' its quadratic part and spends slightly less; SCAD and MCP contribute their
#' own curvature and spend slightly more.
#'
#' The trace over the whole vector does not depend on which coefficients are at
#' the kink, so a criterion built on it assigns the same degrees of freedom to
#' a model that selects nothing and to one that selects everything.
#'
#' # The inverse
#'
#' The inverse is computed through a Cholesky factorization of \eqn{J}, and
#' the result is `NA_real_` when that fails, which a caller reads as an
#' unusable point. Both matrices are densified
#' first: the trace needs the full inverse, which is dense whatever \eqn{J}
#' was.
#'
#' @param J The penalized information \eqn{H + S}, a `p x p` symmetric
#'   matrix, dense or sparse.
#' @param H The likelihood's own information, of the same shape as `J`.
#' @param active A logical vector over the stacked coefficients marking the
#'   ones away from a kink, or `NULL` (the default) for all of them. Both
#'   matrices are subset to it before the trace.
#'
#' @return A single number: the trace, `0` when `active` selects nothing, and
#'   `NA_real_` when `J` has no Cholesky factor.
#'
#' @seealso [aic()] for the criterion that this function feeds,
#'   [statmod_edf()] for the per-term reading, which splits the same trace
#'   by term.
#'
#' @references
#' Zou, H., Hastie, T. and Tibshirani, R. (2007). On the degrees of freedom of
#' the lasso. *The Annals of Statistics* 35(5), 2173--2192.
#'
#' @keywords internal
outer_tau <- function(J, H, active = NULL) {
  J <- as_dense(J)
  H <- as_dense(H)
  if (!is.null(active)) {
    if (!any(active)) return(0)
    J <- J[active, active, drop = FALSE]
    H <- H[active, active, drop = FALSE]
  }
  P <- tryCatch(chol2inv(chol(J)), error = function(e) NULL)
  if (is.null(P)) return(NA_real_)
  sum(P * t(H))
}


#' A Prediction-Error Criterion at Given Coefficients and Hyperparameters
#'
#' @description
#' Evaluates \eqn{-2\ell + \kappa\tau} at one set of coefficients and
#' hyperparameters. The coefficients are taken as given and the function does
#' not fit the model, so the criterion at the penalized mode is obtained by
#' passing the mode.
#'
#' @details
#' Each piece is assembled once: the weighted log-likelihood, the
#' penalized information \eqn{H + S} and the likelihood's own information
#' \eqn{H}, then [outer_tau()] for the trace and [outer_k()] for \eqn{\kappa}.
#'
#' The penalty's own value is reported alongside but does not enter the
#' criterion. A prediction-error criterion prices the fit's complexity
#' through \eqn{\tau}, not through the size of the penalty.
#'
#' Where the penalty's Hessian is exactly zero on the active coordinates, as
#' happens when only a lasso penalizes them, the trace is
#' \eqn{\mathrm{tr}(H_{AA}^{-1}H_{AA}) = |A|}, the count of Zou, Hastie and
#' Tibshirani (2007), and it is returned as that count without forming
#' \eqn{H_{AA}}. The
#' condition is on the whole active block, so a ridge, a random effect or a
#' smooth among the active coordinates puts its \eqn{\lambda} into
#' \eqn{S_{AA}} and the full trace is computed as before.
#'
#' Where the active columns are linearly dependent \eqn{H_{AA}} is singular
#' and the trace does not exist; the degrees of freedom are then
#' \eqn{\mathrm{rank}(X_A)} (Tibshirani and Taylor, 2012), which the count
#' overstates by the deficiency. It happens where a development carries its
#' own intercept beside a lasso over indicators that sum to it, as in
#' `nl(~ a * exp(-r * x), a ~ 1 + lasso(~ g))`. The count is therefore taken
#' only where [design_count_exact()] certifies that every active subset has
#' full rank, and [active_rank()] supplies the degrees of freedom everywhere
#' else.
#'
#' The certificate costs one pivoted decomposition per equation and is read
#' once where the design stands still, so a lasso path pays it at its first
#' evaluation only.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, refreshed at `coef` if any term needs it.
#' @param coef A named list of coefficient vectors, one per distribution
#'   parameter.
#' @param hyper The hyperparameters, per penalized term.
#' @param method An [OuterMethod()] of kind `"aic"` or `"bic"`, read for
#'   `hessian` and for \eqn{\kappa}.
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#' @param active Which coefficients are away from a kink, as
#'   [statmod_active()] reports them, or `NULL` if no penalty has one.
#'
#' @return A list with the components below, or `NULL` when the penalized
#'   information has no Cholesky factor or the trace is not finite:
#'   \describe{
#'     \item{`value`}{the criterion, \eqn{-2\ell + \kappa\tau}.}
#'     \item{`loglik`}{the weighted log-likelihood at `coef`.}
#'     \item{`penalty`}{the penalties' total value there. It is reported for
#'       the caller and does not enter `value`.}
#'     \item{`edf`}{the trace \eqn{\tau}.}
#'   }
#'
#' @seealso [aic()] for the criterion, [outer_tau()] for the trace,
#'   [statmod_marginal()] for the marginal criteria's counterpart.
#'
#' @keywords internal
statmod_pe <- function(spec, design, coef, hyper, method,
                       emethod = NULL, active = NULL) {
  expected <- info_expected(method@hessian)
  ll <- statmod_loglik_at(spec, coef, design)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  S <- zap_nonfinite(S)
  # the trace reads the active coordinates alone, so the information is
  # formed over those columns and not over every coefficient
  if (is.null(active)) {
    H <- statmod_information_at(spec, coef, design, expected, emethod)
    tau <- outer_tau(H + S, H, NULL)
  } else if (!any(active)) {
    tau <- 0
  } else {
    act <- which(active)
    S_aa <- S[act, act, drop = FALSE]
    if (!any(as_dense(S_aa) != 0)) {
      tau <- if (design_count_exact(design)) length(act) else
        active_rank(design, act)
    } else {
      H <- statmod_information_at(spec, coef, design, expected, emethod,
                                  index = act)
      tau <- outer_tau(H + S_aa, H, NULL)
    }
  }
  if (!is.finite(tau)) return(NULL)
  k <- outer_k(method, weight_total(spec))
  list(value = -2 * ll + k * tau, loglik = ll,
       penalty = statmod_penalty_at(spec, coef, hyper, design, "value"),
       edf = tau)
}


#' Is the Count of the Active Coordinates Their Rank?
#'
#' @description
#' Returns `TRUE` where every equation's design has full column rank, which is
#' sufficient for every subset of the active coordinates to have full rank as
#' well, so that \eqn{|A|} is \eqn{\mathrm{rank}(X_A)} whatever the active set
#' turns out to be.
#'
#' @details
#' The rank condition that [statmod_pe()] needs concerns one active set at one
#' point, and checking it there would cost a decomposition per evaluation.
#' Full column rank of the whole design implies the condition for every
#' active set at once. The result is therefore computed once and reused, which
#' keeps the count free on the path for which it was introduced.
#'
#' The answer is kept in the design's `eta_memo` environment, keyed on the
#' column counts, which is safe because that environment is built fresh with
#' each design. A term registering `modelterms7::term_refresh()` recomputes
#' its block as the coefficients move, so the same column counts do not imply
#' the same columns, and a memo cannot detect the refresh and must not be
#' used; for such a design the function returns `FALSE` and [active_rank()]
#' runs instead. That route is a decomposition of \eqn{X_A}, where the trace
#' needs the cross product of the same columns.
#'
#' A design carrying a structural term has no `eta_memo` either, so there the
#' certificate is recomputed at each call and costs one decomposition per
#' equation. The branch that this serves is reached only where a penalty with
#' a kink is present and the penalty's Hessian vanishes on the active
#' coordinates.
#'
#' @param design The design, as [statmod_design()] returns it.
#'
#' @return A single logical.
#'
#' @seealso [active_rank()], which computes the rank for one active set where
#'   this function returns `FALSE`, and [block_column_rank()] for the
#'   decomposition that both read.
#'
#' @keywords internal
design_count_exact <- function(design) {
  rf <- attr(design, "refresh")
  if (!is.null(rf) && length(rf)) return(FALSE)
  mm <- attr(design, "eta_memo")
  key <- vapply(design, function(d) d$npar, integer(1))
  if (!is.null(mm) && !is.null(mm$rank_value) &&
      identical(mm$rank_key, key)) {
    return(mm$rank_value)
  }
  out <- TRUE
  for (d in design) {
    if (is.null(d$X) || !d$npar) next
    r <- block_column_rank(d$X)
    if (is.na(r) || r < ncol(d$X)) { out <- FALSE; break }
  }
  if (!is.null(mm)) { mm$rank_key <- key; mm$rank_value <- out }
  out
}


#' The Rank of the Active Columns
#'
#' @description
#' Computes \eqn{\mathrm{rank}(X_A)}, the degrees of freedom of a fit whose
#' penalty is flat on the active block (Tibshirani and Taylor, 2012), summed
#' over the equations to which the active coordinates belong.
#'
#' @details
#' The information is block diagonal over the equations wherever the family's
#' per-observation information is, and in any case a coordinate is identified
#' only if its own equation's active columns identify it, so the rank is read
#' per equation and added. Each matrix of size \eqn{n \times |A|} is formed
#' once: the columns are subset once and decomposed once.
#'
#' @param design The design, as [statmod_design()] returns it.
#' @param act The active coordinates, as positions in the stacked vector.
#'
#' @return A single number.
#'
#' @seealso [design_count_exact()], whose certificate makes this unnecessary,
#'   and [block_column_rank()].
#'
#' @keywords internal
active_rank <- function(design, act) {
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  r <- 0
  for (a in seq_along(design)) {
    loc <- act[act > offs[a] & act <= offs[a] + npar[a]] - offs[a]
    if (!length(loc)) next
    rk <- block_column_rank(design[[a]]$X[, loc, drop = FALSE])
    if (is.na(rk)) return(NA_real_)
    r <- r + rk
  }
  r
}


#' The Column Rank of One Design Block
#'
#' @description
#' Returns the number of linearly independent columns, computed by the pivoted
#' decomposition that the fit's own solve uses, so that the two agree on which
#' columns a model identifies.
#'
#' @details
#' A dense block goes through `qr()`, at the tolerance of `dqrdc2`, which is
#' the tolerance that [augmented_solve()] uses. A sparse one goes through
#' `Matrix::qr()` and its rank is counted on the JACOBI-EQUILIBRATED diagonal
#' of the triangular factor, which is the correction that
#' [sparse_augmented_solve()] already carries: since \eqn{R'R = A'A}, scaling
#' the columns by their norms scales that diagonal by the same factors, so a
#' block whose columns differ in size is not read as deficient while an exact
#' dependence stays exactly singular.
#'
#' @param X A design block, dense or sparse.
#'
#' @return A single number, or `NA_integer_` where the decomposition fails.
#'
#' @seealso [design_count_exact()], [active_rank()].
#'
#' @keywords internal
block_column_rank <- function(X) {
  if (is.null(X) || !ncol(X)) return(0L)
  if (!inherits(X, "Matrix")) {
    q <- tryCatch(qr(X), error = function(e) NULL)
    return(if (is.null(q)) NA_integer_ else q$rank)
  }
  qrA <- tryCatch(Matrix::qr(X), error = function(e) NULL)
  if (is.null(qrA)) return(NA_integer_)
  Rf <- tryCatch(Matrix::qrR(qrA, backPermute = FALSE), error = function(e) NULL)
  if (is.null(Rf) || nrow(Rf) != ncol(X)) return(NA_integer_)
  dg <- abs(Matrix::diag(Rf))
  q <- qrA@q
  ord <- if (length(q)) q + 1L else seq_len(ncol(X))
  cn <- sqrt(Matrix::colSums(X^2))[ord]
  if (length(cn) != length(dg) || !all(is.finite(cn)) || !all(is.finite(dg))) {
    return(NA_integer_)
  }
  keep <- cn > 0
  if (!any(keep)) return(0L)
  eq <- dg[keep] / cn[keep]
  if (!max(eq) > 0) return(0L)
  sum(eq > ncol(X) * .Machine$double.eps * max(eq))
}


#' The Exact Derivatives of a Prediction-Error Criterion
#'
#' @description
#' Computes the gradient and, at order 2, the Hessian of
#' \eqn{-2\ell + \kappa\tau} over the free scale of the hyperparameters.
#'
#' @details
#' Writing \eqn{P = J^{-1}}, \eqn{A_m = \partial J/\partial\theta_m = S_m +
#' T[\hat\beta_m]} and \eqn{B_m = \partial H/\partial\theta_m = T[\hat\beta_m]},
#' \deqn{\tau_m = \mathrm{tr}(PB_m) - \mathrm{tr}(PA_mPH),}
#' and the term \eqn{-2\ell} contributes
#' \eqn{-2(\partial\rho/\partial\beta)'\hat\beta_m}. Differentiating once
#' more brings in \eqn{A_{ml}}, \eqn{B_{ml}} and \eqn{\hat\beta_{ml}}, which
#' are quantities that [statmod_marginal_hess()] assembles.
#'
#' Where `method@hessian` is `"expected"`, \eqn{\tau} reads the expected
#' information \eqn{H_E}, so its derivative contracts the derivative of
#' \eqn{H_E} in the predictors, read from
#' [distributions7::distrib_dexpected_hessian()] under its own key, while the
#' mode still moves by the observed penalized information. That route has no
#' Hessian, and a request for order 2 signals an error. Beside an unpenalized
#' filter the mode moves jointly with the filter's own parameters, and only
#' the gradient is available.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param idx The outer index.
#' @param order `1` for the gradient alone, `2` for both.
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#'
#' @return A list with `grad` and, at order 2, `hess`; or
#'   `NULL` where the information is not invertible.
#'
#' @seealso [aic()], [statmod_marginal_grad()]
#'
#' @keywords internal
statmod_pe_derivs <- function(spec, design, coef, hyper, method, idx,
                              order = 1L, emethod = NULL) {
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  nh <- nrow(idx)
  n <- spec@n_obs
  kap <- outer_k(method, weight_total(spec))
  expected <- info_expected(method@hessian)
  if (expected && order >= 2L) {
    stop("the expected route of a prediction-error criterion has no exact",
         " Hessian; use hessian = \"observed\" where the second derivative",
         " is needed.", call. = FALSE)
  }

  # the OBSERVED information whatever the criterion reads: the mode is where
  # the penalized likelihood's score vanishes, so this is what it moves by
  H <- statmod_information_at(spec, coef, design, expected = FALSE)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  S <- zap_nonfinite(S)
  # The routes below form a FULL inverse, and the inverse of a sparse matrix
  # is dense, so densifying here gives up nothing sparsity could have kept.
  J <- as_dense(H + S)
  Jfac <- tryCatch(chol(J), error = function(e) NULL)
  if (is.null(Jfac)) return(NULL)
  P <- chol2inv(Jfac)

  th <- statmod_eta(spec, design, coef)$theta
  pieces <- outer_pieces(spec, design, coef, hyper, idx, offs, total, order)
  grho <- statmod_penalty_at(spec, coef, hyper, design, "gradient")
  grho <- unlist(grho[params], use.names = FALSE)

  # ⚠️ AN UNPENALIZED FILTER MOVES WITH THE MODE HERE TOO, for the reason
  # statmod_marginal_grad() records: the mode moves jointly in the
  # coefficients and the filter's own parameters, and the information moves
  # along the total derivative of every predictor. tau reads the information
  # and nothing else, so the same two pieces are all this needs -- the solve
  # with the joint curvature and the joint rows -- and the log-likelihood's
  # term reads the coefficients' part of the movement, the filter's own score
  # vanishing at the mode. Read on the coefficients alone the gradient was out
  # by 9.48 under aic() and 2.22 under bic() on a smooth beside an unpenalized
  # gas(1, 1), relative and FLAT in the step.
  d3j <- NULL
  jm <- NULL
  if (any(vapply(attr(design, "structural"),
                 function(su) identical(su$kind, "filter"), logical(1)))) {
    if (order >= 2L) {
      stop("a prediction-error criterion beside a filter has no exact Hessian;",
           " only the gradient is available there.", call. = FALSE)
    }
    d3j <- distributions7::distrib_deriv3(spec@distrib, spec@response, th,
                                          scale = "link",
                                          threads = spec@threads)
    jm <- filter_joint_movement(spec, design, coef, S,
                                matrix(0, total, total), total, d3j)
    if (is.null(jm)) return(NULL)
  }
  ib <- seq_len(total)
  bhat <- lapply(seq_len(nh), function(m) {
    if (is.null(jm)) -as.numeric(P %*% pieces$c[[m]]) else
      -jm$solve(pieces$c[[m]])
  })
  tv <- lapply(bhat, function(v) {
    if (is.null(jm)) block_predictors(design, params, npar, offs, v) else
      lapply(jm$rows, function(R) as.numeric(R %*% v))
  })
  # how J moves through a penalty whose Hessian moves with the coefficients,
  # which A_m carries beside S_m and T[b_m] -- NULL where every penalty is
  # quadratic in them, and then nothing below changes
  Pm <- lapply(bhat, function(v) statmod_penalty_dbeta(spec, design, coef,
                                                       hyper, v, total))
  add_pm <- function(A, m) if (is.null(Pm[[m]])) A else A + Pm[[m]]

  # ⚠️ THE TRACE ON THE CRITERION'S OWN INFORMATION. aic(information = expected_info())
  # prices tau = tr[(H_E + S)^-1 H_E], so its derivative reads H_E and how H_E
  # moves with the mode, which is dE[l'']/deta and not l''': differentiating an
  # expectation moves the measure as well as the integrand. Only the trace
  # changes matrix; the mode's movement above keeps the observed one. Before
  # this the route differentiated the OBSERVED trace for a criterion reading
  # the expected one. The array and its key are the ones
  # statmod_marginal_grad() reads on the same route through ctx_kmove().
  if (expected) {
    HE <- statmod_information_at(spec, coef, design, expected = TRUE, emethod)
    JEfac <- tryCatch(chol(as_dense(HE + S)), error = function(e) NULL)
    if (is.null(JEfac)) return(NULL)
    PE <- chol2inv(JEfac)
    dE <- distributions7::distrib_dexpected_hessian(
      spec@distrib, spec@response, th, scale = "link",
      method = resolve_emethod(spec@distrib, emethod))
    keyE <- function(a, b, k) distributions7::dexpected_key(params, a, b, k)
    PEH <- PE %*% HE
    g <- numeric(nh)
    for (m in seq_len(nh)) {
      BE <- contract3(spec, design, dE, params, npar, offs, total, tv[[m]],
                      key = keyE)
      PAE <- PE %*% add_pm(pieces$S[[m]] + BE, m)
      tau_m <- sum(PE * t(BE)) - sum(PAE * t(PEH))
      g[m] <- -2 * sum(grho * bhat[[m]][ib]) + kap * tau_m
    }
    return(list(grad = free_scale(g, hyper, idx)))
  }

  d3 <- if (!is.null(d3j)) d3j else
    distributions7::distrib_deriv3(spec@distrib, spec@response, th,
                                   scale = "link", threads = spec@threads)
  d4 <- if (order >= 2L) {
    distributions7::distrib_deriv4(spec@distrib, spec@response, th,
                                   scale = "link", threads = spec@threads)
  } else NULL
  Bm <- lapply(tv, function(t) contract3(spec, design, d3, params, npar, offs,
                                         total, t))
  Am <- lapply(seq_len(nh), function(m) add_pm(pieces$S[[m]] + Bm[[m]], m))
  PH <- P %*% H

  # one per hyperparameter, and the gradient and the Hessian both read them
  PA <- lapply(Am, function(A) P %*% A)

  g <- numeric(nh)
  for (m in seq_len(nh)) {
    tau_m <- sum(P * t(Bm[[m]])) - sum(PA[[m]] * t(PH))
    g[m] <- -2 * sum(grho * bhat[[m]][ib]) + kap * tau_m
  }
  if (order < 2L) return(list(grad = free_scale(g, hyper, idx)))

  # These depend on ONE hyperparameter and were formed inside the PAIR loop,
  # so each dense p by p product was taken nh + 1 times over: measured, 6
  # against 2 at two hyperparameters and 12 against 3 at three. It is the
  # hoist statmod_marginal_hess received in 0.49.0 -- "the two expensive
  # quantities each depend on ONE direction and not on a pair, so they are
  # built once per hyperparameter" -- which the prediction-error route never
  # got, that release having deliberately left aic(), bic() and cv() alone.
  # Nothing about the arithmetic changes: the same product is taken once
  # instead of many times, so the answer is identical to the bit.
  PAPH <- lapply(PA, function(M) M %*% PH)
  PB <- lapply(Bm, function(B) P %*% B)
  Hb <- lapply(bhat, function(v) as.numeric(H %*% v))

  Hm <- matrix(0, nh, nh)
  for (m in seq_len(nh)) {
    for (l in m:nh) {
      key <- pieces$pair[m, l]
      Bml <- contract4(spec, design, d4, params, npar, offs, total, tv[[l]],
                       tv[[m]])
      rhs <- Am[[l]] %*% bhat[[m]] + pieces$S[[m]] %*% bhat[[l]] +
        pieces$c2[[key]]
      bml <- -as.numeric(P %*% rhs)
      Tml <- contract3(spec, design, d3, params, npar, offs, total,
                       block_predictors(design, params, npar, offs, bml))
      Bml <- Bml + Tml
      Aml <- pieces$S2[[key]] + Bml
      # J moves through a penalty whose Hessian depends on the coefficients as
      # well as through H, so its second derivative carries that penalty's
      # second movement; H's does not, so Bml is left as it is
      E2 <- statmod_penalty_second(spec, design, coef, hyper, idx, m, l,
                                   bhat[[m]], bhat[[l]], bml, total)
      if (!is.null(E2)) Aml <- Aml + E2

      # -2 l(b(t)) : the second derivative carries the curvature of the
      # log-likelihood along the two directions and the mode's own second
      # movement, the score at the mode being drho/dbeta
      ll_ml <- -sum(bhat[[l]] * Hb[[m]]) + sum(grho * bml)

      tau_ml <- sum(PA[[l]] * t(PAPH[[m]])) - sum(P * t(Aml %*% PH)) +
        sum(PA[[m]] * t(PAPH[[l]])) - sum(PA[[m]] * t(PB[[l]])) -
        sum(PA[[l]] * t(PB[[m]])) + sum(P * t(Bml))

      v <- -2 * ll_ml + kap * tau_ml
      Hm[m, l] <- v
      Hm[l, m] <- v
    }
  }
  list(grad = free_scale(g, hyper, idx),
       hess = free_scale2(Hm, g, hyper, idx))
}


#' Carry a Hyperparameter Derivative onto the Free Scale
#'
#' @description
#' Carries a hyperparameter derivative onto the free scale by the diagonal
#' chain rule that each hyperparameter's own link induces, at first order for
#' a gradient and at second order for a Hessian.
#'
#' @param g A gradient on the parameter scale.
#' @param M A Hessian on the parameter scale.
#' @param hyper The hyperparameters.
#' @param idx The outer index.
#'
#' @return For `free_scale()` a vector, for `free_scale2()` a matrix, and for
#'   `link_slopes()` a list with the components `h1` and `h2`, the first and
#'   second derivatives of each hyperparameter's inverse link at its value.
#'
#' @keywords internal
free_scale <- function(g, hyper, idx) {
  g * link_slopes(hyper, idx)$h1
}

#' @rdname free_scale
#' @keywords internal
free_scale2 <- function(M, g, hyper, idx) {
  s <- link_slopes(hyper, idx)
  out <- M * outer(s$h1, s$h1)
  diag(out) <- diag(out) + s$h2 * g
  out
}

#' @rdname free_scale
#' @keywords internal
link_slopes <- function(hyper, idx) {
  links <- attr(idx, "links")
  nh <- nrow(idx)
  h1 <- numeric(nh)
  h2 <- numeric(nh)
  for (r in seq_len(nh)) {
    v <- hyper[[idx$parameter[r]]][[idx$term[r]]][[idx$name[r]]]
    e <- linkfunctions7::linkfun(links[[r]], v)
    h1[r] <- linkfunctions7::dlinkinv(links[[r]], e)
    h2[r] <- linkfunctions7::d2linkinv(links[[r]], e)
  }
  list(h1 = h1, h2 = h2)
}
