#' @include outer.R
NULL

# The exact gradient of the marginal criterion.
#
# Writing V(t) = l(b) - rho(b; t) + (q/2) log 2pi - (1/2) log|A'KA| at the
# penalized mode b(t), with K = H + S:
#
#   dV/dt_m = -drho/dt_m - (1/2) [ tr(M dS/dt_m) + u' v_m ]
#
# where M = A(A'KA)^-1 A' (which is K^-1 for REML), v_m = db/dt_m and
# u_c = tr(M dK/db_c).
#
# The first two terms of V contribute only -drho/dt_m, by the envelope
# theorem: their derivative through b carries the factor dl/db - drho/db,
# which is the score of the penalized objective and vanishes at the mode.
# Nothing is left of db/dt there, and it survives only inside the determinant.
#
# v_m comes from differentiating the stationarity condition:
#   (H_obs + S) db/dt_m + d2rho/db dt_m = 0.
#
# u is the contraction of the third derivative of the log-likelihood in the
# link-scale predictors with the diagonal of the "hat" matrix of M, one
# crossprod per distribution parameter, so nothing of order p^3 beyond the
# factorization of K that the criterion has already paid for.

#' Can the Exact Derivatives of the Criterion Be Computed Here?
#'
#' @description
#' Returns `TRUE` when the exact gradient of the marginal criterion in the
#' hyperparameters under estimation (`order = 1`), or the exact gradient and
#' Hessian (`order = 2`), can be assembled for the model.
#'
#' @details
#' The family must carry the derivatives that the order needs
#' ([order_available()]): the third derivative of the log-likelihood for the
#' gradient and the fourth for the Hessian.
#'
#' For the information, \eqn{K} enters the criterion through its determinant,
#' so the gradient needs \eqn{\partial K/\partial\beta}. With the observed
#' information that is the third derivative of the log-likelihood in the
#' link-scale predictors. With the expected information it is the derivative in
#' the predictors of \eqn{-E[\ell'']}, which is not \eqn{-E[\ell''']}. That
#' derivative is supplied by [distributions7::distrib_dexpected_hessian()], and
#' the expected route is admitted when the family returns it
#' ([expected_deriv_ok()]). At order 2 the expected route is admitted only for
#' [reml()] and [ml()], for a model without a block that moves with its
#' coefficients and without a structural term, and for a family that returns
#' [distributions7::distrib_d2expected_hessian()] ([expected_deriv2_ok()]).
#' The expected route is not admitted at either order where a penalty covers
#' the parameters of a structural term, because the criterion is then assembled
#' from the observed joint curvature.
#'
#' For structural terms, a term of the likelihood shape (`regime()` and the
#' marginal break-point term) has no exact derivatives at either order. For a
#' term of the filter shape, order 2 requires a penalty over the term's own
#' parameters and a term that implements [modelterms7::term_fourth()]. Where a
#' penalty covers the term's parameters, the term must also implement
#' [modelterms7::term_third()].
#'
#' For moving blocks, a block that moves with its coefficients (`nl()`,
#' `seg()`, `jump()` and `jseg()`) is admitted at order 1 on either route and
#' at order 2 on the observed route. Its contribution is read from the term
#' ([u_refresh()], [mode_curvature()]) and is never differenced.
#'
#' For penalties, \eqn{\partial S/\partial\theta} and its second
#' derivative are given by the generics of \pkg{penalties7}
#' ([penalties7::penalty_dhessian()], [penalties7::penalty_d2hessian()] and
#' [penalties7::penalty_dcross()]). The function calls the derivative methods
#' of a penalty at probe values and reads an error as a penalty that does not
#' supply them, so that a penalty is admitted according to what its methods
#' return and not according to its class. The check has two parts. Every
#' penalized unit of the model is examined, because the movement of the mode
#' reaches the determinant through each penalty whose Hessian moves with the
#' coefficients (a Student t prior, for instance). A penalty for which
#' [penalties7::beta_quadratic()] returns `TRUE` passes, and so does a penalty
#' with a kink (a lasso, a SCAD or an MCP). Any other penalty must supply
#' [penalties7::penalty_dhessian_beta()], and a mixed covariance class of this
#' kind is not admitted. Beside a structural term such a penalty is not
#' admitted at order 2. The second part examines, with [penalty_answers()], the
#' penalty of every term that owns a hyperparameter under estimation, and for a
#' hyperparameter shared by several terms the penalty of each of them. A
#' penalty with a kink does not remove the exact derivatives of the other
#' hyperparameters, because its own hyperparameters are left out of the outer
#' index.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param idx The outer index, from [outer_hyper_index()].
#' @param method An [OuterMethod()].
#' @param order `1` for the gradient, `2` for the Hessian as well.
#' @param gamma `TRUE` where the outer vector also carries coefficients that the
#'   criterion estimates ([marginal_coords()]), so an empty index is not a
#'   reason to return `FALSE`.
#'
#' @return `TRUE` or `FALSE`.
#'
#' @seealso [statmod_marginal_grad()]
#'
#' @keywords internal
outer_gradient_ok <- function(spec, design, idx, method, order = 1L,
                              gamma = FALSE) {
  if (!nrow(idx) && !isTRUE(gamma)) return(FALSE)
  # THE ORDER THE FAMILY CARRIES.  The exact gradient reads dK/dbeta, which is
  # the family's THIRD derivative, and the Hessian reads the fourth; a family
  # with fewer is not refused here, the search falling back on a difference,
  # so this is a gate and not an error.  See R/order.R for the scale.
  if (!order_available(spec@distrib, 2L + order)) return(FALSE)
  # ⚠️ A TERM OF THE LIKELIHOOD SHAPE -- regime() -- IS REFUSED AT BOTH ORDERS,
  # penalized or not. Its parameters are estimated beside the coefficients and
  # move with the hyperparameter, and the determinant reads the information at
  # predictors weighted by the smoothed posterior, so the exact gradient needs
  # how that posterior moves along the direction the mode moves in, which
  # modelterms7 does not expose. Read on the coefficients alone it was
  # admitted and out: measured on a ridge beside regime(k = 2) at 500
  # observations, 2.17e-04, 2.06e-04 and 2.05e-04 relative at h of 1e-2, 3e-3
  # and 1e-3 against a central difference of the criterion with the mode
  # refitted, FLAT in the step; decomposed, -1.04e-03 from the mode's movement,
  # +1.53e-04 from the posterior weights and +1.24e-05 from the term's own
  # path, against a gap of -8.71e-04.
  #
  # The refusal is not free, and what it costs is stated rather than found:
  # the search this package chooses becomes nelder_mead(), 8.48 s against
  # 5.62 s on that model at the same criterion to 1e-6, and
  # statmod_certificate() answers `unknown` on every model carrying such a
  # term, there being no exact gradient to read the decrement with. It stays
  # until the term supplies the derivative of its posterior.
  for (su in attr(design, "structural")) {
    if (identical(su$kind, "loglik")) return(FALSE)
  }
  # an "auto" that statmod() has not settled reads as "observed", as it does
  # at every other site
  if (info_expected(method@hessian)) {
    # The expected route asks for the same object in the same place -- the
    # movement of K with the coefficients -- but that object is
    # dE[l'']/deta rather than l''', because differentiating an expectation
    # moves the measure as well as the integrand. distributions7 supplies it
    # for the closed form of a family and for each method object, so the
    # derivative is probed with the method that the criterion carries, a family
    # without a closed form being read as expect_opg() (see resolve_emethod()).
    #
    # ORDER 2 reads the next order of the same object,
    # d2E[l'']/deta^2, which a family supplies only where it wrote it out
    # (distributions7::distrib_d2expected_hessian(), from compiled kernels;
    # the rest refuse). It is admitted there for the marginal criteria alone:
    # the prediction-error route's order-2 assembly is written for the
    # observed information, a block that moves with its coefficients would
    # need the expected information's movement through dX/dbeta, which nothing
    # supplies, and a structural term's criterion has no expected twin.
    # Everywhere else statmod_marginal_hess() differences the exact gradient.
    if (order >= 2L) {
      if (!method@kind %in% c("ml", "reml")) return(FALSE)
      if (length(attr(design, "refresh"))) return(FALSE)
      if (length(attr(design, "structural"))) return(FALSE)
      if (!expected_deriv2_ok(spec@distrib,
                              resolve_emethod(spec@distrib,
                                              info_emethod(method@hessian)))) {
        return(FALSE)
      }
    }
    emethod <- resolve_emethod(spec@distrib, info_emethod(method@hessian))
    if (!expected_deriv_ok(spec@distrib, emethod)) return(FALSE)
    # and not where a structural term is penalized: there the criterion is
    # statmod_marginal_full(), which assembles the joint curvature from
    # term_curvature() -- the OBSERVED one -- so that branch has no expected
    # criterion for this to be the derivative of.
    if (structural_penalized(spec, design)) return(FALSE)
  }
  # ⚠️ EVERY PENALTY IS ASKED, not only those whose hyperparameters are under
  # estimation. The mode's movement reaches the determinant through each
  # penalty whose Hessian moves with the coefficients, a Student t prior held at
  # given values among them, so one that cannot supply that movement leaves the
  # gradient without it -- and the order-2 assembly is written for penalties
  # whose Hessian does not move at all. A MIXED class is read by neither route
  # in that case, its prior being multivariate, and is refused.
  for (un in statmod_penalized(spec, design)) {
    pen <- un$penalty
    if (is.null(pen)) next
    ok <- tryCatch({
      th <- as.list(penalty_theta_start(pen))
      k <- as.integer(pen@n_coef)
      if (isTRUE(penalties7::beta_quadratic(pen, th)) ||
          length(penalties7::penalty_kinks(pen, th))) {
        TRUE
      } else if (isTRUE(un$mixed)) {
        FALSE
      } else if (order >= 2L && length(attr(design, "structural"))) {
        # the joint order-2 assembly, statmod_structural_hess(), reads no
        # movement of the penalty's Hessian, so a moving penalty beside a
        # structural term keeps the stencil there
        FALSE
      } else {
        b <- rep(0.37, k)
        v <- rep(0.21, k)
        penalties7::penalty_dhessian_beta(pen, b, th, v)
        if (order >= 2L) {
          # ⚠️ ORDER 2 READS THE SECOND MOVEMENT of that Hessian, along two
          # directions and in each hyperparameter -- statmod_penalty_second().
          # Without it the assembled Hessian on a Student t prior over 30
          # groups was 122 per cent out and flat in the step.
          penalties7::penalty_d2hessian_beta(pen, b, th, v, rep(-0.13, k))
          penalties7::penalty_dhessian_beta_theta(pen, b, th, v)
        }
        TRUE
      }
    }, error = function(e) FALSE)
    if (!isTRUE(ok)) return(FALSE)
  }
  mem <- index_members(idx)
  # A SHARED hyperparameter is exact at BOTH orders. The row's derivative
  # is the sum of its members', and outer_pieces() keys the second order by
  # the index ROW pair and accumulates, so a row standing for members of two
  # terms holds the sum of their tables where a key naming the pair of
  # hyperparameter NAMES could not.
  seen <- unique(paste(mem$parameter, mem$term, sep = "\r"))
  for (s in seen) {
    bits <- strsplit(s, "\r", fixed = TRUE)[[1L]]
    u <- statmod_unit(spec, design, bits[1L], bits[2L])
    if (is.null(u) || !penalty_answers(u$penalty, order)) return(FALSE)
  }
  # Where a penalty covers a STRUCTURAL term's own parameters the determinant
  # spans those parameters too, so the chain term reads how the recursion's
  # own second derivative moves with the mode -- the third derivative of the
  # predictor, contracted in the one direction the mode moves in. A term that
  # answers term_third() supplies it; one that does not leaves the search
  # derivative-free rather than reporting a gradient missing a piece.
  #
  # The order-2 route is NOT extended: the criterion's own second derivative
  # would ask for a fourth order through the recursion, and lbfgs on the exact
  # gradient already buys most of what newton would.
  # ⚠️ A block that MOVES with the coefficients -- nl(), seg(), jump(), jseg()
  # -- is covered at NEITHER order, and this does not refuse it. dK/dbeta gains
  # everything coming from dX/dbeta, which is the term's own SECOND derivative,
  # and the criterion's own second derivative asks for the third. The layer
  # cannot difference the block to get them: measured at h, h/4 and h/16, nl's
  # converges (0.6038 throughout) and seg's break-point column diverges as 1/h
  # (3.6e4, 1.4e5, 5.8e5), being a step function in psi. The term has to supply
  # them, which is piano_nl_derivate.txt.
  #
  # It is left admitted because the measured consequence is small and refusing
  # is not free: the gradient is out by 1.4e-03 (nl) and 1.1e-02 (seg) and the
  # search still reaches the same hyperparameter to 0.01-0.14 per cent, while
  # refusing costs 1.3x on those models alone and 3.6x to 7.4x beside a smooth
  # -- the refusal being model-wide, u is shared by every hyperparameter, so
  # the smooth would lose its gradient too.
  # A MIXED covariance class is answered at order 1 and refused at order 2,
  # which is the structural branch's own rule below rather than one of its
  # own: statmod_structural_grad() assembles on the joint vector already, so
  # what such a class needed was its positions in it, while
  # statmod_marginal_hess() is written over the stacked coefficients and has
  # no joint twin at all. Measured against a central difference of the
  # criterion with the mode refitted, the order-1 route converges O(h^2) --
  # 3.1e-05, 2.8e-06, 3.1e-07 at h of 1e-2, 3e-3, 1e-3 -- where it used to
  # return exactly zero in every coordinate, which reads as stationarity.
  # ⚠️ THE ORDER-2 ROUTE IS REFUSED WHEREVER A STRUCTURAL TERM IS PRESENT,
  # penalized or not, and the unpenalized case is the one this widened. The
  # assembly is over the stacked COEFFICIENTS, and a filter's own parameters
  # are estimated beside them and move with the hyperparameter too, so what it
  # returns is not the criterion's second derivative. Measured against a
  # SECOND difference of the criterion with the mode refitted: on an
  # unpenalized filter beside a smooth, where this used to answer TRUE, the
  # assembled Hessian reads -1.927925 where the criterion's is -1.944612, 0.86
  # per cent out and flat in the step; on a PENALIZED filter it reads -1.1e-08
  # where the criterion's is -1.971456, which is the case below.
  #
  # ⚠️ THE PARAGRAPH ABOVE IS THE STATE BEFORE THE FOURTH ORDER EXISTED, and
  # is kept because what replaced it has to be read against it.
  # statmod_structural_hess() assembles that Hessian on the JOINT vector,
  # reading modelterms7::term_fourth() and the family's fifth derivative, so
  # a term that answers the fourth order is exact at both orders. One that
  # answers only the third -- and regime(), which answers neither -- keeps
  # statmod_hess_stencil(), one central difference of the exact gradient.
  #
  # The SEARCH is still not given newton() here, and that is a measurement
  # rather than a limitation: lbfgs() on the exact gradient reaches the same
  # criterion in 15 evaluations against 86, and 112 against 245, so a second
  # order would have to beat a method starting three times ahead while paying
  # a fourth-order recursion per pair at every evaluation. What the analytic
  # route buys is the REPORTING -- see statmod_structural_hess().
  if (order >= 2L && length(attr(design, "structural"))) {
    # ⚠️ AND NOT WHERE THE TERM IS UNPENALIZED, which this used to answer TRUE
    # for a gas term. statmod_structural_hess() is the second derivative of the
    # criterion whose determinant spans the term's parameters, and with no
    # penalty over them the criterion's does not: measured on a smooth beside
    # an unpenalized gas(1, 1), it read -3.97422 where a second difference of
    # the criterion reads -4.00272, 0.71 per cent out and flat in the step.
    # statmod_hess_stencil() on the exact gradient reads -4.00270, so a
    # certificate there is differenced -- which costs refits the analytic
    # route did not, and reads the right curvature where that one did not.
    # The prediction-error route is refused for the same reason: its order-2
    # assembly is over the stacked coefficients, as the gradient's was.
    if (!structural_penalized(spec, design)) return(FALSE)
    tm <- structural_term_of(spec, design)
    if (is.null(tm) || !answers_term_fourth(tm)) return(FALSE)
  }
  if (structural_penalized(spec, design)) {
    for (u in statmod_penalized(spec, design)) {
      # a MIXED class runs through the same chain term, so the term it reaches
      # has to answer term_third() for the same reason -- and it is named on
      # the class rather than on the unit, the unit's own `term` being the
      # class's key
      tm <- if (isTRUE(u$mixed)) find_term(spec, u$class$sterm)
            else if (isTRUE(u$structural)) spec@terms[[u$param]][[u$term]]
            else next
      if (is.null(tm) || !answers_term_third(tm)) return(FALSE)
    }
  }
  TRUE
}


#' A Term by Its Key in Any Equation
#'
#' @description
#' Returns the term that a key names, searching over every distribution
#' parameter.
#'
#' @details
#' A unit records the parameter under which its penalty is filed. For a
#' covariance class that is the parameter by which the class is keyed, which is
#' not necessarily the parameter carrying the structural term into which the
#' class reaches. A caller that has a term's key and needs the term should use
#' this function instead of assuming that the two agree.
#'
#' @param spec A [StatmodSpec()].
#' @param key A term's key, its call as written.
#'
#' @return The term, or `NULL` if no equation carries that key.
#'
#' @keywords internal
find_term <- function(spec, key) {
  for (p in spec@distrib@params) {
    tm <- spec@terms[[p]][[key]]
    if (!is.null(tm)) return(tm)
  }
  NULL
}


#' Does the Family Supply the Expected Information's Derivative?
#'
#' @description
#' Determines whether [distributions7::distrib_dexpected_hessian()] returns a
#' value for this family. The generic is called once, at a probe value, and the
#' result is not inferred from the family's class.
#'
#' @details
#' The default method in \pkg{distributions7} is one central difference of the
#' family's expected information. It applies where that information is a
#' written-out formula and signals an error where the information is itself
#' approximated, because differencing an approximation would be a difference of
#' a difference. The function returns `FALSE` for a family on which the generic
#' signals an error, and also where a random parameter vector or a response
#' cannot be drawn for the probe.
#'
#' @param distrib A \pkg{distributions7} distribution.
#'
#' @return A single logical.
#'
#' @seealso [outer_gradient_ok()]
#'
#' @keywords internal
expected_deriv_ok <- function(distrib, emethod = NULL) {
  th <- tryCatch(distributions7::generate_random_theta(distrib),
                 error = function(e) NULL)
  if (is.null(th)) return(FALSE)
  y <- tryCatch(distributions7::distrib_rng(distrib, 1L, th),
                error = function(e) NULL)
  if (is.null(y)) return(FALSE)
  ok <- tryCatch({
    distributions7::distrib_dexpected_hessian(distrib, y, th, scale = "link",
                                              method = emethod)
    TRUE
  }, error = function(e) FALSE)
  isTRUE(ok)
}


#' Does the Family Supply the Expected Information's Second Derivative?
#'
#' @description
#' Determines whether [distributions7::distrib_d2expected_hessian()] returns a
#' value for this family. The generic is called once, at a probe value.
#'
#' @details
#' The default method in \pkg{distributions7} signals an error, because
#' differencing the first derivative of the expected information would be a
#' difference of a difference. The function therefore returns `TRUE` only for a
#' family that provides the derivative itself. For any other family the order-2
#' Hessian of the expected route is computed by [statmod_hess_stencil()].
#'
#' @param distrib A \pkg{distributions7} distribution.
#'
#' @return A single logical.
#'
#' @seealso [outer_gradient_ok()], [expected_deriv_ok()]
#'
#' @keywords internal
expected_deriv2_ok <- function(distrib, emethod = NULL) {
  th <- tryCatch(distributions7::generate_random_theta(distrib),
                 error = function(e) NULL)
  if (is.null(th)) return(FALSE)
  y <- tryCatch(distributions7::distrib_rng(distrib, 1L, th),
                error = function(e) NULL)
  if (is.null(y)) return(FALSE)
  ok <- tryCatch({
    distributions7::distrib_d2expected_hessian(distrib, y, th, scale = "link",
                                               method = emethod)
    TRUE
  }, error = function(e) FALSE)
  isTRUE(ok)
}


#' Does a Term Supply Its Third Derivative?
#'
#' @description
#' Determines whether the term implements [modelterms7::term_third()] from the
#' class on which the method is registered, never from a list of class names,
#' so a term written later is covered without an edit here.
#'
#' @details
#' A structural term that has not implemented the method inherits one that
#' signals an error, and an additive term inherits one that returns zero, which
#' is correct for a predictor that is a block of columns. The function
#' therefore returns `FALSE` where the class that owns the method is the base
#' class `structural_term`.
#'
#' @param term A built term.
#'
#' @return A single logical.
#'
#' @seealso [outer_gradient_ok()]
#'
#' @keywords internal
answers_term_third <- function(term) {
  m <- tryCatch(S7::method(modelterms7::term_third, S7::S7_class(term)),
                error = function(e) NULL)
  if (is.null(m)) return(FALSE)
  owner <- attr(m, "signature")[[1]]
  base <- modelterms7::structural_term
  !(identical(attr(owner, "name"), attr(base, "name")) &&
    identical(attr(owner, "package"), attr(base, "package")))
}


#' Does a Term Supply Its Fourth Derivative?
#'
#' @description
#' Determines whether the term implements [modelterms7::term_fourth()], which
#' is the question of [answers_term_third()] one order up, and reads it the
#' same way: from the class on which the method is registered, never from a
#' list of class names.
#'
#' @details
#' The criterion's own second derivative reads a fourth order through the
#' recursion, so a term that implements the third derivative and not the
#' fourth supplies an exact gradient and no exact Hessian. For such a term the
#' function returns `FALSE` and the Hessian is computed by
#' [statmod_hess_stencil()]. A term of the likelihood shape (`regime()` and the
#' marginal break-point term) implements neither derivative, and
#' [outer_gradient_ok()] returns `FALSE` for it at both orders.
#'
#' @param term A built term.
#'
#' @return A single logical.
#'
#' @seealso [outer_gradient_ok()], [statmod_structural_hess()]
#'
#' @keywords internal
answers_term_fourth <- function(term) {
  m <- tryCatch(S7::method(modelterms7::term_fourth, S7::S7_class(term)),
                error = function(e) NULL)
  if (is.null(m)) return(FALSE)
  owner <- attr(m, "signature")[[1]]
  base <- modelterms7::structural_term
  !(identical(attr(owner, "name"), attr(base, "name")) &&
    identical(attr(owner, "package"), attr(base, "package")))
}


#' The Structural Term of a Model, If Any
#'
#' @description
#' Returns the single term of the structural branch, or `NULL` where the model
#' carries none.
#'
#' @details
#' At most one structural term is admitted per formula, so a caller that tests
#' whether the term implements [modelterms7::term_third()] or
#' [modelterms7::term_fourth()] needs one term. The term is found through the
#' design's structural attribute, where the fit records what it built, and not
#' by walking the terms.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#'
#' @return A built term, or `NULL`.
#'
#' @keywords internal
structural_term_of <- function(spec, design) {
  su <- attr(design, "structural")
  if (!length(su)) return(NULL)
  u <- su[[1L]]
  tm <- spec@terms[[u$param]][[u$term]]
  if (is.null(tm)) find_term(spec, u$term) else tm
}


#' Does a Penalty Supply What a Marginal Criterion Needs?
#'
#' @description
#' Evaluates the derivatives of a penalty at a probe value of its
#' hyperparameters and reports whether they could be computed.
#'
#' @details
#' The gradient requires [penalties7::penalty_dhessian()]. A penalty that is
#' not quadratic in the coefficients ([penalties7::beta_quadratic()]) must also
#' supply [penalties7::penalty_dhessian_beta()], because its Hessian moves with
#' the mode. The order-2 route additionally requires
#' [penalties7::penalty_d2hessian()] and [penalties7::penalty_dcross()] and,
#' for a penalty that is not quadratic, [penalties7::penalty_d2hessian_beta()]
#' and [penalties7::penalty_dhessian_beta_theta()]. A penalty without
#' hyperparameters returns `FALSE`.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param order `1` or `2`.
#'
#' @return A single logical.
#'
#' @keywords internal
penalty_answers <- function(pen, order = 1L) {
  th <- as.list(penalty_theta_start(pen))
  k <- as.integer(pen@n_coef)
  if (!length(th) || is.na(k) || k < 1L) return(FALSE)
  b <- rep(0.37, k)
  ok <- tryCatch({
    penalties7::penalty_dhessian(pen, b, th)
    bq <- isTRUE(penalties7::beta_quadratic(pen, th))
    # ⚠️ A PENALTY WHOSE HESSIAN MOVES WITH THE COEFFICIENTS moves the
    # determinant as the mode moves, so the gradient reads dS/dbeta[v] beside
    # the log-likelihood's own third derivative. It used to be admitted here
    # without it: measured on a Student t prior over 30 groups, the gradient was
    # 4.3e-04 out and FLAT in the step, and adding tr(M dS/dbeta[v]) took it to
    # 1.5e-06. A penalty that cannot supply it leaves the search without an
    # exact gradient rather than with one missing that piece.
    if (!bq) penalties7::penalty_dhessian_beta(pen, b, th, rep(0.21, k))
    if (order >= 2L) {
      if (!bq) {
        penalties7::penalty_d2hessian_beta(pen, b, th, rep(0.21, k),
                                           rep(-0.13, k))
        penalties7::penalty_dhessian_beta_theta(pen, b, th, rep(0.21, k))
      }
      penalties7::penalty_d2hessian(pen, b, th)
      penalties7::penalty_dcross(pen, b, th)
    }
    TRUE
  }, error = function(e) FALSE)
  isTRUE(ok)
}


#' The Exact Gradient of the Marginal Criterion
#'
#' @description
#' Computes \eqn{\partial V/\partial\eta} at the penalized mode, over the free
#' scale of the hyperparameters under estimation.
#'
#' @details
#' The gradient is the sum of the envelope term
#' \eqn{-\partial\rho/\partial\theta}, the explicit derivative of the
#' determinant, \eqn{-\tfrac12\mathrm{tr}(M\,\partial S/\partial\theta)}, the
#' implicit one, \eqn{-\tfrac12 u'v}, and, where a penalty's Hessian depends on
#' the coefficients, \eqn{-\tfrac12\mathrm{tr}(M\,\partial S/\partial\beta[v])}.
#' Here \eqn{v = -J^{-1}\,\partial^2\rho/\partial\beta\partial\theta} is how the
#' mode moves, \eqn{u_c = \mathrm{tr}(M\,\partial K/\partial\beta_c)} is how
#' the determinant reads that movement, and \eqn{M} is \eqn{K^{-1}} under
#' [reml()] and the inverse of \eqn{K} on the integrated subspace under
#' [ml()].
#'
#' The matrix \eqn{J} is the curvature of the penalized likelihood at the mode:
#' the observed information plus the penalty's Hessian, plus the second
#' derivative of a block that moves with its coefficients
#' ([mode_curvature()]). It differs from \eqn{K} when the criterion reads the
#' expected information. A coefficient that the criterion estimates, a
#' direction that the model does not identify and a coordinate that a kink of a
#' break-point term holds do not move with \eqn{\theta}, and the function
#' solves on the directions that they leave free ([mode_solver()]). Where a
#' penalty covers the parameters of a structural term, the function delegates
#' to [statmod_structural_grad()].
#'
#' \eqn{u} is assembled without forming any third-derivative array. Writing
#' \eqn{G_{ab,i} = x_{ia}'M_{[a][b]}x_{ib}} for the per-observation diagonal of
#' the block of \eqn{M},
#' \deqn{u_{k} = -X_k'\Big(w \sum_{a,b} \ell'''_{abk}\, G_{ab}\Big),}
#' one crossprod per distribution parameter. The component
#' \eqn{\ell'''_{abk}} is looked up by a name built from the parameter names in
#' the family's own order, never parsed out of one.
#'
#' For each coefficient in `gam`, which the criterion estimates, the entry is
#' the score of the log-likelihood in that coefficient minus half of the
#' movement of the determinant along the direction in which the mode moves
#' when the coefficient changes by one unit.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param idx The outer index.
#' @param basis The integrated subspace, or `NULL`.
#' @param free Whether to carry the result onto the free scale. The Hessian
#'   assembly calls this function with `free = FALSE` and applies its own
#'   second-order chain rule.
#' @param ctx The evaluation context, or `NULL`.
#' @param gam The coefficients that the criterion estimates, from
#'   [marginal_coords()], or `NULL`.
#'
#' @return A numeric vector, one entry per row of `idx` followed by one per
#'   position of `gam$where`, or `NULL` where the determinant does not exist.
#'
#' @seealso [statmod_marginal()], [reml()]
#'
#' @keywords internal
statmod_marginal_grad <- function(spec, design, coef, hyper, method, idx,
                                  basis = NULL, free = TRUE, ctx = NULL,
                                  gam = NULL) {
  if (structural_penalized(spec, design)) {
    return(statmod_structural_grad(spec, design, coef, hyper, method, idx,
                                   basis, free))
  }
  # ⚠️ THE BLOCK AT THE MODE, not the block as it was built.
  # statmod_information_at() refreshes internally, so K -- and therefore M --
  # is assembled on the refreshed block; block_leverage() and u_vector()'s
  # final crossprod were handed this design as it arrived, which for a term
  # registering term_refresh() is the block at the coefficients the fit
  # STARTED from. The two agree for every fixed design, which is why nothing
  # showed until a penalty was put inside nl() or seg(): measured on
  # nl(a ~ 0 + ridge(~grp)), the two blocks differ by 2.07 and u is wrong by
  # 110 -- on the SIGMA row among others, an equation carrying no refreshable
  # term at all, because G_ab reads the mu block whatever row is being formed.
  # Refreshed, that row is -21.902535 against a direct numerical
  # tr(M dK/dbeta) of -21.902535.
  # the design as it arrived carries the structural state, which the joint
  # movement below reads; the refreshed one differs from it only in its blocks
  design0 <- design
  design <- statmod_design_at(spec, coef, design)
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)

  # the information, the penalty, K = H + S and its inverse come from the
  # context when there is one, so a criterion, a gradient and a Hessian read
  # at the same point assemble them once between them
  expected <- info_expected(method@hessian)
  pen <- ctx_penalized(ctx, spec, design, coef, hyper, expected,
                       laplace = TRUE)
  if (is.null(pen)) return(NULL)
  # ⚠️ TWO different matrices, and they coincide only on the observed route.
  # The determinant is of the CRITERION's K, so M and u are read off that one;
  # but the mode is where the penalized LIKELIHOOD's score vanishes, so how it
  # moves is governed by H_obs + S whatever the criterion's own matrix is --
  # which is what the derivation at the top of this file already writes down.
  # Using one for both was invisible while the criterion was always the
  # observed one, and on the expected route it is a systematic error that
  # shrinks with n exactly as H_obs approaches H_exp: measured on a gamma
  # smooth at 300, 1000 and 3000 observations, 1.9e-03, 1.4e-03 and 1.1e-04
  # against a finite difference of the criterion, and FLAT in the inner
  # tolerance, which is what said it was not the mode's location.
  # and a third reason for two: the criterion leaves a kinked penalty's
  # coordinates out of its determinant, while the mode moves in all of them
  mode_pen <- ctx_penalized(ctx, spec, design, coef, hyper, FALSE)
  if (is.null(mode_pen)) return(NULL)
  # how the mode moves: the penalized likelihood's own curvature, which for a
  # refreshable block is not the Gauss-Newton matrix the design gives
  Dm <- mode_curvature(spec, design, coef, params, npar, offs, total)
  # a coordinate pinned in K other than at a boundary (a coefficient the
  # criterion estimates, one the model does not identify), and one that a
  # break-point term's kink holds at the mode, stays where it is when the
  # mode moves; see mode_solver()
  ms <- mode_solver(spec, design, coef, mode_pen$K, mode_pen$inv, Dm, npar,
                    offs)
  Dm <- ms$Dm
  msolve <- ms$solve
  M <- ctx_trace_matrix(ctx, pen, basis, expected)
  if (is.null(M)) return(NULL)
  # ⚠️ THE CRITERION PINS A SATURATED COORDINATE, SO ITS GRADIENT MUST TOO.
  # pin_boundary() replaces the row and column of a coordinate whose curvature
  # is not finite by the identity's, so in the criterion that block is a
  # CONSTANT and differentiates to zero. ctx_trace_matrix() has already zeroed
  # M there; what is left is that the coordinate's own entry of u is zero too.
  # Read through the identity instead, the trace multiplied the family's third
  # derivative at a parameter on its link's clamp -- NaN there, a Student t's
  # nu at double.xmax -- by 1, and the whole gradient came back NaN. Measured
  # on the reference battery's fam-studentt, every route then stopped after
  # ONE evaluation at the starting hyperparameters (criterion -1550.286),
  # where with the gradient of the pinned criterion -- which a central
  # difference of that criterion confirms, -18.09516 against -18.09511 and
  # -18.09516 at h of 3e-3 and 1e-3 -- all three reach -1544.442 and certify.
  held <- union(attr(pen$K, "boundary_at"), attr(mode_pen$K, "boundary_at"))

  # dK/dbeta is the third derivative of the log-likelihood on the observed
  # route and the derivative of the expected information on the expected one.
  # The contraction is the same in both -- one crossprod per distribution
  # parameter against the same leverage diagonal -- so only the array and the
  # key it is read by change.
  km <- ctx_kmove(ctx, spec, design, coef, hyper, method)
  G <- ctx_leverage(ctx, design, M, params, npar, offs, spec@threads)

  # ⚠️ A FILTER THE DETERMINANT DOES NOT SPAN STILL MOVES WITH THE MODE.
  # With no penalty over its own parameters the criterion's K is over the
  # coefficients alone, on the static design, but the predictors it is read
  # at carry the filter's level, which depends on the coefficients through
  # the recursion and on the filter's own parameters -- and those are
  # estimated beside the coefficients and move with the hyperparameter too.
  # So the mode moves JOINTLY, by the joint penalized curvature, and K moves
  # along the TOTAL derivative of every predictor, which is u_vector()'s
  # formula with the joint rows in place of the design. No term_third() is
  # needed: K is on the static design, so the recursion's own second
  # derivative (the E_t of the penalized route) is not in it.
  #
  # Before this the gradient read the coefficients alone for both, and was
  # out by three pieces measured on a smooth beside an unpenalized gas(1, 1)
  # at 400 observations, against a central difference of the criterion with
  # the mode refitted: -8.49e-03 from the mode's movement, +8.90e-04 from the
  # determinant read along X v instead of along the total derivative, and
  # -1.10e-05 from the filter's own path, summing to -7.611e-03 against a gap
  # of -7.610e-03 -- 7.2e-03 relative and FLAT in the step. With both read
  # jointly the gap is 5.6e-05, 5.1e-06 and 1.5e-06 at h of 1e-2, 3e-3 and
  # 1e-3, which is O(h^2) onto the reference's own floor.
  joint <- NULL
  if (any(vapply(attr(design0, "structural"),
                 function(su) identical(su$kind, "filter"), logical(1)))) {
    joint <- filter_joint_movement(
      spec, design0, coef, ctx_penalty(ctx, spec, design, coef, hyper), Dm,
      total, ctx_deriv(ctx, spec, design, coef, hyper, 3L),
      hold = ms$hold, Ak = ms$Ak)
    if (is.null(joint)) return(NULL)
    msolve <- joint$solve
  }
  u <- if (is.null(joint)) {
    u_vector(spec, design, coef, M, params, npar, offs, total,
             d3 = km$deriv, key = km$key, G = G) +
      u_refresh(spec, design, coef, M, params, npar, offs, total, expected,
                ctx_emethod(ctx))
  } else {
    uj <- u_vector(spec, design, coef, M, params, npar, offs, total,
                   d3 = km$deriv, key = km$key, G = G, rows = joint$rows)
    ib <- seq_len(total)
    uj[ib] <- uj[ib] + u_refresh(spec, design, coef, M, params, npar, offs,
                                 total, expected, ctx_emethod(ctx))
    uj
  }
  if (length(held)) u[held] <- 0

  out <- numeric(nrow(idx))
  links <- attr(idx, "links")
  # over the MEMBERS rather than the rows: a shared hyperparameter is one
  # coordinate standing for several penalties, and the criterion's derivative
  # in it is the sum of theirs. With nothing shared there is one member per
  # row and this is the loop that was here before.
  mem <- index_members(idx)
  for (a in seq_along(params)) {
    p <- params[a]
    rows <- which(mem$parameter == p)
    if (!length(rows)) next
    for (nm in unique(mem$term[rows])) {
      # not `u`: this function already uses that name for the contraction of
      # the third derivative, and shadowing it fails several frames down
      un <- statmod_unit(spec, design, p, nm)
      pos <- un$index
      pen <- un$penalty
      # a covariance class spans more than one equation, so its
      # coefficients come from the stacked vector; for an ordinary unit
      # this is what coef[[p]][un$cols] returned
      bt <- unit_beta(un, coef, params)
      th <- as.list(hyper[[p]][[nm]])
      gt <- penalties7::penalty_grad_theta(pen, bt, th)
      cr <- penalties7::penalty_cross(pen, bt, th)
      dS <- penalties7::penalty_dhessian(pen, bt, th)
      for (i in rows[mem$term[rows] == nm]) {
        r <- mem$row[i]
        h <- mem$name[i]
        # the mode moves by -(H+S)^-1 d2rho/dbeta dtheta, with H the penalized
        # LIKELIHOOD's own curvature -- see mode_curvature()
        c_m <- numeric(total)
        c_m[pos] <- as.numeric(cr[[h]])
        v <- -as.numeric(msolve(c_m))
        dS_m <- matrix(0, total, total)
        dS_m[pos, pos] <- as_dense(dS[[h]])
        # the mode's movement reaches the determinant through every penalty
        # whose Hessian moves with the coefficients, not only this one's: see
        # statmod_penalty_dbeta()
        trp <- penalty_dbeta_trace(spec, design, coef, hyper, v, M)
        dtheta <- -as.numeric(gt[[h]]) -
          (sum(M * dS_m) + sum(u * v) + trp) / 2
        # and onto the free scale the search runs on
        out[r] <- out[r] + if (!free) dtheta else {
          eta <- linkfunctions7::linkfun(links[[r]], hyper[[p]][[nm]][[h]])
          dtheta * linkfunctions7::dlinkinv(links[[r]], eta)
        }
      }
    }
  }
  if (is.null(gam) || !length(gam$where)) return(out)
  # THE COEFFICIENTS THE CRITERION ESTIMATES. The inner fit holds each one,
  # so the mode moves with it by b = e_j - P K_m e_j, P the inverse of the
  # mode's curvature over the free coordinates -- which is what the pinned
  # inverse behind msolve() is. The penalty does not depend on such a
  # coefficient, its row of S being zero, so the envelope leaves the
  # log-likelihood's own score, and the determinant moves along b.
  Hm <- as_dense(ctx_information(ctx, spec, design, coef, hyper, FALSE,
                                 ctx_emethod(ctx))) +
    as_dense(ctx_penalty(ctx, spec, design, coef, hyper)) + Dm
  sc <- statmod_score_at(spec, coef, design, index = gam$where)
  gout <- numeric(length(gam$where))
  # every coordinate the mode's curvature pins other than a boundary: the
  # estimated coefficients, the ones the model does not identify and the ones
  # a break-point's kink holds
  hold <- ms$hold
  for (i in seq_along(gam$where)) {
    w <- gam$where[i]
    col <- Hm[, w]
    col[hold] <- 0
    v <- -as.numeric(msolve(col))
    v[hold] <- 0
    b <- v
    b[w] <- 1
    trp <- penalty_dbeta_trace(spec, design, coef, hyper, b, M)
    gout[i] <- sc[i] - (sum(u * b) + trp) / 2
  }
  c(out, gout)
}


#' How the Mode Moves Where an Unpenalized Filter Moves With It
#'
#' @description
#' Computes the pieces that [statmod_marginal_grad()] reads on the joint vector
#' of coefficients and a filter's own parameters when no penalty covers those
#' parameters: a solve with the joint penalized curvature, which gives the
#' movement of the mode, and the derivative of every equation's predictor over
#' that vector, along which the movement of the determinant is read.
#'
#' @details
#' The mode is where the penalized objective's gradient in the coefficients
#' AND in the filter's parameters vanishes, so differentiating that condition
#' in a hyperparameter gives
#' \deqn{(J + S)\,v = -\partial^2\rho/\partial u\,\partial\theta,}
#' with \eqn{J} the joint observed information of
#' [statmod_full_information()] and \eqn{S} the penalty's Hessian placed on the
#' coefficients, the only coordinates that it covers. A block that moves with
#' its coefficients adds its own second derivative on those coordinates, as it
#' does on the coefficient route.
#'
#' The rows are the static design of every equation, except that the equation
#' carrying the filter is replaced by the forward Jacobian of the recursion,
#' which is read from [filter_curvature()] at the same memo slot that the joint
#' information fills, so the recursion runs once at a point. A held level is
#' dropped from both, exactly as the information drops it.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design as it arrived, carrying the structural state.
#' @param coef The coefficients at the penalized mode.
#' @param S The penalty's Hessian over the stacked coefficients.
#' @param Dm The mode's own correction for a block that moves, from
#'   [mode_curvature()].
#' @param total The number of stacked coefficients.
#' @param D3 The family's third derivative on the link scale at the fitted
#'   predictors, which the recursion reads where the memo misses.
#' @param hold The positions of the coefficients that [mode_solver()] holds
#'   (a break-point held at the minimum of its profile), which stay where
#'   they are when the mode moves. The default is `integer(0)`.
#' @param Ak The kink constraints over the stacked coefficients, from
#'   [mode_solver()], or `NULL` for none.
#'
#' @return A list with `rows`, one matrix per distribution parameter over the
#'   estimated coordinates, and `solve`, a function of a vector over the
#'   coefficients returning the joint solve over every estimated coordinate;
#'   or `NULL` where the joint curvature cannot be formed or solved.
#'
#' @seealso [statmod_marginal_grad()], [filter_curvature()],
#'   [statmod_full_information()]
#'
#' @keywords internal
filter_joint_movement <- function(spec, design, coef, S, Dm, total, D3,
                                  hold = integer(0), Ak = NULL) {
  jd <- joint_design_rows(spec, design, coef)
  if (is.null(jd)) return(NULL)
  J <- tryCatch(as_dense(statmod_full_information(spec, coef, design)),
                error = function(e) NULL)
  width <- length(jd$keep)
  if (is.null(J) || nrow(J) != width || !all(is.finite(J))) return(NULL)
  ib <- seq_len(total)
  J[ib, ib] <- J[ib, ib] + as_dense(S) + Dm
  # ⚠️ THE COORDINATES THAT mode_solver() HOLDS STAY HELD ON THE JOINT VECTOR.
  # A break-point that a jump() term holds at the minimum of its profile has
  # a zero row and column in the joint curvature, and a kink is a linear
  # constraint on the coefficients; without either the joint matrix is
  # singular and the gradient was lost. Measured on pig2 with jump(tt) +
  # gas(1, 1) | alpha ~ s(tt): the break-point's row of J is exactly zero,
  # the smallest eigenvalue -8.6e-13, and the certificate read "the outer
  # gradient is not finite". The movement is solved on the directions that
  # the held coordinates and the kink rows leave free, the filter's own
  # parameters among them, as mode_solver() does on the coefficients.
  B <- matrix(0, length(hold), width)
  if (length(hold)) B[cbind(seq_along(hold), hold)] <- 1
  if (!is.null(Ak) && nrow(Ak) > 0L) {
    B <- rbind(cbind(as.matrix(Ak), matrix(0, nrow(Ak), width - total)), B)
  }
  N <- null_basis(B, width)
  R <- crossprod(N, J %*% N)
  fac <- tryCatch(chol(R), error = function(e) NULL)
  rsolve <- if (!is.null(fac)) {
    function(r) backsolve(fac, forwardsolve(t(fac), r))
  } else {
    # the true joint curvature may lose definiteness where it is still
    # invertible, and a mode's movement is worth having from it: the implicit
    # function theorem asks for a nonsingular matrix, not a definite one
    Ri <- tryCatch(solve(R), error = function(e) NULL)
    if (is.null(Ri)) return(NULL)
    function(r) Ri %*% r
  }
  solver <- function(z) N %*% rsolve(crossprod(N, as.numeric(z)))
  d <- spec@distrib
  th <- jd$ev$theta
  gl <- distributions7::distrib_gradient(d, spec@response, th, scale = "link",
                                         threads = spec@threads)
  H <- distributions7::distrib_hessian(d, spec@response, th, scale = "link",
                                       threads = spec@threads)
  V <- jd$V
  if (isTRUE(jd$f$joint)) {
    cv <- joint_filter_curvature(spec, design, jd$f, jd$V, gl, H, D3, th)
    GF <- match(jd$f$params, jd$params)
    for (k in seq_along(GF)) V[[GF[k]]] <- cv$jacobian[[k]]
  } else {
    cv <- filter_curvature(spec, design, jd$f, jd$ap, jd$V, gl, H, D3, th)
    V[[jd$ap]] <- cv$jacobian
  }
  list(rows = lapply(V, function(x) x[, jd$keep, drop = FALSE]),
       solve = function(z) {
         zz <- numeric(width)
         zz[ib] <- z
         as.numeric(solver(zz))
       })
}


#' How a Penalty's Hessian Moves With the Mode
#'
#' @description
#' Computes \eqn{\sum_u \partial S_u/\partial\beta\,[v]} over every penalty on
#' the stacked coefficients whose Hessian depends on them, placed in the
#' stacked coefficient space. `penalty_dbeta_trace()` computes its trace
#' against a matrix, and `penalty_dbeta_blocks()` returns the per-penalty
#' pieces from which both are assembled.
#'
#' @details
#' A marginal criterion's determinant is of \eqn{H + S}, so the mode's movement
#' \eqn{v} reaches it through \eqn{S} wherever \eqn{S} depends on the
#' coefficients, as it does for a heavy-tailed prior on a random effect, and the
#' piece is \eqn{\mathrm{tr}(M\,\partial S/\partial\beta[v])}. A
#' prediction-error criterion reads the matrix itself inside its trace. Every
#' such penalty contributes, not only the one owning the hyperparameter being
#' differentiated, because the mode moves in all of the coefficients at once.
#'
#' A penalty for which [penalties7::beta_quadratic()] returns `TRUE` is
#' skipped (this is every ridge, smooth and Gaussian random effect), so a
#' model carrying only those assembles nothing: the matrix is `NULL` and
#' the trace exactly `0`, which leaves the gradient unchanged. A kinked penalty
#' and a mixed covariance class are skipped as well, and a penalty over a
#' structural term's own parameters is read on the joint vector by
#' [statmod_structural_grad()].
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param v The direction, over the stacked coefficients.
#' @param total The stacked width.
#' @param M The matrix against which the trace is taken, `total` square.
#'
#' @return `statmod_penalty_dbeta()` a `total` square matrix, or `NULL` if
#'   no penalty's Hessian moves; `penalty_dbeta_trace()` a single number;
#'   `penalty_dbeta_blocks()` a list of `pos` and `T` pairs.
#'
#' @seealso [statmod_marginal_grad()], [statmod_pe_derivs()],
#'   [penalties7::penalty_dhessian_beta()]
#'
#' @keywords internal
statmod_penalty_dbeta <- function(spec, design, coef, hyper, v, total) {
  bl <- penalty_dbeta_blocks(spec, design, coef, hyper, v)
  if (!length(bl)) return(NULL)
  out <- matrix(0, total, total)
  for (b in bl) out[b$pos, b$pos] <- out[b$pos, b$pos] + b$T
  out
}

#' @rdname statmod_penalty_dbeta
#' @keywords internal
penalty_dbeta_trace <- function(spec, design, coef, hyper, v, M) {
  s <- 0
  for (b in penalty_dbeta_blocks(spec, design, coef, hyper, v)) {
    s <- s + sum(as_dense(M[b$pos, b$pos, drop = FALSE]) * b$T)
  }
  s
}

#' @rdname statmod_penalty_dbeta
#' @keywords internal
penalty_dbeta_blocks <- function(spec, design, coef, hyper, v) {
  params <- spec@distrib@params
  out <- list()
  for (u in statmod_penalized(spec, design)) {
    if (isTRUE(u$structural) || isTRUE(u$mixed) || is.null(u$penalty) ||
        !length(u$index)) next
    pen <- u$penalty
    th <- as.list(hyper[[u$param]][[u$key]])
    if (isTRUE(penalties7::beta_quadratic(pen, th))) next
    if (length(penalties7::penalty_kinks(pen, th))) next
    pos <- u$index
    bt <- unit_beta(u, coef, params)
    out[[length(out) + 1L]] <- list(pos = pos, T = as_dense(
      penalties7::penalty_dhessian_beta(pen, bt, th, v[pos])))
  }
  out
}


#' How a Moving Penalty Hessian Enters the Criterion's Second Derivative
#'
#' @description
#' Computes the part of \eqn{\partial K_m/\partial t_l} that a penalty whose
#' Hessian depends on the coefficients contributes, placed in the stacked
#' coefficient space, where \eqn{K_m = S_m + T[b_m]} is the determinant's
#' matrix differentiated in hyperparameter \eqn{m}.
#'
#' @details
#' Differentiating \eqn{S_m(\beta, t) + \partial S/\partial\beta\,[b_m]} once
#' more along \eqn{t_l}, while \eqn{\beta} moves by \eqn{b_l}, gives
#' \deqn{\frac{\partial S_m}{\partial\beta}[b_l]
#'   + \frac{\partial S_l}{\partial\beta}[b_m]
#'   + \frac{\partial^2 S}{\partial\beta^2}[b_l, b_m]
#'   + \frac{\partial S}{\partial\beta}[b_{ml}]}
#' beside \eqn{S_{ml}}, which [outer_pieces()] already carries. The first two
#' are the mixed derivative \eqn{\partial^2 S/\partial\beta\,\partial\theta}
#' taken along two directions, \eqn{b_l} for the hyperparameters of row \eqn{m}
#' and \eqn{b_m} for those of row \eqn{l}; each is read from
#' [penalties7::penalty_dhessian_beta_theta()] for the hyperparameters that the
#' penalty owns. The third is [penalties7::penalty_d2hessian_beta()] and the
#' fourth [penalties7::penalty_dhessian_beta()]. A penalty for which
#' [penalties7::beta_quadratic()] returns `TRUE` does not contribute, and a
#' kinked penalty, a mixed covariance class and a penalty over a structural
#' term's own parameters do not contribute either, so for a model carrying
#' only those the function returns `NULL` and the Hessian is left untouched.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param idx The outer index.
#' @param m,l The two rows of the index.
#' @param bm,bl,bml The mode's movement along each and along the pair.
#' @param total The stacked width.
#'
#' @return A `total` square matrix, or `NULL` if no penalty's Hessian
#'   moves with the coefficients.
#'
#' @seealso [statmod_penalty_dbeta()], [statmod_marginal_hess()],
#'   [statmod_pe_derivs()]
#'
#' @keywords internal
statmod_penalty_second <- function(spec, design, coef, hyper, idx, m, l, bm,
                                   bl, bml, total) {
  params <- spec@distrib@params
  mem <- index_members(idx)
  out <- NULL
  for (u in statmod_penalized(spec, design)) {
    if (isTRUE(u$structural) || isTRUE(u$mixed) || is.null(u$penalty) ||
        !length(u$index)) next
    pen <- u$penalty
    th <- as.list(hyper[[u$param]][[u$key]])
    if (isTRUE(penalties7::beta_quadratic(pen, th))) next
    if (length(penalties7::penalty_kinks(pen, th))) next
    pos <- u$index
    bt <- unit_beta(u, coef, params)
    E <- as_dense(penalties7::penalty_d2hessian_beta(pen, bt, th, bl[pos],
                                                     bm[pos])) +
      as_dense(penalties7::penalty_dhessian_beta(pen, bt, th, bml[pos]))
    # the hyperparameters of row r this penalty owns, one member line apiece;
    # a row owned by another penalty moves this one's Hessian by nothing
    owned <- function(r, v) {
      nms <- mem$name[mem$row == r & mem$parameter == u$param &
                        mem$term == u$key]
      if (!length(nms)) return(0)
      D <- penalties7::penalty_dhessian_beta_theta(pen, bt, th, v[pos])
      Reduce(`+`, lapply(nms, function(nm) as_dense(D[[nm]])))
    }
    E <- E + owned(m, bl) + owned(l, bm)
    if (is.null(out)) out <- matrix(0, total, total)
    out[pos, pos] <- out[pos, pos] + E
  }
  out
}


#' The Trace of the Determinant's Movement With the Coefficients
#'
#' @description
#' Computes \eqn{u_c = \mathrm{tr}(M\,\partial K/\partial\beta_c)}, assembling
#' one crossprod per distribution parameter.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param M The matrix against which the trace is taken.
#' @param params The distribution's parameter names.
#' @param npar,offs,total The block sizes, their offsets and the total.
#' @param d3 The array from which \eqn{\partial K/\partial\eta} is read, with
#'   a negative sign: the third derivatives of the log-likelihood on the
#'   observed route, the derivatives in the predictors of \eqn{E[\ell'']} on
#'   the expected one. The value `NULL` selects the family's third derivative
#'   at the fitted predictors.
#' @param G The per-observation diagonals of \eqn{M}, or `NULL`.
#' @param key The builder of `d3`'s component names, or `NULL`.
#' @param rows `NULL` for the contraction against each equation's design, or
#'   one matrix per distribution parameter giving the derivative of that
#'   equation's predictor in a wider vector, as [filter_joint_movement()]
#'   returns them. The result is then as long as their common width.
#'
#' @return A numeric vector as long as the stacked coefficients, or as the
#'   width of `rows`.
#'
#' @keywords internal
u_vector <- function(spec, design, coef, M, params, npar, offs, total,
                     d3 = NULL, G = NULL, key = NULL, rows = NULL) {
  n <- spec@n_obs
  if (is.null(d3)) {
    th <- statmod_eta(spec, design, coef)$theta
    d3 <- distributions7::distrib_deriv3(spec@distrib, spec@response, th,
                                         scale = "link", threads = spec@threads)
  }
  keys <- names(d3)
  # the observed route's array is symmetric in all three indices and keyed by
  # the sorted triple; the expected route's is symmetric in (a, b) only, the
  # measure's own derivative not being, so it carries its own builder
  if (is.null(key)) key <- function(a, b, k) d3_key(params, a, b, k, keys)
  if (is.null(G)) G <- block_leverage(design, M, params, npar, offs,
                                      spec@threads)

  out <- numeric(if (is.null(rows)) total else ncol(rows[[1L]]))
  for (k in seq_along(params)) {
    # an equation with no coefficients has no predictor to move where the
    # operand is the design; where it is the joint rows, the equation
    # carrying a filter moves with the filter's own parameters whatever its
    # design holds
    if (npar[k] == 0L && is.null(rows)) next
    s <- numeric(n)
    for (a in seq_along(params)) {
      if (npar[a] == 0L) next
      for (b in seq_along(params)) {
        if (npar[b] == 0L) next
        # a block of the leverage that is exactly zero contributes nothing,
        # and must not be multiplied by a component that is NaN: that is
        # what a coordinate held at its link's clamp leaves in the family's
        # third derivative, and 0 * NaN is NaN
        if (!any(G[[a]][[b]] != 0)) next
        s <- s + rep_len(d3[[key(a, b, k)]], n) * G[[a]][[b]]
      }
    }
    if (is.null(rows)) {
      out[offs[k] + seq_len(npar[k])] <-
        -as.numeric(crossprod(design[[params[k]]]$X, spec@weights * s))
    } else {
      out <- out - as.numeric(crossprod(rows[[k]], spec@weights * s))
    }
  }
  out
}


#' How the Determinant Reads a Block That Moves With the Coefficients
#'
#' @description
#' Computes the part of \eqn{u_c = \mathrm{tr}(M\,\partial K/\partial\beta_c)}
#' that [u_vector()] does not compute: everything that comes from
#' \eqn{\partial X/\partial\beta} where a term's block depends on its own
#' coefficients.
#'
#' @details
#' With \eqn{H_{(a,j),(b,k)} = -\sum_i w_i \ell_{ab,i}X_a[i,j]X_b[i,k]} and a
#' block that moves, differentiating in \eqn{\beta_c} gives three
#' contributions and [u_vector()] computes one. The other two are
#' transposes under the trace, so with
#' \eqn{R_{ab}[i,j] = \sum_k M_{(a,j),(b,k)}X_b[i,k]} and
#' \eqn{A_a[i,j] = w_i\sum_b \ell_{ab,i}R_{ab}[i,j]},
#' \deqn{\Delta u_c = -2\sum_{i,j}A_a[i,j]\,\partial X_a[i,j]/\partial\beta_c.}
#'
#' The derivative is read from the term, through
#' [modelterms7::term_block_contract()], and is never differenced. A term
#' carries its own chain rule, the links on its parameters and a subformula's
#' design. A break-point column is a step function of its break-point, so a
#' difference quotient of it diverges as the step shrinks. A term that does not
#' implement the contraction inherits zeros, which is correct for a fixed
#' design. A sharp `jump()` or `jseg()` returns zeros as well, because its
#' block is a working linearization with a frozen weight.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, already refreshed at `coef`.
#' @param coef The coefficients at the penalized mode.
#' @param M The matrix against which the trace is taken.
#' @param params,npar,offs,total The block bookkeeping.
#' @param expected Whether the criterion carries the expected information.
#' @param emethod The method that obtains the expected information: `NULL`
#'   (the closed form of the family where it has one, and
#'   [distributions7::expect_opg()] otherwise) or an object that inherits from
#'   [distributions7::ExpectMethod]. It is read only where the expected
#'   information is used.
#' @param units The refreshable terms, from [refresh_units()], or
#'   `NULL` to resolve them here.
#' @param Hl The link-scale curvature, from [refresh_hessian()], or
#'   `NULL` to compute it here.
#'
#' @return A numeric vector as long as the stacked coefficients.
#'
#' @seealso [u_vector()], [refresh_amat()],
#'   [modelterms7::term_block_contract()]
#'
#' @keywords internal
u_refresh <- function(spec, design, coef, M, params, npar, offs, total,
                      expected = FALSE, emethod = NULL, units = NULL,
                      Hl = NULL) {
  out <- numeric(total)
  if (is.null(units)) {
    units <- refresh_units(spec, design, coef, params, npar, offs)
  }
  if (!length(units)) return(out)
  if (is.null(Hl)) Hl <- refresh_hessian(spec, design, coef, expected, emethod)
  for (un in units) {
    cw <- vector("list", length(params))
    for (b in seq_along(params)) {
      if (npar[b] == 0L) next
      cw[[b]] <- Hl[[hess_key(params, un$a, b)]]
    }
    A <- refresh_amat(spec, design, M, params, npar, offs, un$ra, cw)
    # NOT wrapped in a tryCatch: a term that has not written the contraction
    # inherits the base method and gets zeros, which is a legitimate answer and
    # not an error, so anything raised here is a defect and must be seen. A
    # catch-all put around it swallowed "not an exported object" and reported a
    # correction of exactly zero, which reads as "nothing to correct".
    dc <- modelterms7::term_block_contract(un$built, coef = un$bt, A = A)
    if (length(dc) != length(un$cols)) {
      stop(sprintf(paste("term_block_contract() returned %d values for a block",
                         "of %d columns\n  in '%s'."),
                   length(dc), length(un$cols), un$term), call. = FALSE)
    }
    out[un$ra] <- -2 * as.numeric(dc)
  }
  out
}


#' The Curvature by Which the Mode Actually Moves
#'
#' @description
#' Computes the difference between the true Hessian of the penalized
#' log-likelihood and the Gauss-Newton matrix that [statmod_information_at()]
#' returns, where a term's block moves with its coefficients.
#'
#' @details
#' \eqn{v = \partial\hat\beta/\partial\theta} solves
#' \eqn{(\partial^2\rho/\partial\beta^2 - \partial^2\ell/\partial\beta^2)v =
#' -\partial^2\rho/\partial\beta\partial\theta}, and
#' \eqn{\partial^2\ell/\partial\beta^2} is the true second derivative:
#' \deqn{\sum_i\sum_{a,b}\ell_{ab}\frac{\partial\eta_a}{\partial\beta}
#'   \frac{\partial\eta_b}{\partial\beta}
#'   + \sum_i\sum_a \ell_a\frac{\partial^2\eta_a}{\partial\beta^2}.}
#' The design gives the first sum; the second is zero for every fixed
#' block and is not for a refreshable one, where \eqn{X} is the Jacobian and so
#' \eqn{\partial^2\eta_i/\partial\beta_c\partial\beta_d =
#' \partial X_{id}/\partial\beta_c}. It is supported on the term's own block,
#' which keeps the computation cheap: one call of
#' [modelterms7::term_block_contract()] per column, weighted by the
#' score where [u_refresh()] weights by \eqn{M}.
#'
#' This is the matrix by which the mode moves, and it differs from the
#' criterion's matrix. The determinant is of the matrix that
#' [statmod_information_at()] assembles, and the derivative of the determinant
#' is read from that matrix; how the mode moves depends on the penalized
#' likelihood and is read from this one. The two coincide for the observed
#' information when no block moves with its coefficients. With the expected
#' information a single matrix for both would give a systematic error in the
#' gradient, and the two are kept apart.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, already refreshed at `coef`.
#' @param coef The coefficients at the penalized mode.
#' @param params,npar,offs,total The block bookkeeping.
#' @param gl The per-observation weights of each equation's second
#'   derivative, a named list keyed by parameter, or `NULL` (the default) for
#'   the family's score on the link scale at `coef`. The Jeffreys engine
#'   ([jfp_hessian()]) passes its own weights here.
#'
#' @return A symmetric matrix over the stacked coefficients, with nonzero
#'   entries only on the blocks of refreshable terms.
#'
#' @seealso [u_refresh()],
#'   [modelterms7::term_block_contract()]
#'
#' @keywords internal
mode_curvature <- function(spec, design, coef, params, npar, offs, total,
                           gl = NULL) {
  D <- matrix(0, total, total)
  rf <- attr(design, "refresh")
  st <- attr(design, "state")
  if (is.null(rf) || !length(rf) || is.null(st)) return(D)
  n <- spec@n_obs
  w <- spec@weights
  if (is.null(gl)) {
    th <- statmod_eta(spec, design, coef)$theta
    gl <- distributions7::distrib_gradient(spec@distrib, spec@response, th,
                                           scale = "link",
                                           threads = spec@threads)
  }
  for (r in rf) {
    p <- r$param
    a <- match(p, params)
    if (is.na(a) || npar[a] == 0L) next
    cols <- design[[p]]$blocks[[r$term]]
    if (!length(cols)) next
    ra <- offs[a] + cols
    bt <- coef[[p]][cols]
    tm <- modelterms7::term_refresh(st$terms[[p]][[r$term]], bt)
    sc <- w * rep_len(gl[[p]], n)
    for (dd in seq_along(cols)) {
      A <- matrix(0, n, length(cols))
      A[, dd] <- sc
      D[ra, ra[dd]] <- -as.numeric(
        modelterms7::term_block_contract(tm, coef = bt, A = A))
    }
  }
  # symmetric in (c, d) as a second derivative, and symmetrized explicitly:
  # the two orders collect the same terms in a different one, and this file
  # already records that "symmetric by construction" is not a property
  # floating-point addition has
  (D + t(D)) / 2
}


#' The Name of a Third-Derivative Component
#'
#' @description
#' Locates the \eqn{(a, b, k)} entry of a distribution's third derivative,
#' which is keyed by name, never by position.
#'
#' @details
#' The name is built by putting the three parameter names in the family's order
#' and joining them with an underscore, and is then checked against the names
#' that the derivative returned.
#'
#' @param params The parameter names, in the family's order.
#' @param a,b,k Indices into `params`.
#' @param keys The names the derivative actually returned.
#'
#' @return A single string.
#'
#' @keywords internal
d3_key <- function(params, a, b, k, keys) {
  want <- paste(params[sort(c(a, b, k))], collapse = "_")
  if (!want %in% keys) {
    stop(sprintf("The third derivative has no component '%s'.", want),
         call. = FALSE)
  }
  want
}


#' The Rows of the Joint Predictor Derivative
#'
#' @description
#' Computes \eqn{V_{a,t} = \partial\eta_{a,t}/\partial u} over the coefficients
#' of every equation followed by a structural term's own parameters. The
#' function returns each equation's static design placed in its own columns,
#' with zeros in the columns of the structural term's parameters. The callers
#' replace the matrix of the equation that carries the filter by the forward
#' Jacobian of the recursion, which is dense.
#'
#' @details
#' It exists so that the contraction that [u_vector()] carries out is
#' written once. The formula there does not change where a filter is present;
#' only its operand does: \eqn{X} becomes \eqn{[X_p \mid D]}. Everything
#' downstream, the third derivative against the diagonal of \eqn{M} and the
#' movement of the mode, then reads the joint vector with no special case.
#'
#' The rows are returned at the full width, including a held level, and the
#' caller drops that level exactly as [statmod_full_information()] does.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#'
#' @return `NULL` where the model has no filter. Otherwise it is a list with `V`
#'   (one matrix per distribution parameter), `ap` (which equation carries the
#'   filter), `nb`, `np` and `mf` (the numbers of coefficients, of the term's
#'   parameters and the joint width), `keep` (the columns that survive a held
#'   level), `zn` (the term's parameter names), `ev`, `f`, `params` and `n`.
#'
#' @seealso [u_vector()], [statmod_full_information()]
#'
#' @keywords internal
joint_design_rows <- function(spec, design, coef) {
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  nb <- sum(npar)
  n <- spec@n_obs
  ev <- statmod_eta(spec, design, coef)
  if (!length(ev$filters)) return(NULL)
  f <- ev$filters[[1L]]
  ap <- match(f$param, params)
  zn <- modelterms7::term_params(f$tm)
  np <- length(zn)
  mf <- nb + np
  # as_dense() for the reason statmod_full_information() records: an equation
  # carrying a random effect has a sparse design, and writing a dgCMatrix
  # into a slice of a base matrix is a length error
  V <- lapply(seq_along(params), function(a) {
    M <- matrix(0, n, mf)
    if (npar[a] > 0L) {
      M[, offs[a] + seq_len(npar[a])] <- as_dense(design[[params[a]]]$X)
    }
    M
  })
  held <- statmod_structural_state(design)$held[[f$term]]
  keep <- if (length(held)) setdiff(seq_len(mf), nb + match(held, zn))
          else seq_len(mf)
  list(V = V, ap = ap, nb = nb, np = np, mf = mf, keep = keep, zn = zn,
       ev = ev, f = f, params = params, n = n)
}


#' The Exact Gradient Where a Penalty Covers a Filter's Own Parameters
#'
#' @description
#' Computes the exact gradient of the marginal criterion where a penalty covers
#' a structural term's own parameters, so that the determinant spans the joint
#' vector of the coefficients and those parameters.
#'
#' @details
#' The criterion is the same and so are its three pieces, except that
#' \eqn{K} carries the recursion's own second derivative,
#' \deqn{K = -\sum_t w_t\sum_{a,b}\ell_{ab,t}V_{a,t}^\top V_{b,t}
#'   - \sum_t w_t \ell_{p,t}E_t,}
#' so differentiating it along the direction \eqn{v} in which the mode moves
#' gives three contributions instead of one:
#'
#' 1. the family's third derivative against the per-observation diagonal
#'    of \eqn{M}, which is [u_vector()]'s formula with \eqn{V} in
#'    place of \eqn{X};
#' 2. the derivative of \eqn{V_p} itself, which is \eqn{E_t v}, one row
#'    per observation, entering through every \eqn{\ell_{pb}};
#' 3. the derivative of \eqn{\ell_p E_t}, whose first half is \eqn{E}
#'    re-weighted by \eqn{\sum_k \ell_{pk}(V_k\cdot v)} and whose second is
#'    \eqn{\partial^3 e_t/\partial u^3[v]}, the one genuinely new object.
#'
#' Piece 1 is direction-free and is computed once; pieces 2 and 3 are read
#' along each hyperparameter's own direction, which is the trade that keeps
#' the cost at \eqn{O(nm^2)} per hyperparameter instead of \eqn{O(nm^3)} once.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param idx The outer index.
#' @param basis The integrated subspace, or `NULL`.
#' @param free Whether to carry the result onto the free scale.
#'
#' @return A numeric vector, one entry per row of `idx`, or `NULL`
#'   where the determinant does not exist.
#'
#' @seealso [statmod_marginal_grad()],
#'   [modelterms7::term_third()]
#'
#' @keywords internal
statmod_structural_grad <- function(spec, design, coef, hyper, method, idx,
                                    basis = NULL, free = TRUE) {
  jd <- joint_design_rows(spec, design, coef)
  if (is.null(jd)) return(NULL)
  K <- tryCatch(statmod_marginal_full(spec, design, coef, hyper, NULL),
                error = function(e) NULL)
  if (is.null(K)) return(NULL)
  Kfac <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(Kfac)) return(NULL)
  Kinv <- chol2inv(Kfac)
  sst <- statmod_structural_state(design)
  keyt <- jd$f$term
  freep <- setdiff(jd$zn, sst$held[[keyt]])
  M <- if (is.null(basis)) Kinv else {
    A <- structural_joint_basis(spec, design, keyt, freep, jd$nb, basis)
    inner <- tryCatch(chol2inv(chol(crossprod(A, K %*% A))),
                      error = function(e) NULL)
    if (is.null(inner)) return(NULL)
    A %*% inner %*% t(A)
  }

  st <- structural_grad_parts(spec, design, coef, jd, M)
  params <- jd$params
  keep <- jd$keep
  flat <- unlist(coef[params], use.names = FALSE)
  # the joint twin of penalty_dbeta_trace(): every penalty whose Hessian moves
  # with what it covers moves the determinant as the mode moves, whether it
  # covers coefficients or a structural term's own parameters, each read where
  # it lives in the joint vector. A MIXED class is not read here: its prior is
  # multivariate, and outer_gradient_ok() refuses one that is not quadratic.
  jpd <- function(v) {
    s <- 0
    for (u in statmod_penalized(spec, design)) {
      if (isTRUE(u$mixed) || is.null(u$penalty)) next
      th <- as.list(hyper[[u$param]][[u$key]])
      pen <- u$penalty
      if (isTRUE(penalties7::beta_quadratic(pen, th))) next
      if (length(penalties7::penalty_kinks(pen, th))) next
      if (isTRUE(u$structural)) {
        bt <- as.numeric(sst$zeta[[u$term]][u$cols])
        pos <- match(jd$nb + u$cols, keep)
      } else {
        bt <- unit_beta(u, coef, params)
        pos <- match(u$index, keep)
      }
      if (!length(pos) || anyNA(pos)) next
      Tm <- as_dense(penalties7::penalty_dhessian_beta(pen, bt, th, v[pos]))
      s <- s + sum(as_dense(M[pos, pos, drop = FALSE]) * Tm)
    }
    s
  }
  out <- numeric(nrow(idx))
  links <- attr(idx, "links")
  # over the members, and accumulating: see statmod_marginal_grad()
  mem <- index_members(idx)
  for (i in seq_len(nrow(mem))) {
    r <- mem$row[i]
    p <- mem$parameter[i]
    nm <- mem$term[i]
    h <- mem$name[i]
    un <- statmod_unit(spec, design, p, nm)
    pen <- un$penalty
    if (isTRUE(un$mixed)) {
      # A CLASS SPLIT BETWEEN THE TWO HALVES is read where each of its
      # coordinates lives and put back in the class's own interleaved order,
      # which `joint` records: the same read joint_penalty_at() already does
      # for the value and the Hessian, so the ordering is composed in one
      # place. Those positions are ALREADY positions in the vector this
      # function assembles -- class_joint_pieces() numbers a structural
      # coordinate among the FREE parameters, which is what
      # statmod_marginal_full() spans -- so there is nothing to match them
      # into. Measured, matching them into `keep` a second time shifts the
      # structural half one place down and leaves the gradient 8.9e-02 out,
      # flat in the step size, which reads as a missing term and is not one.
      z <- sst$zeta[[un$class$sterm]]
      bt <- numeric(length(un$joint))
      for (pc in un$pieces) {
        bt[match(pc$joint, un$joint)] <- if (isTRUE(pc$structural)) {
          as.numeric(z[pc$zcols])
        } else {
          flat[pc$index]
        }
      }
      pos <- as.integer(un$joint)
    } else if (isTRUE(un$structural)) {
      z <- sst$zeta[[un$term]]
      bt <- as.numeric(z[un$cols])
      pos <- match(jd$nb + un$cols, keep)
    } else {
      bt <- coef[[p]][un$cols]
      pos <- match(un$index, keep)
    }
    # length zero as well as NA: a unit whose positions this vector does not
    # carry contributed an EMPTY scatter rather than being skipped, which is
    # how a mixed class read exactly zero in every coordinate
    if (!length(pos) || anyNA(pos)) next
    th <- as.list(hyper[[p]][[nm]])
    gt <- penalties7::penalty_grad_theta(pen, bt, th)
    cr <- penalties7::penalty_cross(pen, bt, th)
    dS <- penalties7::penalty_dhessian(pen, bt, th)
    c_m <- numeric(length(keep))
    c_m[pos] <- as.numeric(cr[[h]])
    v <- -as.numeric(Kinv %*% c_m)
    dS_m <- matrix(0, length(keep), length(keep))
    dS_m[pos, pos] <- as_dense(dS[[h]])
    chain <- sum(st$u * v) + structural_chain_extra(spec, design, jd, M, st, v)
    dtheta <- -as.numeric(gt[[h]]) - (sum(M * dS_m) + chain + jpd(v)) / 2
    out[r] <- out[r] + if (!free) dtheta else {
      eta <- linkfunctions7::linkfun(links[[r]], hyper[[p]][[nm]][[h]])
      dtheta * linkfunctions7::dlinkinv(links[[r]], eta)
    }
  }
  out
}


#' The Direction-Free Quantities of the Joint Chain Term
#'
#' @description
#' Computes the quantities of [statmod_structural_grad()] that do not depend
#' on which hyperparameter is being differentiated: the family's derivatives,
#' the filter's forward Jacobian, the per-observation diagonal of \eqn{M} and
#' the contraction \eqn{u}.
#'
#' @details
#' Computing them once keeps a model with several hyperparameters affordable,
#' because only the two pieces that read the direction are repeated, and each
#' of those is one pass of the recursion.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param jd The joint rows, from [joint_design_rows()].
#' @param M The matrix against which the trace is taken.
#'
#' @return A list with `u`, `V`, `Vk`, `VM`, `G`, `H`, `D3`, `D4`, `gl`,
#'   `s_at`, `sd_at`, `cd_at`, `seed`, `blocks`, `blocks_data` and `w`: the
#'   contraction \eqn{u}, the rows of the joint predictor derivative with the
#'   filter's Jacobian in place, their products with \eqn{M}, the
#'   per-observation diagonals of \eqn{M}, the family's derivatives, the score
#'   and curvature that drive the recursion, the seed of the recursion, the
#'   `blocks` callback with its data, and the observation weights.
#'
#' @keywords internal
structural_grad_parts <- function(spec, design, coef, jd, M) {
  # an outer search evaluates the criterion and its gradient at the same
  # mode more than once, and everything here depends only on the point and
  # on M; the memo returns the stored object itself, so a hit is
  # bit-identical to recomputing
  sst <- statmod_structural_state(design)
  structural_memo(design, "grad_parts",
                  list(coef = coef,
                       zeta = if (is.null(sst)) NULL else sst$zeta,
                       M = M), function() {
    structural_grad_parts_impl(spec, design, coef, jd, M)
  })
}

#' @rdname structural_grad_parts
#' @keywords internal
structural_grad_parts_impl <- function(spec, design, coef, jd, M) {
  if (isTRUE(jd$f$joint)) return(joint_grad_parts(spec, design, coef, jd, M))
  params <- jd$params
  n <- jd$n
  w <- spec@weights
  keep <- jd$keep
  ap <- jd$ap
  ev <- jd$ev
  f <- jd$f
  d <- spec@distrib
  y <- spec@response
  gl <- distributions7::distrib_gradient(d, y, ev$theta, scale = "link",
                                         threads = spec@threads)
  H <- distributions7::distrib_hessian(d, y, ev$theta, scale = "link",
                                       threads = spec@threads)
  D3 <- distributions7::distrib_deriv3(d, y, ev$theta, scale = "link",
                                       threads = spec@threads)
  D4 <- distributions7::distrib_deriv4(d, y, ev$theta, scale = "link",
                                       threads = spec@threads)
  # what drives the recursion and its derivatives, the log-likelihood's own
  # where gas(scaling = ) is zero; `s_at` stays the log-likelihood's score,
  # the weight every contraction below is taken against
  dr <- filter_driving(spec, ev$theta, ap, filter_scaling(f$tm), gl, H, D3, D4)
  s_at <- rep_len(gl[[f$param]], n)
  sd_at <- rep_len(dr$gl[[f$param]], n)
  cd_at <- rep_len(dr$H[[hess_key(params, ap, ap)]], n)

  # the filter's forward Jacobian, which is V for the equation it sits in.
  # The recursion is re-run at parameters it has already been run at, so the
  # callbacks LOOK UP the derivatives read at that predictor rather than
  # asking the family again -- the same bargain the information makes.
  seed <- jd$V[[ap]]
  mk_blocks <- .structural_blocks(params, ap, jd$V, dr$H, dr$D3, dr$D4, n)
  bd_data <- structural_blocks_data(params, ap, jd$V, dr$H, dr$D3, n,
                                    D4 = dr$D4)
  # the joint information reads the same recursion at the same point (its
  # blocks callback carries D4 as well, which the second order never
  # touches), so the memo's slot is shared with it, the seed in the key
  cv <- structural_memo(design, "curv",
                        list(zeta = f$psi, eta = f$eta_static,
                             g = w * s_at, seed = seed), function() {
    modelterms7::term_curvature(
      f$tm, f$eta_static, y, function(e, i) sd_at[i], function(e, i) cd_at[i],
      f$psi, w * s_at, seed, mk_blocks(NULL),
      score_values = sd_at, curvature_values = cd_at,
      blocks_data = bd_data, threads = spec@threads)
  })
  V <- jd$V
  V[[ap]] <- cv$jacobian
  Vk <- lapply(V, function(x) x[, keep, drop = FALSE])

  VM <- lapply(Vk, function(x) x %*% M)
  G <- vector("list", length(params))
  for (a in seq_along(params)) {
    G[[a]] <- lapply(seq_along(params), function(b) rowSums(VM[[a]] * Vk[[b]]))
  }
  keys3 <- names(D3)
  u <- numeric(length(keep))
  for (k in seq_along(params)) {
    s <- numeric(n)
    for (a in seq_along(params)) {
      for (b in seq_along(params)) {
        s <- s + rep_len(D3[[d3_key(params, a, b, k, keys3)]], n) * G[[a]][[b]]
      }
    }
    u <- u - as.numeric(crossprod(Vk[[k]], w * s))
  }
  # G is the per-observation diagonal of M read between two equations'
  # rows, which the second order needs as well; it is returned rather than
  # recomputed there, both readers being at the same point by construction
  list(u = u, V = V, Vk = Vk, VM = VM, G = G, H = H, D3 = D3, D4 = D4,
       gl = gl, s_at = s_at, sd_at = sd_at, cd_at = cd_at,
       seed = seed, blocks = mk_blocks,
       blocks_data = bd_data, w = w)
}


#' The Chain Term of the Joint Gradient Along a Direction
#'
#' @description
#' Computes the contributions to \eqn{\mathrm{tr}(M\,\partial K/\partial u[v])}
#' that come from the recursion and not from the design: the derivative of the
#' term's curvature along \eqn{v}, the curvature re-weighted by the derivative
#' of the score along \eqn{v}, and the movement of the filter's Jacobian, which
#' enters through every cross derivative of the filter's equation.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param jd The joint rows.
#' @param M The matrix against which the trace is taken.
#' @param st The shared quantities, from [structural_grad_parts()].
#' @param v The direction, over the estimated coordinates.
#'
#' @return A single number.
#'
#' @keywords internal
structural_chain_extra <- function(spec, design, jd, M, st, v) {
  if (isTRUE(jd$f$joint)) {
    return(joint_chain_extra(spec, design, jd, M, st, v))
  }
  params <- jd$params
  n <- jd$n
  ap <- jd$ap
  f <- jd$f
  w <- st$w
  keep <- jd$keep
  # the direction over the FULL width the term indexes by, a held level
  # carrying no component of it
  vfull <- numeric(jd$mf)
  vfull[keep] <- v
  # the derivative of each equation's predictor along v
  dV <- lapply(st$Vk, function(x) as.numeric(x %*% v))

  cv3 <- modelterms7::term_third(
    f$tm, f$eta_static, spec@response,
    function(e, i) st$sd_at[i], function(e, i) st$cd_at[i], f$psi,
    w * st$s_at, st$seed, st$blocks(vfull), vfull,
    score_values = st$sd_at, curvature_values = st$cd_at,
    blocks_data = st$blocks_data, threads = spec@threads)
  W3 <- cv3$curvature[keep, keep, drop = FALSE]
  dphi <- cv3$dphi[, keep, drop = FALSE]

  # (iii) the level's own weight moves too: E re-weighted by the derivative
  # of l_p along v, which is one more pass of the SECOND-order recursion
  kappa <- numeric(n)
  for (k in seq_along(params)) {
    kappa <- kappa + rep_len(st$H[[hess_key(params, ap, k)]], n) * dV[[k]]
  }
  cvk <- modelterms7::term_curvature(
    f$tm, f$eta_static, spec@response,
    function(e, i) st$sd_at[i], function(e, i) st$cd_at[i], f$psi,
    w * kappa, st$seed, st$blocks(NULL),
    score_values = st$sd_at, curvature_values = st$cd_at,
    blocks_data = st$blocks_data, threads = spec@threads)
  Wk <- cvk$curvature[keep, keep, drop = FALSE]

  # (ii) V_p itself moves, by E_t v, and it enters through every l_pb. The
  # two orderings of the pair are transposes and contribute equally
  b2 <- numeric(n)
  for (b in seq_along(params)) {
    b2 <- b2 + rep_len(st$H[[hess_key(params, ap, b)]], n) *
      rowSums(st$VM[[b]] * dphi)
  }
  -2 * sum(w * b2) - sum(M * Wk) - sum(M * W3)
}


#' The Model's Derivative Pieces for a Filter's Recursion
#'
#' @description
#' Builds the `blocks` callback that [modelterms7::term_curvature()],
#' [modelterms7::term_third()] and [modelterms7::term_fourth()] take, at a
#' direction, at two directions or without one.
#'
#' @details
#' The pieces are built on the active set that the term passes to the callback,
#' so a panel's outer products have the same size whatever the number of
#' groups.
#' `dcurv` serves twice, as the derivative of the curvature along the
#' direction and as the factor multiplying the movement of \eqn{V_p}, and
#' `N` is where the family's fourth derivative enters: each order of
#' differentiating the predictor through the recursion pulls in one more order
#' of the family.
#'
#' With TWO directions the callback serves [modelterms7::term_fourth()]
#' instead, and carries three quantities more: `Q`, the fourth
#' derivative with two of its indices on the filter's own equation;
#' `P`, the FIFTH derivative contracted against both directions, which
#' is the only place that order enters; and `cppp`, the scalar by which the
#' level's own curvature moves, which a contraction cannot recover.
#' `N` is then a list of two, one per direction.
#'
#' @param params The distribution's parameter names.
#' @param ap Which of them carries the filter.
#' @param Vs The static rows.
#' @param H,D3,D4,D5 The family's derivatives at the fitted predictors.
#'   `D5` is needed only where two directions are given.
#' @param n The number of observations.
#'
#' @return A function of the direction (`NULL`, one vector, or a list
#'   of two) returning a `blocks` callback.
#'
#' @keywords internal
.structural_blocks <- function(params, ap, Vs, H, D3, D4, n, D5 = NULL) {
  # The components are recycled ONCE, out here, and their keys built once
  # with them. The first version read them through
  # at <- function(x, i) rep_len(x, n)[i] inside the closure below, which
  # allocates an n-long vector -- and rebuilds a sort + paste key -- per
  # observation per parameter pair: O(n^2) work against a total of O(n),
  # and measured at 17.5% of a panel fit at 60 groups, growing with the
  # groups. The same hoist structural_callbacks() received in 0.18.0.
  np <- length(params)
  Hr <- vector("list", np)
  for (q in seq_len(np)) {
    if (q != ap) Hr[[q]] <- rep_len(H[[hess_key(params, ap, q)]], n)
  }
  D3r <- vector("list", np)
  for (r in seq_len(np)) {
    D3r[[r]] <- vector("list", np)
    for (r2 in seq_len(np)) {
      D3r[[r]][[r2]] <- rep_len(D3[[deriv3_key(params, ap, r, r2)]], n)
    }
  }
  D4r <- NULL
  if (!is.null(D4)) {
    D4r <- vector("list", np)
    for (r in seq_len(np)) {
      D4r[[r]] <- vector("list", np)
      for (r2 in seq_len(np)) {
        D4r[[r]][[r2]] <- vector("list", np)
        for (r3 in seq_len(np)) {
          D4r[[r]][[r2]][[r3]] <-
            rep_len(D4[[deriv4_key(params, ap, r, r2, r3)]], n)
        }
      }
    }
  }
  D5r <- NULL
  if (!is.null(D5)) {
    D5r <- vector("list", np)
    for (r in seq_len(np)) {
      D5r[[r]] <- vector("list", np)
      for (r2 in seq_len(np)) {
        D5r[[r]][[r2]] <- vector("list", np)
        for (r3 in seq_len(np)) {
          D5r[[r]][[r2]][[r3]] <- vector("list", np)
          for (r4 in seq_len(np)) {
            D5r[[r]][[r2]][[r3]][[r4]] <-
              rep_len(D5[[deriv5_key(params, ap, r, r2, r3, r4)]], n)
          }
        }
      }
    }
  }
  function(vfull) {
    dirs <- if (is.null(vfull)) list() else
      if (is.list(vfull)) vfull else list(vfull)
    third <- length(dirs) >= 1L
    fourth <- length(dirs) >= 2L
    function(e, i, D, act = NULL) {
      if (is.null(act)) act <- seq_len(ncol(Vs[[1L]]))
      mk <- length(act)
      cross <- numeric(mk)
      for (q in seq_along(params)) {
        if (q == ap) next
        cross <- cross + Hr[[q]][i] * Vs[[q]][i, act]
      }
      vr <- lapply(seq_along(params), function(r)
        if (r == ap) D else Vs[[r]][i, act])
      M <- matrix(0, mk, mk)
      for (r in seq_along(params)) {
        for (r2 in seq_along(params)) {
          M <- M + D3r[[r]][[r2]][i] * outer(vr[[r]], vr[[r2]])
        }
      }
      if (!third) return(list(cross = cross, M = M))
      # the predictor of every equation differentiated along each direction;
      # for the filter's own it is the CURRENT jacobian row, which only the
      # recursion has
      dvs <- lapply(dirs, function(vf)
        vapply(seq_along(params), function(r) sum(vr[[r]] * vf[act]),
               numeric(1)))
      dcurv <- numeric(mk)
      for (r in seq_along(params)) {
        dcurv <- dcurv + D3r[[ap]][[r]][i] * vr[[r]]
      }
      Nof <- function(dv) {
        N <- matrix(0, mk, mk)
        for (r in seq_along(params)) {
          for (r2 in seq_along(params)) {
            co <- 0
            for (r3 in seq_along(params)) {
              co <- co + D4r[[r]][[r2]][[r3]][i] * dv[r3]
            }
            if (co != 0) N <- N + co * outer(vr[[r]], vr[[r2]])
          }
        }
        N
      }
      if (!fourth) {
        return(list(cross = cross, M = M, dcurv = dcurv,
                    N = Nof(dvs[[1L]])))
      }
      # the fourth order's three further pieces. `Q` carries two of the
      # fourth derivative's indices on the filter's own equation, `P` is the
      # fifth contracted against both directions, and `cppp` is the scalar
      # the level's own curvature moves by.
      Q <- matrix(0, mk, mk)
      P <- matrix(0, mk, mk)
      for (r in seq_along(params)) {
        for (r2 in seq_along(params)) {
          qc <- D4r[[ap]][[r]][[r2]][i]
          if (qc != 0) Q <- Q + qc * outer(vr[[r]], vr[[r2]])
          pc <- 0
          for (r3 in seq_along(params)) {
            for (r4 in seq_along(params)) {
              pc <- pc + D5r[[r]][[r2]][[r3]][[r4]][i] *
                dvs[[1L]][r3] * dvs[[2L]][r4]
            }
          }
          if (pc != 0) P <- P + pc * outer(vr[[r]], vr[[r2]])
        }
      }
      list(cross = cross, M = M, dcurv = dcurv,
           N = list(Nof(dvs[[1L]]), Nof(dvs[[2L]])), Q = Q, P = P,
           cppp = D3r[[ap]][[ap]][i])
    }
  }
}

#' The Same Pieces as Data for the Compiled Recursion
#'
#' @description
#' Lays out the quantities that the callback of `.structural_blocks()` reads
#' as matrices, so that the compiled routes of `modelterms7` can read
#' them without calling back into R. The mixed second derivatives have one
#' column per distribution parameter (the filter's own column is zero and the
#' loop skips it). The third derivatives have one column per parameter pair,
#' and the pair \eqn{(r, r_2)} is at column \eqn{(r-1)\,n_p + r_2}. The static
#' jacobian rows are densified, and the filter's parameter index is included.
#'
#' @param params The distribution's parameter names.
#' @param ap Which of them carries the filter.
#' @param Vs The static rows.
#' @param H,D3 The family's derivatives at the fitted predictors.
#' @param n The number of observations.
#' @param D4,D5 The fourth and fifth derivatives, which the compiled third
#'   and fourth orders read: one column per index triple \eqn{(r, r_2, r_3)}
#'   after the filter's own, at \eqn{((r-1) n_p + r_2 - 1) n_p + r_3}, and
#'   likewise for the quadruples of the fifth. `NULL` leaves them out.
#'
#' @return A list with `H`, `D3`, `Vs` and `ap`, and `D4` and `D5` where
#'   they were given.
#'
#' @keywords internal
structural_blocks_data <- function(params, ap, Vs, H, D3, n, D4 = NULL,
                                   D5 = NULL) {
  np <- length(params)
  Hc <- matrix(0, n, np)
  for (q in seq_len(np)) {
    if (q != ap) Hc[, q] <- rep_len(H[[hess_key(params, ap, q)]], n)
  }
  D3m <- matrix(0, n, np * np)
  for (r in seq_len(np)) {
    for (r2 in seq_len(np)) {
      D3m[, (r - 1L) * np + r2] <-
        rep_len(D3[[deriv3_key(params, ap, r, r2)]], n)
    }
  }
  out <- list(H = Hc, D3 = D3m,
              Vs = lapply(Vs, function(x) if (is.matrix(x)) x else as_dense(x)),
              ap = ap)
  if (!is.null(D4)) {
    D4m <- matrix(0, n, np^3)
    for (r in seq_len(np)) for (r2 in seq_len(np)) for (r3 in seq_len(np)) {
      D4m[, ((r - 1L) * np + r2 - 1L) * np + r3] <-
        rep_len(D4[[deriv4_key(params, ap, r, r2, r3)]], n)
    }
    out$D4 <- D4m
  }
  if (!is.null(D5)) {
    D5m <- matrix(0, n, np^4)
    for (r in seq_len(np)) for (r2 in seq_len(np)) for (r3 in seq_len(np)) {
      for (r4 in seq_len(np)) {
        D5m[, (((r - 1L) * np + r2 - 1L) * np + r3 - 1L) * np + r4] <-
          rep_len(D5[[deriv5_key(params, ap, r, r2, r3, r4)]], n)
      }
    }
    out$D5 <- D5m
  }
  out
}


#' The Joint Subspace over Which a Marginal Criterion Integrates
#'
#' @description
#' Returns the basis onto which [ml()] projects the joint curvature: the
#' coefficients' own range basis, to which one column per penalized parameter
#' of the structural term is appended.
#'
#' @details
#' It is written once because [statmod_marginal_full()] and
#' [statmod_structural_grad()] must project onto the same subspace,
#' and a separate composition in each of the two callers would agree only by
#' accident.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param key The structural term's name.
#' @param free Its free parameters, in order.
#' @param nb The number of coefficients.
#' @param basis The coefficients' basis.
#'
#' @return A matrix with `nb` plus `length(free)` rows.
#'
#' @keywords internal
structural_joint_basis <- function(spec, design, key, free, nb, basis) {
  cols <- structural_range_cols(spec, design, key, free)
  A <- matrix(0, nb + length(free), ncol(basis) + length(cols))
  A[seq_len(nrow(basis)), seq_len(ncol(basis))] <- basis
  if (length(cols)) {
    A[cbind(nb + cols, ncol(basis) + seq_along(cols))] <- 1
  }
  A
}


#' The Per-Observation Diagonal of Each Block of a Matrix
#'
#' @description
#' Computes \eqn{G_{ab,i} = x_{ia}'M_{[a][b]}x_{ib}}, the diagonal of
#' \eqn{X_a M_{ab} X_b'} taken one observation at a time.
#'
#' @details
#' This is the one quantity to which every trace against \eqn{M} reduces. For a
#' matrix of the design's own form,
#' \deqn{\mathrm{tr}(M\,X'WX) = \sum_i w_i\,(XMX')_{ii},}
#' so a contraction of a third or fourth derivative never has to be assembled
#' as a \eqn{p\times p} matrix in order to be traced against \eqn{M}.
#'
#' It is computed once per evaluation point and read by the gradient's
#' contraction and by every pair of the Hessian.
#'
#' Where both blocks are sparse and the pairs of nonzeros in a row number
#' fewer than a quarter of \eqn{p_a p_b}, the diagonal is computed by
#' [leverage_pairs()] from the nonzeros of each row. Otherwise it is read off
#' the dense product that has the fewer columns.
#'
#' @param design The design.
#' @param M The matrix against which the traces are taken.
#' @param params The distribution's parameter names.
#' @param npar,offs The block sizes and their offsets.
#' @param threads How many threads the sparse route's kernel may use.
#'
#' @return A list of lists, `G[[a]][[b]]` a vector as long as the sample, or
#'   `NULL` where equation `a` or equation `b` has no coefficients.
#'
#' @seealso [u_vector()], [trace_design_form()]
#'
#' @keywords internal
block_leverage <- function(design, M, params, npar, offs, threads = 1L) {
  np <- length(params)
  tri <- vector("list", np)
  for (a in seq_len(np)) {
    if (npar[a] == 0L) next
    X <- design[[params[a]]]$X
    if (methods::is(X, "sparseMatrix")) {
      tri[[a]] <- row_nonzeros(X)
    }
  }
  G <- vector("list", np)
  for (a in seq_len(np)) G[[a]] <- vector("list", np)
  for (a in seq_len(np)) {
    if (npar[a] == 0L) next
    Xa <- design[[params[a]]]$X
    # x_ai' M_ab x_bi is a SCALAR, so the (b, a) block is the same number: it
    # is mirrored rather than computed twice, which used to double the work on
    # every off-diagonal pair.
    for (b in a:np) {
      if (npar[b] == 0L) next
      Xb <- design[[params[b]]]$X
      Mab <- M[offs[a] + seq_len(npar[a]), offs[b] + seq_len(npar[b]),
               drop = FALSE]
      # The gate compares the work of the two routes: the pairs of nonzeros
      # the kernel visits against the n pa pb products the dense route forms.
      # It used to require a product of densities below 1e-3 and fewer than
      # 2e7 pairs, which kept a smooth beside 200 random-effect levels (11
      # nonzeros a row of 210, density 0.052) on the dense route: 36 ms
      # against 3 ms for the same diagonal to 2.5e-15, and 524 ms against 3
      # ms at 1000 levels (2026-10-26, note/prestazioni.md). The cap dated
      # from a route that materialized the pairs; the kernel reads the
      # nonzeros row by row and needs no memory beyond them.
      pairs <- if (is.null(tri[[a]]) || is.null(tri[[b]])) Inf else
        length(tri[[a]]$i) / nrow(Xa) * length(tri[[b]]$i)
      v <- if (!is.null(tri[[a]]) && !is.null(tri[[b]]) && is.matrix(Mab) &&
               4 * pairs < as.double(nrow(Xa)) * npar[a] * npar[b]) {
        leverage_pairs(tri[[a]], tri[[b]], Mab, nrow(Xa), threads)
      } else if (npar[a] <= npar[b]) {
        # the intermediate is n x npar[a] this way round and n x npar[b] the
        # other, and the two give the same diagonal. Read the cheaper one: a
        # one-column block beside a five-hundred-column one was materializing
        # a dense n x 500 product to keep n numbers out of it.
        rowSums(Xa * (Xb %*% t(Mab)))
      } else {
        rowSums((Xa %*% Mab) * Xb)
      }
      G[[a]][[b]] <- v
      if (a != b) G[[b]][[a]] <- v
    }
  }
  G
}


#' A Design's Nonzeros, Ordered by Row
#'
#' @param X A sparse design block.
#'
#' @return A list with `i`, `j` and `v`, sorted by row.
#'
#' @seealso [leverage_pairs()]
#'
#' @keywords internal
row_nonzeros <- function(X) {
  tt <- methods::as(X, "TsparseMatrix")
  o <- order(tt@i)
  list(i = tt@i[o] + 1L, j = tt@j[o] + 1L, v = tt@x[o])
}


#' The Leverage Diagonal Over the Nonzeros of Two Rows
#'
#' @description
#' Computes
#' \eqn{G_i = \sum_{j\in J_i}\sum_{k\in K_i} X_a[i,j]X_b[i,k]M_{ab}[j,k]}, the
#' same quantity that [block_leverage()] otherwise reads off a dense
#' \eqn{n\times p_b} product.
#'
#' @details
#' Where a design is built from grouping indicators a row has one nonzero per
#' block, so the quadratic form is over a handful of entries, while the dense
#' product computes \eqn{p_a p_b} of them per row to keep one. The compiled
#' kernel `leverage_pairs_cpp()` reads the nonzeros of the two designs row by
#' row and sums them without expanding the pairs, and each row is written by
#' one thread. [block_leverage()] takes this route where both blocks are
#' sparse and the pairs of nonzeros in a row number fewer than a quarter of
#' the products that the dense route forms.
#'
#' @param ta,tb The two designs' nonzeros, from [row_nonzeros()].
#' @param Mab The block of \eqn{M}.
#' @param n The number of observations.
#' @param threads How many threads the kernel may use.
#'
#' @return A numeric vector as long as the sample.
#'
#' @seealso [block_leverage()]
#'
#' @keywords internal
leverage_pairs <- function(ta, tb, Mab, n, threads = 1L) {
  na <- tabulate(ta$i, n)
  nb <- tabulate(tb$i, n)
  if (!length(ta$i) || !length(tb$i)) return(numeric(n))
  # The expansion this used to build -- seven vectors as long as the PAIR
  # count, about 190 MB at the 3.38 million pairs a fit over 500 random-effect
  # levels reaches -- is what the kernel replaces: both designs arrive ordered
  # by row, so a row's two ranges are contiguous and the sum is a short double
  # loop per observation. Row i is written in full by one thread.
  astart <- if (n > 1L) c(0L, cumsum(na)[-n]) else 0L
  bstart <- if (n > 1L) c(0L, cumsum(nb)[-n]) else 0L
  leverage_pairs_cpp(as.integer(astart), as.integer(na),
                     as.integer(ta$j), as.numeric(ta$v),
                     as.integer(bstart), as.integer(nb),
                     as.integer(tb$j), as.numeric(tb$v),
                     Mab, as.integer(n), as.integer(threads))
}


#' The Trace Against a Contraction, Without Forming It
#'
#' @description
#' Computes \eqn{\mathrm{tr}(M\,T)} where \eqn{T} is the third derivative of
#' the penalized objective contracted once, or the fourth contracted twice.
#'
#' @details
#' [contract3()] and [contract4()] assemble
#' \eqn{-X_a'\,\mathrm{diag}(\omega_i w_i)\,X_b} block by block, and where the
#' result is only ever traced against \eqn{M} the assembly is waste: the trace
#' of that block against \eqn{M_{ab}} is \eqn{-\sum_i \omega_i w_i G_{ab,i}}.
#' The block and its transpose both contribute, so an off-diagonal pair counts
#' twice.
#'
#' The contraction of the derivative with the direction is formed exactly as
#' those two functions form it, which keeps the two routes on the same
#' arithmetic.
#'
#' @param spec A [StatmodSpec()].
#' @param G The per-observation diagonals, from [block_leverage()].
#' @param deriv The third or fourth derivatives on the link scale.
#' @param params,npar The parameter names and block sizes.
#' @param tv The direction in which the derivative is contracted.
#' @param tu A second direction, for a fourth derivative; `NULL` for a
#'   third.
#' @param key A function of the three or four parameter positions returning
#'   the component of `deriv` to read, or `NULL` for the observed route's keys,
#'   which are symmetric in every position. On the expected route `deriv` is
#'   the derivative of the expected information, symmetric in its first two
#'   positions only, and carries its own builder.
#'
#' @return A single number.
#'
#' @seealso [contract3()], [contract4()]
#'
#' @keywords internal
trace_design_form <- function(spec, G, deriv, params, npar, tv, tu = NULL,
                              key = NULL) {
  n <- spec@n_obs
  keys <- names(deriv)
  if (is.null(key)) {
    key <- if (is.null(tu)) function(a, b, k) d3_key(params, a, b, k, keys)
           else function(a, b, k, q) d4_key(params, a, b, k, q, keys)
  }
  total <- 0
  for (a in seq_along(params)) {
    if (npar[a] == 0L) next
    for (b in a:length(params)) {
      if (npar[b] == 0L) next
      w <- numeric(n)
      for (k in seq_along(params)) {
        if (npar[k] == 0L) next
        if (is.null(tu)) {
          w <- w + rep_len(deriv[[key(a, b, k)]], n) * tv[[k]]
        } else {
          for (q in seq_along(params)) {
            if (npar[q] == 0L) next
            w <- w + rep_len(deriv[[key(a, b, k, q)]], n) *
              tv[[k]] * tu[[q]]
          }
        }
      }
      mult <- if (a == b) 1 else 2
      total <- total - mult * sum(spec@weights * w * G[[a]][[b]])
    }
  }
  total
}
