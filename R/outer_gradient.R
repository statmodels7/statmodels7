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

#' Can the Exact Gradient Be Computed Here?
#'
#' @description
#' `TRUE` when every hyperparameter under estimation belongs to a penalty
#' whose second derivative in the coefficients is linear in the
#' hyperparameters and free of the coefficients, and the criterion uses the
#' observed information.
#'
#' @details
#' **Why the observed information.** \eqn{K} enters the criterion through
#' its determinant, so the gradient needs \eqn{\partial K/\partial\beta}. With
#' the observed information that is the third derivative of the log-likelihood
#' in the link-scale predictors, which every family of \pkg{distributions7}
#' carries in closed form. With the expected information it would be the
#' derivative in \eqn{\beta} of \eqn{-E[\ell'']}, which is not
#' \eqn{-E[\ell''']} and is not one of that package's generics. So
#' `reml(hessian = "observed")` is what the exact route asks for, and
#' `"expected"` keeps the derivative-free search.
#'
#' **Why the penalty is asked.** \eqn{\partial S/\partial\theta} and its
#' second derivative are generics of \pkg{penalties7}
#' ([penalties7::penalty_dhessian()],
#' [penalties7::penalty_d2hessian()],
#' [penalties7::penalty_dcross()]). A penalty that answers them is
#' estimable by a marginal criterion whatever its shape: the quadratic, the
#' additive, the structured and the separable branches all do, so a ridge, a
#' random effect and a heavy-tailed prior are covered as well as a spline. One
#' that does not, a SCAD or an MCP or anything with a kink, rejects, and the
#' search stays derivative-free. Nothing here tests a penalty's behavior to
#' find out what it is; it is asked.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param idx The outer index, from [outer_hyper_index()].
#' @param method An [OuterMethod()].
#' @param order `1` for the gradient, `2` for the Hessian as well.
#'
#' @return `TRUE` or `FALSE`.
#'
#' @seealso [statmod_marginal_grad()]
#'
#' @keywords internal
outer_gradient_ok <- function(spec, design, idx, method, order = 1L) {
  if (!nrow(idx)) return(FALSE)
  if (!identical(method@hessian, "observed")) {
    # The expected route asks for the same object in the same place -- the
    # movement of K with the coefficients -- but that object is
    # dE[l'']/deta rather than l''', because differentiating an expectation
    # moves the measure as well as the integrand. distributions7 supplies it
    # wherever the family writes its expected information out and refuses
    # elsewhere, so the family is ASKED rather than tested.
    #
    # Order 2 is not extended: the criterion's own second derivative would
    # want the next order of the same object, and lbfgs on the exact gradient
    # already buys most of what newton would -- which is the same judgement
    # the structural branch below records.
    if (order >= 2L) return(FALSE)
    if (!expected_deriv_ok(spec@distrib)) return(FALSE)
    # and not where a structural term is penalized: there the criterion is
    # statmod_marginal_full(), which assembles the joint curvature from
    # term_curvature() -- the OBSERVED one -- so that branch has no expected
    # criterion for this to be the derivative of.
    if (structural_penalized(spec, design)) return(FALSE)
  }
  seen <- unique(paste(idx$parameter, idx$term, sep = "\r"))
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
  if (structural_penalized(spec, design)) {
    if (order >= 2L) return(FALSE)
    for (u in statmod_penalized(spec, design)) {
      if (!isTRUE(u$structural)) next
      if (!answers_term_third(spec@terms[[u$param]][[u$term]])) return(FALSE)
    }
  }
  TRUE
}


#' Does the Family Supply the Expected Information's Derivative?
#'
#' @description
#' Whether [distributions7::distrib_dexpected_hessian()] answers for
#' this family, asked at a probe and never inferred from its class.
#'
#' @details
#' The default route in \pkg{distributions7} is one central difference of the
#' family's own expected information, which is a single stencil on an analytic
#' quantity wherever that information is a written-out formula and refuses
#' where it is itself an integral. Six of forty univariate families refuse, and
#' the reason is cost, never accuracy: measured at 100 observations they
#' cost 1880 to 147300 ms against a median of 0.183 ms for the others, so 2p
#' of those calls per criterion evaluation is not a slower route but an
#' unusable one.
#'
#' @param distrib A \pkg{distributions7} distribution.
#'
#' @return A single logical.
#'
#' @seealso [outer_gradient_ok()]
#'
#' @keywords internal
expected_deriv_ok <- function(distrib) {
  th <- tryCatch(distributions7::generate_random_theta(distrib),
                 error = function(e) NULL)
  if (is.null(th)) return(FALSE)
  y <- tryCatch(distributions7::distrib_rng(distrib, 1L, th),
                error = function(e) NULL)
  if (is.null(y)) return(FALSE)
  ok <- tryCatch({
    distributions7::distrib_dexpected_hessian(distrib, y, th, scale = "link")
    TRUE
  }, error = function(e) FALSE)
  isTRUE(ok)
}


#' Does a Term Supply Its Third Derivative?
#'
#' @description
#' Whether the term implements [modelterms7::term_third()], read
#' from the class the method is registered on, never from a list of
#' class names, so a term written later is covered without an edit here.
#'
#' @details
#' A structural term that has not written it inherits a method that signals an
#' error, and an additive term inherits one that returns zero, which is the
#' right answer for a predictor that is a block of columns. The question is
#' therefore whether the owning class is the refusing base.
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


#' Does a Penalty Supply What a Marginal Criterion Needs?
#'
#' @description
#' Asks the penalty for the derivatives, at a probe value of its
#' hyperparameters, and reports whether it answered.
#'
#' @details
#' The order-2 route additionally requires the penalty to be quadratic in the
#' coefficients ([penalties7::beta_quadratic()]), since otherwise the
#' third and fourth derivatives of the penalty in \eqn{\beta} enter the
#' criterion's own second derivative and \pkg{penalties7} does not carry them.
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
    if (order >= 2L) {
      if (!isTRUE(penalties7::beta_quadratic(pen, th))) return(FALSE)
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
#' \eqn{\partial V/\partial\eta} at the penalized mode, over the free scale of
#' the hyperparameters under estimation.
#'
#' @details
#' The three pieces are the envelope term \eqn{-\partial\rho/\partial\theta},
#' the explicit derivative of the determinant \eqn{\mathrm{tr}(M\,\partial
#' S/\partial\theta)}, and the implicit one \eqn{u'v}, where
#' \eqn{v = -(H+S)^{-1}\partial^2\rho/\partial\beta\partial\theta} is how the
#' mode moves and \eqn{u_c = \mathrm{tr}(M\,\partial K/\partial\beta_c)} is how
#' the determinant reads that movement.
#'
#' \eqn{u} is assembled without forming any third-derivative array. Writing
#' \eqn{G_{ab,i} = x_{ia}'M_{[a][b]}x_{ib}} for the per-observation diagonal of
#' the block of \eqn{M},
#' \deqn{u_{k} = -X_k'\Big(w \sum_{a,b} \ell'''_{abk}\, G_{ab}\Big),}
#' one crossprod per distribution parameter. The component
#' \eqn{\ell'''_{abk}} is looked up by a name built from the parameter names in
#' the family's own order, never parsed out of one.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param idx The outer index.
#' @param basis The integrated subspace, or `NULL`.
#' @param free Whether to carry the result onto the free scale. The Hessian
#'   asks for the parameter scale, having its own second-order chain rule to
#'   apply.
#'
#' @return A numeric vector, one entry per row of `idx`, or `NULL`
#'   where the determinant does not exist.
#'
#' @seealso [statmod_marginal()], [reml()]
#'
#' @keywords internal
statmod_marginal_grad <- function(spec, design, coef, hyper, method, idx,
                                  basis = NULL, free = TRUE, ctx = NULL) {
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
  design <- statmod_design_at(spec, coef, design)
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)

  # the information, the penalty, K = H + S and its inverse come from the
  # context when there is one, so a criterion, a gradient and a Hessian read
  # at the same point assemble them once between them
  expected <- identical(method@hessian, "expected")
  pen <- ctx_penalized(ctx, spec, design, coef, hyper, expected)
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
  mode_pen <- if (!expected) pen else
    ctx_penalized(ctx, spec, design, coef, hyper, FALSE)
  if (is.null(mode_pen)) return(NULL)
  # how the mode moves: the penalized likelihood's own curvature, which for a
  # refreshable block is not the Gauss-Newton matrix the design gives
  Dm <- mode_curvature(spec, design, coef, params, npar, offs, total)
  msolve <- NULL
  if (any(Dm != 0)) {
    Km <- as_dense(mode_pen$K) + Dm
    fac <- tryCatch(chol(Km), error = function(e) NULL)
    # the TRUE Hessian may lose definiteness where Gauss-Newton cannot, and a
    # mode's movement is worth having approximately rather than not at all
    if (!is.null(fac)) {
      msolve <- function(z) backsolve(fac, forwardsolve(t(fac), z))
    }
  }
  if (is.null(msolve)) {
    Kinv <- mode_pen$inv
    msolve <- function(z) as.numeric(Kinv %*% z)
  }
  M <- ctx_trace_matrix(ctx, pen, basis, expected)
  if (is.null(M)) return(NULL)

  # dK/dbeta is the third derivative of the log-likelihood on the observed
  # route and the derivative of the expected information on the expected one.
  # The contraction is the same in both -- one crossprod per distribution
  # parameter against the same leverage diagonal -- so only the array and the
  # key it is read by change.
  km <- ctx_kmove(ctx, spec, design, coef, hyper, method)
  u <- u_vector(spec, design, coef, M, params, npar, offs, total,
                d3 = km$deriv, key = km$key,
                G = ctx_leverage(ctx, design, M, params, npar, offs,
                                 spec@threads)) +
    u_refresh(spec, design, coef, M, params, npar, offs, total, expected,
              ctx_approx(ctx))

  out <- numeric(nrow(idx))
  links <- attr(idx, "links")
  for (a in seq_along(params)) {
    p <- params[a]
    rows <- which(idx$parameter == p)
    if (!length(rows)) next
    for (nm in unique(idx$term[rows])) {
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
      for (r in rows[idx$term[rows] == nm]) {
        h <- idx$name[r]
        # the mode moves by -(H+S)^-1 d2rho/dbeta dtheta, with H the penalized
        # LIKELIHOOD's own curvature -- see mode_curvature()
        c_m <- numeric(total)
        c_m[pos] <- as.numeric(cr[[h]])
        v <- -as.numeric(msolve(c_m))
        dS_m <- matrix(0, total, total)
        dS_m[pos, pos] <- as_dense(dS[[h]])
        dtheta <- -as.numeric(gt[[h]]) -
          (sum(M * dS_m) + sum(u * v)) / 2
        # and onto the free scale the search runs on
        out[r] <- if (!free) dtheta else {
          eta <- linkfunctions7::linkfun(links[[r]], hyper[[p]][[nm]][[h]])
          dtheta * linkfunctions7::dlinkinv(links[[r]], eta)
        }
      }
    }
  }
  out
}


#' The Trace of the Determinant's Movement With the Coefficients
#'
#' @description
#' \eqn{u_c = \mathrm{tr}(M\,\partial K/\partial\beta_c)}, assembled one
#' crossprod per distribution parameter.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param M The matrix the trace is taken against.
#' @param params The distribution's parameter names.
#' @param npar,offs,total The block sizes, their offsets and the total.
#'
#' @return A numeric vector as long as the stacked coefficients.
#'
#' @keywords internal
u_vector <- function(spec, design, coef, M, params, npar, offs, total,
                     d3 = NULL, G = NULL, key = NULL) {
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

  out <- numeric(total)
  for (k in seq_along(params)) {
    if (npar[k] == 0L) next
    s <- numeric(n)
    for (a in seq_along(params)) {
      if (npar[a] == 0L) next
      for (b in seq_along(params)) {
        if (npar[b] == 0L) next
        s <- s + rep_len(d3[[key(a, b, k)]], n) * G[[a]][[b]]
      }
    }
    out[offs[k] + seq_len(npar[k])] <-
      -as.numeric(crossprod(design[[params[k]]]$X, spec@weights * s))
  }
  out
}


#' How the Determinant Reads a Block That Moves With the Coefficients
#'
#' @description
#' The part of \eqn{u_c = \mathrm{tr}(M\,\partial K/\partial\beta_c)} that
#' [u_vector()] does not compute: everything coming from
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
#' **The derivative is asked of the term**, through
#' [modelterms7::term_block_contract()], and never differenced here.
#' Two reasons, both measured. A term carries its own chain rule, the links
#' on its parameters and a subformula's design; and a break-point column is a
#' step function in its break-point, so a difference quotient of it diverges as
#' the step shrinks and never converges. A term that does not implement the
#' contraction inherits zeros, which is exactly right for a fixed design.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, already refreshed at `coef`.
#' @param coef The coefficients at the penalized mode.
#' @param M The matrix the trace is taken against.
#' @param params,npar,offs,total The block bookkeeping.
#' @param expected Whether the criterion carries the expected information.
#' @param approx The approximation for the expected information.
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
                      expected = FALSE, approx = "bartlett", units = NULL,
                      Hl = NULL) {
  out <- numeric(total)
  if (is.null(units)) {
    units <- refresh_units(spec, design, coef, params, npar, offs)
  }
  if (!length(units)) return(out)
  if (is.null(Hl)) Hl <- refresh_hessian(spec, design, coef, expected, approx)
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


#' The Curvature the Mode Actually Moves By
#'
#' @description
#' What separates the true Hessian of the penalized log-likelihood from the
#' Gauss-Newton matrix [statmod_information_at()] returns, where a
#' term's block moves with its coefficients.
#'
#' @details
#' \eqn{v = \partial\hat\beta/\partial\theta} solves
#' \eqn{(\partial^2\rho/\partial\beta^2 - \partial^2\ell/\partial\beta^2)v =
#' \partial^2\rho/\partial\beta\partial\theta}, and
#' \eqn{\partial^2\ell/\partial\beta^2} is the TRUE second derivative:
#' \deqn{\sum_i\sum_{a,b}\ell_{ab}\frac{\partial\eta_a}{\partial\beta}
#'   \frac{\partial\eta_b}{\partial\beta}
#'   + \sum_i\sum_a \ell_a\frac{\partial^2\eta_a}{\partial\beta^2}.}
#' The first sum is what the design gives; the second is zero for every fixed
#' block and is not for a refreshable one, where \eqn{X} is the Jacobian and so
#' \eqn{\partial^2\eta_i/\partial\beta_c\partial\beta_d =
#' \partial X_{id}/\partial\beta_c}. It is supported on the term's own block,
#' and that is what keeps this cheap: one call of
#' [modelterms7::term_block_contract()] per column, weighted by the
#' score where [u_refresh()] weights by \eqn{M}.
#'
#' Note that it is the **mode's** matrix and not the criterion's. The
#' determinant is of
#' whatever [statmod_information_at()] assembles, and its derivative
#' reads that one; how the mode moves is a fact about the penalized likelihood
#' and reads this one. Confusing the two is the defect this file already
#' records for the expected information, in a second place.
#'
#' Measured on `nl(a ~ 0 + ridge(~grp))` against a finite difference of
#' the criterion with the mode refitted: the gradient goes from 1.5e-04 to
#' 8.8e-09 at 320 observations and from 1.6e-05 to 1.6e-09 at 960, on the
#' observed route and the expected one alike.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, already refreshed at `coef`.
#' @param coef The coefficients at the penalized mode.
#' @param params,npar,offs,total The block bookkeeping.
#'
#' @return A square matrix, zero everywhere no refreshable term reaches.
#'
#' @seealso [u_refresh()],
#'   [modelterms7::term_block_contract()]
#'
#' @keywords internal
mode_curvature <- function(spec, design, coef, params, npar, offs, total) {
  D <- matrix(0, total, total)
  rf <- attr(design, "refresh")
  st <- attr(design, "state")
  if (is.null(rf) || !length(rf) || is.null(st)) return(D)
  n <- spec@n_obs
  w <- spec@weights
  th <- statmod_eta(spec, design, coef)$theta
  gl <- distributions7::distrib_gradient(spec@distrib, spec@response, th,
                                         scale = "link", threads = spec@threads)
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
#' The name is built by putting the three parameter names in the family's own
#' order and joining them, the direction \pkg{distributions7} sanctions, and
#' then checked against the enumeration, never trusted.
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
    stop(sprintf("No third-derivative component '%s'.", want), call. = FALSE)
  }
  want
}


#' The Rows of the Joint Predictor Derivative
#'
#' @description
#' \eqn{V_{a,t} = \partial\eta_{a,t}/\partial u} over the coefficients of
#' every equation followed by a structural term's own parameters: the
#' equation's design placed in its own columns, and for the equation carrying
#' a filter the forward Jacobian of the recursion, which is dense.
#'
#' @details
#' It exists so that the contraction [u_vector()] performs is
#' written once. The formula there does not change where a filter is present;
#' only its operand does, \eqn{X} becoming \eqn{[X_p \mid D]}. Everything
#' downstream, the third derivative against the diagonal of \eqn{M} and the
#' movement of the mode, then reads the joint vector with no special case.
#'
#' The rows are returned at the full width, a held level included, and the
#' caller drops it exactly as [statmod_full_information()] does.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#'
#' @return A list with `V` (one matrix per distribution parameter),
#'   `ap` (which equation carries the filter), `keep` (the columns
#'   that survive a held level), `ev`, `f` and the sizes.
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
#' [statmod_marginal_grad()] over the joint vector of coefficients
#' and a structural term's parameters, over which the determinant spans
#' there.
#'
#' @details
#' The criterion is the same and so are its three pieces; what changes is that
#' \eqn{K} carries the recursion's own second derivative,
#' \deqn{K = -\sum_t w_t\sum_{a,b}\ell_{ab,t}V_{a,t}^\top V_{b,t}
#'   - \sum_t w_t \ell_{p,t}E_t,}
#' so differentiating it along the direction \eqn{v} the mode moves in gives
#' three contributions instead of one:
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
  out <- numeric(nrow(idx))
  links <- attr(idx, "links")
  for (r in seq_len(nrow(idx))) {
    p <- idx$parameter[r]
    nm <- idx$term[r]
    h <- idx$name[r]
    un <- statmod_unit(spec, design, p, nm)
    pen <- un$penalty
    if (isTRUE(un$structural)) {
      z <- sst$zeta[[un$term]]
      bt <- as.numeric(z[un$cols])
      pos <- match(jd$nb + un$cols, keep)
    } else {
      bt <- coef[[p]][un$cols]
      pos <- match(un$index, keep)
    }
    if (anyNA(pos)) next
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
    dtheta <- -as.numeric(gt[[h]]) - (sum(M * dS_m) + chain) / 2
    out[r] <- if (!free) dtheta else {
      eta <- linkfunctions7::linkfun(links[[r]], hyper[[p]][[nm]][[h]])
      dtheta * linkfunctions7::dlinkinv(links[[r]], eta)
    }
  }
  out
}


#' What the Joint Chain Term Needs Before a Direction Is Known
#'
#' @description
#' The quantities of [statmod_structural_grad()] that do not depend
#' on which hyperparameter is being differentiated: the family's derivatives,
#' the filter's forward Jacobian, the per-observation diagonal of \eqn{M} and
#' the contraction \eqn{u}.
#'
#' @details
#' Computing them once is what keeps a model with several hyperparameters
#' affordable: only the two pieces that read the direction are repeated, and
#' each of those is one pass of the recursion.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param jd The joint rows, from [joint_design_rows()].
#' @param M The matrix the trace is taken against.
#'
#' @return A list of the shared quantities.
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
  s_at <- rep_len(gl[[f$param]], n)
  c_at <- rep_len(H[[hess_key(params, ap, ap)]], n)

  # the filter's forward Jacobian, which is V for the equation it sits in.
  # The recursion is re-run at parameters it has already been run at, so the
  # callbacks LOOK UP the derivatives read at that predictor rather than
  # asking the family again -- the same bargain the information makes.
  seed <- jd$V[[ap]]
  mk_blocks <- .structural_blocks(params, ap, jd$V, H, D3, D4, n)
  bd_data <- structural_blocks_data(params, ap, jd$V, H, D3, n)
  # the joint information reads the same recursion at the same point (its
  # blocks callback carries D4 as well, which the second order never
  # touches), so the memo's slot is shared with it, the seed in the key
  cv <- structural_memo(design, "curv",
                        list(zeta = f$psi, eta = f$eta_static,
                             g = w * s_at, seed = seed), function() {
    modelterms7::term_curvature(
      f$tm, f$eta_static, y, function(e, i) s_at[i], function(e, i) c_at[i],
      f$psi, w * s_at, seed, mk_blocks(NULL),
      score_values = s_at, curvature_values = c_at,
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
  list(u = u, V = V, Vk = Vk, VM = VM, H = H, D3 = D3, D4 = D4,
       s_at = s_at, c_at = c_at, seed = seed, blocks = mk_blocks,
       blocks_data = bd_data, w = w)
}


#' The Two Pieces of the Chain Term That Read the Direction
#'
#' @description
#' The contributions to \eqn{\mathrm{tr}(M\,\partial K/\partial u[v])} that
#' come from the recursion, never from the design: the derivative of the
#' filter's own Jacobian, and the derivative of the term the level contributes
#' to the information.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param jd The joint rows.
#' @param M The matrix the trace is taken against.
#' @param st The shared quantities, from [structural_grad_parts()].
#' @param v The direction, over the estimated coordinates.
#'
#' @return A single number.
#'
#' @keywords internal
structural_chain_extra <- function(spec, design, jd, M, st, v) {
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
    function(e, i) st$s_at[i], function(e, i) st$c_at[i], f$psi,
    w * st$s_at, st$seed, st$blocks(vfull), vfull)
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
    function(e, i) st$s_at[i], function(e, i) st$c_at[i], f$psi,
    w * kappa, st$seed, st$blocks(NULL),
    score_values = st$s_at, curvature_values = st$c_at,
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
#' Builds the `blocks` callback [modelterms7::term_curvature()]
#' and [modelterms7::term_third()] take, at a direction or without
#' one.
#'
#' @details
#' The pieces are built on the active set the term asks for, so a panel's
#' outer products are of the same size whatever the number of groups.
#' `dcurv` serves twice, as the derivative of the curvature along the
#' direction and as the factor multiplying the movement of \eqn{V_p}, and
#' `N` is where the family's fourth derivative enters: each order of
#' differentiating the predictor through the recursion pulls in one more order
#' of the family.
#'
#' @param params The distribution's parameter names.
#' @param ap Which of them carries the filter.
#' @param Vs The static rows.
#' @param H,D3,D4 The family's derivatives at the fitted predictors.
#' @param n The number of observations.
#'
#' @return A function of the direction returning a `blocks` callback.
#'
#' @keywords internal
.structural_blocks <- function(params, ap, Vs, H, D3, D4, n) {
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
  function(vfull) {
    third <- !is.null(vfull)
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
      # the predictor of every equation differentiated along the direction;
      # for the filter's own it is the CURRENT jacobian row, which only the
      # recursion has
      dv <- vapply(seq_along(params), function(r) sum(vr[[r]] * vfull[act]),
                   numeric(1))
      dcurv <- numeric(mk)
      for (r in seq_along(params)) {
        dcurv <- dcurv + D3r[[ap]][[r]][i] * vr[[r]]
      }
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
      list(cross = cross, M = M, dcurv = dcurv, N = N)
    }
  }
}

#' The Same Pieces as Data for the Compiled Recursion
#'
#' @description
#' The quantities the callback of `.structural_blocks()` reads, laid out
#' as matrices so `modelterms7`'s compiled second-order route can read
#' them without calling back into R: the mixed second derivatives one column
#' per distribution parameter (the filter's own column zero, the loop skips
#' it), the third derivatives one column per parameter pair with pair
#' \eqn{(r, r_2)} at column \eqn{(r-1)\,n_p + r_2}, the static jacobian rows
#' densified, and the filter's parameter index.
#'
#' @param params The distribution's parameter names.
#' @param ap Which of them carries the filter.
#' @param Vs The static rows.
#' @param H,D3 The family's derivatives at the fitted predictors.
#' @param n The number of observations.
#'
#' @return A list with `H`, `D3`, `Vs` and `ap`.
#'
#' @keywords internal
structural_blocks_data <- function(params, ap, Vs, H, D3, n) {
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
  list(H = Hc, D3 = D3m,
       Vs = lapply(Vs, function(x) if (is.matrix(x)) x else as_dense(x)),
       ap = ap)
}


#' The Subspace a Marginal Criterion Integrates Over, Jointly
#'
#' @description
#' The basis [ml()] projects the joint curvature onto: the
#' coefficients' own range basis, with one column per penalized parameter of
#' the structural term appended.
#'
#' @details
#' It is written once because [statmod_marginal_full()] and
#' [statmod_structural_grad()] must project onto the same subspace,
#' and two callers composing it separately would agree only by accident.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param key The structural term's name.
#' @param free Its free parameters, in order.
#' @param nb The number of coefficients.
#' @param basis The coefficients' basis.
#'
#' @return A matrix.
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
#' \eqn{G_{ab,i} = x_{ia}'M_{[a][b]}x_{ib}}, the diagonal of \eqn{X_a M_{ab}
#' X_b'} taken one observation at a time.
#'
#' @details
#' This is the one quantity every trace against \eqn{M} reduces to. For a
#' matrix of the design's own form,
#' \deqn{\mathrm{tr}(M\,X'WX) = \sum_i w_i\,(XMX')_{ii},}
#' so a contraction of a third or fourth derivative never has to be assembled
#' as a \eqn{p\times p} matrix in order to be traced against \eqn{M}: measured
#' at 8000 observations and 69 coefficients, forming it and taking the trace
#' costs 25.5 ms where the weighted sum costs 0.031 ms, and at 20000
#' observations and 503 coefficients 3480 ms against 0.066 ms.
#'
#' It is computed once per evaluation point and read by the gradient's
#' contraction and by every pair of the Hessian.
#'
#' @param design The design.
#' @param M The matrix the traces are taken against.
#' @param params The distribution's parameter names.
#' @param npar,offs The block sizes and their offsets.
#' @param threads How many threads the sparse route's kernel may use.
#'
#' @return A list of lists, `G[[a]][[b]]` a vector as long as the sample.
#'
#' @seealso [u_vector()], [trace_design_form()]
#'
#' @keywords internal
block_leverage <- function(design, M, params, npar, offs, threads = 1L) {
  np <- length(params)
  dens <- rep(1, np)
  tri <- vector("list", np)
  for (a in seq_len(np)) {
    if (npar[a] == 0L) next
    X <- design[[params[a]]]$X
    if (methods::is(X, "sparseMatrix")) {
      dens[a] <- Matrix::nnzero(X) / max(1, prod(dim(X)))
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
      # The density gate bounds the RATIO of work and says nothing about its
      # size, and the expansion is materialized: a design carrying many dense
      # columns beside its indicators can pass the ratio and still ask for
      # hundreds of megabytes of pairs. The cap is on the absolute count.
      pairs <- if (is.null(tri[[a]]) || is.null(tri[[b]])) Inf else
        length(tri[[a]]$i) / nrow(Xa) * length(tri[[b]]$i)
      v <- if (!is.null(tri[[a]]) && !is.null(tri[[b]]) && is.matrix(Mab) &&
               dens[a] * dens[b] < 1e-3 && pairs < 2e7) {
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
#' \eqn{G_i = \sum_{j\in J_i}\sum_{k\in K_i} X_a[i,j]X_b[i,k]M_{ab}[j,k]}, the
#' same quantity [block_leverage()] otherwise reads off a dense
#' \eqn{n\times p_b} product.
#'
#' @details
#' Where a design is built from grouping indicators a row has one nonzero per
#' block, so the quadratic form is over a handful of entries and the dense
#' product computes \eqn{p_a p_b} of them per row to keep one. The pairs are
#' expanded once and the whole sum is vectorized.
#'
#' **It is taken only where it wins**, and the threshold is measured
#' and nothing about it is assumed. At a combined density of 3.6e-05 (a
#' random intercept
#' over 500 levels) it is 14.2 times the dense route; at 0.18 (three smooths
#' and a random effect) it is 50 times slower, R's per-element indexing being
#' far dearer than a BLAS flop, and at a dense block 50 times slower again.
#' Interpolating the two measurements puts the crossover at a combined density
#' near 1.1e-03, which is the gate: below it the route is taken, at it the two
#' cost the same, and above it the dense product stands.
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
#' \eqn{\mathrm{tr}(M\,T)} where \eqn{T} is the third derivative of the
#' penalized objective contracted once, or the fourth contracted twice.
#'
#' @details
#' [contract3()] and [contract4()] assemble
#' \eqn{-X_a'\,\mathrm{diag}(\omega_i w_i)\,X_b} block by block, and where the
#' result is only ever traced against \eqn{M} the assembly is waste: the trace
#' of that block against \eqn{M_{ab}} is \eqn{-\sum_i \omega_i w_i G_{ab,i}}.
#' The block and its transpose both contribute, so an off-diagonal pair counts
#' twice.
#'
#' The weights \eqn{w_i} are built exactly as those two functions build them,
#' and that keeps the two routes the same arithmetic in the same order
#' where it matters.
#'
#' @param spec A [StatmodSpec()].
#' @param G The per-observation diagonals, from [block_leverage()].
#' @param deriv The third or fourth derivatives on the link scale.
#' @param params,npar The parameter names and block sizes.
#' @param tv The direction the derivative is contracted in.
#' @param tu A second direction, for a fourth derivative; `NULL` for a
#'   third.
#'
#' @return A single number.
#'
#' @seealso [contract3()], [contract4()]
#'
#' @keywords internal
trace_design_form <- function(spec, G, deriv, params, npar, tv, tu = NULL) {
  n <- spec@n_obs
  keys <- names(deriv)
  total <- 0
  for (a in seq_along(params)) {
    if (npar[a] == 0L) next
    for (b in a:length(params)) {
      if (npar[b] == 0L) next
      w <- numeric(n)
      for (k in seq_along(params)) {
        if (npar[k] == 0L) next
        if (is.null(tu)) {
          w <- w + rep_len(deriv[[d3_key(params, a, b, k, keys)]], n) * tv[[k]]
        } else {
          for (q in seq_along(params)) {
            if (npar[q] == 0L) next
            w <- w + rep_len(deriv[[d4_key(params, a, b, k, q, keys)]], n) *
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
