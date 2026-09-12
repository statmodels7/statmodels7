#' @include outer_refresh.R
NULL

# The exact Hessian of the marginal criterion.
#
# Differentiating dV/dt_m = -rho_m - (1/2) tr(M K_m) once more:
#
#   d2V/dt_m dt_l = -rho_ml + b_m' J b_l
#                   + (1/2) tr(M K_l M K_m) - (1/2) tr(M dK_m/dt_l)
#
# The first two terms come from the envelope form of the gradient: rho_m is
# read at the mode, so its total derivative picks up c_m' b_l, and c_m = -J b_m
# turns that into b_m' J b_l.
#
# The determinant contributes twice. M itself moves, dM/dt_l = -M K_l M, which
# is the same expression whether M is K^-1 or A(A'KA)^-1 A'; and K_m moves,
#
#   dK_m/dt_l = S_ml + U[b_l, b_m] + T[b_ml],
#
# with T[v] the third derivative of the penalized objective in beta contracted
# once, U[v, u] the fourth contracted twice, and
#
#   b_ml = -J^-1 ( (S_l + T[b_l]) b_m + S_m b_l + c_ml ).
#
# Every term of dK_m/dt_l that would carry a third or fourth derivative of the
# PENALTY in beta is absent because the route requires beta_quadratic(); what
# is left of T and U is the log-likelihood's own third and fourth derivatives
# in the link-scale predictors, which distributions7 carries in closed form for
# every family.

#' The Exact Hessian of the Marginal Criterion
#'
#' @description
#' \eqn{\partial^2 V/\partial\eta^2} at the penalized mode, over the free scale
#' of the hyperparameters under estimation.
#'
#' @details
#' The pieces are those of the gradient differentiated once more: the
#' hyperparameter Hessian of the penalty, the movement of the mode read
#' through the penalized curvature, and two contributions from the
#' determinant, one from \eqn{M} moving and one from \eqn{K_m} moving.
#'
#' **Everything the penalty contributes is asked of the penalty**
#' ([penalties7::penalty_hess_theta()],
#' [penalties7::penalty_dhessian()],
#' [penalties7::penalty_d2hessian()],
#' [penalties7::penalty_dcross()]), so a penalty that is not
#' quadratic in its hyperparameters is covered by the same assembly with no
#' branch here: a ridge, a random effect, a structured prior.
#'
#' **Onto the free scale** the chain rule is second order and diagonal,
#' each hyperparameter having its own link: with \eqn{\theta = h(\eta)},
#' \deqn{\partial^2 V/\partial\eta_m\partial\eta_l =
#'   h_m'h_l'\,\partial^2 V/\partial\theta_m\partial\theta_l
#'   + \delta_{ml}\,h_m''\,\partial V/\partial\theta_m.}
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param idx The outer index.
#' @param basis The integrated subspace, or `NULL`.
#'
#' @return A square matrix, one row per row of `idx`, or `NULL` where
#'   the determinant does not exist.
#'
#' @seealso [statmod_marginal_grad()], [reml()]
#'
#' @keywords internal
statmod_marginal_hess <- function(spec, design, coef, hyper, method, idx,
                                  basis = NULL, ctx = NULL,
                                  inner = NULL) {
  # ⚠️ A MODEL CARRYING A STRUCTURAL TERM IS NOT ANSWERED BY THE ASSEMBLY
  # BELOW, which is written over the stacked coefficients while a filter's
  # own parameters move with the hyperparameter beside them -- measured, 0.86
  # per cent out on an unpenalized filter beside a smooth and essentially
  # zero on a penalized one. What answers it is the same assembly written on
  # the JOINT vector, statmod_structural_hess(), which needs a fourth order
  # through the recursion and the family's fifth derivative with it. A term
  # that does not supply the fourth -- regime() supplies neither it nor the
  # third -- keeps statmod_hess_stencil(), one central difference of the
  # exact gradient.
  # ⚠️ AND ONLY WHERE A PENALTY COVERS THE TERM'S OWN PARAMETERS. The joint
  # assembly differentiates the criterion whose determinant spans them, which
  # is the criterion only there; with no such penalty the determinant is over
  # the coefficients alone, and the joint assembly returns the second
  # derivative of a different function. Measured on a smooth beside an
  # unpenalized gas(1, 1): -3.97422 where a second difference of the
  # criterion with the mode refitted reads -4.00295, -4.00272 and -4.00275 at
  # h of 3e-2, 1e-2 and 5e-3, 0.71 per cent out and flat; the stencil of the
  # exact gradient reads -4.00270.
  if (length(attr(design, "structural"))) {
    tm <- structural_term_of(spec, design)
    if (!is.null(tm) && answers_term_fourth(tm) &&
        structural_penalized(spec, design)) {
      h <- tryCatch(statmod_structural_hess(spec, design, coef, hyper, method,
                                            idx, basis),
                    error = function(e) NULL)
      if (!is.null(h) && all(is.finite(h))) return(h)
    }
    return(statmod_hess_stencil(spec, design, coef, hyper, method, idx,
                                basis, inner))
  }
  # the block AT THE MODE, for the reason statmod_marginal_grad() records: the
  # leverage diagonal and the contractions must read the same block K was
  # assembled on, and for a refreshable term the design as it arrives is the
  # one built at the starting coefficients
  design <- statmod_design_at(spec, coef, design)
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  nh <- nrow(idx)

  # one assembly and one factorization for the criterion, this Hessian and the
  # gradient it reads at the end, where there used to be three and two
  # order 2 is refused on the expected route by outer_gradient_ok(), so this
  # reads the observed matrix; it is passed rather than assumed so that a
  # caller reaching here directly gets the criterion's own.
  pen <- ctx_penalized(ctx, spec, design, coef, hyper,
                       identical(method@hessian, "expected"))
  if (is.null(pen)) return(NULL)
  K <- pen$K
  Kinv <- pen$inv
  M <- ctx_trace_matrix(ctx, pen, basis,
                        identical(method@hessian, "expected"))
  if (is.null(M)) return(NULL)

  d3 <- ctx_deriv(ctx, spec, design, coef, hyper, 3L)
  d4 <- ctx_deriv(ctx, spec, design, coef, hyper, 4L)
  G <- ctx_leverage(ctx, design, M, params, npar, offs, spec@threads)

  # ⚠️ dX/dbeta reaches this assembly in THREE places -- the matrix dK/dt_m
  # below, the trace of dK_m/dt_l against M, and the twice-contracted fourth
  # derivative -- and correcting one alone is measurably worse than correcting
  # none: a weakly identified nl went from 2.02e-01 to 6.35e+02 with a standard
  # error of NaN when only the first was written, because a corrected dK/dt was
  # then contracted against a direction the next comment shows is 66 per cent
  # wrong and mixed with two uncorrected traces. The three are resolved from ONE
  # list of units so they cannot disagree about which blocks move; an empty list
  # leaves every expression below untouched.
  units <- refresh_units(spec, design, coef, params, npar, offs)
  Hl <- if (length(units)) {
    refresh_hessian(spec, design, coef, identical(method@hessian, "expected"),
                    ctx_approx(ctx))
  } else NULL
  # the trace of the refresh part of T[v] against M is v'u, with u the same
  # contraction the gradient already forms: term_block_contract() is the
  # adjoint of term_block_deriv(), so no second route to it is needed
  uref <- if (length(units)) {
    u_refresh(spec, design, coef, M, params, npar, offs, total,
              identical(method@hessian, "expected"), ctx_approx(ctx),
              units = units, Hl = Hl)
  } else NULL
  # the weight the block's SECOND derivative is paired with depends on neither
  # direction, so it is built once here rather than once per pair below
  acurv <- if (length(units)) {
    refresh_curv_amat(spec, design, M, params, npar, offs, units, Hl)
  } else NULL
  # the score the mode's own third derivative is weighted by, read once
  glref <- if (length(units)) {
    distributions7::distrib_gradient(
      spec@distrib, spec@response,
      statmod_eta(spec, design, coef)$theta, scale = "link",
      threads = spec@threads)
  } else NULL

  # ⚠️ HOW THE MODE MOVES is governed by the penalized LIKELIHOOD's own
  # curvature, which for a block that moves with its coefficients is not the
  # matrix the determinant is of -- the two that statmod_marginal_grad() keeps
  # apart through mode_curvature(), and that this assembly did not. Measured on
  # a weakly identified nl against two refits of the mode, b_m read off K is
  # 6.6e-01 wrong and 5.2e-05 read off K + Dm; everything below inherits it,
  # which is why correcting dK/dbeta alone moved the Hessian by nothing.
  Dm <- if (length(units)) {
    mode_curvature(spec, design, coef, params, npar, offs, total)
  } else NULL
  Jmat <- K
  msolve <- function(z) as.numeric(Kinv %*% z)
  if (!is.null(Dm) && any(Dm != 0)) {
    Jt <- as_dense(K) + Dm
    fac <- tryCatch(chol(Jt), error = function(e) NULL)
    # the true Hessian may lose definiteness where Gauss-Newton cannot, and a
    # mode's movement is worth having approximately rather than not at all
    if (!is.null(fac)) {
      Jmat <- Jt
      # as.numeric because one caller passes the result of a `%*%`, which would
      # otherwise carry a one-column matrix through every consumer below
      msolve <- function(z) as.numeric(backsolve(fac, forwardsolve(t(fac), z)))
    }
  }

  # the pieces each hyperparameter owns, in the stacked coefficient space
  pieces <- outer_pieces(spec, design, coef, hyper, idx, offs, total)
  bhat <- lapply(seq_len(nh), function(m) -msolve(pieces$c[[m]]))
  tv <- lapply(bhat, function(v) block_predictors(design, params, npar, offs,
                                                 v))
  # what each direction costs a moving block, once per hyperparameter: both
  # quantities depend on ONE direction, so building them inside the pair loop
  # would repeat an O(np^2) product for every pair
  dref <- if (length(units)) {
    lapply(seq_len(nh), function(m)
      refresh_direction(spec, design, M, params, npar, offs, d3, tv[[m]],
                        bhat[[m]], units))
  } else NULL
  Tm <- lapply(seq_len(nh), function(m) {
    Tb <- contract3(spec, design, d3, params, npar, offs, total, tv[[m]])
    if (!length(units)) return(Tb)
    R <- contract3_refresh(spec, design, params, npar, offs, total, dref[[m]],
                           Hl, units)
    Tb + R + t(R)
  })
  Km <- lapply(seq_len(nh), function(m) pieces$S[[m]] + Tm[[m]])

  out <- matrix(0, nh, nh)
  for (m in seq_len(nh)) {
    for (l in m:nh) {
      key <- pieces$pair[m, l]
      Sml <- pieces$S2[[key]]
      cml <- pieces$c2[[key]]
      # the block's SECOND derivative in the pair's two directions, read once
      # and consumed twice: the mode's second movement below and the
      # twice-contracted fourth derivative further down
      f2 <- if (length(units)) {
        lapply(units, function(un)
          refresh_dblock2(un, bhat[[l]][un$ra], bhat[[m]][un$ra], spec@n_obs))
      } else NULL
      rhs <- (pieces$S[[l]] + Tm[[l]]) %*% bhat[[m]] +
        pieces$S[[m]] %*% bhat[[l]] + cml
      # ⚠️ T[b_l] is the third derivative of the PENALIZED OBJECTIVE, whose
      # second is K + D and not K, so its own third derivative is not
      # dK/dbeta alone. The missing piece is dD/dbeta, and leaving it out is
      # not a small error: measured against the mode refitted and differenced
      # twice, b_ml is wrong by 7.6 to 9.0 per cent without it.
      if (length(units)) {
        rhs <- rhs + refresh_mode_third(spec, params, npar, units, Hl, glref,
                                        tv[[l]], dref[[m]], f2, total)
      }
      bml <- -msolve(rhs)
      # dK_m/dt_l enters ONLY through its trace against M, so the two
      # contractions are never assembled: tr(M X'WX) is a weighted sum of the
      # per-observation diagonal G. Measured at 8000 observations and 69
      # coefficients, forming the matrix and tracing costs 25.5 ms where the
      # sum costs 0.031 ms, and the pair loop does it twice per pair.
      tr_dKm <- sum(M * Sml) +
        trace_design_form(spec, G, d4, params, npar, tv[[l]], tv[[m]]) +
        trace_design_form(spec, G, d3, params, npar,
                          block_predictors(design, params, npar, offs, bml))
      if (length(units)) {
        tr_dKm <- tr_dKm + sum(uref * bml) +
          trace_refresh4(spec, M, params, npar, Hl, dref[[l]], dref[[m]],
                         units, G, d3, bhat[[l]], f2, acurv)
      }
      v <- -pieces$rho2[m, l] +
        sum(bhat[[m]] * as.numeric(Jmat %*% bhat[[l]])) +
        sum((M %*% Km[[l]]) * t(M %*% Km[[m]])) / 2 -
        tr_dKm / 2
      out[m, l] <- v
      out[l, m] <- v
    }
  }

  # and onto the free scale the search runs on
  links <- attr(idx, "links")
  g <- statmod_marginal_grad(spec, design, coef, hyper, method, idx, basis,
                             free = FALSE, ctx = ctx)
  h1 <- numeric(nh)
  h2 <- numeric(nh)
  for (r in seq_len(nh)) {
    v <- hyper[[idx$parameter[r]]][[idx$term[r]]][[idx$name[r]]]
    e <- linkfunctions7::linkfun(links[[r]], v)
    h1[r] <- linkfunctions7::dlinkinv(links[[r]], e)
    h2[r] <- linkfunctions7::d2linkinv(links[[r]], e)
  }
  out <- out * outer(h1, h1)
  diag(out) <- diag(out) + h2 * g
  out
}


#' The Per-Hyperparameter Pieces of the Outer Derivatives
#'
#' @description
#' The penalty's contributions to the criterion's gradient and Hessian, placed
#' in the stacked coefficient space.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param idx The outer index.
#' @param offs,total The block offsets and the total width.
#' @param order `1` for the first-order pieces alone, `2` for the
#'   second-order ones as well.
#'
#' @details
#' At order 1 the second-order generics are not called at all. That is what
#' lets a penalty supplying only `penalty_dhessian()` give an exact
#' gradient: asking it for a derivative the gradient does not use would have
#' rejected it for a quantity nobody wanted.
#'
#' @return A list with `S` (one matrix per hyperparameter) and `c`
#'   (one vector), and at order 2 also `S2` and `c2` (one per pair,
#'   keyed), `rho2` (the hyperparameter Hessian) and `pair` (the key
#'   of each pair).
#'
#' @keywords internal
outer_pieces <- function(spec, design, coef, hyper, idx, offs, total,
                         order = 2L) {
  params <- spec@distrib@params
  nh <- nrow(idx)
  Sm <- vector("list", nh)
  cm <- vector("list", nh)
  S2 <- list()
  c2 <- list()
  rho2 <- matrix(0, nh, nh)
  pair <- matrix("", nh, nh)

  # over the MEMBERS: a shared row stands for several penalties, and the
  # quantity it needs is the sum of theirs, exactly as the gradient's is.
  # Where nothing is shared there is one member per row and this is the loop
  # over rows that was here before.
  mem <- index_members(idx)
  for (r in seq_len(nh)) {
    Sm[[r]] <- matrix(0, total, total)
    cm[[r]] <- numeric(total)
  }
  terms <- unique(paste(mem$parameter, mem$term, sep = "\r"))
  for (s in terms) {
    bits <- strsplit(s, "\r", fixed = TRUE)[[1L]]
    p <- bits[1L]
    nm <- bits[2L]
    a <- match(p, params)
    un <- statmod_unit(spec, design, p, nm)
    # A PENALTY OVER A STRUCTURAL TERM'S OWN PARAMETERS acts on none of these
    # coordinates. Everything here is placed in a matrix over the stacked
    # COEFFICIENTS, and such a penalty has no position in that vector: its
    # coordinates are the term's own parameters, which the marginal criterion
    # spans in a joint matrix this function does not build. It was already
    # contributing nothing -- the writes below land at NULL positions and are
    # no-ops -- but it reached the penalty first, at an empty coefficient
    # vector, and a multivariate prior asked for a derivative at no rows warns
    # where a univariate one returns empty in silence.
    if (isTRUE(un$structural) || isTRUE(un$mixed)) next
    pos <- un$index
    pen <- un$penalty
    bt <- unit_beta(un, coef, params)
    th <- as.list(hyper[[p]][[nm]])
    dS <- penalties7::penalty_dhessian(pen, bt, th)
    cr <- penalties7::penalty_cross(pen, bt, th)
    lines <- which(mem$parameter == p & mem$term == nm)
    for (i in lines) {
      r <- mem$row[i]
      Sm[[r]][pos, pos] <- Sm[[r]][pos, pos] + as_dense(dS[[mem$name[i]]])
      cm[[r]][pos] <- cm[[r]][pos] + as.numeric(cr[[mem$name[i]]])
    }
    if (order < 2L) next
    # THE SECOND ORDER IS KEYED BY THE INDEX ROW PAIR, not by the pair of
    # hyperparameter NAMES within one term, and every entry ACCUMULATES. A
    # shared row stands for the hyperparameters of several penalties held at
    # one value, so its second derivative is the sum of theirs exactly as its
    # first is; and two names of ONE unit may sit in two different rows, which
    # is why the loop runs over that unit's member LINES and reads the
    # penalty's own table by name. Nothing outside the unit contributes to the
    # pairs it writes: the penalty is a sum, so no second derivative mixes two
    # units. Where nothing is shared each unit's lines map onto rows one for
    # one and this is the loop that was here before.
    d2S <- penalties7::penalty_d2hessian(pen, bt, th)
    dcr <- penalties7::penalty_dcross(pen, bt, th)
    ht <- penalties7::penalty_hess_theta(pen, bt, th)
    for (i in lines) {
      for (j in lines) {
        r <- mem$row[i]
        q <- mem$row[j]
        nk <- pair_key(mem$name[i], mem$name[j], names(ht))
        # the ROW pair, so that two units of one group land on one entry. The
        # (r, q) and (q, r) keys are distinct and each is filled, where a key
        # naming the pair of names served both at once.
        key <- paste0(r, "_", q)
        pair[r, q] <- key
        rho2[r, q] <- rho2[r, q] + as.numeric(ht[[nk]])[1L]
        if (is.null(S2[[key]])) {
          S2[[key]] <- matrix(0, total, total)
          c2[[key]] <- numeric(total)
        }
        S2[[key]][pos, pos] <- S2[[key]][pos, pos] + as_dense(d2S[[nk]])
        c2[[key]][pos] <- c2[[key]][pos] + as.numeric(dcr[[nk]])
      }
    }
  }
  if (order < 2L) return(list(S = Sm, c = cm))
  # two hyperparameters of DIFFERENT terms are independent: the penalty is a
  # sum, so no second derivative mixes them
  zeroS <- matrix(0, total, total)
  zeroc <- numeric(total)
  for (m in seq_len(nh)) {
    for (l in seq_len(nh)) {
      if (nzchar(pair[m, l])) next
      key <- paste0("\r", m, "_", l)
      pair[m, l] <- key
      S2[[key]] <- zeroS
      c2[[key]] <- zeroc
    }
  }
  list(S = Sm, c = cm, S2 = S2, c2 = c2, rho2 = rho2, pair = pair)
}


#' The Outer Hessian by One Difference of the Exact Gradient
#'
#' @description
#' The criterion's second derivative in the hyperparameters where no analytic
#' route exists: one central difference of [statmod_marginal_grad()], with the
#' coefficients refitted at every probe.
#'
#' @details
#' This is what a model carrying a structural term gets, and the reason is
#' that the analytic assembly of [statmod_marginal_hess()] spans the stacked
#' coefficients while a filter's own parameters are estimated beside them and
#' move with the hyperparameter as well. Extending it is not a matter of
#' bookkeeping: each order of differentiation through the recursion pulls in
#' one more order of the response's family, so the first derivative reads the
#' family's fourth through [modelterms7::term_third()] and the second would
#' read a fifth, which does not exist.
#'
#' Differencing an ANALYTIC quantity once is the licence this toolkit already
#' grants itself for the Student t's degrees of freedom and for the marginal
#' break-point's prior rows. What it forbids is a difference of a difference,
#' and there is one layer here.
#'
#' # The two probes start from the same place
#'
#' Both refit from the coefficients given and from the structural state as it
#' stands, restored before each probe, so the two differ in the hyperparameter
#' alone. That is what makes the result stable: the mode's own location error
#' is nearly the same at \eqn{+h} and \eqn{-h} and cancels in the difference
#' rather than being amplified by \eqn{1/h}. Measured on a penalized filter
#' against a second difference of the criterion, the result is FLAT over four
#' decades of the step -- -1.971435 at every \eqn{h} from 1e-2 down to 3e-5,
#' against a criterion second difference of -1.971456 -- where the assembled
#' Hessian reads -1.1e-08.
#'
#' # Where it refuses
#'
#' A stencil inherits the reproducibility of the quantity it differences,
#' divided by the step, so where the mode is poorly located the answer is
#' noise rather than a curvature. It is therefore computed TWICE, at `h` and
#' at `3 * h`, and refused where the two disagree by more than
#' [hess_stencil_tol()]. The two regimes are six orders apart and nothing
#' sits between them: measured, 5.5e-08 on a penalized filter, 9.3e-07 on an
#' unpenalized one beside a smooth and 1.1e-06 on an ordinary smooth, against
#' 6.0e-01 and 5.7e-01 on two mixed covariance classes whose correlation the
#' search left at \eqn{|z| = 8.53}, where the chart's conditioning is 1e10 and
#' the exact gradient itself reads 1e-3. A refusal there is the answer, and
#' the consumers already handle `NULL`.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients the probes start from.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param idx The hyperparameter index, from [outer_hyper_index()].
#' @param basis The integrated basis, or `NULL`.
#' @param inner The inner optimizer the probes refit with; `iwls()` where
#'   none is given.
#' @param h The step, on the free scale the search runs on.
#'
#' @return The Hessian on the free scale, or `NULL` where a probe could not
#'   be evaluated or the two steps disagree.
#'
#' @seealso [statmod_marginal_hess()], [statmod_marginal_grad()],
#'   [hess_stencil_step()]
#'
#' @keywords internal
statmod_hess_stencil <- function(spec, design, coef, hyper, method, idx,
                                 basis = NULL, inner = NULL,
                                 h = hess_stencil_step()) {
  if (!nrow(idx)) return(NULL)
  if (is.null(inner)) inner <- iwls()
  if (S7::S7_inherits(inner, Iwls)) inner <- iwls_resolve(inner, spec@distrib)
  cfg <- inner_settings(inner)
  blocks <- tryCatch(statmod_blocks(spec, design), error = function(e) NULL)
  if (is.null(blocks)) return(NULL)
  eta0 <- hyper_to_eta(hyper, idx)
  nh <- length(eta0)
  beta0 <- unlist(coef[spec@distrib@params], use.names = FALSE)
  # a structural term's own parameters live in an ENVIRONMENT the inner fit
  # writes into as it goes, so a probe moves them and the next one would
  # start from wherever the last left them. That is the ratchet outer_fit()
  # already guards against on a step the search takes back, and it is the
  # same guard here: the state is restored before every probe and once more
  # on the way out, so the caller's design is left as it was found.
  sst <- statmod_structural_state(design)
  z0 <- if (is.null(sst)) NULL else sst$zeta
  restore <- function() if (!is.null(z0)) sst$zeta <- z0
  on.exit(restore(), add = TRUE)
  grad_at <- function(eta) {
    restore()
    hy <- eta_to_hyper(eta, idx, hyper)
    r <- tryCatch(statmod_alternate(spec, design, blocks, hy, inner, beta0,
                                    cfg$expected, cfg$approx, cfg$maxit,
                                    cfg$tol, verbosity(0),
                                    hold_refresh = TRUE),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    v <- tryCatch(statmod_marginal_grad(spec, design, r$obj$split(r$par), hy,
                                        method, idx, basis),
                  error = function(e) NULL)
    if (is.null(v) || !all(is.finite(v))) NULL else v
  }
  at_step <- function(step) {
    H <- matrix(0, nh, nh)
    for (m in seq_len(nh)) {
      ep <- eta0; ep[[m]] <- ep[[m]] + step
      em <- eta0; em[[m]] <- em[[m]] - step
      gp <- grad_at(ep)
      gm <- grad_at(em)
      if (is.null(gp) || is.null(gm)) return(NULL)
      H[, m] <- (gp - gm) / (2 * step)
    }
    # the criterion's Hessian is symmetric; the two columns of a pair are
    # computed from different probes and agree only up to what the probes
    # resolve, so the average is taken rather than one of the two kept
    (H + t(H)) / 2
  }
  A <- at_step(h)
  if (is.null(A)) return(NULL)
  B <- at_step(3 * h)
  if (is.null(B)) return(NULL)
  sc <- max(abs(A))
  if (!is.finite(sc) || sc <= 0) return(NULL)
  if (max(abs(A - B)) / sc > hess_stencil_tol()) return(NULL)
  A
}


#' The Step and the Tolerance of the Outer Hessian's Stencil
#'
#' @description
#' The step [statmod_hess_stencil()] differences at, and how far its two
#' readings may disagree before it refuses.
#'
#' @details
#' Both are measured rather than taken from a library rule.
#' [numericals7::fd_step()] would give \eqn{\epsilon^{1/3}}, about 6e-6, which
#' is right for differencing a function evaluated to machine precision and
#' wrong here: the gradient is computed by refitting a mode, so what bounds
#' the step below is that reproducibility and not the rounding of a double.
#'
#' Swept on a penalized filter against a second difference of the criterion,
#' the result is 8.0e-05 out at \eqn{h = 0.1}, where the truncation still
#' shows, and then FLAT at 1.06e-05 -- the reference's own error -- for every
#' \eqn{h} from 1e-2 down to 3e-5. `1e-3` is two decades below where
#' truncation matters and two above where anything else begins.
#'
#' The tolerance separates two regimes that are SIX ORDERS apart, with
#' nothing between them: 5.5e-08, 9.3e-07 and 1.1e-06 where the mode is well
#' located, against 6.0e-01 and 5.7e-01 where it is not. `1e-3` sits three
#' orders above the worst resolved reading and three below the best
#' unresolved one, and a relative error of that size in the curvature is
#' 5e-4 in a standard error, under the fourth significant figure a summary
#' prints.
#'
#' @return A single number.
#'
#' @seealso [statmod_hess_stencil()]
#'
#' @keywords internal
hess_stencil_step <- function() 1e-3

#' @rdname hess_stencil_step
#' @keywords internal
hess_stencil_tol <- function() 1e-3


#' The Key of a Hyperparameter Pair
#'
#' @description
#' Locates the entry of a penalty's hyperparameter Hessian, which is keyed by
#' name, never by position.
#'
#' @param a,b The two hyperparameter names.
#' @param keys The names the penalty actually returned.
#'
#' @return A single string.
#'
#' @keywords internal
pair_key <- function(a, b, keys) {
  k <- paste(a, b, sep = "_")
  if (k %in% keys) return(k)
  k <- paste(b, a, sep = "_")
  if (k %in% keys) return(k)
  stop(sprintf("No hyperparameter component for '%s' and '%s'.", a, b),
       call. = FALSE)
}


#' The Predictors a Coefficient Direction Induces
#'
#' @description
#' \eqn{(X_k v_k)_i} for each distribution parameter, which is how a movement
#' of the coefficients is felt by the log-density.
#'
#' @param design The design.
#' @param params The parameter names.
#' @param npar,offs The block sizes and offsets.
#' @param v A stacked coefficient vector.
#'
#' @return A list of numeric vectors, one per parameter.
#'
#' @keywords internal
block_predictors <- function(design, params, npar, offs, v) {
  lapply(seq_along(params), function(k) {
    if (npar[k] == 0L) return(NULL)
    as.numeric(design[[params[k]]]$X %*% v[offs[k] + seq_len(npar[k])])
  })
}


#' The Third Derivative of the Objective Contracted Once
#'
#' @description
#' \eqn{T[v] = (\partial K/\partial\beta)\cdot v}, a matrix over the stacked
#' coefficients.
#'
#' @details
#' Each block is a weighted crossproduct, the weight being
#' \eqn{w_i\sum_k \ell'''_{abk}(X_k v_k)_i}: the third derivative never
#' appears as an array.
#'
#' On the expected route `d3` is the derivative of the expected information in
#' the predictors, [distributions7::distrib_dexpected_hessian()], and the same
#' contraction gives \eqn{(\partial H_E/\partial\beta)\cdot v}; that array is
#' symmetric in its first two positions only, so it is read through its own
#' `key`.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param d3 The third derivatives in the link scale, or on the expected route
#'   the derivative of the expected information.
#' @param params,npar,offs,total The block bookkeeping.
#' @param tv The predictors of the direction.
#' @param key A function of three parameter positions returning the name of the
#'   component of `d3` to read, or `NULL` for the observed route's key, which
#'   is symmetric in all three positions.
#'
#' @return A square matrix.
#'
#' @keywords internal
contract3 <- function(spec, design, d3, params, npar, offs, total, tv,
                      key = NULL) {
  keys <- names(d3)
  if (is.null(key)) key <- function(a, b, k) d3_key(params, a, b, k, keys)
  n <- spec@n_obs
  out <- zero_information(design, total)
  for (a in seq_along(params)) {
    if (npar[a] == 0L) next
    for (b in a:length(params)) {
      if (npar[b] == 0L) next
      w <- numeric(n)
      for (k in seq_along(params)) {
        # skipped where the direction does not move that equation's predictor,
        # which for a direction over the coefficients alone is where the
        # equation has none; over the joint vector a filter moves its
        # equation's predictor whatever its design holds
        if (is.null(tv[[k]])) next
        w <- w + rep_len(d3[[key(a, b, k)]], n) * tv[[k]]
      }
      blk <- -wcrossprod(design[[params[a]]]$X, spec@weights * w,
                         design[[params[b]]]$X, spec@threads)
      ra <- offs[a] + seq_len(npar[a])
      rb <- offs[b] + seq_len(npar[b])
      out[ra, rb] <- blk
      if (a != b) out[rb, ra] <- t(blk)
    }
  }
  out
}


#' The Fourth Derivative of the Objective Contracted Twice
#'
#' @description
#' \eqn{U[v, u] = (\partial^2 K/\partial\beta^2)\cdot(v, u)}, a matrix over the
#' stacked coefficients.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param d4 The fourth derivatives in the link scale.
#' @param params,npar,offs,total The block bookkeeping.
#' @param tv,tu The predictors of the two directions.
#'
#' @return A square matrix.
#'
#' @keywords internal
contract4 <- function(spec, design, d4, params, npar, offs, total, tv, tu) {
  keys <- names(d4)
  n <- spec@n_obs
  out <- zero_information(design, total)
  for (a in seq_along(params)) {
    if (npar[a] == 0L) next
    for (b in a:length(params)) {
      if (npar[b] == 0L) next
      w <- numeric(n)
      for (k in seq_along(params)) {
        if (npar[k] == 0L) next
        for (q in seq_along(params)) {
          if (npar[q] == 0L) next
          w <- w + rep_len(d4[[d4_key(params, a, b, k, q, keys)]], n) *
            tv[[k]] * tu[[q]]
        }
      }
      blk <- -wcrossprod(design[[params[a]]]$X, spec@weights * w,
                         design[[params[b]]]$X, spec@threads)
      ra <- offs[a] + seq_len(npar[a])
      rb <- offs[b] + seq_len(npar[b])
      out[ra, rb] <- blk
      if (a != b) out[rb, ra] <- t(blk)
    }
  }
  out
}


#' The Name of a Fourth-Derivative Component
#'
#' @description
#' Locates the \eqn{(a, b, k, q)} entry, by a name built from the parameter
#' names in the family's own order.
#'
#' @param params The parameter names.
#' @param a,b,k,q Indices into `params`.
#' @param keys The names the derivative returned.
#'
#' @return A single string.
#'
#' @keywords internal
d4_key <- function(params, a, b, k, q, keys) {
  want <- paste(params[sort(c(a, b, k, q))], collapse = "_")
  if (!want %in% keys) {
    stop(sprintf("No fourth-derivative component '%s'.", want), call. = FALSE)
  }
  want
}


#' The Smoothing-Parameter Correction to the Effective Degrees of Freedom
#'
#' @description
#' The amount by which \eqn{\mathrm{tr}[(H+S)^{-1}H]} understates the
#' complexity of a fit whose hyperparameters were themselves estimated.
#'
#' @details
#' The ordinary effective degrees of freedom read the smoothing parameters
#' as though they were known, and they were not: they were chosen from the
#' same data. Propagating their uncertainty into the coefficients gives the
#' corrected Bayesian covariance
#'
#' \deqn{V' = V_\beta + J V_\theta J^\top, \qquad
#'   J = \partial\hat\beta/\partial\theta,}
#'
#' and the corrected count is \eqn{\mathrm{tr}(V' H)}. \eqn{J} comes from
#' the implicit function theorem at the penalized mode: differentiating
#' \eqn{\partial(-\ell + \rho)/\partial\beta = 0} in the hyperparameter
#' gives \eqn{(H+S)J_k = -\partial^2\rho/\partial\beta\,\partial\theta_k},
#' whose right-hand side is \pkg{penalties7}'s `penalty_cross()`.
#' \eqn{V_\theta} is the inverse of the outer criterion's own Hessian,
#' which [statmod_marginal_hess()] returns, negated because that
#' criterion is a maximand.
#'
#' Everything is on the hyperparameter's link scale, which is where the
#' outer criterion optimizes and therefore the only scale on which its
#' Hessian is a variance.
#'
#' **Against mgcv.** This is the first of the two terms mgcv sums into
#' `edf2`, and it agrees with mgcv's to about 3e-4 on a univariate
#' smooth once the difference between the two bases is allowed for. mgcv
#' adds a second term for the Gaussian scale, which it profiles out of the
#' fit and whose uncertainty it must therefore add back; here every
#' distribution parameter carries its own equation and its own coefficients,
#' so that uncertainty is already inside \eqn{H}. The residual difference is
#' measured at 0.16, 0.10 and 0.05 effective parameters at n = 200, 400 and
#' 2000, falling with the sample size.
#'
#' **Where it does not apply.** A kinked penalty has no hyperparameter
#' the outer criterion estimates, [outer_hyper_index()] skipping it, so
#' there is nothing to propagate and the correction is zero. That is
#' not an approximation: the map from the hyperparameter to the penalized
#' mode turns a corner whenever a coefficient joins or leaves the active
#' set, and a delta method needs a derivative that does not exist there.
#'
#' **A model carrying a filter.** There the mode moves in the JOINT vector,
#' the coefficients followed by the structural term's own parameters, which
#' is what the criterion's determinant spans. A penalty over those
#' parameters is a column of no design, so read on the coefficients alone
#' [hyper_mode_cross()] skips it and the mode moves by nothing in exactly
#' the coordinates that penalty shrinks. Measured on a converged filter
#' with a penalized loading over ten groups of forty, the correction on the
#' coefficients was EXACTLY 0 against 1.108889 on the joint vector, a
#' quarter of that model's whole effective count of 4.4976, moving cAIC by
#' 2.22 and cBIC by 6.64.
#'
#' The count this corrects was ALREADY read on that vector.
#' [statmod_edf()] takes a structural model's effective degrees of freedom
#' from `joint_smoother_diag()`, which reads the very same two matrices,
#' so the base count was joint while its correction was on the
#' coefficients -- two halves of one number read on two different vectors.
#'
#' The two matrices are the ones [statmod_marginal_full()] and
#' [statmod_full_information()] already build, so there is no second
#' assembly to disagree with the criterion's. That route reads the
#' OBSERVED information, a filter having no expected one to offer, so
#' `expected` and `approx` do not reach it. A model with no structural
#' term of the filter shape is untouched by construction rather than by
#' tolerance, [statmod_marginal_full()] returning `NULL` there.
#'
#' @param spec A [StatmodSpec()].
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param design The design.
#' @param method The outer method that estimated them, or `NULL`.
#' @param expected Whether the information is the expected one. Not read
#'   where the model carries a filter, whose information is observed.
#' @param approx The approximation for the expected information. Not read
#'   where the model carries a filter.
#'
#' @return A list with `total`, the scalar correction, `per`, one entry per
#'   penalty key, and `n_hyper`, how many hyperparameters were estimated.
#'   Zero throughout where none was; a zero `total` beside a positive
#'   `n_hyper` means the curvature could not be read, and a caller
#'   reporting to a reader has to tell the two apart.
#'
#' @references
#' Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and model
#' selection for general smooth models. *Journal of the American
#' Statistical Association*, 111(516), 1548--1563.
#'
#' @seealso [statmod_marginal_hess()],
#'   [penalties7::penalty_cross()]
#'
#' @keywords internal
statmod_edf_correction <- function(spec, coef, hyper, design, method,
                                   expected = TRUE, approx = "opg") {
  # `n_hyper` is what tells a zero correction from an unavailable one: with
  # no estimated hyperparameter there is nothing to propagate and zero is the
  # answer, while with one there is something and zero means the curvature
  # could not be read. A reader told the wrong reason is worse off than
  # one told nothing.
  zero <- list(total = 0, per = numeric(0), n_hyper = 0L)
  if (is.null(method) || !method@kind %in% c("ml", "reml")) return(zero)
  params <- spec@distrib@params

  blocks <- statmod_blocks(spec, design)
  idx <- outer_hyper_index(spec, blocks)
  if (!nrow(idx)) return(zero)
  zero$n_hyper <- nrow(idx)

  # THE VECTOR THE MODE MOVES IN is the one the criterion's determinant
  # spans, and for a model carrying a filter that is the JOINT vector: the
  # coefficients followed by the term's own parameters. A penalty over those
  # parameters is a column of no design, so in coefficient space
  # hyper_mode_cross() skips it and the mode moves by nothing in exactly the
  # coordinates such a penalty shrinks. Measured on a converged filter with a
  # penalized loading over ten groups of forty, the correction came back
  # EXACTLY 0 -- the cross matrix identically zero, one penalty skipped --
  # where the joint vector gives 1.108889 against a total edf of 4.4976, so
  # cAIC moved 2.22 and cBIC 6.64. It was not a lower bound slightly low; it
  # was nothing at all.
  #
  # NEITHER MATRIX IS ASSEMBLED HERE. statmod_marginal_full() is the one
  # place K + S is built on that vector and statmod_full_information() the
  # one place K is, so this reads the two the criterion and vcov() read --
  # measured, an assembly written out here is identical() to the first.
  # The joint route reads the OBSERVED information, there being no expected
  # one for a filter, so `expected` and `approx` do not reach it.
  #
  # A model with NO filter takes the branch that was here before, by
  # construction rather than by tolerance: statmod_marginal_full() returns
  # NULL where the design carries no structural term of the filter shape,
  # and a term of the LIKELIHOOD shape -- regime() -- carries no penalty
  # over its own parameters for this to have skipped.
  M <- tryCatch(statmod_marginal_full(spec, design, coef, hyper),
                error = function(e) NULL)
  H <- if (is.null(M)) NULL else
    tryCatch(statmod_full_information(spec, coef, design),
             error = function(e) NULL)
  joint <- !is.null(M) && !is.null(H)
  if (joint) {
    Vb <- tryCatch(solve(as.matrix(M)), error = function(e) NULL)
  } else {
    H <- statmod_information_at(spec, coef, design, expected, approx)
    S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
    S <- zap_nonfinite(S)
    Vb <- tryCatch(solve(H + S), error = function(e) NULL)
  }
  if (is.null(Vb)) return(zero)

  # J = -Vb %*% d2rho/dbeta dtheta, one column per estimated hyperparameter
  J <- -Vb %*% hyper_mode_cross(spec, design, coef, hyper, idx,
                                nrow(as.matrix(H)), joint = joint)$cross

  Ho <- tryCatch(statmod_marginal_hess(spec, design, coef, hyper, method,
                                       idx, NULL),
                 error = function(e) NULL)
  if (is.null(Ho)) return(zero)
  # the criterion is a maximand, so its Hessian is negative definite at the
  # optimum and the variance is the inverse of its negative
  Vth <- tryCatch(solve(-as.matrix(Ho)), error = function(e) NULL)
  if (is.null(Vth)) return(zero)

  # per hyperparameter, so a summary can say which penalty the extra
  # complexity belongs to, and in total
  contrib <- vapply(seq_len(nrow(idx)), function(k) {
    v <- J[, k, drop = FALSE] %*% Vth[k, k, drop = FALSE] %*% t(J[, k, drop = FALSE])
    sum(v * t(H))
  }, numeric(1))
  per <- tapply(contrib, paste(idx$parameter, idx$term, sep = "\r"), sum)
  total <- sum((J %*% Vth %*% t(J)) * t(H))
  list(total = total, per = per, n_hyper = nrow(idx))
}


#' The Mixed Derivative of the Penalty in the Coefficients and the
#' Hyperparameters
#'
#' @description
#' \eqn{\partial^2\rho / \partial\beta \partial\theta}, written into the
#' stacked coefficient vector with one column per estimated hyperparameter.
#'
#' @details
#' This is the one ingredient of the penalized mode's movement that nothing
#' else computes, and both consumers of that movement read it here rather than
#' assembling it each: [statmod_edf_correction()], which contracts it against
#' the information to price what estimating a hyperparameter cost, and
#' [hyper_correction()], which keeps the matrix and adds it to a variance.
#'
#' A shared hyperparameter is ONE column standing for several penalties, so
#' each member writes into the group's column and they accumulate. Where
#' nothing is shared each member is its own row and the lookup is what was
#' here before the groups existed.
#'
#' A penalty over a STRUCTURAL term's own parameters covers positions among
#' those parameters rather than columns of a design, so where the caller's
#' matrix spans the coefficients alone there is nowhere to write it and it is
#' skipped. It is counted rather than passed over in silence: a correction
#' assembled without it is incomplete, and a caller reporting to a reader has
#' to be able to say so.
#'
#' With `joint` the matrix spans the vector the mode really moves in for such
#' a model, the coefficients followed by the term's own free parameters, and
#' those penalties have rows after all. What they lacked was an address and
#' not a derivative: [unit_joint_positions()] says where each unit's
#' coordinates live in that vector and [unit_joint_beta()] reads their values,
#' so nothing here is derived. [hyper_correction()] asks for it because
#' [vcov.StatmodFit()] inverts the joint penalized information; the
#' coefficient-space consumer does not, and its answer is unchanged.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param idx The hyperparameter index, from [outer_hyper_index()].
#' @param n How many rows the matrix carries: the stacked coefficients, and
#'   with `joint` a structural term's free parameters after them.
#' @param joint Whether those rows include that tail. `FALSE`, the default,
#'   is the coefficient-only matrix, which is what
#'   [statmod_edf_correction()] contracts.
#'
#' @return A list with `cross`, an `n` by `nrow(idx)` matrix, and `skipped`,
#'   how many penalties contributed nothing to it.
#'
#' @seealso [statmod_edf_correction()], [hyper_correction()],
#'   [penalties7::penalty_cross()]
#'
#' @keywords internal
hyper_mode_cross <- function(spec, design, coef, hyper, idx, n,
                             joint = FALSE) {
  cross <- matrix(0, n, nrow(idx))
  mem <- index_members(idx)
  skipped <- 0L
  for (un in statmod_penalized(spec, design)) {
    # A structural term's penalty indexes the term's OWN parameters, which
    # are not columns of any design. Where this matrix spans the coefficients
    # alone, writing it would land on whichever coefficients happen to occupy
    # those positions; where it spans the joint vector those parameters are
    # its tail and the positions are theirs.
    own <- isTRUE(un$structural) || isTRUE(un$mixed)
    pos <- if (own && !joint) integer(0)
           else unit_joint_positions(un, spec, design)
    if (!length(pos) || max(pos) > n) {
      skipped <- skipped + 1L
      next
    }
    cr <- tryCatch(
      penalties7::penalty_cross(un$penalty,
                                unit_joint_beta(un, spec, design, coef),
                                as.list(hyper[[un$param]][[un$key]]),
                                scale = "link"),
      error = function(e) NULL)
    if (is.null(cr)) {
      skipped <- skipped + 1L
      next
    }
    for (h in names(cr)) {
      # through the MEMBER table: a shared hyperparameter is one column
      # standing for several penalties, so each member writes into the
      # group's column and they accumulate.
      k <- mem$row[mem$parameter == un$param & mem$term == un$key &
                     mem$name == h]
      if (!length(k)) next
      cross[pos, k] <- cross[pos, k] + as.numeric(cr[[h]])
    }
  }
  list(cross = cross, skipped = skipped)
}


#' What a Hyperparameter's Own Uncertainty Adds to a Variance
#'
#' @description
#' \eqn{J V_\theta J'}, the movement of the penalized mode under the estimated
#' hyperparameters, as a matrix to be added to \eqn{V_b}.
#'
#' @details
#' The two matrices [vcov.StatmodFit()] reports by default are conditional on
#' the hyperparameters: both are read at the value the outer search stopped
#' at, as though it had been known. It was estimated from the same data, and
#' what that costs is
#' \deqn{V' = V_b + J V_\theta J', \qquad
#'   J = -(H + S)^{-1} \frac{\partial^2 \rho}{\partial\beta \partial\theta},}
#' the delta method applied to the map from the hyperparameter to the mode
#' (Wood, Pya and Safken, 2016). It is the same quantity
#' [statmod_edf_correction()] contracts against the information to obtain a
#' count of parameters; here the matrix itself is what is wanted.
#'
#' \eqn{V_b} is PASSED IN rather than recomputed, and that is what keeps the
#' two halves of the sum describing one model. The caller has already settled
#' which information \eqn{H} is, which coordinates are held and which are
#' aliased; a correction built on a second inverse, regularized differently,
#' would not be the movement of the mode whose variance it is added to.
#'
#' # Where there is nothing to add, and where it cannot be read
#'
#' The correction is exactly zero where no hyperparameter was estimated by a
#' differentiable criterion. A kinked penalty's is the argument of a minimum
#' over a grid, which [outer_hyper_index()] skips, and the map from it to the
#' mode turns a corner whenever a coefficient joins or leaves the active set,
#' so there is no derivative to propagate. That is a property of the model and
#' not a failure, and `n_hyper` is zero.
#'
#' It is UNAVAILABLE, with `n_hyper` positive and `C` `NULL`, where the
#' criterion's own Hessian cannot be read, which is where the search left a
#' coordinate at the edge of its range.
#'
#' It is PARTIAL, with `complete` false, where some of it could be read and
#' some could not: a hyperparameter [hyper_variance()] held contributes
#' nothing. The matrix returned is then a lower bound on the correction rather
#' than the whole of it.
#'
#' A penalty over a STRUCTURAL term's own parameters no longer makes it
#' partial. The mode of such a model moves in the joint vector, coefficients
#' and the term's free parameters together, which is what `Vb` already spans
#' here, and [hyper_mode_cross()] is asked for the same vector rather than for
#' the coefficients with a tail of zeros after them.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param method The outer method that estimated them, or `NULL`.
#' @param Vb The bayesian variance the correction is to be added to, over the
#'   coordinates `keep` names.
#' @param keep A logical vector over the joint vector, the stacked
#'   coefficients followed by any structural tail, saying which coordinates
#'   `Vb` spans. It is the caller's whole vector and not its head: a tail
#'   coordinate dropped as a flat direction leaves `Vb` narrower than the
#'   count alone would say.
#' @param nz How many of those coordinates are the structural tail.
#'
#' @return A list with `C`, the correction or `NULL`, `n_hyper`, how many
#'   hyperparameters a differentiable criterion estimated, and `complete`,
#'   whether every one of them contributed.
#'
#' @references
#' Wood, S. N., Pya, N. and Safken, B. (2016). Smoothing parameter and model
#' selection for general smooth models. *Journal of the American
#' Statistical Association*, 111(516), 1548--1563.
#'
#' @seealso [vcov.StatmodFit()], [statmod_edf_correction()],
#'   [statmod_hyper_vcov()]
#'
#' @keywords internal
hyper_correction <- function(spec, design, coef, hyper, method, Vb, keep,
                             nz = 0L) {
  none <- list(C = NULL, n_hyper = 0L, complete = TRUE)
  if (is.null(method) || !method@kind %in% c("ml", "reml")) return(none)
  idx <- outer_hyper_index(spec, statmod_blocks(spec, design))
  if (!nrow(idx)) return(none)
  out <- list(C = NULL, n_hyper = nrow(idx), complete = FALSE)

  # THE MATRIX SPANS WHAT Vb SPANS. For a model carrying a structural term
  # that is the joint vector, and the tail used to be padded with zeros while
  # hyper_mode_cross() skipped the penalty belonging in it, so the mode moved
  # by nothing in exactly the coordinates such a penalty shrinks -- measured,
  # the correction came back identical to Vb with the warning raised.
  cr <- hyper_mode_cross(spec, design, coef, hyper, idx, length(keep),
                         joint = nz > 0L)
  X <- cr$cross[keep, , drop = FALSE]
  Ho <- tryCatch(statmod_marginal_hess(spec, design, coef, hyper, method,
                                       idx, NULL),
                 error = function(e) NULL)
  if (is.null(Ho)) return(out)
  # the criterion is a maximand, so its Hessian is negative definite at the
  # optimum and the variance is the inverse of its negative
  Vth <- hyper_variance(-as.matrix(Ho))
  if (is.null(Vth)) return(out)
  # a coordinate hyper_variance() held has no variance to propagate. Its
  # column contributes nothing, which makes the correction partial rather
  # than wrong, and NA would make the whole product missing.
  ok <- is.finite(diag(Vth))
  if (!all(ok)) {
    Vth[!ok, ] <- 0
    Vth[, !ok] <- 0
  }
  J <- -Vb %*% X
  out$C <- as.matrix(J %*% Vth %*% t(J))
  out$complete <- all(ok) && cr$skipped == 0L
  out
}


#' The Variance of the Estimated Hyperparameters
#'
#' @description
#' Inverts the negative of the outer criterion's Hessian, holding any
#' coordinate whose own curvature cannot produce a variance.
#'
#' @details
#' At a maximum the criterion's Hessian is negative definite and its negative
#' inverts to a variance. A hyperparameter driven to the edge of its range,
#' or one the search left before reaching a maximum, has a curvature there
#' that is zero or of the wrong sign, and no variance follows from it.
#'
#' Such a coordinate used to cost every other one its standard error, the
#' whole matrix being refused. It is held instead and the rest is inverted,
#' which is the variance CONDITIONAL on it — the same reading `vcov()` gives
#' a coefficient the information carries nothing about. That is the marginal
#' variance only where the coupling contributes nothing to the kept
#' curvature, so the Schur correction \eqn{A_{kb}A_{bb}^{-1}A_{bk}} is
#' computed and compared against the kept diagonal rather than assumed
#' negligible; above `schur` the whole matrix is refused as before. The
#' default is the size at which the correction cannot move the four
#' significant digits the summary prints.
#'
#' @param A The negative of the outer Hessian, with dimnames.
#' @param schur The largest relative Schur correction a held coordinate may
#'   contribute to a kept one's curvature.
#'
#' @return A matrix of the same shape as `A` with the variance in the kept
#'   rows and columns and `NA` elsewhere, or `NULL` when no coordinate is
#'   usable or the coupling is too large to ignore.
#'
#' @seealso [statmod_hyper_vcov()], its only caller.
#'
#' @keywords internal
hyper_variance <- function(A, schur = 1e-4) {
  usable <- function(M) !is.null(M) && all(is.finite(M)) && all(diag(M) > 0)
  V <- tryCatch(solve(A), error = function(e) NULL)
  if (usable(V)) return(V)
  p <- ncol(A)
  d <- diag(A)
  bad <- which(!is.finite(d) | d <= 0 |
                 apply(!is.finite(A), 1L, any))
  keep <- setdiff(seq_len(p), bad)
  if (!length(bad) || !length(keep)) return(NULL)
  W <- tryCatch(solve(A[keep, keep, drop = FALSE]), error = function(e) NULL)
  if (!usable(W)) return(NULL)
  Bi <- tryCatch(solve(A[bad, bad, drop = FALSE]), error = function(e) NULL)
  if (is.null(Bi) || any(!is.finite(Bi))) return(NULL)
  Akb <- A[keep, bad, drop = FALSE]
  corr <- diag(as.matrix(Akb %*% Bi %*% t(Akb)))
  if (any(abs(corr) > schur * d[keep])) return(NULL)
  out <- matrix(NA_real_, p, p, dimnames = dimnames(A))
  out[keep, keep] <- W
  out
}


#' The Variance of the Hyperparameters a Marginal Criterion Estimated
#'
#' @description
#' The asymptotic variance matrix of the estimated hyperparameters, on the free
#' scale their links carry them onto.
#'
#' @details
#' A hyperparameter estimated by [reml()] or [ml()] is
#' the maximizer of a criterion that is twice differentiable in it, so it has a
#' variance like any other maximum-likelihood estimate: the inverse of the
#' negative Hessian of that criterion at the point reached, which
#' [statmod_marginal_hess()] already computes exactly. It is read on
#' the free scale because that is where the criterion was maximized and where
#' the quadratic approximation behind it is reasonable; a variance for a
#' smoothing parameter on its own scale, where the estimate is often several
#' orders of magnitude from zero and the criterion far from symmetric, would
#' describe a shape the criterion does not have.
#'
#' **Where it does not apply.** A hyperparameter chosen along a path, which
#' is what [aic()], [bic()] and [cv()] do to a kinked penalty's, is the
#' argument of a minimum over a grid. It is no root of a derivative, so there
#' is no Hessian to invert and no standard error follows. Its uncertainty is
#' a resampling question, and `NULL` is returned in place of a number of
#' another kind.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method The outer method that estimated them, or `NULL`.
#' @param inner The inner optimizer the fit used, which the stencil route
#'   refits its probes with; `iwls()` where none is given.
#'
#' @return A square matrix, one row per estimated hyperparameter, whose
#'   dimnames join the distribution parameter, the term and the
#'   hyperparameter's own name with a carriage return, a character no name
#'   can contain, so the three stay recoverable from the key. The index
#'   rides on the attribute `"idx"`. `NULL` where there is nothing to
#'   report.
#'
#' @seealso [statmod_marginal_hess()], [summary.StatmodFit()]
#'
#' @keywords internal
statmod_hyper_vcov <- function(spec, design, coef, hyper, method,
                               inner = NULL) {
  if (is.null(method) || !method@kind %in% c("ml", "reml")) return(NULL)
  blocks <- statmod_blocks(spec, design)
  idx <- outer_hyper_index(spec, blocks)
  if (!nrow(idx)) return(NULL)
  # A MIXED CLASS and a penalty over a STRUCTURAL term's own parameters were
  # both refused here, and neither is refused now: such a model carries a
  # structural term, so statmod_marginal_hess() answers it with
  # statmod_hess_stencil(), whose own two-step check refuses where the
  # curvature is not resolved. What that replaces is a wrong number rather
  # than a missing one -- measured on a panel of twenty groups whose level is
  # developed over a random effect, the assembled curvature read 1.09e-06
  # where the criterion's is 28.096, so the standard error came out 669.2 on
  # the free scale against 0.18866 and the interval covered the whole
  # positive line.
  basis <- integrated_basis(spec, design, method@kind)
  Ho <- tryCatch(statmod_marginal_hess(spec, design, coef, hyper, method, idx,
                                       basis, inner = inner),
                 error = function(e) NULL)
  if (is.null(Ho)) return(NULL)
  # the criterion is a maximand, so its Hessian is negative definite at the
  # optimum and the variance is the inverse of its negative. A hyperparameter
  # driven to the edge of its range leaves a curvature that is zero or of the
  # wrong sign there, and no interval follows from it.
  A <- -as.matrix(Ho)
  k <- paste(idx$parameter, idx$term, idx$name, sep = "\r")
  dimnames(A) <- list(k, k)
  V <- hyper_variance(A)
  if (is.null(V)) return(NULL)
  structure(V, idx = idx)
}


#' The Joint Derivative of the Penalized Information Along One Direction
#'
#' @description
#' \eqn{\partial K/\partial u\,[v]} assembled as a matrix over the
#' coefficients and the filter's own parameters, where
#' [structural_chain_extra()] returns only its trace against \eqn{M}.
#'
#' @details
#' The gradient needs the trace and nothing else, so it never forms this. The
#' HESSIAN needs the matrix twice over: in \eqn{\mathrm{tr}(M K_l M K_m)},
#' which no contraction reduces, and as the operator carrying the mode's
#' second movement. Writing \eqn{V_a} for each equation's rows,
#' \eqn{\mathrm{d}\varphi = E v} for the filter's own moving, and \eqn{E} for
#' the second derivative of the predictor the filter produces,
#'
#' \deqn{\frac{\partial K}{\partial u}[v] =
#'   -\sum_i w_i\Big[\sum_{a,b}\Big(\sum_k\ell_{abk}(V_k\cdot v)\Big)
#'     V_a^\top V_b
#'   + \sum_b \ell_{pb}\big(\mathrm{d}\varphi\otimes V_b
#'     + V_b\otimes\mathrm{d}\varphi\big)\Big]
#'   - W(\kappa_v) - W_3[v],}
#'
#' with \eqn{\kappa_v = \sum_k \ell_{pk}(V_k\cdot v)} the weight the level's
#' own term is re-read at and \eqn{W_3[v]} what
#' [modelterms7::term_third()] returns. Traced against \eqn{M} it
#' reproduces [structural_chain_extra()] exactly, which is what a test
#' asserts of it.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param jd The joint rows, from [joint_design_rows()].
#' @param st The shared quantities, from [structural_grad_parts()].
#' @param v The direction, over the estimated coordinates.
#'
#' @return A square symmetric matrix over the estimated coordinates.
#'
#' @seealso [structural_chain_extra()], [statmod_structural_hess()]
#'
#' @keywords internal
structural_dk_matrix <- function(spec, design, jd, st, v) {
  params <- jd$params
  n <- jd$n
  ap <- jd$ap
  f <- jd$f
  w <- st$w
  keep <- jd$keep
  nk <- length(keep)
  vfull <- numeric(jd$mf)
  vfull[keep] <- v
  dV <- lapply(st$Vk, function(x) as.numeric(x %*% v))

  cv3 <- modelterms7::term_third(
    f$tm, f$eta_static, spec@response,
    function(e, i) st$s_at[i], function(e, i) st$c_at[i], f$psi,
    w * st$s_at, st$seed, st$blocks(vfull), vfull)
  dphi <- cv3$dphi[, keep, drop = FALSE]
  out <- -cv3$curvature[keep, keep, drop = FALSE]

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
  out <- out - cvk$curvature[keep, keep, drop = FALSE]

  keys3 <- names(st$D3)
  for (a in seq_along(params)) {
    for (b in seq_along(params)) {
      cab <- numeric(n)
      for (k in seq_along(params)) {
        cab <- cab + rep_len(st$D3[[d3_key(params, a, b, k, keys3)]], n) *
          dV[[k]]
      }
      out <- out - crossprod(st$Vk[[a]] * (w * cab), st$Vk[[b]])
    }
  }
  for (b in seq_along(params)) {
    hb <- rep_len(st$H[[hess_key(params, ap, b)]], n)
    Z <- crossprod(dphi * (w * hb), st$Vk[[b]])
    out <- out - Z - t(Z)
  }
  (out + t(out)) / 2
}


#' The Joint Second Derivative of the Penalized Information, Traced
#'
#' @description
#' \eqn{\mathrm{tr}(M\,\partial^2 K/\partial u^2[v, w])}: the quantity
#' [structural_chain_extra()] gives at the first order, one order up
#' and contracted against a second direction.
#'
#' @details
#' Nine terms, and none of them assembles a matrix over the coefficients: a
#' term in \eqn{V_a^\top V_b} traces as a weighted sum of the
#' per-observation diagonal \eqn{G}, a term carrying
#' \eqn{\mathrm{d}\varphi} traces against the rows \eqn{M} has already been
#' applied to, and the three terms carrying the recursion's own derivatives
#' trace against what [modelterms7::term_curvature()],
#' [modelterms7::term_third()] and [modelterms7::term_fourth()]
#' return at the right weights.
#'
#' The last of those is the fourth derivative of the predictor through the
#' recursion, which is the object this whole route exists for, and the only
#' place the family's FIFTH derivative enters is `P` inside its
#' `blocks` callback.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param jd The joint rows.
#' @param M The matrix the trace is taken against.
#' @param st The shared quantities, from [structural_grad_parts()].
#' @param blk4 The `blocks` factory carrying the family's fifth
#'   derivative, as [.structural_blocks()] builds it.
#' @param v,w The two directions, over the estimated coordinates.
#'
#' @return A single number.
#'
#' @seealso [structural_chain_extra()], [statmod_structural_hess()]
#'
#' @keywords internal
structural_chain_extra2 <- function(spec, design, jd, M, st, blk4, v, w) {
  params <- jd$params
  np <- length(params)
  n <- jd$n
  ap <- jd$ap
  f <- jd$f
  wt <- st$w
  keep <- jd$keep
  vfull <- numeric(jd$mf)
  vfull[keep] <- v
  wfull <- numeric(jd$mf)
  wfull[keep] <- w
  dVv <- lapply(st$Vk, function(x) as.numeric(x %*% v))
  dVw <- lapply(st$Vk, function(x) as.numeric(x %*% w))
  H <- st$H
  D3 <- st$D3
  D4 <- st$D4
  keys3 <- names(D3)

  # the recursion's own quantities, along both directions and against both
  cv4 <- modelterms7::term_fourth(
    f$tm, f$eta_static, spec@response,
    function(e, i) st$s_at[i], function(e, i) st$c_at[i], f$psi,
    wt * st$s_at, st$seed, blk4(list(vfull, wfull)), list(vfull, wfull))
  dphi_v <- cv4$dphi[[1L]][, keep, drop = FALSE]
  dphi_w <- cv4$dphi[[2L]][, keep, drop = FALSE]
  dpsi <- cv4$dpsi[, keep, drop = FALSE]
  # v'E w, one per observation: the predictor's own second derivative read
  # along the two directions
  phi_vw <- as.numeric(dphi_w %*% v)

  # (ix) the level's term at the fourth order
  tot <- sum(M * cv4$curvature[keep, keep, drop = FALSE])

  # (vi) the level's term re-weighted by the second derivative of l_p along
  # the two directions
  gk <- numeric(n)
  for (k in seq_along(params)) {
    for (k2 in seq_along(params)) {
      gk <- gk + rep_len(D3[[d3_key(params, ap, k, k2, keys3)]], n) *
        dVv[[k]] * dVw[[k2]]
    }
  }
  gk <- gk + rep_len(H[[hess_key(params, ap, ap)]], n) * phi_vw
  cv0 <- modelterms7::term_curvature(
    f$tm, f$eta_static, spec@response,
    function(e, i) st$s_at[i], function(e, i) st$c_at[i], f$psi,
    wt * gk, st$seed, st$blocks(NULL),
    score_values = st$s_at, curvature_values = st$c_at,
    blocks_data = st$blocks_data, threads = spec@threads)
  tot <- tot + sum(M * cv0$curvature[keep, keep, drop = FALSE])

  # (vii) and (viii): the level's third derivative along one direction,
  # weighted by how l_p moves along the other
  kap <- function(dv) {
    out <- numeric(n)
    for (k in seq_along(params)) {
      out <- out + rep_len(H[[hess_key(params, ap, k)]], n) * dv[[k]]
    }
    out
  }
  kv <- kap(dVv)
  kw <- kap(dVw)
  for (pair in list(list(kv, wfull), list(kw, vfull))) {
    cv3 <- modelterms7::term_third(
      f$tm, f$eta_static, spec@response,
      function(e, i) st$s_at[i], function(e, i) st$c_at[i], f$psi,
      wt * pair[[1L]], st$seed, st$blocks(pair[[2L]]), pair[[2L]])
    tot <- tot + sum(M * cv3$curvature[keep, keep, drop = FALSE])
  }

  # (i) the family's fifth and fourth derivatives against the leverage
  # diagonal, which is the u_vector() identity one order up
  for (a in seq_along(params)) {
    for (b in seq_along(params)) {
      co <- numeric(n)
      for (k in seq_along(params)) {
        for (k2 in seq_along(params)) {
          co <- co + rep_len(D4[[deriv4_key(params, a, b, k, k2)]], n) *
            dVv[[k]] * dVw[[k2]]
        }
      }
      co <- co + rep_len(D3[[d3_key(params, a, b, ap, keys3)]], n) * phi_vw
      tot <- tot + sum(wt * co * st$G[[a]][[b]])
    }
  }

  # (ii) and (iii): V_p moves along one direction inside a third derivative
  # read along the other
  for (dd in list(list(dVv, dphi_w), list(dVw, dphi_v))) {
    for (b in seq_along(params)) {
      co <- numeric(n)
      for (k in seq_along(params)) {
        co <- co + rep_len(D3[[d3_key(params, ap, b, k, keys3)]], n) *
          dd[[1L]][[k]]
      }
      tot <- tot + 2 * sum(wt * co * rowSums(st$VM[[b]] * dd[[2L]]))
    }
  }

  # (iv) V_p's own SECOND movement, which is the third derivative of the
  # predictor contracted against both directions
  for (b in seq_along(params)) {
    hb <- rep_len(H[[hess_key(params, ap, b)]], n)
    tot <- tot + 2 * sum(wt * hb * rowSums(st$VM[[b]] * dpsi))
  }

  # (v) both movements of V_p at once
  tot <- tot + 2 * sum(wt * rep_len(H[[hess_key(params, ap, ap)]], n) *
                         rowSums((dphi_v %*% M) * dphi_w))
  -tot
}


#' The Penalty's Pieces on the Joint Vector
#'
#' @description
#' [outer_pieces()] over the coefficients AND a filter's own
#' parameters, which is the vector a marginal criterion's determinant spans
#' where a penalty covers those parameters.
#'
#' @details
#' The arithmetic is [outer_pieces()]'s and the difference is where each
#' unit's coordinates live: an ordinary unit is matched into the kept
#' coefficients, a penalty over a structural term's own parameters is placed
#' among the filter's free ones, and a MIXED covariance class is read where
#' each of its coordinates lives and put back in the class's own interleaved
#' order. Those three readings are the gradient's own, in
#' [statmod_structural_grad()], and are composed here so the two
#' cannot disagree about a position.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param idx The outer index.
#' @param jd The joint rows.
#'
#' @return A list with `S`, `c`, `S2`, `c2`, `rho2`
#'   and `pair`, in the shape [outer_pieces()] returns them.
#'
#' @seealso [statmod_structural_hess()], [outer_pieces()]
#'
#' @keywords internal
structural_outer_pieces <- function(spec, design, coef, hyper, idx, jd) {
  params <- jd$params
  keep <- jd$keep
  nk <- length(keep)
  nh <- nrow(idx)
  sst <- statmod_structural_state(design)
  flat <- unlist(coef[params], use.names = FALSE)
  Sm <- vector("list", nh)
  cm <- vector("list", nh)
  for (r in seq_len(nh)) {
    Sm[[r]] <- matrix(0, nk, nk)
    cm[[r]] <- numeric(nk)
  }
  S2 <- list()
  c2 <- list()
  rho2 <- matrix(0, nh, nh)
  pair <- matrix("", nh, nh)

  mem <- index_members(idx)
  terms <- unique(paste(mem$parameter, mem$term, sep = "\r"))
  for (s in terms) {
    bits <- strsplit(s, "\r", fixed = TRUE)[[1L]]
    p <- bits[1L]
    nm <- bits[2L]
    un <- statmod_unit(spec, design, p, nm)
    if (is.null(un)) next
    pen <- un$penalty
    bt <- unit_joint_beta(un, spec, design, coef)
    pos <- unit_joint_positions(un, spec, design)
    if (!length(pos) || anyNA(pos) || max(pos) > nk) next
    th <- as.list(hyper[[p]][[nm]])
    dS <- penalties7::penalty_dhessian(pen, bt, th)
    cr <- penalties7::penalty_cross(pen, bt, th)
    lines <- which(mem$parameter == p & mem$term == nm)
    for (i in lines) {
      r <- mem$row[i]
      Sm[[r]][pos, pos] <- Sm[[r]][pos, pos] + as_dense(dS[[mem$name[i]]])
      cm[[r]][pos] <- cm[[r]][pos] + as.numeric(cr[[mem$name[i]]])
    }
    d2S <- penalties7::penalty_d2hessian(pen, bt, th)
    dcr <- penalties7::penalty_dcross(pen, bt, th)
    ht <- penalties7::penalty_hess_theta(pen, bt, th)
    for (i in lines) {
      for (j in lines) {
        r <- mem$row[i]
        q <- mem$row[j]
        nk2 <- pair_key(mem$name[i], mem$name[j], names(ht))
        key <- paste0(r, "_", q)
        pair[r, q] <- key
        rho2[r, q] <- rho2[r, q] + as.numeric(ht[[nk2]])[1L]
        if (is.null(S2[[key]])) {
          S2[[key]] <- matrix(0, nk, nk)
          c2[[key]] <- numeric(nk)
        }
        S2[[key]][pos, pos] <- S2[[key]][pos, pos] + as_dense(d2S[[nk2]])
        c2[[key]][pos] <- c2[[key]][pos] + as.numeric(dcr[[nk2]])
      }
    }
  }
  zeroS <- matrix(0, nk, nk)
  zeroc <- numeric(nk)
  for (m in seq_len(nh)) {
    for (l in seq_len(nh)) {
      if (nzchar(pair[m, l])) next
      key <- paste0("\r", m, "_", l)
      pair[m, l] <- key
      S2[[key]] <- zeroS
      c2[[key]] <- zeroc
    }
  }
  list(S = Sm, c = cm, S2 = S2, c2 = c2, rho2 = rho2, pair = pair)
}


#' The Exact Outer Hessian of a Model Carrying a Structural Term
#'
#' @description
#' [statmod_marginal_hess()] over the joint vector of coefficients and
#' a filter's own parameters, which is what the determinant spans there.
#'
#' @details
#' # Why it needs a fourth order
#'
#' Each order of differentiating the predictor through the recursion pulls in
#' one more order of the response's family, the score the recursion is driven
#' by being read at the predictor it produces. The gradient reads
#' \eqn{\partial^3 e/\partial u^3} in one direction through
#' [modelterms7::term_third()]; the criterion's own second derivative
#' reads \eqn{\partial^4 e/\partial u^4} in two, through
#' [modelterms7::term_fourth()], and the family's FIFTH derivative with
#' it.
#'
#' # The shape is [statmod_marginal_hess()]'s
#'
#' With \eqn{u} the joint vector, \eqn{K} the penalized information over it
#' and \eqn{M} the matrix the trace is taken against,
#'
#' \deqn{\frac{\partial^2 V}{\partial t_m\partial t_l}
#'   = -\rho_{ml} + \hat b_m^\top K \hat b_l
#'   + \tfrac{1}{2}\mathrm{tr}(MK_lMK_m)
#'   - \tfrac{1}{2}\mathrm{tr}\Big(M\frac{\partial K_m}{\partial t_l}\Big),}
#'
#' with \eqn{\hat b_m = -K^{-1}c_m} the mode's movement, \eqn{K_m = S_m +
#' \partial K/\partial u[\hat b_m]} and the last trace carrying the penalty's
#' second derivative, the twice-contracted second derivative of \eqn{K} and
#' the once-contracted one at the mode's second movement. Every piece is the
#' one the coefficient-space assembly uses, read on the joint vector.
#'
#' # What it does not carry
#'
#' A block that MOVES with its coefficients -- `nl()`, `seg()` --
#' beside the filter contributes nothing here, exactly as it contributes
#' nothing to [statmod_structural_grad()]: that correction is written
#' in the coefficient-space assembly and has no joint twin. Such a model is
#' admitted at both orders and the approximation is the gradient's own.
#'
#' # Cost
#'
#' One [modelterms7::term_fourth()] and three lower recursions per PAIR
#' of hyperparameters, against the stencil's four refits per hyperparameter.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param idx The outer index.
#' @param basis The integrated subspace, or `NULL`.
#'
#' @return A square matrix on the free scale, one row per row of `idx`,
#'   or `NULL` where the joint matrix could not be formed.
#'
#' @seealso [statmod_structural_grad()], [statmod_marginal_hess()],
#'   [statmod_hess_stencil()] for the route it replaces.
#'
#' @keywords internal
statmod_structural_hess <- function(spec, design, coef, hyper, method, idx,
                                    basis = NULL) {
  if (!nrow(idx)) return(NULL)
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
  # the fifth derivative is read HERE and not in structural_grad_parts(),
  # whose result the gradient shares: the gradient never needs it, and it
  # costs 2p evaluations of the fourth
  # no `threads`: the fifth order is ONE central difference of the analytic
  # fourth, and the count reaches the family's own kernels through those two
  # calls rather than through this one
  D5 <- tryCatch(
    distributions7::distrib_deriv5(spec@distrib, spec@response,
                                   jd$ev$theta, scale = "link"),
    error = function(e) NULL)
  if (is.null(D5)) return(NULL)
  blk4 <- .structural_blocks(jd$params, jd$ap, jd$V, st$H, st$D3, st$D4,
                             jd$n, D5)

  nh <- nrow(idx)
  pieces <- structural_outer_pieces(spec, design, coef, hyper, idx, jd)
  msolve <- function(z) as.numeric(Kinv %*% z)
  bhat <- lapply(seq_len(nh), function(m) -msolve(pieces$c[[m]]))
  Tm <- lapply(bhat, function(v) structural_dk_matrix(spec, design, jd, st, v))
  Km <- lapply(seq_len(nh), function(m) pieces$S[[m]] + Tm[[m]])

  out <- matrix(0, nh, nh)
  for (m in seq_len(nh)) {
    for (l in m:nh) {
      key <- pieces$pair[m, l]
      Sml <- pieces$S2[[key]]
      cml <- pieces$c2[[key]]
      rhs <- as.numeric(Km[[l]] %*% bhat[[m]]) +
        as.numeric(pieces$S[[m]] %*% bhat[[l]]) + cml
      bml <- -msolve(rhs)
      tr_dKm <- sum(M * Sml) +
        structural_chain_extra2(spec, design, jd, M, st, blk4,
                                bhat[[m]], bhat[[l]]) +
        sum(st$u * bml) +
        structural_chain_extra(spec, design, jd, M, st, bml)
      v <- -pieces$rho2[m, l] +
        sum(bhat[[m]] * as.numeric(K %*% bhat[[l]])) +
        sum((M %*% Km[[l]]) * t(M %*% Km[[m]])) / 2 -
        tr_dKm / 2
      out[m, l] <- v
      out[l, m] <- v
    }
  }

  # and onto the free scale the search runs on
  links <- attr(idx, "links")
  g <- statmod_structural_grad(spec, design, coef, hyper, method, idx, basis,
                               free = FALSE)
  h1 <- numeric(nh)
  h2 <- numeric(nh)
  for (r in seq_len(nh)) {
    val <- hyper[[idx$parameter[r]]][[idx$term[r]]][[idx$name[r]]]
    e <- linkfunctions7::linkfun(links[[r]], val)
    h1[r] <- linkfunctions7::dlinkinv(links[[r]], e)
    h2[r] <- linkfunctions7::d2linkinv(links[[r]], e)
  }
  out <- out * outer(h1, h1)
  diag(out) <- diag(out) + h2 * g
  out
}
