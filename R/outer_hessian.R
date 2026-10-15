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
# T and U are the log-likelihood's own third and fourth derivatives in the
# link-scale predictors, which distributions7 carries in closed form for every
# family, plus -- for a penalty whose Hessian depends on the coefficients --
# the penalty's: dS/dbeta[v] inside T, and in dK_m/dt_l the second movement
# d2S/dbeta2[b_l, b_m] and dS_m/dbeta[b_l] + dS_l/dbeta[b_m], which
# statmod_penalty_second() assembles.

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
#' **Contributions of the penalty.** Every contribution of the penalty comes
#' from a generic of \pkg{penalties7}: [penalties7::penalty_hess_theta()],
#' [penalties7::penalty_cross()], [penalties7::penalty_dhessian()],
#' [penalties7::penalty_dcross()] and [penalties7::penalty_d2hessian()], with
#' [penalties7::penalty_dhessian_beta()] and its second-order companions for a
#' penalty whose Hessian depends on the coefficients. The assembly therefore
#' contains no branch for a particular penalty class (a ridge, a random
#' effect, a structured prior).
#'
#' **Onto the free scale** the chain rule is second order and diagonal,
#' each hyperparameter having its own link: with \eqn{\theta = h(\eta)},
#' \deqn{\partial^2 V/\partial\eta_m\partial\eta_l =
#'   h_m'h_l'\,\partial^2 V/\partial\theta_m\partial\theta_l
#'   + \delta_{ml}\,h_m''\,\partial V/\partial\theta_m.}
#'
#' **Routes.** For a model with a structural term the function returns
#' [statmod_structural_hess()] where the term implements
#' [modelterms7::term_fourth()] and a penalty covers the term's own parameters,
#' and [statmod_hess_stencil()] otherwise, including the case where that
#' assembly returns no finite matrix. On the expected information the
#' assembly runs only where [outer_gradient_ok()] admits order 2, and the
#' stencil is used elsewhere. Any other model is answered by the assembly
#' described above.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param idx The outer index.
#' @param basis The integrated subspace, or `NULL`.
#' @param ctx The evaluation context, or `NULL`.
#' @param inner The inner optimizer, read by the differenced route.
#' @param gam The coefficients the criterion estimates, from
#'   [marginal_coords()], or `NULL`. Their rows and columns follow those of
#'   the hyperparameters.
#'
#' @return A square matrix with one row for each row of `idx`, followed by one
#'   row for each coefficient in `gam`, or `NULL` where the penalized
#'   information or the matrix of the determinant cannot be formed or, on a
#'   differenced route, where the two steps disagree.
#'
#' @seealso [statmod_marginal_grad()], [reml()]
#'
#' @keywords internal
statmod_marginal_hess <- function(spec, design, coef, hyper, method, idx,
                                  basis = NULL, ctx = NULL,
                                  inner = NULL, gam = NULL) {
  ng <- if (is.null(gam)) 0L else length(gam$where)
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
                                basis, inner, gam = gam))
  }
  # ⚠️ ON THE EXPECTED ROUTE THE ASSEMBLY BELOW IS THE CRITERION'S SECOND
  # DERIVATIVE ONLY WHERE outer_gradient_ok() ADMITS ORDER 2, which asks the
  # family for d2E[l'']/deta^2. It used to run ungated for the three callers
  # that reach it without that gate -- statmod_hyper_vcov(),
  # hyper_correction() and statmod_edf_correction() -- and there it traced the
  # OBSERVED third and fourth derivatives against the EXPECTED penalized
  # matrix. Measured on one smooth at 4000 observations, 1e-4 to 3e-4 relative
  # against a second difference of the criterion (gamma1 -2.90522 against
  # -2.90568, beta1 -3.24139 against -3.24222), exact only on the Poisson,
  # whose log link is canonical. Where order 2 is not admitted the exact
  # gradient is differenced instead.
  expected <- identical(method@hessian, "expected")
  if (expected && !outer_gradient_ok(spec, design, idx, method, 2L,
                                     gamma = ng > 0L)) {
    return(statmod_hess_stencil(spec, design, coef, hyper, method, idx,
                                basis, inner, gam = gam))
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
  # ⚠️ TWO MATRICES ON THE EXPECTED ROUTE, and this assembly reads each in the
  # places it belongs. The determinant is of the criterion's K = H_exp + S, so
  # M, the traces tr(M K_l M K_m) and tr(M dK_m/dt_l), and the leverage
  # diagonal are read off it, with dK/dbeta the derivative of the expected
  # information (dE, d2E). The mode is where the penalized LIKELIHOOD's score
  # vanishes, so how it moves -- b_m, b_ml, the right-hand side of b_ml with
  # the observed third derivative, and b_m' J b_l -- is governed by
  # J = H_obs + S whatever the criterion's matrix is. On the observed route
  # the two coincide, which is why one T served both.
  pen <- ctx_penalized(ctx, spec, design, coef, hyper, expected,
                       laplace = TRUE)
  if (is.null(pen)) return(NULL)
  mode_pen <- ctx_penalized(ctx, spec, design, coef, hyper, FALSE)
  if (is.null(mode_pen)) return(NULL)
  K <- mode_pen$K
  Kinv <- mode_pen$inv
  M <- ctx_trace_matrix(ctx, pen, basis, expected)
  if (is.null(M)) return(NULL)

  d3 <- ctx_deriv(ctx, spec, design, coef, hyper, 3L)
  d4 <- if (expected) NULL else ctx_deriv(ctx, spec, design, coef, hyper, 4L)
  kE <- if (expected) ctx_kmove(ctx, spec, design, coef, hyper, method) else NULL
  kE2 <- if (expected) ctx_kmove2(ctx, spec, design, coef, hyper) else NULL
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
  # the coordinates pinned in K other than at a boundary (the coefficients
  # the criterion estimates and those the model does not identify) and those
  # a break-point's kink holds at the mode. The mode does not move in them,
  # so every solve below leaves them at zero; see mode_solver()
  ms <- mode_solver(spec, design, coef, K, Kinv, Dm, npar, offs)
  hold <- ms$hold
  Dm <- ms$Dm
  Jmat <- ms$J
  msolve <- ms$solve

  # the pieces each hyperparameter owns, in the stacked coefficient space
  pieces <- outer_pieces(spec, design, coef, hyper, idx, offs, total)
  bhat <- lapply(seq_len(nh), function(m) -msolve(pieces$c[[m]]))
  # THE COEFFICIENTS THE CRITERION ESTIMATES are outer coordinates with no
  # penalty of their own: S_j, c_j and every second-order penalty piece are
  # zero, and the direction the mode moves in carries the coordinate itself,
  # b_j = e_j - P K_m e_j, with K_m the mode's full curvature.
  Hm <- NULL
  if (ng) {
    Hm <- as_dense(ctx_information(ctx, spec, design, coef, hyper, FALSE,
                                   ctx_approx(ctx))) +
      as_dense(ctx_penalty(ctx, spec, design, coef, hyper))
    if (!is.null(Dm)) Hm <- Hm + Dm
    for (w in gam$where) {
      col <- Hm[, w]
      col[hold] <- 0
      b <- -msolve(col)
      b[hold] <- 0
      b[w] <- 1
      bhat[[length(bhat) + 1L]] <- b
    }
    pieces$S <- c(pieces$S, rep(list(matrix(0, total, total)), ng))
  }
  nv <- nh + ng
  tv <- lapply(bhat, function(v) block_predictors(design, params, npar, offs,
                                                 v))
  # what each direction costs a moving block, once per hyperparameter: both
  # quantities depend on ONE direction, so building them inside the pair loop
  # would repeat an O(np^2) product for every pair
  dref <- if (length(units)) {
    lapply(seq_len(nv), function(m)
      refresh_direction(spec, design, M, params, npar, offs, d3, tv[[m]],
                        bhat[[m]], units))
  } else NULL
  Tm <- lapply(seq_len(nv), function(m) {
    Tb <- contract3(spec, design, d3, params, npar, offs, total, tv[[m]])
    if (!length(units)) return(Tb)
    R <- contract3_refresh(spec, design, params, npar, offs, total, dref[[m]],
                           Hl, units)
    Tb + R + t(R)
  })
  # a penalty whose Hessian moves with the coefficients moves K along b_m too;
  # NULL, and nothing added, where every penalty is quadratic in them
  Pm <- lapply(seq_len(nv), function(m)
    statmod_penalty_dbeta(spec, design, coef, hyper, bhat[[m]], total))
  Tm <- lapply(seq_len(nv), function(m) {
    if (is.null(Pm[[m]])) Tm[[m]] else Tm[[m]] + Pm[[m]]
  })
  # the determinant's matrix moves by the derivative of ITS information, which
  # on the expected route is dE and not the observed third derivative Tm holds
  Km <- if (expected) {
    lapply(seq_len(nv), function(m) {
      TE <- contract3(spec, design, kE$deriv, params, npar, offs, total,
                      tv[[m]], key = kE$key)
      if (!is.null(Pm[[m]])) TE <- TE + Pm[[m]]
      pieces$S[[m]] + TE
    })
  } else {
    lapply(seq_len(nv), function(m) pieces$S[[m]] + Tm[[m]])
  }

  out <- matrix(0, nv, nv)
  for (m in seq_len(nv)) {
    for (l in m:nv) {
      # a pair touching an estimated coefficient has no penalty piece
      both_hyper <- m <= nh && l <= nh
      if (both_hyper) {
        key <- pieces$pair[m, l]
        Sml <- pieces$S2[[key]]
        cml <- pieces$c2[[key]]
        rho2 <- pieces$rho2[m, l]
      } else {
        Sml <- 0
        cml <- 0
        rho2 <- 0
      }
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
      if (length(hold)) rhs[hold] <- 0
      bml <- -msolve(rhs)
      if (length(hold)) bml[hold] <- 0
      # dK_m/dt_l enters ONLY through its trace against M, so the two
      # contractions are never assembled: tr(M X'WX) is a weighted sum of the
      # per-observation diagonal G. Measured at 8000 observations and 69
      # coefficients, forming the matrix and tracing costs 25.5 ms where the
      # sum costs 0.031 ms, and the pair loop does it twice per pair.
      tr_dKm <- if (expected) {
        sum(M * Sml) +
          trace_design_form(spec, G, kE2$deriv, params, npar, tv[[l]], tv[[m]],
                            key = kE2$key) +
          trace_design_form(spec, G, kE$deriv, params, npar,
                            block_predictors(design, params, npar, offs, bml),
                            key = kE$key)
      } else {
        sum(M * Sml) +
          trace_design_form(spec, G, d4, params, npar, tv[[l]], tv[[m]]) +
          trace_design_form(spec, G, d3, params, npar,
                            block_predictors(design, params, npar, offs, bml))
      }
      if (length(units)) {
        tr_dKm <- tr_dKm + sum(uref * bml) +
          trace_refresh4(spec, M, params, npar, Hl, dref[[l]], dref[[m]],
                         units, G, d3, bhat[[l]], f2, acurv)
      }
      E2 <- statmod_penalty_second(spec, design, coef, hyper, idx, m, l,
                                   bhat[[m]], bhat[[l]], bml, total)
      if (!is.null(E2)) tr_dKm <- tr_dKm + sum(M * E2)
      # the second derivative of l_p at the mode is -E_a'K_m b_b - r_a'b_b
      # - rho_ab; for two hyperparameters E vanishes and -r_m'b_l is
      # b_m'J b_l, the form kept here
      first <- if (both_hyper) {
        sum(bhat[[m]] * as.numeric(Jmat %*% bhat[[l]]))
      } else if (m > nh) {
        -sum(Hm[gam$where[m - nh], ] * bhat[[l]])
      } else {
        -sum(Hm[gam$where[l - nh], ] * bhat[[m]])
      }
      v <- -rho2 + first +
        sum((M %*% Km[[l]]) * t(M %*% Km[[m]])) / 2 -
        tr_dKm / 2
      out[m, l] <- v
      out[l, m] <- v
    }
  }

  # and onto the free scale the search runs on
  links <- attr(idx, "links")
  g <- statmod_marginal_grad(spec, design, coef, hyper, method, idx, basis,
                             free = FALSE, ctx = ctx, gam = gam)
  # an estimated coefficient is on its own scale already
  h1 <- rep(1, nv)
  h2 <- numeric(nv)
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
#' At order 1 the function calls [penalties7::penalty_dhessian()] and
#' [penalties7::penalty_cross()] only. The second-order generics
#' [penalties7::penalty_d2hessian()], [penalties7::penalty_dcross()] and
#' [penalties7::penalty_hess_theta()] are called at order 2, so a penalty that
#' supplies the first two generics gives an exact gradient without
#' implementing the others.
#'
#' A shared hyperparameter occupies one row of `idx` and receives the sum of
#' the contributions of its members. A penalty over the own parameters of a
#' structural term and a mixed covariance class contribute nothing here;
#' [structural_outer_pieces()] places them on the joint vector.
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
#' [statmod_marginal_hess()] uses this function for a model with a structural
#' term for which the joint assembly [statmod_structural_hess()] does not
#' apply: a filter without a penalty over its own parameters, a term of the
#' likelihood shape such as `regime()`, or a term that does not implement
#' [modelterms7::term_fourth()]. It also uses it on the expected information
#' where [outer_gradient_ok()] does not admit order 2. The analytic assembly
#' over the coefficients alone is not the criterion's Hessian when a filter's
#' own parameters move with the hyperparameter. [statmod_structural_hess()]
#' covers that case with the fourth derivative of the recursion from
#' [modelterms7::term_fourth()] and the fifth derivative of the family from
#' [distributions7::distrib_deriv5()].
#'
#' The function differences the exact gradient once, so the result contains a
#' single layer of finite differences. The matrix is symmetrized by averaging
#' it with its transpose.
#'
#' # Probes
#'
#' Each probe refits the coefficients from the coefficients supplied and from
#' the structural state as it stands, which is restored before each probe, so
#' that the two probes of a column differ in the hyperparameter alone. The
#' location error of the mode is nearly the same at \eqn{+h} and \eqn{-h} and
#' cancels in the difference. Where the inner optimizer is [iwls()], the probes
#' use its tolerance `min(tol, hess_stencil_inner_tol())`, a budget of at
#' least 500 iterations and no stopping criterion of its own; any other inner
#' optimizer is used as given.
#'
#' # Rejection of an unresolved stencil
#'
#' A stencil inherits the reproducibility of the quantity that it differences,
#' divided by the step, so where the mode is poorly located the result is
#' noise. The matrix is therefore computed at `h` and at `3 * h`, and the
#' function returns `NULL` where the largest absolute difference between the
#' two exceeds [hess_stencil_tol()] times the largest absolute entry of the
#' first. Callers treat a `NULL` result as an unavailable curvature.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients from which the probes start.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param idx The hyperparameter index, from [outer_hyper_index()].
#' @param basis The integrated basis, or `NULL`.
#' @param inner The inner optimizer with which the probes refit; `iwls()` where
#'   none is given.
#' @param h The step, on the free scale on which the search runs. For a
#'   coefficient that the criterion estimates the step is
#'   \eqn{h \max(1, |\gamma_j|)}.
#' @param gam The coefficients that the criterion estimates, from
#'   [marginal_coords()], or `NULL`. Each probe holds them at its own values
#'   and refits the rest, as [outer_fit()] does.
#'
#' @return The Hessian on the free scale over the hyperparameters followed
#'   by the coefficients in `gam`, or `NULL` where a probe could not be
#'   evaluated or the two steps disagree.
#'
#' @seealso [statmod_marginal_hess()], [statmod_marginal_grad()],
#'   [hess_stencil_step()]
#'
#' @keywords internal
statmod_hess_stencil <- function(spec, design, coef, hyper, method, idx,
                                 basis = NULL, inner = NULL,
                                 h = hess_stencil_step(), gam = NULL) {
  ng <- if (is.null(gam)) 0L else length(gam$where)
  if (!nrow(idx) && !ng) return(NULL)
  if (is.null(inner)) inner <- iwls()
  if (S7::S7_inherits(inner, Iwls)) {
    inner <- iwls_resolve(inner, spec@distrib)
    # the probes locate the mode far more tightly than a fit does. At the
    # fit's own tolerance a probe can stop where it started, the score there
    # already under the threshold, and the difference then misses the mode's
    # movement: a bias CONSTANT in the step, which the h-against-3h check
    # passes wherever both steps are small enough to share it. See
    # hess_stencil_inner_tol().
    inner@tol <- min(inner@tol, hess_stencil_inner_tol())
    inner@maxit <- max(inner@maxit, 500)
    inner@criterion <- NULL
  }
  cfg <- inner_settings(inner)
  blocks <- tryCatch(statmod_blocks(spec, design), error = function(e) NULL)
  if (is.null(blocks)) return(NULL)
  beta0 <- unlist(coef[spec@distrib@params], use.names = FALSE)
  # the outer vector is (eta, gamma): the coefficients the criterion
  # estimates are held at the probe's values, as outer_fit() holds them, and
  # their step is relative to their own size
  nh <- nrow(idx)
  eta0 <- c(hyper_to_eta(hyper, idx), if (ng) beta0[gam$where])
  nv <- length(eta0)
  scale <- c(rep(1, nh), if (ng) pmax(1, abs(beta0[gam$where])))
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
  grad_at <- function(v) {
    restore()
    hy <- eta_to_hyper(v[seq_len(nh)], idx, hyper)
    sp <- spec
    b0 <- beta0
    if (ng) {
      g <- v[nh + seq_len(ng)]
      sp <- statmod_hold(spec, gam, g)
      b0[gam$where] <- g
    }
    r <- tryCatch(statmod_alternate(sp, design, blocks, hy, inner, b0,
                                    cfg$expected, cfg$approx, cfg$maxit,
                                    cfg$tol, verbosity(0),
                                    hold_refresh = TRUE),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    u <- tryCatch(statmod_marginal_grad(sp, design, r$obj$split(r$par), hy,
                                        method, idx, basis, gam = gam),
                  error = function(e) NULL)
    if (is.null(u) || length(u) != nv || !all(is.finite(u))) NULL else u
  }
  at_step <- function(step) {
    H <- matrix(0, nv, nv)
    for (m in seq_len(nv)) {
      s <- step * scale[[m]]
      ep <- eta0; ep[[m]] <- ep[[m]] + s
      em <- eta0; em[[m]] <- em[[m]] - s
      gp <- grad_at(ep)
      gm <- grad_at(em)
      if (is.null(gp) || is.null(gm)) return(NULL)
      H[, m] <- (gp - gm) / (2 * s)
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
#' The step at which [statmod_hess_stencil()] differences, the largest
#' relative disagreement between its two readings, and the tolerance of the
#' probes' inner fit.
#'
#' @details
#' `hess_stencil_step()` is `1e-3`, on the free scale of the hyperparameters.
#' The default of [numericals7::fd_step()], about 6e-6, is intended for a
#' function evaluated to machine precision. Here the gradient is computed by
#' refitting a mode, and the reproducibility of that refit bounds the step from
#' below.
#'
#' `hess_stencil_tol()` is `1e-3`. It is compared with the largest absolute
#' difference between the matrices at `h` and `3 * h`, divided by the largest
#' absolute entry of the first. A relative error of `1e-3` in the curvature
#' corresponds to a relative error of about `5e-4` in a standard error.
#'
#' `hess_stencil_inner_tol()` is `1e-10`, the largest tolerance that the probes
#' use when the inner optimizer is [iwls()]. At the tolerance of a fit a probe can
#' stop at its starting point, and the difference then misses the movement of
#' the mode. The resulting bias does not grow as \eqn{1/h}, so the comparison
#' of the two steps does not detect it.
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

#' @rdname hess_stencil_step
#' @keywords internal
hess_stencil_inner_tol <- function() 1e-10


#' The Key of a Hyperparameter Pair
#'
#' @description
#' Locates the entry of a penalty's hyperparameter Hessian, which is keyed by
#' name, never by position.
#'
#' @param a,b The two hyperparameter names.
#' @param keys The names that the penalty returned.
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
#' @return A list with one element for each distribution parameter: a numeric
#'   vector with one value per observation, or `NULL` for a parameter whose
#'   equation has no coefficients.
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
#' Each block is minus a weighted crossproduct of the two designs, the weight
#' being \eqn{w_i\sum_k \ell'''_{abk}(X_k v_k)_i}, so that the third
#' derivative is never formed as an array over the coefficients. The sign makes
#' \eqn{T[v]} the derivative of \eqn{K = -\ell'' + S}.
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
#' @details
#' The blocks follow [contract3()]: each is minus a weighted crossproduct of the
#' two designs, with the weight formed from the fourth derivatives and the
#' predictors of the two directions.
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
#' @param keys The names that the derivative returned.
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
#' The ordinary effective degrees of freedom treat the smoothing parameters as
#' known, although they were chosen from the same data. Propagating their
#' uncertainty into the coefficients gives the corrected Bayesian covariance
#'
#' \deqn{V' = V_\beta + J V_\theta J^\top, \qquad
#'   J = \partial\hat\beta/\partial\theta,}
#'
#' and the corrected count is \eqn{\mathrm{tr}(V' H)}. \eqn{J} comes from
#' the implicit function theorem at the penalized mode: differentiating
#' \eqn{\partial(-\ell + \rho)/\partial\beta = 0} in the hyperparameter
#' gives \eqn{(H+S)J_k = -\partial^2\rho/\partial\beta\,\partial\theta_k},
#' whose right-hand side is \pkg{penalties7}'s `penalty_cross()`.
#' \eqn{V_\theta} is the inverse of the negative of the outer criterion's
#' Hessian, which [statmod_marginal_hess()] returns; the negative is taken
#' because the criterion is maximized. The Hessian is that of the criterion
#' the fit maximized, so for [ml()] the determinant spans the range space of
#' the penalty ([integrated_basis()]).
#'
#' Everything is on the hyperparameter's link scale, which is the scale on
#' which the outer criterion is optimized and therefore the only scale on which
#' its Hessian is a variance.
#'
#' **Relation to mgcv.** The correction corresponds to the term of the
#' effective degrees of freedom `edf2` of mgcv that arises from the uncertainty
#' of the smoothing parameters (Wood, Pya and Safken, 2016). The scale of a
#' Gaussian response is a distribution parameter with its own coefficients here,
#' so its uncertainty is part of \eqn{H} and needs no separate term.
#'
#' **Where it does not apply.** The function returns zero, with `n_hyper` zero,
#' for every outer criterion other than [reml()] and [ml()], including [aic()]
#' and [bic()] used for a smooth penalty. A kinked penalty has no hyperparameter
#' that the outer criterion estimates, because [outer_hyper_index()] skips it,
#' so there is nothing to propagate and the correction is zero. The map from
#' the hyperparameter to the penalized mode turns a corner whenever a
#' coefficient joins or leaves the active set, and a delta method needs a
#' derivative that does not exist there.
#'
#' **Coefficients estimated on the criterion.** Where the criterion also
#' estimates coefficients ([marginal_coords()]), the function returns
#' [marginal_edf_correction()].
#'
#' **A model carrying a filter.** The mode of such a model moves in the joint
#' vector of the coefficients and the own parameters of the structural term,
#' which is the vector that the determinant of the criterion spans. A penalty
#' over those parameters is not a column of any design, so [hyper_mode_cross()]
#' omits it when it is called for the coefficients alone. The correction is
#' therefore computed on the joint vector.
#' [statmod_edf()] takes the effective degrees of freedom of a structural
#' model from `joint_smoother_diag()`, which reads the same two matrices, so
#' the base count and the correction are both computed on the joint vector.
#' The two matrices are those that [statmod_marginal_full()] and
#' [statmod_full_information()] build, so that this function and the
#' criterion use one assembly. The joint route reads the observed information,
#' and `expected` and `approx` do not reach it. A model with no structural term
#' of the filter shape takes the coefficient route, because
#' [statmod_marginal_full()] returns `NULL` there.
#'
#' @param spec A [StatmodSpec()].
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param design The design.
#' @param method The outer method that estimated them, or `NULL`.
#' @param expected Whether the information \eqn{H} is the expected one. Not
#'   read where the model carries a filter, whose information is observed. The
#'   caller passes the information that the fit used, so that \eqn{V_\beta},
#'   \eqn{H} and the criterion's Hessian read one information.
#' @param approx The approximation for the expected information. Not read
#'   where the model carries a filter.
#' @param pinned The coordinates that the fit left out of the criterion's
#'   determinant, as it records them in `methods$pinned`.
#'
#' @return A list with `total`, the scalar correction, `per`, one entry for
#'   each parameter and term, and `n_hyper`, the number of hyperparameters that
#'   were estimated. `per` uses the diagonal of \eqn{V_\theta} alone, so its
#'   sum differs from `total` when the estimated hyperparameters are
#'   correlated. The total and `n_hyper` are zero and `per` is empty where no
#'   hyperparameter was estimated. A zero `total` beside a positive `n_hyper`
#'   means that the curvature could not be read.
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
                                   expected = TRUE, approx = "opg",
                                   pinned = NULL) {
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
  # where the criterion estimated coefficients, the mode moves with a
  # hyperparameter at those coefficients held, and the hyperparameter's
  # variance is read over (eta, gamma)
  mc <- fit_marginal_context(spec, design, coef, method, pinned)
  if (length(mc$gam$where)) {
    return(marginal_edf_correction(mc, coef, hyper, method, idx, expected,
                                   approx, zero))
  }

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

  # the curvature of the criterion the fit maximized: under ml() the
  # determinant spans the penalty's range space, and NULL is reml()'s
  Ho <- tryCatch(statmod_marginal_hess(spec, design, coef, hyper, method,
                                       idx, integrated_basis(spec, design,
                                                             method@kind)),
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
#' This function computes the cross derivative of the penalty that the
#' movement of the penalized mode requires. Three functions read it here
#' instead of assembling it separately: [statmod_edf_correction()], which
#' contracts the movement with the information to obtain the additional
#' effective degrees of freedom, [hyper_correction()], which adds the resulting
#' matrix to a variance, and [marginal_blocks()], which builds the movement of
#' the integrated coefficients.
#'
#' A shared hyperparameter is one column that stands for several penalties, so
#' each member writes into the column of its group and the contributions
#' accumulate. Where nothing is shared, each member has its own column.
#'
#' A penalty over the own parameters of a structural term covers positions among
#' those parameters, which are not columns of a design. Where the matrix spans
#' the coefficients alone, there is no position to write to and the penalty is
#' skipped. The skipped penalties are counted in `skipped`, because a
#' correction assembled without them is incomplete.
#'
#' With `joint` the matrix spans the vector in which the mode of such a model
#' moves, the coefficients followed by the free parameters of the term, and
#' the penalties over those parameters have rows. [unit_joint_positions()]
#' returns the positions of the coordinates of each unit in that vector and
#' [unit_joint_beta()] returns their values. [hyper_correction()] and
#' [statmod_edf_correction()] pass `joint = TRUE` for a model that carries a
#' filter, because [vcov.StatmodFit()] and [statmod_edf()] use the joint
#' penalized information. [marginal_blocks()] uses the matrix over the
#' coefficients alone.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param idx The hyperparameter index, from [outer_hyper_index()].
#' @param n How many rows the matrix carries: the stacked coefficients, and
#'   with `joint` a structural term's free parameters after them.
#' @param joint Whether those rows include that tail. `FALSE`, the default,
#'   gives the matrix over the coefficients alone.
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
#' The matrices that [vcov.StatmodFit()] returns for `type = "bayesian"` and
#' `type = "frequentist"` are conditional on the hyperparameters: each is read
#' at the value at which the outer search stopped, as though that value were
#' known. The value was estimated from the same data, and the cost of that is
#' \deqn{V' = V_b + J V_\theta J', \qquad
#'   J = -(H + S)^{-1} \frac{\partial^2 \rho}{\partial\beta \partial\theta},}
#' the delta method applied to the map from the hyperparameter to the mode
#' (Wood, Pya and Safken, 2016). The same quantity is contracted against the
#' information by [statmod_edf_correction()] to obtain a count of parameters;
#' this function returns the matrix.
#'
#' \eqn{V_b} is passed in and is not recomputed. The caller has already settled
#' which information \eqn{H} is, which coordinates are held and which are
#' aliased, and a correction built on a second inverse, regularized
#' differently, would not describe the movement of the mode to whose variance
#' it is added.
#'
#' # Availability of the correction
#'
#' The function returns `C = NULL` with `n_hyper = 0` where no hyperparameter
#' was estimated by [reml()] or [ml()]. That includes a fit whose outer
#' criterion is [aic()] or [bic()], and a kinked penalty, whose hyperparameter
#' is the argument of a minimum over a grid that [outer_hyper_index()] skips.
#' The map from that hyperparameter to the mode turns a corner whenever a
#' coefficient joins or leaves the active set, so no derivative exists to
#' propagate.
#'
#' The correction is unavailable, with `n_hyper` positive and `C` `NULL`, where
#' the criterion's own Hessian cannot be read, or where [hyper_variance()]
#' cannot invert it because the coupling of a held coordinate to the kept ones
#' is too large to ignore.
#'
#' The correction is partial, with `complete` false, where some of it could be
#' read and some could not: a hyperparameter that [hyper_variance()] held
#' contributes nothing, and a penalty for which [hyper_mode_cross()] could not
#' return a column (counted in its `skipped`) contributes nothing. The matrix
#' returned then omits those contributions.
#'
#' A penalty over the own parameters of a structural term does not make the
#' correction partial. The mode of such a model moves in the joint vector of
#' the coefficients and the free parameters of the term, which `Vb` spans, and
#' [hyper_mode_cross()] is called with `joint = TRUE` to return rows for the
#' whole vector.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param hyper The hyperparameters.
#' @param method The outer method that estimated them, or `NULL`.
#' @param Vb The bayesian variance to which the correction is added, over the
#'   coordinates that `keep` marks.
#' @param keep A logical vector over the joint vector (the stacked
#'   coefficients followed by any structural tail) that marks the coordinates
#'   `Vb` spans. It is the whole vector of the caller and not its head,
#'   because a tail coordinate dropped as a flat direction leaves `Vb` narrower
#'   than the count alone would indicate.
#' @param nz How many of those coordinates are the structural tail.
#'
#' @return A list with `C`, the correction or `NULL`, `n_hyper`, the number of
#'   hyperparameters that [reml()] or [ml()] estimated, and `complete`,
#'   whether every one of them contributed and no penalty was skipped.
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
                                       idx, integrated_basis(spec, design,
                                                             method@kind)),
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
#' Such a coordinate is held, and the remaining coordinates are inverted. The
#' result is the variance conditional on the held coordinates, which is the
#' reading that [vcov.StatmodFit()] gives to the remaining coefficients when
#' it leaves out a coefficient about which the information is flat. The
#' conditional variance equals the marginal variance only where the coupling
#' contributes nothing to the curvature of the kept coordinates. The Schur
#' correction \eqn{A_{kb}A_{bb}^{-1}A_{bk}} is therefore computed and compared
#' with the kept diagonal. Above `schur` times the kept diagonal, the function
#' returns `NULL`.
#'
#' A variance is returned only where the inverted matrix is positive definite,
#' which is decided from the smallest eigenvalue of the matrix equilibrated to
#' a unit diagonal. A positive diagonal of the inverse does not suffice,
#' because an indefinite matrix can have one.
#'
#' A coordinate that the certificate places at the edge of its range is held
#' whatever its curvature, and the Schur test is not applied to it, because
#' the criterion is flat in that coordinate and its curvature and its coupling
#' to the other coordinates are below the resolution of the criterion.
#'
#' @param A The negative of the outer Hessian, with dimnames.
#' @param schur The largest ratio of the Schur correction to the diagonal entry
#'   of a kept coordinate. The default is `1e-4`.
#' @param at_edge Positions of coordinates at the edge of their range, held
#'   without the Schur test.
#'
#' @return A matrix of the same shape as `A` with the variance in the kept
#'   rows and columns and `NA` elsewhere, or `NULL` when the matrix is not
#'   positive definite and no coordinate can be held to make the kept part so,
#'   when no coordinate remains, or when the Schur test fails.
#'
#' @seealso [statmod_hyper_vcov()], [hyper_correction()], [marginal_vcov()] and
#'   [marginal_edf_correction()], which also call it.
#'
#' @keywords internal
hyper_variance <- function(A, schur = 1e-4, at_edge = integer(0)) {
  # positive definite on the equilibrated scale, so the test does not turn on
  # how different the hyperparameters' scales are
  pd <- function(M) {
    if (!all(is.finite(M)) || any(diag(M) <= 0)) return(FALSE)
    s <- 1 / sqrt(diag(M))
    ev <- eigen(s * t(s * M), symmetric = TRUE, only.values = TRUE)$values
    min(ev) > ncol(M) * .Machine$double.eps
  }
  usable <- function(M) !is.null(M) && all(is.finite(M)) && all(diag(M) > 0)
  p <- ncol(A)
  at_edge <- intersect(as.integer(at_edge), seq_len(p))
  if (!length(at_edge)) {
    V <- if (pd(A)) tryCatch(solve(A), error = function(e) NULL) else NULL
    if (usable(V)) return(V)
  }
  d <- diag(A)
  bad <- setdiff(which(!is.finite(d) | d <= 0 |
                         apply(!is.finite(A), 1L, any)), at_edge)
  keep <- setdiff(seq_len(p), c(bad, at_edge))
  if (!length(c(bad, at_edge)) || !length(keep)) return(NULL)
  if (!pd(A[keep, keep, drop = FALSE])) return(NULL)
  W <- tryCatch(solve(A[keep, keep, drop = FALSE]), error = function(e) NULL)
  if (!usable(W)) return(NULL)
  if (length(bad)) {
    # the Schur test over the coordinates held for their curvature alone; a
    # coordinate at the edge has no coupling the criterion can resolve
    Bi <- tryCatch(solve(A[bad, bad, drop = FALSE]), error = function(e) NULL)
    if (is.null(Bi) || any(!is.finite(Bi))) return(NULL)
    Akb <- A[keep, bad, drop = FALSE]
    corr <- diag(as.matrix(Akb %*% Bi %*% t(Akb)))
    if (any(abs(corr) > schur * d[keep])) return(NULL)
  }
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
#' A hyperparameter estimated by [reml()] or [ml()] is the maximizer of a
#' criterion that is twice differentiable in it, so its variance is the
#' inverse of the negative Hessian of that criterion at the point reached,
#' which [statmod_marginal_hess()] returns. If the criterion also estimates
#' coefficients, the variance is the block of the hyperparameters in the
#' inverse of the negative Hessian over the hyperparameters and those
#' coefficients. [statmod_marginal_hess()] computes the Hessian analytically
#' where a route applies and by a difference of the exact gradient otherwise.
#' The matrix is read on the free scale because the criterion was maximized on
#' that scale, where the quadratic approximation is reasonable. A variance on
#' the scale of the hyperparameter itself would describe a shape that the
#' criterion does not have where the estimate is far from zero and the
#' criterion is asymmetric.
#'
#' **Where it does not apply.** The function returns `NULL` for every outer
#' criterion other than [reml()] and [ml()], including [aic()] and [bic()]
#' used as the outer criterion for a smooth penalty. A hyperparameter that a
#' path selects (the penalty of a kinked term, chosen by [aic()], [bic()] or
#' [cv()] over a grid) is the argument of a minimum, not a root of a
#' derivative, so no Hessian exists to invert. Its uncertainty can be assessed
#' by resampling.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method The outer method that estimated them, or `NULL`.
#' @param inner The inner optimizer that the fit used, with which the stencil
#'   route refits its probes; `iwls()` where none is given.
#' @param pinned The coordinates that the fit left out of the criterion's
#'   determinant, as it records them in `methods$pinned`.
#' @param at_edge The `boundary_key` of [statmod_certificate()]: the
#'   hyperparameters at the edge of their range, held by [hyper_variance()]
#'   so that the others carry the variance conditional on them.
#'
#' @return A square matrix, one row per estimated hyperparameter, whose
#'   dimnames join the distribution parameter, the term and the
#'   hyperparameter's own name with a carriage return, a character that no
#'   name can contain, so that the three are recoverable from the key. The
#'   index is stored in the attribute `"idx"`. The rows and columns of held
#'   coordinates are `NA` (see [hyper_variance()]). The value is `NULL` where
#'   there is nothing to report.
#'
#' @seealso [statmod_marginal_hess()], [summary.StatmodFit()]
#'
#' @keywords internal
statmod_hyper_vcov <- function(spec, design, coef, hyper, method,
                               inner = NULL, pinned = NULL,
                               at_edge = character(0)) {
  if (is.null(method) || !method@kind %in% c("ml", "reml")) return(NULL)
  blocks <- statmod_blocks(spec, design)
  idx <- outer_hyper_index(spec, blocks)
  if (!nrow(idx)) return(NULL)
  # on the criterion the fit maximized: where it estimated coefficients too,
  # the hyperparameters' variance is the eta block of the inverse over
  # (eta, gamma), so the coefficients' own uncertainty reaches it
  mc <- fit_marginal_context(spec, design, coef, method, pinned)
  ng <- length(mc$gam$where)
  if (ng) {
    basis <- integrated_basis(mc$spec, mc$design, method@kind, gamma = TRUE)
    Ho <- tryCatch(statmod_marginal_hess(mc$spec, mc$design, coef, hyper,
                                         method, idx, basis, inner = inner,
                                         gam = mc$gam),
                   error = function(e) NULL)
    if (is.null(Ho) || !all(is.finite(Ho))) return(NULL)
    k <- paste(idx$parameter, idx$term, idx$name, sep = "\r")
    Vv <- hyper_variance(-as.matrix(Ho), at_edge = which(k %in% at_edge))
    if (is.null(Vv)) return(NULL)
    V <- Vv[seq_len(nrow(idx)), seq_len(nrow(idx)), drop = FALSE]
    dimnames(V) <- list(k, k)
    return(structure(V, idx = idx))
  }
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
  V <- hyper_variance(A, at_edge = which(k %in% at_edge))
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
#' The gradient uses only the trace of this matrix against \eqn{M} and does
#' not form the matrix. The Hessian uses the matrix in two places: in
#' \eqn{\mathrm{tr}(M K_l M K_m)}, which no contraction reduces, and as the
#' operator that carries the second movement of the mode. Write \eqn{V_a} for
#' the rows of each equation, \eqn{\mathrm{d}\varphi = E v} for the movement of
#' the filter's own parameters, and \eqn{E} for the second derivative of the
#' predictor that the filter produces. Then
#'
#' \deqn{\frac{\partial K}{\partial u}[v] =
#'   -\sum_i w_i\Big[\sum_{a,b}\Big(\sum_k\ell_{abk}(V_k\cdot v)\Big)
#'     V_a^\top V_b
#'   + \sum_b \ell_{pb}\big(\mathrm{d}\varphi\otimes V_b
#'     + V_b\otimes\mathrm{d}\varphi\big)\Big]
#'   - W(\kappa_v) - W_3[v],}
#'
#' with \eqn{\kappa_v = \sum_k \ell_{pk}(V_k\cdot v)} the weight at which the
#' term of the level is read again and \eqn{W_3[v]} the curvature that
#' [modelterms7::term_third()] returns. The trace of the matrix against
#' \eqn{M} equals the value of [structural_chain_extra()] up to rounding.
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
    function(e, i) st$sd_at[i], function(e, i) st$cd_at[i], f$psi,
    w * st$s_at, st$seed, st$blocks(vfull), vfull,
    score_values = st$sd_at, curvature_values = st$cd_at,
    blocks_data = st$blocks_data, threads = spec@threads)
  dphi <- cv3$dphi[, keep, drop = FALSE]
  out <- -cv3$curvature[keep, keep, drop = FALSE]

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
#' The terms are listed in the comments of the code. None of them assembles a
#' matrix over the coefficients: a term in \eqn{V_a^\top V_b} traces as a
#' weighted sum of the per-observation diagonal \eqn{G}, a term carrying
#' \eqn{\mathrm{d}\varphi} traces against the rows to which \eqn{M} has
#' already been applied, and the terms that carry the derivatives of the
#' recursion trace against what [modelterms7::term_curvature()],
#' [modelterms7::term_third()] and [modelterms7::term_fourth()] return at the
#' appropriate weights.
#'
#' The last of those is the fourth derivative of the predictor through the
#' recursion. The fifth derivative of the family enters only through `P` in
#' the `blocks` callback of that call.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param jd The joint rows.
#' @param M The matrix against which the trace is taken.
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
    function(e, i) st$sd_at[i], function(e, i) st$cd_at[i], f$psi,
    wt * st$s_at, st$seed, blk4(list(vfull, wfull)), list(vfull, wfull),
    score_values = st$sd_at, curvature_values = st$cd_at,
    blocks_data = attr(blk4, "data"), threads = spec@threads)
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
    function(e, i) st$sd_at[i], function(e, i) st$cd_at[i], f$psi,
    wt * gk, st$seed, st$blocks(NULL),
    score_values = st$sd_at, curvature_values = st$cd_at,
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
      function(e, i) st$sd_at[i], function(e, i) st$cd_at[i], f$psi,
      wt * pair[[1L]], st$seed, st$blocks(pair[[2L]]), pair[[2L]],
      score_values = st$sd_at, curvature_values = st$cd_at,
      blocks_data = st$blocks_data, threads = spec@threads)
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
#' [outer_pieces()] for the vector that consists of the coefficients and the
#' own parameters of a filter. The determinant of a marginal criterion spans
#' that vector where a penalty covers those parameters.
#'
#' @details
#' The arithmetic is that of [outer_pieces()]; the difference is where the
#' coordinates of each unit are placed. An ordinary unit is matched into the
#' kept coefficients, a penalty over the own parameters of a structural term is
#' placed among the free parameters of the filter, and a mixed covariance
#' class is read where each of its coordinates lives and put back in the
#' interleaved order of the class. These three readings are those of
#' [statmod_structural_grad()], and they are reused here so that the gradient
#' and the Hessian place each coordinate at the same position.
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
#' # Orders of differentiation
#'
#' Each order of differentiating the predictor through the recursion requires
#' one more order of the response's family, because the score that drives the
#' recursion is read at the predictor that the recursion produces. The gradient
#' reads \eqn{\partial^3 e/\partial u^3} in one direction through
#' [modelterms7::term_third()]. The second derivative of the criterion reads
#' \eqn{\partial^4 e/\partial u^4} in two directions through
#' [modelterms7::term_fourth()], and the fifth derivative of the family, from
#' [distributions7::distrib_deriv5()], enters with it.
#'
#' # Form of the assembly
#'
#' The assembly has the form of [statmod_marginal_hess()]. With \eqn{u} the
#' joint vector, \eqn{K} the penalized information over it and \eqn{M} the
#' matrix against which the trace is taken,
#'
#' \deqn{\frac{\partial^2 V}{\partial t_m\partial t_l}
#'   = -\rho_{ml} + \hat b_m^\top K \hat b_l
#'   + \tfrac{1}{2}\mathrm{tr}(MK_lMK_m)
#'   - \tfrac{1}{2}\mathrm{tr}\Big(M\frac{\partial K_m}{\partial t_l}\Big),}
#'
#' with \eqn{\hat b_m = -K^{-1}c_m} the mode's movement, \eqn{K_m = S_m +
#' \partial K/\partial u[\hat b_m]} and the last trace carrying the penalty's
#' second derivative, the twice-contracted second derivative of \eqn{K} and
#' the once-contracted one at the second movement of the mode. Every piece is
#' the one that the coefficient-space assembly uses, read on the joint vector.
#'
#' # Blocks that move with their coefficients
#'
#' A block that moves with its coefficients (`nl()`, `seg()`) beside the filter
#' contributes nothing here, as it contributes nothing to
#' [statmod_structural_grad()]: the correction for such a block is written in
#' the coefficient-space assembly and has no counterpart on the joint vector.
#' Such a model is admitted at both orders, and the result has the
#' approximation of the gradient.
#'
#' # Cost
#'
#' The assembly makes one call of [modelterms7::term_fourth()], together with
#' calls of the lower-order recursions, for each pair of hyperparameters. The
#' stencil refits the model four times for each hyperparameter.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param idx The outer index.
#' @param basis The integrated subspace, or `NULL`.
#'
#' @return A square matrix on the free scale with one row for each row of
#'   `idx`, or `NULL` where the joint matrix cannot be formed or factorized or
#'   the fifth derivative of the family is unavailable.
#'
#' @seealso [statmod_structural_grad()], [statmod_marginal_hess()],
#'   [statmod_hess_stencil()] for the route that it replaces.
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
  dr4 <- filter_driving(spec, jd$ev$theta, jd$ap, filter_scaling(jd$f$tm),
                        st$gl, st$H, st$D3, st$D4, D5)
  blk4 <- .structural_blocks(jd$params, jd$ap, jd$V, dr4$H, dr4$D3, dr4$D4,
                             jd$n, dr4$D5)
  # the same pieces as data, for the compiled fourth order
  attr(blk4, "data") <- structural_blocks_data(jd$params, jd$ap, jd$V, dr4$H,
                                               dr4$D3, jd$n, D4 = dr4$D4,
                                               D5 = dr4$D5)

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
