#' @include outer_gradient.R
NULL

# How a block that moves with its coefficients enters the outer derivatives.
#
# With H_{(a,j),(b,k)} = -sum_i w_i l_ab,i X_a[i,j] X_b[i,k] and a block whose
# X is the Jacobian of the term's own contribution, differentiating in beta
# gives three contributions rather than one: the log-likelihood's own third
# derivative with X held, and two from dX/dbeta that are transposes of each
# other under the trace. The gradient consumes them through `u_refresh()`; the
# HESSIAN needs the same derivative in three further places, and the three are
# written here because correcting one of them alone is worse than correcting
# none -- measured, a weakly identified nl goes from 2.02e-01 to 6.35e+02 with
# a standard error of NaN when only dK/dt is corrected.
#
# The HESSIAN's twice-contracted term needs two further pieces, and they are
# written together for the same reason: measured against a mixed second
# difference of H, adding either one alone leaves a curved nl WORSE or barely
# better (2.62e-02 today, 2.88e-02 with the first alone, 2.65e-03 with the
# second alone, 5.19e-08 with both).
#
# The first is d2X/dbeta2, the term's own third derivative, which
# term_block_deriv2() supplies. The second reads only dX/dbeta and had been
# missing since these corrections were written: the third-derivative
# contraction carries the DIRECTION'S PREDICTOR (X_m v_m), and where X_m is a
# Jacobian that predictor moves with beta too. On a bilinear f, where
# d2X/dbeta2 is exactly zero, it is the whole of the gap -- 2.32e-02 against
# an arbitrary trace matrix, 2.19e-08 once it is in.


#' The Terms Whose Block Moves With Their Coefficients
#'
#' @description
#' Returns one entry per refreshable term, carrying everything that its
#' derivative needs: the distribution parameter of the equation in which the
#' term sits, its columns within that equation and in the stacked coefficient
#' vector, its own coefficients and the term rebuilt at them.
#'
#' @details
#' The entries are resolved once and passed to every consumer, so the
#' gradient's contraction and the Hessian's three corrections cannot disagree
#' about which terms move or where their columns are.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, already refreshed at `coef`.
#' @param coef The coefficients.
#' @param params,npar,offs The block bookkeeping.
#'
#' @return A list, which is empty if the design has no block that moves. Each
#'   element is a list with `a`, the index of the term's equation; `param`, the
#'   distribution parameter of that equation; `term`, the name of the term
#'   within the equation; `cols`, its columns within the equation; `ra`, the
#'   same columns in the stacked coefficient vector; `bt`, the term's
#'   coefficients; and `built`, the term refreshed at `bt`.
#'
#' @seealso [u_refresh()], [contract3_refresh()]
#'
#' @keywords internal
refresh_units <- function(spec, design, coef, params, npar, offs) {
  rf <- attr(design, "refresh")
  st <- attr(design, "state")
  if (is.null(rf) || !length(rf) || is.null(st)) return(list())
  out <- list()
  for (r in rf) {
    p <- r$param
    a <- match(p, params)
    if (is.na(a) || npar[a] == 0L) next
    cols <- design[[p]]$blocks[[r$term]]
    if (!length(cols)) next
    tm <- st$terms[[p]][[r$term]]
    if (is.null(tm)) next
    bt <- coef[[p]][cols]
    out[[length(out) + 1L]] <- list(
      a = a, param = p, term = r$term, cols = cols,
      ra = offs[a] + cols, bt = bt,
      built = modelterms7::term_refresh(tm, bt))
  }
  out
}


#' The Constraints That the Kinks Place on the Stacked Coefficients
#'
#' @description
#' Reads, through [modelterms7::term_kink_constraints()], the linear
#' constraints that hold the objective on its kinks at the coefficients given,
#' for every term that recomputes its own block, and places their rows in the
#' columns of the stacked vector.
#'
#' @details
#' A break-point that sits on an observation is a kink of the objective. Where
#' its position is a single coefficient the row is that coefficient's unit
#' vector; where the position is developed over a design the row is the
#' development's row of that observation, and the coefficients can still move
#' along the hyperplane on which the break-point stays on the observation.
#' [iwls_fit()] solves on that subspace where a line search rejects every step
#' length, and [mode_solver()] moves the mode in it.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients, a named list with one numeric vector per
#'   parameter.
#' @param split The stacked positions of each parameter's coefficients, a
#'   named list as the objective's `split` returns it for `seq_along()` of the
#'   stacked vector.
#' @param total The length of the stacked vector.
#'
#' @return A numeric matrix with `total` columns and one row per constraint,
#'   with no rows if the design has no term with a kink.
#'
#' @seealso [iwls_fit()], [mode_solver()], [refresh_units()].
#'
#' @keywords internal
kink_constraints <- function(spec, design, coef, split, total) {
  rf <- attr(design, "refresh")
  st <- attr(design, "state")
  out <- matrix(0, 0, total)
  if (is.null(rf) || !length(rf) || is.null(st)) return(out)
  for (r in rf) {
    p <- r$param
    cols <- design[[p]]$blocks[[r$term]]
    tm <- st$terms[[p]][[r$term]]
    if (!length(cols) || is.null(tm)) next
    A <- modelterms7::term_kink_constraints(tm, coef = coef[[p]][cols])
    if (!nrow(A)) next
    M <- matrix(0, nrow(A), total)
    M[, split[[p]][cols]] <- A
    out <- rbind(out, M)
  }
  out
}


#' A Basis of the Directions That a Set of Linear Constraints Leaves Free
#'
#' @description
#' Returns a `p x r` matrix \eqn{N} with \eqn{B N = 0} and
#' \eqn{N^\top N = I}, whose columns span the null space of `B`.
#'
#' @details
#' Where every row of `B` is a multiple of a unit vector, the constraints hold
#' coordinates, and \eqn{N} is the identity's columns of the other
#' coordinates, so the reduced problem is the problem without those
#' coordinates exactly. Otherwise \eqn{N} is read off the complete QR
#' decomposition of \eqn{B^\top}.
#'
#' @param B A numeric matrix with `p` columns.
#' @param p The number of columns.
#'
#' @return The matrix \eqn{N}, with zero columns where `B` has full column
#'   rank.
#'
#' @keywords internal
null_basis <- function(B, p = ncol(B)) {
  B <- as.matrix(B)
  if (!nrow(B)) return(diag(p))
  nz <- rowSums(B != 0)
  B <- B[nz > 0, , drop = FALSE]
  if (!nrow(B)) return(diag(p))
  if (all(rowSums(B != 0) == 1L)) {
    held <- unique(apply(B != 0, 1L, which))
    return(diag(p)[, setdiff(seq_len(p), held), drop = FALSE])
  }
  q <- qr(t(B))
  if (q$rank >= p) return(matrix(0, p, 0L))
  Q <- qr.Q(q, complete = TRUE)
  Q[, (q$rank + 1L):p, drop = FALSE]
}


#' The Movement of the Mode Under Held Coordinates and Kinks
#'
#' @description
#' Builds the solve that gives the movement of the penalized mode under a
#' perturbation \eqn{c}: \eqn{v = -N (N^\top (K + D_m) N)^{-1} N^\top c}, with
#' \eqn{N} a basis of the directions that the held coordinates and the kinks
#' leave free.
#'
#' @details
#' The held coordinates are those pinned in \eqn{K} other than at a boundary
#' (a coefficient that the criterion estimates, a direction that the model
#' does not identify, a coefficient that a kinked penalty holds at zero). The
#' kinks are the break-points that sit on an observation at the mode
#' ([kink_constraints()]): there the one-sided slopes have opposite signs, so
#' a small perturbation leaves the break-point on the observation, as a lasso
#' leaves a coefficient at zero. Where the break-point's position is a single
#' coefficient, that coefficient does not move; where it is developed over a
#' design, the mode moves along the hyperplane on which the break-point stays
#' on the observation. The solve is the bordered system of the equality
#' constrained problem, written in the basis of its null space. The
#' criterion's own determinant is not affected: it keeps the side convention
#' of the break-point term.
#'
#' Where \eqn{N^\top (K + D_m) N} is not positive definite, the solve uses
#' \eqn{K} in place of \eqn{K + D_m}.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design at the mode.
#' @param coef The coefficients at the mode, a named list.
#' @param K The mode's penalized information, as [ctx_penalized()] returns
#'   it with `laplace = FALSE`.
#' @param Kinv Its inverse.
#' @param Dm The block's own curvature from [mode_curvature()], or `NULL`.
#' @param npar,offs The coefficient counts and offsets per parameter.
#'
#' @return A list with `solve`, a function of \eqn{c} that returns \eqn{v}
#'   without the sign; `J`, the curvature used in the solve; `Dm`, the
#'   curvature in which the held rows and columns are set to zero; and
#'   `hold`, the held positions.
#'
#' @seealso [statmod_marginal_grad()], [statmod_marginal_hess()],
#'   [null_basis()]
#'
#' @keywords internal
mode_solver <- function(spec, design, coef, K, Kinv, Dm, npar, offs) {
  total <- sum(npar)
  split <- stats::setNames(lapply(seq_along(npar), function(a)
    offs[a] + seq_len(npar[a])), names(design))
  Ak <- kink_constraints(spec, design, coef, split, total)
  hold <- setdiff(attr(K, "held_at"), attr(K, "boundary_at"))
  if (!is.null(Dm) && length(hold)) {
    Dm[hold, ] <- 0
    Dm[, hold] <- 0
  }
  J <- K
  solve <- function(z) as.numeric(Kinv %*% z)
  kinked <- nrow(Ak) > 0L
  if ((!is.null(Dm) && any(Dm != 0)) || kinked) {
    Jt <- as_dense(K) + if (is.null(Dm)) 0 else Dm
    if (kinked) {
      E <- matrix(0, length(hold), total)
      if (length(hold)) E[cbind(seq_along(hold), hold)] <- 1
      N <- null_basis(rbind(Ak, E), total)
      reduced <- function(Jm) {
        fac <- tryCatch(chol(crossprod(N, Jm %*% N)), error = function(e) NULL)
        if (is.null(fac)) return(NULL)
        function(z) {
          r <- crossprod(N, as.numeric(z))
          as.numeric(N %*% backsolve(fac, forwardsolve(t(fac), r)))
        }
      }
      f <- reduced(Jt)
      if (is.null(f)) {
        Jt <- as_dense(K)
        f <- reduced(Jt)
      }
      if (!is.null(f)) {
        J <- Jt
        solve <- f
      }
    } else {
      fac <- tryCatch(chol(Jt), error = function(e) NULL)
      # the true Hessian may lose definiteness where Gauss-Newton cannot, and
      # the mode's movement is then read off K
      if (!is.null(fac)) {
        J <- Jt
        solve <- function(z) {
          as.numeric(backsolve(fac, forwardsolve(t(fac), as.numeric(z))))
        }
      }
    }
  }
  list(solve = solve, J = J, Dm = Dm, hold = hold)
}


#' The Log-Density's Curvature Used by the Refresh Corrections
#'
#' @description
#' Returns \eqn{\ell_{ab}} on the link scale, from whichever route was used to
#' assemble the criterion's \eqn{K}.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients.
#' @param expected Whether the information is the expected one.
#' @param approx The approximation for the expected information.
#'
#' @return A named list of numeric vectors.
#'
#' @keywords internal
refresh_hessian <- function(spec, design, coef, expected = FALSE,
                            approx = "opg") {
  th <- statmod_eta(spec, design, coef)$theta
  if (expected) {
    distributions7::distrib_expected_hessian(spec@distrib, spec@response, th,
                                             scale = "link", approx = approx, threads = spec@threads)
  } else {
    distributions7::distrib_hessian(spec@distrib, spec@response, th,
                                    scale = "link", threads = spec@threads)
  }
}


#' The Weight for Contracting a Refreshable Block
#'
#' @description
#' Computes \eqn{A_a[i,j] = w_i\sum_b c_{ab}(i)\,(X_b M_{ab}^\top)[i,j]}, the
#' matrix with which a term's own derivative is paired.
#'
#' @details
#' The same shape serves three quantities and differs only in \eqn{c_{ab}}: the
#' curvature \eqn{\ell_{ab}} for the gradient's contraction, and the third
#' derivative already contracted in a direction for each of the two mixed terms
#' of the Hessian.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param M The matrix against which the trace is taken.
#' @param params,npar,offs The block bookkeeping.
#' @param ra The term's rows in the stacked coefficient vector.
#' @param cw One length-\eqn{n} vector per distribution parameter, or
#'   `NULL` where that parameter carries no coefficient.
#'
#' @return A numeric matrix, one row per observation and one column per
#'   coefficient of the term.
#'
#' @keywords internal
refresh_amat <- function(spec, design, M, params, npar, offs, ra, cw) {
  n <- spec@n_obs
  A <- matrix(0, n, length(ra))
  for (b in seq_along(params)) {
    if (npar[b] == 0L) next
    rb <- offs[b] + seq_len(npar[b])
    Mab <- as_dense(M[ra, rb, drop = FALSE])
    A <- A + rep_len(cw[[b]], n) * as_dense(design[[params[b]]]$X %*% t(Mab))
  }
  A * spec@weights
}


#' A Refreshable Block's Derivative Along One Direction
#'
#' @description
#' Calls [modelterms7::term_block_deriv()] on one unit and checks the shape of
#' the result. A block of another shape signals an error that names the term.
#'
#' @param un One entry of [refresh_units()].
#' @param v The direction, as long as the term's coefficients.
#' @param n The number of observations.
#'
#' @return A numeric matrix, `n` by the term's coefficient count.
#'
#' @keywords internal
refresh_dblock <- function(un, v, n) {
  D <- as_dense(modelterms7::term_block_deriv(un$built, coef = un$bt, v = v))
  if (!is.matrix(D) || nrow(D) != n || ncol(D) != length(un$cols)) {
    stop(sprintf(paste("term_block_deriv() returned a %s block where %d by %d",
                       "was wanted\n  in '%s'."),
                 paste(dim(as.matrix(D)), collapse = " by "), n,
                 length(un$cols), un$term), call. = FALSE)
  }
  D
}


#' A Refreshable Block's Second Derivative Along Two Directions
#'
#' @description
#' Calls [modelterms7::term_block_deriv2()] on one unit and checks the shape of
#' the result. A block of another shape signals an error that names the term.
#'
#' @details
#' The base method on `model_term` returns zeros, which is exact for a design
#' that does not move. A sharp [modelterms7::seg()] returns zeros because its
#' second derivative vanishes away from the break-points, and a sharp
#' [modelterms7::jump()] or [modelterms7::jseg()] returns zeros because its
#' block is a working linearization with a frozen weight and not a Jacobian.
#' The function does not difference in any case.
#'
#' @param un One entry of [refresh_units()].
#' @param v,u The two directions, each of which is as long as the term's
#'   coefficients.
#' @param n The number of observations.
#'
#' @return A numeric matrix, `n` by the term's coefficient count.
#'
#' @seealso [refresh_dblock()], [trace_refresh4()]
#'
#' @keywords internal
refresh_dblock2 <- function(un, v, u, n) {
  D <- as_dense(modelterms7::term_block_deriv2(un$built, coef = un$bt, v = v,
                                               u = u))
  if (!is.matrix(D) || nrow(D) != n || ncol(D) != length(un$cols)) {
    stop(sprintf(paste("term_block_deriv2() returned a %s block where %d by %d",
                       "was wanted\n  in '%s'."),
                 paste(dim(as.matrix(D)), collapse = " by "), n,
                 length(un$cols), un$term), call. = FALSE)
  }
  D
}


#' The Weight for the Second Derivative of a Block
#'
#' @description
#' Computes [refresh_amat()] with the curvature, one matrix per unit.
#'
#' @details
#' It is the same weight that [u_refresh()] builds for the gradient, and
#' it depends on neither direction, so it is built once per Hessian rather
#' than once per pair of hyperparameters.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param M The matrix against which the trace is taken.
#' @param params,npar,offs The block bookkeeping.
#' @param units The refreshable terms, from [refresh_units()].
#' @param Hl The link-scale curvature, from [refresh_hessian()].
#'
#' @return A list of matrices, one per unit.
#'
#' @seealso [trace_refresh4()], [u_refresh()]
#'
#' @keywords internal
refresh_curv_amat <- function(spec, design, M, params, npar, offs, units, Hl) {
  lapply(units, function(un) {
    cw <- vector("list", length(params))
    for (b in seq_along(params)) {
      if (npar[b] == 0L) next
      cw[[b]] <- Hl[[hess_key(params, un$a, b)]]
    }
    refresh_amat(spec, design, M, params, npar, offs, un$ra, cw)
  })
}


#' The Third Derivative Already Contracted in One Direction
#'
#' @description
#' Computes \eqn{\sum_k \ell_{abk}(X_k v_k)_i}, the per-observation weight
#' that [contract3()] builds for the \eqn{(a,b)} block.
#'
#' @param spec A [StatmodSpec()].
#' @param d3 The third derivatives on the link scale.
#' @param params,npar The parameter names and block sizes.
#' @param a,b The block indices.
#' @param tv The predictors of the direction.
#'
#' @return A numeric vector as long as the sample.
#'
#' @keywords internal
d3_direction <- function(spec, d3, params, npar, a, b, tv) {
  n <- spec@n_obs
  keys <- names(d3)
  w <- numeric(n)
  for (k in seq_along(params)) {
    if (npar[k] == 0L) next
    w <- w + rep_len(d3[[d3_key(params, a, b, k, keys)]], n) * tv[[k]]
  }
  w
}


#' How a Moving Block Enters the Contracted Third Derivative
#'
#' @description
#' Returns the part of \eqn{T[v] = (\partial K/\partial\beta)\cdot v} that
#' [contract3()] does not compute, as a matrix over the stacked
#' coefficients.
#'
#' @details
#' Differentiating \eqn{H} in \eqn{\beta} and contracting with \eqn{v} leaves,
#' beside the third derivative that [contract3()] carries, two terms in
#' \eqn{\partial X/\partial\beta}: with \eqn{D_a = (\partial X_a/\partial\beta)v},
#' \deqn{R[(a,j),(b,k)] = -\sum_i w_i\,\ell_{ab,i}\,D_a[i,j]\,X_b[i,k],}
#' and the other is its transpose, so the correction is \eqn{R + R^\top}. The
#' derivative is read from the term through
#' [modelterms7::term_block_deriv()] and is never differenced here, for
#' the reason that [u_refresh()] records: a break-point column is a step
#' function in its break-point.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param params,npar,offs,total The block bookkeeping.
#' @param dir One entry per unit, from [refresh_direction()].
#' @param Hl The link-scale curvature, from [refresh_hessian()].
#' @param units The refreshable terms, from [refresh_units()].
#'
#' @return A square matrix, whose transpose completes the correction.
#'
#' @seealso [contract3()], [u_refresh()]
#'
#' @keywords internal
contract3_refresh <- function(spec, design, params, npar, offs, total, dir, Hl,
                              units) {
  R <- zero_information(design, total)
  if (!length(units)) return(R)
  n <- spec@n_obs
  w <- spec@weights
  for (s in seq_along(units)) {
    un <- units[[s]]
    D <- dir[[s]]$D
    for (b in seq_along(params)) {
      if (npar[b] == 0L) next
      rb <- offs[b] + seq_len(npar[b])
      lab <- rep_len(Hl[[hess_key(params, un$a, b)]], n)
      R[un$ra, rb] <- R[un$ra, rb] -
        as_dense(crossprod(D * (w * lab), design[[params[b]]]$X))
    }
  }
  R
}


#' The Derivatives of a Moving Block Along One Direction
#'
#' @description
#' Returns the two quantities that every consumer of
#' \eqn{\partial X/\partial\beta} reads for a given direction: the block's own
#' derivative along it, and the third derivative already contracted in it and
#' carried onto the term's columns.
#'
#' @details
#' Both depend on a single direction, never on a pair, so they are built
#' once per hyperparameter and combined in the Hessian's pair loop: computing
#' them inside that loop would repeat an \eqn{O(np^2)} product for every pair.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param M The matrix against which the trace is taken.
#' @param params,npar,offs The block bookkeeping.
#' @param d3 The third derivatives on the link scale.
#' @param tv The predictors of the direction.
#' @param v The direction, over the stacked coefficients.
#' @param units The refreshable terms, from [refresh_units()].
#'
#' @return A list, one entry per unit, with `D` and `A`.
#'
#' @seealso [contract3_refresh()], [trace_refresh4()]
#'
#' @keywords internal
refresh_direction <- function(spec, design, M, params, npar, offs, d3, tv, v,
                              units) {
  lapply(units, function(un) {
    cw <- vector("list", length(params))
    for (b in seq_along(params)) {
      if (npar[b] == 0L) next
      cw[[b]] <- d3_direction(spec, d3, params, npar, un$a, b, tv)
    }
    list(D = refresh_dblock(un, v[un$ra], spec@n_obs),
         A = refresh_amat(spec, design, M, params, npar, offs, un$ra, cw))
  })
}


#' How a Moving Block Enters the Twice-Contracted Fourth Derivative
#'
#' @description
#' Returns the part of \eqn{\mathrm{tr}(M\,U[v,u])} that [trace_design_form()]
#' does not compute, where \eqn{U} is the second derivative of \eqn{K} in the
#' coefficients contracted in two directions.
#'
#' @details
#' Differentiating the three contributions of \eqn{\partial K/\partial\beta}
#' once more leaves five terms. With \eqn{D_a = (\partial X_a/\partial\beta)v}
#' and \eqn{E_a = (\partial X_a/\partial\beta)u}, three of them read only the
#' first derivative,
#' \deqn{N_1[(a,j),(b,k)] = -\sum_i w_i\Big(\sum_m \ell_{abm}(X_mv_m)_i\Big)
#'     E_a[i,j]X_b[i,k],}
#' \eqn{N_2} is the same expression when the two directions are exchanged, and
#' \eqn{N_3[(a,j),(b,k)] = -\sum_i w_i\ell_{ab,i}D_a[i,j]E_b[i,k]}, which is
#' supported where both blocks move. A fourth reads the second,
#' \eqn{N_4[(a,j),(b,k)] = -\sum_i w_i\ell_{ab,i}F_a[i,j]X_b[i,k]} with
#' \eqn{F_a = (\partial^2X_a/\partial\beta^2)[v,u]} from
#' [modelterms7::term_block_deriv2()]. Those four enter as
#' \eqn{\mathrm{tr}(M(N + N^\top)) = 2\sum M\odot N}. \eqn{N_1}, \eqn{N_2} and
#' \eqn{N_4} are read off [refresh_amat()], with the third derivative for the
#' first two and the curvature for the fourth, and \eqn{N_3} is summed from
#' the two blocks' derivatives with the rows of \eqn{M} between them, so that
#' the function assembles neither \eqn{N} nor any contraction.
#'
#' The fifth is of another shape and is added separately. The
#' third-derivative contraction carries the direction's predictor
#' \eqn{(X_mv_m)_i}, which moves with \eqn{\beta} as everything else does where
#' \eqn{X_m} is a Jacobian, and differentiating it leaves
#' \deqn{-\sum_i w_i\Big(\sum_m \ell_{abm}\,\dot v_m(i)\Big)X_a[i,j]X_b[i,k],
#'   \qquad \dot v_m = (\partial X_m/\partial\beta[u])\,v_m.}
#' That is the shape that [trace_design_form()] already computes, so the
#' piece is one further call of it with \eqn{\dot v} in place of the
#' predictors. \eqn{\dot v} is the predictor's own second derivative
#' contracted in both directions and is therefore symmetric in them, which is
#' a free check on the assembly.
#'
#' @param spec A [StatmodSpec()].
#' @param M The matrix against which the trace is taken.
#' @param params,npar The parameter names and block sizes.
#' @param Hl The link-scale curvature, from [refresh_hessian()].
#' @param dv,du The two directions, from [refresh_direction()].
#' @param units The refreshable terms, from [refresh_units()].
#' @param G The leverage diagonal, from [block_leverage()].
#' @param d3 The third derivatives on the link scale.
#' @param v The first direction over the stacked coefficients.
#' @param f2 The block's second derivative in the two directions, one matrix
#'   per unit, from [refresh_dblock2()].
#' @param acurv The curvature weights, from [refresh_curv_amat()].
#'
#' @return A single number.
#'
#' @seealso [trace_design_form()], [contract3_refresh()],
#'   [refresh_dblock2()]
#'
#' @keywords internal
trace_refresh4 <- function(spec, M, params, npar, Hl, dv, du, units, G, d3,
                           v, f2, acurv) {
  if (!length(units)) return(0)
  n <- spec@n_obs
  w <- spec@weights
  out <- 0
  vdot <- lapply(seq_along(params), function(m) numeric(n))
  for (s in seq_along(units)) {
    un <- units[[s]]
    out <- out - sum(du[[s]]$D * dv[[s]]$A) - sum(dv[[s]]$D * du[[s]]$A)
    for (t2 in seq_along(units)) {
      Mab <- as_dense(M[un$ra, units[[t2]]$ra, drop = FALSE])
      lab <- rep_len(Hl[[hess_key(params, un$a, units[[t2]]$a)]], n)
      out <- out - sum(w * lab * rowSums((dv[[s]]$D %*% Mab) * du[[t2]]$D))
    }
    out <- out - sum(f2[[s]] * acurv[[s]])
    vdot[[un$a]] <- vdot[[un$a]] + as.numeric(du[[s]]$D %*% v[un$ra])
  }
  2 * out + trace_design_form(spec, G, d3, params, npar, vdot)
}


#' How a Moving Block Enters the Mode's Second Movement
#'
#' @description
#' Returns the part of \eqn{\partial^3L/\partial\beta^3[v,u]} that
#' [contract3()] and [contract3_refresh()] do not carry,
#' as a vector over the stacked coefficients.
#'
#' @details
#' \eqn{b_{ml}} solves \eqn{J b_{ml} = -[(S_l + T[b_l])b_m + S_m b_l + c_{ml}]}
#' with \eqn{T} the third derivative of the penalized objective in \eqn{\beta}.
#' Where a block moves, the objective's second derivative is not \eqn{K} but
#' \eqn{K + D}, with \eqn{D} the term that [mode_curvature()] builds, so
#' its third derivative is not \eqn{\partial K/\partial\beta} either. The term
#' \eqn{\partial D/\partial\beta} is missing, and differentiating
#' \eqn{D_{cd} = -\sum_i w_i\ell_a\,\Xi_a[i,c,d]} once and contracting gives two
#' pieces,
#' \deqn{-\sum_i w_i\Big(\sum_k \ell_{ak}(X_kv_k)_i\Big)\,
#'   \big[(\partial X_a/\partial\beta)u\big][i,c]
#'   \;-\;\sum_i w_i\,\ell_a(i)\,F_a[i,c],}
#' of which the first reads the block's first derivative and the second reads
#' its second derivative \eqn{F_a = (\partial^2X_a/\partial\beta^2)[v,u]}. Both
#' are one `crossprod` against a per-observation weight, so the function forms
#' neither the third-derivative array of the predictor nor any contraction of
#' it.
#'
#' @param spec A [StatmodSpec()].
#' @param params,npar,total The block bookkeeping.
#' @param units The refreshable terms, from [refresh_units()].
#' @param Hl The link-scale curvature, from [refresh_hessian()].
#' @param gl The link-scale score.
#' @param tv The predictors of the first direction.
#' @param du The second direction, from [refresh_direction()].
#' @param f2 The block's second derivative, one matrix per unit.
#'
#' @return A numeric vector as long as the stacked coefficients.
#'
#' @seealso [mode_curvature()], [trace_refresh4()]
#'
#' @keywords internal
refresh_mode_third <- function(spec, params, npar, units, Hl, gl, tv, du, f2,
                               total) {
  out <- numeric(total)
  if (!length(units)) return(out)
  n <- spec@n_obs
  w <- spec@weights
  for (s in seq_along(units)) {
    un <- units[[s]]
    hk <- numeric(n)
    for (k in seq_along(params)) {
      if (npar[k] == 0L) next
      hk <- hk + rep_len(Hl[[hess_key(params, un$a, k)]], n) * tv[[k]]
    }
    out[un$ra] <- out[un$ra] -
      as.numeric(crossprod(du[[s]]$D, w * hk)) -
      as.numeric(crossprod(f2[[s]], w * rep_len(gl[[un$param]], n)))
  }
  out
}
