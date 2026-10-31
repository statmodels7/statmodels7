#' @include jfp_engine.R outer_refresh.R
NULL

# The pieces of the third and fourth derivatives of log|I_tilde| when the
# predictors are not linear in the joint vector (a term whose block is
# recomputed at the coefficients, nl() or a smoothed seg(), sits in some
# equation) and when the block's own rows move (the block lies in the
# subformula of such a term). A block has rows in one equation a, so
#
#   I_tilde = sum_t w_t J_t V_t' V_t,   J_t = I_aa(eta_t),
#
# and a derivative of I_tilde along a set of directions distributes them over
# the three factors by the Leibniz rule. In whitened rows W_S = V_S U^-1 a
# derivative along a set R of directions, written D_R(f; P), is
#
#   sum over R = S1 + S2 + S3 of  d_S1 f  P(S2, S3),
#
# written out case by case in jfp_leibniz() for |R| = 0, 1, 2. The scalar f is
# J or one of its derivatives in the predictors (J_c, J_cd), and its derivative
# along S1 is the Faa di Bruno sum in the predictors, which for one or two
# directions is
#
#   d_v f  = sum_c f_c eta_cv,
#   d_vw f = sum_cd f_cd eta_cv eta_dw + sum_c f_c eta_cvw,
#
# the second term existing only where a predictor is not linear (jfp_mv_J()).
# A free index k of the joint vector enters the information through
# J_k = sum_c J_c eta_ck, so that d_S J_k = sum over S' in S of
# d_{S - S'} J_c eta_{c k S'}: eta_ck is the Jacobian column, and eta_{ckS'}
# with S' not empty is a column of the derivative of a moving term's block
# (term_block_deriv(), term_block_deriv2()). It enters the rows only for k in
# the term that holds the block, through V_{kS}.
#
# The three objects that R/jfp_outer.R reads are then
#
#   A_R    = sum_t w D_R(J; Pr(W, W)),                   p x p
#   B_R[k] = sum_t w [ sum_c sum_S' eta_{ckS'} D_{R-S'}(J_c; Pr(W, W))
#                      + D_R(J; Pr(W_k, W) + Pr(W, W_k)) ],   one row per k
#   T_R(Y)_kl = sum_t w D_{R + {k, l}}(J; Y(W, W)),
#
# with Pr(A, B) the per-observation vec(A_t' B_t) and Y(A, B) = A_t Y B_t'.
# T_R splits by where k and l land: both on the information (the cross
# product of the designs, plus the derivatives of the predictor in k or l, or
# in both), one on the information and one on the rows, or both on the rows
# (on two factors, or on one, where the row's own second derivative in k and
# l is contracted against the block's columns, jfp_mv_cv()).


#' The Movement of the Third and Fourth Derivatives of a Jeffreys Block
#'
#' @description
#' `jfp_mv_ctx()` sets up, for one block and up to two directions, the
#' context from which `jfp_mv_white()`, `jfp_mv_B()` and `jfp_mv_pair()`
#' compute the three objects of the third and fourth derivatives of
#' \eqn{\log|\tilde I|}: the whitened derivative of the block's matrix along
#' the directions, the rows \eqn{B_R} of the low-rank terms, and the trace
#' terms \eqn{T_R(Y)_{kl} = \partial_R\,\mathrm{tr}(\tilde Y\,\partial_{kl}
#' \tilde I)} with \eqn{\tilde Y = U^{-1}YU^{-\top}}.
#'
#' @details
#' With the block's rows in one equation, \eqn{\tilde I = \sum_t w_tJ_t
#' \tilde V_t^\top\tilde V_t} with \eqn{J_t} the family's expected
#' information of that equation. Each derivative distributes the directions
#' over \eqn{J} and the two factors \eqn{\tilde V} by the Leibniz rule
#' ([jfp_leibniz()]), and a derivative of \eqn{J} is the Faa di Bruno sum in
#' the predictors, with the second derivative of a predictor where a term in
#' its equation recomputes its block at the coefficients ([jfp_mv_J()]). A
#' free index enters the information through the Jacobian column and through
#' the derivative of the moving block along the directions, and it enters the
#' rows only in the term that holds the block. For a model whose predictors
#' are linear and a block whose rows do not move, every share beyond the
#' design's cross products is zero and the objects are those of the static
#' derivatives.
#'
#' @param st A state, from [jfp_state()].
#' @param bs A block's entry of the state.
#' @param ds A named list of directions (`v`, and `w` for the fourth
#'   derivative), each from [jfp_dir()].
#' @param cx A context, from `jfp_mv_ctx()`.
#' @param R A character vector naming the directions taken, a subset of
#'   `names(ds)`.
#' @param Y A symmetric \eqn{\tilde p\times\tilde p} matrix.
#'
#' @return `jfp_mv_ctx()` an environment; `jfp_mv_white()` a
#'   \eqn{\tilde p\times\tilde p} matrix; `jfp_mv_B()` a matrix over the joint
#'   vector and the \eqn{\tilde p^2} entries; `jfp_mv_pair()` a list with
#'   `sparse` (per-observation weights keyed as [statmod_information_at()]
#'   reads its `H`) and `dense` (a list of entries `idx`, `X`, each adding
#'   `X` in the rows `idx` and its transpose in the columns `idx`).
#'
#' @keywords internal
jfp_mv_ctx <- function(st, bs, ds) {
  live <- jfp_live(bs)
  if (length(live) != 1L) {
    stop("A Jeffreys block must have rows in exactly one equation.",
         call. = FALSE)
  }
  if (is.null(st$units)) {
    st$units <- refresh_units(st$spec, st$design, st$coef, st$params,
                              st$npar, st$offs)
  }
  cx <- new.env(parent = emptyenv())
  cx$st <- st
  cx$bs <- bs
  cx$a <- live
  cx$ds <- ds
  cx$units <- st$units
  cx$moving <- length(st$units) > 0L
  cx$host <- identical(bs$block$kind, "refresh")
  cx$p <- ncol(bs$Z[[live]])
  cx$m <- if (cx$host) length(bs$tcols) else 0L
  cx$Ui <- backsolve(bs$U, diag(cx$p))
  cx$memo <- list()
  p <- cx$p
  ii <- rep(seq_len(p), times = p)
  jj <- rep(seq_len(p), each = p)
  cx$pr <- function(x, y) x[, ii, drop = FALSE] * y[, jj, drop = FALSE]
  cx$WW <- jfp_mv_pfun(cx$pr, function(S) jfp_mv_W(cx, S),
                       function(S) jfp_mv_W(cx, S))
  cx
}

#' @rdname jfp_mv_ctx
#' @keywords internal
jfp_mv_white <- function(cx, R) {
  F <- jfp_leibniz(function(S) jfp_mv_J(cx, integer(0), S), R, cx$WW)
  out <- matrix(colSums(cx$st$spec@weights * F), cx$p, cx$p)
  (out + t(out)) / 2
}

#' @rdname jfp_mv_ctx
#' @keywords internal
jfp_mv_B <- function(cx, R) {
  st <- cx$st
  w <- st$spec@weights
  B <- matrix(0, st$total, cx$p * cx$p)
  for (cc in seq_along(st$params)) {
    if (st$npar[cc] == 0L) next
    jq <- function(S) jfp_mv_J(cx, cc, S)
    F <- jfp_leibniz(jq, R, cx$WW)
    rows <- st$offs[cc] + seq_len(st$npar[cc])
    if (!is.null(F)) {
      B[rows, ] <- B[rows, ] +
        as.matrix(crossprod(st$design[[st$params[cc]]]$X, w * F))
    }
    for (i in seq_along(cx$units)) {
      un <- cx$units[[i]]
      if (un$a != cc) next
      for (S in jfp_nonempty(R)) {
        F1 <- jfp_leibniz(jq, setdiff(R, S), cx$WW)
        if (is.null(F1)) next
        B[un$ra, ] <- B[un$ra, ] + crossprod(jfp_mv_etak(cx, i, S), w * F1)
      }
    }
  }
  if (cx$host) {
    jJ <- function(S) jfp_mv_J(cx, integer(0), S)
    Wf <- function(S) jfp_mv_W(cx, S)
    ti <- cx$bs$tidx
    for (k in seq_len(cx$m)) {
      Wk <- function(S) jfp_mv_Wk(cx, S)[[k]]
      F <- jfp_leibniz(jJ, R, jfp_mv_pfun(cx$pr, Wk, Wf)) +
        jfp_leibniz(jJ, R, jfp_mv_pfun(cx$pr, Wf, Wk))
      B[ti[k], ] <- B[ti[k], ] + colSums(w * F)
    }
  }
  B
}

#' @rdname jfp_mv_ctx
#' @keywords internal
jfp_mv_pair <- function(cx, Y, R) {
  st <- cx$st
  params <- st$params
  np <- length(params)
  n <- st$spec@n_obs
  w <- st$spec@weights
  total <- st$total
  e <- character(0)
  prY <- function(x, y) rowSums((x %*% Y) * y)
  Wf <- function(S) jfp_mv_W(cx, S)
  WWy <- jfp_mv_pfun(prY, Wf, Wf)
  sparse <- list()
  dense <- list()
  # both indices on the information through the Jacobian columns: the
  # weighted cross product of the designs
  for (cc in seq_len(np)) {
    for (dd in cc:np) {
      f <- jfp_leibniz(function(S) jfp_mv_J(cx, c(cc, dd), S), R, WWy)
      sparse[[hess_key(params, cc, dd)]] <- if (is.null(f)) numeric(n) else f
    }
  }
  if (cx$moving) {
    units <- cx$units
    for (i in seq_along(units)) {
      un <- units[[i]]
      cc <- un$a
      # k through the derivative of a moving block along the directions S, l
      # through the Jacobian column; the transpose carries the other order
      for (S in jfp_nonempty(R)) {
        Ek <- jfp_mv_etak(cx, i, S)
        X <- matrix(0, length(un$ra), total)
        for (dd in seq_len(np)) {
          if (st$npar[dd] == 0L) next
          f <- jfp_leibniz(function(S0) jfp_mv_J(cx, c(cc, dd), S0),
                           setdiff(R, S), WWy)
          if (is.null(f)) next
          rows <- st$offs[dd] + seq_len(st$npar[dd])
          X[, rows] <- X[, rows] + as.matrix(
            crossprod(Ek, (w * f) * st$design[[params[dd]]]$X))
        }
        dense[[length(dense) + 1L]] <- list(idx = un$ra, X = X)
      }
      # k and l both through the derivatives of moving blocks, one direction
      # each; the pair (w, v) is the transpose of (v, w)
      if (length(R) == 2L) {
        Ev <- jfp_mv_etak(cx, i, R[1L])
        X <- matrix(0, length(un$ra), total)
        for (j in seq_along(units)) {
          u2 <- units[[j]]
          f <- jfp_leibniz(function(S0) jfp_mv_J(cx, c(cc, u2$a), S0), e, WWy)
          X[, u2$ra] <- X[, u2$ra] +
            crossprod(Ev, (w * f) * jfp_mv_etak(cx, j, R[2L]))
        }
        dense[[length(dense) + 1L]] <- list(idx = un$ra, X = X)
      }
      # the second derivative of the predictor in k and l, moved along S
      for (S in c(list(e), jfp_nonempty(R))) {
        f <- jfp_leibniz(function(S0) jfp_mv_J(cx, cc, S0), setdiff(R, S),
                         WWy)
        if (is.null(f)) next
        X <- matrix(0, length(un$ra), total)
        X[, un$ra] <- jfp_mv_etakk(cx, i, S, w * f) / 2
        dense[[length(dense) + 1L]] <- list(idx = un$ra, X = X)
      }
    }
  }
  if (cx$host) {
    bs <- cx$bs
    m <- cx$m
    ti <- bs$tidx
    Wkf <- function(k) function(S) jfp_mv_Wk(cx, S)[[k]]
    # k on the information, l on the rows (a factor 2 for the two factors)
    for (cc in seq_len(np)) {
      if (st$npar[cc] == 0L) next
      jq <- function(S0) jfp_mv_J(cx, cc, S0)
      for (S in c(list(e), jfp_nonempty(R))) {
        Rl <- matrix(0, n, m)
        for (l in seq_len(m)) {
          f <- jfp_leibniz(jq, setdiff(R, S), jfp_mv_pfun(prY, Wkf(l), Wf))
          if (!is.null(f)) Rl[, l] <- 2 * f
        }
        if (!length(S)) {
          X <- matrix(0, m, total)
          rows <- st$offs[cc] + seq_len(st$npar[cc])
          X[, rows] <- t(as.matrix(crossprod(st$design[[params[cc]]]$X,
                                             w * Rl)))
          dense[[length(dense) + 1L]] <- list(idx = ti, X = X)
        } else {
          for (i in seq_along(cx$units)) {
            un <- cx$units[[i]]
            if (un$a != cc) next
            X <- matrix(0, m, total)
            X[, un$ra] <- t(crossprod(jfp_mv_etak(cx, i, S), w * Rl))
            dense[[length(dense) + 1L]] <- list(idx = ti, X = X)
          }
        }
      }
    }
    # k and l both on the rows: on the two factors, and on one factor through
    # the rows' second derivative in k and l
    jJ <- function(S0) jfp_mv_J(cx, integer(0), S0)
    G <- matrix(0, m, m)
    for (k in seq_len(m)) {
      for (l in k:m) {
        f <- jfp_leibniz(jJ, R, jfp_mv_pfun(prY, Wkf(k), Wkf(l)))
        if (!is.null(f)) {
          G[k, l] <- 2 * sum(w * f)
          G[l, k] <- G[k, l]
        }
      }
    }
    for (S in c(list(e), jfp_nonempty(R))) {
      H1 <- jfp_leibniz(jJ, setdiff(R, S), function(S2, S3)
        if (length(S2)) NULL else jfp_mv_W(cx, S3))
      if (is.null(H1)) next
      Q <- (w * H1) %*% Y %*% t(cx$Ui) %*% t(bs$block$C)
      G <- G + 2 * jfp_mv_cv(cx, S, Q)
    }
    X <- matrix(0, m, total)
    X[, ti] <- (G + t(G)) / 4
    dense[[length(dense) + 1L]] <- list(idx = ti, X = X)
  }
  list(sparse = sparse, dense = dense)
}


#' The Leibniz Rule Over Three Factors
#'
#' @description
#' Returns \eqn{\sum_{R = S_1\sqcup S_2\sqcup S_3}\partial_{S_1}f\;
#' P(S_2, S_3)} for a set \eqn{R} of at most two directions, written out case
#' by case: one term for no direction, three for one, nine for two.
#'
#' @param jq A function of a set of directions returning the derivative of the
#'   scalar factor along them, a numeric vector over the observations.
#' @param R A character vector of at most two direction names.
#' @param P A function of two sets of directions returning the product of the
#'   other two factors differentiated along them, or `NULL` where it is zero.
#'
#' @return The sum, or `NULL` where every term is zero.
#'
#' @keywords internal
jfp_leibniz <- function(jq, R, P) {
  acc <- NULL
  add <- function(S1, S2, S3) {
    x <- P(S2, S3)
    if (is.null(x)) return(invisible(NULL))
    j <- jq(S1)
    if (is.null(j)) return(invisible(NULL))
    acc <<- if (is.null(acc)) j * x else acc + j * x
    invisible(NULL)
  }
  e <- character(0)
  if (length(R) == 0L) {
    add(e, e, e)
  } else if (length(R) == 1L) {
    add(R, e, e)
    add(e, R, e)
    add(e, e, R)
  } else {
    v <- R[1L]
    w <- R[2L]
    add(R, e, e)
    add(v, w, e)
    add(v, e, w)
    add(w, v, e)
    add(w, e, v)
    add(e, R, e)
    add(e, v, w)
    add(e, w, v)
    add(e, e, R)
  }
  acc
}


#' The Building Blocks of the Moving Shares
#'
#' @description
#' `jfp_dir()` carries a direction over the joint vector with the movement of
#' every predictor along it. `jfp_mv_J()` returns the derivative of the
#' expected information of the block's equation, with the predictor indices
#' `fixed`, along the directions `R`. `jfp_mv_e2()` returns the second
#' derivative of every predictor along the two directions. `jfp_mv_W()` and
#' `jfp_mv_Wk()` return the whitened derivatives of the block's rows along
#' `R`, the second one also in each coefficient of the term that holds the
#' block. `jfp_mv_etak()` returns the derivative of a moving term's block
#' along `R`; `jfp_mv_etakk()` its derivative in two of the term's
#' coefficients along `R`, contracted against per-observation weights;
#' `jfp_mv_cv()` the derivative of the block's rows in two coefficients of
#' its term along `R`, contracted against weights on the block's columns.
#' `jfp_mv_pfun()` wraps a product of two factors with a memory, and
#' `jfp_nonempty()` lists the non-empty subsets of at most two directions.
#'
#' @param st A state, from [jfp_state()].
#' @param v A direction over the joint vector.
#' @param cx A context, from [jfp_mv_ctx()].
#' @param fixed Equation indices of the derivative in the predictors.
#' @param R,S A character vector naming directions.
#' @param i The position of a moving term among the context's units.
#' @param g Per-observation weights.
#' @param Q Weights over the observations and the block's columns.
#' @param pr A product of two matrices.
#' @param A,B Functions of a set of directions returning a factor.
#'
#' @return As described; `NULL` where a factor is zero.
#'
#' @keywords internal
jfp_dir <- function(st, v) {
  list(u = as.numeric(v), e1 = jfp_eta_dir(st, v))
}

#' @rdname jfp_dir
#' @keywords internal
jfp_mv_J <- function(cx, fixed, R) {
  key <- paste0("J", paste(fixed, collapse = "."), "|",
                paste(R, collapse = ""))
  if (is.null(cx$memo[[key]])) {
    st <- cx$st
    a <- cx$a
    out <- switch(length(R) + 1L,
      jfp_dI(st, a, a, fixed, list()),
      jfp_dI(st, a, a, fixed, list(cx$ds[[R]]$e1)),
      {
        o <- jfp_dI(st, a, a, fixed, list(cx$ds[[R[1L]]]$e1,
                                          cx$ds[[R[2L]]]$e1))
        if (cx$moving) o <- o + jfp_dI(st, a, a, fixed, list(jfp_mv_e2(cx)))
        o
      })
    cx$memo[[key]] <- out
  }
  cx$memo[[key]]
}

#' @rdname jfp_dir
#' @keywords internal
jfp_mv_e2 <- function(cx) {
  if (is.null(cx$memo$e2)) {
    st <- cx$st
    n <- st$spec@n_obs
    out <- lapply(st$params, function(p) numeric(n))
    w <- cx$ds[[2L]]$u
    for (i in seq_along(cx$units)) {
      un <- cx$units[[i]]
      out[[un$a]] <- out[[un$a]] +
        as.numeric(jfp_mv_etak(cx, i, names(cx$ds)[1L]) %*% w[un$ra])
    }
    cx$memo$e2 <- out
  }
  cx$memo$e2
}

#' @rdname jfp_dir
#' @keywords internal
jfp_mv_W <- function(cx, S) {
  bs <- cx$bs
  if (!length(S)) return(bs$Z[[cx$a]])
  if (!cx$host) return(NULL)
  key <- paste0("W|", paste(S, collapse = ""))
  if (is.null(cx$memo[[key]])) {
    dv <- lapply(S, function(d) cx$ds[[d]]$u[bs$tidx])
    D <- if (length(S) == 1L) {
      modelterms7::term_block_deriv(bs$tm, coef = bs$bt, v = dv[[1L]])
    } else {
      modelterms7::term_block_deriv2(bs$tm, coef = bs$bt, v = dv[[1L]],
                                     u = dv[[2L]])
    }
    cx$memo[[key]] <- as.matrix(D)[, bs$pos, drop = FALSE] %*%
      bs$block$C %*% cx$Ui
  }
  cx$memo[[key]]
}

#' @rdname jfp_dir
#' @keywords internal
jfp_mv_Wk <- function(cx, S) {
  key <- paste0("Wk|", paste(S, collapse = ""))
  if (is.null(cx$memo[[key]])) {
    bs <- cx$bs
    m <- cx$m
    dv <- lapply(S, function(d) cx$ds[[d]]$u[bs$tidx])
    cx$memo[[key]] <- lapply(seq_len(m), function(k) {
      if (!length(S)) return(jfp_vk(bs, k) %*% cx$Ui)
      ek <- numeric(m)
      ek[k] <- 1
      D <- if (length(S) == 1L) {
        modelterms7::term_block_deriv2(bs$tm, coef = bs$bt, v = dv[[1L]],
                                       u = ek)
      } else {
        modelterms7::term_block_deriv3(bs$tm, coef = bs$bt, v = dv[[1L]],
                                       u = dv[[2L]], s = ek)
      }
      as.matrix(D)[, bs$pos, drop = FALSE] %*% bs$block$C %*% cx$Ui
    })
  }
  cx$memo[[key]]
}

#' @rdname jfp_dir
#' @keywords internal
jfp_mv_etak <- function(cx, i, S) {
  key <- paste0("E", i, "|", paste(S, collapse = ""))
  if (is.null(cx$memo[[key]])) {
    un <- cx$units[[i]]
    dv <- lapply(S, function(d) cx$ds[[d]]$u[un$ra])
    D <- if (length(S) == 1L) {
      modelterms7::term_block_deriv(un$built, coef = un$bt, v = dv[[1L]])
    } else {
      modelterms7::term_block_deriv2(un$built, coef = un$bt, v = dv[[1L]],
                                     u = dv[[2L]])
    }
    cx$memo[[key]] <- as.matrix(D)
  }
  cx$memo[[key]]
}

#' @rdname jfp_dir
#' @keywords internal
jfp_mv_etakk <- function(cx, i, S, g) {
  un <- cx$units[[i]]
  m <- length(un$ra)
  dv <- lapply(S, function(d) cx$ds[[d]]$u[un$ra])
  G <- matrix(0, m, m)
  for (k in seq_len(m)) {
    ek <- numeric(m)
    ek[k] <- 1
    D <- switch(length(S) + 1L,
      modelterms7::term_block_deriv(un$built, coef = un$bt, v = ek),
      modelterms7::term_block_deriv2(un$built, coef = un$bt, v = dv[[1L]],
                                     u = ek),
      modelterms7::term_block_deriv3(un$built, coef = un$bt, v = dv[[1L]],
                                     u = dv[[2L]], s = ek))
    G[k, ] <- as.numeric(crossprod(as.matrix(D), g))
  }
  (G + t(G)) / 2
}

#' @rdname jfp_dir
#' @keywords internal
jfp_mv_cv <- function(cx, S, Q) {
  bs <- cx$bs
  m <- cx$m
  dv <- lapply(S, function(d) cx$ds[[d]]$u[bs$tidx])
  G <- matrix(0, m, m)
  for (j in seq_along(bs$pos)) {
    ej <- numeric(m)
    ej[bs$pos[j]] <- 1
    for (k in seq_len(m)) {
      ek <- numeric(m)
      ek[k] <- 1
      D <- switch(length(S) + 1L,
        modelterms7::term_block_deriv2(bs$tm, coef = bs$bt, v = ej, u = ek),
        modelterms7::term_block_deriv3(bs$tm, coef = bs$bt, v = ej, u = ek,
                                       s = dv[[1L]]),
        modelterms7::term_block_deriv4(bs$tm, coef = bs$bt, v = ej, u = ek,
                                       s = dv[[1L]], r = dv[[2L]]))
      G[k, ] <- G[k, ] + as.numeric(crossprod(as.matrix(D), Q[, j]))
    }
  }
  G
}

#' @rdname jfp_dir
#' @keywords internal
jfp_mv_pfun <- function(pr, A, B) {
  memo <- new.env(parent = emptyenv())
  function(S2, S3) {
    key <- paste0(paste(S2, collapse = ""), "|", paste(S3, collapse = ""))
    if (!exists(key, envir = memo, inherits = FALSE)) {
      a <- A(S2)
      b <- if (is.null(a)) NULL else B(S3)
      assign(key, if (is.null(b)) NULL else pr(a, b), envir = memo)
    }
    get(key, envir = memo, inherits = FALSE)
  }
}

#' @rdname jfp_dir
#' @keywords internal
jfp_nonempty <- function(R) {
  switch(length(R) + 1L, list(), list(R), list(R[1L], R[2L], R))
}
