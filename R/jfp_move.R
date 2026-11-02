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
#' In a model with a filter the Jacobian columns are the joint rows of every
#' equation ([jfp_jac()]), and each equation that the filter drives is a
#' moving unit whose derivatives along the directions come from the filter's
#' kernels ([jfp_mv_rows()]). A block in the subformula of a filter's
#' parameter has rows in every equation that the filter drives, so
#' \eqn{\tilde I = \sum_t w_t\sum_{q,r}I_{qr,t}\tilde V_{q,t}^\top
#' \tilde V_{r,t}} and the three objects are summed over the ordered pairs
#' \eqn{(q, r)} ([jfp_mv_setpair()]); the rows' derivatives are read along
#' the block's whitened directions \eqn{d_i}, one kernel call per direction.
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
  filt <- identical(bs$block$kind, "filter")
  if (length(live) != 1L && !filt) {
    stop("A Jeffreys block must have rows in exactly one equation.",
         call. = FALSE)
  }
  cx <- new.env(parent = emptyenv())
  cx$st <- st
  cx$bs <- bs
  cx$live <- live
  cx$a <- live[1L]
  cx$ds <- ds
  if (!is.null(st$joint)) {
    # a filter bends the predictors of the equations it drives in every
    # coordinate of the joint vector: one unit per driven equation, read
    # from the filter's kernels
    if (is.null(st$units)) {
      st$units <- refresh_units(st$spec, st$design, st$coef, st$params,
                                st$npar, st$offs)
    }
    driven <- match(st$joint$driven, st$params)
    if (length(st$units) &&
        (filt || !all(vapply(st$units, function(u) u$a %in% driven,
                             logical(1))))) {
      stop("The third derivative of a Jeffreys penalty beside a filter and ",
           "a term whose block moves is implemented only for a block that ",
           "moves in an equation
  that the filter drives, with the ",
           "Jeffreys block outside the filter.", call. = FALSE)
    }
    # a block that moves in a driven equation enters through the filter
    # (jfp_static_rows(), jfp_static_curv())
    cx$static <- length(st$units) > 0L
    cx$units <- lapply(driven, function(a)
      list(a = a, ra = seq_len(st$total), filter = TRUE))
    cx$keep <- st$joint$keep
    cx$mf <- st$joint$mf
    cx$zero <- stats::setNames(lapply(st$params, function(p)
      numeric(st$spec@n_obs)), st$params)
  } else {
    if (is.null(st$units)) {
      st$units <- refresh_units(st$spec, st$design, st$coef, st$params,
                                st$npar, st$offs)
    }
    cx$units <- st$units
  }
  cx$moving <- length(cx$units) > 0L
  cx$host <- bs$block$kind %in% c("refresh", "filter")
  cx$filt <- filt
  cx$p <- ncol(bs$Z[[cx$a]])
  cx$m <- if (filt) st$total else if (cx$host) length(bs$tcols) else 0L
  cx$ti <- if (filt) seq_len(st$total) else bs$tidx
  cx$Ui <- backsolve(bs$U, diag(cx$p))
  cx$memo <- list()
  p <- cx$p
  ii <- rep(seq_len(p), times = p)
  jj <- rep(seq_len(p), each = p)
  cx$pr <- function(x, y) x[, ii, drop = FALSE] * y[, jj, drop = FALSE]
  # the ordered pairs (q, r) of the equations in which the block has rows:
  # I_tilde = sum_t w sum_qr I_qr V_q' V_r
  cx$pairs <- as.matrix(expand.grid(q = live, r = live))
  cx$WWl <- list()
  cx$q <- cx$a
  cx$r <- cx$a
  cx
}

#' @rdname jfp_mv_ctx
#' @keywords internal
jfp_mv_white <- function(cx, R) {
  out <- matrix(0, cx$p, cx$p)
  for (k in seq_len(nrow(cx$pairs))) {
    jfp_mv_setpair(cx, k)
    F <- jfp_leibniz(function(S) jfp_mv_J(cx, integer(0), S), R, cx$WW)
    if (!is.null(F)) {
      out <- out + matrix(colSums(cx$st$spec@weights * F), cx$p, cx$p)
    }
  }
  (out + t(out)) / 2
}

#' @rdname jfp_mv_ctx
#' @keywords internal
jfp_mv_B <- function(cx, R) {
  st <- cx$st
  w <- st$spec@weights
  B <- matrix(0, st$total, cx$p * cx$p)
  for (pk in seq_len(nrow(cx$pairs))) {
    jfp_mv_setpair(cx, pk)
    q <- cx$q
    r <- cx$r
    for (cc in seq_along(st$params)) {
      jc <- jfp_jac(st, cc)
      if (is.null(jc)) next
      jq <- function(S) jfp_mv_J(cx, cc, S)
      F <- jfp_leibniz(jq, R, cx$WW)
      if (!is.null(F)) {
        B[jc$idx, ] <- B[jc$idx, ] + as.matrix(crossprod(jc$X, w * F))
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
      Wq <- function(S) jfp_mv_W(cx, S, q)
      Wr <- function(S) jfp_mv_W(cx, S, r)
      ti <- cx$ti
      for (k in seq_len(cx$m)) {
        Wkq <- function(S) jfp_mv_Wk(cx, S, q)[[k]]
        Wkr <- function(S) jfp_mv_Wk(cx, S, r)[[k]]
        F <- jfp_leibniz(jJ, R, jfp_mv_pfun(cx$pr, Wkq, Wr)) +
          jfp_leibniz(jJ, R, jfp_mv_pfun(cx$pr, Wq, Wkr))
        B[ti[k], ] <- B[ti[k], ] + colSums(w * F)
      }
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
  sparse <- list()
  dense <- list()
  addw <- function(key, f) {
    sparse[[key]] <<- if (is.null(sparse[[key]])) f else sparse[[key]] + f
  }
  Sall <- c(list(e), jfp_nonempty(R))
  # the weights against which a row's own second derivative in k and l is
  # contracted, per equation of that row and per subset S of the directions
  Qacc <- list()
  G <- if (cx$host) matrix(0, cx$m, cx$m) else NULL
  for (pk in seq_len(nrow(cx$pairs))) {
    jfp_mv_setpair(cx, pk)
    q <- cx$q
    r <- cx$r
    Wq <- function(S) jfp_mv_W(cx, S, q)
    Wr <- function(S) jfp_mv_W(cx, S, r)
    WWy <- jfp_mv_pfun(prY, Wq, Wr)
    # both indices on the information through the Jacobian columns: the
    # weighted cross product of the rows
    for (cc in seq_len(np)) {
      for (dd in cc:np) {
        f <- jfp_leibniz(function(S) jfp_mv_J(cx, c(cc, dd), S), R, WWy)
        addw(hess_key(params, cc, dd), if (is.null(f)) numeric(n) else f)
      }
    }
    if (cx$moving) {
      units <- cx$units
      for (i in seq_along(units)) {
        un <- units[[i]]
        cc <- un$a
        # k through the derivative of a moving block along the directions
        # S, l through the Jacobian column; the transpose carries the other
        # order
        for (S in jfp_nonempty(R)) {
          Ek <- jfp_mv_etak(cx, i, S)
          X <- matrix(0, length(un$ra), total)
          for (dd in seq_len(np)) {
            jd <- jfp_jac(st, dd)
            if (is.null(jd)) next
            f <- jfp_leibniz(function(S0) jfp_mv_J(cx, c(cc, dd), S0),
                             setdiff(R, S), WWy)
            if (is.null(f)) next
            X[, jd$idx] <- X[, jd$idx] + as.matrix(
              crossprod(Ek, (w * f) * jd$X))
          }
          dense[[length(dense) + 1L]] <- list(idx = un$ra, X = X)
        }
        # k and l both through the derivatives of moving blocks, one
        # direction each; the pair (w, v) is the transpose of (v, w)
        if (length(R) == 2L) {
          Ev <- jfp_mv_etak(cx, i, R[1L])
          X <- matrix(0, length(un$ra), total)
          for (j in seq_along(units)) {
            u2 <- units[[j]]
            f <- jfp_leibniz(function(S0) jfp_mv_J(cx, c(cc, u2$a), S0), e,
                             WWy)
            X[, u2$ra] <- X[, u2$ra] +
              crossprod(Ev, (w * f) * jfp_mv_etak(cx, j, R[2L]))
          }
          dense[[length(dense) + 1L]] <- list(idx = un$ra, X = X)
        }
        # the second derivative of the predictor in k and l, moved along S
        for (S in Sall) {
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
      m <- cx$m
      ti <- cx$ti
      Wkq <- function(k) function(S) jfp_mv_Wk(cx, S, q)[[k]]
      Wkr <- function(k) function(S) jfp_mv_Wk(cx, S, r)[[k]]
      # k on the information, l on the rows of the first factor (a factor 2
      # for the two factors, the sum running over the ordered pairs)
      for (cc in seq_len(np)) {
        jc <- jfp_jac(st, cc)
        if (is.null(jc)) next
        jq <- function(S0) jfp_mv_J(cx, cc, S0)
        for (S in Sall) {
          Rl <- matrix(0, n, m)
          for (l in seq_len(m)) {
            f <- jfp_leibniz(jq, setdiff(R, S), jfp_mv_pfun(prY, Wkq(l), Wr))
            if (!is.null(f)) Rl[, l] <- 2 * f
          }
          if (!length(S)) {
            X <- matrix(0, m, total)
            X[, jc$idx] <- t(as.matrix(crossprod(jc$X, w * Rl)))
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
      # k and l both on the rows, on the two factors
      jJ <- function(S0) jfp_mv_J(cx, integer(0), S0)
      for (k in seq_len(m)) {
        for (l in k:m) {
          f <- jfp_leibniz(jJ, R, jfp_mv_pfun(prY, Wkq(k), Wkr(l)))
          if (!is.null(f)) {
            G[k, l] <- G[k, l] + 2 * sum(w * f)
            if (l != k) G[l, k] <- G[l, k] + 2 * sum(w * f)
          }
        }
      }
      # and on the first factor through its second derivative in k and l
      for (S in Sall) {
        H1 <- jfp_leibniz(jJ, setdiff(R, S), function(S2, S3)
          if (length(S2)) NULL else Wr(S3))
        if (is.null(H1)) next
        key <- paste(q, paste(S, collapse = ""))
        Q <- (w * H1) %*% Y
        Qacc[[key]] <- if (is.null(Qacc[[key]])) list(q = q, S = S, Q = Q)
                       else list(q = q, S = S, Q = Qacc[[key]]$Q + Q)
      }
    }
  }
  if (cx$host) {
    for (z in Qacc) G <- G + 2 * jfp_mv_cv(cx, z$S, z$Q, z$q)
    X <- matrix(0, cx$m, total)
    X[, cx$ti] <- (G + t(G)) / 4
    dense[[length(dense) + 1L]] <- list(idx = cx$ti, X = X)
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
#' In a model with a filter, `jfp_mv_setpair()` selects the ordered pair of
#' equations whose information and rows the other functions read,
#' `jfp_mv_full()` places directions over the filter's full width,
#' `jfp_mv_rows()` returns the predictors' derivative along the directions
#' (and along the block's whitened direction `i`, when given), one row per
#' observation and per equation, and `jfp_jac()` returns an equation's
#' Jacobian columns with their positions in the joint vector; the weights
#' `Q` of `jfp_mv_cv()` are then on the whitened directions.
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
#' @param eq The equation of a factor of the block's rows.
#' @param k The position of an ordered pair of equations in the context.
#' @param cc An equation.
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
  key <- paste0("J", cx$q, ":", cx$r, "|", paste(fixed, collapse = "."), "|",
                paste(R, collapse = ""))
  if (is.null(cx$memo[[key]])) {
    st <- cx$st
    q <- cx$q
    r <- cx$r
    out <- switch(length(R) + 1L,
      jfp_dI(st, q, r, fixed, list()),
      jfp_dI(st, q, r, fixed, list(cx$ds[[R]]$e1)),
      {
        o <- jfp_dI(st, q, r, fixed, list(cx$ds[[R[1L]]]$e1,
                                          cx$ds[[R[2L]]]$e1))
        if (cx$moving) o <- o + jfp_dI(st, q, r, fixed, list(jfp_mv_e2(cx)))
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
jfp_mv_W <- function(cx, S, eq = cx$a) {
  bs <- cx$bs
  if (!length(S)) return(bs$Z[[eq]])
  if (!cx$host) return(NULL)
  key <- paste0("W", eq, "|", paste(S, collapse = ""))
  if (is.null(cx$memo[[key]])) {
    cx$memo[[key]] <- if (cx$filt) {
      # the rows' derivative along S: the predictors' derivative of order
      # 1 + |S| along S and the whitened directions d_i
      jfp_mv_rows(cx, NULL, S)[[eq]] %*% bs$D
    } else {
      dv <- lapply(S, function(d) cx$ds[[d]]$u[bs$tidx])
      D <- if (length(S) == 1L) {
        modelterms7::term_block_deriv(bs$tm, coef = bs$bt, v = dv[[1L]])
      } else {
        modelterms7::term_block_deriv2(bs$tm, coef = bs$bt, v = dv[[1L]],
                                       u = dv[[2L]])
      }
      as.matrix(D)[, bs$pos, drop = FALSE] %*% bs$block$C %*% cx$Ui
    }
  }
  cx$memo[[key]]
}

#' @rdname jfp_dir
#' @keywords internal
jfp_mv_Wk <- function(cx, S, eq = cx$a) {
  key <- paste0("Wk", eq, "|", paste(S, collapse = ""))
  if (is.null(cx$memo[[key]])) {
    bs <- cx$bs
    m <- cx$m
    if (cx$filt) {
      # column i of the k-th matrix: the derivative of order 2 + |S| of the
      # predictor along d_i, S and the coordinate k
      keep <- cx$keep
      n <- cx$st$spec@n_obs
      Rs <- lapply(seq_len(cx$p), function(i)
        jfp_mv_rows(cx, i, S)[[eq]][, keep, drop = FALSE])
      cx$memo[[key]] <- lapply(seq_len(m), function(k)
        matrix(vapply(Rs, function(Ri) Ri[, k], numeric(n)), n, cx$p))
    } else {
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
  }
  cx$memo[[key]]
}

#' @rdname jfp_dir
#' @keywords internal
jfp_mv_etak <- function(cx, i, S) {
  un <- cx$units[[i]]
  if (isTRUE(un$filter)) {
    out <- jfp_mv_rows(cx, NULL, S)[[un$a]]
    if (isTRUE(cx$static)) {
      key <- paste0("SR|", paste(S, collapse = ""))
      if (is.null(cx$memo[[key]])) {
        cx$memo[[key]] <- jfp_static_rows(cx$st, jfp_mv_full(cx, S))
      }
      out <- out + cx$memo[[key]][[un$a]]
    }
    return(out[, cx$keep, drop = FALSE])
  }
  key <- paste0("E", i, "|", paste(S, collapse = ""))
  if (is.null(cx$memo[[key]])) {
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
  if (isTRUE(un$filter)) {
    gg <- cx$zero
    gg[[un$a]] <- g
    keep <- cx$keep
    G <- if (!length(S)) {
      cx$st$joint$curv_full(gg, scale = FALSE)
    } else {
      cx$st$joint$higher(gg, jfp_mv_full(cx, S), scale = FALSE)$curv
    }
    if (isTRUE(cx$static)) {
      G <- G + jfp_static_curv(cx$st, gg, jfp_mv_full(cx, S))
    }
    G <- G[keep, keep, drop = FALSE]
    return((G + t(G)) / 2)
  }
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
jfp_mv_cv <- function(cx, S, Q, eq = cx$a) {
  bs <- cx$bs
  m <- cx$m
  if (cx$filt) {
    # sum_t sum_i Q_ti d^{2+|S|} eta_eq [d_i, S, k, l], one kernel call per
    # whitened direction
    keep <- cx$keep
    G <- matrix(0, m, m)
    for (i in seq_len(cx$p)) {
      gg <- cx$zero
      gg[[eq]] <- Q[, i]
      dirs <- c(list(bs$D[, i]), jfp_mv_full(cx, S))
      G <- G + cx$st$joint$higher(gg, dirs,
                                  scale = FALSE)$curv[keep, keep, drop = FALSE]
    }
    return(G)
  }
  Q <- Q %*% t(cx$Ui) %*% t(bs$block$C)
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


#' @rdname jfp_dir
#' @keywords internal
jfp_mv_setpair <- function(cx, k) {
  cx$q <- cx$pairs[k, 1L]
  cx$r <- cx$pairs[k, 2L]
  key <- paste(cx$q, cx$r)
  if (is.null(cx$WWl[[key]])) {
    q <- cx$q
    r <- cx$r
    cx$WWl[[key]] <- jfp_mv_pfun(cx$pr, function(S) jfp_mv_W(cx, S, q),
                                 function(S) jfp_mv_W(cx, S, r))
  }
  cx$WW <- cx$WWl[[key]]
  invisible(cx)
}


#' @rdname jfp_dir
#' @keywords internal
jfp_mv_full <- function(cx, S) {
  lapply(S, function(d) {
    x <- numeric(cx$mf)
    x[cx$keep] <- cx$ds[[d]]$u
    x
  })
}


#' @rdname jfp_dir
#' @keywords internal
jfp_mv_rows <- function(cx, i, S) {
  key <- paste0("R", if (is.null(i)) "" else i, "|", paste(S, collapse = ""))
  if (is.null(cx$memo[[key]])) {
    dirs <- c(if (!is.null(i)) list(cx$bs$D[, i]), jfp_mv_full(cx, S))
    cx$memo[[key]] <- cx$st$joint$higher(cx$zero, dirs)$rows
  }
  cx$memo[[key]]
}


#' @rdname jfp_dir
#' @keywords internal
jfp_jac <- function(st, cc) {
  if (!is.null(st$joint)) {
    return(list(idx = seq_len(st$total), X = st$joint$V[[cc]]))
  }
  if (st$npar[cc] == 0L) return(NULL)
  list(idx = st$offs[cc] + seq_len(st$npar[cc]),
       X = st$design[[st$params[cc]]]$X)
}


#' A Recomputed Block Beside a Filter in the Jeffreys Derivatives
#'
#' @description
#' Returns the parts of the derivatives of the predictors that come from a
#' term whose block is recomputed at the coefficients ([modelterms7::nl()],
#' a smoothed [modelterms7::seg()]) in a model with a filter:
#' `jfp_static_rows()` the part of
#' \eqn{\partial^{1+|S|}\eta/\partial u\,\partial S} for one or two
#' directions, one matrix per equation over every coordinate, and
#' `jfp_static_curv()` the part of \eqn{\sum_t\sum_c g_{c,t}\,
#' \partial^{2+|S|}\eta_{c,t}/\partial u^2\partial S} for none, one or two
#' directions. `jfp_movers()` returns the terms, and `jfp_static_aug()` the
#' predictors' derivatives with further unknowns whose static rows are
#' given.
#'
#' @details
#' A filter writes \eqn{\eta = F(s, \psi)} from the static predictors
#' \eqn{s(u)} of every equation and its own parameters. The kernels of the
#' filter ([jfp_joint_rows()]) differentiate \eqn{F} with the static rows
#' \eqn{R = \partial s/\partial u} held at their value, which is the whole
#' derivative only where \eqn{s} is linear in \eqn{u}. Where a block moves,
#' the derivatives of \eqn{s} of order two and more are
#' \eqn{s''[a, b] = \partial R[a]\,b} and so on
#' ([modelterms7::term_block_deriv()] and its siblings), and the derivative
#' of \eqn{\eta} along a set of indices is the sum over the partitions of the
#' set into blocks, a block of one index entering \eqn{F} as a full
#' direction and a block of \eqn{j \ge 2} indices as the static perturbation
#' \eqn{s^{(j)}} along them. The kernels return every such term once the
#' perturbations are placed as further unknowns: a column with the static
#' rows \eqn{x} and no part in \eqn{\psi} is a direction along which \eqn{F}
#' moves by \eqn{x} alone, so its entry in the kernels' result is the
#' derivative of \eqn{F} with \eqn{x} in one argument. A free index that
#' lies in a block of two or more indices is carried by one such column per
#' value of the index.
#'
#' The partitions written out are those of \eqn{\{k, v\}} and
#' \eqn{\{k, v, w\}} for the rows, and of \eqn{\{k, l\}},
#' \eqn{\{k, l, v\}} and \eqn{\{k, l, v, w\}} for the contraction, without
#' the partition into singletons, which is the kernels' own result. In an
#' equation that no filter drives \eqn{F} is the identity in that
#' equation's static predictor, so the unknown's column enters there only
#' through the contraction of the Jacobian.
#'
#' @param st A state, from [jfp_state()], of a model with a filter.
#' @param dirs A list of directions over every coordinate (the full width of
#'   the filter's rows).
#' @param g A named list of per-observation weights, already multiplied by
#'   the observation weights.
#' @param cols A list of column blocks, each a list with `name` and `rows`
#'   (one matrix or `NULL` per distribution parameter).
#'
#' @return `jfp_static_rows()` one matrix per distribution parameter over
#'   every coordinate; `jfp_static_curv()` a square matrix over every
#'   coordinate; `jfp_movers()` the units of [refresh_units()];
#'   `jfp_static_aug()` what [jfp_joint_rows()] returns, with `pos` (the
#'   columns of each block), `oi` (the columns of the original coordinates),
#'   `pad()` (a direction over the original coordinates placed on the wider
#'   width) and `unit()` (the direction of one added column).
#'
#' @keywords internal
jfp_movers <- function(st) {
  if (is.null(st$units)) {
    st$units <- refresh_units(st$spec, st$design, st$coef, st$params,
                              st$npar, st$offs)
  }
  st$units
}

#' @rdname jfp_movers
#' @keywords internal
jfp_static_aug <- function(st, cols) {
  n <- st$spec@n_obs
  np <- length(st$params)
  pos <- list()
  ext <- lapply(seq_len(np), function(cc) matrix(0, n, 0L))
  q <- 0L
  for (b in cols) {
    j <- max(vapply(b$rows, function(x) if (is.null(x)) 0L else ncol(x),
                    integer(1)))
    for (cc in seq_len(np)) {
      x <- b$rows[[cc]]
      ext[[cc]] <- cbind(ext[[cc]], if (is.null(x)) matrix(0, n, j) else x)
    }
    pos[[b$name]] <- st$joint$nb + q + seq_len(j)
    q <- q + j
  }
  out <- jfp_joint_rows(st$spec, st$design0, st$coef, st$ev, extra = ext)
  nb <- st$joint$nb
  mf <- st$joint$mf
  oi <- c(seq_len(nb), nb + q + seq_len(mf - nb))
  out$pos <- pos
  out$oi <- oi
  out$pad <- function(d) {
    x <- numeric(mf + q)
    x[oi] <- d
    x
  }
  out$unit <- function(j) {
    x <- numeric(mf + q)
    x[j] <- 1
    x
  }
  out
}

#' @rdname jfp_movers
#' @param name The name of a block of columns.
#' @param f A function of a unit returning the block's static rows in the
#'   unit's equation.
#' @keywords internal
jfp_static_block <- function(st, name, f) {
  units <- jfp_movers(st)
  np <- length(st$params)
  lapply(seq_along(units), function(i) {
    un <- units[[i]]
    rows <- vector("list", np)
    rows[un$a] <- list(as_dense(f(un)))
    list(name = paste0(name, i), rows = rows)
  })
}

#' @rdname jfp_movers
#' @keywords internal
jfp_static_pairs <- function(st, name, f) {
  units <- jfp_movers(st)
  np <- length(st$params)
  lapply(seq_along(units), function(i) {
    un <- units[[i]]
    m <- length(un$ra)
    # column (k - 1) m + l: the static predictor's derivative in the unit's
    # coordinates k and l, with f's further directions
    X <- do.call(cbind, lapply(seq_len(m), function(k) {
      ek <- numeric(m)
      ek[k] <- 1
      as_dense(f(un, ek))
    }))
    rows <- vector("list", np)
    rows[un$a] <- list(X)
    list(name = paste0(name, i), rows = rows)
  })
}

#' @rdname jfp_movers
#' @param v,w Directions over every coordinate.
#' @keywords internal
jfp_static_xvw <- function(st, v, w) {
  n <- st$spec@n_obs
  rows <- vector("list", length(st$params))
  for (un in jfp_movers(st)) {
    x <- as.numeric(as_dense(modelterms7::term_block_deriv(
      un$built, coef = un$bt, v = v[un$ra])) %*% w[un$ra])
    rows[[un$a]] <- if (is.null(rows[[un$a]])) matrix(x, n, 1L) else
      rows[[un$a]] + x
  }
  list(name = "xvw", rows = rows)
}

#' @rdname jfp_movers
#' @keywords internal
jfp_static_rows <- function(st, dirs) {
  units <- jfp_movers(st)
  np <- length(st$params)
  n <- st$spec@n_obs
  mf <- st$joint$mf
  out <- lapply(seq_len(np), function(cc) matrix(0, n, mf))
  v <- dirs[[1L]]
  D1 <- function(d) function(un) modelterms7::term_block_deriv(
    un$built, coef = un$bt, v = d[un$ra])
  if (length(dirs) == 1L) {
    A <- jfp_static_aug(st, jfp_static_block(st, "Dv", D1(v)))
    for (cc in seq_len(np)) {
      for (i in seq_along(units)) {
        ra <- units[[i]]$ra
        out[[cc]][, ra] <- out[[cc]][, ra] +
          A$Vfull[[cc]][, A$pos[[paste0("Dv", i)]], drop = FALSE]
      }
    }
    return(out)
  }
  w <- dirs[[2L]]
  A <- jfp_static_aug(st, c(
    list(jfp_static_xvw(st, v, w)),
    jfp_static_block(st, "Dv", D1(v)), jfp_static_block(st, "Dw", D1(w)),
    jfp_static_block(st, "D2", function(un) modelterms7::term_block_deriv2(
      un$built, coef = un$bt, v = v[un$ra], u = w[un$ra]))))
  zero <- stats::setNames(lapply(st$params, function(p) numeric(n)),
                          st$params)
  # {v, w} as one block; {k, w} with v; {k, v} with w; {k, v, w}
  rx <- A$higher(zero, list(A$unit(A$pos$xvw)), scale = FALSE)$rows
  rv <- A$higher(zero, list(A$pad(v)), scale = FALSE)$rows
  rw <- A$higher(zero, list(A$pad(w)), scale = FALSE)$rows
  for (cc in seq_len(np)) {
    out[[cc]] <- out[[cc]] + rx[[cc]][, A$oi, drop = FALSE]
    for (i in seq_along(units)) {
      ra <- units[[i]]$ra
      out[[cc]][, ra] <- out[[cc]][, ra] +
        rv[[cc]][, A$pos[[paste0("Dw", i)]], drop = FALSE] +
        rw[[cc]][, A$pos[[paste0("Dv", i)]], drop = FALSE] +
        A$Vfull[[cc]][, A$pos[[paste0("D2", i)]], drop = FALSE]
    }
  }
  out
}

#' @rdname jfp_movers
#' @keywords internal
jfp_static_curv <- function(st, g, dirs = list()) {
  units <- jfp_movers(st)
  np <- length(st$params)
  n <- st$spec@n_obs
  mf <- st$joint$mf
  out <- matrix(0, mf, mf)
  # the Jacobian's added columns contracted against the weights
  gjac <- function(A, nm) {
    s <- 0
    for (cc in seq_len(np)) {
      gc <- g[[st$params[cc]]]
      if (is.null(gc) || !any(gc != 0)) next
      s <- s + as.numeric(crossprod(rep_len(gc, n),
                                    A$Vfull[[cc]][, A$pos[[nm]],
                                                  drop = FALSE]))
    }
    s
  }
  # a vector over a unit's coordinate pairs, added to its square block
  addsq <- function(x, i) {
    ra <- units[[i]]$ra
    m <- length(ra)
    X <- matrix(x, m, m, byrow = TRUE)
    out[ra, ra] <<- out[ra, ra] + (X + t(X)) / 2
  }
  D1 <- function(d) function(un) modelterms7::term_block_deriv(
    un$built, coef = un$bt, v = d[un$ra])
  Xkl <- function(un, ek) modelterms7::term_block_deriv(un$built,
                                                         coef = un$bt, v = ek)
  X3 <- function(d) function(un, ek) modelterms7::term_block_deriv2(
    un$built, coef = un$bt, v = d[un$ra], u = ek)
  if (!length(dirs)) {
    # {k, l}
    A <- jfp_static_aug(st, jfp_static_pairs(st, "X", Xkl))
    for (i in seq_along(units)) addsq(gjac(A, paste0("X", i)), i)
    return((out + t(out)) / 2)
  }
  v <- dirs[[1L]]
  A1 <- matrix(0, mf, mf)
  if (length(dirs) == 1L) {
    A <- jfp_static_aug(st, c(jfp_static_block(st, "Dv", D1(v)),
                              jfp_static_pairs(st, "X", Xkl),
                              jfp_static_pairs(st, "X3v", X3(v))))
    cA <- A$curv_full(g, scale = FALSE)
    pv <- A$pad(v)
    for (i in seq_along(units)) {
      ra <- units[[i]]$ra
      # {k, v} with l, and {l, v} with k by the transpose
      A1[ra, ] <- A1[ra, ] + cA[A$pos[[paste0("Dv", i)]], A$oi, drop = FALSE]
      # {k, l} with v, and {k, l, v}
      addsq(as.numeric(cA[A$pos[[paste0("X", i)]], , drop = FALSE] %*% pv) +
              gjac(A, paste0("X3v", i)), i)
    }
    out <- out + A1 + t(A1)
    return((out + t(out)) / 2)
  }
  w <- dirs[[2L]]
  X4 <- function(un, ek) modelterms7::term_block_deriv3(
    un$built, coef = un$bt, v = v[un$ra], u = w[un$ra], s = ek)
  A <- jfp_static_aug(st, c(
    list(jfp_static_xvw(st, v, w)),
    jfp_static_block(st, "Dv", D1(v)), jfp_static_block(st, "Dw", D1(w)),
    jfp_static_block(st, "D2", function(un) modelterms7::term_block_deriv2(
      un$built, coef = un$bt, v = v[un$ra], u = w[un$ra])),
    jfp_static_pairs(st, "X", Xkl), jfp_static_pairs(st, "X3v", X3(v)),
    jfp_static_pairs(st, "X3w", X3(w)), jfp_static_pairs(st, "X4", X4)))
  pv <- A$pad(v)
  pw <- A$pad(w)
  hx <- A$higher(g, list(A$unit(A$pos$xvw)), scale = FALSE)$curv
  hv <- A$higher(g, list(pv), scale = FALSE)$curv
  hw <- A$higher(g, list(pw), scale = FALSE)$curv
  cA <- A$curv_full(g, scale = FALSE)
  # {v, w} with k and l apart
  out <- out + hx[A$oi, A$oi]
  ra_all <- integer(0)
  pv_all <- integer(0)
  pw_all <- integer(0)
  for (i in seq_along(units)) {
    ra <- units[[i]]$ra
    # {k, v} with l and w, {k, w} with l and v, {k, v, w} with l, and the
    # same with k and l exchanged by the transpose
    A1[ra, ] <- A1[ra, ] +
      hw[A$pos[[paste0("Dv", i)]], A$oi, drop = FALSE] +
      hv[A$pos[[paste0("Dw", i)]], A$oi, drop = FALSE] +
      cA[A$pos[[paste0("D2", i)]], A$oi, drop = FALSE]
    ra_all <- c(ra_all, ra)
    pv_all <- c(pv_all, A$pos[[paste0("Dv", i)]])
    pw_all <- c(pw_all, A$pos[[paste0("Dw", i)]])
    # k and l in one block: {k, l} with v and w, {k, l} with {v, w},
    # {k, l, v} with w, {k, l, w} with v, {k, l, v, w}
    addsq(as.numeric(hw[A$pos[[paste0("X", i)]], , drop = FALSE] %*% pv) +
            as.numeric(cA[A$pos[[paste0("X", i)]], A$pos$xvw]) +
            as.numeric(cA[A$pos[[paste0("X3v", i)]], , drop = FALSE] %*% pw) +
            as.numeric(cA[A$pos[[paste0("X3w", i)]], , drop = FALSE] %*% pv) +
            gjac(A, paste0("X4", i)), i)
  }
  out <- out + A1 + t(A1)
  # {k, v} with {l, w}, and {k, w} with {l, v}
  C <- cA[pv_all, pw_all, drop = FALSE]
  out[ra_all, ra_all] <- out[ra_all, ra_all] + C + t(C)
  (out + t(out)) / 2
}
