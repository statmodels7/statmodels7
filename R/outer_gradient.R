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
#' \code{TRUE} when every hyperparameter under estimation belongs to a penalty
#' whose second derivative in the coefficients is linear in the
#' hyperparameters and free of the coefficients, and the criterion uses the
#' observed information.
#'
#' @details
#' \strong{Why the observed information.} \eqn{K} enters the criterion through
#' its determinant, so the gradient needs \eqn{\partial K/\partial\beta}. With
#' the observed information that is the third derivative of the log-likelihood
#' in the link-scale predictors, which every family of \pkg{distributions7}
#' carries in closed form. With the expected information it would be the
#' derivative in \eqn{\beta} of \eqn{-E[\ell'']}, which is not
#' \eqn{-E[\ell''']} and is not one of that package's generics. So
#' \code{reml(hessian = "observed")} is what the exact route asks for, and
#' \code{"expected"} keeps the derivative-free search.
#'
#' \strong{Why the penalty is asked.} \eqn{\partial S/\partial\theta} and its
#' second derivative are generics of \pkg{penalties7}
#' (\code{\link[penalties7]{penalty_dhessian}},
#' \code{\link[penalties7]{penalty_d2hessian}},
#' \code{\link[penalties7]{penalty_dcross}}). A penalty that answers them is
#' estimable by a marginal criterion whatever its shape: the quadratic, the
#' additive, the structured and the separable branches all do, so a ridge, a
#' random effect and a heavy-tailed prior are covered as well as a spline. One
#' that does not -- a SCAD, an MCP, anything with a kink -- rejects, and the
#' search stays derivative-free. Nothing here tests a penalty's behaviour to
#' find out what it is; it is asked.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param idx The outer index, from \code{\link{outer_hyper_index}}.
#' @param method An \code{\link{OuterMethod}}.
#' @param order \code{1} for the gradient, \code{2} for the Hessian as well.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @seealso \code{\link{statmod_marginal_grad}}
#'
#' @keywords internal
outer_gradient_ok <- function(spec, design, idx, method, order = 1L) {
  if (!identical(method@hessian, "observed")) return(FALSE)
  if (!nrow(idx)) return(FALSE)
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
  if (structural_penalized(spec, design)) {
    if (order >= 2L) return(FALSE)
    for (u in statmod_penalized(spec, design)) {
      if (!isTRUE(u$structural)) next
      if (!answers_term_third(spec@terms[[u$param]][[u$term]])) return(FALSE)
    }
  }
  TRUE
}


#' Does a Term Supply Its Third Derivative?
#'
#' @description
#' Whether the term implements \code{\link[modelterms7]{term_third}}, read
#' from the class the method is registered on rather than from a list of
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
#' @seealso \code{\link{outer_gradient_ok}}
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
#' coefficients (\code{\link[penalties7]{beta_quadratic}}), since otherwise the
#' third and fourth derivatives of the penalty in \eqn{\beta} enter the
#' criterion's own second derivative and \pkg{penalties7} does not carry them.
#'
#' @param pen A \pkg{penalties7} penalty.
#' @param order \code{1} or \code{2}.
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
#' \eqn{\ell'''_{abk}} is looked up by a name BUILT from the parameter names in
#' the family's own order, never parsed out of one.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An \code{\link{OuterMethod}}.
#' @param idx The outer index.
#' @param basis The integrated subspace, or \code{NULL}.
#' @param free Whether to carry the result onto the free scale. The Hessian
#'   asks for the parameter scale, having its own second-order chain rule to
#'   apply.
#'
#' @return A numeric vector, one entry per row of \code{idx}, or \code{NULL}
#'   where the determinant does not exist.
#'
#' @seealso \code{\link{statmod_marginal}}, \code{\link{reml}}
#'
#' @keywords internal
statmod_marginal_grad <- function(spec, design, coef, hyper, method, idx,
                                  basis = NULL, free = TRUE) {
  if (structural_penalized(spec, design)) {
    return(statmod_structural_grad(spec, design, coef, hyper, method, idx,
                                   basis, free))
  }
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)

  H <- statmod_information_at(spec, coef, design, expected = FALSE)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian")
  S[!is.finite(S)] <- 0
  # The routes below form a FULL inverse, and the inverse of a sparse matrix
  # is dense, so densifying here gives up nothing sparsity could have kept.
  K <- as_dense(H + S)
  Kfac <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(Kfac)) return(NULL)
  Kinv <- chol2inv(Kfac)
  M <- if (is.null(basis)) Kinv else {
    A <- basis
    inner <- tryCatch(chol2inv(chol(crossprod(A, K %*% A))),
                      error = function(e) NULL)
    if (is.null(inner)) return(NULL)
    A %*% inner %*% t(A)
  }

  u <- u_vector(spec, design, coef, M, params, npar, offs, total)

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
      cols <- un$cols
      pos <- un$index
      pen <- un$penalty
      bt <- coef[[p]][cols]
      th <- as.list(hyper[[p]][[nm]])
      gt <- penalties7::penalty_grad_theta(pen, bt, th)
      cr <- penalties7::penalty_cross(pen, bt, th)
      dS <- penalties7::penalty_dhessian(pen, bt, th)
      for (r in rows[idx$term[rows] == nm]) {
        h <- idx$name[r]
        # the mode moves by -(H+S)^-1 d2rho/dbeta dtheta
        c_m <- numeric(total)
        c_m[pos] <- as.numeric(cr[[h]])
        v <- -as.numeric(Kinv %*% c_m)
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
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients.
#' @param M The matrix the trace is taken against.
#' @param params The distribution's parameter names.
#' @param npar,offs,total The block sizes, their offsets and the total.
#'
#' @return A numeric vector as long as the stacked coefficients.
#'
#' @keywords internal
u_vector <- function(spec, design, coef, M, params, npar, offs, total) {
  n <- spec@n_obs
  th <- statmod_eta(spec, design, coef)$theta
  d3 <- distributions7::distrib_deriv3(spec@distrib, spec@response, th,
                                       scale = "link")
  keys <- names(d3)

  # the per-observation diagonal of each block of M
  G <- vector("list", length(params))
  for (a in seq_along(params)) {
    G[[a]] <- vector("list", length(params))
    if (npar[a] == 0L) next
    Xa <- design[[params[a]]]$X
    for (b in seq_along(params)) {
      if (npar[b] == 0L) next
      Mab <- M[offs[a] + seq_len(npar[a]), offs[b] + seq_len(npar[b]),
               drop = FALSE]
      G[[a]][[b]] <- rowSums((Xa %*% Mab) * design[[params[b]]]$X)
    }
  }

  out <- numeric(total)
  for (k in seq_along(params)) {
    if (npar[k] == 0L) next
    s <- numeric(n)
    for (a in seq_along(params)) {
      if (npar[a] == 0L) next
      for (b in seq_along(params)) {
        if (npar[b] == 0L) next
        s <- s + rep_len(d3[[d3_key(params, a, b, k, keys)]], n) * G[[a]][[b]]
      }
    }
    out[offs[k] + seq_len(npar[k])] <-
      -as.numeric(crossprod(design[[params[k]]]$X, spec@weights * s))
  }
  out
}


#' The Name of a Third-Derivative Component
#'
#' @description
#' Locates the \eqn{(a, b, k)} entry of a distribution's third derivative,
#' which is keyed by name and not by position.
#'
#' @details
#' The name is BUILT by putting the three parameter names in the family's own
#' order and joining them, the direction \pkg{distributions7} sanctions, and
#' then checked against the enumeration rather than trusted.
#'
#' @param params The parameter names, in the family's order.
#' @param a,b,k Indices into \code{params}.
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
#' It exists so that the contraction \code{\link{u_vector}} performs is
#' written once. The formula there does not change where a filter is present;
#' only its operand does, \eqn{X} becoming \eqn{[X_p \mid D]}. Everything
#' downstream -- the third derivative against the diagonal of \eqn{M}, the
#' movement of the mode -- then reads the joint vector without a special case.
#'
#' The rows are returned at the FULL width, a held level included, and the
#' caller drops it exactly as \code{\link{statmod_full_information}} does.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients.
#'
#' @return A list with \code{V} (one matrix per distribution parameter),
#'   \code{ap} (which equation carries the filter), \code{keep} (the columns
#'   that survive a held level), \code{ev}, \code{f} and the sizes.
#'
#' @seealso \code{\link{u_vector}}, \code{\link{statmod_full_information}}
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
#' \code{\link{statmod_marginal_grad}} over the joint vector of coefficients
#' and a structural term's parameters, which is what the determinant spans
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
#' \enumerate{
#'   \item the family's third derivative against the per-observation diagonal
#'     of \eqn{M}, which is \code{\link{u_vector}}'s formula with \eqn{V} in
#'     place of \eqn{X};
#'   \item the derivative of \eqn{V_p} itself, which is \eqn{E_t v}, one row
#'     per observation, entering through every \eqn{\ell_{pb}};
#'   \item the derivative of \eqn{\ell_p E_t}, whose first half is \eqn{E}
#'     re-weighted by \eqn{\sum_k \ell_{pk}(V_k\cdot v)} and whose second is
#'     \eqn{\partial^3 e_t/\partial u^3[v]}, the one genuinely new object.
#' }
#'
#' Piece 1 is direction-free and is computed once; pieces 2 and 3 are read
#' along each hyperparameter's own direction, which is the trade that keeps
#' the cost at \eqn{O(nm^2)} per hyperparameter instead of \eqn{O(nm^3)} once.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An \code{\link{OuterMethod}}.
#' @param idx The outer index.
#' @param basis The integrated subspace, or \code{NULL}.
#' @param free Whether to carry the result onto the free scale.
#'
#' @return A numeric vector, one entry per row of \code{idx}, or \code{NULL}
#'   where the determinant does not exist.
#'
#' @seealso \code{\link{statmod_marginal_grad}},
#'   \code{\link[modelterms7]{term_third}}
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
#' The quantities of \code{\link{statmod_structural_grad}} that do not depend
#' on which hyperparameter is being differentiated: the family's derivatives,
#' the filter's forward Jacobian, the per-observation diagonal of \eqn{M} and
#' the contraction \eqn{u}.
#'
#' @details
#' Computing them once is what keeps a model with several hyperparameters
#' affordable: only the two pieces that read the direction are repeated, and
#' each of those is one pass of the recursion.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param coef The coefficients.
#' @param jd The joint rows, from \code{\link{joint_design_rows}}.
#' @param M The matrix the trace is taken against.
#'
#' @return A list of the shared quantities.
#'
#' @keywords internal
structural_grad_parts <- function(spec, design, coef, jd, M) {
  params <- jd$params
  n <- jd$n
  w <- spec@weights
  keep <- jd$keep
  ap <- jd$ap
  ev <- jd$ev
  f <- jd$f
  d <- spec@distrib
  y <- spec@response
  gl <- distributions7::distrib_gradient(d, y, ev$theta, scale = "link")
  H <- distributions7::distrib_hessian(d, y, ev$theta, scale = "link")
  D3 <- distributions7::distrib_deriv3(d, y, ev$theta, scale = "link")
  D4 <- distributions7::distrib_deriv4(d, y, ev$theta, scale = "link")
  s_at <- rep_len(gl[[f$param]], n)
  c_at <- rep_len(H[[hess_key(params, ap, ap)]], n)

  # the filter's forward Jacobian, which is V for the equation it sits in.
  # The recursion is re-run at parameters it has already been run at, so the
  # callbacks LOOK UP the derivatives read at that predictor rather than
  # asking the family again -- the same bargain the information makes.
  seed <- jd$V[[ap]]
  mk_blocks <- .structural_blocks(params, ap, jd$V, H, D3, D4, n)
  cv <- modelterms7::term_curvature(
    f$tm, f$eta_static, y, function(e, i) s_at[i], function(e, i) c_at[i],
    f$psi, w * s_at, seed, mk_blocks(NULL))
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
       s_at = s_at, c_at = c_at, seed = seed, blocks = mk_blocks, w = w)
}


#' The Two Pieces of the Chain Term That Read the Direction
#'
#' @description
#' The contributions to \eqn{\mathrm{tr}(M\,\partial K/\partial u[v])} that
#' come from the recursion rather than from the design: the derivative of the
#' filter's own Jacobian, and the derivative of the term the level contributes
#' to the information.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param jd The joint rows.
#' @param M The matrix the trace is taken against.
#' @param st The shared quantities, from \code{\link{structural_grad_parts}}.
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
    w * kappa, st$seed, st$blocks(NULL))
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
#' Builds the \code{blocks} callback \code{\link[modelterms7]{term_curvature}}
#' and \code{\link[modelterms7]{term_third}} take, at a direction or without
#' one.
#'
#' @details
#' The pieces are built on the ACTIVE SET the term asks for, so a panel's
#' outer products are of the same size whatever the number of groups.
#' \code{dcurv} serves twice, as the derivative of the curvature along the
#' direction and as the factor multiplying the movement of \eqn{V_p}, and
#' \code{N} is where the family's fourth derivative enters: each order of
#' differentiating the predictor through the recursion pulls in one more order
#' of the family.
#'
#' @param params The distribution's parameter names.
#' @param ap Which of them carries the filter.
#' @param Vs The static rows.
#' @param H,D3,D4 The family's derivatives at the fitted predictors.
#' @param n The number of observations.
#'
#' @return A function of the direction returning a \code{blocks} callback.
#'
#' @keywords internal
.structural_blocks <- function(params, ap, Vs, H, D3, D4, n) {
  at <- function(x, i) rep_len(x, n)[i]
  function(vfull) {
    third <- !is.null(vfull)
    function(e, i, D, act = NULL) {
      if (is.null(act)) act <- seq_len(ncol(Vs[[1L]]))
      mk <- length(act)
      cross <- numeric(mk)
      for (q in seq_along(params)) {
        if (q == ap) next
        cross <- cross + at(H[[hess_key(params, ap, q)]], i) * Vs[[q]][i, act]
      }
      vr <- lapply(seq_along(params), function(r)
        if (r == ap) D else Vs[[r]][i, act])
      M <- matrix(0, mk, mk)
      for (r in seq_along(params)) {
        for (r2 in seq_along(params)) {
          M <- M + at(D3[[deriv3_key(params, ap, r, r2)]], i) *
            outer(vr[[r]], vr[[r2]])
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
        dcurv <- dcurv + at(D3[[deriv3_key(params, ap, ap, r)]], i) * vr[[r]]
      }
      N <- matrix(0, mk, mk)
      for (r in seq_along(params)) {
        for (r2 in seq_along(params)) {
          co <- 0
          for (r3 in seq_along(params)) {
            co <- co + at(D4[[deriv4_key(params, ap, r, r2, r3)]], i) * dv[r3]
          }
          if (co != 0) N <- N + co * outer(vr[[r]], vr[[r2]])
        }
      }
      list(cross = cross, M = M, dcurv = dcurv, N = N)
    }
  }
}


#' The Subspace a Marginal Criterion Integrates Over, Jointly
#'
#' @description
#' The basis \code{\link{ml}()} projects the joint curvature onto: the
#' coefficients' own range basis, with one column per penalized parameter of
#' the structural term appended.
#'
#' @details
#' It is written once because \code{\link{statmod_marginal_full}} and
#' \code{\link{statmod_structural_grad}} must project onto the SAME subspace,
#' and two callers composing it separately would agree only by accident.
#'
#' @param spec A \code{\link{StatmodSpec}}.
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
