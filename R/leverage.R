#' @include sandwich.R
NULL

#' The Per-Observation Pieces of the Hat Matrix of a Fit
#'
#' @description
#' Assembles, over the coordinates `cols` of the coefficient vector, the
#' stacked design of the observations, the information of each observation
#' and its symmetric square root. [leverage_blocks()] and the corrected
#' sandwich variances of [hc2_sandwich()] are computed from these pieces.
#'
#' @details
#' Let \eqn{d} be the number of distribution parameters that have at least
#' one coefficient. Observation \eqn{i} contributes the \eqn{d \times P}
#' matrix \eqn{Z_i}, whose row \eqn{a} holds the row of the design of
#' equation \eqn{a} in the columns of that equation and zeros elsewhere, so
#' that the score of the observation in the coefficients is
#' \eqn{u_i = Z_i' s_i}, with \eqn{s_i} the vector of its derivatives of the
#' log-likelihood in the \eqn{d} link-scale predictors. The matrix
#' \eqn{W_i} is the \eqn{d \times d} information of the observation in those
#' predictors (the negative second derivative, or its expectation), times the
#' prior weight. It follows `expected`: the expected information (with the
#' method of the fit's information when that is expected, and the closed form
#' of the family otherwise) or the observed one. The expected information of
#' one observation is the covariance of its score under the model and is
#' positive semi-definite; the observed information of one observation can be
#' indefinite. Its symmetric square root \eqn{W_i^{1/2}} is taken from the
#' eigen decomposition of the positive part (eigenvalues below
#' \eqn{10^{-12}} times the largest in absolute value are set to zero), and the
#' number of observations with a negative eigenvalue is returned as `clipped`.
#' The pseudo-inverse of the root is kept beside it, since an observation can
#' have a singular \eqn{W_i}.
#'
#' @param fit The [StatmodFit()].
#' @param design The design of the fit.
#' @param cols The positions, among the coefficients stacked by equation, of
#'   the coordinates that the variance covers.
#' @param expected,emethod As in [statmod_information_at()]: whether the
#'   information is the expected one, and the method that obtains it.
#'
#' @return A list with `params` (the names of the \eqn{d} equations that have
#'   coefficients), `n`, `d`, `Z` (a sparse matrix with \eqn{nd} rows, ordered
#'   by observation and then by equation, and one column for each element of
#'   `cols`), `W`, `root` and `root_pinv` (arrays of dimension
#'   \eqn{n \times d \times d}) and `clipped` (how many observations had a
#'   negative eigenvalue set to zero).
#'
#' @keywords internal
leverage_parts <- function(fit, design, cols, expected, emethod) {
  spec <- fit@spec
  params <- spec@distrib@params
  n <- spec@n_obs
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  act <- which(npar > 0L)
  d <- length(act)
  ev <- statmod_eta(spec, design, fit@coefficients)
  # the roots follow the information of the fit: the expected one (the working
  # covariance of the scores, positive semi-definite) or the observed one, whose
  # per-observation matrix can be indefinite and is then replaced by its
  # positive part in psd_roots(), the count being returned as `clipped`
  H <- if (expected) statmod_family_hessian(spec, ev$theta, TRUE, emethod) else
    statmod_family_hessian(spec, ev$theta, FALSE, NULL)
  w <- spec@weights
  W <- array(0, c(n, d, d))
  for (a in seq_len(d)) {
    for (b in a:d) {
      v <- -w * rep_len(H[[hess_key(params, act[a], act[b])]], n)
      W[, a, b] <- v
      W[, b, a] <- v
    }
  }
  rt <- psd_roots(W)
  ii <- jj <- xx <- vector("list", d)
  for (a in seq_len(d)) {
    q <- act[a]
    stack <- offs[[q]] + seq_len(npar[[q]])
    loc <- which(stack %in% cols)
    gj <- match(stack[loc], cols)
    X <- design[[params[q]]]$X[, loc, drop = FALSE]
    Xs <- if (isS4(X)) X else Matrix::Matrix(X, sparse = TRUE)
    sm <- Matrix::summary(methods::as(methods::as(Xs, "dMatrix"),
                                      "generalMatrix"))
    ii[[a]] <- (sm$i - 1L) * d + a
    jj[[a]] <- gj[sm$j]
    xx[[a]] <- sm$x
  }
  Z <- Matrix::sparseMatrix(i = unlist(ii), j = unlist(jj), x = unlist(xx),
                            dims = c(n * d, length(cols)), repr = "C")
  list(params = params[act], n = n, d = d, Z = Z, W = W,
       root = rt$root, root_pinv = rt$root_pinv, clipped = rt$clipped)
}


#' Symmetric Square Roots of a Stack of Small Information Matrices
#'
#' @description
#' For each slice of an array of symmetric matrices, computes the symmetric
#' square root and its pseudo-inverse from the eigen decomposition.
#'
#' @details
#' Eigenvalues below \eqn{10^{-12}} times the largest eigenvalue of the slice
#' are set to zero and have a zero pseudo-inverse. This covers the negative
#' eigenvalues of an observed information and the zero information of an
#' observation with weight zero.
#'
#' @param W An array of dimension \eqn{n \times d \times d}.
#'
#' @return A list with `root` and `root_pinv`, arrays of the same dimension,
#'   and `clipped`, the number of slices with a negative eigenvalue.
#'
#' @keywords internal
psd_roots <- function(W) {
  n <- dim(W)[[1L]]
  d <- dim(W)[[2L]]
  if (d == 1L) {
    v <- W[, 1L, 1L]
    top <- max(abs(v), 0)
    ok <- v > 1e-12 * top
    r <- ifelse(ok, sqrt(pmax(v, 0)), 0)
    ri <- ifelse(ok, 1 / sqrt(ifelse(ok, v, 1)), 0)
    return(list(root = array(r, c(n, 1L, 1L)),
                root_pinv = array(ri, c(n, 1L, 1L)),
                clipped = sum(v < -1e-12 * top)))
  }
  R <- Rp <- array(0, dim(W))
  clipped <- 0L
  for (i in seq_len(n)) {
    e <- eigen(matrix(W[i, , ], d, d), symmetric = TRUE)
    lam <- e$values
    ok <- lam > 1e-12 * max(abs(lam), 0)
    if (any(lam < -1e-12 * max(abs(lam), 0))) clipped <- clipped + 1L
    sq <- ifelse(ok, sqrt(pmax(lam, 0)), 0)
    sqi <- ifelse(ok, 1 / sqrt(ifelse(ok, lam, 1)), 0)
    R[i, , ] <- e$vectors %*% (sq * t(e$vectors))
    Rp[i, , ] <- e$vectors %*% (sqi * t(e$vectors))
  }
  list(root = R, root_pinv = Rp, clipped = clipped)
}


#' The Leverage Block of Each Observation
#'
#' @description
#' Computes \eqn{H_{ii} = W_i^{1/2} Z_i V_b Z_i' W_i^{1/2}} for every
#' observation, from the pieces of [leverage_parts()] and the inverse
#' penalized information \eqn{V_b}.
#'
#' @param parts The list that [leverage_parts()] returns.
#' @param Vb The inverse penalized information over the coordinates of
#'   `parts$Z`, a dense matrix.
#' @param chunk How many observations to process at a time.
#'
#' @return An array of dimension \eqn{n \times d \times d}.
#'
#' @keywords internal
leverage_blocks <- function(parts, Vb, chunk = 4000L) {
  n <- parts$n
  d <- parts$d
  Q <- array(0, c(n, d, d))
  starts <- seq(1L, n, by = chunk)
  for (s0 in starts) {
    obs <- s0:min(n, s0 + chunk - 1L)
    m <- length(obs)
    ridx <- as.vector(t(outer((obs - 1L) * d, seq_len(d), "+")))
    Zc <- parts$Z[ridx, , drop = FALSE]
    ZV <- as.matrix(Zc %*% Vb)
    Zd <- as.matrix(Zc)
    for (a in seq_len(d)) {
      ra <- seq(a, m * d, by = d)
      for (b in a:d) {
        rb <- seq(b, m * d, by = d)
        v <- rowSums(ZV[ra, , drop = FALSE] * Zd[rb, , drop = FALSE])
        Q[obs, a, b] <- v
        Q[obs, b, a] <- v
      }
    }
  }
  Rt <- parts$root
  out <- array(0, c(n, d, d))
  for (a in seq_len(d)) for (b in seq_len(d)) {
    acc <- numeric(n)
    for (cc in seq_len(d)) for (e in seq_len(d)) {
      acc <- acc + Rt[, a, cc] * Q[, cc, e] * Rt[, e, b]
    }
    out[, a, b] <- acc
  }
  out
}


#' Leverage of the Observations in a Distributional Regression
#'
#' @description
#' Returns, for each observation \eqn{i}, the block \eqn{H_{ii}} of the hat
#' matrix of the fit over the coefficients of all the equations. For a model
#' with one equation it is the usual leverage; for a model with \eqn{d}
#' equations it is a \eqn{d \times d} matrix, and its trace is the number of
#' effective parameters that observation \eqn{i} uses.
#'
#' @details
#' Let \eqn{\hat\beta} be the stacked coefficients of the equations that have
#' coefficients, of length \eqn{P}. Write \eqn{Z_i} for the \eqn{d \times P}
#' matrix whose row \eqn{a} holds the row of the design of equation \eqn{a}
#' for observation \eqn{i}, in the columns of that equation. The score of
#' observation \eqn{i} in the coefficients is \eqn{u_i = Z_i' s_i}, where
#' \eqn{s_i} is the vector of the derivatives of its log-likelihood in the
#' \eqn{d} predictors (on the link scale of each parameter), multiplied by the
#' prior weight. Let \eqn{W_i} be the \eqn{d \times d} information of the
#' observation in those predictors, times the prior weight, and \eqn{K} the
#' penalized information \eqn{\sum_i Z_i' W_i Z_i + S} that the Bayesian
#' variance \eqn{V_b = K^{-1}} inverts, with \eqn{S} the Hessian of the
#' penalty. The linearized fit satisfies
#' \eqn{W^{-1/2}\hat s \approx (I - H)\, W^{-1/2} s}, with the hat matrix
#' \deqn{H = W^{1/2} Z K^{-1} Z' W^{1/2},}
#' where \eqn{Z} stacks the \eqn{Z_i}, \eqn{W} is block diagonal in the
#' \eqn{W_i} and \eqn{W^{1/2}} is the symmetric square root. The diagonal
#' block of observation \eqn{i} is
#' \deqn{H_{ii} = W_i^{1/2} Z_i K^{-1} Z_i' W_i^{1/2}.}
#' It is symmetric, and for an unpenalized fit with \eqn{W} the information
#' at the estimate it is a projection, so that \eqn{\sum_i \mathrm{tr}\,
#' H_{ii} = P}. With penalties the sum is the effective degrees of freedom.
#'
#' For a Gaussian linear model with a constant scale and the expected
#' information, the block of the mean equation is the scalar
#' \eqn{h_i = x_i'(X'X)^{-1}x_i} of [stats::hatvalues()]. For one equation
#' of a larger model, the diagonal element \eqn{(a, a)} of the block is the
#' leverage that [plot_influence()] computes for that equation, provided that
#' `information` is [distributions7::expected_info()]. Both read
#' \eqn{K^{-1}} and not the Bayesian variance of [vcov.StatmodFit()], which
#' for a coordinate that the marginal criterion estimates holds the
#' curvature of that criterion.
#'
#' The matrix \eqn{W_i} and the matrix \eqn{K} are both read on `information`.
#' With [distributions7::expected_info()] the hat matrix is a projection and
#' its trace is \eqn{P} when there is no penalty. With
#' [distributions7::observed_info()] the matrix \eqn{W_i} is the observed
#' information of the observation, which can be indefinite; its root is then
#' the symmetric square root of the positive part (the eigenvalues that are
#' not positive are set to zero, with a zero pseudo-inverse), `clipped`
#' counts the observations where that happened, and the trace is close to
#' \eqn{P} but not equal to it. With prior weights \eqn{w_i} the block is the
#' one of the row as a unit, \eqn{w_i} times the block of one of its copies.
#'
#' @param fit A [StatmodFit()].
#' @param information An object from [distributions7::observed_info()] or
#'   [distributions7::expected_info()], or `NULL` for the information that
#'   [vcov.StatmodFit()] resolves for the family.
#'
#' @return A list with `blocks`, an array of dimension \eqn{n \times d \times
#'   d} of the matrices \eqn{H_{ii}}, `leverage`, the \eqn{n \times d} matrix
#'   of their diagonals with the equations as columns, and `trace`, the
#'   vector of the traces, and `clipped`, the number of observations whose
#'   information had a negative eigenvalue set to zero (always zero with the
#'   expected information). The \eqn{d} equations are those that have
#'   coefficients. Coordinates whose variance is missing (aliased, or held at
#'   zero by a penalty) are not part of \eqn{Z_i}.
#'
#' @references
#' Hoaglin, D. C. and Welsch, R. E. (1978). The hat matrix in regression and
#' ANOVA. *The American Statistician*, 32(1), 17--22.
#'
#' @seealso [hc2_sandwich()], which uses the blocks; [plot_influence()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.3)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' lv <- statmod_leverage(fit, information = distributions7::expected_info())
#' range(lv$leverage[, "mu"] - hatvalues(lm(y ~ x, dd)))
#' @export
statmod_leverage <- function(fit, information = NULL) {
  if (!S7::S7_inherits(fit, StatmodFit)) {
    stop("'fit' must be a fitted statmod() model.", call. = FALSE)
  }
  vt <- bayesian_vcov(information = information)
  spec <- fit@spec
  design <- statmod_design(spec)
  if (length(attr(design, "structural")) > 0L ||
      length(attr(design, "refresh")) > 0L) {
    stop(paste0("The leverage is not available for a model with a structural ",
                "term or with a term whose block moves with its coefficients."),
         call. = FALSE)
  }
  V <- vcov(fit, vcov_type = vt, readable = FALSE)
  cols <- which(is.finite(diag(V)))
  info <- vcov_information(vt, spec@distrib)
  parts <- leverage_parts(fit, design, cols, info_expected(info),
                          info_emethod(info))
  # K^-1 itself: the Bayesian variance above can carry the curvature of the
  # marginal criterion in place of K^-1 for a coordinate that criterion
  # estimates (the scale of a Gaussian under REML), and the hat matrix is the
  # projection that K^-1 defines
  coef <- fit@coefficients
  K <- statmod_information_at(spec, coef, design, info_expected(info),
                              info_emethod(info)) +
    statmod_penalty_at(spec, coef, fit@hyper, design, "hessian")
  K <- as_dense(K)[cols, cols, drop = FALSE]
  Vb <- as_dense(solve_pd(K, "the penalized information", rownames(V)[cols]))
  blocks <- leverage_blocks(parts, Vb)
  dimnames(blocks) <- list(NULL, parts$params, parts$params)
  lev <- vapply(seq_len(parts$d), function(a) blocks[, a, a], numeric(parts$n))
  lev <- matrix(lev, parts$n, parts$d, dimnames = list(NULL, parts$params))
  tr <- rowSums(lev)
  list(blocks = blocks, leverage = lev, trace = tr, clipped = parts$clipped)
}
