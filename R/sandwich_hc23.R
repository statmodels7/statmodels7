#' @include sandwich.R leverage.R
NULL

#' @title Sandwich Variances with a Leverage Correction
#'
#' @description
#' Select the meat of the frequentist variance \eqn{V = V_b B V_b} with the
#' scores of the observations rescaled by their leverage: `hc2_sandwich()`
#' and `hc3_sandwich()`, which with a `cluster` are the cluster-robust
#' estimators CR2 and CR3, and `hc4_sandwich()`. Pass the object as the
#' `sandwich` argument of [frequentist_vcov()].
#'
#' @details
#' Let \eqn{s_i} be the vector of the derivatives of the log-likelihood of
#' observation \eqn{i} in the \eqn{d} predictors of the model, \eqn{Z_i} the
#' \eqn{d \times P} matrix that places the design rows of the observation in
#' the stacked coefficients, \eqn{W_i} the \eqn{d \times d} information of the
#' observation, and \eqn{H_{ii}} its leverage block, as defined at
#' [statmod_leverage()]. The score of the observation is
#' \eqn{u_i = Z_i' s_i}. Each estimator replaces \eqn{s_i} by
#' \deqn{s_i^* = W_i^{1/2}\, (I - H_{ii})^{-p}\, W_i^{-1/2}\, s_i,}
#' and then forms \eqn{B = \sum_i u_i^* u_i^{*\prime}} with
#' \eqn{u_i^* = Z_i' s_i^*}. The estimators differ in the exponent \eqn{p}:
#' \itemize{
#'   \item HC2 (`hc2_sandwich()`): \eqn{p = 1/2}, as in MacKinnon and White
#'     (1985);
#'   \item HC3 (`hc3_sandwich()`): \eqn{p = 1}, as in MacKinnon and White
#'     (1985);
#'   \item HC4 (`hc4_sandwich()`): \eqn{p = \delta_i/2} with
#'     \eqn{\delta_i = \min\{4, \mathrm{tr}(H_{ii})/\bar h\}}, where
#'     \eqn{\bar h} is the mean of \eqn{\mathrm{tr}(H_{ii})} over the
#'     observations, as in Cribari-Neto (2004).
#' }
#' The matrix power is taken on the eigenvalues of the symmetric matrix
#' \eqn{H_{ii}}. When a model has one equation, \eqn{W_i} is a scalar and
#' \eqn{s_i^* = s_i (1 - h_i)^{-p}}, the textbook estimators; for a Gaussian
#' linear model with a constant scale they equal the estimators built from
#' [stats::hatvalues()]. An eigenvalue of \eqn{H_{ii}} equal to one (a
#' direction that the fit reproduces exactly) is given the factor zero, which
#' is the pseudo-inverse, and a warning reports it.
#'
#' With `cluster = ~ id`, the observations of a cluster \eqn{g} are corrected
#' together: \eqn{s_g^* = W_g^{1/2}\,(I - H_{gg})^{-p}\,W_g^{-1/2}\,s_g},
#' where \eqn{s_g} stacks the vectors \eqn{s_i} of its observations,
#' \eqn{W_g} is block diagonal and \eqn{H_{gg}} is the diagonal block of the
#' hat matrix for the cluster, of dimension equal to the number of
#' observations in the cluster times \eqn{d}. The exponent \eqn{1/2} is
#' CR2 (Bell and McCaffrey, 2002; Pustejovsky and Tipton, 2018) and the
#' exponent \eqn{1} is CR3. A cluster with one observation each reduces them
#' to HC2 and HC3. `hc4_sandwich()` has no `cluster` argument, since no
#' version of HC4 for clusters is established. The estimators apply no
#' degrees-of-freedom factor, so that `hc2_sandwich()` and `hc3_sandwich()`
#' coincide with [hc0_sandwich()] when the leverage is zero.
#'
#' The correction is the exact bias removal for the linearized fit when the
#' working model is the model itself (the observations independent, with
#' variance \eqn{W_i^{-1}} for their scores), and it uses the information
#' selected by the `information` argument of [frequentist_vcov()] when it is
#' given, for the roots \eqn{W_i^{1/2}} as well as for the inverse
#' \eqn{K^{-1}}. When that argument is `NULL`, the leverage is computed on the
#' expected information, whatever information the bread \eqn{V_b} resolves
#' for the family. With the observed information, an observation whose matrix
#' \eqn{W_i} is indefinite contributes the positive part, with the
#' pseudo-inverse root ([statmod_leverage()] returns how many observations
#' that concerns). With penalties the leverage is that of
#' the penalized fit, so that the sum of the traces is the effective degrees of
#' freedom of the model. Prior weights, random effects, aliased coordinates
#' and the models that are rejected are treated as at [hc0_sandwich()]. Prior
#' weights enter the leverage through \eqn{W_i}; without a cluster the leverage
#' of one of the \eqn{w_i} copies of a row is \eqn{1/w_i} times that of the
#' row, which is what the estimator uses, and with a cluster the copies of a
#' row are corrected within the cluster. With `weights = "sampling"` and no
#' cluster the leverage is that of the row, not divided by \eqn{w_i}, and the
#' mean leverage of HC4 is the mean over the rows with positive weight.
#'
#' **Few clusters.** With CR2 and `reference = "t"` (the default),
#' [confint.StatmodFit()] and the Wald rows of [summary.StatmodFit()] and
#' [statmod_test()] use a \eqn{t} distribution whose degrees of freedom are
#' computed for each coefficient by the Satterthwaite approximation of
#' Pustejovsky and Tipton (2018). For a coefficient \eqn{c'\beta}, with
#' \eqn{t_g = (I - H_{gg})^{-1/2} W_g^{1/2} Z_g V_b c}, the variance estimate
#' is \eqn{\sum_g (t_g' \hat r_g)^2} in the standardized scores
#' \eqn{\hat r = W^{-1/2}\hat s}. If the true standardized scores \eqn{r} had
#' identity covariance, then \eqn{\hat r = (I - H) r} and the estimate is a
#' quadratic form of \eqn{r} whose degrees of freedom are
#' \eqn{(\sum_g q_g'q_g)^2 / \sum_{g,h} (q_g'q_h)^2}, where \eqn{q_g} is the
#' vector \eqn{(I - H)\, t_g} with \eqn{t_g} placed in the rows of cluster
#' \eqn{g}. `reference = "normal"` keeps the normal distribution. The other
#' sandwich types always use the normal distribution.
#'
#' @param cluster `NULL` (the default) or a one-sided formula with a single
#'   variable of the data, such as `~ id`.
#' @param data `NULL` (the default) or the data frame that holds the cluster
#'   variable, with one row for each observation of the fit.
#' @param reference For `hc2_sandwich()` with a `cluster`: `"t"` (the default)
#'   or `"normal"`, the reference distribution of the Wald intervals and
#'   tests.
#' @param weights `"frequency"` (the default) or `"sampling"`: whether the
#'   prior weights of the fit are frequencies or sampling weights, as at
#'   [hc0_sandwich()].
#'
#' @return A `SandwichMethod` object with the properties `type` (`"HC2"`,
#'   `"HC3"` or `"HC4"`), `cluster`, `data`, `reference` and `weights`.
#'
#' @references
#' MacKinnon, J. G. and White, H. (1985). Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' *Journal of Econometrics*, 29(3), 305--325.
#'
#' Bell, R. M. and McCaffrey, D. F. (2002). Bias reduction in standard errors
#' for linear regression with multi-stage samples. *Survey Methodology*,
#' 28(2), 169--181.
#'
#' Cribari-Neto, F. (2004). Asymptotic inference under heteroskedasticity of
#' unknown form. *Computational Statistics & Data Analysis*, 45(2), 215--233.
#'
#' Pustejovsky, J. E. and Tipton, E. (2018). Small-sample methods for
#' cluster-robust variance estimation and hypothesis testing in fixed effects
#' models. *Journal of Business & Economic Statistics*, 36(4), 672--683.
#'
#' @seealso [hc0_sandwich()], [statmod_leverage()], [frequentist_vcov()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150), id = rep(1:30, each = 5))
#' dd$y <- 1 + 2 * dd$x + rnorm(150, sd = 0.2 + dd$x)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' sqrt(diag(vcov(fit, vcov_type = frequentist_vcov(
#'   sandwich = hc3_sandwich()), param = "mu")))
#' vcov(fit, vcov_type = frequentist_vcov(
#'   sandwich = hc2_sandwich(cluster = ~ id)), param = "mu")
#' @name hc2_sandwich
#' @export
hc2_sandwich <- function(cluster = NULL, data = NULL,
                         reference = c("t", "normal"),
                         weights = c("frequency", "sampling")) {
  SandwichMethod(type = "HC2", cluster = cluster, data = data,
                 reference = match.arg(reference),
                 weights = match.arg(weights))
}

#' @rdname hc2_sandwich
#' @export
hc3_sandwich <- function(cluster = NULL, data = NULL,
                         weights = c("frequency", "sampling")) {
  SandwichMethod(type = "HC3", cluster = cluster, data = data,
                 weights = match.arg(weights))
}

#' @rdname hc2_sandwich
#' @export
hc4_sandwich <- function(weights = c("frequency", "sampling")) {
  SandwichMethod(type = "HC4", weights = match.arg(weights))
}


#' The Exponent of the Leverage Correction
#'
#' @param type `"HC2"`, `"HC3"` or `"HC4"`.
#'
#' @return `-1/2` for HC2, `-1` for HC3, `NA` for HC4, whose exponent depends
#'   on the observation.
#'
#' @keywords internal
sandwich_power <- function(type) {
  switch(type, HC2 = -0.5, HC3 = -1, HC4 = NA_real_)
}


#' Apply a Function of a Symmetric Block Through Its Eigenvalues
#'
#' @description
#' Computes \eqn{(I - H)^{p}} for a symmetric \eqn{H}, with the factor zero on
#' the directions where an eigenvalue of \eqn{H} reaches one.
#'
#' @param H A symmetric matrix.
#' @param power The exponent, negative.
#'
#' @return A list with the matrix `F` and `capped`, the number of eigenvalues
#'   treated as one.
#'
#' @keywords internal
leverage_power <- function(H, power) {
  e <- eigen((H + t(H)) / 2, symmetric = TRUE)
  lam <- e$values
  cap <- lam >= 1 - 1e-8
  fv <- ifelse(cap, 0, (1 - lam)^power)
  list(F = e$vectors %*% (fv * t(e$vectors)), capped = sum(cap))
}


#' Rescale the Scores of a Fit by the Leverage
#'
#' @description
#' Computes the scores \eqn{s_i^*} of [hc2_sandwich()] (or of the CR2, CR3 and
#' HC4 variants) and, for CR2, the closure that evaluates the Satterthwaite
#' degrees of freedom.
#'
#' @param fit The [StatmodFit()].
#' @param design The design of the fit.
#' @param Vb The inverse penalized information over the coordinates `cols`.
#' @param cols The positions, among the coefficients stacked by equation, of
#'   the coordinates that `Vb` covers.
#' @param g The list of weighted scores of [statmod_score_obs()].
#' @param id `NULL` or the cluster numbers.
#' @param sw The [SandwichMethod] object.
#' @param expected,emethod The information that the bread was built from.
#' @param lev_info The information of the leverage: the roots \eqn{W_i^{1/2}}
#'   and the inverse \eqn{K^{-1}} of the hat matrix are read on it. `NULL`
#'   takes the information of the bread.
#'
#' @return A list with `g`, the list `g` with the rescaled scores in the place
#'   of the originals, and `df_fn`, `NULL` unless the type is CR2 with a
#'   `t` reference, otherwise a function of the positions (among `cols`) that
#'   returns the degrees of freedom there.
#'
#' @keywords internal
sandwich_adjusted <- function(fit, design, Vb, cols, g, id, sw, expected,
                              emethod, lev_info = NULL) {
  spec <- fit@spec
  w <- spec@weights
  if (is.null(lev_info)) {
    lev_info <- if (expected) distributions7::expected_info(method = emethod)
    else distributions7::observed_info()
  }
  lev_exp <- info_expected(lev_info)
  lev_meth <- info_emethod(lev_info)
  parts <- leverage_parts(fit, design, cols, lev_exp, lev_meth)
  # the hat matrix is the projection of the leverage's information: when that
  # differs from the bread's, K^-1 is rebuilt on it
  Vl <- if (identical(lev_exp, expected) && identical(lev_meth, emethod)) Vb else
    leverage_inverse(fit, design, cols, lev_exp, lev_meth, colnames(Vb))
  n <- parts$n
  d <- parts$d
  S <- vapply(parts$params, function(p) g[[p]], numeric(n))
  S <- matrix(S, n, d)
  Rt <- parts$root
  Rp <- parts$root_pinv
  capped <- 0L
  Sst <- S
  df_fn <- NULL

  if (is.null(id)) {
    Hb <- leverage_blocks(parts, Vl)
    pos <- w > 0
    Hc <- Hb
    sampling <- identical(sw@weights, "sampling")
    if (!sampling) for (a in seq_len(d)) for (b in seq_len(d)) {
      Hc[, a, b] <- ifelse(pos, Hb[, a, b] / ifelse(pos, w, 1), 0)
    }
    tr <- vapply(seq_len(d), function(a) Hc[, a, a], numeric(n))
    tr <- rowSums(matrix(tr, n, d))
    pw <- sandwich_power(sw@type)
    if (is.na(pw)) {
      hbar <- if (sampling) mean(tr[pos]) else sum(w * tr) / sum(w)
      delta <- pmin(4, tr / hbar)
      pw_i <- -delta / 2
    } else pw_i <- rep(pw, n)
    if (d == 1L) {
      h <- Hc[, 1L, 1L]
      cap <- h >= 1 - 1e-8
      fac <- ifelse(cap, 0, (1 - pmin(h, 1 - 1e-8))^pw_i)
      fac[!pos] <- 1
      capped <- sum(cap & pos)
      Sst[, 1L] <- S[, 1L] * fac
    } else {
      for (i in which(pos)) {
        fp <- leverage_power(matrix(Hc[i, , ], d, d), pw_i[[i]])
        capped <- capped + fp$capped
        Wh <- matrix(Rt[i, , ], d, d)
        Whp <- matrix(Rp[i, , ], d, d)
        Sst[i, ] <- S[i, ] + drop(Wh %*% ((fp$F - diag(d)) %*% (Whp %*% S[i, ])))
      }
    }
  } else {
    pw <- sandwich_power(sw@type)
    groups <- split(seq_len(n), id)
    blk <- vector("list", length(groups))
    for (k in seq_along(groups)) {
      idx <- groups[[k]]
      zr <- as.vector(t(outer((idx - 1L) * d, seq_len(d), "+")))
      Zg <- parts$Z[zr, , drop = FALSE]
      m <- length(idx)
      Whg <- Whgp <- matrix(0, m * d, m * d)
      for (j in seq_len(m)) {
        rr <- (j - 1L) * d + seq_len(d)
        Whg[rr, rr] <- matrix(Rt[idx[[j]], , ], d, d)
        Whgp[rr, rr] <- matrix(Rp[idx[[j]], , ], d, d)
      }
      ZV <- as.matrix(Zg %*% Vl)
      Mg <- ZV %*% t(as.matrix(Zg))
      Hg <- Whg %*% Mg %*% Whg
      fp <- leverage_power(Hg, pw)
      capped <- capped + fp$capped
      sg <- as.vector(t(S[idx, , drop = FALSE]))
      sstar <- sg + drop(Whg %*% ((fp$F - diag(m * d)) %*% (Whgp %*% sg)))
      Sst[idx, ] <- matrix(sstar, m, d, byrow = TRUE)
      blk[[k]] <- list(zr = zr, Zg = Zg, Whg = Whg, F = fp$F)
    }
    if (identical(sw@type, "HC2") && identical(sw@reference, "t")) {
      df_fn <- cr2_df_closure(parts, blk, Vl)
    }
  }
  if (capped > 0L) {
    warning(sprintf(paste0("The leverage of %d direction%s of the hat matrix ",
                           "reaches one; the correction uses the ",
                           "pseudo-inverse there."),
                    capped, if (capped == 1L) "" else "s"), call. = FALSE)
  }
  for (a in seq_len(d)) g[[parts$params[[a]]]] <- Sst[, a]
  list(g = g, df_fn = df_fn)
}


#' Satterthwaite Degrees of Freedom of the CR2 Variance
#'
#' @description
#' Returns a function that computes the Satterthwaite degrees of freedom of
#' Pustejovsky and Tipton (2018) for the coefficients of the model, as
#' described at [hc2_sandwich()].
#'
#' @param parts The list of [leverage_parts()].
#' @param blk The list of clusters built by [sandwich_adjusted()], each with
#'   the rows `zr`, the design `Zg`, the roots `Whg` and the matrix
#'   \eqn{(I - H_{gg})^{-1/2}} as `F`.
#' @param Vb The inverse penalized information.
#'
#' @return A function of `k`, positions among the coordinates of `Vb`, and
#'   `contrast`, that returns a numeric vector of degrees of freedom. With
#'   `contrast` (a matrix with one row for each coordinate of `Vb` and one
#'   column for each quantity), the degrees of freedom are those of the
#'   contrasts \eqn{c'\beta} given by the columns, and `k` is ignored.
#'
#' @keywords internal
cr2_df_closure <- function(parts, blk, Vb) {
  force(parts); force(blk); force(Vb)
  Jm <- NULL
  function(k = seq_len(ncol(Vb)), chunk = 40L, contrast = NULL) {
    if (!is.null(contrast)) {
      contrast <- as.matrix(contrast)
      VCall <- Vb %*% contrast
      k <- seq_len(ncol(contrast))
    }
    if (is.null(Jm)) {
      n <- parts$n; d <- parts$d
      ij <- expand.grid(a = seq_len(d), b = seq_len(d))
      oi <- rep(seq_len(n), each = nrow(ij))
      Wh <- Matrix::sparseMatrix(
        i = (oi - 1L) * d + rep(ij$a, n), j = (oi - 1L) * d + rep(ij$b, n),
        x = as.vector(t(vapply(seq_len(nrow(ij)), function(r)
          parts$root[, ij$a[[r]], ij$b[[r]]], numeric(n)))),
        dims = c(n * d, n * d))
      Y <- Wh %*% parts$Z
      Jm <<- as.matrix(Matrix::crossprod(Y))
    }
    G <- length(blk)
    out <- rep(NA_real_, length(k))
    if (!length(k)) return(out)
    for (s0 in seq(1L, length(k), by = chunk)) {
      sel <- s0:min(length(k), s0 + chunk - 1L)
      VC <- if (is.null(contrast)) Vb[, k[sel], drop = FALSE] else
        VCall[, sel, drop = FALSE]
      A <- array(0, c(nrow(Vb), G, length(sel)))
      tau <- matrix(0, G, length(sel))
      for (gi in seq_len(G)) {
        b <- blk[[gi]]
        Tg <- b$F %*% (b$Whg %*% as.matrix(b$Zg %*% VC))
        A[, gi, ] <- as.matrix(Matrix::crossprod(b$Zg, b$Whg %*% Tg))
        tau[gi, ] <- colSums(Tg^2)
      }
      for (j in seq_along(sel)) {
        Ak <- A[, , j]
        if (G == 1L) Ak <- matrix(Ak, ncol = 1L)
        Bk <- Vb %*% Ak
        Gm <- diag(tau[, j], G) - 2 * crossprod(Ak, Bk) +
          crossprod(Bk, Jm %*% Bk)
        tr <- sum(diag(Gm))
        den <- sum(Gm^2)
        out[sel[[j]]] <- if (is.finite(den) && den > 0) tr^2 / den else NA_real_
      }
    }
    out
  }
}


#' The Reference Quantile of a Wald Interval
#'
#' @description
#' Returns the critical value for each coefficient of a table: the normal
#' quantile, or the \eqn{t} quantile where the variance carries Satterthwaite
#' degrees of freedom.
#'
#' @param V A variance matrix, possibly with the attribute `"df_fn"` that
#'   [vcov.StatmodFit()] copies from the sandwich.
#' @param nm The names of the coefficients, rows of `V`.
#' @param se Their standard errors; a coordinate without a finite one gets no
#'   degrees of freedom computed.
#' @param level The confidence level.
#'
#' @return A list with `z` (a vector of critical values) and `df` (the
#'   degrees of freedom, `Inf` where the reference is the normal).
#'
#' @keywords internal
wald_reference <- function(V, nm, se, level) {
  fn <- attr(V, "df_fn")
  df <- rep(Inf, length(nm))
  if (!is.null(fn)) {
    use <- which(is.finite(se))
    if (length(use)) {
      got <- fn(nm[use])
      df[use[!is.na(got)]] <- got[!is.na(got)]
    }
  }
  list(z = stats::qt((1 + level) / 2, df), df = df)
}


#' The Two-Sided P-Value of a Wald Statistic
#'
#' @param stat The statistics \eqn{z}.
#' @param df `NULL` for the normal distribution, otherwise the degrees of
#'   freedom, one for each statistic (`Inf` for the normal).
#'
#' @return A numeric vector.
#'
#' @keywords internal
wald_pvalue <- function(stat, df = NULL) {
  if (is.null(df)) return(2 * stats::pnorm(-abs(stat)))
  2 * stats::pt(-abs(stat), df)
}


#' Degrees of Freedom Addressed by Coordinate Name
#'
#' @param fn The function of positions that [cr2_df_closure()] returns.
#' @param names The names of the coordinates that the sandwich covers, in the
#'   order of its columns.
#'
#' @return A function of a character vector of names that returns the degrees
#'   of freedom, `NA` for a name the variance does not cover. Its second
#'   argument `contrast` is a matrix whose rows are named by coordinates and
#'   whose columns are contrasts; the degrees of freedom of those contrasts
#'   are returned, `NA` for a contrast that reads a coordinate the variance
#'   does not cover.
#'
#' @keywords internal
sandwich_df_by_name <- function(fn, names) {
  force(fn); force(names)
  structure(function(nm = NULL, contrast = NULL) {
    if (!is.null(contrast)) {
      pos <- match(rownames(contrast), names)
      miss <- is.na(pos)
      bad <- if (any(miss)) colSums(contrast[miss, , drop = FALSE] != 0) > 0 else
        rep(FALSE, ncol(contrast))
      Cf <- matrix(0, length(names), ncol(contrast))
      Cf[pos[!miss], ] <- contrast[!miss, , drop = FALSE]
      res <- rep(NA_real_, ncol(contrast))
      if (any(!bad)) res[!bad] <- fn(contrast = Cf[, !bad, drop = FALSE])
      return(res)
    }
    pos <- match(nm, names)
    res <- rep(NA_real_, length(nm))
    ok <- !is.na(pos)
    if (any(ok)) res[ok] <- fn(pos[ok])
    res
  }, class = c("statmod_df_fn", "function"))
}


#' Degrees of Freedom of One Coefficient Under a CR2 Variance
#'
#' @param fit The [StatmodFit()].
#' @param vcov_type The [vcov_method] object.
#' @param name The name of the coefficient, `"parameter:coefficient"`.
#'
#' @return The Satterthwaite degrees of freedom when `vcov_type` is a
#'   frequentist variance with [hc2_sandwich()] clustered and
#'   `reference = "t"`, and `Inf` otherwise.
#'
#' @keywords internal
sandwich_wald_df <- function(fit, vcov_type, name) {
  sw <- if (S7::prop_exists(vcov_type, "sandwich")) vcov_type@sandwich else NULL
  if (is.null(sw) || !identical(sw@type, "HC2") || is.null(sw@cluster) ||
      !identical(sw@reference, "t")) {
    return(Inf)
  }
  V <- stats::vcov(fit, vcov_type = vcov_type, readable = FALSE)
  fn <- attr(V, "df_fn")
  if (is.null(fn)) return(Inf)
  r <- fn(name)
  if (is.na(r)) Inf else r
}


#' @param x The function of coordinate names.
#' @param ... Unused.
#' @rdname sandwich_df_by_name
#' @export
print.statmod_df_fn <- function(x, ...) {
  cat("<Satterthwaite degrees of freedom, as a function of coefficient names>
")
  invisible(x)
}
