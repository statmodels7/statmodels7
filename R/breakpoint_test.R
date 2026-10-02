# The test that a break-point exists, after Davies (1987): under the null the
# position is not defined, so the statistic is a process in the position and
# the test reads its supremum.


#' Test That a Break-Point Exists
#'
#' @description
#' Tests the hypothesis that a model has no break-point against the
#' alternative of the break-point term it carries: no change of slope for a
#' [modelterms7::seg()], no change of level for a [modelterms7::jump()], and
#' neither for a [modelterms7::jseg()]. The term can sit in the equation of
#' any distribution parameter, beside any other term, and inside a parameter
#' of [modelterms7::nl()].
#'
#' @details
#' Under the null hypothesis the position \eqn{\psi} of the break-point does
#' not enter the model, so the classical tests of the change against zero
#' have no fixed reference distribution: the position is a nuisance parameter
#' present only under the alternative (Davies, 1987). The test reads Rao's
#' score statistic for the change at a fixed position \eqn{q},
#' \deqn{S(q) = U_c(q)' \, [K(q)^{-1}]_{cc} \, U_c(q),}
#' where \eqn{U_c(q)} is the score of the change coefficients at the fit of
#' the null model and \eqn{K(q)} the penalized information of the model with
#' the break-point held at \eqn{q}. One fit of the null model serves every
#' \eqn{q}: the process costs one gradient and one Hessian per position. The
#' statistic of the test is \eqn{M = \sup_q S(q)}.
#'
#' Where the position moves the model continuously (a change of slope, or a
#' smoothed step), \eqn{S(q)} is evaluated at `k` positions equally spaced
#' between the confinement limits of the term, and the p-value is the upper
#' bound of Davies (1987) for a \eqn{\chi^2_s} process,
#' \deqn{P(\sup_q S(q) > M) \le P(\chi^2_s > M) +
#'   V M^{(s-1)/2} e^{-M/2} 2^{-s/2} / \Gamma(s/2),}
#' with \eqn{V = \sum_j |S(q_{j+1})^{1/2} - S(q_j)^{1/2}|} the total
#' variation of the root of the process and \eqn{s} the number of change
#' coefficients. This is the test of `segmented::davies.test()`, read on the
#' score of the model at hand.
#'
#' Where the step is sharp the process is constant between consecutive
#' values of the covariate and jumps at each of them, and the bound, which
#' counts the variation over every interval, is too wide to be useful. The
#' supremum is then taken over every interval inside the confinement limits,
#' and its distribution under the null is estimated by a parametric
#' bootstrap: `n_boot` responses are drawn from the null fit, the null model
#' is refitted to each, and the p-value is the share of replicates whose
#' supremum reaches the observed one, \eqn{(1 + \#\{M^* \ge M\})/(1 +
#' B)}. The draws are conditional on the fitted random effects, if any.
#'
#' The hyperparameters are held at the values of the fit, so the test is
#' conditional on the smoothing the data chose. The term must carry one
#' break-point whose position has no development; for a choice between one
#' and several break-points, compare the fits by an information criterion.
#'
#' @param fit A [StatmodFit()] carrying the break-point term.
#' @param param The distribution parameter whose equation carries the term,
#'   or `NULL` to find it.
#' @param term The label of the term (`"seg"`, `"jump"`, `"jseg"` by
#'   default), or `NULL` to find it. Needed only where the model carries
#'   more than one break-point term.
#' @param k The number of positions at which a continuous process is read.
#' @param n_boot The number of bootstrap replicates for a sharp step.
#' @param seed A seed for the bootstrap, or `NULL`.
#'
#' @return An object of class `"htest"` with the supremum of the score
#'   process as `statistic`, the number of change coefficients as
#'   `parameter`, the `p.value`, the position of the supremum as
#'   `estimate`, and the process itself as `process`, a data frame with
#'   columns `psi` and `score`.
#'
#' @references
#' Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is
#' present only under the alternative. *Biometrika*, 74, 33--43.
#'
#' Muggeo, V. M. R. (2003). Estimating regression models with unknown
#' break-points. *Statistics in Medicine*, 22, 3055--3071.
#'
#' @seealso [statmod_test()] for a test of one coefficient against one value.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(120, 0, 10))
#' dd$y <- 1 + 0.5 * dd$x - 0.8 * pmax(dd$x - 6, 0) + rnorm(120, sd = 0.5)
#' fit <- statmod(y ~ seg(x), distributions7::gaussian1_distrib(), dd)
#' statmod_breakpoint_test(fit)
#'
#' @export
statmod_breakpoint_test <- function(fit, param = NULL, term = NULL, k = 10L,
                                    n_boot = 199L, seed = NULL) {
  if (!S7::S7_inherits(fit, StatmodFit)) {
    stop("'fit' must be a statmod fit.", call. = FALSE)
  }
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 3L) {
    stop("'k' must be a whole number of at least 3.", call. = FALSE)
  }
  n_boot <- as.integer(n_boot)
  if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 1L) {
    stop("'n_boot' must be a positive whole number.", call. = FALSE)
  }
  spec0 <- fit@spec
  if (length(statmod_structural(spec0))) {
    stop("A model carrying a structural term has no break-point test.",
         call. = FALSE)
  }
  loc <- bp_test_locate(spec0, param, term)
  p <- loc$param
  design0 <- statmod_design(spec0)
  nms <- design0[[p]]$coef_names
  pos_name <- paste0(loc$prefix, "psi1")
  chg <- intersect(paste0(loc$prefix, c("gamma1", "delta1")), nms)
  params <- spec0@distrib@params
  npar <- vapply(design0, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  at_chg <- offs[[match(p, params)]] + match(chg, nms)
  at_pos <- offs[[match(p, params)]] + match(pos_name, nms)
  tm <- loc$term
  bp <- tm@blueprint
  sharp <- loc$sharp
  lim <- bp$lim

  # the positions: every interval of a sharp step, k points of a continuous
  # process
  if (sharp) {
    u <- sort(unique(bp$xv))
    qs <- (u[-1L] + u[-length(u)]) / 2
    qs <- qs[qs > lim[1L] & qs < lim[2L]]
  } else {
    qs <- seq(lim[1L], lim[2L], length.out = k + 2L)[-c(1L, k + 2L)]
  }
  if (length(qs) < 2L) {
    stop("The confinement limits leave fewer than two positions to test.",
         call. = FALSE)
  }

  method <- fit@methods$smooth
  cfg <- inner_settings(method)
  hyper <- fit@hyper

  # the model with the break-point at q and its changes held at zero
  spec_at <- function(q, sp = spec0) {
    hc <- sp@held_coef
    hv <- hc[[p]]
    if (is.null(hv)) hv <- numeric(0)
    hv[chg] <- 0
    terms <- sp@terms
    if (sharp) {
      t0 <- terms[[p]][[loc$name]]
      if (isTRUE(t0@blueprint$held)) t0 <- modelterms7::seg_hold(t0, FALSE)
      t1 <- modelterms7::seg_hold(modelterms7::seg_relocate(t0, q))
      terms[[p]][[loc$name]] <- t1
    } else {
      hv[pos_name] <- q
    }
    hc[[p]] <- hv
    S7::set_props(sp, terms = terms, held_coef = hc)
  }

  # the null fit: the changes at zero, at any position
  null_fit <- function(sp, beta = NULL) {
    s0 <- spec_at(qs[[1L]], sp)
    d0 <- statmod_design(s0)
    b0 <- statmod_blocks(s0, d0)
    # from the package's own start: the coefficients of the alternative with
    # the change set to zero can put the mean far from the data, and the
    # dispersion then runs to the edge of its chart
    if (is.null(beta)) {
      beta <- statmod_start(s0, d0, statmod_objective(s0, hyper, d0),
                            start_intercepts())
    }
    beta[at_chg] <- 0
    if (!sharp) beta[at_pos] <- qs[[1L]]
    res <- statmod_alternate(s0, d0, b0, hyper, method, beta, cfg$expected,
                             cfg$approx, cfg$maxit, cfg$tol, verbosity(0))
    list(par = res$par, value = res$value, spec = s0, design = d0,
         converged = isTRUE(res$converged))
  }

  # the score process at the null estimates
  process <- function(sp, par0) {
    vapply(qs, function(q) {
      sq <- spec_at(q, sp)
      dq <- statmod_design(sq)
      obj <- statmod_objective(sq, hyper, dq, cfg$expected, cfg$approx)
      v <- par0
      if (!sharp) v[at_pos] <- q
      g <- -obj$gr(v)
      K <- tryCatch(as_dense(obj$he(v)), error = function(e) NULL)
      if (is.null(K) || !all(is.finite(g))) return(NA_real_)
      held <- held_positions(sq, dq, obj, v)$where
      held <- setdiff(held, at_chg)
      free <- setdiff(seq_along(v), held)
      Kf <- K[free, free, drop = FALSE]
      ci <- match(at_chg, free)
      Vi <- tryCatch(solve(Kf), error = function(e) NULL)
      if (is.null(Vi)) return(NA_real_)
      uc <- g[at_chg]
      as.numeric(crossprod(uc, Vi[ci, ci, drop = FALSE] %*% uc))
    }, numeric(1))
  }

  f0 <- null_fit(spec0)
  S <- process(spec0, f0$par)
  if (!any(is.finite(S))) {
    stop("The score process could not be evaluated at any position.",
         call. = FALSE)
  }
  M <- max(S, na.rm = TRUE)
  s <- length(chg)
  if (!sharp) {
    r <- sqrt(pmax(S[is.finite(S)], 0))
    V <- sum(abs(diff(r)))
    pv <- stats::pchisq(M, s, lower.tail = FALSE) +
      V * M^((s - 1) / 2) * exp(-M / 2) * 2^(-s / 2) / gamma(s / 2)
    pv <- min(1, pv)
    meth <- "Davies' test for a break-point (score process, upper bound)"
  } else {
    if (!is.null(seed)) {
      if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        stats::runif(1)
      }
      old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
      set.seed(seed)
    }
    cf0 <- statmod_objective(f0$spec, hyper, f0$design)$split(f0$par)
    th <- statmod_eta(f0$spec, f0$design, cf0)$theta
    Mb <- rep(NA_real_, n_boot)
    for (b in seq_len(n_boot)) {
      yb <- distributions7::distrib_rng(spec0@distrib, spec0@n_obs, th)
      sb <- S7::set_props(spec0, response = yb)
      fb <- tryCatch(null_fit(sb, f0$par), error = function(e) NULL)
      if (is.null(fb)) next
      Sb <- tryCatch(process(sb, fb$par), error = function(e) NULL)
      if (!is.null(Sb) && any(is.finite(Sb))) Mb[b] <- max(Sb, na.rm = TRUE)
    }
    ok <- is.finite(Mb)
    if (!any(ok)) {
      stop("No bootstrap replicate could be refitted.", call. = FALSE)
    }
    pv <- (1 + sum(Mb[ok] >= M)) / (1 + sum(ok))
    meth <- sprintf("Test for a break-point (score process, parametric bootstrap, %d replicates)",
                    sum(ok))
  }
  kind <- tm@kind
  out <- list(statistic = c("sup S" = M), parameter = c(df = s),
              p.value = pv, estimate = c(psi = qs[which.max(S)]),
              method = meth,
              data.name = sprintf("%s in the equation of %s", tm@label, p),
              alternative = switch(kind,
                seg = "a change of slope",
                jump = "a change of level",
                jseg = "a change of level and of slope"),
              process = data.frame(psi = qs, score = S))
  class(out) <- "htest"
  out
}


#' Where the Break-Point Term of a Fit Is
#'
#' @description
#' Finds the break-point term a test reads: its equation, the prefix of its
#' coefficient names, the built term, and whether its step is sharp. The term
#' is looked for among the terms of every equation and among the sub-terms of
#' the parameters of a [modelterms7::nl()] term.
#'
#' @param spec A [StatmodSpec()].
#' @param param,term As [statmod_breakpoint_test()].
#'
#' @return A list with `param`, `name` (the term's name in its equation),
#'   `prefix`, `term` and `sharp`.
#'
#' @keywords internal
bp_test_locate <- function(spec, param = NULL, term = NULL) {
  found <- list()
  add <- function(p, nm, prefix, tm) {
    found[[length(found) + 1L]] <<- list(param = p, name = nm, prefix = prefix,
                                         term = tm)
  }
  for (p in spec@distrib@params) {
    if (!is.null(param) && !identical(p, param)) next
    for (nm in names(spec@terms[[p]])) {
      tm <- spec@terms[[p]][[nm]]
      if (S7::S7_inherits(tm, modelterms7::SegTerm)) {
        add(p, nm, paste0(tm@label, "."), tm)
      } else if (S7::S7_inherits(tm, modelterms7::NlTerm)) {
        subs <- tm@blueprint$subs
        for (a in names(subs)) for (st in subs[[a]]) {
          if (S7::S7_inherits(st, modelterms7::SegTerm)) {
            add(p, nm, paste0(tm@label, ".", a, ".", st@label, "."), st)
          }
        }
      }
    }
  }
  if (!is.null(term)) {
    found <- Filter(function(f) identical(f$term@label, term), found)
  }
  if (!length(found)) {
    stop("The fit carries no seg(), jump() or jseg() term to test.",
         call. = FALSE)
  }
  if (length(found) > 1L) {
    stop("The fit carries several break-point terms; name one with 'param' ",
         "and 'term'.", call. = FALSE)
  }
  f <- found[[1L]]
  bp <- f$term@blueprint
  if (bp$npsi != 1L) {
    stop("The test reads a term with one break-point (npsi = 1).",
         call. = FALSE)
  }
  if (any(vapply(bp$Z, Negate(is.null), logical(1)))) {
    stop("The test reads a term whose coefficients have no development.",
         call. = FALSE)
  }
  f$sharp <- !identical(f$term@kind, "seg") && is.null(bp$smooth)
  if (f$sharp && !grepl("^[^.]+\\.$", f$prefix)) {
    stop("A sharp step inside nl() is not a model statmod() fits.",
         call. = FALSE)
  }
  f
}
