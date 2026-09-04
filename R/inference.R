#' @include predict.R
#' @importFrom stats vcov confint
NULL

# What a fitted model reports about its own uncertainty.
#
# The coefficients live on the LINK scale by construction -- they are the
# coefficients of a linear predictor, free on the whole line -- so a symmetric
# Wald interval is built where the quantity actually ranges and needs no
# mapping back. That is the one place this differs from
# distributions7::fit_distrib(), which estimates a bounded parameter and has to
# carry its interval through the link.

#' Where Each Stacked Coefficient Comes From
#'
#' @description
#' One row per stacked coefficient, naming the distribution parameter, the term
#' and the coefficient, and saying whether the term carries a penalty and
#' whether that penalty has a kink.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#'
#' @return A data frame with as many rows as there are coefficients.
#'
#' @keywords internal
coef_labels <- function(spec, design) {
  params <- spec@distrib@params
  units <- statmod_penalized(spec, design)
  rows <- list()
  for (p in params) {
    d <- design[[p]]
    if (d$npar == 0L) next
    term <- rep(NA_character_, d$npar)
    pen <- rep(FALSE, d$npar)
    kink <- rep(FALSE, d$npar)
    for (nm in names(d$blocks)) term[d$blocks[[nm]]] <- nm
    for (u in units) {
      # a class is considered under EVERY equation one of its members sits
      # in; `param` on it is a convention for the hyperparameter store
      here <- if (is.null(u$pieces)) identical(u$param, p) else
        any(vapply(u$pieces, function(z) identical(z$param, p), TRUE))
      if (!here) next
      # a structural unit's penalty sits on the term's OWN parameters,
      # which contribute no column: its cols index the term's parameter
      # vector, and writing them here grew pen past the design and
      # recycled the labels into duplicate rows
      if (isTRUE(u$structural)) next
      # a covariance class has no columns of its own: its members do, each in
      # its own equation, so the marks go on the pieces that live in this one
      cs <- if (is.null(u$pieces)) u$cols else
        unlist(lapply(Filter(function(z) identical(z$param, p), u$pieces),
                      function(z) z$cols), use.names = FALSE)
      pen[cs] <- TRUE
      kink[cs] <- penalty_has_kink(u$penalty)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = p, term = term, coefficient = d$coef_names,
      penalized = pen, kinked = kink, stringsAsFactors = FALSE)
  }
  if (!length(rows)) {
    return(data.frame(parameter = character(0), term = character(0),
                      coefficient = character(0), penalized = logical(0),
                      kinked = logical(0), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- paste(out$parameter, out$coefficient, sep = ":")
  out
}


#' Which Information Matrix a Fit Reports
#'
#' @description
#' `TRUE` when [vcov.StatmodFit()] should invert the expected information:
#' when the fit itself inverted it, which only [iwls()] does and only unless
#' asked otherwise, AND the family writes that information out in closed
#' form.
#'
#' @details
#' The default follows the fit rather than choosing for itself, so a standard
#' error comes from the same matrix the step did, and a caller who wants the
#' other one asks for it.
#'
#' The second condition is the part that is not about the fit. Where a family
#' has no closed expected information the scoring step is driven by an
#' approximation of it -- the outer product of the observed scores, which
#' costs one gradient and is positive semidefinite by construction -- and that
#' is a good matrix to take a step with, the score being exact, but not a good
#' one to read a standard error off. Measured on a Poisson-inverse gaussian
#' regression at \eqn{n = 500}: the outer product gives standard errors 5.7
#' per cent from those of the exact expectation where the observed Hessian
#' gives 0.6 per cent, for coefficients agreeing to \eqn{10^{-6}}. So the
#' report falls back to the observed information there, which every family has
#' and which is exact.
#'
#' @param object A [StatmodFit()].
#'
#' @return A single logical.
#'
#' @seealso [distributions7::expected_hessian_exact()], the predicate this
#'   reads, and [vcov.StatmodFit()], which follows it.
#'
#' @keywords internal
fit_expected <- function(object) {
  m <- object@methods$smooth
  # An optimizers7 optimizer is given the observed information by
  # inner_settings(), so it did not invert the expected one and the report
  # must not say it did: the rule above is that a standard error comes from
  # the matrix the step used.
  used <- S7::S7_inherits(m, Iwls) && identical(m@hessian, "expected")
  used && distributions7::expected_hessian_exact(object@spec@distrib)
}


#' @title The Variance Matrix of a Fit
#' @name vcov.StatmodFit
#' @description
#' The variance of the estimated coefficients, over every distribution
#' parameter's block at once.
#' @details
#' **Two matrices, and they differ only when something is penalized.**
#' Writing \eqn{H} for the information of the log-likelihood and \eqn{S} for
#' the second derivative of the penalty,
#' \deqn{V_b = (H + S)^{-1}, \qquad V_f = (H+S)^{-1} H (H+S)^{-1}.}
#' The first is the posterior variance under the prior the penalty is the
#' negative logarithm of, and it is what an interval around a penalized term
#' should be built from: it carries the smoothing bias as though it were
#' variance, and that is why such intervals cover at about their nominal
#' rate. The second is the sampling variance of the penalized estimator at a
#' fixed penalty, which is smaller and covers less. With no penalty \eqn{S = 0}
#' and both are \eqn{H^{-1}}.
#'
#' **A coefficient a kinked penalty has set to zero has no row.** At zero
#' the penalty is not twice differentiable, so \eqn{S} does not exist there and
#' no curvature can be read; the entry is `NA`. The coefficients a lasso
#' or an MCP left non-zero do get a variance, and it is conditional on that
#' selection, which [summary.StatmodFit()] says in a note instead of leaving
#' a reader to assume otherwise.
#' **Which information.** `expected` says which matrix \eqn{H} is. Its
#' default is the expected information where the fit inverted it AND the
#' family writes it out in closed form, and the observed Hessian otherwise
#' ([fit_expected()]). The two agree asymptotically and not in a sample: one
#' is the information averaged over the model, the other the curvature of the
#' likelihood at the data in hand.
#'
#' Where the family has no closed form, `expected = TRUE` reaches an
#' approximation and `approx` says which. `"opg"`, the default, is the outer
#' product of the observed scores and costs one gradient; `"bartlett"`
#' evaluates the expectation itself, a sum over the support for a discrete
#' family and a quadrature for a continuous one, and is orders of magnitude
#' dearer -- 89.06 s against 0.64 s on a Poisson-inverse gaussian regression
#' at \eqn{n = 500}. The expensive route is reachable and is not the default.
#'
#' @param object A [StatmodFit()].
#' @param type `"bayesian"` or `"frequentist"`.
#' @param expected Whether the expected information is used. Defaults to
#'   [fit_expected()]: the expected one where the fit inverted it and the
#'   family writes it out, the observed Hessian otherwise.
#' @param approx How the expected information is approximated for a family
#'   with no closed form: `"opg"` (the default), `"bartlett"`, `"integrate"`
#'   or `"mc"`. Read only when `expected` is `TRUE` and the family has no
#'   closed form.
#' @param ... Unused.
#' @return A square matrix over the stacked coefficients, with dimnames
#'   `parameter:coefficient`.
#' @seealso [confint.StatmodFit()], [summary.StatmodFit()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' sqrt(diag(vcov(fit)))
#' @keywords internal
vcov.StatmodFit <- function(object, type = c("bayesian", "frequentist"),
                            expected = NULL,
                            approx = c("opg", "bartlett", "integrate", "mc"),
                            readable = TRUE,
                            parameter = NULL, ...) {
  type <- match.arg(type)
  approx <- match.arg(approx)
  if (is.null(expected)) expected <- fit_expected(object)
  spec <- object@spec
  design <- statmod_design(spec)
  coef <- object@coefficients
  lab <- coef_labels(spec, design)
  total <- nrow(lab)
  nm <- rownames(lab)

  # A model carrying a structural term is inverted over the coefficients AND
  # the term's own parameters together, and the coefficient block of that
  # inverse is taken. Inverting the coefficient block alone would report the
  # variance that would hold with the term's parameters known, which is not
  # what was estimated. The information there is the observed one: neither a
  # filter nor a mixture over states has an expected information to offer.
  fil <- length(attr(design, "structural")) > 0L
  H <- if (fil) statmod_full_information(spec, coef, design) else
    statmod_information_at(spec, coef, design, expected, approx)
  S <- statmod_penalty_at(spec, coef, object@hyper, design, "hessian")
  nz <- nrow(H) - total
  if (nz > 0L) {
    # A structural term's own information is assembled dense, so the padded
    # penalty is dense too rather than a sum of two kinds whose result would
    # depend on which side of the bind it arrived on.
    S <- as_dense(S)
    S <- rbind(cbind(S, matrix(0, total, nz)),
               matrix(0, nz, total + nz))
    # A penalty over the term's OWN parameters belongs in that tail. Leaving
    # it out is not a conservative choice: the deviations of a panel are not
    # identified without it -- a constant added to a population value and
    # taken off every deviation leaves the filter unchanged -- so the matrix
    # is singular along exactly that direction and nothing can be reported.
    ps <- structural_penalty_block(spec, design, object@hyper, nz)
    if (!is.null(ps)) S[total + seq_len(nz), total + seq_len(nz)] <- ps
  }

  keep <- rep(TRUE, total)
  beta <- unlist(coef[spec@distrib@params], use.names = FALSE)
  keep[lab$kinked & beta == 0] <- FALSE
  frz <- frozen_block(spec, lab)
  keep[frz] <- FALSE
  # a kinked penalty contributes no curvature away from its kink either, and
  # any non-finite entry would be the kink itself reached by a hair
  S <- zap_nonfinite(S)

  out <- matrix(NA_real_, total, total, dimnames = list(nm, nm))
  if (!any(keep)) return(out)
  keep_full <- c(keep, rep(TRUE, nz))
  A <- (H + S)[keep_full, keep_full, drop = FALSE]
  lb <- c(nm[keep], rep("", nz))
  # what a SMALL eigenvalue means is decided inside solve_pd(), on the
  # equilibrated matrix: a smoothing parameter at 1e15 and a break-point
  # term's annealed working columns both separate the scales without
  # flattening any direction, and per-direction scaling forgives both
  Vb <- tryCatch(solve_pd(A, "the penalized information", lb),
                 error = function(e) e)
  # A COORDINATE THE INFORMATION CARRIES NOTHING ABOUT IS HELD, and the rest
  # is reported. A combination c'beta has a finite variance exactly where c
  # is orthogonal to the null space, so a coefficient the null space does
  # not touch keeps its variance whatever happens to the others: refusing
  # the whole matrix loses those, which on a count model whose dispersion
  # ran to its Poisson limit in one province was 56 standard errors thrown
  # away for one that does not exist. It runs ONLY where solve_pd() has
  # already refused, so a fit whose matrix is invertible is untouched, and
  # where nothing can be held the original message stands.
  if (inherits(Vb, "error")) {
    flat <- uninformative_coords(A)
    if (!length(flat)) stop(conditionMessage(Vb), call. = FALSE)
    warning(held_condition(sprintf(paste0(
      "The penalized information carries nothing about %s, so",
      " %s no
  variance and %s row%s missing here. The rest of the matrix",
      " stands:
  a coefficient is estimable exactly where it is orthogonal",
      " to the flat
  direction, and these coordinates carry the whole of",
      " it."),
      paste(ifelse(nzchar(lb[flat]), lb[flat], "a structural parameter"),
            collapse = ", "),
      if (length(flat) == 1L) "it has" else "they have",
      if (length(flat) == 1L) "its" else "their",
      if (length(flat) == 1L) " is" else "s are")))
    keep_full[which(keep_full)[flat]] <- FALSE
    A <- (H + S)[keep_full, keep_full, drop = FALSE]
    Vb <- solve_pd(A, "the penalized information", lb[-flat])
  }
  V <- if (type == "bayesian") Vb else
    Vb %*% H[keep_full, keep_full, drop = FALSE] %*% Vb
  # THE JOINT INVERSE IS KEPT WHOLE. Its coefficient block is not the
  # inverse of the coefficient block wherever the two are correlated, which
  # is why the matrix is built jointly; and the structural block is the
  # variance of the term's own parameters, which used to be computed here
  # and dropped, leaving a model whose predictor is a filter with no
  # variance to report at all.
  if (any(frz)) {
    warning(frozen_condition(sprintf(paste0(
      "The block of %s is a working linearization with",
                           " a frozen weight,\n  not a Jacobian, so its",
                           " curvature is not the model's: every one of its",
                           "\n  coefficients would be reported with a",
                           " standard error of zero. They are\n  missing",
                           " here instead. The rest of the matrix stands and",
                           " is conditional\n  on the break-points where",
                           " they are; resample for their own",
      " uncertainty."),
      paste(unique(lab$term[frz]), collapse = ", "))))
  }
  out <- matrix(NA_real_, total + nz, total + nz)
  out[keep_full, keep_full] <- V
  rn <- c(nm, if (nz) rep("", nz) else character(0))
  who <- c(lab$parameter, if (nz) rep(NA_character_, nz) else character(0))
  if (nz) {
    tl <- structural_tail_names(spec, design)
    if (length(tl) == nz) {
      rn[total + seq_len(nz)] <- unname(tl)
      who[total + seq_len(nz)] <- names(tl)
    }
  }
  dimnames(out) <- list(rn, rn)
  if (isTRUE(readable)) {
    out <- readable_vcov(spec, design, object, out)
  } else if (nz) {
    # BY EQUATION, as the readable route and `coef()` both are. The joint
    # order puts every design block first and the structural tail last, so a
    # model whose predictor is a filter listed the scale's intercept before
    # the mean's own parameters and did not line up with `coef()`.
    ord <- order(match(who, spec@distrib@params))
    out <- out[ord, ord, drop = FALSE]
  }
  if (!is.null(parameter)) out <- select_parameter(out, parameter, spec)
  out
}

#' The Names of the Structural Tail of the Joint Information
#'
#' @description
#' The free parameters of every structural term, under the names a
#' coefficient of the same term would carry.
#'
#' @details
#' A level an intercept in the same equation carries is held and is not in
#' the joint information, so it is not here either; the tail is the free ones
#' in the term's own order, which is the order the information was assembled
#' in.
#'
#' @param spec The fitted specification.
#' @param design The design.
#'
#' @return A character vector, empty where there is no structural term.
#'
#' @keywords internal
structural_tail_names <- function(spec, design) {
  su <- attr(design, "structural")
  if (!length(su)) return(character(0))
  sst <- statmod_structural_state(design)
  out <- character(0)
  for (u in su) {
    tm <- spec@terms[[u$param]][[u$term]]
    free <- setdiff(modelterms7::term_params(tm), sst$held[[u$term]])
    lb <- tryCatch(tm@label, error = function(e) "")
    if (length(lb) == 1L && nzchar(lb)) free <- paste(lb, free, sep = ".")
    out <- c(out, stats::setNames(paste(u$param, free, sep = ":"),
                                  rep(u$param, length(free))))
  }
  out
}

#' The Variance of the Quantities a Fit Reports
#'
#' @description
#' The delta method over the joint variance, with the Jacobian
#' [readable_joint()] supplies.
#'
#' @details
#' A quantity that reads a coordinate whose variance is missing has no
#' variance either, and its row and column are left missing. Such a
#' coordinate is one a kinked penalty set to zero, or a parameter held by an
#' intercept. Computing the entry from the coordinates that do have a
#' variance would report a number for something the fit did not estimate.
#'
#' @param spec The fitted specification.
#' @param design The design.
#' @param fit The fit.
#' @param V The joint variance over the coordinates.
#'
#' @return The variance over the quantities.
#'
#' @keywords internal
readable_vcov <- function(spec, design, fit, V) {
  rj <- tryCatch(readable_joint(spec, design, fit), error = function(e) NULL)
  if (is.null(rj) || nrow(rj$jacobian) == 0L ||
      ncol(rj$jacobian) != nrow(V)) {
    return(V)
  }
  J <- rj$jacobian
  ok <- !is.na(diag(V))
  bad <- rj$held | apply(J[, !ok, drop = FALSE] != 0, 1L, any)
  Vz <- V
  Vz[is.na(Vz)] <- 0
  out <- J %*% Vz %*% t(J)
  # SYMMETRIZED, and not symmetric by construction. An entry and its
  # transpose collect the same products in a different ORDER, and
  # floating-point addition is not associative, so the two can differ
  # in the last bit -- measured, they did on Linux and not on Windows
  # or macOS, which is what a bit-for-bit assertion turns into a red
  # job on one platform of five. A variance matrix is symmetric, so
  # the halving is the contract rather than a repair.
  out <- (out + t(out)) / 2
  out[bad, ] <- NA_real_
  out[, bad] <- NA_real_
  nm <- paste(rj$parameter, rj$name, sep = ":")
  dimnames(out) <- list(nm, nm)
  out
}

#' One Equation's Submatrix of a Variance Matrix
#'
#' @description
#' The rows and columns of the named distribution parameters.
#'
#' @param V A variance matrix whose names are `parameter:name`.
#' @param parameter The distribution parameters to keep.
#' @param spec The fitted specification, for the message.
#'
#' @return The submatrix.
#'
#' @keywords internal
select_parameter <- function(V, parameter, spec) {
  who <- sub(":.*$", "", rownames(V))
  keep <- who %in% parameter
  if (!any(keep)) {
    stop(sprintf(paste0("'parameter' matched nothing. The model's ",
                        "distribution parameters are %s."),
                 paste(spec@distrib@params, collapse = ", ")), call. = FALSE)
  }
  V[keep, keep, drop = FALSE]
}
S7::method(vcov, StatmodFit) <- vcov.StatmodFit

#' The Quantities a Fit Reports, and the Map From Its Coordinates
#'
#' @description
#' One row per quantity the model is about, over the coefficients of every
#' equation and the own parameters of a structural term, with the Jacobian
#' from the coordinates those were estimated on.
#'
#' @details
#' This is the one place the readable view is built, so that
#' [coef.StatmodFit()], [vcov.StatmodFit()] and
#' [confint.StatmodFit()] cannot report a quantity under one name
#' and index it under another.
#'
#' A term says what it is about through
#' [modelterms7::term_readable()], which gives the quantities and
#' the Jacobian. The coordinates that Jacobian touches are replaced by the
#' quantities; a coordinate no quantity reads stands where it is, with a unit
#' row of its own. That is what leaves a developed parameter intact: its
#' development is a vector of coefficients over covariates with no single
#' value, so nothing is declared for it and nothing is taken away.
#'
#' The joint coordinate vector is the design coefficients of every equation
#' in order, then the **free** parameters of each structural term. Free,
#' because a level an intercept in the same equation carries is held and is
#' absent from the information the variance comes from. A quantity that reads a held
#' parameter is marked: its value stands, and its variance would be that of
#' the rest alone, so it is not reported.
#'
#' @param spec The fitted specification.
#' @param design The design.
#' @param fit The fit.
#'
#' @return A list with `name`, `value`, `parameter`,
#'   `term`, `scale` (one link per quantity), `held` (whether
#'   the quantity reads a parameter that is not estimated), `jacobian`
#'   (quantities by joint coordinates) and `n_design`.
#'
#' @seealso [modelterms7::term_readable()],
#'   [statmod_structural_table()]
#'
#' @keywords internal
readable_joint <- function(spec, design, fit) {
  params <- spec@distrib@params
  npar <- vapply(params, function(p) design[[p]]$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  su <- attr(design, "structural")
  sst <- if (length(su)) statmod_structural_state(design) else NULL

  tails <- list()
  at <- total
  for (u in su) {
    tm <- spec@terms[[u$param]][[u$term]]
    nmp <- modelterms7::term_params(tm)
    hl <- if (is.null(sst)) character(0) else sst$held[[u$term]]
    free <- setdiff(nmp, hl)
    tails[[length(tails) + 1L]] <- list(
      u = u, tm = tm, nm = nmp, held = hl, free = free, at = at)
    at <- at + length(free)
  }
  njoint <- at
  ident <- linkfunctions7::identity_link()

  rows <- list()
  add <- function(name, value, parameter, term, scale, held, j) {
    rows[[length(rows) + 1L]] <<- list(
      name = name, value = value, parameter = parameter, term = term,
      scale = scale, held = held, j = j)
  }

  for (pi in seq_along(params)) {
    p <- params[[pi]]
    d <- design[[p]]
    if (!d$npar) next
    cf <- fit@coefficients[[p]]
    term_of <- rep(NA_character_, d$npar)
    for (nm in names(d$blocks)) term_of[d$blocks[[nm]]] <- nm
    # which coordinates a term's quantities cover, and what they are
    swap <- vector("list", d$npar)
    drop <- rep(FALSE, d$npar)
    for (nm in names(d$blocks)) {
      idx <- d$blocks[[nm]]
      if (!length(idx)) next
      term <- spec@terms[[p]][[nm]]
      rd <- tryCatch(modelterms7::term_readable(term, cf[idx]),
                     error = function(e) NULL)
      if (is.null(rd) || !length(rd$name)) next
      sup <- which(apply(rd$jacobian != 0, 2L, any))
      if (!length(sup)) next
      lb <- tryCatch(term@label, error = function(e) "")
      drop[idx[sup]] <- TRUE
      swap[[idx[sup][1L]]] <- list(rd = rd, idx = idx, term = nm, lb = lb)
    }
    for (i in seq_len(d$npar)) {
      s <- swap[[i]]
      if (!is.null(s)) {
        for (k in seq_along(s$rd$name)) {
          j <- numeric(njoint)
          j[offs[[pi]] + s$idx] <- s$rd$jacobian[k, ]
          nmk <- if (nzchar(s$lb)) paste(s$lb, s$rd$name[[k]], sep = ".") else
            s$rd$name[[k]]
          sc <- if (is.null(s$rd$scale)) ident else s$rd$scale[[k]]
          add(nmk, s$rd$value[[k]], p, s$term, sc, FALSE, j)
        }
      }
      if (drop[[i]]) next
      j <- numeric(njoint)
      j[offs[[pi]] + i] <- 1
      add(d$coef_names[[i]], cf[[i]], p, term_of[[i]], ident, FALSE, j)
    }
  }

  for (tl in tails) {
    z <- sst$zeta[[tl$u$term]]
    rd <- tryCatch(modelterms7::term_readable(tl$tm, z),
                   error = function(e) NULL)
    lb <- tryCatch(tl$tm@label, error = function(e) "")
    if (is.null(rd) || !length(rd$name)) {
      vals <- unlist(z)[tl$free]
      for (k in seq_along(tl$free)) {
        j <- numeric(njoint)
        j[tl$at + k] <- 1
        nmk <- if (nzchar(lb)) paste(lb, tl$free[[k]], sep = ".") else
          tl$free[[k]]
        add(nmk, vals[[k]], tl$u$param, tl$u$term, ident, FALSE, j)
      }
      next
    }
    # the Jacobian is over the term's WHOLE parameter vector and the joint
    # matrix carries the free ones alone, so the held columns are left out
    # and a quantity that reads one is marked rather than given the variance
    # of the rest
    fi <- match(tl$free, tl$nm)
    for (k in seq_along(rd$name)) {
      j <- numeric(njoint)
      j[tl$at + seq_along(fi)] <- rd$jacobian[k, fi]
      touches <- length(tl$held) > 0L &&
        any(rd$jacobian[k, match(tl$held, tl$nm)] != 0)
      nmk <- if (nzchar(lb)) paste(lb, rd$name[[k]], sep = ".") else
        rd$name[[k]]
      sc <- if (is.null(rd$scale)) ident else rd$scale[[k]]
      add(nmk, rd$value[[k]], tl$u$param, tl$u$term, sc, touches, j)
    }
  }

  # BY EQUATION, and stably. The joint coordinate order puts every design
  # block first and the structural tail last, so a model whose predictor is a
  # filter listed the scale's intercept before the mean's own parameters. The
  # Jacobian's COLUMNS keep the joint order, which is the order the variance
  # matrix is in; only the rows move.
  ord <- order(match(vapply(rows, function(z) z$parameter, character(1)),
                     params))
  rows <- rows[ord]
  list(name = vapply(rows, function(z) z$name, character(1)),
       value = vapply(rows, function(z) as.numeric(z$value), numeric(1)),
       parameter = vapply(rows, function(z) z$parameter, character(1)),
       term = vapply(rows, function(z) as.character(z$term)[[1L]],
                     character(1)),
       scale = lapply(rows, function(z) z$scale),
       held = vapply(rows, function(z) z$held, logical(1)),
       jacobian = do.call(rbind, lapply(rows, function(z) z$j)),
       n_design = total)
}


#' Invert a Matrix That Ought to Be Positive Definite
#'
#' @description
#' An inverse through the Cholesky factor, signaling an error naming the
#' matrix when a direction is flat.
#'
#' @details
#' A failure here is a statement about the fit and not about the arithmetic:
#' at a maximum the penalized information is positive definite, so a matrix
#' that is not says something about where the run stopped.
#'
#' The test is `lmin > tol * ref` on the smallest eigenvalue, and never
#' whether `chol()` raised. On an exactly singular matrix the pivot that
#' should be zero comes out positive or negative according to rounding, so
#' the second answer differs between platforms for one matrix.
#'
#' `ref` is the matrix's own scale, or the scale of the unpenalized
#' information where the caller holds it. That is what separates a flat
#' direction from the scale separation a large smoothing parameter
#' legitimately produces. Returning a pseudo-inverse instead would give a
#' standard error for a direction the data does not identify.
#'
#' The smallest eigenvalue is estimated, never computed, from LAPACK's
#' condition estimator (`dpocon`) read on the Cholesky factor the
#' inverse needs anyway: `rcond` is
#' \eqn{1/(\lVert A\rVert_1\lVert A^{-1}\rVert_1)}, so
#' `rcond * ||A||_1` is \eqn{1/\lVert A^{-1}\rVert_1}, which for a
#' symmetric matrix lies between \eqn{\lambda_{\min}/\sqrt{p}} and
#' \eqn{\lambda_{\min}}. The estimate therefore errs on the small side and
#' the test is conservative by at most a factor \eqn{\sqrt{p}}, plus
#' whatever the estimator's own slack is; the two cases it has to keep
#' apart are separated by some fifty orders of magnitude, so neither
#' reaches the other. It replaced a full eigendecomposition, which answers
#' the same question exactly and costs \eqn{O(p^3)} with a large constant:
#' measured at \eqn{p = 1022}, 1.18 s against the Cholesky's 0.25.
#'
#' The message names the directions, never the causes. A first version
#' offered two causes, the run not having reached a maximum and two columns
#' of the design carrying the same information, and on a Student t fitted to
#' `iris` neither was right: the design was full rank and the score was
#' 4e-5. What had
#' happened is the third and commonest case, a parameter drifting to where its
#' information vanishes, and no list of guesses would have said so. The
#' eigenvector of the smallest eigenvalue does: it is read off and the
#' coefficients that load on it are printed.
#'
#' @param A A square matrix.
#' @param what What the matrix is, for the message.
#' @param labels The names of the coefficients `A` is indexed by.
#'
#' @return The inverse.
#'

#' @keywords internal
solve_pd <- function(A, what, labels = NULL) {
  # The verdict comes from the smallest eigenvalue and not from whether
  # chol() raised. On a matrix with an exactly zero eigenvalue -- two columns
  # of the design carrying the same information is the ordinary way to get
  # one -- the pivot that should be zero comes out positive or negative
  # according to rounding, so the same fit was refused on Windows and
  # accepted on Linux and macOS. Reading the condition estimate off the
  # factor is a statement about the matrix instead: where the factorization
  # succeeds on a singular matrix by luck, the estimate is at the rounding
  # scale and the answer is the same on every platform.
  #
  # The test runs on the JACOBI-EQUILIBRATED matrix, D^-1/2 A D^-1/2 with
  # D its diagonal, which is what tells a flat direction from scale
  # separation whatever produced the separation. A reference scale passed
  # by the caller (the unpenalized information's largest diagonal) used to
  # do that, and it covered one source only: a smoothing parameter at 1e15
  # separates the scales and the reference forgave it, but the committed
  # working block of a break-point term carries auxiliary columns near
  # 1/(2 c d) -- 1e8 at the annealed floor, 1e16 on the information's
  # diagonal -- so the DESIGN itself set the reference and ordinary
  # curvature at 242 in the ordinary-scale directions read as flat
  # (measured, on the flagship jseg's own summary). Equilibration makes
  # the test per-direction: unit diagonal, so the smallest eigenvalue is
  # small only where a direction is flat AGAINST ITS OWN SCALE, and an
  # exact collinearity stays exactly singular. The inverse is recovered
  # through the same scaling, D^-1/2 (D^-1/2 A D^-1/2)^-1 D^-1/2.
  #
  # A non-positive diagonal entry means that coordinate carries no
  # information at all, which is the flat case stated directly.
  # A non-finite entry is what a parameter run out of its range leaves behind,
  # and it is not a matrix to decompose: the factorization raises its own
  # error there, which would replace the message below with one naming
  # neither the fit nor the direction. flat_directions() reports those rows
  # instead.
  # The test and the message's flat direction both need a dense matrix, and
  # there is no cheap sparse counterpart to replace them with. Densifying
  # HERE is deliberate rather than a lapse: this runs once, at vcov(), on a
  # p by p matrix, where the sparsity that matters is in the n by p design
  # and the per-iteration products taken against it. The same judgement is
  # recorded for the observed Hessian of a regime mixture, which is also
  # computed once and left in R.
  A <- as_dense(A)
  if (ncol(A) > 0L && all(is.finite(A)) && all(diag(A) > 0)) {
    s <- 1 / sqrt(diag(A))
    Ae <- A * tcrossprod(s)
    ch <- tryCatch(chol(Ae), error = function(e) NULL)
    if (!is.null(ch)) {
      # the 1-norm, which is what the condition estimator is expressed in
      anorm <- max(colSums(abs(Ae)))
      rc <- chol_rcond_cpp(ch, anorm)
      lmin <- if (is.na(rc)) 0 else rc * anorm
      if (is.finite(lmin) && lmin > 1e-12 * anorm) {
        # Inverted through the factor already in hand. A condition number of
        # 1e15 born of scale separation costs the well-determined directions
        # nothing here and the shrunk ones simply report variances near zero.
        out <- chol2inv(ch) * tcrossprod(s)
        return((out + t(out)) / 2)
      }
    }
  }
  stop(sprintf(paste0("%s is not positive definite, so there is no variance",
                      "\n  matrix at this point.%s\n  A fit can reach such a",
                      " point without failing: a parameter that runs\n  to",
                      " the edge of its space leaves no information behind,",
                      " and the score\n  is then small because the surface is",
                      " flat, not because it is a maximum.\n  Compare the",
                      " fitted parameters against the data before reading",
                      " anything\n  else."),
               what, flat_directions(A, labels)), call. = FALSE)
}


#' The Condition a Frozen Block Raises
#'
#' @description
#' A warning of its own class, so that a caller reporting the same thing in
#' its own channel can muffle it in place of repeating it.
#'
#' @details
#' [summary.StatmodFit()] calls [vcov()] more than
#' once and would raise the warning once per call; it says it once, as a
#' note.
#'
#' @param msg The message.
#'
#' @return A condition.
#'
#' @keywords internal
frozen_condition <- function(msg) {
  structure(class = c("statmod_frozen_block", "warning", "condition"),
            list(message = msg, call = NULL))
}


#' The Condition a Held Coordinate Warns Through
#'
#' @description
#' A classed warning, so a caller can silence or catch the one raised when
#' [vcov.StatmodFit()] holds a coordinate the information carries nothing
#' about, without silencing every other warning of a fit.
#'
#' @param msg The message.
#'
#' @return A condition of class `statmod_held_coord`.
#'
#' @seealso [uninformative_coords()], which finds them, and
#'   [vcov.StatmodFit()], which raises this.
#'
#' @keywords internal
held_condition <- function(msg) {
  structure(class = c("statmod_held_coord", "warning", "condition"),
            list(message = msg, call = NULL))
}

#' Which Coefficients Belong to a Block That Is Not a Jacobian
#'
#' @description
#' The positions of every coefficient of a term whose design block is a
#' working linearization with a frozen weight.
#'
#' @details
#' A discontinuous break-point term is fitted through a block whose weight is
#' held at the previous iterate, so the curvature that block carries is the
#' working model's, never the model's. Measured on a jump at 400
#' observations against a bootstrap of 200 resamples: the working information
#' gives the change of level and the auxiliary coordinate a standard error of
#' **exactly zero**, where the resamples give 0.063 and 0.540, and the
#' position read off them 1.8e-05 against 0.090. A zero looks like a number and is
#' worse than a gap, so those coefficients are left missing.
#'
#' The question is asked of the term through
#' [modelterms7::term_jacobian_block()] and never of its class,
#' so a construction whose block is a Jacobian keeps its inference: a
#' continuous [modelterms7::seg()], and a discontinuous one
#' smoothed by an [penalties7::abs_smoother()], both answer yes.
#'
#' @param spec The fitted specification.
#' @param lab The coefficient labels.
#'
#' @return A logical vector over the coefficients.
#'
#' @keywords internal
frozen_block <- function(spec, lab) {
  out <- rep(FALSE, nrow(lab))
  for (p in names(spec@terms)) {
    for (nm in names(spec@terms[[p]])) {
      ok <- tryCatch(modelterms7::term_jacobian_block(spec@terms[[p]][[nm]]),
                     error = function(e) TRUE)
      if (isTRUE(ok)) next
      out[lab$parameter == p & !is.na(lab$term) & lab$term == nm] <- TRUE
    }
  }
  out
}


#' Is a Matrix Worth Factorizing Sparsely?
#'
#' @description
#' Whether a symmetric matrix is large enough and sparse enough that a sparse
#' Cholesky beats a dense one. The question is put to the matrix and to
#' nothing else: no term, no model and no caller's preference enters.
#'
#' @details
#' The two conditions are the measured crossover. On the penalized
#' information of a random intercept over \eqn{m} levels at 20000
#' observations, the sparse route against the dense one, with coercion,
#' factorization, log-determinant and full inverse each timed with the
#' repetition loop sized by elapsed time:
#'
#' \tabular{rrrrr}{
#'   **m** \tab **p** \tab **density** \tab **whole route**
#'     \tab **inverse** \cr
#'   20 \tab 23 \tab 0.282 \tab 1.08x \tab **0.13x** \cr
#'   50 \tab 53 \tab 0.128 \tab 1.06x \tab **0.33x** \cr
#'   100 \tab 103 \tab 0.067 \tab 1.18x \tab 1.14x \cr
#'   200 \tab 203 \tab 0.034 \tab 1.74x \tab 2.9x \cr
#'   500 \tab 503 \tab 0.014 \tab 5.28x \tab 7.6x \cr
#'   1000 \tab 1003 \tab 0.007 \tab 2.50x \tab 11.2x \cr
#'   2000 \tab 2003 \tab 0.003 \tab 4.20x \tab 7.0x
#' }
#'
#' The whole-route column builds the matrix afresh on each repetition and is
#' the one the thresholds are read from. The inverse column is the like-for-
#' like comparison, each route carrying its own factorization; an earlier
#' version of it timed the sparse solves against a factor built once outside
#' the loop, which flattered the sparse side without changing where it loses.
#'
#' Below about a hundred coefficients the sparse route loses badly:
#' its fixed cost is the coercion and the S4 dispatch around it, which does not
#' shrink with the matrix. On the fully dense penalized information of a single
#' smooth (\eqn{p = 16}, density 1) it measures 0.01x, a hundred times
#' slower, and preventing that is what the size condition is for.
#'
#' **Both quantities are read off the matrix, and the first is its storage.**
#' A matrix held as a base matrix is refused whatever its zeros,
#' which reads like a test of the container instead of the mathematics, so
#' it is worth saying why it is neither an oversight nor a term test.
#' [statmod_information_at()] accumulates into the design's own kind,
#' so the penalized matrix is stored sparsely exactly when the design is, and
#' \pkg{modelterms7} builds a block sparse only when asked
#' (`sparse = TRUE`, whose default is `FALSE`). Measured on
#' `y ~ 0 + g + s(x)` over 400 levels at 20000 observations, whose
#' penalized matrix is 5 per cent nonzero either way: built dense the fit takes
#' 104.24 s and this factorization is **0.16 per cent** of it, the time
#' being in the \eqn{O(np^2)} products against a dense design
#' (`statmod_information_at` 48.8 per cent, `crossprod` 57.2 per cent
#' of self time); built sparse the same fit takes 2.19 s. So where the storage
#' is dense the factorization is not what a fit is spending its time on, and
#' coercing a dense \eqn{p \times p} matrix here to save a share of that size
#' would cost more than it returns.
#'
#' **The like-for-like comparison is the one that says this is not a term
#' test**, and it is the check `piano_lme4.txt` section 5 asks for. With
#' every design built the same way, this route is worth 1.38x on
#' `0 + g + s(x)` over 400 levels, 1.33x on `random(~1|g)` over 500
#' and 1.07x on `s(x, by = g)` over 60: an unpenalized indicator block,
#' a random effect and a factor-`by` smooth, gaining together and in the
#' order their sizes predict. Nothing here asks which term or which family
#' produced the matrix.
#'
#' @param M A symmetric matrix.
#' @param min_dim The smallest order worth the fixed cost.
#' @param max_density The largest fraction of nonzeros worth it.
#'
#' @return A single logical.
#'
#' @seealso [pd_factor()]
#'
#' @keywords internal
worth_sparse <- function(M, min_dim = 100L, max_density = 0.10) {
  if (!isS4(M)) return(FALSE)
  p <- ncol(M)
  if (!length(p) || p < min_dim) return(FALSE)
  nz <- tryCatch(length(M@x), error = function(e) NA_integer_)
  if (is.na(nz)) return(FALSE)
  nz / (as.numeric(p) * p) <= max_density
}


#' The Smallest Eigenvalue of a Sparse Factor's Matrix, Estimated
#'
#' @description
#' \eqn{1/\lVert A^{-1}\rVert_1} from a sparse Cholesky factor, the same
#' quantity LAPACK's `dpocon` produces from a dense one.
#'
#' @details
#' The sparse route needs a condition estimate of its own, and it cannot
#' borrow the dense one: `chol_rcond_cpp` reads a dense triangular
#' factor. `Matrix::rcond` is not the answer either -- measured, it costs
#' 10.3 ms at p = 503 and 500 ms at p = 2003, more than the dense
#' factorization the sparse route exists to replace. Higham's one-norm
#' estimator applied to the factor's own solves costs 0.58 ms at p = 53 and
#' 0.80 ms at p = 1003, nearly flat, because it is a handful of triangular
#' solves and an R loop around them.
#'
#' For a symmetric matrix \eqn{\lVert A^{-1}\rVert_1 \ge \lVert
#' A^{-1}\rVert_2 = 1/\lambda_{\min}}, so the quantity returned is at or below
#' the smallest eigenvalue, and the estimator's own error is a further
#' underestimate of the norm in the other direction. It is used exactly as the
#' dense estimate is: to separate a matrix comfortably positive definite from
#' one that is not, two situations that differ by some fifteen orders of
#' magnitude here (measured on a design with two identical columns, 1.4e-14
#' relative to the matrix's scale, against 1.5e3 for a hyperparameter driven
#' to 1e15). A factor of two either way cannot move that verdict, which is the
#' argument already recorded for the dense estimator.
#'
#' @param L A `CHMfactor`.
#' @param p The order of the matrix.
#'
#' @return A single number, or `NA_real_` where the estimate failed.
#'
#' @seealso [pd_factor()], [pd_logdet()]
#'
#' @keywords internal
sparse_lmin <- function(L, p) {
  ax <- function(x) as.matrix(Matrix::solve(L, as.matrix(x)))
  est <- tryCatch(
    Matrix::onenormest(A.x = ax, At.x = ax, n = p, t = 1, silent = TRUE)$est,
    error = function(e) NA_real_)
  if (!is.finite(est) || est <= 0) return(NA_real_)
  1 / est
}


#' Factorize a Penalized Information Once
#'
#' @description
#' The Cholesky factor of a matrix a Laplace approximation needs positive
#' definite, together with its log-determinant, in whichever storage the
#' matrix itself calls for.
#'
#' @details
#' This is the one place the penalized matrix is factorized. The criterion
#' wants its log-determinant, the gradient wants the mode's movement and the
#' Hessian wants both plus the inverse; before this existed the criterion and
#' [ctx_penalized()] each factorized the same matrix at the same point,
#' which at \eqn{p = 503} was 12.4 ms spent twice.
#'
#' **The verdict is unchanged and so is its property.** Whether the
#' matrix is accepted never turns on whether a factorization raised: where the
#' cheap test is inconclusive the eigendecomposition answers about the matrix.
#' The sparse route carries its own condition estimate
#' ([sparse_lmin()]) rather than the dense one, and falls back to
#' the dense route where that estimate cannot be formed, so a refusal is
#' reached by the same reasoning on either storage.
#'
#' @param M A symmetric matrix, sparse or dense.
#' @param scale A reference magnitude, as [pd_logdet()] takes.
#'
#' @return A list with `logdet`, `ok`, `factor` and
#'   `sparse`. The factor is `NULL` where the answer came from the
#'   eigendecomposition.
#'
#' @seealso [pd_logdet()], [ctx_penalized()]
#'
#' @keywords internal
pd_factor <- function(M, scale = NULL) {
  no <- function() list(logdet = NA_real_, ok = FALSE, factor = NULL,
                        sparse = FALSE)
  if (!length(ncol(M)) || !ncol(M)) return(no())

  if (worth_sparse(M)) {
    p <- ncol(M)
    if (!all(is.finite(M@x))) return(no())
    anorm <- max(Matrix::colSums(abs(M)))
    ref <- if (is.null(scale) || !is.finite(scale) || scale <= 0) anorm
           else min(scale, anorm)
    Ms <- tryCatch(Matrix::forceSymmetric(M), error = function(e) NULL)
    # CHOLMOD reports a matrix it cannot factorize as a WARNING and not as an
    # error, where the dense route raises and is caught. A hyperparameter the
    # search should step away from is not something to warn a caller about --
    # the criterion returns NULL there and the search backtracks -- so the
    # warning is suppressed and the verdict is read off the return value, as
    # it is on the dense route.
    L <- if (is.null(Ms)) NULL else
      suppressWarnings(tryCatch(Matrix::Cholesky(Ms, LDL = FALSE, super = NA),
                                error = function(e) NULL))
    if (!is.null(L)) {
      lmin <- sparse_lmin(L, p)
      if (is.finite(lmin) && lmin > sqrt(.Machine$double.eps) * ref) {
        ld <- tryCatch(2 * as.numeric(Matrix::determinant(
          L, logarithm = TRUE, sqrt = TRUE)$modulus),
          error = function(e) NA_real_)
        if (is.finite(ld)) {
          return(list(logdet = ld, ok = TRUE, factor = L, sparse = TRUE))
        }
      }
    }
    # inconclusive: the dense route asks the matrix, which is what settles it
  }

  d <- pd_logdet_dense(as_dense(M), scale)
  c(d, list(sparse = FALSE))
}


#' The Log-Determinant of a Penalized Information, Robustly
#'
#' @description
#' \eqn{\log|M|} for a matrix a Laplace approximation needs to be positive
#' definite, by the cheap route where that is safe and a costlier one where
#' it is not.
#'
#' @details
#' **Why not `chol()` alone.** A marginal criterion read the
#' determinant off `chol(M)` and reported the criterion as nonexistent
#' whenever the factorization raised. At a condition number near the rounding
#' floor, whether it raises is decided by the arithmetic and never by the
#' matrix:
#' measured on a hierarchical score-driven panel, \eqn{K+S} had a smallest
#' eigenvalue of 4.3e-11 against a condition number of 8.0e15, and the outer
#' search then backtracked through a dozen points reported unavailable towards
#' one that had been available a moment earlier. The same doubt this package
#' already records for [solve_pd()] and for basis7's rank tests.
#'
#' **The three routes.** The factorization is tried first, being O(p^3/3)
#' and the common case. Where it succeeds, LAPACK's condition estimator reads
#' the smallest eigenvalue off the factor already in hand for O(p^2), and a
#' matrix comfortably away from the floor is accepted with the determinant the
#' factor gives. The eigendecomposition is computed only where that test is
#' inconclusive, or where the factorization raised at all. It costs more and
#' answers about the matrix: a factorization that failed by rounding luck on
#' a
#' matrix that is in fact positive definite is recovered there, and one that
#' is genuinely rank deficient is refused deterministically rather than
#' according to the platform.
#'
#' @param M A symmetric matrix.
#' @param scale A reference magnitude, as [solve_pd()] takes: the
#'   unpenalized information's own scale, so that a hyperparameter legitimately
#'   sent to 1e15 is told apart from a flat direction.
#'
#' @return A list with `logdet` and `ok`, or `ok = FALSE` where
#'   the matrix is not positive definite.
#'
#' @seealso [solve_pd()], [statmod_marginal()]
#'
#' @keywords internal
pd_logdet <- function(M, scale = NULL) {
  r <- pd_factor(M, scale)
  r$factor <- NULL
  r$sparse <- NULL
  r
}


#' The Dense Route of pd_logdet
#'
#' @description
#' The three routes described at [pd_logdet()], on a dense matrix.
#'
#' @details
#' Split out so that [pd_factor()] can reach it as the fallback of
#' the sparse route without restating the verdict: there is one place that
#' decides whether a matrix is positive definite, and one set of thresholds.
#'
#' @param M A dense symmetric matrix.
#' @param scale A reference magnitude.
#'
#' @return A list with `logdet`, `ok` and, on a refusal reached
#'   through the eigendecomposition, `min_ev` and `max_ev`.
#'
#' @seealso [pd_logdet()]
#'
#' @keywords internal
pd_logdet_dense <- function(M, scale = NULL) {
  if (!ncol(M) || !all(is.finite(M))) return(list(logdet = NA_real_, ok = FALSE))
  anorm <- max(colSums(abs(M)))
  ref <- if (is.null(scale) || !is.finite(scale) || scale <= 0) anorm
         else min(scale, anorm)
  # ⚠️ The two thresholds are NOT solve_pd()'s, and the difference is the
  # quantity. Inverting a matrix at a condition number of 1e14 loses most of
  # the answer, so solve_pd() refuses there; a log-DETERMINANT is a sum of
  # logarithms and survives it -- measured, chol and eigen agree to six
  # significant figures at 1e14 (-644.725631 against -644.725128). What
  # matters here is only that every eigenvalue is genuinely POSITIVE, and an
  # eigenvalue below eps times the largest is not distinguishable from zero
  # in double precision whatever the factorization reports.
  ch <- tryCatch(chol(M), error = function(e) NULL)
  if (!is.null(ch)) {
    rc <- tryCatch(chol_rcond_cpp(ch, anorm), error = function(e) NA_real_)
    lmin <- if (is.na(rc)) 0 else rc * anorm
    # comfortably clear of the floor: chol's answer is its own evidence, and
    # this is the common case, so the eigendecomposition is never computed
    if (is.finite(lmin) && lmin > sqrt(.Machine$double.eps) * ref) {
      # the factor travels with the answer so that a caller wanting the
      # inverse as well does not compute a second one of the same matrix
      return(list(logdet = 2 * sum(log(diag(ch))), ok = TRUE, factor = ch))
    }
  }
  # inconclusive or refused: ask the MATRIX rather than the arithmetic
  ev <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev) || !all(is.finite(ev))) {
    return(list(logdet = NA_real_, ok = FALSE))
  }
  if (min(ev) > .Machine$double.eps * max(ev)) {
    return(list(logdet = sum(log(ev)), ok = TRUE))
  }
  list(logdet = NA_real_, ok = FALSE, min_ev = min(ev), max_ev = max(ev))
}




#' Which Coefficients a Singular Curvature Is Flat In
#'
#' @description
#' The eigenvector of the smallest eigenvalue, reported as the coefficients
#' that load on it.
#'
#' @param A A square matrix.
#' @param labels Its coefficient names.
#'
#' @return A single string, empty when nothing can be said.
#'
#' @keywords internal
flat_directions <- function(A, labels) {
  if (is.null(labels) || !length(labels)) return("")
  # a parameter far enough out makes its own row non-finite rather than merely
  # small, and then there is no eigenvector to read: the rows themselves are
  # the answer, and they are the more direct one
  bad <- labels[apply(!is.finite(A), 1L, any)]
  if (length(bad)) {
    return(sprintf(paste0("\n  Its entries are not finite in the rows of: %s",
                          "\n  which is what a parameter that has run out of",
                          " its range looks like here."),
                   paste(bad, collapse = ", ")))
  }
  e <- tryCatch(eigen(A, symmetric = TRUE), error = function(e) NULL)
  if (is.null(e)) return("")
  k <- which.min(e$values)
  v <- abs(e$vectors[, k])
  hit <- labels[v > 0.2 * max(v)]
  if (!length(hit)) return("")
  sprintf(paste0("\n  It is flat along a direction carried by: %s\n  (its",
                 " smallest eigenvalue is %s, against %s at the largest)."),
          paste(hit, collapse = ", "), format(signif(e$values[k], 3)),
          format(signif(max(e$values), 3)))
}


#' Which Coordinates a Singular Information Carries Nothing About
#'
#' @description
#' The positions whose row and column of a penalized information are empty,
#' so that holding them and inverting the rest reports the variance of every
#' other coordinate exactly.
#'
#' @details
#' A combination \eqn{c'\hat\beta} has a finite asymptotic variance exactly
#' where \eqn{c} is orthogonal to the null space \eqn{\mathcal{N}} of
#' \eqn{K}, and then \eqn{\mathrm{Var}(c'\hat\beta) = c'K^{+}c} with
#' \eqn{K^{+}} the Moore-Penrose inverse. Reading `diag(K^+)` is therefore
#' not enough on its own: for a \eqn{c} that is NOT orthogonal it returns a
#' number where the truth is infinite. The rule is in two parts, invert away
#' from \eqn{\mathcal{N}} and report nothing on it, and this function is the
#' second part.
#'
#' A candidate is a coordinate with no curvature of its own. Its diagonal is
#' not finite, which is what a parameter run out of its range leaves behind;
#' or it is not positive, which is [solve_pd()]'s own reading of that entry
#' and is a statement about the coordinate whatever scale it is on; or, among
#' the coordinates whose diagonal IS positive, it carries the whole of a null
#' direction of the JACOBI-EQUILIBRATED matrix. The equilibration is what
#' tells a flat direction from scale separation, and requiring one coordinate
#' to carry the whole direction is what excludes the collinear case, where
#' the null vector is \eqn{(e_1 - e_2)/\sqrt{2}} and each of the two carries
#' half of it.
#'
#' **A candidate is held only where holding it disturbs nothing**, and the
#' amount it disturbs is the Schur correction it removes from the others,
#' \eqn{K_{Aj}K_{jj}^{-1}K_{jA}}. That test has two branches, because it has
#' two scales.
#'
#' Where \eqn{K_{jj}} is POSITIVE the coordinate has a scale of its own and
#' the correction, read against \eqn{K_{AA}} in equilibrated units where that
#' block has a unit diagonal, is
#' \deqn{\max_k \Bigl(\frac{|K_{jk}|}{\sqrt{K_{kk}}}\Bigr)^{2}
#'       \Big/ K_{jj} \;\le\; \texttt{schur}.}
#' It is dimensionless and needs no reference: a coordinate measured a
#' thousand times smaller than its neighbours has a small row AND a small
#' diagonal, and the ratio sees through both. Scaling one down by
#' \eqn{10^{-14}} leaves the ratio at one, and nothing is held.
#'
#' Where \eqn{K_{jj}} is NOT positive the coordinate has no scale of its own
#' -- the information in that direction is zero rather than small -- so that
#' ratio is \eqn{0/0} and says nothing. The only scale left is the matrix's,
#' and the question becomes whether the row sits at the level its entries
#' were accumulated to: \eqn{\max_k |K_{jk}| \le \texttt{row\_tol}\max|K|}.
#' Measured on a count model whose dispersion left its range, \eqn{K_{jj}} is
#' exactly zero, the row's largest entry is \eqn{1.2\times10^{-15}} against a
#' matrix whose largest is 2562, and \eqn{\lVert K_{j\cdot}\rVert/\lVert
#' K\rVert} is \eqn{4.4\times10^{-19}}.
#'
#' Where neither branch passes, dropping the coordinate would leave
#' \eqn{K_{AA}^{-1}}, the variance CONDITIONAL on that coefficient being
#' known, which is smaller than the marginal one; nothing is held and the
#' caller's refusal stands.
#'
#' @param A A penalized information, over the coefficients and any structural
#'   tail.
#' @param tol The relative eigenvalue below which a direction is flat.
#' @param share How much of the null space one coordinate must carry to be
#'   held, as a share of one.
#' @param schur The largest Schur correction, in equilibrated units, that
#'   counts as no disturbance to the coordinates that are kept, where the
#'   candidate has a positive diagonal.
#' @param row_tol The share of the matrix's largest entry below which a row
#'   is read as accumulated rounding, where the candidate does not.
#'
#' @return An integer vector of positions in `A`, possibly empty.
#'
#' @seealso [solve_pd()], which refuses the matrix these come from, and
#'   [vcov.StatmodFit()], which holds them.
#'
#' @keywords internal
uninformative_coords <- function(A, tol = 1e-10, share = 1 - 1e-6,
                                 schur = 1e-8, row_tol = 1e-12) {
  A <- as_dense(A)
  p <- ncol(A)
  if (!p) return(integer(0))
  # BY THE DIAGONAL, which is boundary_coords()' rule: a coordinate at a
  # bound makes its whole row non-finite, cross terms included, so a row
  # test marks its neighbours as well. Whatever is left non-finite among
  # the survivors goes too, there being no matrix there either.
  bad <- boundary_coords(A)
  rest <- setdiff(seq_len(p), bad)
  if (length(rest)) {
    extra <- rest[apply(!is.finite(A[rest, rest, drop = FALSE]), 1L, any)]
    bad <- c(bad, extra)
    rest <- setdiff(rest, extra)
  }
  if (!length(rest)) return(sort(bad))
  d <- diag(A)
  # A NON-POSITIVE DIAGONAL is the flat case stated directly: that
  # coordinate has no curvature of its own, whatever scale it is on, and
  # the sign of an entry at the rounding level is a coin toss. A positive
  # one equilibrates and is read as a direction of the correlation matrix.
  cand <- rest[d[rest] <= 0]
  rest <- setdiff(rest, cand)
  if (length(rest) >= 2L) {
    B <- A[rest, rest, drop = FALSE]
    s <- 1 / sqrt(diag(B))
    Be <- B * tcrossprod(s)
    e <- tryCatch(eigen(Be, symmetric = TRUE), error = function(e) NULL)
    if (!is.null(e)) {
      ev <- abs(e$values)
      flat <- ev <= tol * max(ev)
      if (any(flat)) {
        # each coordinate's share of the null space: one where the
        # direction IS that coordinate, a half where two of them share it
        load <- rowSums(e$vectors[, flat, drop = FALSE]^2)
        cand <- c(cand, rest[load > share])
        rest <- setdiff(rest, rest[load > share])
      }
    }
  }
  if (!length(cand)) return(sort(bad))
  # HOLDING IS EXACT ONLY WHERE IT DISTURBS NOTHING, and the amount it
  # disturbs is the Schur correction A_Aj A_jj^-1 A_jA it removes from the
  # others. Read in equilibrated units, where the kept block has a unit
  # diagonal, that is max_k (|A_jk|/sqrt(A_kk))^2 / |A_jj|, and it needs no
  # reference scale: a coordinate measured a thousand times smaller than
  # its neighbours has a small row AND a small diagonal, and the ratio sees
  # through both. Where the correction is not negligible, dropping the
  # coordinate would report the variance CONDITIONAL on it as if it were
  # the marginal one, so nothing is held and the caller's refusal stands.
  if (length(rest)) {
    sk <- sqrt(d[rest])
    scale <- max(abs(A[rest, rest]))
    ok <- vapply(cand, function(j) {
      row <- abs(A[j, rest])
      if (!all(is.finite(row))) return(FALSE)
      if (d[j] > 0) {
        # a scale of its own: the correction is r^2/A_jj, dimensionless
        r <- max(row / sk)
        r * r <= schur * d[j]
      } else {
        # NO scale of its own, so r^2/A_jj is 0/0 and says nothing. The
        # only scale left is the matrix's, and the question is whether the
        # row sits at the level its entries were accumulated to.
        max(row) <= row_tol * scale
      }
    }, logical(1))
    cand <- cand[ok]
  }
  sort(c(bad, cand))
}


#' @title Confidence Intervals for a Fit
#' @name confint.StatmodFit
#' @description
#' Wald intervals for the coefficients of every distribution parameter.
#' @details
#' The interval is symmetric about the estimate and needs no mapping back. A
#' coefficient of a linear predictor is unbounded whatever the distribution
#' parameter it belongs to, the link having already carried that parameter onto
#' the whole line, so the scale the interval is built on is the scale the
#' quantity lives on. What the interval does not do is respect a bound on the
#' parameter itself; for that, map an interval for the predictor through the
#' inverse link at the covariate values of interest.
#'
#' The variance comes from [vcov.StatmodFit()], so the same two
#' conventions apply, and a coefficient a kinked penalty set to zero has
#' `NA` in place of an interval.
#' @param object A [StatmodFit()].
#' @param parm Which coefficients: a distribution parameter's name, a vector of
#'   `parameter:coefficient` labels, or `NULL` for all of them.
#' @param level The confidence level.
#' @param type Passed to [vcov.StatmodFit()].
#' @param ... Passed to [vcov.StatmodFit()].
#' @return A data frame with the parameter, the term, the coefficient, the
#'   estimate, its standard error and the two limits.
#' @seealso [vcov.StatmodFit()], [summary.StatmodFit()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(80))
#' dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' confint(fit)
#' confint(fit, "sigma")
#' @keywords internal
confint.StatmodFit <- function(object, parm = NULL, level = 0.95,
                               type = c("bayesian", "frequentist"),
                               readable = TRUE, ...) {
  type <- match.arg(type)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  spec <- object@spec
  design <- statmod_design(spec)
  z <- stats::qnorm((1 + level) / 2)
  V <- vcov(object, readable = readable, type = type, ...)
  se <- sqrt(diag(V))
  if (isTRUE(readable)) {
    rj <- readable_joint(spec, design, object)
    est <- rj$value
    # EACH INTERVAL ON THE SCALE THAT KEEPS ITS QUANTITY IN ITS OWN SET, and
    # mapped back, as every interval in the toolkit is: a loading is
    # positive and rides a log, an autoregressive coefficient is confined to
    # no interval a scalar link expresses and rides the identity. Read on
    # the raw scale a positive quantity routinely gets a negative lower end.
    lo <- hi <- rep(NA_real_, length(est))
    for (i in seq_along(est)) {
      if (!is.finite(se[[i]])) next
      g <- rj$scale[[i]]
      t0 <- linkfunctions7::linkfun(g, est[[i]])
      st <- se[[i]] * abs(linkfunctions7::dlinkfun(g, est[[i]]))
      ends <- sort(c(linkfunctions7::linkinv(g, t0 - z * st),
                     linkfunctions7::linkinv(g, t0 + z * st)))
      lo[[i]] <- ends[[1L]]
      hi[[i]] <- ends[[2L]]
    }
    out <- data.frame(parameter = rj$parameter, term = rj$term,
                      coefficient = rj$name, estimate = est, se = se,
                      lower = lo, upper = hi, stringsAsFactors = FALSE)
    rownames(out) <- rownames(V)
    if (is.null(parm)) return(out)
    keep <- rownames(out) %in% parm | out$parameter %in% parm |
      out$term %in% parm | out$coefficient %in% parm
    if (!any(keep)) {
      stop(sprintf(paste0("'parm' matched nothing. It takes a distribution",
                          " parameter\n  (%s), a term's name, a quantity's",
                          " name, or a label of the\n  form",
                          " 'parameter:name'."),
                   paste(spec@distrib@params, collapse = ", ")),
           call. = FALSE)
    }
    return(out[keep, , drop = FALSE])
  }
  lab <- coef_labels(spec, design)
  est <- unlist(object@coefficients[spec@distrib@params], use.names = FALSE)
  se <- se[seq_len(nrow(lab))]
  out <- data.frame(lab[, c("parameter", "term", "coefficient")],
                    estimate = est, se = se,
                    lower = est - z * se, upper = est + z * se,
                    stringsAsFactors = FALSE)
  rownames(out) <- rownames(lab)
  if (is.null(parm)) return(out)

  keep <- rownames(out) %in% parm | out$parameter %in% parm |
    out$term %in% parm
  if (!any(keep)) {
    stop(sprintf(paste0("'parm' matched nothing. It takes a distribution",
                        " parameter\n  (%s), a term's name, or a label of",
                        " the form 'parameter:coefficient'."),
                 paste(spec@distrib@params, collapse = ", ")), call. = FALSE)
  }
  out[keep, , drop = FALSE]
}
S7::method(confint, StatmodFit) <- confint.StatmodFit


#' What Kind of Block a Term Reports As
#'
#' @description
#' Which of the four readings a term gets in a summary: its coefficients, a
#' smooth's linear part and smoothing parameter, a random effect's variance
#' parameters, or a selection's survivors.
#'
#' @details
#' The classification is by the term's class and by its penalties, not by its
#' label, so a term given a name of its own is read the same way. The
#' penalties are the ones the term declares through
#' [modelterms7::term_penalties()], so a term penalized over part of
#' its parameters, a segmented term's changes or a filter's deviations, is
#' read as penalized, never as parametric, and is a selection when any
#' of its penalties has a kink.
#'
#' @param term A built term.
#'
#' @return One of `"structural"`, `"breakpoint"`,
#'   `"parametric"`, `"smooth"`, `"random"`,
#'   `"selection"`, `"penalized"`.
#'
#' @keywords internal
term_block_kind <- function(term) {
  # a structural term is asked about FIRST OF ALL: it contributes no design
  # columns, so what it is has nothing to do with whether it carries a
  # penalty. Asked later it came back "parametric" whenever it carried none,
  # which put a term with no columns among the terms whose columns are one
  # unpenalized block.
  if (S7::S7_inherits(term, modelterms7::structural_term)) {
    return("structural")
  }
  # a break-point term is asked about next, before its penalties are
  # looked at: its block is a working linearization whose coefficients are
  # not the quantities of the model, so it wants a section of its own
  # whether or not a development of its coefficients carries a penalty
  if (S7::S7_inherits(term, modelterms7::SegTerm)) return("breakpoint")
  # a LABELLED random effect is asked about before its penalties too, and
  # for the same reason: it declares none, its coefficients being covered
  # by the covariance class's, and read after the penalties it came back
  # "parametric" -- forty grouping indicators printed one per line under
  # a heading that says they are an unpenalized block
  if (S7::S7_inherits(term, modelterms7::RandomTerm) &&
      !is.na(modelterms7::term_tag(term))) {
    return("random")
  }
  ent <- modelterms7::term_penalties(term)
  if (!length(ent)) return("parametric")
  if (S7::S7_inherits(term, modelterms7::RandomTerm)) return("random")
  if (S7::S7_inherits(term, modelterms7::SmoothTerm)) return("smooth")
  if (any(vapply(ent, function(e) penalty_has_kink(e$penalty), logical(1)))) {
    return("selection")
  }
  "penalized"
}


#' Which Coefficients of a Smooth Are the Linear Part
#'
#' @description
#' `TRUE` for the columns a Demmler-Reinsch smooth carries its linear
#' effect in, which are the ones worth printing.
#'
#' @details
#' The rest of the block are coefficients of an orthonormal basis of the
#' wiggly part; individually they say nothing, and what they say jointly is the
#' effective degrees of freedom, which the block header reports instead.
#'
#' The question is asked of the term's own specification (`spec$linear`)
#' and never of a suffix in a coefficient's name, a name being a label and
#' this is a fact about the construction.
#'
#' @param term A built smooth term.
#' @param k The number of columns in its block.
#'
#' @return A logical vector of length `k`.
#'
#' @keywords internal
smooth_linear_cols <- function(term, k) {
  out <- rep(FALSE, k)
  sp <- tryCatch(term@spec, error = function(e) NULL)
  if (is.list(sp) && isTRUE(sp$linear) && k > 0L) out[1L] <- TRUE
  out
}


#' A Summary of a Fitted Model
#'
#' @description
#' What [summary.StatmodFit()] returns: the blocks of each
#' distribution parameter, the degrees of freedom, the information criteria and
#' whatever has to be said about how the numbers should be read.
#'
#' @param call The fit's call.
#' @param distrib_name The distribution's name.
#' @param n_obs The number of observations.
#' @param tables A named list, one entry per distribution parameter, each a
#'   list of block records.
#' @param edf The per-term degrees of freedom.
#' @param loglik The maximized log-likelihood.
#' @param df The effective degrees of freedom in total.
#' @param aic,bic The information criteria.
#' @param converged Whether every loop stopped on its own rule.
#' @param elapsed The elapsed time, in seconds.
#' @param level The confidence level the intervals were built at.
#' @param type Which variance convention was used.
#' @param notes Character vector of things the reader has to know.
#'
#' @return An object of class `StatmodSummary`.
#'
#' @seealso [summary.StatmodFit()]
#'
#' @examples
#' dd <- data.frame(y = rnorm(30), x = runif(30))
#' S7::S7_inherits(summary(statmod(y ~ x,
#'                                 distributions7::gaussian1_distrib(), dd)),
#'                 StatmodSummary)
#'
#' @name StatmodSummary-class
#' @aliases StatmodSummary
#' @keywords internal
#' @export
StatmodSummary <- S7::new_class("StatmodSummary",
  properties = list(
    call = S7::class_any,
    distrib_name = S7::class_character,
    n_obs = S7::class_numeric,
    tables = S7::class_list,
    # the link each equation is written on, named by parameter. Every
    # coefficient in a block is a coefficient of the LINEAR PREDICTOR, so
    # what it means for the parameter depends on the link, and a summary
    # that does not say which one leaves the reader to remember the
    # family's defaults.
    links = S7::new_property(S7::class_character, default = character(0)),
    edf = S7::class_any,
    structural = S7::class_any,
    loglik = S7::class_numeric,
    df = S7::class_numeric,
    aic = S7::class_numeric,
    bic = S7::class_numeric,
    converged = S7::class_logical,
    elapsed = S7::class_numeric,
    level = S7::class_numeric,
    type = S7::class_character,
    notes = S7::class_character,
    # the reading of statmod_certificate(), or NULL where it could not be
    # taken. A list rather than columns: its four readings are of different
    # kinds and one of them is a set of names.
    certificate = S7::new_property(S7::class_any, default = NULL)
  )
)


#' @title Summarize a Fitted Model
#' @name summary.StatmodFit
#' @description
#' One coefficient table per distribution parameter, carrying the estimate,
#' its standard error, the Wald statistic, its p-value and an interval,
#' together with the degrees of freedom, the information criteria and the
#' qualifications the numbers carry.
#' @details
#' **Each distribution parameter is read as blocks, not as one list of
#' coefficients**, because most of a fitted model's coefficients are not
#' quantities anybody reads. The blocks are
#' \describe{
#'   \item{the parametric terms}{every unpenalized term together, one row per
#'     coefficient, which is the ordinary table.}
#'   \item{one block per smooth}{the linear component's coefficient where the
#'     construction carries one, the smoothing parameter, and the effective
#'     degrees of freedom. The coefficients of the wiggly part are not shown:
#'     they are coordinates in an orthonormal basis and say nothing one at a
#'     time, while what they say together is the edf.}
#'   \item{one block per random effect}{the parameters of the effects'
#'     distribution, usually called the variance component, and the
#'     edf. Not the effects themselves, of which there is one per level.}
#'   \item{one block per selection}{a lasso, a SCAD or an MCP: its
#'     hyperparameters, how many coefficients survived, and those coefficients.
#'     The ones set exactly to zero are counted, not listed.}
#'   \item{one block per other penalized term}{its coefficients, which stay
#'     interpretable under a ridge, together with its hyperparameters.}
#' }
#'
#' **A hyperparameter is the first row of its block**, since it governs
#' every coefficient under it, and the cell where its standard error would be
#' says what put the value there. One estimated by [reml()] or
#' [ml()] maximizes a twice differentiable criterion, so it carries
#' a standard error and an interval, both read on the free scale its link
#' defines and mapped back ([statmod_hyper_vcov()]). One chosen by
#' [aic()], [bic()] or [cv()] over a kinked
#' penalty is the argument of a minimum over a grid, so the row names the
#' criterion and leaves the remaining columns empty: there is no curvature at
#' such a point to read a standard error from. One the caller set is marked
#' fixed.
#'
#' **What a Wald p-value means here depends on the row**, and the summary
#' says which is which, in place of printing one column and leaving it at
#' that.
#' For an unpenalized coefficient it is the usual thing. For a coefficient in a
#' penalized block it is conditional on the smoothing parameter, which was not
#' estimated jointly with it, and it does not account for the shrinkage of the
#' estimate towards zero. For a coefficient a kinked penalty selected, the row
#' exists only because that coefficient survived the selection, and a naive
#' interval there under-covers.
#'
#' **The degrees of freedom** are the effective ones, summed over the
#' terms, so a penalized term counts what it spends instead of how many
#' columns it has. The information criteria are built on that count.
#' @param object A [StatmodFit()].
#' @param level The confidence level.
#' @param type Which variance matrix: passed to [vcov.StatmodFit()].
#' @param expected Which information the standard errors are read off, passed
#'   to [vcov.StatmodFit()]. `NULL`, the default, takes the expected one where
#'   the fit inverted it and the family writes it out, and the observed
#'   Hessian otherwise; `TRUE` or `FALSE` names one.
#' @param approx How the expected information is approximated for a family
#'   with no closed form, passed to [vcov.StatmodFit()]: `"opg"` (the
#'   default), `"bartlett"`, `"integrate"` or `"mc"`.
#' @param correct Whether the degrees of freedom carry what the estimation
#'   of the hyperparameters cost. The ordinary count reads them as known,
#'   and they were chosen from the same data, so a criterion built on it is
#'   too generous. See [statmod_edf_correction()]. Defaults to
#'   `FALSE` because it changes a number a reader may be comparing with
#'   an earlier fit; it is zero where no hyperparameter was estimated.
#' @param ... Passed to [vcov.StatmodFit()].
#' @return A [StatmodSummary()].
#' @seealso [vcov.StatmodFit()], [confint.StatmodFit()]
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(120))
#' dd$y <- 1 + 2 * dd$x + rnorm(120, sd = 0.4)
#' summary(statmod(y ~ x | sigma ~ x,
#'                 distributions7::gaussian1_distrib(), dd))
#' @keywords internal
summary.StatmodFit <- function(object, level = 0.95,
                               type = c("bayesian", "frequentist"),
                               expected = NULL,
                               approx = c("opg", "bartlett", "integrate", "mc"),
                               correct = FALSE, ...) {
  type <- match.arg(type)
  approx <- match.arg(approx)
  # A TERM WHOSE BLOCK IS A WORKING LINEARIZATION says so ONCE, as a note,
  # rather than once per call to vcov() -- of which this function makes more
  # than one. The condition carries a class of its own, so muffling it cannot
  # swallow another warning raised on the way.
  frozen_msg <- character(0)
  catch_frozen <- function(expr) {
    withCallingHandlers(expr, statmod_frozen_block = function(cnd) {
      frozen_msg <<- unique(c(frozen_msg, conditionMessage(cnd)))
      invokeRestart("muffleWarning")
    })
  }
  # ON THE RAW COORDINATES, both of them. This function does its own reading
  # of what a term is about -- readable_rows() for a design block, the
  # structural table for a term with none -- and it INDEXES the variance
  # matrix by the design's own coefficient names. Given the readable matrix
  # those keys are not in it, and every delta-method standard error a
  # break-point term reports came back missing.
  ci <- catch_frozen(confint(object, level = level, type = type,
                             expected = expected, approx = approx,
                             readable = FALSE, ...))
  spec <- object@spec
  design <- statmod_design(spec)
  ci$statistic <- ci$estimate / ci$se
  ci$p_value <- 2 * stats::pnorm(-abs(ci$statistic))
  lab <- coef_labels(spec, design)

  # one variance matrix for every block that needs one: a term reported by
  # the quantities it is about carries them across by the delta method
  V <- catch_frozen(tryCatch(vcov(object, readable = FALSE, type = type,
                                  expected = expected, approx = approx, ...),
                             error = function(e) NULL))
  strc <- tryCatch(statmod_structural_table(object, level),
                   error = function(e) NULL)
  tables <- lapply(spec@distrib@params, function(p)
    summary_blocks(object, spec, design, p, ci, level, V, strc))
  names(tables) <- spec@distrib@params

  ll <- logLik.StatmodFit(object)
  df <- attr(ll, "df")
  notes <- c(character(0), frozen_msg,
             tryCatch(class_notes(spec, design), error = function(e) character(0)))
  # WHERE THE PARAMETERS ENDED UP, for a fit that did not converge. It
  # qualifies the fit rather than describing it, so it is a note and not a
  # line of print(): the measured case is a scale that ran to 1e-15 while
  # 380 of 400 coefficients survived, which is read once when something
  # looks wrong and never otherwise.
  if (!object@converged) {
    r <- tryCatch(fitted_ranges(object), error = function(e) "")
    if (nzchar(r)) notes <- c(notes, r)
  }
  # read HERE and not at the end, because a note below asks whether the point
  # is a maximum and that is the certificate's answer. It costs one outer
  # gradient and one solve, which is nothing beside a summary that already
  # inverts the penalized information, and a fold of cv() or a path point
  # that never prints itself pays nothing at all.
  cert <- tryCatch(statmod_certificate(object), error = function(e) NULL)

  # The count above reads the hyperparameters as though they were known,
  # and they were estimated from the same data. Adding what that costs is
  # off by default because it changes a number a reader may be comparing
  # with an earlier fit.
  corr <- 0
  if (isTRUE(correct)) {
    cc <- tryCatch(statmod_edf_correction(spec, object@coefficients,
                                          object@hyper, design,
                                          object@methods$outer),
                   error = function(e)
                     list(total = 0, per = numeric(0), n_hyper = 0L))
    corr <- cc$total
    df <- df + corr
    nh <- if (is.null(cc$n_hyper)) 0L else cc$n_hyper
    notes <- c(notes, if (corr > 0) sprintf(paste0(
      "The degrees of freedom carry %.3f for the hyperparameters, which ",
      "were\n  estimated rather than given. Without it the criteria are ",
      "too generous."), corr) else if (nh > 0L) paste0(
      "The correction for the estimated hyperparameters could not be ",
      "computed\n  here, so the degrees of freedom read them as though they ",
      "were given.\n  A shared hyperparameter is the case: the criterion's ",
      "curvature is not\n  defined over a group.") else paste0(
      "No hyperparameter here was estimated by a marginal criterion, so ",
      "there is\n  nothing for the correction to propagate and it is zero."))
  }
  if (any(lab$penalized)) {
    # WHICH criteria were in force, so a reader knows whether the number at
    # the head of a penalized block was chosen or given. The two kinds carry
    # different guarantees and are named separately.
    src <- unique(unlist(lapply(tables, function(bl)
      unlist(lapply(bl, function(b) b$table$source[b$table$role %in%
                                                     c("fixed", "estimated")]),
             use.names = FALSE)), use.names = FALSE))
    marg <- setdiff(intersect(src, c("reml", "ml")), NA)
    path <- setdiff(intersect(src, c("aic", "bic", "cv")), NA)
    # WHETHER THE STANDARD ERROR IS THERE, and not merely whether the
    # criterion was marginal. A hyperparameter the search left before
    # reaching a maximum has a curvature of the wrong sign in its own
    # direction and no variance follows from it, so the row carries the
    # estimate and nothing else -- and a note promising an interval beside
    # it would be false of the fit in front of the reader.
    mse <- unlist(lapply(tables, function(bl) unlist(lapply(bl, function(b) {
      r <- b$table$role %in% "estimated" & b$table$source %in% c("reml", "ml")
      if (!any(r)) return(logical(0))
      is.finite(b$table$se[r])
    }), use.names = FALSE)), use.names = FALSE)
    if (length(marg) && any(mse)) {
      notes <- c(notes, sprintf(paste0(
        "A hyperparameter marked %s was estimated by that criterion, and its",
        " standard\n  error and interval are read on the free scale its link ",
        "defines and mapped\n  back. Every coefficient beside it is still ",
        "conditional on the value reached."),
        paste(toupper(marg), collapse = " or ")))
    }
    if (length(marg) && any(!mse)) {
      notes <- c(notes, sprintf(paste0(
        "A hyperparameter marked %s with no standard error was estimated at ",
        "a point\n  where the criterion's curvature in its own direction is ",
        "not negative: the\n  search did not reach a maximum there, so no ",
        "variance follows from it. Its\n  estimate stands and every ",
        "coefficient beside it is conditional on it."),
        paste(toupper(marg), collapse = " or ")))
    }
    if (length(path)) {
      notes <- c(notes, sprintf(paste0(
        "A hyperparameter marked %s was chosen by a path over its own values,",
        " not held.\n  It is the argument of a minimum over a grid rather ",
        "than the root of a\n  derivative, so no standard error follows from ",
        "it; its uncertainty is a\n  resampling question."),
        paste(toupper(path), collapse = " or ")))
    }
    if ("fixed" %in% src) {
      notes <- c(notes, paste0(
        "A hyperparameter marked fixed is held at the value it was given, not",
        " estimated,\n  so it has no standard error and every interval beside",
        " it is conditional on it."))
    }
  }
  if (any(lab$kinked)) {
    nz <- sum(lab$kinked &
                unlist(object@coefficients[spec@distrib@params],
                       use.names = FALSE) != 0)
    notes <- c(notes, sprintf(paste0(
      "%d of %d coefficients under a kinked penalty survived. Their rows ",
      "are\n  conditional on that selection, and the ones at zero have no ",
      "variance at\n  all: at the kink there is no curvature to read."),
      nz, sum(lab$kinked)))
  }
  # ⚠️ WHETHER THE POINT IS A MAXIMUM is the CERTIFICATE's answer and not the
  # search's flag, and this note used to read the flag: a fit whose search
  # stopped on its own rule at a point located to 2e-11 printed "the fit did
  # not converge, so everything above is read at a point that is not a
  # maximum" directly under a certificate saying CONVERGED. The two are
  # different questions, so the note is emitted only where the question it
  # answers -- is this a maximum -- has actually been answered no.
  bad <- if (is.null(cert) || is.null(cert$state)) !object@converged else
    !cert$state %in% c("converged", "boundary")
  if (bad) {
    notes <- c(notes, paste0(
      "The point reported is not certified as a maximum, so every estimate ",
      "and\n  interval above is read where the surface is still moving."))
  }

  # A smoothed break-point term: the smoother and its width are part of the
  # model -- the transition's width, the bent-cable reading -- so they are
  # reported; and where a break-point carries a random development, a
  # smoother declaring a scale correction (the probit's convolution
  # identity) gets the corrected scale printed beside the apparent one.
  notes <- c(notes, tryCatch(smoothed_notes(spec, object),
                             error = function(e) character(0)))

  # A structural term contributes no columns, so its block is built from
  # what it reports rather than from a design; `strc` is computed above,
  # where summary_blocks() can be given it.
  if (!is.null(strc) && any(strc$held)) {
    notes <- c(notes, paste0(
      "A level marked held is carried by an intercept in the same equation ",
      "and is\n  not estimated: the two are exactly confounded, so only one ",
      "of them can be."))
  }

  StatmodSummary(
    call = object@call, distrib_name = spec@distrib@distrib_name,
    n_obs = spec@n_obs, tables = tables, edf = object@edf,
    links = vapply(spec@distrib@link_params,
                   function(g) g@link_name, character(1)),
    structural = strc,
    loglik = object@loglik, df = df,
    aic = -2 * object@loglik + 2 * df,
    bic = -2 * object@loglik + log(spec@n_obs) * df,
    converged = object@converged, elapsed = object@elapsed,
    level = level, type = type, notes = notes,
    certificate = cert)
}
S7::method(summary, StatmodFit) <- summary.StatmodFit


#' The Quantities a Penalty's Hyperparameters Are About
#'
#' @description
#' Replaces the coordinate rows of a penalty whose hyperparameters are a chart
#' with the quantities it declares through
#' [penalties7::penalty_readable()]: the standard deviations and
#' correlations of a correlated random effect, in place of the logarithms of
#' a
#' Cholesky diagonal and the entries below it.
#'
#' @details
#' The standard error is the delta method, and it composes two Jacobians: the
#' penalty's, which is in the parameter scale of its hyperparameters, and the
#' link's, the variance matrix being on the free scale the outer criterion was
#' maximized on. Each interval is built on the scale the quantity declares,
#' log for a standard deviation and Fisher's z for a correlation, and
#' mapped
#' back, so a standard deviation cannot be given a negative lower end and a
#' correlation cannot be given an interval that leaves \eqn{(-1, 1)}. That is
#' the rule every other interval in the toolkit follows.
#'
#' No test is printed, for the reason the coordinate rows print none: the null
#' a \eqn{z} of value over standard error reports on is that the quantity is
#' zero, which for a standard deviation is the edge of its range.
#'
#' @param rd The result of [penalties7::penalty_readable()].
#' @param th The penalty's hyperparameters, as fitted.
#' @param Vh The hyperparameter variance matrix, or `NULL`.
#' @param p The parameter the term sits in.
#' @param key The penalty's key.
#' @param level The confidence level.
#' @param role,src What the coordinate rows reported.
#' @param cols The column names of a summary block.
#'
#' @return A data frame of rows, in the shape of a summary block.
#'
#' @keywords internal
readable_hyper_rows <- function(rd, th, Vh, p, key, level, role, src, cols) {
  nm <- names(th)
  k <- length(rd$value)
  se <- rep(NA_real_, k)
  if (!is.null(Vh)) {
    j <- match(paste(p, key, nm, sep = "\r"), rownames(Vh))
    if (!anyNA(j)) {
      vb <- as.matrix(Vh[j, j, drop = FALSE])
      lk <- attr(attr(Vh, "idx"), "links")[j]
      d <- vapply(seq_along(nm), function(i) {
        eta <- linkfunctions7::linkfun(lk[[i]], th[[i]])
        linkfunctions7::dlinkinv(lk[[i]], eta)
      }, 0)
      jac <- rd$jacobian * rep(d, each = k)
      if (all(is.finite(vb)) && all(is.finite(jac))) {
        se <- sqrt(pmax(diag(jac %*% vb %*% t(jac)), 0))
      }
    }
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  v <- as.numeric(rd$value)
  tr <- as.character(rd$transform)
  lo <- hi <- rep(NA_real_, k)
  for (i in seq_len(k)) {
    if (!is.finite(se[i])) next
    ends <- switch(tr[i],
      log = exp(log(v[i]) + c(-1, 1) * z * se[i] / v[i]),
      atanh = tanh(atanh(v[i]) + c(-1, 1) * z * se[i] / (1 - v[i]^2)),
      v[i] + c(-1, 1) * z * se[i])
    ends <- sort(ends)
    lo[i] <- ends[[1L]]
    hi[i] <- ends[[2L]]
  }
  out <- data.frame(name = names(rd$value), estimate = v, se = se,
                    statistic = NA_real_, p_value = NA_real_,
                    lower = lo, upper = hi,
                    role = rep(role[[1L]], k), source = rep(src[[1L]], k),
                    stringsAsFactors = FALSE)
  stats::setNames(out, c(cols, "source"))
}

#' The Blocks of One Distribution Parameter
#'
#' @description
#' Groups a parameter's terms into the readings a summary prints: the
#' parametric terms together, and one block per penalized term.
#'
#' @param fit A [StatmodFit()].
#' @param spec The specification.
#' @param design The design.
#' @param p The distribution parameter.
#' @param ci The flat interval table, as [confint.StatmodFit()]
#'   returns it with the statistic and the p-value added.
#' @param level The confidence level the intervals are built at.
#' @param V The variance matrix over the stacked coefficients, or
#'   `NULL`. It is needed only by a term reported through
#'   [modelterms7::term_readable()], whose quantities are
#'   functions of several coefficients at once and whose standard errors
#'   are therefore the delta method and no diagonal entry.
#'
#' @return A list of block records, each with `kind`, `label`,
#'   `n_coef`, `edf`, `n_zero` and `table`, together
#'   with `head` and `components`: a term written in parameters
#'   of its own that develops one of them over covariates reports that
#'   parameter as a compartment of its own, carrying its hyperparameter and
#'   its sub-terms' rows, and `table` keeps only what is left.
#'
#' @param st The structural table, or `NULL`. A structural term has no
#'   design columns, so its block is built from what it reports, never from
#'   a block of the design, and its hyperparameter is reported there instead
#'   of in a block of its own carrying nothing else.
#'
#' @keywords internal
summary_blocks <- function(fit, spec, design, p, ci, level = 0.95,
                           V = NULL, st = NULL) {
  rows <- ci[ci$parameter == p, , drop = FALSE]
  cols <- c("name", "estimate", "se", "statistic", "p_value", "lower",
            "upper", "role")
  # `source` says what put the number there -- a criterion by name, or the
  # caller -- which `role` alone cannot, every criterion answering "estimated"
  all_cols <- c(cols, "source")
  empty <- stats::setNames(
    data.frame(character(0), numeric(0), numeric(0), numeric(0), numeric(0),
               numeric(0), numeric(0), character(0), character(0),
               stringsAsFactors = FALSE), all_cols)

  coef_rows <- function(nm) {
    r <- rows[rows$term == nm, , drop = FALSE]
    if (!nrow(r)) return(empty)
    out <- data.frame(name = r$coefficient, estimate = r$estimate, se = r$se,
                      statistic = r$statistic, p_value = r$p_value,
                      lower = r$lower, upper = r$upper, role = "coefficient",
                      source = "", stringsAsFactors = FALSE)
    stats::setNames(out, all_cols)
  }
  # WHICH hyperparameter was estimated and by what. A marginal criterion --
  # reml(), ml() -- maximizes a twice differentiable function of it, so it
  # carries a standard error and an interval, both read on the free scale its
  # link defines and mapped back (statmod_hyper_vcov). A path -- aic(), bic(),
  # cv() over a kinked penalty -- chooses the argument of a minimum over a
  # grid, so there is no curvature to read and the row reports the value and
  # the criterion that chose it. Only a value the caller set is held.
  outer_ran <- !is.null(fit@methods$outer)
  Vh <- if (outer_ran) tryCatch(
    statmod_hyper_vcov(spec, design, fit@coefficients, fit@hyper,
                       fit@methods$outer), error = function(e) NULL) else NULL
  spc <- fit@methods$sparse_criterion
  spc_keys <- fit@methods$sparse_hyper
  if (is.null(spc_keys)) spc_keys <- character(0)
  # WHAT WAS HELD is the terms' answer and nobody else's: a hyperparameter
  # a term fixed is fixed whatever criterion ran beside it.
  held <- tryCatch(statmod_held(spec, design), error = function(e) character(0))
  # A term may carry more than one penalty, each filed under a key of its
  # own, so the rows of a term are those of every key belonging to it. Where
  # there are several the hyperparameter is named for the penalty as well:
  # two lambdas in one block, one on the slope changes and one on the jumps,
  # are not the same number and cannot appear under the same name.
  hyper_parts <- function(nm) {
    ent <- modelterms7::term_penalties(spec@terms[[p]][[nm]])
    des <- statmod_design(spec)
    # A LABELLED term declares no penalty of its own: the covariance class
    # carries it. The class's hyperparameters are reported ONCE, under its
    # first member, and a note names the equations the block spans -- a
    # class is one penalty and printing it under each member would say
    # there are several.
    cl <- class_unit_of(spec, des, p, nm)
    if (!is.null(cl)) {
      first <- cl$pieces[[1L]]
      # under the FIRST member only, so one penalty prints once
      if (identical(first$param, p) && identical(first$term, nm)) {
        # ADDED to the term's own entries, never in place of them: a term may
        # carry a penalty of its own beside a labelled sub-term, and replacing
        # them would drop it from the page
        ent <- c(ent, list(list(name = "", penalty = cl$penalty,
                                index = cl$index, class_key = cl$key)))
      }
    }
    if (!length(ent)) return(list(rows = list(), index = list()))
    own <- length(ent) - as.integer(!is.null(cl) &&
                                    identical(cl$pieces[[1L]]$param, p) &&
                                    identical(cl$pieces[[1L]]$term, nm))
    out <- lapply(ent, function(e) {
      key <- if (is.null(e$class_key)) {
        statmod_entry_key(nm, utils::head(ent, own), e)
      } else e$class_key
      th <- fit@hyper[[p]][[key]]
      if (is.null(th) || !length(th)) return(empty)
      u <- statmod_unit(spec, des, p, key)
      if (is.null(u)) return(empty)
      lab <- if (length(ent) > 1L && nzchar(e$name)) {
        paste(e$name, names(th), sep = ".")
      } else {
        names(th)
      }
      marginal <- outer_ran && !penalty_has_kink(u$penalty)
      hk <- paste(p, key, names(th), sep = "\r")
      # PER HYPERPARAMETER, and the test has to be the vector one: `ifelse`
      # returns a result the length of its TEST, so a scalar `marginal`
      # returned one answer and recycled it over a penalty carrying two
      # hyperparameters -- the elastic net's alpha, which no path varies,
      # was reported as chosen by the criterion that had chosen its lambda
      chosen <- (rep(marginal, length(hk)) | hk %in% spc_keys) &
        !(hk %in% held)
      role <- ifelse(chosen, "estimated", "fixed")
      by <- if (marginal) fit@methods$outer@kind else
        if (is.null(spc)) "estimated" else spc@kind
      src <- ifelse(chosen, by, "fixed")
      r <- data.frame(name = lab, estimate = as.numeric(th),
                      se = NA_real_, statistic = NA_real_, p_value = NA_real_,
                      lower = NA_real_, upper = NA_real_, role = role,
                      source = src, stringsAsFactors = FALSE)
      # the interval is built where the criterion was maximized and mapped
      # back, as every other interval in the toolkit is, so a positive
      # hyperparameter cannot be given a negative lower end; the standard
      # error printed beside it is the delta method onto its own scale. No
      # test accompanies it: the null a z of value/se reports on is that the
      # hyperparameter is zero, which for a smoothing parameter is the edge
      # of its range and not an interior hypothesis.
      if (!is.null(Vh) && marginal) {
        z <- stats::qnorm(1 - (1 - level) / 2)
        lk <- attr(Vh, "idx")
        for (i in seq_along(th)) {
          k <- paste(p, key, names(th)[[i]], sep = "\r")
          j <- match(k, rownames(Vh))
          if (is.na(j)) next
          link <- attr(lk, "links")[[j]]
          # a coordinate the criterion's curvature cannot produce a
          # variance for is NA here and keeps the row's empty se and
          # interval; sort() drops an NA, so reaching the ends below with
          # one would ask for an element that is not there
          if (!is.finite(Vh[j, j])) next
          se_eta <- sqrt(Vh[j, j])
          eta <- linkfunctions7::linkfun(link, r$estimate[[i]])
          r$se[[i]] <- abs(linkfunctions7::dlinkinv(link, eta)) * se_eta
          ends <- sort(c(linkfunctions7::linkinv(link, eta - z * se_eta),
                         linkfunctions7::linkinv(link, eta + z * se_eta)))
          r$lower[[i]] <- ends[[1L]]
          r$upper[[i]] <- ends[[2L]]
        }
      }
      # WHERE THE HYPERPARAMETERS ARE A CHART and not the quantities, the
      # quantities are reported instead: a correlated random effect is about
      # the standard deviations and correlations of its effects, and the
      # logarithms of a Cholesky diagonal are the coordinates that produce
      # them. The penalty declares which case it is; every other branch
      # answers NULL and its rows stand.
      rd <- tryCatch(penalties7::penalty_readable(u$penalty, th),
                     error = function(e) NULL)
      if (!is.null(rd) && length(rd$value)) {
        rr <- readable_hyper_rows(rd, th, if (marginal) Vh else NULL, p, key,
                                  level, role, src, cols)
        # A hyperparameter the readable block does not DESCRIBE keeps its own
        # row: a multivariate Student t is about the standard deviations and
        # the correlations of its scale matrix, and its degrees of freedom are
        # none of those. The question is asked of the Jacobian -- a column that
        # is zero throughout is a coordinate no quantity depends on -- so a
        # family that declares more later is covered without an edit.
        keep <- apply(rd$jacobian, 2L, function(z) all(z == 0))
        if (any(keep)) {
          rr <- rbind(rr, stats::setNames(r[keep, , drop = FALSE],
                                          c(cols, "source")))
        }
        return(rr)
      }
      stats::setNames(r, c(cols, "source"))
    })
    # PER ENTRY and not stacked, because an entry belongs where its
    # coefficients do: a penalty covering a developed parameter's columns is
    # that parameter's hyperparameter and is reported inside its
    # compartment, while one covering the term's own columns stays with the
    # term. Which is which is read off the columns rather than off the
    # entry's name, so a term that names its entries differently is covered
    # without an edit.
    list(rows = out, index = lapply(ent, `[[`, "index"))
  }
  hyper_rows <- function(nm) {
    hp <- hyper_parts(nm)
    if (!length(hp$rows)) return(empty)
    do.call(rbind, hp$rows)
  }
  # A term whose block is a working linearization is reported by what it is
  # ABOUT and not by the coefficients it is fitted through: a break-point
  # term's auxiliary pair carries the position as -g/delta, which is a
  # number no reader wants to compute. term_readable() gives the quantities
  # and the Jacobian from the coefficients, so the variance comes across by
  # the delta method -- the same route segmented reports a break-point's
  # standard error by. A term that answers NULL, which is what a developed
  # one does, falls back to its coefficients.
  readable_rows <- function(nm) {
    term <- spec@terms[[p]][[nm]]
    idx <- design[[p]]$blocks[[nm]]
    rd <- tryCatch(modelterms7::term_readable(term, fit@coefficients[[p]][idx]),
                   error = function(e) NULL)
    if (is.null(rd) || !length(rd$name)) return(NULL)
    # V is indexed by NAME rather than by position: coef_labels() skips a
    # parameter with no coefficients, so a stacked offset computed here
    # would not be the one the matrix was built with
    key <- paste(p, design[[p]]$coef_names[idx], sep = ":")
    se <- rep(NA_real_, length(rd$name))
    if (!is.null(V) && all(key %in% rownames(V))) {
      Vb <- as.matrix(V[key, key, drop = FALSE])
      if (all(is.finite(Vb))) {
        se <- sqrt(pmax(diag(rd$jacobian %*% Vb %*% t(rd$jacobian)), 0))
      }
    }
    z <- stats::qnorm(1 - (1 - level) / 2)
    st <- rd$value / se
    # A break-point gets an estimate and an interval and NO test. The null
    # a z of value/se would report on is that the position is zero, which
    # is not a hypothesis anyone holds; the one a reader wants is that
    # there is no break-point at all, and under it the position is a
    # nuisance parameter that vanishes, so the classical p-value is wrong
    # by a factor of three to five (Davies' problem). segmented prints the
    # estimate and the standard error alone for the same reason.
    pos <- grepl("^psi[0-9]*$", rd$name)
    st[pos] <- NA_real_
    out <- data.frame(name = rd$name, estimate = rd$value, se = se,
                      statistic = st, p_value = 2 * stats::pnorm(-abs(st)),
                      lower = rd$value - z * se, upper = rd$value + z * se,
                      role = "coefficient", source = "",
                      stringsAsFactors = FALSE)
    stats::setNames(out, all_cols)
  }
  term_edf <- function(nm) {
    if (is.null(fit@edf)) return(NA_real_)
    e <- fit@edf
    v <- e$edf[e$parameter == p & e$term == nm]
    if (length(v)) v[1L] else NA_real_
  }

  # WHAT A SUB-TERM IS CALLED where a compartment header names it. The
  # term's own label is used where it has one; a parametric development is
  # named for what it carries, an intercept alone being the population
  # value of the parameter and anything else a set of covariates.
  sub_label <- function(s) {
    lb <- tryCatch(s@label, error = function(e) "")
    if (length(lb) == 1L && nzchar(lb)) return(lb)
    nms <- tryCatch(modelterms7::term_coef_names(s),
                    error = function(e) character(0))
    if (length(nms) == 1L && endsWith(nms, "(Intercept)")) "intercept" else
      "covariates"
  }
  # WHAT A HYPERPARAMETER IS, rather than which coordinate carries it. A
  # compartment's `sigma` sits under a term of a model whose distribution
  # has a `sigma` of its own, and the two are different quantities; a
  # gaussian prior's is the scale of the effects it shrinks and a quadratic
  # penalty's `lambda` is a precision. Any other name is left as it stands.
  hyper_label <- function(nm, pen) {
    if (identical(nm, "sigma") &&
        grepl("gaussian", pen@penalty_name, fixed = TRUE)) {
      return("effect sd")
    }
    if (identical(nm, "lambda") &&
        isTRUE(tryCatch(penalties7::is_quadratic(pen),
                        error = function(e) FALSE))) {
      return("precision")
    }
    nm
  }
  # WHAT A SUB-TERM CONTRIBUTES to its compartment: its rows filtered the way
  # a block of that kind is filtered at the top level, so a smooth shows its
  # unpenalized part, a kinked penalty shows what it kept, and a random
  # effect shows no coefficients at all. What a reader wants of a set of
  # predictions is how many there are and how far they spread, which is one
  # line, where the predictions themselves are a column of numbers nobody
  # reads and are still in `coef()`.
  #
  # `rows_at` is how the rows of a set of the term's own coefficients are
  # obtained. An additive term reads its design block's rows by column; a
  # structural one has no design and reads what it REPORTS by position in
  # its parameter vector. Everything below is common to the two.
  sub_block <- function(s, ix, rows_at) {
    sk <- term_block_kind(s)
    rr <- rows_at(ix)
    nms <- tryCatch(modelterms7::term_coef_names(s),
                    error = function(e) character(0))
    if (length(nms) == nrow(rr)) rr$name <- nms
    keep <- switch(sk,
      smooth = smooth_linear_cols(s, nrow(rr)),
      selection = rr$estimate != 0,
      rep(TRUE, nrow(rr)))
    if (identical(sk, "random")) keep <- rep(FALSE, nrow(rr))
    lines <- character(0)
    if (identical(sk, "random") && nrow(rr)) {
      lines <- sprintf("%d predictions, sd %s, range %s to %s", nrow(rr),
                       format(signif(stats::sd(rr$estimate), 3L)),
                       format(signif(min(rr$estimate), 3L)),
                       format(signif(max(rr$estimate), 3L)))
      names(lines) <- "effects"
    }
    list(kind = sk, table = rr[keep, , drop = FALSE], lines = lines)
  }
  # ONE COMPARTMENT PER DEVELOPED PARAMETER, carrying its own hyperparameter
  # first and then each sub-term's rows in the order the block binds them.
  compartments <- function(term, rows_at, dev, hp) {
    ent <- modelterms7::term_penalties(term)
    lapply(dev, function(cp) {
      mine <- vapply(hp$index, function(ii) {
        length(ii) > 0L && all(ii %in% cp$index)
      }, logical(1))
      hr <- hp$rows[mine]
      pens <- lapply(ent[mine], function(e) e$penalty)
      for (i in seq_along(hr)) {
        if (!nrow(hr[[i]])) next
        hr[[i]]$name <- vapply(hr[[i]]$name, hyper_label, character(1),
                               pen = pens[[i]], USE.NAMES = FALSE)
      }
      parts <- lapply(seq_along(cp$subs), function(i) {
        sub_block(cp$subs[[i]], cp$sub_index[[i]], rows_at)
      })
      tb <- do.call(rbind, c(hr, lapply(parts, function(z) z$table)))
      ln <- unlist(lapply(parts, function(z) z$lines))
      list(name = cp$name,
           header = sprintf("%s  ~ %s", cp$name,
                            paste(vapply(cp$subs, sub_label, character(1)),
                                  collapse = " + ")),
           table = if (is.null(tb)) empty else tb,
           lines = if (is.null(ln)) character(0) else ln,
           n_coef = length(cp$index))
    })
  }
  # THE DEVELOPED PARAMETERS READ AT A GLANCE, one line each: the population
  # value of the development, and what develops it. A parameter that is a
  # number of its own is one row of the table below and needs no line above
  # it; a developed one is spread over a compartment where its population
  # value is labeled by the development's intercept, so this is the only
  # place the parameter's own name appears beside a number.
  head_rows <- function(rows_at, dev) {
    out <- lapply(dev, function(cp) {
      r <- NULL
      note <- paste("~", paste(vapply(cp$subs, sub_label, character(1)),
                               collapse = " + "))
      for (i in seq_along(cp$subs)) {
        j <- which(tryCatch(modelterms7::term_coef_names(cp$subs[[i]]),
                            error = function(e) character(0)) ==
                   "(Intercept)")
        if (length(j) == 1L) {
          r <- rows_at(cp$sub_index[[i]][[j]])
          break
        }
      }
      if (!is.null(r) && !nrow(r)) r <- NULL
      data.frame(name = cp$name,
                 estimate = if (is.null(r)) NA_real_ else r$estimate[[1L]],
                 lower = if (is.null(r)) NA_real_ else r$lower[[1L]],
                 upper = if (is.null(r)) NA_real_ else r$upper[[1L]],
                 note = note, stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
  }
  # WHICH COLUMNS ARE THE TERM'S OWN, and which hyperparameters go with
  # them: everything a developed parameter does not claim.
  own_hyper <- function(hp, dev) {
    hp$rows[!vapply(hp$index, function(ii) {
      any(vapply(dev, function(cp) length(ii) > 0L && all(ii %in% cp$index),
                 logical(1)))
    }, logical(1))]
  }
  developed <- function(term, n_own) {
    cp <- tryCatch(modelterms7::term_components(term),
                   error = function(e) list())
    if (!length(cp)) return(list())
    idx <- unlist(lapply(cp, function(z) z$index), use.names = FALSE)
    if (!length(idx) || max(idx) > n_own) return(list())
    cp
  }
  # A STRUCTURAL TERM contributes no design columns, so its block is built
  # from what it REPORTS -- the quantities term_readable() names, with the
  # variance of the joint information behind them -- and its hyperparameter
  # sits with the coordinates it shrinks rather than in a block of its own
  # carrying nothing else.
  structural_rows <- function(r) {
    stats::setNames(data.frame(
      name = ifelse(r$held, paste0(r$name, " (held)"), r$name),
      estimate = r$estimate, se = r$se, statistic = NA_real_,
      p_value = NA_real_, lower = r$lower, upper = r$upper,
      role = "coefficient", source = "", stringsAsFactors = FALSE), all_cols)
  }
  structural_block <- function(nm, term) {
    r <- st[st$parameter == p & st$term == nm, , drop = FALSE]
    hp <- hyper_parts(nm)
    cp_all <- developed(term, length(modelterms7::term_params(term)))
    dev <- cp_all[vapply(cp_all, function(z) length(z$subs) > 0L, logical(1))]
    rows_at <- function(ix) {
      structural_rows(r[!is.na(r$position) & r$position %in% ix, ,
                        drop = FALSE])
    }
    own <- structural_rows(r[!(r$component %in% names(dev)), , drop = FALSE])
    tb <- do.call(rbind, c(own_hyper(hp, dev), list(own)))
    list(kind = "structural", label = block_label("structural"), term = nm,
         n_coef = nrow(r), edf = term_edf(nm), n_zero = 0L,
         table = if (is.null(tb)) empty else tb,
         head = if (length(dev)) head_rows(rows_at, dev) else NULL,
         components = compartments(term, rows_at, dev, hp))
  }

  blocks <- list()
  para <- character(0)
  for (nm in names(spec@terms[[p]])) {
    if (identical(term_block_kind(spec@terms[[p]][[nm]]), "parametric")) {
      para <- c(para, nm)
    }
  }
  if (length(para)) {
    tb <- do.call(rbind, lapply(para, coef_rows))
    blocks[[length(blocks) + 1L]] <- list(
      kind = "parametric", label = "Parametric terms", term = NA_character_,
      n_coef = nrow(tb), edf = sum(vapply(para, term_edf, numeric(1))),
      n_zero = 0L, table = tb, head = NULL, components = list())
  }

  for (nm in names(spec@terms[[p]])) {
    term <- spec@terms[[p]][[nm]]
    kind <- term_block_kind(term)
    if (identical(kind, "parametric")) next
    k <- length(design[[p]]$blocks[[nm]])
    if (!k) {
      # no design columns: a structural term, reported by what it says it
      # reports. Where the structural table could not be built there is
      # nothing to show and the term is passed over rather than given an
      # empty block.
      if (is.null(st) || !any(st$parameter == p & st$term == nm)) next
      blocks[[length(blocks) + 1L]] <- structural_block(nm, term)
      next
    }
    cr <- coef_rows(nm)
    hp <- hyper_parts(nm)
    rows_at <- function(ix) cr[ix, , drop = FALSE]
    # a component whose parameter is a number of its own needs no
    # compartment: it is one coefficient and belongs in the term's own table
    cp_all <- if (nrow(cr) == k) developed(term, k) else list()
    dev <- cp_all[vapply(cp_all, function(z) length(z$subs) > 0L, logical(1))]
    if (length(dev)) {
      own <- setdiff(seq_len(k),
                     unlist(lapply(dev, function(z) z$index),
                            use.names = FALSE))
      # the own table reports the same quantities a term without a
      # development reports: the position of a break-point and not the pair
      # of working coefficients it is read off
      rr <- if (identical(kind, "breakpoint")) readable_rows(nm) else NULL
      tb <- do.call(rbind, c(own_hyper(hp, dev), list(
        if (is.null(rr)) cr[own, , drop = FALSE] else rr)))
      blocks[[length(blocks) + 1L]] <- list(
        kind = kind, label = block_label(kind), term = nm, n_coef = k,
        edf = term_edf(nm), n_zero = 0L,
        table = if (is.null(tb)) empty else tb,
        head = head_rows(rows_at, dev),
        components = compartments(term, rows_at, dev, hp))
      next
    }
    keep <- switch(kind,
      smooth = smooth_linear_cols(term, nrow(cr)),
      # a coefficient a kinked penalty set to zero is counted, not listed
      selection = cr$estimate != 0,
      rep(TRUE, nrow(cr)))
    if (identical(kind, "random")) keep <- rep(FALSE, nrow(cr))
    body <- if (identical(kind, "breakpoint")) {
      rr <- readable_rows(nm)
      if (is.null(rr)) cr else rr
    } else {
      cr[keep, , drop = FALSE]
    }
    # the hyperparameters come FIRST in every penalized block: they govern
    # everything below them, and a table that opens with a hundred selected
    # coefficients buries the one number that produced that selection
    tb <- rbind(hyper_rows(nm), body)
    blocks[[length(blocks) + 1L]] <- list(
      kind = kind, label = block_label(kind), term = nm, n_coef = k,
      edf = term_edf(nm),
      n_zero = if (identical(kind, "selection")) sum(cr$estimate == 0) else 0L,
      table = tb, head = NULL, components = list())
  }
  blocks
}


#' The Covariance Class One Term Belongs To
#'
#' @description
#' The penalized unit of the covariance class that covers a term's
#' coefficients, or `NULL` where the term is not labelled.
#'
#' @details
#' A labelled term declares no penalty of its own, the class carrying it, so a
#' reader asking a term for its hyperparameters has to ask the class instead.
#' The lookup is by the pair of the equation and the term's name, which is what
#' a class's pieces record.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param param The distribution parameter the term sits in.
#' @param nm The term's name in the formula.
#'
#' @return One unit, as [statmod_penalized()] returns it, or `NULL`.
#'
#' @seealso [statmod_classes()] for the classes, [summary_blocks()] for the
#'   reader that needs this.
#'
#' @keywords internal
class_unit_of <- function(spec, design, param, nm) {
  for (u in statmod_penalized(spec, design)) {
    if (is.null(u$pieces)) next
    hit <- vapply(u$pieces, function(z)
      identical(z$param, param) && identical(z$term, nm), TRUE)
    if (any(hit)) return(u)
  }
  NULL
}


#' What a Summary Says About the Covariance Classes
#'
#' @description
#' One note per class spanning more than one term, naming the label, the
#' grouping and the terms whose coefficients share the block.
#'
#' @details
#' A class's hyperparameters are printed once, under its first member, so
#' without the note a reader sees a covariance of four coordinates under a
#' term carrying two columns and nothing saying where the other two came from.
#'
#' A class of one member gets no note: there is nothing shared to report, and
#' its block is the random effect it would have been without a label.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#'
#' @return A character vector, possibly empty.
#'
#' @seealso [summary.StatmodFit()], which collects it.
#'
#' @keywords internal
class_notes <- function(spec, design) {
  out <- character(0)
  for (u in tryCatch(statmod_penalized(spec, design),
                     error = function(e) list())) {
    if (is.null(u$pieces) || length(u$pieces) < 2L) next
    who <- vapply(u$pieces, function(z) sprintf("'%s' in '%s'", z$term, z$param),
                  "")
    out <- c(out, sprintf(paste0(
      "The covariance block '%s' is shared by %s: one %d-variate prior over ",
      "the effects of one level of %s, with the hyperparameters reported ",
      "under the first of them."),
      u$key, paste(who, collapse = " and "), u$class$dim, u$class$group))
  }
  out
}


#' The Heading a Block of Each Kind Is Printed Under
#'
#' @description
#' The name a summary gives a block, from the kind
#' [term_block_kind()] answered with.
#'
#' @param kind The kind of the block.
#'
#' @return A single string.
#'
#' @keywords internal
block_label <- function(kind) {
  switch(kind, smooth = "Smooth", random = "Random effect",
         selection = "Selection", breakpoint = "Break-points",
         parametric = "Parametric terms", structural = "Structural",
         "Penalized")
}


#' The Notes a Smoothed Break-Point Term Adds to a Summary
#'
#' @description
#' One note per smoothed term, naming the smoother and the width the build
#' resolved, the width of the transition being part of the model, and, where
#' a break-point carries a random development
#' under a Gaussian precision and the smoother declares a scale correction,
#' the corrected scale beside the apparent one.
#'
#' @details
#' The correction is the smoother's own: the probit satisfies the exact
#' convolution identity \eqn{\tau^2_{\mathrm{apparent}} = \tau^2 + h^2}, so
#' the corrected scale is \eqn{\sqrt{\tau^2 - h^2}}; a smoother declaring
#' none (the hyperbolic, the quintic) gets the apparent scale alone, which
#' the random effect's own block already reports. The apparent scale is
#' read off the ridge precision as \eqn{1/\sqrt{\lambda}}, which is only a
#' scale where the penalty is the quadratic branch with that
#' hyperparameter; any other development is left without the note rather
#' than given a number of the wrong meaning.
#'
#' @param spec The fitted specification, whose terms are the ones the fit
#'   left.
#' @param object The fit.
#'
#' @return A character vector, possibly empty.
#'
#' @seealso [penalties7::abs_smoother()]
#'
#' @keywords internal
smoothed_notes <- function(spec, object) {
  out <- character(0)
  for (p in names(spec@terms)) {
    for (nm in names(spec@terms[[p]])) {
      tm <- spec@terms[[p]][[nm]]
      if (!S7::S7_inherits(tm, modelterms7::SegTerm)) next
      bp <- tryCatch(tm@blueprint, error = function(e) NULL)
      smx <- if (is.list(bp)) bp$smooth else NULL
      if (is.null(smx)) next
      sm <- smx$sm
      out <- c(out, sprintf(paste0(
        "'%s' in '%s' is smoothed (%s, %s = %s%s): its break-points are ",
        "ordinary\n  parameters and that is the width of the transition, ",
        "so its rows are the\n  smoothed model's."),
        nm, p, sm@smoother_name, sm@width_name,
        format(smx$width, digits = 3),
        if (!is.null(smx$w_group)) ", per group" else ""))
      if (is.null(sm@tau_correction)) next
      ent <- tryCatch(modelterms7::term_penalties(tm),
                      error = function(e) list())
      for (e in ent) {
        if (!grepl("^psi", e$name)) next
        pen <- e$penalty
        key <- statmod_entry_key(nm, ent, e)
        th <- tryCatch(object@hyper[[p]][[key]], error = function(err) NULL)
        if (is.null(th)) next
        # The apparent scale of the break-point deviations. A gaussian
        # prior written by its scale carries it as sigma, which is what
        # random() declares; the quadratic ridge carries the precision, so
        # the scale is 1/sqrt(lambda). Any other prior is left without the
        # note rather than given a number of the wrong meaning: the
        # convolution identity composes GAUSSIAN variances.
        tau <- NULL
        if ("sigma" %in% pen@params &&
            grepl("gaussian", pen@penalty_name, fixed = TRUE)) {
          tau <- suppressWarnings(as.numeric(th[["sigma"]]))
        } else if (identical(pen@params, "lambda") &&
                   isTRUE(tryCatch(penalties7::is_quadratic(pen),
                                   error = function(err) FALSE))) {
          lam <- suppressWarnings(as.numeric(th[["lambda"]]))
          if (length(lam) == 1L && is.finite(lam) && lam > 0) {
            tau <- 1 / sqrt(lam)
          }
        }
        if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) ||
            tau <= 0) next
        tauc <- sm@tau_correction(tau, smx$width)
        out <- c(out, sprintf(paste0(
          "The scale of the random break-points of '%s' composes with the ",
          "smoothing\n  width: apparent tau %s, corrected ",
          "sqrt(tau^2 - %s^2) = %s."),
          nm, format(tau, digits = 4), sm@width_name,
          format(tauc, digits = 4)))
      }
    }
  }
  out
}

#' @title Print a Model Summary
#' @name print.StatmodSummary
#' @description
#' The call, then each distribution parameter's blocks, then the degrees of
#' freedom, the criteria and the notes.
#' @param x A [StatmodSummary()].
#' @param digits Significant digits in the tables.
#' @param notes Whether to print the qualifications the numbers carry.
#'   `FALSE` by default, when the foot says how many there are: they
#'   state conventions, never facts of the fit, so they read the same
#'   under every model. They are on the summary's `notes` property
#'   either way.
#' @param n How many coefficient rows a block shows, `Inf` or `NA` for all
#'   of them. Defaults to the option `statmodels7.summary_rows`, and to 10
#'   where that is unset. A hyperparameter row is shown whatever `n` is: it
#'   governs the coefficients under it, and every one of them is
#'   conditional on the value it reached. A block short enough to fit in
#'   twelve rows is never abridged, so raising `n` changes nothing for a
#'   parametric block of ordinary size.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @seealso [summary.StatmodFit()]
#' @keywords internal
print.StatmodSummary <- function(x, digits = 4L, notes = FALSE,
                                 n = NULL, ...) {
  cat("A statmod fit\n\n")
  cat("Call:  ", paste(deparse(x@call), collapse = "\n        "), "\n\n",
      sep = "")
  cat("Distribution: ", x@distrib_name, "     Observations: ", x@n_obs,
      "\n", sep = "")

  for (p in names(x@tables)) {
    lk <- if (p %in% names(x@links)) x@links[[p]] else ""
    cat("\n", strrep("=", 3L), " ", p,
        if (nzchar(lk)) paste0("   [", lk, " link]") else "", "\n", sep = "")
    blocks <- x@tables[[p]]
    if (!length(blocks)) {
      cat("  (no coefficients)\n")
      next
    }
    for (b in blocks) print_block(b, digits, n)
  }

  cat(sprintf("\n%.0f%% intervals, %s variance\n", 100 * x@level, x@type))
  # WHICH log-likelihood and WHICH degrees of freedom, because the pairing
  # is what makes the criterion mean anything and the two conventions in
  # common use are not comparable. This one is conditional: the likelihood
  # is read at the fitted coefficients, a penalized coefficient among them,
  # and the count is the effective degrees of freedom. A mixed-model package
  # reporting a MARGINAL likelihood integrates its random effects out and
  # counts variance components instead, and its AIC is a different number
  # answering a different question (Vaida and Blanchard, 2005).
  cat(sprintf(paste0("conditional log-likelihood %.6f    effective df %.2f",
                     "\ncAIC %.3f    cBIC %.3f\n"),
              x@loglik, x@df, x@aic, x@bic))
  # ⚠️ THE TWO LINES BELOW ANSWER DIFFERENT QUESTIONS and used to read as a
  # contradiction: a fit could print DID NOT CONVERGE in capitals and
  # certificate: CONVERGED two lines under it, with nothing on screen saying
  # why both could be true. The first is a property of the SEARCH -- whether
  # it met its own stopping rule -- and the second of the POINT it stopped
  # at, which is the question a reader has. So the point gets the capitals
  # and the search says what it is.
  cat(sprintf("fitted in %s   search: %s\n", format_duration(x@elapsed),
              if (x@converged) "converged" else "DID NOT CONVERGE"))
  if (!is.null(x@certificate)) {
    ct <- x@certificate
    cat(sprintf("certificate: %s", toupper(ct$state)))
    if (is.finite(ct$gradient)) {
      cat(sprintf("   outer gradient %.3g", ct$gradient))
    }
    if (is.finite(ct$mode_error)) {
      cat(sprintf("   %.3g above the mode", ct$mode_error))
    }
    cat("\n")
    if (length(ct$boundary)) {
      cat("  at a boundary: ", paste(ct$boundary, collapse = ", "), "\n",
          sep = "")
    }
    for (r in ct$reason) {
      cat(paste0("  ", strwrap(r, width = 74L)), sep = "\n")
      cat("\n")
    }
  }
  # THE NOTES ARE COUNTED AND NOT PRINTED. Almost every one of them states a
  # CONVENTION rather than a fact of this fit -- why a hyperparameter's test
  # columns are empty, what a coefficient beside an estimated one is
  # conditional on -- so it reads the same under every model and fills the
  # foot of every summary with a paragraph nobody reads twice. Nothing is
  # dropped in silence: the count says they are there and how to see them,
  # and `summary(fit)@notes` has carried them all along.
  if (length(x@notes)) {
    if (isTRUE(notes)) {
      cat("\n")
      for (n in x@notes) {
        cat(paste0("  ", strwrap(n, width = 74L)), sep = "\n")
        cat("\n")
      }
    } else {
      cat(sprintf("%d note%s: print(summary(fit), notes = TRUE)\n",
                  length(x@notes), if (length(x@notes) > 1L) "s" else ""))
    }
  }
  invisible(x)
}
S7::method(print, StatmodSummary) <- print.StatmodSummary


#' The Rows of a Block Formatted for Printing
#'
#' @description
#' The six numeric columns of a summary table rendered as strings, with the
#' cells a hyperparameter row has no number for left empty and the name of
#' each such row carrying what put the value there.
#'
#' @details
#' A hyperparameter row prints numbers where there are any: one estimated by
#' a marginal criterion carries a standard error and an interval. Where there
#' is none the columns are blank. What put the value there goes in the name,
#' on every hyperparameter row, never only on the ones with nothing else
#' in them: written into the column where a standard error would have been it
#' marked a held or path-chosen row and never a REML one, whose column is
#' occupied, so the note at the foot spoke of a mark that was never printed.
#'
#' The whole block is formatted in one call, its own rows and every
#' compartment's together, so that the widths a column is padded to are the
#' same throughout and the compartments line up under the table they sit
#' beneath.
#'
#' @param tb A summary table.
#' @param digits Significant digits.
#'
#' @return A list with `cells`, a character matrix of six columns, and
#'   `name`, the row labels.
#'
#' @keywords internal
format_block_cells <- function(tb, digits = 4L) {
  hyp <- tb$role %in% c("fixed", "estimated")
  fixed <- hyp & !is.finite(tb$se)
  num <- function(v) ifelse(is.na(v), "", format(signif(v, digits)))
  out <- cbind(
    estimate = format(signif(tb$estimate, digits)),
    se = num(tb$se),
    z = num(tb$statistic),
    p = ifelse(is.na(tb$p_value), "",
               format.pval(tb$p_value, digits = digits, eps = 1e-16)),
    lower = num(tb$lower),
    upper = num(tb$upper))
  out[fixed, c("se", "z", "p", "lower", "upper")] <- ""
  # an estimated one that does carry an interval still has no test: the null
  # a z would report on is that the hyperparameter is zero, the edge of its
  # range rather than an interior hypothesis
  out[hyp & !fixed, c("z", "p")] <- ""
  src <- if (is.null(tb$source)) tb$role else
    ifelse(nzchar(tb$source), tb$source, tb$role)
  nm <- tb$name
  nm[hyp] <- sprintf("%s [%s]", nm[hyp], src[hyp])
  list(cells = out, name = nm)
}

#' Which Rows of a Table a Summary Prints
#'
#' @description
#' Every row of a short table, and of a long one the hyperparameters
#' together with the first few coefficients, the rest reported as a count.
#'
#' @details
#' A block of a few coefficients is printed whole: a threshold that cut it
#' would hide the very numbers a reader opened the summary for. A block of
#' many is a column of numbers nobody reads to the end, and what is dropped
#' is still in [coef()]. The hyperparameters are never dropped,
#' whatever the length: they govern every coefficient under them.
#'
#' @param tb A summary table.
#' @param cap The length above which a table is cut.
#' @param show How many coefficient rows a cut table keeps.
#'
#' @return An integer vector of row positions.
#'
#' @keywords internal
block_rows_shown <- function(tb, n = NULL, cap = 12L, show = 10L) {
  if (is.null(n)) n <- getOption("statmodels7.summary_rows", show)
  if (is.na(n) || !is.finite(n)) return(seq_len(nrow(tb)))
  n <- max(0L, as.integer(n))
  # A HYPERPARAMETER ROW IS NEVER HIDDEN. It governs every coefficient under
  # it, and a block that showed ten coefficients and dropped the smoothing
  # parameter that produced them would bury the one number the rest are
  # conditional on.
  hyp <- which(tb$role %in% c("fixed", "estimated"))
  rest <- setdiff(seq_len(nrow(tb)), hyp)
  if (nrow(tb) <= max(cap, n + length(hyp))) return(seq_len(nrow(tb)))
  sort(c(hyp, utils::head(rest, n)))
}

#' Print One Block of a Model Summary
#'
#' @description
#' The heading of a block, the term read at a glance where it is written in
#' parameters of its own, its own coefficients, and one indented compartment
#' per parameter developed over covariates.
#'
#' @details
#' A term that develops one of its own parameters carries columns that mean
#' different things: a break-point's population value and its per-group
#' deviations are not comparable quantities, and a table that stacks them
#' reads as a list of numbers and no longer as a model. Each developed
#' parameter is therefore printed as a compartment of its own, headed by what
#' develops it, opening with its hyperparameter under a name that says what
#' the hyperparameter is, and rendering each sub-term the way a block of that
#' kind is rendered at the top level. A random development reports the scale
#' of its effects and one line saying how many predictions there are and how
#' far they spread, the predictions themselves being in [coef()].
#'
#' @param b A block record from [summary_blocks()].
#' @param digits Significant digits.
#'
#' @return `NULL`, invisibly. Called for the printing.
#'
#' @keywords internal
print_block <- function(b, digits = 4L, n = NULL) {
  head <- if (is.na(b$term)) b$label else b$term
  bits <- character(0)
  if (!identical(b$kind, "parametric")) {
    bits <- c(bits, sprintf("%d %s", b$n_coef,
                            if (identical(b$kind, "structural")) "parameters"
                            else "coefficients"))
    if (is.finite(b$edf)) bits <- c(bits, sprintf("edf %.2f", b$edf))
  }
  if (identical(b$kind, "selection")) {
    bits <- c(bits, sprintf("%d selected, %d at zero",
                            b$n_coef - b$n_zero, b$n_zero))
  }
  cat("\n", head, sep = "")
  if (length(bits)) cat("   [", paste(bits, collapse = ", "), "]", sep = "")
  cat("\n")

  if (!is.null(b$head) && nrow(b$head)) {
    print_block_head(b$head, digits)
    cat("\n")
  }

  comp <- if (is.null(b$components)) list() else b$components
  # the term's own table first, then one section per compartment, each
  # carrying the rows it keeps and the free-text lines a random development
  # reports instead of its predictions
  secs <- list(list(header = NULL, indent = 2L, tb = b$table,
                    lines = character(0)))
  for (cp in comp) {
    secs[[length(secs) + 1L]] <- list(header = cp$header, indent = 4L,
                                      tb = cp$table, lines = cp$lines)
  }
  for (i in seq_along(secs)) {
    keep <- block_rows_shown(secs[[i]]$tb, n)
    secs[[i]]$hidden <- nrow(secs[[i]]$tb) - length(keep)
    secs[[i]]$tb <- secs[[i]]$tb[keep, , drop = FALSE]
  }
  ns <- vapply(secs, function(z) nrow(z$tb), integer(1))
  if (!sum(ns)) {
    cat("  (nothing to report on its own)\n")
    return(invisible(NULL))
  }
  # THE WHOLE BLOCK IS FORMATTED AT ONCE, so a compartment's numbers are
  # padded to the same widths as the table above it and the columns line up
  # down the block rather than restarting at every section
  fm <- format_block_cells(do.call(rbind, lapply(secs, function(z) z$tb)),
                           digits)
  used <- apply(fm$cells, 2L, function(z) any(nzchar(z)))
  used[[1L]] <- TRUE
  fm$cells <- fm$cells[, used, drop = FALSE]
  # a name inside a compartment is the sub-term's own, and the term's prefix
  # is dropped from its own rows: the heading already says which term this is
  nm <- fm$name
  # the COEFFICIENT rows alone: a hyperparameter's name is not built from
  # the term's, so including it in the question leaves the names no prefix in
  # common and the stripping fires nowhere
  at <- 0L
  for (i in seq_along(secs)) {
    if (ns[[i]] && (i > 1L || !is.na(b$term))) {
      ic <- at + which(secs[[i]]$tb$role == "coefficient")
      if (length(ic)) nm[ic] <- drop_common_prefix(nm[ic])
    }
    at <- at + ns[[i]]
  }
  labs <- character(0)
  ind <- integer(0)
  at <- 0L
  rows <- list()
  for (i in seq_along(secs)) {
    ii <- if (ns[[i]]) at + seq_len(ns[[i]]) else integer(0)
    at <- at + ns[[i]]
    rows[[i]] <- ii
    labs <- c(labs, nm[ii], names(secs[[i]]$lines))
    ind <- c(ind, rep(secs[[i]]$indent, ns[[i]] + length(secs[[i]]$lines)))
  }
  w <- max(nchar(labs) + ind)
  cw <- pmax(nchar(colnames(fm$cells)),
             apply(nchar(fm$cells), 2L, max))
  pad <- function(x, k) {
    mapply(function(v, j) formatC(v, width = j, flag = " "), x, k,
           USE.NAMES = FALSE)
  }
  cat(sub("[ ]+$", "", paste0(strrep(" ", w + 2L),
      paste(pad(colnames(fm$cells), cw), collapse = " "))), "\n", sep = "")
  for (i in seq_along(secs)) {
    s <- secs[[i]]
    if (!is.null(s$header)) cat("\n", strrep(" ", s$indent - 2L), s$header,
                                "\n", sep = "")
    for (r in rows[[i]]) {
      cat(sub("[ ]+$", "", paste0(
        strrep(" ", s$indent),
        formatC(nm[[r]], width = -(w - s$indent + 2L), flag = " "),
        paste(pad(fm$cells[r, ], cw), collapse = " "))), "\n", sep = "")
    }
    for (k in seq_along(s$lines)) {
      cat(sub("[ ]+$", "", paste0(
        strrep(" ", s$indent),
        formatC(names(s$lines)[[k]], width = -(w - s$indent + 2L),
                flag = " "), s$lines[[k]])), "\n", sep = "")
    }
    if (s$hidden > 0L) {
      cat(strrep(" ", s$indent),
          sprintf("... %d more, in coef()\n", s$hidden), sep = "")
    }
  }
  invisible(NULL)
}

#' The Term Read at a Glance
#'
#' @description
#' One line per parameter a term is written in, giving its value and interval
#' where it is a number and, where it is developed over covariates, the
#' population value of that development beside what develops it.
#'
#' @param hd The head record of a block.
#' @param digits Significant digits.
#'
#' @return `NULL`, invisibly. Called for the printing.
#'
#' @keywords internal
print_block_head <- function(hd, digits = 4L) {
  val <- ifelse(is.na(hd$estimate), "", format(signif(hd$estimate, digits)))
  ci <- ifelse(is.na(hd$lower) | is.na(hd$upper), "",
               sprintf("[%s, %s]", format(signif(hd$lower, digits)),
                       format(signif(hd$upper, digits))))
  w <- max(nchar(hd$name))
  for (i in seq_len(nrow(hd))) {
    cat(sub("[ ]+$", "", paste0(
      "  ", formatC(hd$name[[i]], width = -(w + 2L), flag = " "),
      formatC(val[[i]], width = max(nchar(val)), flag = " "), "  ",
      formatC(ci[[i]], width = -max(nchar(ci)), flag = " "), "  ",
      hd$note[[i]])), "\n", sep = "")
  }
  invisible(NULL)
}

#' Drop the Prefix a Set of Coefficient Names Share
#'
#' @description
#' The first dotted piece, where every name carries the same one and
#' removing it leaves every name non-empty.
#'
#' @details
#' A term composes its coefficients' names from its own and its parameters',
#' so inside the block of one term the leading piece repeats on every row
#' and says what the heading has already said. It is dropped for the
#' printing alone; [coef()] and the summary's own tables keep the
#' names the fit was built with, which are the ones another call can be
#' indexed by.
#'
#' **one piece** and never every piece they share: the term's own name is
#' one, and a set of coefficients that happen to agree further along, `r.1`
#' and `r.2` of a matrix column, would otherwise be left as the bare
#' numbers.
#'
#' @param nms The names.
#'
#' @return The names, shortened where they share a prefix.
#'
#' @keywords internal
drop_common_prefix <- function(nms) {
  if (length(nms) < 1L) return(nms)
  sp <- strsplit(nms, ".", fixed = TRUE)
  if (min(lengths(sp)) < 2L) return(nms)
  if (length(unique(vapply(sp, function(z) z[[1L]], character(1)))) != 1L) {
    return(nms)
  }
  vapply(sp, function(z) paste(z[-1L], collapse = "."), character(1),
         USE.NAMES = FALSE)
}


#' What the Fit Certifies About the Point It Reports
#'
#' @description
#' Three readings taken at the reported point and independent of the path the
#' search took: the outer criterion's gradient, how far the coefficients sit
#' above the penalized mode, and which hyperparameters have run to a boundary.
#'
#' @details
#' **Why a certificate and not the optimizer's flag.** The flag says
#' whether a search stopped on its own rule, which is a statement about the
#' search. Measured across shapes, it does not order fits by quality: on one
#' model the default reported success at a criterion of -1783.47 while the same
#' data under [optimizers7::lbfgs()] reached -1664.43 and reported
#' failure. What a reader wants is a property of the point.
#'
#' **The state comes from the gradient and the mode error is reported
#' beside it, not folded into it.** Measured at the reported point over six
#' shapes, the outer gradient separates by five orders, 4.7e-07, 7.8e-07,
#' 5.8e-05, 7.7e-05 and 3.0e-04 on fits that are right, against 28.8 on one
#' that is not, while the mode error does not: it reads 1.8e-16 to 6.1e-12
#' on
#' four of them, 22.8 on the failing one, and 0.114 on a random-changepoint
#' `seg` whose answer is right to a correlation of 0.9932. A number that
#' does not separate cannot decide a state, and a certificate that says how far
#' from the mode is worth more than a boolean that hides it.
#'
#' `tol` is 1e-2 instead of the geometric middle of the two groups: the
#' two ways of being wrong are not symmetric, and a certificate that says not
#' converged at a good point is visible and checkable where one that certifies
#' a bad point is the failure this exists to remove.
#'
#' **What it costs** is one outer gradient and one solve, once, at a point
#' the fit already holds. Nothing is refitted: measured, the criterion
#' reconstructed from `fit@spec` equals the one the fit reports exactly on
#' every shape, so the reading is of the fitted model and of no other.
#'
#' **Where there is no outer gradient there are two cases, and they get
#' different answers.** A model with **no penalty**, which covers `linpar`,
#' `nl`, `seg`, `jump` and `jseg`, has no hyperparameter for a
#' gradient to be about, so the only question left is whether the inner fit
#' reached its mode, and the mode error answers it: measured over the
#' reference battery it reads 5.2e-11 to 7.9e-05 on fits that are right
#' against 1.215 on a `jump` fitted to data carrying a slope and a slope
#' change it has no term for. A model whose only hyperparameters are
#' **kinked**,
#' `lasso`, `scad`, `mcp`, swept along a path because a Laplace
#' a Laplace approximation at a mode sitting on the kink having no meaning,
#' gets neither
#' reading and stays `"unknown"`: at a coefficient the penalty has set to
#' zero the score does not vanish but lies in the subdifferential, so the mode
#' error is not a statement about being at a mode. Measured on a lasso, its
#' 4.7e-03 is carried by a coordinate whose coefficient is exactly 0 and whose
#' score is -0.715.
#'
#' A form whose criterion has no exact gradient
#' ([outer_gradient_ok()]) is also `"unknown"` and never
#' approximated: 2p refits to difference it would cost more than the fit.
#'
#' **The boundary label, and why its threshold needs no derivation.**
#' A hyperparameter may run to an edge and belong there: on a covariate that
#' is pure noise the smoothing parameter reaches 9.2e+08, the criterion is
#' genuinely flat, and calling that fit unconverged would be wrong. A
#' coordinate is reported as sitting at a boundary when its free value
#' exceeds `edge` AND its own gradient component has already met
#' `tol`. Because of that second condition the threshold cannot change
#' the verdict: a coordinate it moves out of the interior set had already
#' passed the test, so the maximum that decides the state is unaffected, and
#' both `"converged"` and `"boundary"` are certified. What
#' `edge` decides is how the point is described. The default separates
#' the measured cases with room on both sides: coordinates that ran to an
#' edge sit at 9.3, 10.5 and 20.6 on the free scale against 0.13, 0.30 and
#' 2.01 for the ones that did not.
#'
#' @param fit A [StatmodFit()].
#' @param tol The largest outer gradient a certified point may carry.
#' @param edge The free value beyond which a hyperparameter whose gradient
#'   has already met `tol` is reported as sitting at a boundary. It
#'   decides the label alone and never the verdict; see the details.
#'
#' @return A list with `state` (`"converged"`, `"boundary"`,
#'   `"not converged"` or `"unknown"`), `gradient`,
#'   `mode_error`, `boundary` and `reason`.
#'
#' @seealso [statmod()], [mode_error_limit()],
#'   [criterion_resolution()]
#'
#' @examples
#' dd <- data.frame(x = runif(120))
#' dd$y <- sin(4 * dd$x) + rnorm(120, 0, 0.3)
#' statmod_certificate(statmod(y ~ s(x, k = 8),
#'                             distributions7::gaussian1_distrib(), dd))$state
#'
#' @export
statmod_certificate <- function(fit, tol = 1e-2, edge = 8) {
  out <- list(state = "unknown", gradient = NA_real_, mode_error = NA_real_,
              boundary = character(0), reason = character(0))
  method <- fit@methods$outer
  spec <- fit@spec
  design <- tryCatch(statmod_design(spec), error = function(e) NULL)
  if (is.null(design)) {
    out$reason <- "the design could not be rebuilt from the fit"
    return(out)
  }
  cf <- fit@coefficients
  hy <- fit@hyper
  ctx <- tryCatch(outer_context(spec, design, cf, hy, "bartlett"),
                  error = function(e) NULL)

  # HOW FAR ABOVE THE MODE, in log-likelihood units: the decrease the mode's
  # own Newton correction predicts. Reported whether or not a criterion ran.
  if (!is.null(ctx)) {
    pen <- tryCatch(ctx_penalized(ctx, spec, design, cf, hy, FALSE),
                    error = function(e) NULL)
    if (!is.null(pen)) {
      obj <- tryCatch(statmod_objective(spec, hy, design, FALSE, "bartlett"),
                      error = function(e) NULL)
      sc <- if (is.null(obj)) NULL else
        tryCatch(obj$gr(obj$stack(cf)), error = function(e) NULL)
      if (!is.null(sc) && all(is.finite(sc))) {
        db <- tryCatch(as.numeric(as.matrix(pen$inv) %*% sc),
                       error = function(e) NULL)
        if (!is.null(db) && all(is.finite(db))) {
          out$mode_error <- 0.5 * sum(sc * db)
        }
      }
    }
  }

  blocks <- tryCatch(statmod_blocks(spec, design), error = function(e) NULL)
  idx <- if (is.null(blocks)) NULL else outer_hyper_index(spec, blocks)
  no_outer <- is.null(method) || !method@kind %in% c("ml", "reml") ||
    is.null(idx) || !nrow(idx)
  if (no_outer) {
    # ⚠️ TWO QUITE DIFFERENT WAYS TO HAVE NO OUTER GRADIENT, and they call for
    # different answers. Measured over the reference battery, seven cases of
    # twenty-nine land here:
    #
    #   NO PENALTY AT ALL -- linpar, nl, seg, jump, jseg sharp and smoothed.
    #     There is no hyperparameter to estimate, so there is nothing an outer
    #     gradient could say, and the only question left is whether the inner
    #     fit reached its mode. The mode error answers exactly that, and it
    #     answers well: 5.2e-11, 1.5e-10, 1.9e-10, 2.0e-06 and 7.9e-05 on fits
    #     that are right, against 1.215 on a `jump` fitted to data carrying a
    #     slope and a slope change it has no term for -- misspecified, so its
    #     break-point iteration never settles, its annealing runs to the floor
    #     and its block reaches 3.0e+13 on the information's diagonal. On data
    #     where the jump IS the truth the same term reads 4.2e-04.
    #
    #   A KINKED HYPERPARAMETER -- lasso, scad, mcp -- chosen by a PATH
    #     because a Laplace approximation at a mode sitting on the kink is
    #     arithmetic without a meaning. There is a hyperparameter, but it is
    #     the argmin over a grid rather than the root of a derivative. AND THE
    #     MODE ERROR IS NOT A READING HERE EITHER: at a coefficient the
    #     penalty has set to zero the score does not vanish, it lies in the
    #     subdifferential. Measured on a lasso, the mode error of 4.7e-03 is
    #     carried by a coordinate whose coefficient is exactly 0 and whose
    #     score is -0.715. So this stays unknown, and says why.
    pen_units <- tryCatch(statmod_penalized(spec, design),
                          error = function(e) list())
    if (length(pen_units)) {
      out$reason <- paste0(
        "the only hyperparameters here belong to a penalty with a kink, ",
        "chosen along a path rather than by a criterion with a derivative; ",
        "and at a coefficient the penalty has set to zero the score does not ",
        "vanish, so the mode error is not a reading either")
      return(out)
    }
    if (!is.finite(out$mode_error)) {
      out$reason <- "the model carries no penalty, and the mode error could not be read"
      return(out)
    }
    out$state <- if (out$mode_error <= mode_error_limit()) "converged" else
      "not converged"
    out$reason <- sprintf(paste0(
      "the model carries no penalty, so there is no outer gradient; the ",
      "reading is the inner fit's own, %.4g log-likelihood units above its ",
      "mode against %g"), out$mode_error, mode_error_limit())
    return(out)
  }
  if (!outer_gradient_ok(spec, design, idx, method, 1L)) {
    out$reason <- "this form has no exact outer gradient, and differencing it would cost more than the fit"
    return(out)
  }
  basis <- integrated_basis(spec, design, method@kind)
  g <- tryCatch(statmod_marginal_grad(spec, design, cf, hy, method, idx, basis,
                                      ctx = ctx),
                error = function(e) NULL)
  if (is.null(g) || !all(is.finite(g))) {
    # ⚠️ A NON-FINITE OUTER GRADIENT AT A BOUNDARY IS NOT THE SAME COMPLAINT,
    # and saying which it is costs one read. Where a coefficient sits at the
    # clamp its link keeps it strictly inside, the family's THIRD and FOURTH
    # derivatives there are not finite -- measured on the Student t, all ten
    # components of the third go from nu = 1e150, which is sqrt(double.xmax)
    # and so the signature of a product of two quantities of order nu -- and
    # the outer gradient reads exactly those. The fit itself may be perfectly
    # well located: on the reference battery's `fam-studentt` the mode error
    # is 1.3e-10 against a limit of 1e-3 and the criterion is the best of the
    # three routes, so what is missing is the reading and not the answer.
    #
    # The state stays `not converged` rather than `boundary`, deliberately:
    # certifying a point whose hyperparameters were never verified would
    # claim more than has been checked. What it can honestly do is name the
    # coordinate, so the reader is not left with "not finite".
    H0 <- tryCatch(ctx_information(ctx, spec, design, cf, hy, FALSE,
                                   "bartlett"), error = function(e) NULL)
    frozen <- if (is.null(H0)) integer(0) else boundary_coords(H0)
    if (length(frozen)) {
      npar <- vapply(design, function(d) d$npar, integer(1))
      ends <- cumsum(npar)
      who <- unique(names(npar)[vapply(frozen, function(j)
        which(j <= ends)[1], integer(1))])
      out$boundary <- who
      out$reason <- sprintf(paste0(
        "the outer gradient is not finite at the reported point, because a ",
        "coefficient of %s sits at the bound its link keeps it strictly ",
        "inside and the family's third and fourth derivatives are not finite ",
        "there. The fit itself is %.4g log-likelihood units above its mode"),
        paste(who, collapse = ", "), out$mode_error)
    } else {
      out$reason <- "the outer gradient is not finite at the reported point"
    }
    out$state <- "not converged"
    return(out)
  }
  # A COORDINATE AT A BOUNDARY is one whose criterion has stopped moving in it
  # while its value has run far from where it started. Both halves are needed:
  # a small gradient component alone is what convergence looks like, and a
  # large value alone is an ordinary answer on a wide scale.
  #
  # `edge` IS NOT DERIVED FROM ANYTHING, and it does not have to be, because
  # THE CONJUNCTION BELOW IS WHAT MAKES IT SAFE: a coordinate is called an
  # edge only if it has ALREADY met `tol`, so removing it from `interior`
  # cannot raise the maximum that decides the verdict. Move the threshold in
  # either direction and the only thing that changes is whether the state
  # reads `converged` or `boundary`, both of which are certified. Delete the
  # `abs(g) <= tol` conjunct, however, and the threshold starts excusing
  # coordinates from the gradient test, at which point this paragraph is
  # false and the number needs an argument of its own.
  #
  # What the default separates, measured: the coordinates that ran to an edge
  # sit at |eta| of 9.3, 10.5 and 20.6 -- a smoothing parameter of 9.2e+08 on
  # pure noise, prior scales of 9.2e-05 and 2.8e-05 -- against 0.13, 0.30 and
  # 2.01 for the ones that did not, so 8 sits in a wide gap rather than on a
  # boundary of its own.
  eta <- hyper_to_eta(hy, idx)
  at_edge <- which(abs(eta) > edge & abs(g) <= tol)
  if (length(at_edge)) {
    out$boundary <- vapply(at_edge, function(k)
      paste(idx$parameter[k], idx$term[k], idx$name[k], sep = "/"),
      character(1))
  }
  interior <- setdiff(seq_along(g), at_edge)
  out$gradient <- if (length(interior)) max(abs(g[interior])) else 0
  if (out$gradient <= tol) {
    out$state <- if (length(at_edge)) "boundary" else "converged"
  } else {
    out$state <- "not converged"
    out$reason <- sprintf(
      "the outer criterion's gradient is %.4g at the reported point, against %g",
      out$gradient, tol)
  }
  if (is.finite(out$mode_error) && out$mode_error > mode_error_limit()) {
    out$reason <- c(out$reason, sprintf(
      paste0("the coefficients sit %.4g log-likelihood units above the",
             " penalized mode, so everything read there -- the criterion,",
             " its gradient, vcov() -- carries that"),
      out$mode_error))
  }
  out
}
