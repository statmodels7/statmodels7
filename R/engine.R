#' The Compiled Inner Fit
#'
#' @description
#' `engine_inner_fit()` runs the inner penalized fit of [fit_smooth()] in
#' compiled code (`src/engine_iwls.cpp`) when every component of the model
#' has a compiled route, and returns `NULL` otherwise, in which case
#' [fit_smooth()] calls [iwls_fit()]. The compiled fit performs the
#' operations of [iwls_fit()] in the same order, with the families and links
#' read through the registries of distributions7 and linkfunctions7.
#' `engine_route()` decides whether a model is covered and assembles the
#' inputs.
#'
#' @details
#' A model is covered when `spec@engine` is `"auto"`; the method is an
#' [iwls()] without a `criterion` whose expected information is the closed
#' form (`method = NULL`), [distributions7::expect_opg()],
#' [distributions7::expect_mc()] or [distributions7::expect_integrate()] (the
#' same object for the information and its backup); the block is the
#' whole coefficient vector, with no break-point kink, no moving block and no
#' nonlinear curvature unit; the family is a univariate family of the
#' registry without a wrapper and every link has a compiled route; no term
#' carries a structural state, a Jeffreys penalty or a mixed penalty unit;
#' and every penalized unit is one that `engine_penalty_units()` describes.
#'
#' Where the R route would signal an error the compiled fit stops and
#' `engine_inner_fit()` returns `NULL`, so that [iwls_fit()] reproduces it.
#' With [distributions7::expect_mc()] the random stream is then restored to
#' its state before the compiled fit, so that [iwls_fit()] makes the same
#' draws as when it runs alone.
#'
#' @param spec,design,obj,hyper As in [fit_smooth()].
#' @param beta The full starting vector, with the held values written in.
#' @param method A settled [iwls()].
#' @param groups,frozen As [fit_smooth()] passes them to [iwls_fit()].
#' @param verbose Whether to print the iterations.
#' @param damp_on_reject,guard As in [iwls_fit()].
#'
#' @return `engine_inner_fit()`: the list that [iwls_fit()] returns, or
#'   `NULL`. `engine_route()`: `NULL`, or a list with the elements `model`,
#'   `pen` and `method` that the compiled fit reads.
#'
#' @seealso [fit_smooth()], [iwls_fit()], [statmod()] (argument `engine`)
#' @keywords internal
engine_inner_fit <- function(spec, design, obj, beta, hyper, method, groups,
                             frozen, verbose, damp_on_reject, guard) {
  rt <- engine_route(spec, design, obj, beta, hyper, method)
  if (is.null(rt)) return(NULL)
  rt$method$damp_on_reject <- isTRUE(damp_on_reject)
  rt$method$guard <- isTRUE(guard)
  # a simulated expected information draws from R's stream: where the
  # compiled fit stops part way, the stream is put back where it was, so
  # that iwls_fit() makes the draws that the R route makes
  mc <- identical(rt$method$emethod, 2L)
  if (mc) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  res <- engine_iwls_fit(rt$model, rt$pen, rt$method, as.numeric(beta),
                         as.integer(frozen), lapply(groups, as.integer),
                         isTRUE(verbose))
  if (!identical(res$status, 0L)) {
    if (mc) {
      if (had_seed) assign(".Random.seed", seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    }
    return(NULL)
  }
  note <- switch(as.character(res$note_kind),
    "1" = paste0("the score is not finite at iteration ", res$note_it,
                 ", so the run stopped at the last usable point"),
    "2" = paste0("the line search found no acceptable step at iteration ",
                 res$note_it, if (isTRUE(damp_on_reject)) ", even with damping"),
    "3" = paste0("the objective stopped moving at iteration ", res$note_it),
    "4" = paste0("converged by the dimensionless rule (score ",
                 format(res$dimless, digits = 3), " against ",
                 format(res$score, digits = 3),
                 " absolute; the response's scale is the gap)"),
    NULL)
  if (isTRUE(verbose)) {
    cat(sprintf("  %-5s %14.6f %12.3e %10s\n", "end", res$value, res$score,
                if (res$converged) "converged" else "stopped"))
    if (!is.null(note)) cat("  ", note, "\n", sep = "")
  }
  hist <- res$history
  if (!is.null(hist)) {
    hist <- data.frame(iteration = hist$iteration, objective = hist$objective,
                       score = hist$score, step = hist$step, rank = hist$rank,
                       route = hist$route, held = hist$held, damp = hist$damp)
  }
  list(par = res$par, value = res$value, converged = res$converged,
       iterations = res$iterations, score = res$score, note = note,
       aliased = res$aliased, fallback = res$fallback,
       switched = res$switched, curvature = c(full = 0L, search = 0L),
       history = hist)
}

#' @rdname engine_inner_fit
#' @keywords internal
engine_route <- function(spec, design, obj, beta, hyper, method) {
  if (!identical(spec@engine, "auto")) return(NULL)
  if (!S7::S7_inherits(method, Iwls) || !is.null(method@criterion)) {
    return(NULL)
  }
  # the expectation methods with a compiled route: the closed form (NULL),
  # the outer product of the scores, simulation and integration; the
  # information and its backup must use the same one
  em <- c(list(info_emethod(method@hessian)),
          if (!is.null(method@backup)) list(info_emethod(method@backup)))
  ecode <- vapply(em, function(m) {
    if (is.null(m)) 0L
    else if (S7::S7_inherits(m, distributions7::ExpectOPG)) 1L
    else if (S7::S7_inherits(m, distributions7::ExpectMC)) 2L
    else if (S7::S7_inherits(m, distributions7::ExpectIntegrate)) 3L
    else NA_integer_
  }, integer(1))
  if (anyNA(ecode) || length(unique(ecode)) > 1L) return(NULL)
  if (length(em) > 1L && !identical(em[[1L]], em[[2L]])) return(NULL)
  how <- solver_how(method@solver)
  em_ <- engine_model(spec, design, obj$split(beta), hyper, how)
  if (is.null(em_)) return(NULL)
  model <- em_$model
  d <- spec@distrib
  st <- solver_tol(method@solver)
  sw <- method@switch_after
  mcfg <- list(
    how = how, stol = if (is.null(st)) NA_real_ else as.numeric(st),
    floor = solver_floor(method@solver),
    maxit = as.integer(method@maxit), tol = method@tol,
    step_halving = as.integer(method@step_halving),
    decrement_limit = method@decrement_limit,
    switch_after = if (length(sw)) as.numeric(sw) else Inf,
    switch_rate = if (length(method@switch_rate)) as.numeric(method@switch_rate)
                  else NA_real_,
    tol_def = iwls_definite_tol(method),
    expected = info_expected(method@hessian),
    backup = !is.null(method@backup),
    switching = isTRUE(method@switching) && length(sw) && is.finite(sw),
    polish = isTRUE(iwls_polish()),
    emethod = ecode[[1L]],
    escore = ecode[[1L]] >= 2L && identical(em[[1L]]@identity, "score"),
    nsim = if (ecode[[1L]] == 2L) as.numeric(em[[1L]]@nsim) else NA_real_,
    etol = if (ecode[[1L]] == 3L) {
      m <- em[[1L]]
      as.numeric(c(m@quad_atol, m@quad_rtol, m@max_panels, m@series_atol,
                   m@series_rtol, m@max_terms))
    } else rep(NA_real_, 6L),
    support = as.numeric(d@bounds)
  )
  list(model = model, pen = em_$pen, method = mcfg)
}

#' @rdname engine_inner_fit
#' @param coef The coefficients, split by distribution parameter.
#' @param how The solver of the inner fit, `"qr"`, `"svd"` or `"chol"`,
#'   which decides whether the square root of the penalty is formed.
#' @details
#' `engine_model()` returns the data of the model as the compiled routes
#' read them, `list(model, pen, covered)`, or `NULL` where the model is not
#' covered: the family, its links and the design blocks, and the penalty at
#' `hyper` (its units, its Hessian and, for the QR and SVD solvers, its
#' square root). `covered` is the highest derivative order of the family in
#' the registry.
#' @keywords internal
engine_model <- function(spec, design, coef, hyper, how) {
  if (has_sharp_breakpoint(spec) || length(attr(design, "refresh")) ||
      length(nl_curvature_units(design)) ||
      !is.null(statmod_structural_state(design))) {
    return(NULL)
  }
  d <- spec@distrib
  params <- d@params
  if (S7::S7_inherits(d, distributions7::multivariate_distrib)) return(NULL)
  route <- distributions7::distrib_scalar_route(d)
  cls <- attr(S7::S7_class(d), "name")
  if (is.null(route) || !identical(route$name, cls)) return(NULL)
  lk <- lapply(params, function(p) {
    r <- linkfunctions7::link_scalar_route(d@link_params[[p]])
    if (!is.null(r)) r$bounds <- as.numeric(d@link_params[[p]]@link_bounds)
    r
  })
  if (any(vapply(lk, is.null, logical(1)))) return(NULL)
  pr <- engine_family_probe(cls, vapply(lk, `[[`, character(1), "name"))
  if (pr$id < 0L || pr$id >= 100000L || pr$covered < 2L ||
      any(pr$links < 0L)) {
    return(NULL)
  }
  units <- statmod_penalized(spec, design)
  if (any(vapply(units, function(u) {
    isTRUE(u$mixed) || isTRUE(u$jeffreys) || isTRUE(u$structural)
  }, logical(1)))) {
    return(NULL)
  }
  if (!is.null(jfp_penalty(spec, design, coef, "parts"))) return(NULL)
  n <- spec@n_obs
  yv <- spec@response
  if (!is.numeric(yv) || !is.null(dim(yv)) || length(yv) != n) return(NULL)
  X <- lapply(params, function(p) {
    blk <- design[[p]]
    if (!is.null(blk$adj)) return(FALSE)
    if (blk$npar == 0L) return(NULL)
    x <- blk$X
    if (isS4(x)) {
      x <- tryCatch(methods::as(methods::as(x, "generalMatrix"),
                                "CsparseMatrix"), error = function(e) FALSE)
      if (!methods::is(x, "dgCMatrix")) return(FALSE)
    } else if (!is.matrix(x) || !is.double(x)) {
      return(FALSE)
    }
    x
  })
  if (any(vapply(X, isFALSE, logical(1)))) return(NULL)
  sparse <- design_sparse(design)

  # the penalty at these hyperparameters: its units, for the value and the
  # gradient (D17), and its Hessian S, for the curvature
  pu <- engine_penalty_units(spec, design, units, hyper)
  if (is.null(pu)) return(NULL)
  S <- statmod_penalty_at(spec, coef, hyper, design, "hessian",
                          jeffreys = FALSE)
  Sd <- as_dense(S)
  C <- if (how %in% c("qr", "svd")) penalty_sqrt(S) else NULL
  if (sparse) {
    S <- methods::as(methods::as(S, "generalMatrix"), "CsparseMatrix")
    if (!is.null(C) && isS4(C)) {
      C <- methods::as(methods::as(C, "generalMatrix"), "CsparseMatrix")
    }
  } else {
    S <- Sd
  }
  if (!is.null(C) && isS4(C) && !methods::is(C, "dgCMatrix")) C <- as_dense(C)

  bounds <- d@params_bounds
  cons <- route$constants
  cm <- if (length(cons)) {
    do.call(cbind, lapply(cons, function(v) rep_len(as.numeric(v), n)))
  } else {
    matrix(0, n, 0L)
  }
  offsets <- lapply(params, function(pp) {
    o <- spec@offsets[[pp]]
    if (is.null(o)) NULL else rep_len(as.numeric(o), n)
  })
  model <- list(
    X = X, y = as.numeric(spec@response), w = rep_len(as.numeric(spec@weights), n),
    offsets = offsets, sparse = sparse, family = cls, constants = cm,
    lo = vapply(params, function(pp) as.numeric(bounds[[pp]][1]), numeric(1)),
    hi = vapply(params, function(pp) as.numeric(bounds[[pp]][2]), numeric(1)),
    links = lk, threads = as.integer(spec@threads)
  )
  list(model = model, pen = list(S = S, C = C, units = pu), covered = pr$covered)
}

#' @rdname engine_inner_fit
#' @param units The penalized units, [statmod_penalized()].
#' @details
#' `engine_penalty_units()` describes each penalized unit to the compiled
#' fit by the quantities that its `penalty_value()` and `penalty_gradient()`
#' methods read at the given hyperparameters, so that the compiled value and
#' gradient repeat those methods operation by operation: for a
#' `QuadraticPenalty` the map, `P`, `penalty_matrix` and the three constants;
#' for an `AdditivePenalty` the summed matrix of `additive_sum()`, the
#' matrices and their smoothing parameters; for a `DistribPenalty` with a
#' gaussian prior of fixed mean, the mean and the standard deviation, or, for
#' a multivariate one, the mean, the factor of `parameters7::param_factor()`,
#' the log-determinant and the inverse covariance. It returns `NULL` for any
#' other unit.
#' @keywords internal
engine_penalty_units <- function(spec, design, units, hyper) {
  out <- vector("list", length(units))
  for (k in seq_along(units)) {
    # a unit with pieces (a covariance class) is read and updated at its
    # stacked positions u$index, as statmod_penalty_at() reads and updates it
    u <- units[[k]]
    pen <- u$penalty
    th <- as.list(hyper[[u$param]][[u$key]])
    idx <- as.integer(u$index)
    mp <- if ("map" %in% S7::prop_names(pen)) pen@map else NULL
    if (!is.null(mp) && !is.matrix(mp)) mp <- as_dense(mp)
    if (S7::S7_inherits(pen, penalties7::QuadraticPenalty)) {
      lam <- th[[1]]
      r <- pen@rank
      out[[k]] <- list(kind = 1L, idx = idx, map = mp,
                       P = engine_matrix(pen@P), M = engine_matrix(pen@penalty_matrix),
                       lam = lam, a = lam / 2, c1 = r / 2 * log(lam), rank = r,
                       c2 = r / 2 * log(2 * pi), c3 = pen@penalty_logpdet / 2)
    } else if (S7::S7_inherits(pen, penalties7::AdditivePenalty)) {
      a <- penalties7::additive_sum(pen, th)
      out[[k]] <- list(kind = 2L, idx = idx, map = NULL, S = engine_matrix(a$S),
                       c1 = 0.5 * a$logpdet, c2 = pen@rank / 2 * log(2 * pi),
                       matrices = lapply(pen@matrices, engine_matrix),
                       lams = as.numeric(unlist(th[pen@params])))
    } else if (S7::S7_inherits(pen, penalties7::DistribPenalty)) {
      par <- pen@parent
      # a fixed() wrapper, by its two properties (its classes are internal)
      if (!all(c("parent_distrib", "fixed_params") %in% S7::prop_names(par))) {
        return(NULL)
      }
      full <- distributions7::fixed_full_theta(par, th)
      lj <- if (is.null(mp)) NA_real_ else penalties7::dp_log_jacobian(pen)
      parent <- par@parent_distrib
      if (identical(attr(S7::S7_class(parent), "name"), "Gaussian1Distrib") &&
          pen@block == 1L) {
        out[[k]] <- list(kind = 3L, idx = idx, map = mp,
                         mu = as.numeric(full[[1]]), sigma = as.numeric(full[[2]]),
                         lj = as.numeric(lj))
        if (length(out[[k]]$mu) != 1L || length(out[[k]]$sigma) != 1L) return(NULL)
      } else if (S7::S7_inherits(parent, distributions7::MvGaussianDistrib) &&
                 pen@block == parent@n_dim) {
        pc <- distributions7::mvg_pieces(parent, full)
        l <- unname(unclass(parameters7::param_factor(pc$s, pc$eta)))
        out[[k]] <- list(kind = 4L, idx = idx, map = mp, block = as.integer(pen@block),
                         mu = as.numeric(pc$mu), L = as.matrix(l),
                         inverted = isTRUE(parent@inverted),
                         logdet = as.numeric(pc$logdet), p2pi = pc$p * base::log(2 * pi),
                         sinv = unname(as.matrix(pc$sigma_inv)), lj = as.numeric(lj))
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  }
  out
}

#' @rdname engine_inner_fit
#' @param x A base matrix or a `Matrix` object.
#' @details
#' `engine_matrix()` returns a penalty matrix as the compiled fit reads it:
#' a base matrix as it is, any `Matrix` object as a `dgCMatrix`, since the R
#' route multiplies by it as a sparse matrix.
#' @keywords internal
engine_matrix <- function(x) {
  if (!isS4(x)) return(as.matrix(x))
  methods::as(methods::as(methods::as(x, "dMatrix"), "generalMatrix"), "CsparseMatrix")
}
