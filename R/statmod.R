#' @include blocks.R
NULL

#' A Fitted Model
#'
#' @description
#' What \code{\link{statmod}} returns: the specification kept whole, the
#' coefficients and hyperparameters it reached, and the record of how it got
#' there.
#'
#' @param spec The \code{\link{StatmodSpec}} that was fitted.
#' @param coefficients A named list, one vector per distribution parameter.
#' @param hyper The hyperparameters, per penalized term.
#' @param loglik The maximized weighted log-likelihood.
#' @param objective The value of the penalized objective.
#' @param edf Effective degrees of freedom, per term and total.
#' @param fitted The fitted parameters, per distribution parameter.
#' @param converged Whether every loop stopped on its own rule.
#' @param elapsed The elapsed time, in seconds.
#' @param criterion The marginal criterion at the estimated hyperparameters,
#'   \code{NA} when none was used.
#' @param history A list of data frames: \code{outer}, \code{blocks},
#'   \code{inner}.
#' @param methods What fitted each block.
#'
#' @return An object of class \code{StatmodFit}.
#'
#' @seealso \code{\link{statmod}}
#'
#' @examples
#' dd <- data.frame(y = rnorm(30), x = runif(30))
#' S7::S7_inherits(statmod(y ~ x, distributions7::gaussian1_distrib(), dd),
#'                 StatmodFit)
#'
#' @name StatmodFit-class
#' @aliases StatmodFit
#' @keywords internal
#' @export
StatmodFit <- S7::new_class("StatmodFit",
  properties = list(
    spec = S7::class_any,
    coefficients = S7::class_list,
    hyper = S7::class_list,
    loglik = S7::class_numeric,
    objective = S7::class_numeric,
    edf = S7::class_any,
    fitted = S7::class_list,
    converged = S7::class_logical,
    elapsed = S7::class_numeric,
    criterion = S7::class_numeric,
    history = S7::class_list,
    methods = S7::class_list,
    call = S7::class_any
  )
)


#' Fit a Model
#'
#' @description
#' Reads one formula carrying every parameter of a distribution, assembles the
#' terms it names into a penalized likelihood, and fits it.
#'
#' @details
#' \strong{The formula.} The equations of the distribution's parameters are
#' separated by \code{|}, the first carrying the response:
#' \preformatted{    y ~ x1 + ridge(R) + lasso(L)  |  sigma ~ z  |  nu ~ 1}
#' A parameter with no equation gets an intercept. See
#' \code{\link{statmod_equations}}, whose recovery is not the obvious one.
#'
#' \strong{The fitting scheme.} The terms split in two by a property each one
#' already reports. Every term whose penalty is twice differentiable in its
#' coefficients -- an unpenalized block, a ridge, a spline, a random effect --
#' is estimated in ONE system by \code{inner_method}, because their joint
#' curvature exists and using it is what makes a fit converge in a handful of
#' iterations. A term whose penalty has a kink -- lasso, scad, mcp -- is
#' estimated by a method of its own with everything else held fixed. The fit
#' alternates between the two until the objective and every block stop moving.
#'
#' \strong{The objective is unaveraged}: minus the weighted log-likelihood
#' plus the penalties at full size, since a penalty is a negative log-prior
#' and a posterior adds the two at full size. What is scaled instead is the
#' stopping rule, so that a threshold means the same thing at \eqn{n = 10} and
#' at \eqn{n = 10^7}.
#'
#' \strong{The budget and the stopping rule belong to the method.} There is no
#' \code{maxit} and no \code{tol} here: they are set on \code{inner_method},
#' which is \code{\link{iwls}(maxit =, tol =)} or an optimizer with its own
#' \code{maxit} and \code{criterion}, and the alternation reads them from there
#' (see \code{\link{method_budget}}). Carrying a second copy would let a caller
#' set both and be obeyed by neither.
#'
#' \strong{The hyperparameters.} With \code{outer_method = NULL}, the default,
#' each one sits where \code{hyper} put it, or at the probe value of its bounds
#' otherwise -- a placeholder rather than a choice, and it matters, since a
#' lasso at \eqn{\lambda = 1} against an unaveraged log-likelihood of a few
#' hundred observations selects nothing at all. With
#' \code{outer_method = \link{reml}()} or \code{\link{ml}()} they are estimated
#' by a marginal criterion, \code{outer_optimizer} searching over them and the
#' coefficients being refitted at each. Only a twice differentiable penalty
#' takes part: a lasso, a SCAD or an MCP keeps the value it was given.
#'
#' \strong{Verbosity} has three levels, naming the loops rather than counting
#' them: \code{1} the outer search and the alternation, \code{2} the inner
#' method's own iterations, \code{3} the optimizers' traces as well. A named
#' form is accepted too, as \code{verbose = c(outer = TRUE, blocks = FALSE)},
#' since watching the hyperparameters move while silencing a chatty inner
#' optimizer is the common case.
#'
#' @param formula The model formula.
#' @param distrib A \pkg{distributions7} distribution object.
#' @param data A data frame.
#' @param weights Optional prior weights, taken as given and not normalized.
#' @param offsets Optional named list of offsets, one per parameter.
#' @param inner_method How the smooth block is fitted: \code{\link{iwls}()} or
#'   an \pkg{optimizers7} optimizer.
#' @param outer_method How the hyperparameters are estimated:
#'   \code{\link{reml}()}, \code{\link{ml}()}, or \code{NULL} to hold them.
#' @param outer_optimizer The optimizer that searches over them.
#' @param hyper Optional hyperparameters, a named list of named lists as
#'   \code{list(mu = list(lasso = c(lambda = 5)))}. They are held at these
#'   values.
#' @param start Optional starting coefficients, a named list.
#' @param verbose A level from 0 to 3, or a named logical vector.
#'
#' @return An object of class \code{\link{StatmodFit}}.
#'
#' @seealso \code{\link{statmod_spec}}, \code{\link{iwls}},
#'   \code{\link{loglik}}
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.5)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' fit
#'
#' # every parameter can be modelled
#' statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#'
#' @export
statmod <- function(formula, distrib, data, weights = NULL, offsets = NULL,
                    inner_method = iwls(), outer_method = NULL,
                    outer_optimizer = optimizers7::nelder_mead(),
                    hyper = NULL, start = NULL, verbose = 0) {
  t0 <- proc.time()[["elapsed"]]
  cl <- match.call()
  vb <- verbosity(verbose)
  budget <- method_budget(inner_method)
  maxit <- budget$maxit
  tol <- budget$tol

  spec <- statmod_spec(formula, distrib, data, weights, offsets)
  design <- statmod_design(spec)
  hyper <- statmod_hyper_merge(spec, statmod_hyper_start(spec), hyper)
  blocks <- statmod_blocks(spec, design)

  expected <- !S7::S7_inherits(inner_method, Iwls) ||
    identical(inner_method@hessian, "expected")
  approx <- if (S7::S7_inherits(inner_method, Iwls)) inner_method@approx else
    "bartlett"
  obj <- statmod_objective(spec, hyper, design, expected, approx)
  beta <- statmod_start(spec, design, obj, start)

  if (is.null(outer_method)) {
    res <- statmod_alternate(spec, design, blocks, hyper, inner_method, beta,
                             expected, approx, maxit, tol, vb)
    crit <- NA_real_
  } else {
    if (!S7::S7_inherits(outer_method, OuterMethod)) {
      stop("'outer_method' must be reml(), ml() or NULL.", call. = FALSE)
    }
    res <- outer_fit(spec, design, blocks, hyper, inner_method, outer_method,
                     outer_optimizer, beta, approx, maxit, tol, vb)
    hyper <- res$hyper
    crit <- res$criterion
  }

  coef <- res$obj$split(res$par)
  fitted <- statmod_eta(spec, design, coef)$theta
  StatmodFit(
    spec = spec, coefficients = coef, hyper = hyper,
    loglik = statmod_loglik_at(spec, coef, design),
    objective = res$value,
    edf = statmod_edf(spec, coef, design, hyper, expected, approx),
    fitted = fitted, converged = res$converged,
    elapsed = proc.time()[["elapsed"]] - t0,
    criterion = crit,
    history = list(
      outer = res$hist_outer,
      blocks = if (length(res$hist_blocks)) do.call(rbind, res$hist_blocks)
        else NULL,
      inner = if (length(res$hist_inner)) do.call(rbind, res$hist_inner)
        else NULL
    ),
    methods = list(smooth = inner_method, outer = outer_method,
                   sparse = vapply(blocks$sparse, function(b)
                     paste(b$param, b$term, sep = "/"), character(1))),
    call = cl
  )
}


#' The Alternation Between the Smooth Block and the Rest
#'
#' @description
#' Fits the terms whose penalties are twice differentiable in one system and
#' each remaining block by a method of its own, sweeping until the objective
#' stops moving.
#'
#' @details
#' It is a function of its own because the outer search calls it once per
#' hyperparameter it tries, warm-started from the previous coefficients.
#'
#' @param spec The specification.
#' @param design The design.
#' @param blocks The block split.
#' @param hyper The hyperparameters.
#' @param inner_method How the smooth block is fitted.
#' @param beta The starting coefficients, stacked.
#' @param expected Whether the information is the expected one.
#' @param approx The approximation for the expected information.
#' @param maxit,tol The budget and the tolerance.
#' @param vb The resolved verbosity.
#'
#' @return A list with \code{par}, \code{value}, \code{converged}, \code{obj},
#'   \code{hist_blocks} and \code{hist_inner}.
#'
#' @seealso \code{\link{statmod}}
#'
#' @keywords internal
statmod_alternate <- function(spec, design, blocks, hyper, inner_method, beta,
                              expected, approx, maxit, tol, vb) {
  obj <- statmod_objective(spec, hyper, design, expected, approx)
  value <- obj$fn(beta)
  hist_blocks <- list()
  hist_inner <- list()
  converged <- FALSE

  for (sweep in seq_len(as.integer(maxit))) {
    before <- value

    # the smooth block, all of it at once
    if (length(blocks$smooth)) {
      if (vb$blocks) {
        cat(sprintf("[sweep %d] smooth block: %d coefficients\n", sweep,
                    length(blocks$smooth)))
      }
      res <- fit_smooth(obj, beta, blocks$smooth, spec, design, hyper,
                        inner_method, vb)
      beta <- res$par
      value <- res$value
      hist_blocks[[length(hist_blocks) + 1L]] <- data.frame(
        sweep = sweep, block = "smooth", objective = value,
        change = before - value, iterations = res$iterations,
        converged = res$converged
      )
      if (!is.null(res$history)) {
        res$history$sweep <- sweep
        res$history$block <- "smooth"
        hist_inner[[length(hist_inner) + 1L]] <- res$history
      }
    }

    # each non-smooth block in turn, the others held fixed
    for (bl in blocks$sparse) {
      v0 <- value
      if (vb$blocks) {
        cat(sprintf("[sweep %d] %s in %s: %d coefficients\n", sweep,
                    bl$term, bl$param, length(bl$index)))
      }
      res <- sparse_fit(obj, beta, bl, hyper, verbose = vb$optimizer)
      beta <- res$par
      value <- res$value
      hist_blocks[[length(hist_blocks) + 1L]] <- data.frame(
        sweep = sweep, block = paste(bl$param, bl$term, sep = "/"),
        objective = value, change = v0 - value,
        iterations = res$iterations, converged = res$converged
      )
    }

    rel <- abs(before - value) / max(1, abs(value))
    if (vb$blocks) {
      cat(sprintf("[sweep %d] objective %.8f, relative change %.3e\n",
                  sweep, value, rel))
    }
    if (!length(blocks$sparse) || rel < tol) {
      converged <- TRUE
      break
    }
  }
  list(par = beta, value = value, converged = converged, obj = obj,
       hist_blocks = hist_blocks, hist_inner = hist_inner)
}


#' The Budget and the Stopping Rule of the Alternation
#'
#' @description
#' Reads the iteration budget and the tolerance off the method that fits the
#' smooth block, which is where they are set.
#'
#' @details
#' \code{\link{statmod}} carries neither a \code{maxit} nor a \code{tol} of its
#' own. An argument accepted and ignored is worse than one that signals an
#' error, and that is what a second copy would be: a caller setting
#' \code{iwls(maxit = 20)} and a loose \code{maxit = 100} would get one of them
#' with nothing said about the other. \pkg{distributions7}'s
#' \code{fit_distrib()} shed the same pair for the same reason.
#'
#' \code{\link{iwls}} carries both directly. An \pkg{optimizers7} optimizer
#' carries \code{maxit} and a \code{criterion}; the tolerance is the largest
#' one the criterion tree contains, since a combined rule stops at whichever of
#' its parts fires first and the alternation should not ask for more precision
#' than the loop inside it can deliver. A criterion carrying no tolerance at
#' all leaves the default of \code{\link[optimizers7]{crit_grad}}, read from
#' that function rather than copied as a number.
#'
#' @param method \code{\link{iwls}()} or an \pkg{optimizers7} optimizer.
#'
#' @return A list with \code{maxit} and \code{tol}.
#'
#' @seealso \code{\link{statmod}}, \code{\link{iwls}}
#'
#' @keywords internal
method_budget <- function(method) {
  if (S7::S7_inherits(method, Iwls)) {
    return(list(maxit = as.integer(method@maxit), tol = method@tol))
  }
  if (!S7::S7_inherits(method, optimizers7::optimizer)) {
    stop(paste0("'inner_method' must be iwls() or an optimizers7 optimizer.\n",
                "  The budget and the stopping rule are read off it, so there",
                " is\n  nowhere else for them to come from."), call. = FALSE)
  }
  list(maxit = as.integer(method@maxit),
       tol = criterion_tol(method@criterion))
}


#' The Tolerance a Criterion Asks For
#'
#' @description
#' The largest \code{tol} in a criterion, walking a combined one into its
#' parts.
#'
#' @param crit An \pkg{optimizers7} criterion.
#'
#' @return A single number.
#'
#' @keywords internal
criterion_tol <- function(crit) {
  nms <- S7::prop_names(crit)
  if ("criteria" %in% nms) {
    parts <- vapply(crit@criteria, criterion_tol, numeric(1))
    return(max(parts))
  }
  if ("tol" %in% nms) return(as.numeric(crit@tol)[1L])
  eval(formals(optimizers7::crit_grad)$tol)
}


#' Fit the Smooth Block
#'
#' @description
#' Runs \code{inner_method} on the jointly fitted coefficients, the others
#' held fixed.
#'
#' @param obj The objective.
#' @param beta The current stacked coefficients.
#' @param idx The smooth block's indices.
#' @param spec The specification.
#' @param design The design.
#' @param hyper The hyperparameters.
#' @param method \code{\link{iwls}()} or an optimizer.
#' @param vb The resolved verbosity.
#'
#' @return A list with \code{par}, \code{value}, \code{converged},
#'   \code{iterations} and \code{history}.
#'
#' @keywords internal
fit_smooth <- function(obj, beta, idx, spec, design, hyper, method, vb) {
  whole <- length(idx) == length(beta)
  fn <- function(b) {
    v <- beta
    v[idx] <- b
    obj$fn(v)
  }
  gr <- function(b) {
    v <- beta
    v[idx] <- b
    obj$gr(v)[idx]
  }

  if (S7::S7_inherits(method, Iwls)) {
    sub <- list(fn = fn, gr = gr, npar = length(idx),
                split = obj$split, stack = obj$stack)
    pieces_at <- function(b) {
      v <- beta
      v[idx] <- b
      p <- iwls_pieces(spec, design, obj$split(v), hyper, method)
      if (whole) return(p)
      # a subset of the coefficients takes the corresponding submatrix; the
      # square-root routes need the whole design, so the assembled one is used
      A <- if (is.null(p$A)) crossprod(p$R) + crossprod(p$C) else p$A
      list(R = NULL, C = NULL, A = A[idx, idx, drop = FALSE])
    }
    res <- iwls_fit(sub, beta[idx], method, spec@n_obs, pieces_at,
                    verbose = vb$inner)
    out <- beta
    out[idx] <- res$par
    return(list(par = out, value = obj$fn(out), converged = res$converged,
                iterations = res$iterations, history = res$history))
  }

  res <- optimizers7::minimize(method, fn, beta[idx], gr = gr)
  out <- beta
  out[idx] <- res@par
  list(par = out, value = obj$fn(out), converged = res@converged,
       iterations = res@iterations, history = NULL)
}


#' Starting Coefficients
#'
#' @description
#' An intercept-only fit per distribution parameter, with every other
#' coefficient at zero.
#'
#' @details
#' \code{\link[distributions7]{distrib_start}} gives a starting value computed
#' from the data rather than from the parameter's domain, which is what makes
#' a fit arrive in a handful of iterations instead of spending its budget
#' travelling. The penalized blocks start at zero, where their penalty is
#' smallest.
#'
#' @param spec The specification.
#' @param design The design.
#' @param obj The objective.
#' @param start Optional user starting values, a named list.
#'
#' @return A stacked numeric vector.
#'
#' @keywords internal
statmod_start <- function(spec, design, obj, start = NULL) {
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  out <- lapply(design, function(d) numeric(d$npar))
  names(out) <- params

  th0 <- tryCatch(
    distributions7::distrib_start(spec@distrib, spec@response, 1L),
    error = function(e) NULL)
  for (p in params) {
    if (design[[p]]$npar == 0L) next
    v <- if (!is.null(th0) && !is.null(th0[[p]])) th0[[p]][[1L]] else NULL
    if (is.null(v) || !is.finite(v)) next
    # the intercept carries the whole of it when there is one
    if (identical(design[[p]]$coef_names[1L], "(Intercept)") ||
        grepl("\\(Intercept\\)$", design[[p]]$coef_names[1L])) {
      out[[p]][1L] <- linkfunctions7::linkfun(links[[p]], v)
    }
  }
  if (!is.null(start)) {
    for (p in names(start)) {
      if (!p %in% params) {
        stop(sprintf("'start' names '%s', which is not a parameter.", p),
             call. = FALSE)
      }
      if (length(start[[p]]) != design[[p]]$npar) {
        stop(sprintf("'start$%s' has length %d but '%s' has %d coefficients.",
                     p, length(start[[p]]), p, design[[p]]$npar),
             call. = FALSE)
      }
      out[[p]] <- as.numeric(start[[p]])
    }
  }
  obj$stack(out)
}


#' Effective Degrees of Freedom, Per Term
#'
#' @description
#' Asks each term what it spends, through \pkg{modelterms7}'s \code{edf()}.
#'
#' @details
#' A smooth penalized term counts \eqn{\mathrm{tr}[(H+S)^{-1}H]} over its own
#' block, so it needs the unpenalized curvature there and not only its
#' coefficients and its hyperparameters. That block is cut out of the
#' likelihood's information, which is computed once for every term rather than
#' per term.
#'
#' The arguments are passed BY NAME. \code{edf()}'s third argument is the
#' curvature and its fourth the hyperparameters, and a positional call put the
#' hyperparameters where the curvature belongs: every smooth term then reported
#' \code{NA}, the total degrees of freedom counted the unpenalized terms alone,
#' and AIC and BIC were built on that count.
#'
#' @param spec The specification.
#' @param coef The coefficients.
#' @param design The design.
#' @param hyper The hyperparameters.
#' @param expected Whether the curvature is the expected information.
#' @param approx The approximation for the expected information.
#'
#' @return A data frame with one row per term.
#'
#' @keywords internal
statmod_edf <- function(spec, coef, design, hyper, expected = TRUE,
                        approx = "bartlett") {
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  H <- NULL
  rows <- list()
  for (a in seq_along(params)) {
    p <- params[a]
    for (nm in names(spec@terms[[p]])) {
      cols <- design[[p]]$blocks[[nm]]
      pen <- modelterms7::term_penalty(spec@terms[[p]][[nm]])
      v <- if (is.null(pen)) {
        as.numeric(length(cols))
      } else {
        if (is.null(H)) {
          H <- statmod_information_at(spec, coef, design, expected, approx)
        }
        idx <- offs[a] + cols
        tryCatch(
          modelterms7::edf(spec@terms[[p]][[nm]],
                           coef = coef[[p]][cols],
                           hessian = H[idx, idx, drop = FALSE],
                           theta = as.list(hyper[[p]][[nm]])),
          error = function(e) {
            # the count is not worth abandoning a fit for, but a message that
            # named none of this would leave the reader to guess why a column
            # went missing
            warning(sprintf("edf() failed for '%s' in '%s': %s", nm, p,
                            conditionMessage(e)), call. = FALSE)
            NA_real_
          })
      }
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = p, term = nm, coefficients = length(cols),
        edf = as.numeric(v)[1L]
      )
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}


#' Resolve the Verbosity Setting
#'
#' @description
#' Turns a level or a named logical vector into the three switches the fit
#' reads.
#'
#' @details
#' The levels name the loops rather than counting them: \code{1} shows the
#' alternation, \code{2} the inner method's own iterations, \code{3} the
#' optimizers' traces. The named form exists because the three are genuinely
#' independent -- watching the alternation while silencing a chatty inner
#' optimizer is the common case.
#'
#' @param verbose A number from 0 to 3, or a named logical vector with any of
#'   \code{blocks}, \code{inner}, \code{optimizer}.
#'
#' @return A list of three logicals.
#'
#' @keywords internal
verbosity <- function(verbose) {
  switches <- c("outer", "blocks", "inner", "optimizer")
  if (is.logical(verbose) && !is.null(names(verbose))) {
    bad <- setdiff(names(verbose), switches)
    if (length(bad)) {
      stop(sprintf(paste0("'verbose' names '%s'. The switches are %s."),
                   bad[1L], paste0("'", switches, "'", collapse = ", ")),
           call. = FALSE)
    }
    # a switch left out is off, so that naming one of the four is a complete
    # request rather than an error
    on <- function(nm) nm %in% names(verbose) && isTRUE(verbose[[nm]])
    return(list(outer = on("outer"), blocks = on("blocks"),
                inner = on("inner"), optimizer = on("optimizer")))
  }
  if (!is.numeric(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be a number from 0 to 3 or a named logical vector.",
         call. = FALSE)
  }
  lv <- as.integer(verbose)
  list(outer = lv >= 1L, blocks = lv >= 1L, inner = lv >= 2L,
       optimizer = lv >= 3L)
}
