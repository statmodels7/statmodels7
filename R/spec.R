#' @include formula.R
NULL

#' The Specification of a Model, Before It Is Fitted
#'
#' @description
#' Everything the formula, the data and the distribution produce: one
#' equation per distribution parameter, the terms each equation names built
#' against the data, the response, and the prior weights and offsets.
#'
#' @param formula The model formula, as given.
#' @param distrib The \pkg{distributions7} object.
#' @param equations A named list of one-sided formulas, one per parameter, in
#'   the family's order.
#' @param terms A named list, one entry per parameter, each a named list of
#'   built \pkg{modelterms7} terms.
#' @param response The evaluated left-hand side.
#' @param n_obs The number of observations.
#' @param weights Prior weights, one per observation.
#' @param offsets A named list of offsets, one per parameter or \code{NULL}.
#' @param intercepts A named logical, whether each equation carried one.
#'
#' @return An object of class \code{StatmodSpec}.
#'
#' @seealso \code{\link{statmod_spec}}
#'
#' @examples
#' dd <- data.frame(y = rnorm(10), x = runif(10))
#' S7::S7_inherits(statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd),
#'                 StatmodSpec)
#'
#' @name StatmodSpec-class
#' @aliases StatmodSpec
#' @keywords internal
#' @export
StatmodSpec <- S7::new_class("StatmodSpec",
  properties = list(
    formula = S7::class_any,
    distrib = S7::class_any,
    equations = S7::class_list,
    terms = S7::class_list,
    response = S7::class_any,
    n_obs = S7::class_integer,
    weights = S7::class_numeric,
    offsets = S7::class_list,
    intercepts = S7::class_logical,
    newdata = S7::class_any
  )
)


#' Build a Model Specification
#'
#' @description
#' Splits the formula into one equation per distribution parameter, interprets
#' each with \pkg{modelterms7} and builds its terms against the data.
#'
#' @details
#' The equations are interpreted in an environment where \pkg{modelterms7}'s
#' term constructors shadow whatever the user has attached, so that \code{s()}
#' means ours even with \pkg{mgcv} on the search path. See
#' \code{\link{statmod_equations}} for the split itself, which is not the
#' obvious one.
#'
#' Prior weights enter the log-likelihood as \eqn{\sum_i w_i \ell_i} and are
#' taken as given. They are deliberately NOT normalized: dividing by their sum
#' would turn the log-likelihood into a mean, shrinking every standard error
#' by \eqn{\sqrt{n}} and making the information criteria incomparable with an
#' unweighted fit of the same model.
#'
#' @param formula The model formula.
#' @param distrib A \pkg{distributions7} distribution object.
#' @param data A data frame.
#' @param weights Optional prior weights, one per observation.
#' @param offsets Optional named list of offsets, one per parameter.
#' @param need_response Whether the left-hand side must evaluate. A likelihood
#'   needs it; a prediction does not, and new data routinely has no response
#'   column.
#'
#' @return An object of class \code{\link{StatmodSpec}}.
#'
#' @seealso \code{\link{statmod_equations}}, \code{\link{statmod}}
#'
#' @examples
#' dd <- data.frame(y = rnorm(20), x = runif(20), z = runif(20))
#' spec <- statmod_spec(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd)
#' names(spec@terms)
#'
#' @export
statmod_spec <- function(formula, distrib, data, weights = NULL,
                         offsets = NULL, need_response = TRUE) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!S7::S7_inherits(distrib, distributions7::distrib)) {
    stop("'distrib' must be a distributions7 distribution object.",
         call. = FALSE)
  }
  params <- distrib@params
  split <- statmod_equations(formula, params)

  # the response is evaluated once, in the data, exactly as
  # interpret_formula() would evaluate it for a single-parameter model
  env <- environment(formula)
  if (is.null(env)) env <- baseenv()
  # prediction needs the design and not the response, and new data routinely
  # has no response column; a likelihood needs it and says so
  response <- tryCatch(eval(split$response, data, env),
                       error = function(e) if (need_response) stop(e) else NULL)
  if (is.null(response)) {
    response <- rep(NA_real_, nrow(data))
  }
  n <- if (is.matrix(response)) nrow(response) else length(response)
  if (n == 0L) stop("The response is empty.", call. = FALSE)

  built <- statmod_terms(split$equations, data, env)
  terms_by_param <- built$terms
  intercepts <- built$intercepts

  weights <- check_weights(weights, n)
  offsets <- check_offsets(offsets, params, n)

  StatmodSpec(
    formula = formula, distrib = distrib,
    equations = split$equations, terms = terms_by_param,
    response = response, n_obs = as.integer(n),
    weights = weights, offsets = offsets, intercepts = intercepts,
    newdata = NULL
  )
}


#' The Same Model Read on Other Rows
#'
#' @description
#' A specification carrying the fitted terms and a new data frame, so that
#' every block is reapplied to those rows rather than rebuilt from them.
#'
#' @details
#' A term records how its block was made -- a factor's levels and contrasts, a
#' spline's knots, a basis reparametrization -- and
#' \code{\link[modelterms7]{term_predict}} reapplies that record. Rebuilding
#' instead gives a block of the same shape, multiplying the same coefficients,
#' that means something else: measured on \code{y ~ s(x, k = 10)} at 200
#' observations, predicting on 40 of the rows the model was fitted to differed
#' from the fitted values there by 0.237, and on the 51 rows with
#' \eqn{|x| < 0.5}, where the rebuilt knots move furthest, by 1.19. The whole
#' data handed back agrees exactly, which is why nothing noticed.
#'
#' @param spec The fitted \code{\link{StatmodSpec}}.
#' @param data The rows to read it on.
#' @param need_response Whether the response has to be there.
#'
#' @return A \code{\link{StatmodSpec}} whose \code{newdata} is set.
#'
#' @seealso \code{\link{statmod_design}}
#'
#' @keywords internal
statmod_respec <- function(spec, data, need_response = TRUE) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  env <- environment(spec@formula)
  if (is.null(env)) env <- baseenv()
  split <- statmod_equations(spec@formula, spec@distrib@params)
  response <- tryCatch(eval(split$response, data, env),
                       error = function(e) if (need_response) stop(e) else NULL)
  if (is.null(response)) response <- rep(NA_real_, nrow(data))
  n <- if (is.matrix(response)) nrow(response) else length(response)
  if (n == 0L) stop("The response is empty.", call. = FALSE)
  S7::set_props(spec, response = response, n_obs = as.integer(n),
                weights = rep(1, n),
                offsets = check_offsets(NULL, spec@distrib@params, n),
                newdata = data)
}


#' Interpret and Build Each Parameter's Terms
#'
#' @description
#' Runs \pkg{modelterms7}'s interpreter on every equation and builds the terms
#' it names against the data.
#'
#' @details
#' The equations are interpreted with \pkg{modelterms7}'s constructors in front
#' of the search path, so that \code{s()} means ours whatever the user has
#' attached. A factor covariate needs no special handling: the interpreter
#' collects bare covariates into one \code{linpar()}, whose block comes from
#' \code{model.matrix} and therefore carries the contrasts.
#'
#' @param equations A named list of one-sided formulas.
#' @param data A data frame.
#' @param env The environment the original formula carried.
#'
#' @return A list with \code{terms} (a named list per parameter) and
#'   \code{intercepts} (a named logical).
#'
#' @keywords internal
statmod_terms <- function(equations, data, env) {
  shim <- terms_first(env)
  params <- names(equations)
  out_terms <- stats::setNames(vector("list", length(params)), params)
  intercepts <- stats::setNames(logical(length(params)), params)
  for (p in params) {
    eq <- equations[[p]]
    environment(eq) <- shim
    out <- modelterms7::interpret_formula(eq, data)
    intercepts[[p]] <- out$intercept
    out_terms[[p]] <- lapply(out$terms, function(tm)
      modelterms7::term_build(tm, data))
    names(out_terms[[p]]) <- names(out$terms)
  }
  reject_unfittable(out_terms)
  list(terms = out_terms, intercepts = intercepts)
}


#' Reject a Term the Fitting Scheme Does Not Cover
#'
#' @description
#' Signals an error naming any term whose block is not a fixed design, which
#' is what the alternation of \code{\link{statmod}} assembles.
#'
#' @details
#' Two shapes are outside that assembly, and both are read off the term rather
#' than from a list of class names, so a term written later is covered without
#' an edit here.
#'
#' A \strong{structural} term rewrites the likelihood instead of contributing a
#' predictor, so it has no design block at all and answers neither
#' \code{term_matrix()} nor \code{term_npar()}. Reaching it through the design
#' produced an error naming one of those generics, which says nothing about the
#' cause.
#'
#' A term whose block \strong{depends on its own coefficients} registers a
#' \code{\link[modelterms7]{term_refresh}} method of its own, the base method
#' on \code{model_term} being the identity; the class a method was registered
#' on is \code{attr(m, "signature")[[1]]}. For those the block is a Jacobian
#' and the working solution is an increment, so assembling it once and solving
#' for the coefficients estimates something else while reporting convergence:
#' measured on \code{seg()}, the break-point stays at its starting value and
#' the fitted mean of a continuous construction carries a step. Rejecting them
#' is what keeps that out of a returned object until the alternation refreshes
#' a block between inner fits.
#'
#' Every equation is examined before the error is raised, so a model carrying
#' one such term in the mean and another in the scale reports both rather than
#' the first.
#'
#' @param terms The built terms, a named list of named lists, one per
#'   distribution parameter.
#'
#' @return \code{NULL}, invisibly; called for the error.
#'
#' @seealso \code{\link{statmod_terms}}, \code{\link{statmod}}
#'
#' @keywords internal
reject_unfittable <- function(terms) {
  found <- character(0)
  for (p in names(terms)) {
    for (nm in names(terms[[p]])) {
      why <- unfittable_reason(terms[[p]][[nm]])
      if (nzchar(why)) {
        found <- c(found, sprintf("'%s' in '%s': %s", nm, p, why))
      }
    }
  }
  if (length(found) == 0L) return(invisible(NULL))
  stop(if (length(found) == 1L) {
    paste0("statmod() cannot fit ", found)
  } else {
    paste0("statmod() cannot fit these terms:\n",
           paste0("  ", found, collapse = "\n"))
  }, call. = FALSE)
}


#' Why a Term Is Outside the Fitting Scheme
#'
#' @description
#' Returns the reason a term cannot be assembled as a fixed design block, or
#' the empty string when it can.
#'
#' @param term One built term.
#'
#' @return A single string.
#'
#' @seealso \code{\link{reject_unfittable}}
#'
#' @keywords internal
unfittable_reason <- function(term) {
  if (S7::S7_inherits(term, modelterms7::structural_term)) {
    return(paste("a structural term rewrites the likelihood instead of",
                 "contributing a predictor, and the fitting scheme has no",
                 "route for one yet."))
  }
  if (refreshes_own_block(term)) {
    return(paste("its block depends on its own coefficients and the",
                 "alternation does not refresh a block between inner fits",
                 "yet, so fitting it here would hold the break-point or the",
                 "nonlinear parameters at their starting values and report",
                 "convergence. Iterate it with modelterms7::term_refresh()",
                 "instead; see its help page."))
  }
  ""
}


#' Does a Term Recompute Its Own Block?
#'
#' @description
#' \code{TRUE} when the term registers a
#' \code{\link[modelterms7]{term_refresh}} method of its own rather than
#' inheriting the identity registered on \code{model_term}.
#'
#' @details
#' The owning class of a method is \code{attr(m, "signature")[[1]]}, and it is
#' compared by name and package rather than by \code{identical()}: an S7 class
#' re-created from the same definition is not identical to the original, which
#' is what happens whenever a package's code is re-evaluated rather than
#' loaded.
#'
#' @param term One built term.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{unfittable_reason}}
#'
#' @keywords internal
refreshes_own_block <- function(term) {
  m <- tryCatch(S7::method(modelterms7::term_refresh, S7::S7_class(term)),
                error = function(e) NULL)
  if (is.null(m)) return(FALSE)
  owner <- attr(m, "signature")[[1]]
  base <- modelterms7::model_term
  !(identical(attr(owner, "name"), attr(base, "name")) &&
    identical(attr(owner, "package"), attr(base, "package")))
}


#' Validate Prior Weights
#'
#' @description
#' Returns a vector of prior weights of the right length, defaulting to one
#' per observation.
#'
#' @details
#' They are not normalized. Making them sum to one would turn
#' \eqn{\sum_i w_i \ell_i} into a mean, which is the averaging trap the
#' objective's own scale already has to avoid: every standard error would
#' shrink by \eqn{\sqrt{n}} and the information criteria would stop being
#' comparable.
#'
#' @param weights The weights, or \code{NULL}.
#' @param n The number of observations.
#'
#' @return A numeric vector of length \code{n}.
#'
#' @keywords internal
check_weights <- function(weights, n) {
  if (is.null(weights)) return(rep(1, n))
  if (!is.numeric(weights)) {
    stop("'weights' must be numeric.", call. = FALSE)
  }
  if (length(weights) != n) {
    stop(sprintf("'weights' has length %d but there are %d observations.",
                 length(weights), n), call. = FALSE)
  }
  if (anyNA(weights) || any(weights < 0)) {
    stop("'weights' must be non-negative and complete.", call. = FALSE)
  }
  if (all(weights == 0)) {
    stop("'weights' are all zero: there is nothing to fit.", call. = FALSE)
  }
  as.numeric(weights)
}


#' Validate Offsets
#'
#' @description
#' Returns a named list with one offset per parameter, \code{NULL} where none
#' was given.
#'
#' @param offsets A named list, or \code{NULL}.
#' @param params The parameter names.
#' @param n The number of observations.
#'
#' @return A named list of length \code{length(params)}.
#'
#' @keywords internal
check_offsets <- function(offsets, params, n) {
  out <- stats::setNames(vector("list", length(params)), params)
  if (is.null(offsets)) return(out)
  if (!is.list(offsets) || is.null(names(offsets))) {
    stop("'offsets' must be a named list, one entry per parameter.",
         call. = FALSE)
  }
  bad <- setdiff(names(offsets), params)
  if (length(bad)) {
    stop(sprintf(paste0("'%s' is not a parameter of this distribution.\n",
                        "  Its parameters are: %s."),
                 bad[1L], paste(params, collapse = ", ")), call. = FALSE)
  }
  for (p in names(offsets)) {
    o <- offsets[[p]]
    if (!is.numeric(o) || (length(o) != 1L && length(o) != n)) {
      stop(sprintf("The offset for '%s' must be numeric of length 1 or %d.",
                   p, n), call. = FALSE)
    }
    out[[p]] <- rep_len(as.numeric(o), n)
  }
  out
}


#' The Design of a Specification
#'
#' @description
#' Returns, per parameter, the terms' blocks side by side and the names of the
#' coefficients they carry.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#'
#' @return A named list with one entry per parameter, each a list with
#'   \code{X}, \code{coef_names}, \code{npar} and \code{blocks} (the column
#'   range each term occupies).
#'
#' @seealso \code{\link{statmod_spec}}
#'
#' @examples
#' dd <- data.frame(y = rnorm(20), x = runif(20))
#' d <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
#' vapply(statmod_design(d), function(z) z$npar, integer(1))
#'
#' @export
statmod_design <- function(spec) {
  lapply(spec@terms, function(tms) {
    if (!length(tms)) {
      return(list(X = matrix(0, spec@n_obs, 0L), coef_names = character(0),
                  npar = 0L, blocks = list()))
    }
    # on other rows every block is REAPPLIED, never rebuilt: a term records
    # how its block was made and term_predict() replays that record, where
    # building again would give a block of the same shape multiplying the same
    # coefficients and meaning something else
    mats <- if (is.null(spec@newdata)) lapply(tms, modelterms7::term_matrix)
      else lapply(tms, modelterms7::term_predict, newdata = spec@newdata)
    nms <- lapply(tms, modelterms7::term_coef_names)
    widths <- vapply(mats, ncol, integer(1))
    ends <- cumsum(widths)
    starts <- ends - widths + 1L
    list(
      X = do.call(cbind, mats),
      coef_names = unlist(nms, use.names = FALSE),
      npar = as.integer(sum(widths)),
      blocks = stats::setNames(
        lapply(seq_along(tms), function(j)
          if (widths[j] == 0L) integer(0) else seq.int(starts[j], ends[j])),
        names(tms))
    )
  })
}
