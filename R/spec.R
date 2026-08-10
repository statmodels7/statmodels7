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
    intercepts = S7::class_logical
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
                         offsets = NULL) {
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
  response <- eval(split$response, data, env)
  n <- if (is.matrix(response)) nrow(response) else length(response)
  if (n == 0L) stop("The response is empty.", call. = FALSE)

  # each equation is interpreted with our terms in front of the search path
  shim <- terms_first(env)
  terms_by_param <- stats::setNames(vector("list", length(params)), params)
  intercepts <- stats::setNames(logical(length(params)), params)
  for (p in params) {
    eq <- split$equations[[p]]
    environment(eq) <- shim
    out <- modelterms7::interpret_formula(eq, data)
    intercepts[[p]] <- out$intercept
    terms_by_param[[p]] <- lapply(out$terms, function(tm)
      modelterms7::term_build(tm, data))
    names(terms_by_param[[p]]) <- names(out$terms)
  }

  weights <- check_weights(weights, n)
  offsets <- check_offsets(offsets, params, n)

  StatmodSpec(
    formula = formula, distrib = distrib,
    equations = split$equations, terms = terms_by_param,
    response = response, n_obs = as.integer(n),
    weights = weights, offsets = offsets, intercepts = intercepts
  )
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
    mats <- lapply(tms, modelterms7::term_matrix)
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
