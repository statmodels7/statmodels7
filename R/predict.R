#' @include report.R
#' @importFrom stats predict fitted
NULL

#' @title Predict From a Fitted Model
#' @name predict.StatmodFit
#' @description
#' The linear predictors, the distribution's parameters, or the response's
#' mean, at the fitting data or at new data.
#' @details
#' New data goes through the terms' blueprints, so a factor keeps the levels
#' and the contrasts it was fitted with rather than being rebuilt from
#' whatever the new frame happens to contain -- which is the whole reason a
#' term records a blueprint.
#'
#' \code{type = "link"} gives the linear predictors, one vector per
#' distribution parameter; \code{"parameter"} maps each through its inverse
#' link, which is what the density takes; \code{"response"} gives the mean of
#' the distribution at those parameters, and is refused for a family that has
#' none.
#' @param object A \code{\link{StatmodFit}}.
#' @param newdata A data frame, or \code{NULL} for the fitting data.
#' @param type One of \code{"parameter"}, \code{"link"} or \code{"response"}.
#' @param ... Unused.
#' @return A named list of vectors for \code{"link"} and \code{"parameter"},
#'   and a vector for \code{"response"}.
#' @seealso \code{\link{statmod}}
#' @keywords internal
predict.StatmodFit <- function(object, newdata = NULL,
                               type = c("parameter", "link", "response"),
                               ...) {
  type <- match.arg(type)
  spec <- spec_at(object, newdata, need_response = FALSE)
  design <- statmod_design(spec)
  ep <- statmod_eta(spec, design, object@coefficients)
  if (type == "link") return(ep$eta)
  if (type == "parameter") return(ep$theta)
  # `mean` is the base generic distributions7 registers a method on, so it is
  # called unqualified. Only a missing method becomes the friendly message:
  # a catch-all here would report any failure as "this family has no mean",
  # which is the shape of error message the toolkit records as worse than
  # none.
  m <- withCallingHandlers(
    tryCatch(mean(spec@distrib, ep$theta),
             error = function(e) {
               if (grepl("method", conditionMessage(e), fixed = TRUE)) {
                 stop(sprintf(paste0(
                   "'%s' has no mean, so type = \"response\" has no answer.\n",
                   "  Ask for \"parameter\" instead."),
                   spec@distrib@distrib_name), call. = FALSE)
               }
               stop(e)
             }),
    warning = function(w) invokeRestart("muffleWarning"))
  rep_len(m, spec@n_obs)
}
S7::method(predict, StatmodFit) <- predict.StatmodFit


#' @title The Fitted Values of a Model
#' @name fitted.StatmodFit
#' @description The distribution's parameters at the fitting data.
#' @param object A \code{\link{StatmodFit}}.
#' @param ... Unused.
#' @return A named list of vectors.
#' @seealso \code{\link{predict.StatmodFit}}
#' @keywords internal
fitted.StatmodFit <- function(object, ...) object@fitted
S7::method(fitted, StatmodFit) <- fitted.StatmodFit


#' Simulate a Response From a Written Model
#'
#' @description
#' Takes a formula, a distribution and a data frame of covariates, and draws a
#' response from that model -- at coefficients the caller supplies, or at
#' coefficients drawn at random.
#'
#' @details
#' This is not \code{\link[stats]{simulate}}, which draws from a model that has
#' already been fitted. The \code{r} prefix is R's own for a random draw, so
#' the two names cannot be confused.
#'
#' The point of it is to have data whose truth is known: write the model, draw
#' from it, fit it back, and see whether the fit recovers what was put in. A
#' covariate needs no declaring -- a factor becomes its contrasts and a
#' numeric stays itself, because the design comes from the same interpreter a
#' fit uses.
#'
#' \strong{The coefficients.} \code{par = NULL} draws them, each independently
#' from \code{rnorm(1, 0, sd)}, which on the link scale gives predictors of
#' order one. A named list fixes them instead, one vector per distribution
#' parameter in the design's order; a parameter left out of that list is
#' drawn. \code{coef()} on the result reports what was used, drawn or given,
#' so a simulation is reproducible from its own output.
#'
#' \strong{The response's name} is the formula's left-hand side when it is a
#' symbol, and \code{"y"} otherwise.
#'
#' @param formula The model formula, as \code{\link{statmod}} takes it.
#' @param distrib A \pkg{distributions7} distribution object.
#' @param data A data frame of covariates.
#' @param par Optional named list of coefficient vectors.
#' @param sd The standard deviation of the drawn coefficients.
#' @param offsets Optional named list of offsets.
#'
#' @return The data frame with the response added, carrying the attributes
#'   \code{"par"} (the coefficients used) and \code{"theta"} (the parameters
#'   they gave).
#'
#' @seealso \code{\link{statmod}}, \code{\link{predict.StatmodFit}}
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(50), g = factor(rep(c("a", "b"), 25)))
#'
#' # coefficients drawn
#' sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd)
#' attr(sim, "par")
#'
#' # or given, and recovered by a fit
#' sim2 <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
#'                  par = list(mu = c(1, 2), sigma = log(0.3)))
#' statmod(y ~ x, distributions7::gaussian1_distrib(), sim2)
#'
#' @export
rstatmod <- function(formula, distrib, data, par = NULL, sd = 1,
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
  env <- environment(formula)
  if (is.null(env)) env <- baseenv()

  # the response does not exist yet, which is the point, so the terms are
  # built from the covariates alone
  built <- statmod_terms(split$equations, data, env)
  n <- nrow(data)
  design <- lapply(built$terms, function(tms) {
    if (!length(tms)) {
      return(list(X = matrix(0, n, 0L), coef_names = character(0),
                  npar = 0L, blocks = list()))
    }
    mats <- lapply(tms, modelterms7::term_matrix)
    list(X = do.call(cbind, mats),
         coef_names = unlist(lapply(tms, modelterms7::term_coef_names),
                             use.names = FALSE),
         npar = as.integer(sum(vapply(mats, ncol, integer(1)))),
         blocks = list())
  })

  coef <- draw_coefficients(design, params, par, sd)

  off <- check_offsets(offsets, params, n)
  links <- distrib@link_params
  theta <- stats::setNames(vector("list", length(params)), params)
  for (p in params) {
    e <- if (design[[p]]$npar == 0L) rep(0, n) else
      as.numeric(design[[p]]$X %*% coef[[p]])
    if (!is.null(off[[p]])) e <- e + off[[p]]
    theta[[p]] <- linkfunctions7::linkinv(links[[p]], e)
  }

  y <- distributions7::distrib_rng(distrib, n, theta)
  nm <- if (is.name(split$response)) as.character(split$response) else "y"
  out <- data
  out[[nm]] <- y
  attr(out, "par") <- Map(stats::setNames, coef,
                          lapply(design, `[[`, "coef_names"))
  attr(out, "theta") <- theta
  out
}


#' Draw or Validate the Coefficients of a Simulation
#'
#' @description
#' Returns one coefficient vector per distribution parameter, drawn from a
#' normal where the caller gave none.
#'
#' @param design The design blocks.
#' @param params The parameter names.
#' @param par A named list, or \code{NULL}.
#' @param sd The standard deviation of the drawn coefficients.
#'
#' @return A named list of numeric vectors.
#'
#' @keywords internal
draw_coefficients <- function(design, params, par, sd) {
  if (!is.null(par)) {
    if (!is.list(par) || is.null(names(par))) {
      stop("'par' must be a named list, one entry per parameter.",
           call. = FALSE)
    }
    bad <- setdiff(names(par), params)
    if (length(bad)) {
      stop(sprintf(paste0("'par' names '%s', which is not a parameter.\n",
                          "  They are: %s."),
                   bad[1L], paste(params, collapse = ", ")), call. = FALSE)
    }
  }
  stats::setNames(lapply(params, function(p) {
    k <- design[[p]]$npar
    if (!is.null(par) && !is.null(par[[p]])) {
      v <- as.numeric(par[[p]])
      if (length(v) != k) {
        stop(sprintf("'par$%s' has length %d but '%s' has %d coefficients.",
                     p, length(v), p, k), call. = FALSE)
      }
      return(v)
    }
    stats::rnorm(k, 0, sd)
  }), params)
}
