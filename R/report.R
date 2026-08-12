#' @include statmod.R
#' @importFrom stats logLik coef
NULL

#' Format a Duration in the Unit It Deserves
#'
#' @description
#' Turns a number of seconds into a string in microseconds, milliseconds,
#' seconds, minutes, hours or days, whichever keeps it readable.
#'
#' @details
#' A fit that took 340 microseconds and one that took 2.6 hours must both read
#' at a glance, which a single unit cannot do.
#'
#' Zero is a reading and not the absence of one: on a coarse clock a fast fit
#' measures exactly zero seconds, so it is reported as such rather than
#' suppressed. A guard of the form \code{if (elapsed > 0)} would be a claim
#' that zero cannot be measured, and for a duration that claim is false.
#'
#' @param seconds A number of seconds.
#' @param digits Significant digits.
#'
#' @return A single string.
#'
#' @examples
#' format_duration(c(3.4e-4, 0.25, 90, 7200))
#'
#' @export
format_duration <- function(seconds, digits = 3L) {
  vapply(seconds, function(s) {
    if (!is.finite(s)) return("NA")
    if (s < 0) return("negative")
    if (s == 0) return("0 s")
    if (s < 1e-3) return(paste0(signif(s * 1e6, digits), " us"))
    if (s < 1) return(paste0(signif(s * 1e3, digits), " ms"))
    if (s < 60) return(paste0(signif(s, digits), " s"))
    if (s < 3600) return(paste0(signif(s / 60, digits), " min"))
    if (s < 86400) return(paste0(signif(s / 3600, digits), " h"))
    paste0(signif(s / 86400, digits), " d")
  }, character(1))
}


#' @title Print a Fitted Model
#' @name print.StatmodFit
#' @description
#' The call, the distribution, one line per equation with its terms and what
#' they spend, the log-likelihood, the elapsed time and whether every loop
#' stopped on its own rule.
#' @param x A \code{\link{StatmodFit}}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{statmod}}
#' @keywords internal
print.StatmodFit <- function(x, ...) {
  spec <- x@spec
  cat("A statmod fit\n\n")
  cat("Call:  ", paste(deparse(x@call), collapse = "\n        "), "\n\n",
      sep = "")
  cat("Distribution: ", spec@distrib@distrib_name, "\n", sep = "")
  cat("Observations: ", spec@n_obs, sep = "")
  sw <- sum(spec@weights)
  if (!isTRUE(all.equal(sw, spec@n_obs))) {
    cat(sprintf("  (prior weights summing to %s, not %d)",
                format(signif(sw, 6)), spec@n_obs))
  }
  cat("\n\n")

  for (p in spec@distrib@params) {
    rhs <- paste(deparse(spec@equations[[p]][[2L]]), collapse = " ")
    cat(sprintf("  %-10s ~ %s\n", p, rhs))
    if (!is.null(x@edf)) {
      rows <- x@edf[x@edf$parameter == p, , drop = FALSE]
      for (i in seq_len(nrow(rows))) {
        cat(sprintf("  %-10s   %-14s %3d coef", "", rows$term[i],
                    rows$coefficients[i]))
        if (is.finite(rows$edf[i]) &&
            !isTRUE(all.equal(rows$edf[i], rows$coefficients[i]))) {
          cat(sprintf(", edf %.2f", rows$edf[i]))
        }
        cat("\n")
      }
    }
  }

  cat("\n")
  cat(sprintf("log-likelihood %.6f    objective %.6f\n",
              x@loglik, x@objective))
  if (!is.null(x@methods$outer) && length(x@criterion) &&
      is.finite(x@criterion)) {
    cat(sprintf("%s %.6f over %d hyperparameter evaluation(s)\n",
                toupper(x@methods$outer@kind), x@criterion,
                if (is.null(x@history$outer)) 0L else nrow(x@history$outer)))
  }
  cat(sprintf("fitted in %s, %s\n", format_duration(x@elapsed),
              if (x@converged) "converged" else "DID NOT CONVERGE"))
  if (!x@converged) cat(fitted_ranges(x))
  if (!is.null(x@history$blocks) && nrow(x@history$blocks) > 1L) {
    cat(sprintf("%d sweeps over %d block(s)\n",
                max(x@history$blocks$sweep),
                length(unique(x@history$blocks$block))))
  }
  invisible(x)
}
S7::method(print, StatmodFit) <- print.StatmodFit


#' The Model as a Function of Parameters and Data
#'
#' @description
#' \code{loglik()}, \code{gradient()} and \code{hessian()} evaluate a fitted
#' model at parameters and data of the caller's choosing.
#'
#' @details
#' The name is \code{loglik} and not \code{logLik}: R's \code{logLik()}
#' returns the maximized value of the fitted model and carries \code{df} and
#' \code{nobs}, and overloading it would give one name two behaviours. Both
#' exist, and \code{loglik(fit)} with no arguments must equal
#' \code{logLik(fit)} to the last digit, which is the cheapest check that the
#' callable route and the fitting route are the same model.
#'
#' They are generics rather than closures stored in the fit. A closure
#' captures its environment, which means the data: a fit would then carry the
#' frame twice and keep a stale copy after the data changed. These rebuild
#' from the specification the fit keeps whole, running the terms' blueprints
#' against new data by the same path \code{predict()} takes.
#'
#' @param object A \code{\link{StatmodFit}}.
#' @param par A named list of coefficient vectors, one per distribution
#'   parameter. Defaults to the fitted ones.
#' @param data A data frame. Defaults to the data the model was fitted to.
#' @param ... Passed to methods.
#'
#' @return A number for \code{loglik}, a named list of vectors for
#'   \code{gradient}, a matrix for \code{hessian}.
#'
#' @seealso \code{\link{statmod}}
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.5)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' loglik(fit)
#' loglik(fit, par = list(mu = c(0, 0), sigma = 0))
#'
#' @export
loglik <- S7::new_generic("loglik", "object",
  function(object, par = NULL, data = NULL, ...) S7::S7_dispatch())

#' @rdname loglik
#' @export
gradient <- S7::new_generic("gradient", "object",
  function(object, par = NULL, data = NULL, ...) S7::S7_dispatch())

#' @rdname loglik
#' @param expected Whether the expected information is wanted.
#' @export
hessian <- S7::new_generic("hessian", "object",
  function(object, par = NULL, data = NULL, expected = FALSE, ...)
    S7::S7_dispatch())


#' Rebuild a Specification Against New Data
#'
#' @description
#' Returns the fit's specification when \code{data} is \code{NULL}, and one
#' built against the new data otherwise.
#'
#' @param fit A \code{\link{StatmodFit}}.
#' @param data A data frame or \code{NULL}.
#' @param need_response Whether the response must be there. A likelihood needs
#'   it; a prediction does not, and new data routinely has no response column.
#'
#' @return A \code{\link{StatmodSpec}}.
#'
#' @keywords internal
spec_at <- function(fit, data, need_response = TRUE) {
  if (is.null(data)) return(fit@spec)
  statmod_respec(fit@spec, data, need_response = need_response)
}

#' Resolve a Parameter Structure
#'
#' @description
#' Returns the fitted coefficients when \code{par} is \code{NULL}, and
#' validates a supplied structure against the design otherwise.
#'
#' @param fit A \code{\link{StatmodFit}}.
#' @param par A named list or \code{NULL}.
#' @param design The design to validate against.
#'
#' @return A named list of coefficient vectors.
#'
#' @keywords internal
par_at <- function(fit, par, design) {
  if (is.null(par)) return(fit@coefficients)
  params <- fit@spec@distrib@params
  bad <- setdiff(names(par), params)
  if (length(bad)) {
    stop(sprintf(paste0("'par' names '%s', which is not a parameter.\n",
                        "  They are: %s."),
                 bad[1L], paste(params, collapse = ", ")), call. = FALSE)
  }
  out <- fit@coefficients
  for (p in names(par)) {
    if (length(par[[p]]) != design[[p]]$npar) {
      stop(sprintf("'par$%s' has length %d but '%s' has %d coefficients.",
                   p, length(par[[p]]), p, design[[p]]$npar), call. = FALSE)
    }
    out[[p]] <- as.numeric(par[[p]])
  }
  out
}

#' @rdname loglik
#' @name loglik.StatmodFit
#' @keywords internal
S7::method(loglik, StatmodFit) <- function(object, par = NULL, data = NULL,
                                           ...) {
  spec <- spec_at(object, data)
  design <- statmod_design(spec)
  statmod_loglik_at(spec, par_at(object, par, design), design)
}

#' @rdname loglik
#' @name gradient.StatmodFit
#' @keywords internal
S7::method(gradient, StatmodFit) <- function(object, par = NULL, data = NULL,
                                             ...) {
  spec <- spec_at(object, data)
  design <- statmod_design(spec)
  statmod_score_at(spec, par_at(object, par, design), design)
}

#' @rdname loglik
#' @name hessian.StatmodFit
#' @keywords internal
S7::method(hessian, StatmodFit) <- function(object, par = NULL, data = NULL,
                                            expected = FALSE, ...) {
  spec <- spec_at(object, data)
  design <- statmod_design(spec)
  # the information is minus the Hessian, and what is asked for here is the
  # Hessian of the log-likelihood itself
  -statmod_information_at(spec, par_at(object, par, design), design, expected)
}


#' @title The Maximized Log-Likelihood of a Fit
#' @name logLik.StatmodFit
#' @description
#' The value R's convention expects, carrying the degrees of freedom and the
#' number of observations. \code{\link{loglik}} is the other thing: the model
#' evaluated at parameters and data of the caller's choosing.
#' @details
#' \strong{Which likelihood, and it matters.} The default is the
#' CONDITIONAL one: the log-density at the fitted coefficients, a penalized
#' coefficient among them, paired with the effective degrees of freedom
#' \eqn{\mathrm{tr}[(H+S)^{-1}H]}. A criterion built on the pair is the
#' conditional AIC, and its question is how well the model describes the
#' groups, curves and states actually observed.
#'
#' A mixed-model package reports the MARGINAL likelihood instead: the random
#' effects are integrated out and the count is the number of estimated
#' parameters, variance components among them. Its question is about the
#' population, and its AIC is a different number that is not comparable with
#' the conditional one. Mixing the two -- a marginal likelihood against an
#' effective count, or the reverse -- is what neither convention allows
#' (Vaida and Blanchard, 2005).
#'
#' \code{type = "marginal"} returns the value the outer criterion evaluated
#' while choosing the hyperparameters, with the number of estimated
#' parameters as its degrees of freedom. It is available only where a
#' marginal criterion actually ran: \code{\link{ml}} or \code{\link{reml}}.
#' Where the hyperparameters were held, or found by a prediction criterion,
#' there is no marginal likelihood to report and asking for one is an error
#' rather than a number that would look like one.
#' @param object A \code{\link{StatmodFit}}.
#' @param type \code{"conditional"} (default) or \code{"marginal"}.
#' @param ... Unused.
#' @return A \code{logLik} object.
#' @references
#' Vaida, F. and Blanchard, S. (2005). Conditional Akaike information for
#' mixed-effects models. \emph{Biometrika}, 92(2), 351--370.
#' @seealso \code{\link{loglik}}
#' @keywords internal
logLik.StatmodFit <- function(object,
                              type = c("conditional", "marginal"), ...) {
  type <- match.arg(type)
  if (identical(type, "marginal")) {
    outer <- object@methods$outer
    if (is.null(outer) || !outer@kind %in% c("ml", "reml")) {
      stop(paste0(
        "the marginal log-likelihood is the value ml() or reml() optimizes,",
        "\n  and this fit used ",
        if (is.null(outer)) "held hyperparameters" else
          sprintf("%s()", outer@kind),
        ". There is no marginal likelihood to\n  report; logLik(type = ",
        "\"conditional\") is what this fit has."), call. = FALSE)
    }
    if (!length(object@criterion) || !is.finite(object@criterion)) {
      stop("the outer criterion was not recorded for this fit.",
           call. = FALSE)
    }
    # The parameters of the MARGINAL model, which is a different count from
    # the conditional one and not an effective anything: the coefficients it
    # did not integrate away, one per hyperparameter it estimated, and a
    # structural term's own parameters, which are estimated rather than
    # integrated. A penalized coefficient is a random effect under another
    # name and is integrated out, so it does not appear.
    design <- statmod_design(object@spec)
    lab <- coef_labels(object@spec, design)
    nh <- length(unlist(object@hyper, use.names = FALSE))
    ns <- sum(vapply(statmod_structural_par(object@spec, design),
                     function(u) length(u$parameter) - length(u$held),
                     numeric(1)))
    return(structure(as.numeric(object@criterion),
                     df = sum(!lab$penalized) + nh + ns,
                     nobs = object@spec@n_obs, class = "logLik"))
  }
  df <- if (is.null(object@edf)) {
    sum(lengths(object@coefficients))
  } else {
    # a term whose count could not be obtained falls back to its number of
    # columns, which is the upper bound: dropping it would understate the
    # dimension of the model and flatter every criterion built on it
    e <- object@edf
    sum(ifelse(is.na(e$edf), e$coefficients, e$edf))
  }
  structure(object@loglik, df = df, nobs = object@spec@n_obs,
            class = "logLik")
}
S7::method(logLik, StatmodFit) <- logLik.StatmodFit

#' @title The Coefficients of a Fit
#' @name coef.StatmodFit
#' @description One vector per distribution parameter.
#' @param object A \code{\link{StatmodFit}}.
#' @param ... Unused.
#' @return A named list of numeric vectors.
#' @seealso \code{\link{statmod}}
#' @keywords internal
coef.StatmodFit <- function(object, ...) {
  params <- object@spec@distrib@params
  design <- statmod_design(object@spec)
  stats::setNames(lapply(params, function(p)
    stats::setNames(object@coefficients[[p]], design[[p]]$coef_names)), params)
}
S7::method(coef, StatmodFit) <- coef.StatmodFit


#' What Each Distribution Parameter Reached
#'
#' @description
#' One line per parameter of the distribution, giving the range of its fitted
#' values, for a fit that did not converge.
#'
#' @details
#' The measured case this exists for: a lasso at a fixed hyperparameter, with a
#' free scale and a design the model can interpolate. Fitting the coefficients
#' shrinks the residuals, which shrinks the scale, which raises the working
#' weights, which makes the penalty count for relatively less, which lets more
#' coefficients in. At 200 observations and 400 columns the scale reached
#' 3.8e-15 and 380 of the 400 coefficients survived, where the same block
#' fitted at a held scale kept the five that were real.
#'
#' Nothing here diagnoses that. It reports where the parameters ended up, which
#' is a fact, and a scale at 1e-15 says the rest on its own. Naming a cause
#' would mean picking a threshold for what counts as running away, and the same
#' fit at 100 columns converges to a scale of 0.77 that is nothing of the kind.
#'
#' @param x A \code{\link{StatmodFit}}.
#'
#' @return A single string, empty when the parameters cannot be read.
#'
#' @seealso \code{\link{statmod}}
#'
#' @keywords internal
fitted_ranges <- function(x) {
  th <- tryCatch(x@fitted, error = function(e) NULL)
  if (!length(th)) return("")
  rows <- vapply(names(th), function(p) {
    v <- as.numeric(th[[p]])
    v <- v[is.finite(v)]
    if (!length(v)) return(sprintf("  %-10s not finite", p))
    if (length(unique(v)) == 1L) {
      sprintf("  %-10s %s", p, format(signif(v[1L], 4)))
    } else {
      sprintf("  %-10s %s to %s", p, format(signif(min(v), 4)),
              format(signif(max(v), 4)))
    }
  }, character(1))
  paste0("  the parameters it reached:\n", paste(rows, collapse = "\n"), "\n")
}
