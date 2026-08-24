#' @include statmod.R
#' @importFrom stats logLik coef residuals
NULL

#' Format a Duration in the Unit It Deserves
#'
#' @description
#' Turns a number of seconds into a string in microseconds, milliseconds,
#' seconds, minutes, hours or days, whichever keeps it readable.
#'
#' @details
#' A fit that took 340 microseconds and one that took 2.6 hours must both
#' read at a glance, and no single unit does that.
#'
#' The thresholds are the obvious ones: below a millisecond it reads in
#' microseconds, below a second in milliseconds, below a minute in seconds,
#' then minutes, hours and days.
#'
#' Zero is a reading. On a coarse clock a fast fit measures exactly zero
#' seconds, and it is reported as `"0 s"` rather than suppressed. A guard of
#' the form `if (elapsed > 0)` would claim that zero cannot be measured,
#' which for a duration is false, and would make what the object prints
#' depend on the platform's timer resolution.
#'
#' @param seconds A number of seconds. Vectorized: a vector in gives a vector
#'   out, each element in its own unit.
#' @param digits How many significant digits to keep, passed to [signif()].
#'   Defaults to 3.
#'
#' @return A character vector the length of `seconds`, each element a number
#'   and a unit separated by a space, as `"340 us"` or `"2.6 h"`.
#'
#' @seealso [print.StatmodFit()], which reports a fit's elapsed time with
#'   this.
#'
#' @examples
#' format_duration(c(3.4e-4, 0.25, 90, 7200))
#'
#' # Zero is a reading, not a missing value.
#' format_duration(0)
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
#' Prints a compact report of a fitted model: the call, the distribution, one
#' line per equation naming its terms and the effective degrees of freedom
#' each spends, the conditional log-likelihood, the elapsed time and whether
#' every loop stopped on its own rule.
#'
#' [summary.StatmodFit()] is the long form, with coefficients, standard
#' errors and the hyperparameters.
#' @param x A [StatmodFit()].
#' @param ... Unused.
#' @return `x`, invisibly, as a print method should.
#' @seealso [summary.StatmodFit()] for the full report,
#'   [statmod_certificate()] for the verdict behind the convergence line.
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

  w <- if (is.null(x@edf) || !nrow(x@edf)) 0L else max(nchar(x@edf$term))
  for (p in spec@distrib@params) {
    rhs <- paste(deparse(spec@equations[[p]][[2L]]), collapse = " ")
    cat(sprintf("  %-10s ~ %s\n", p, rhs))
    if (!is.null(x@edf)) {
      rows <- x@edf[x@edf$parameter == p, , drop = FALSE]
      for (i in seq_len(nrow(rows))) {
        tm <- spec@terms[[p]][[rows$term[i]]]
        st <- !is.null(tm) &&
          identical(term_block_kind(tm), "structural")
        k <- if (st) length(modelterms7::term_params(tm)) else
          rows$coefficients[i]
        cat(sprintf("  %-10s   %s %3d %s", "",
                    formatC(rows$term[i], width = -w, flag = " "), k,
                    if (st) "param" else "coef"))
        if (is.finite(rows$edf[i]) && !isTRUE(all.equal(rows$edf[i], k))) {
          cat(sprintf(", edf %.2f", rows$edf[i]))
        }
        cat("\n")
      }
    }
  }

  cat("\n")
  # NO LIKELIHOOD HERE. There are two of them and they answer different
  # questions -- the conditional one the criteria are built on, and the
  # penalized one the inner fit minimizes -- so a single line either carried
  # a number whose meaning changed with the model or invited the two views
  # to be read against each other. `summary()` reports the conditional one
  # beside the effective degrees of freedom and the criteria, which is where
  # the pairing means something; this view says what the model IS.
  if (!is.null(x@methods$outer) && length(x@criterion) &&
      is.finite(x@criterion)) {
    cat(sprintf("%s %.6f over %d hyperparameter evaluation(s)\n",
                toupper(x@methods$outer@kind), x@criterion,
                if (is.null(x@history$outer)) 0L else nrow(x@history$outer)))
  }
  # the same distinction summary() draws: this is the SEARCH's own report and
  # not a verdict on the point, which is what statmod_certificate() answers.
  # The word is kept rather than replaced by "met its stopping rule", which
  # is more precise and which nobody greps for.
  cat(sprintf("fitted in %s   search: %s\n", format_duration(x@elapsed),
              if (x@converged) "converged" else "DID NOT CONVERGE"))
  if (!is.null(x@history$blocks) && nrow(x@history$blocks) > 1L) {
    cat(sprintf("%d pass(es) over %d block(s)\n",
                max(x@history$blocks$pass),
                length(unique(x@history$blocks$block))))
  }
  invisible(x)
}
S7::method(print, StatmodFit) <- print.StatmodFit


#' The Model as a Function of Parameters and Data
#'
#' @description
#' Turn a fitted model back into a function. `loglik()`, `gradient()` and
#' `hessian()` evaluate the model's log-likelihood and its first two
#' derivatives at coefficients and data of the caller's choosing, so a fit
#' can be profiled, differenced or handed to an optimizer of one's own.
#'
#' With no arguments beyond the fit they read at the fitted coefficients and
#' the fitting data.
#'
#' @details
#' # The name
#'
#' `loglik`, not `logLik`. R's [stats::logLik()] returns the maximized value
#' of the fitted model and carries `df` and `nobs` with it, and overloading
#' it would give one name two behaviors. Both exist here, and
#' `loglik(fit)` with no further arguments equals `logLik(fit)` to the last
#' digit, which is the cheapest check that the callable route and the fitting
#' route are the same model.
#'
#' # Generics, not closures on the fit
#'
#' A closure captures its environment, which here means the data. A fit
#' carrying one would hold the frame twice and would keep a stale copy after
#' the data changed.
#'
#' These rebuild from the specification the fit keeps whole, running each
#' term's blueprint against `data` by the same path [predict.StatmodFit()]
#' takes. A factor's levels, a spline's knots and a basis reparametrization
#' are therefore reapplied, never relearned.
#'
#' # No penalty enters
#'
#' All three are the likelihood alone. A model's penalized objective at given
#' coefficients is not what a caller profiling a likelihood wants, and the
#' hyperparameters are not arguments here.
#'
#' @param object A [StatmodFit()].
#' @param par A named list of coefficient vectors, one per distribution
#'   parameter, each as long as that equation's design is wide. Defaults to
#'   the fitted coefficients. Validated against the design, so a wrong length
#'   is an error, never a recycling.
#' @param data A data frame with the columns the model names. Defaults to the
#'   data the model was fitted to. New data must carry the response, since a
#'   likelihood needs one.
#' @param expected `hessian()` only: `TRUE` for the expected information,
#'   `FALSE` (the default) for the observed Hessian.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return `loglik()` gives a single number. `gradient()` gives a named list
#'   of numeric vectors, one per distribution parameter. `hessian()` gives a
#'   symmetric `p x p` matrix over the stacked coefficients, `p` being their
#'   total count.
#'
#' @seealso [logLik.StatmodFit()] for the maximized value with its degrees of
#'   freedom, [vcov.StatmodFit()] for the variance matrix,
#'   [predict.StatmodFit()] for the same reapplication on new data.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + dd$x + rnorm(40, sd = 0.5)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#'
#' # At the fitted coefficients it is the maximized value.
#' all.equal(loglik(fit), as.numeric(logLik(fit)))
#'
#' # Anywhere else it is lower, and the gradient says which way to go.
#' loglik(fit, par = list(mu = c(0, 0), sigma = 0))
#' gradient(fit, par = list(mu = c(0, 0), sigma = 0))
#'
#' # At the optimum the gradient vanishes.
#' max(abs(unlist(gradient(fit))))
#'
#' # And the observed information is positive definite there.
#' eigen(-hessian(fit), only.values = TRUE)$values
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
#' Returns the fit's own specification when `data` is `NULL`, and one carrying
#' the fitted terms reapplied to the new rows otherwise. This is the one place
#' the two cases are told apart, so every route that accepts a `data`
#' argument treats `NULL` alike.
#'
#' @param fit A [StatmodFit()].
#' @param data A data frame, or `NULL` for the data the model was fitted to.
#' @param need_response `TRUE` where the response must be present, as for a
#'   log-likelihood; `FALSE` for a prediction, new data routinely having no
#'   response column.
#'
#' @return A [StatmodSpec()]: `fit@spec` itself when `data` is `NULL`, and
#'   [statmod_respec()]'s result otherwise.
#'
#' @keywords internal
spec_at <- function(fit, data, need_response = TRUE) {
  if (is.null(data)) return(fit@spec)
  statmod_respec(fit@spec, data, need_response = need_response)
}

#' Resolve a Parameter Structure
#'
#' @description
#' Returns the fitted coefficients when `par` is `NULL`, and otherwise checks
#' a supplied structure against the design before returning it.
#'
#' @details
#' Three things are checked: that `par` is a named list, that it names every
#' distribution parameter, and that each vector is as long as that equation's
#' design is wide. Each failure signals an error naming the parameter, which
#' is the alternative to R recycling a short vector into a different model
#' without a word.
#'
#' @param fit A [StatmodFit()], read for its fitted coefficients.
#' @param par A named list of coefficient vectors, or `NULL`.
#' @param design The design to validate the lengths against.
#'
#' @return A named list of coefficient vectors, one per distribution
#'   parameter in the family's order.
#'
#' @seealso [loglik()], whose three methods all pass through this.
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
#' number of observations. [loglik()] is the other thing: the model
#' evaluated at parameters and data of the caller's choosing.
#' @details
#' # Which likelihood, and why it matters
#'
#' The default is the **conditional** one: the log-density at the fitted
#' coefficients, a penalized coefficient among them, paired with the
#' effective degrees of freedom \eqn{\mathrm{tr}[(H+S)^{-1}H]}. A criterion
#' built on that pair is the conditional AIC, and it asks how well the model
#' describes the groups, curves and states actually observed.
#'
#' A mixed-model package reports the **marginal** likelihood instead: the
#' random effects are integrated out and the count is the number of estimated
#' parameters, variance components among them. That asks about the
#' population, and its AIC is a different number, not comparable with the
#' conditional one.
#'
#' Neither convention allows the halves to be mixed: a marginal likelihood
#' against an effective count, or a conditional one against a parameter count
#' (Vaida and Blanchard, 2005).
#'
#' Measured on a random intercept over 120 groups, the two readings of the
#' same fit are -4148.59 on 115.72 effective degrees of freedom and -4372.79
#' on 4, against `lme4::lmer`'s marginal -4371.71 on 4.
#'
#' # When the marginal one is available
#'
#' `type = "marginal"` returns the value the outer criterion evaluated while
#' choosing the hyperparameters, with the number of estimated parameters as
#' its degrees of freedom. It exists only where a marginal criterion actually
#' ran, which is [ml()] or [reml()].
#'
#' Where the hyperparameters were held, or chosen by [aic()], [bic()] or
#' [cv()], there is no marginal likelihood to report and asking for one
#' signals an error instead of returning a number that would look like one.
#' @param object A [StatmodFit()].
#' @param type `"conditional"` (the default) or `"marginal"`.
#' @param ... Unused.
#' @return A `logLik` object: a single number with attributes `df`, the
#'   degrees of freedom, and `nobs`. [stats::AIC()] and [stats::BIC()] read
#'   it, so `AIC(fit)` is the conditional AIC.
#' @references
#' Vaida, F. and Blanchard, S. (2005). Conditional Akaike information for
#' mixed-effects models. *Biometrika*, 92(2), 351--370.
#' @seealso [loglik()]
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

#' @title The Coefficients of a Fitted Model
#'
#' @description
#' One named vector per distribution parameter: the quantities the model is
#' written in, or the coordinates it was estimated on.
#'
#' @details
#' # What `readable` moves, and what it does not
#'
#' Most coefficients are the same under either reading. A coefficient of a
#' linear predictor is what it is. Two kinds of parameter are reported under
#' a different name from the one they are carried under, and those are the
#' ones this argument moves.
#'
#' A **break-point** term is fitted through a working pair and its position
#' is read off that pair. A discontinuous term carries \eqn{g} and reports
#' \eqn{\psi = -g/\delta}, so at `readable = FALSE` the vector holds a number
#' that is no quantity of the model at all.
#'
#' A **score-driven** term's persistence rides a partial autocorrelation, the
#' stationary region not being a box, and what the literature calls
#' \eqn{\beta_j} is the autoregressive coefficient the whole chart produces.
#' At \eqn{q = 2} a fit reporting \eqn{\beta_1 = 0.761} has a free coordinate
#' of \eqn{\mathrm{pacf}_1 = 0.857}.
#'
#' Where a term declares no quantities of its own the coordinates stand. So
#' do the coefficients of a parameter developed over covariates: a
#' development is a vector with no single value to report, so the term
#' declares nothing for it.
#'
#' # A structural term is present under both readings
#'
#' It contributes no design columns, so its parameters are in no block. They
#' used to be in neither reading, and a model whose whole predictor is a
#' score-driven filter answered `numeric(0)`.
#'
#' # When to ask for `FALSE`
#'
#' `readable = FALSE` gives the vector the fit was estimated on, in the order
#' and under the names [vcov.StatmodFit()] is indexed by, with a structural
#' term's part on the unconstrained scale its charts define. That is what a
#' caller feeding a fit back to [loglik()] or to an optimizer needs.
#'
#' Hyperparameters are not coefficients and are not here. [hyper()] reports
#' them.
#'
#' @param object A [StatmodFit()].
#' @param readable `TRUE`, the default, reports the quantities the model is
#'   written in. `FALSE` reports the coordinates it was estimated on.
#' @param ... Unused.
#'
#' @return A named list with one entry per distribution parameter, in the
#'   family's order, each a named numeric vector. The names are the
#'   coefficient labels, composed from each term's own label, so two terms of
#'   one kind in one formula stay apart.
#'
#' @seealso [hyper()] for the hyperparameters, [vcov.StatmodFit()] for the
#'   variance matrix indexed by the `readable = FALSE` names,
#'   [summary.StatmodFit()] for the two together,
#'   [modelterms7::term_readable()] for what a term declares.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(200, 0, 10)))
#' d$y <- 0.3 * d$x + 1.5 * pmax(d$x - 6, 0) + rnorm(200, 0, 0.4)
#' fit <- statmod(y ~ modelterms7::seg(x, psi = 4),
#'                distributions7::gaussian1_distrib(), d)
#'
#' # The break-point is a quantity of the model, near the true 6.
#' coef(fit)$mu
#'
#' # The coordinates the fit ran on are the same numbers here, seg() being
#' # continuous, and vcov() is indexed by these names.
#' coef(fit, readable = FALSE)$mu
#' rownames(vcov(fit))
#'
#' @keywords internal
coef.StatmodFit <- function(object, readable = TRUE, ...) {
  spec <- object@spec
  params <- spec@distrib@params
  design <- statmod_design(spec)
  stats::setNames(lapply(params, function(p) {
    v <- stats::setNames(object@coefficients[[p]], design[[p]]$coef_names)
    if (isTRUE(readable)) v <- coef_readable(spec, design, object, p, v)
    c(v, coef_structural(spec, object, p, readable))
  }), params)
}

#' The Quantities a Term Reports in Place of Its Coordinates
#'
#' @description
#' The coefficient vector of one equation with each term's declared
#' quantities put where the coordinates they are read from were.
#'
#' @details
#' A term says what it is about through [modelterms7::term_readable()], which
#' returns the quantities and the Jacobian from the coefficients they are
#' read from. The columns that Jacobian touches are exactly the coordinates
#' to replace; a coordinate no quantity reads stays where it is.
#'
#' That rule keeps a developed parameter intact. A development is a vector of
#' coefficients over covariates with no single value to report, so the term
#' declares nothing for it and nothing is taken away.
#'
#' The names are composed as the term composes its coefficients', from its
#' own label, so two terms of one kind in one formula stay apart.
#'
#' @param spec The fitted specification.
#' @param design The design.
#' @param fit The [StatmodFit()].
#' @param p The distribution parameter naming the equation, a string.
#' @param v The named coefficient vector of that equation.
#'
#' @return `v`, with each term's declared quantities in place of the
#'   coordinates they are read from. `v` unchanged where no term of the
#'   equation declares any.
#'
#' @seealso [coef.StatmodFit()], the caller,
#'   [modelterms7::term_readable()] for what a term declares.
#'
#' @keywords internal
coef_readable <- function(spec, design, fit, p, v) {
  keep <- rep(TRUE, length(v))
  ins <- list()
  for (nm in names(spec@terms[[p]])) {
    idx <- design[[p]]$blocks[[nm]]
    if (!length(idx)) next
    term <- spec@terms[[p]][[nm]]
    rd <- tryCatch(modelterms7::term_readable(term,
                                              fit@coefficients[[p]][idx]),
                   error = function(e) NULL)
    if (is.null(rd) || !length(rd$name)) next
    sup <- which(apply(rd$jacobian != 0, 2L, any))
    if (!length(sup)) next
    keep[idx[sup]] <- FALSE
    ins[[as.character(idx[sup][1L])]] <- stats::setNames(
      as.numeric(rd$value), paste(term@label, rd$name, sep = "."))
  }
  if (!length(ins)) return(v)
  out <- list()
  for (i in seq_along(v)) {
    k <- as.character(i)
    if (!is.null(ins[[k]])) out[[length(out) + 1L]] <- ins[[k]]
    if (keep[[i]]) out[[length(out) + 1L]] <- v[i]
  }
  unlist(out)
}

#' What a Structural Term Contributes to the Coefficients
#'
#' @description
#' The own parameters of a structural term of one equation, as quantities or
#' as the coordinates they were estimated on.
#'
#' @details
#' A structural term rewrites the likelihood instead of adding columns to a
#' design, so its parameters sit in no block and appeared in no reading of
#' [coef.StatmodFit()]: a model whose whole predictor is a score-driven
#' filter answered with an empty vector. They are named from the term's
#' label, as every other coefficient of that term is.
#'
#' A **held** parameter is one an intercept in the same equation carries, so
#' it is not estimated. It is reported under either reading, at the value it
#' is held at. Leaving it out would make the vector shorter than the term's
#' own parameter count and break the correspondence with [vcov.StatmodFit()].
#'
#' @param spec The fitted specification.
#' @param fit The [StatmodFit()].
#' @param p The distribution parameter naming the equation, a string.
#' @param readable `TRUE` for the quantities the term declares, `FALSE` for
#'   the coordinates on the unconstrained scale.
#'
#' @return A named numeric vector of that term's own parameters.
#'   `numeric(0)` where the equation carries no structural term.
#'
#' @seealso [coef.StatmodFit()], the caller,
#'   [statmod_latent()] for a latent-state term's smoothed states.
#'
#' @keywords internal
coef_structural <- function(spec, fit, p, readable) {
  out <- numeric(0)
  sp <- fit@structural
  if (!length(sp)) return(out)
  for (nm in names(sp)) {
    if (!identical(sp[[nm]]$equation, p)) next
    term <- spec@terms[[p]][[nm]]
    z <- sp[[nm]]$unconstrained
    if (isTRUE(readable)) {
      rd <- tryCatch(modelterms7::term_readable(term, z),
                     error = function(e) NULL)
      w <- if (is.null(rd) || !length(rd$name)) unlist(sp[[nm]]$parameter) else
        stats::setNames(as.numeric(rd$value), rd$name)
    } else {
      w <- unlist(z)
    }
    lb <- tryCatch(term@label, error = function(e) "")
    if (length(lb) == 1L && nzchar(lb)) names(w) <- paste(lb, names(w),
                                                          sep = ".")
    out <- c(out, w)
  }
  out
}
S7::method(coef, StatmodFit) <- coef.StatmodFit

#' @title The Residuals of a Fitted Model
#' @name residuals.StatmodFit
#' @description
#' One residual per observation, comparing it with the whole distribution the
#' model puts on it rather than with any one of that distribution's
#' parameters.
#'
#' @details
#' # One residual per observation
#'
#' A residual asks whether an observation is consistent with the law its row
#' was given, and that law is one object carrying every parameter at once. So
#' there is one residual per observation, however many parameters the model
#' develops over covariates.
#'
#' A per-parameter quantity exists and is a different thing: the contribution
#' \eqn{\partial \ell_i / \partial \eta_{ip}} says which equation an
#' observation strains, and the partial residuals of one equation are what a
#' term's effect is drawn against.
#'
#' # The quantile residual, and why it is the default
#'
#' \deqn{r_i = \Phi^{-1}(F(y_i; \hat\theta_i))}
#'
#' Under a correct model \eqn{F(y_i; \theta_i)} is exactly uniform, so
#' \eqn{r_i} is exactly standard normal: whatever the family, and whichever
#' of its parameters are modeled. It privileges no parameter, and its
#' reference distribution needs no asymptotics (Dunn and Smyth, 1996).
#'
#' # Where the distribution function jumps
#'
#' The construction is randomized: \eqn{u_i} is drawn uniformly on
#' \eqn{(F(y_i^-), F(y_i))} and the residual is \eqn{\Phi^{-1}(u_i)}. Exact
#' again, at the price of being random, so two calls give two answers.
#'
#' This applies to every discrete family, and at the atom alone to a mixed
#' one, which is the zero-adjusted wrapper of a continuous parent. `seed`
#' makes a call reproducible without disturbing the caller's stream; left
#' `NULL`, the ambient state is used and nothing is set.
#'
#' # Pearson and response residuals
#'
#' \eqn{(y_i - \mathbb{E}[Y_i]) / \mathrm{sd}(Y_i)} and its numerator alone.
#' Both are defined against the mean, and for a skewed family the Pearson
#' residual is not standard normal even where the model is right, so its
#' quantile-quantile plot misleads in exactly the case a distributional model
#' is for. They are here because they are familiar.
#'
#' @param object A [StatmodFit()].
#' @param type `"quantile"` (the default), `"pearson"` or `"response"`.
#'   Matched with [match.arg()].
#' @param seed An integer to seed the randomization with, or `NULL` for the
#'   ambient state. Read only where the distribution function jumps, so it
#'   changes nothing for a continuous family.
#' @param ... Unused.
#'
#' @return A numeric vector with one entry per observation. Standard normal
#'   under a correct model for `"quantile"`; approximately standardized for
#'   `"pearson"`; on the response's own scale for `"response"`.
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
#' *Journal of Computational and Graphical Statistics* 5(3), 236--244.
#'
#' @seealso [fitted.StatmodFit()] and [predict.StatmodFit()] for the fitted
#'   parameters these are read against.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(200, -2, 2))
#' d$y <- 1 + 0.8 * d$x + rnorm(200, 0, 0.5)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)
#'
#' # Under a correct model the quantile residuals are standard normal.
#' r <- residuals(fit)
#' c(mean = mean(r), sd = stats::sd(r))
#' stats::shapiro.test(r)$p.value
#'
#' # For a discrete family the same residual is randomized, so two calls
#' # differ unless a seed is given.
#' dp <- data.frame(x = runif(200, -1, 1))
#' dp$y <- rpois(200, exp(1 + dp$x))
#' fp <- statmod(y ~ x, distributions7::poisson_distrib(), dp)
#' identical(residuals(fp, seed = 1), residuals(fp, seed = 1))
#' identical(residuals(fp, seed = 1), residuals(fp, seed = 2))
#'
#' @keywords internal
residuals.StatmodFit <- function(object,
                                 type = c("quantile", "pearson", "response"),
                                 seed = NULL, ...) {
  type <- match.arg(type)
  spec <- object@spec
  d <- spec@distrib
  y <- spec@response
  th <- object@fitted
  if (!length(th)) {
    stop("The fit carries no fitted parameters to compare the response with.",
         call. = FALSE)
  }
  if (identical(type, "response") || identical(type, "pearson")) {
    m <- mean(d, th)
    r <- as.numeric(y) - as.numeric(m)
    if (identical(type, "response")) return(r)
    s <- as.numeric(distributions7::std_dev(d, th))
    return(r / s)
  }
  fy <- as.numeric(distributions7::distrib_cdf(d, y, th))
  # WHERE THE DISTRIBUTION FUNCTION JUMPS, and by how much. A discrete family
  # jumps at every observation and its mass IS what distrib_pdf returns
  # there; a mixed one jumps at its declared atoms alone, where the same call
  # returns the atom's probability. Asking the density rather than
  # differencing the distribution function at the previous support point is
  # what keeps this right for a family whose support is not the integers.
  jump <- rep(FALSE, length(fy))
  if (S7::S7_inherits(d, distributions7::discrete_distrib)) {
    jump[] <- TRUE
  } else {
    at <- tryCatch(distributions7::distrib_atoms(d, th),
                   error = function(e) list(y = numeric(0)))
    if (length(at$y)) jump <- as.numeric(y) %in% at$y
  }
  u <- fy
  if (any(jump)) {
    py <- as.numeric(distributions7::distrib_pdf(d, y, th))
    lo <- pmax(fy - py, 0)
    if (!is.null(seed)) {
      # the caller's stream is left where it was: a residual is not a reason
      # to move it
      old <- if (exists(".Random.seed", envir = globalenv())) {
        get(".Random.seed", envir = globalenv())
      } else NULL
      set.seed(seed)
      on.exit({
        if (is.null(old)) {
          suppressWarnings(rm(".Random.seed", envir = globalenv()))
        } else {
          assign(".Random.seed", old, envir = globalenv())
        }
      }, add = TRUE)
    }
    u[jump] <- stats::runif(sum(jump), lo[jump], fy[jump])
  }
  stats::qnorm(pmin(pmax(u, 0), 1))
}
S7::method(residuals, StatmodFit) <- residuals.StatmodFit


#' @title The Hyperparameters of a Fitted Model
#'
#' @description
#' Reports every hyperparameter of every penalty the model carries: one row
#' each, with the value, whether it was held, and what put it there. A
#' smoothing parameter, a prior scale, a lasso's \eqn{\lambda} and an elastic
#' net's \eqn{\alpha} all appear here.
#'
#' @details
#' # A hyperparameter is not a coefficient
#'
#' It governs the coefficients under it instead of sitting beside them, and
#' the two are estimated by different routes and reported with different
#' qualifications. [coef.StatmodFit()] holds the coefficients and this holds
#' the hyperparameters; neither holds the other.
#'
#' # The two scales
#'
#' `scale = "parameter"`, the default, is the scale the penalty is written
#' on: a smoothing parameter is a positive number, and a gaussian prior's
#' `sigma` is a scale. That is what a reader wants.
#'
#' `scale = "link"` is the free scale the outer search runs on, through each
#' hyperparameter's own link. That is what a caller comparing two fits'
#' searches wants. Where a hyperparameter carries no link the two coincide.
#'
#' # What `source` says that `held` cannot
#'
#' `held` is a logical and says only whether the value moved. `source` says
#' what put it there:
#'
#' - `"fixed"` for one the term itself held, as `s(x, lambda = 2)`.
#' - the criterion's name, `"reml"` or `"ml"`, for one a marginal criterion
#'   maximized.
#' - the criterion that scored the path, `"bic"` and so on, for one chosen
#'   along a path over its own values.
#'
#' The distinction has a consequence a reader needs. A value chosen along a
#' path is the argument of a minimum over a grid, not the root of a
#' derivative, so no standard error follows from it. One a marginal criterion
#' reached carries one, from the curvature of that criterion, and
#' [summary.StatmodFit()] prints it.
#'
#' @param fit A [StatmodFit()].
#' @param scale `"parameter"` (the default) or `"link"`. Matched with
#'   [match.arg()].
#'
#' @return A data frame with one row per hyperparameter and six columns:
#'   \describe{
#'     \item{`parameter`}{the distribution parameter whose equation the
#'       penalized term sits in.}
#'     \item{`term`}{the term's key, its call as written.}
#'     \item{`name`}{the hyperparameter's own name, as the penalty names it.}
#'     \item{`estimate`}{its value, on the scale asked for.}
#'     \item{`held`}{a logical: whether the term fixed it.}
#'     \item{`source`}{what put the value there, as above.}
#'   }
#'   A data frame of no rows, with those columns, where the model carries no
#'   penalty at all.
#'
#' @seealso [coef.StatmodFit()] for the coefficients,
#'   [summary.StatmodFit()] for the two printed together with standard
#'   errors, [statmod_held()] for which are held.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(80, 0, 1))
#' d$y <- sin(3 * d$x) + rnorm(80, 0, 0.3)
#'
#' # Estimated by REML, which is what source says.
#' fit <- statmod(y ~ s(x, k = 6), distributions7::gaussian1_distrib(), d)
#' hyper(fit)
#'
#' # The same value on the scale the outer search ran on.
#' hyper(fit, scale = "link")
#'
#' # Held by the term instead, and reported as fixed.
#' held <- statmod(y ~ s(x, k = 6, lambda = 2),
#'                 distributions7::gaussian1_distrib(), d)
#' hyper(held)[, c("name", "estimate", "held", "source")]
#'
#' @export
hyper <- function(fit, scale = c("parameter", "link")) {
  scale <- match.arg(scale)
  spec <- fit@spec
  held <- tryCatch(statmod_held(spec, statmod_design(spec)),
                   error = function(e) character(0))
  outer_kind <- if (is.null(fit@methods$outer)) NA_character_ else
    fit@methods$outer@kind
  spc <- fit@methods$sparse_criterion
  spc_kind <- if (is.null(spc)) NA_character_ else spc@kind
  rows <- list()
  for (u in statmod_penalty_keys(spec)) {
    th <- fit@hyper[[u$param]][[u$key]]
    if (is.null(th) || !length(th)) next
    kink <- isTRUE(tryCatch(penalty_has_kink(u$penalty),
                            error = function(e) FALSE))
    for (h in names(th)) {
      v <- as.numeric(th[[h]])
      lk <- u$penalty@link_params[[h]]
      if (identical(scale, "link") && !is.null(lk)) {
        v <- linkfunctions7::linkfun(lk, v)
      }
      # HELD is the term's answer and nobody else's, and it is asked of the
      # same enumeration the outer index asks: a hyperparameter a term fixed
      # is fixed whatever criterion ran beside it.
      hk <- paste(u$param, u$key, h, sep = "\r")
      is_held <- hk %in% held || h %in% names(u$fixed)
      src <- if (is_held) "fixed" else if (kink) spc_kind else outer_kind
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = u$param, term = u$key, name = h, estimate = v,
        held = is_held, source = if (is.na(src)) "" else src,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) {
    return(data.frame(parameter = character(0), term = character(0),
                      name = character(0), estimate = numeric(0),
                      held = logical(0), source = character(0),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


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
#' @param x A [StatmodFit()].
#'
#' @return A single string, empty when the parameters cannot be read. It is
#'   a note of [summary.StatmodFit()] rather than a line of
#'   [print.StatmodFit()]: it qualifies the fit rather than
#'   describing it, and it is read once when something looks wrong.
#'
#' @seealso [statmod()]
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
  paste0("The fit did not converge. The parameters it reached: ",
         paste(trimws(gsub("[ ]{2,}", " ", rows)), collapse = "; "),
         ". A scale at 1e-15 says the rest on its own.")
}


#' A Term Key Shortened for Display
#'
#' @description
#' The key with the arguments of its leading call reduced to the first one,
#' so that a trace line names the term instead of repeating its whole
#' specification.
#'
#' @details
#' A key is the call that produced the term, and for a penalty over a
#' sub-term it is that call followed by `::parameter::sub-term`. The
#' call is what grows: printed in full, one line of an outer trace carried
#' the deparsed `gas(p = 1, q = 1, time = t, by = ~ridge(~id), links =
#' list(...))` three times over, which is a line no reader can use. Only the
#' leading call is shortened, and only past its first argument, so
#' `s(x, k = 20)` and `s(z, k = 8)` stay apart; everything after
#' `::` is kept whole, that being what distinguishes one entry of a term
#' from another.
#'
#' Where shortening would make two labels the same the FULL ones are
#' returned, all of them: a shorter label that is ambiguous is worse than a
#' long one, and deciding per label would leave a reader unable to tell which
#' convention a given line follows.
#'
#' @param x A character vector of keys.
#'
#' @return A character vector the same length.
#'
#' @keywords internal
short_keys <- function(x) {
  one <- function(k) {
    parts <- strsplit(k, "::", fixed = TRUE)[[1L]]
    call <- parts[1L]
    op <- regexpr("(", call, fixed = TRUE)
    if (op < 1L || !endsWith(call, ")")) return(k)
    head_ <- substr(call, 1L, op)
    args <- substr(call, op + 1L, nchar(call) - 1L)
    if (!nzchar(args)) return(k)
    # the first argument, found by scanning for a comma at nesting depth zero
    depth <- 0L
    cut <- nchar(args) + 1L
    ch <- strsplit(args, "", fixed = TRUE)[[1L]]
    for (i in seq_along(ch)) {
      c_ <- ch[i]
      if (c_ %in% c("(", "[", "{")) depth <- depth + 1L
      else if (c_ %in% c(")", "]", "}")) depth <- depth - 1L
      else if (c_ == "," && depth == 0L) { cut <- i; break }
    }
    if (cut > nchar(args)) return(k)
    parts[1L] <- paste0(head_, substr(args, 1L, cut - 1L), ", ...)")
    paste(parts, collapse = "::")
  }
  out <- vapply(x, one, character(1), USE.NAMES = FALSE)
  if (anyDuplicated(out) && !anyDuplicated(x)) x else out
}


#' A Titled Rule for a Verbose Trace
#'
#' @description
#' One blank line, then a titled rule with the method that will do the work
#' named on the right, so that a reader of a running fit can see where one
#' step ends and the next begins.
#'
#' @details
#' The trace has three nested things to say -- which outer step, which pass
#' of the alternation inside it, and what each block did -- and printed as
#' undifferentiated lines they are unreadable, as a panel fit with
#' three hyperparameters and 130 outer evaluations demonstrated. Naming the
#' method on every rule answers the question a reader of a slow fit actually
#' has, and that is the one running now.
#'
#' @param title The step, e.g. `"outer 3"`.
#' @param method What runs it, or `NULL`.
#' @param indent How far in, in spaces.
#' @param char The rule's character.
#'
#' @return Invisibly `NULL`; prints.
#'
#' @keywords internal
vb_rule <- function(title, method = NULL, indent = 0L, char = "-") {
  pad <- strrep(" ", indent)
  left <- paste0(pad, strrep(char, 2L), " ", title, " ")
  right <- if (is.null(method) || !nzchar(method)) "" else
    paste0(" [", method, "]")
  fill <- max(0L, 72L - nchar(left) - nchar(right))
  cat("\n", left, strrep(char, fill), right, "\n", sep = "")
  invisible(NULL)
}

#' A Detail Line of a Verbose Trace
#'
#' @param ... Passed to `sprintf`.
#' @param indent How far in, in spaces.
#'
#' @return Invisibly `NULL`; prints.
#'
#' @keywords internal
vb_say <- function(..., indent = 5L) {
  cat(strrep(" ", indent), sprintf(...), "\n", sep = "")
  invisible(NULL)
}

#' The Name of Whatever Is About to Run
#'
#' @description
#' An optimizer's own name, a criterion's kind, or a short description of a
#' step that is not an optimizer.
#'
#' @param x An optimizer, a criterion, or `NULL`.
#' @param default What to say when there is no object to ask.
#'
#' @return A single string.
#'
#' @keywords internal
vb_name <- function(x, default = "") {
  if (is.null(x)) return(default)
  if (S7::S7_inherits(x) && "name" %in% names(S7::props(x))) return(x@name)
  if (S7::S7_inherits(x, Iwls)) return("iwls")
  if (S7::S7_inherits(x) && "kind" %in% names(S7::props(x))) {
    return(toupper(x@kind))
  }
  default
}
