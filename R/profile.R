#' @include inference_objects.R restricted.R marginal.R path.R
NULL

# The profile likelihood of one or several coefficients.
#
# A profile is a table: for each combination of values of the coefficients
# that are held, the maximum of the penalized likelihood over every other
# coefficient. This file builds the table from the restricted refit of
# restricted.R, with the two differences that the table needs. The
# coefficients that the outer criterion estimates (marginal_coords()) are
# estimated by the criterion at every point and are not carried to the joint
# mode, and the hyperparameters are either kept at the values of the fit or
# searched again at every point.


#' @title Values of a Profile
#' @name profile_range
#'
#' @description
#' Selects the values at which [statmod_profile()] evaluates the profile of a
#' coefficient: `n` equally spaced values from `width` standard errors below
#' the estimate to `width` standard errors above it, on the scale of the
#' coefficient in the linear predictor.
#'
#' @details
#' The standard error is the one of `vcov(fit)`, so it follows the variance
#' convention that `vcov()` uses by default. The grid is symmetric about the
#' estimate, and it contains the estimate itself when `n` is odd.
#'
#' The scale is the unconstrained one. A coefficient of a linear predictor
#' lives there, and so do the parameters of a nonlinear term, the positions
#' of the break-points and the parameters of a covariance structure; the
#' scale of the coefficients is described at [statmod_profile()].
#'
#' @param width A single positive number: the half-width of the grid in
#'   standard errors.
#' @param n A single integer, at least 2: the number of values.
#'
#' @return An object of class `ProfileRange`, with the properties `width` and
#'   `n`.
#'
#' @seealso [statmod_profile()], which takes it in its argument `values`.
#'
#' @examples
#' profile_range()
#' profile_range(width = 2, n = 9)
#' @export
profile_range <- function(width = 3, n = 21) {
  ProfileRange(width = width, n = as.integer(n))
}

#' @rdname profile_range
#' @aliases ProfileRange
ProfileRange <- S7::new_class(
  "ProfileRange",
  properties = list(width = S7::class_numeric, n = S7::class_integer),
  validator = function(self) {
    if (length(self@width) != 1L || !is.finite(self@width) || self@width <= 0) {
      return("'width' must be a single positive number.")
    }
    if (length(self@n) != 1L || is.na(self@n) || self@n < 2L) {
      return("'n' must be a single integer of at least 2.")
    }
    NULL
  }
)


#' @title Methods of Profiling
#' @name profile_method
#'
#' @description
#' `profile_method` is the abstract parent of the objects that select how
#' [statmod_profile()] obtains the profile at each value of the grid. It
#' serves to test membership with `S7::S7_inherits()`. [profile_grid()] is
#' the method of this version.
#'
#' @param name A single string, the short name of the method.
#'
#' @return Nothing: the class is abstract and its constructor signals an
#'   error.
#'
#' @seealso [profile_grid()], [statmod_profile()].
#'
#' @examples
#' S7::S7_inherits(profile_grid(), profile_method)
#' try(profile_method(name = "mine"))
#' @export
profile_method <- S7::new_class(
  "profile_method",
  properties = list(name = S7::class_character),
  abstract = TRUE
)

#' @rdname profile_grid
#' @aliases ProfileGrid
ProfileGrid <- S7::new_class(
  "ProfileGrid",
  parent = profile_method,
  properties = list(start = S7::class_character),
  validator = function(self) {
    if (length(self@start) != 1L ||
        !self@start %in% c("neighbor", "fit", "both")) {
      return("'start' must be one of \"neighbor\", \"fit\" and \"both\".")
    }
    NULL
  }
)


#' @title Profile by Refitting on a Grid
#' @name profile_grid
#'
#' @description
#' Selects the exact method of [statmod_profile()]: at every value of the grid
#' the model is fitted again with the coefficients held at that value, with
#' the tolerances, the inner optimizer and the outer criterion of the fit.
#'
#' @details
#' The three starts differ in the coefficients from which each refit begins.
#' `"neighbor"` starts the refit at a value from the refit of the adjacent
#' value that lies nearer to the estimate, going outward from the estimate in
#' both directions, so that the coefficients move little between two refits.
#' The first value on each side starts from the coefficients of the fit.
#' With several held coefficients the chains run along the first one, one
#' pair of chains (below and above the estimate) for each combination of the
#' values of the others. `"fit"` starts every refit at the coefficients of the
#' fit, with the held ones set to the value of the point; the points are then
#' independent of each other. `"both"` runs the two and keeps, at every point,
#' the refit with the larger penalized log-likelihood, at twice the cost. A
#' profile with several local maxima, as that of a break-point has, can be
#' reached from one start and missed from the other.
#'
#' With `cores` greater than one in [statmod_profile()], the chains of
#' `"neighbor"` and the points of `"fit"` are distributed over worker
#' processes.
#'
#' @param start One of `"neighbor"` (the default), `"fit"` and `"both"`.
#'
#' @return A `ProfileGrid` object, inheriting from [profile_method], with the
#'   properties `name` (`"grid"`) and `start`.
#'
#' @seealso [statmod_profile()], [profile_method].
#'
#' @examples
#' profile_grid()
#' profile_grid(start = "both")
#' @export
profile_grid <- function(start = c("neighbor", "fit", "both")) {
  ProfileGrid(name = "grid", start = match.arg(start))
}


#' Compute a Profile
#'
#' @description
#' The generic through which [statmod_profile()] reaches a method: it
#' receives the context built for the call and the grid of values, and returns
#' the list of points that [profile_assemble()] turns into a table.
#'
#' @param method A [profile_method].
#' @param ctx The context of the call, as `profile_context()` builds it.
#' @param grid The data frame of held values, one column per coefficient.
#' @param cores The number of worker processes.
#'
#' @return A list with one element per row of `grid`, each a list of
#'   `par`, `hyper`, `objective`, `loglik`, `criterion` and `converged`.
#'
#' @keywords internal
profile_compute <- S7::new_generic("profile_compute", "method")


#' @title The Profile Likelihood of Coefficients
#' @name statmod_profile
#'
#' @description
#' Returns, for each combination of values of one or several coefficients, the
#' maximum of the penalized log-likelihood over all the other coefficients
#' when those are held at the values. The result is a data frame, one row per
#' combination, which serves a plot, an interval or a diagnostic; the
#' function computes none of those.
#'
#' @details
#' # What the columns report
#'
#' The columns `profile_loglik`, `profile_penalty` and `profile_pen_loglik`
#' belong to the INNER criterion, the penalized log-likelihood of the
#' coefficients at the hyperparameters in force. At the point where the
#' coefficients are held at the values of the row, `profile_loglik` is the
#' log-likelihood without the penalty, `profile_penalty` is the penalty with
#' the sign that makes the two add up, that is the logarithm of the prior
#' density of the coefficients, \eqn{-\tfrac12 P_\lambda(\beta)} with its
#' normalizing constant, and `profile_pen_loglik` is their sum, which is the
#' quantity that the refit maximizes. The constant depends on the
#' hyperparameters and not on the coefficients, so the penalty is not
#' necessarily negative, and with the hyperparameters held its changes along
#' the grid are those of the quadratic form. A model without penalty has
#' `profile_penalty` equal to zero, and the two log-likelihoods then
#' coincide.
#'
#' The column `outer_criterion` belongs to the OUTER criterion ([reml()],
#' [ml()], [aic()] or [bic()]), a different function of a different object:
#' the value, on the sign and the scale that the fit reports in
#' `fit@criterion`, of the criterion that chooses the hyperparameters or the
#' coefficients that the criterion estimates. It is not the profile of that
#' criterion over the held coefficients in any sense that the first three
#' columns have for the inner one, and it is `NA` when the point involved no
#' outer search: when the fit has no outer criterion, or when the
#' hyperparameters are held (`hyper = "held"`) and the criterion estimates no
#' coefficient.
#'
#' The column `converged` states that the point is at the maximum of its
#' refit. It is `TRUE` when the stopping rule of the refit fired, or when the
#' Newton decrement at the point is below the tolerance of
#' [statmod_certificate()], and in both cases the final inner fit converged.
#' The second condition is needed because the line search of the outer
#' optimizer can exhaust its backtracks at a point that is already the
#' optimum, without recording that its stopping rule fired. Such a point was
#' measured with the same criterion, to \eqn{10^{-6}}, as the neighboring
#' points whose flag was `TRUE`. The row of the attribute `maximum` reports
#' instead `fit@converged`, which is only the flag of the stopping rule. A
#' refit that signals an error gives a row of `NA` with `converged = FALSE`,
#' and the function warns.
#'
#' # The coefficients to hold: `which`
#'
#' A label has the form `"parameter:coefficient"`: the distribution
#' parameter whose equation carries the coefficient, a colon, and the name of
#' the coefficient as the second column of the table of [coef_labels()] gives
#' it, which is the name that `coef()` and `summary()` print. The labels of a
#' fit are the row names of `summary(fit)`, for example `"mu:(Intercept)"`,
#' `"mu:x"` or `"sigma:s(z).1"`. Labels must be unique.
#'
#' The coefficient of the effect of a covariate `z` on the parameter `theta`
#' of a [modelterms7::nl()] term, written `nl(~ a * exp(-theta * x),
#' theta ~ 1 + z)` in the equation of `mu`, has the label
#' `"mu:nl.theta.z"`. The name after the colon is `nl.`, the name of the
#' parameter of the term, a dot and the name of the coefficient of its
#' subformula, so that the intercept of `theta` is
#' `"mu:nl.theta.(Intercept)"`. A parameter of the term without a subformula
#' is a coefficient by itself, as in `"mu:nl.a"`.
#'
#' The position of a break-point of [modelterms7::seg()] is the coefficient
#' `psi1` (then `psi2`, and so on) of the term, as in `"mu:seg.psi1"`. When
#' the position has a subformula, as in `seg(w, psi ~ x)`, its coefficients
#' are `"mu:seg.psi1.(Intercept)"` and `"mu:seg.psi1.x"`. The names of all
#' the coefficients of a fit are the row names of
#' `statmodels7:::coef_labels(fit@spec, statmodels7:::statmod_design(fit@spec))`.
#'
#' When two or more positions of the same `seg()` term without subformulas
#' are in `which`, the positions of a model are ordered and distinct: the
#' combinations of the grid in which a position is not larger than the
#' position before it are not evaluated and are not in the table.
#'
#' # The scale
#'
#' The values are on the scale of the coefficient in the stacked vector of the
#' fit, which is unconstrained: a coefficient of a linear predictor is on
#' the scale of the link, and so are the coefficients of the predictor of a
#' parameter of a nonlinear term, the positions of the break-points and the
#' parameters of a covariance structure. Values on another scale must be
#' mapped to it by the caller. A monotone map relabels the grid and leaves the
#' columns of the profile unchanged.
#'
#' # The hyperparameters, and the coefficients of the criterion
#'
#' With `hyper = "held"` the hyperparameters (the smoothing parameters and the
#' scales of the random effects) stay at the values of the fit, and with
#' `hyper = "refit"` the outer search is run again at every point. In both
#' cases, when the outer criterion estimates coefficients ([reml()] and
#' [ml()] estimate the unpenalized coefficients of every parameter except the
#' first), those coefficients are estimated by the criterion at every point
#' and are not carried to the joint mode of the penalized likelihood, which
#' [statmod_restrict()] does. A held coefficient that the criterion would
#' have estimated is held, and the criterion estimates the rest. The cost of
#' a point is then that of an outer search over those coordinates, and it
#' exceeds the cost of a restricted fit.
#'
#' The outer optimizer is the default one of [statmod()], and the stopping
#' rules, the inner method and the information are those of the fit.
#'
#' # Coefficients that cannot be held yet
#'
#' A coefficient under a penalty with a kink (a lasso, a SCAD or an MCP), a
#' coefficient of a model with a structural term ([modelterms7::gas()],
#' [modelterms7::regime()], the marginal break-point term) and any
#' coefficient with `hyper = "refit"` in a model whose kinked penalties are
#' tuned by a path signal an error. An aliased coefficient has no estimate to
#' profile and signals an error.
#'
#' @param fit A [StatmodFit()].
#' @param which A character vector of unique labels `"parameter:coefficient"`,
#'   of one equation or of several.
#' @param values Either one [profile_range()], which applies to every
#'   coefficient of `which`, or a list. The named elements of a list carry
#'   names of `which` and are numeric vectors (the values, on the scale
#'   described above) or [profile_range()] objects. At most one element
#'   without a name sets the default for the coefficients that no element
#'   names, and without it that default is `profile_range()`. Two unnamed
#'   elements, or a name that is not in `which`, are errors. The values of
#'   each coefficient are sorted increasingly and repeated values are dropped.
#'   With several coefficients the table has one row for each combination of
#'   values, and the first coefficient varies fastest.
#' @param method A [profile_method]: [profile_grid()].
#' @param hyper `"held"` or `"refit"`: whether the hyperparameters stay at the
#'   values of the fit or are searched again at every point.
#' @param cores The number of worker processes, a positive integer. The
#'   default runs in this process. Under `pkgload` the run stays sequential.
#'
#' @return A data frame (`check.names = FALSE`) with one column for each
#'   coefficient of `which`, named by its label, then `profile_loglik`,
#'   `profile_penalty`, `profile_pen_loglik`, `outer_criterion` and
#'   `converged`. Its attributes are `maximum`, a one-row data frame with the
#'   same columns at the unrestricted fit (the numbers of the fit itself,
#'   which under [reml()] or [ml()] can lie slightly below the nearest rows
#'   of the table, because at a row the criterion does not integrate the held
#'   coefficients), `method`, the [profile_method]
#'   object, `hyper`, and `criterion`, the outer criterion object of the fit
#'   (`NULL` if it has none). The hyperparameters and the other coefficients
#'   of each point are not returned.
#'
#' @seealso [profile_grid()], [profile_range()], [statmod_test()] for the
#'   likelihood ratio at one value.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#'
#' # the profile of the slope, on five standard errors
#' prof <- statmod_profile(fit, which = "mu:x",
#'                         values = profile_range(width = 2, n = 9))
#' prof
#'
#' # the log-likelihood falls by z^2 / 2 at z standard errors
#' attr(prof, "maximum")
#'
#' # two coefficients, with values chosen by the caller for the second
#' statmod_profile(fit, which = c("mu:x", "mu:(Intercept)"),
#'                 values = list("mu:x" = profile_range(width = 1, n = 3),
#'                               "mu:(Intercept)" = c(0.8, 1, 1.2)))
#' @export
statmod_profile <- function(fit, which,
                            values = profile_range(width = 3, n = 21),
                            method = profile_grid(),
                            hyper = c("held", "refit"), cores = 1) {
  if (!S7::S7_inherits(fit, StatmodFit)) {
    stop("'fit' must be a StatmodFit.", call. = FALSE)
  }
  hyper <- match.arg(hyper)
  if (!S7::S7_inherits(method, profile_method)) {
    stop("'method' must be a profile method, as profile_grid() returns.",
         call. = FALSE)
  }
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1) {
    stop("'cores' must be a single positive integer.", call. = FALSE)
  }
  cores <- as.integer(cores)
  ctx <- profile_context(fit, which, hyper)
  vals <- profile_values(ctx, values)
  grid <- profile_prune(ctx, expand.grid(vals, KEEP.OUT.ATTRS = FALSE,
                                         stringsAsFactors = FALSE))
  pts <- profile_compute(method, ctx, grid, cores)
  profile_assemble(ctx, grid, pts, method)
}


#' The Context of a Profile
#'
#' @description
#' Resolves the labels of `which` against the fit, rejects the coefficients
#' that cannot be held, and collects what every refit of the profile reads:
#' the specification, the stacked estimates, the positions of the held
#' coefficients, the settings of the inner fit and the outer criterion.
#'
#' @param fit A [StatmodFit()].
#' @param which The labels.
#' @param hyper `"held"` or `"refit"`.
#'
#' @return A list.
#'
#' @keywords internal
profile_context <- function(fit, which, hyper) {
  if (!is.character(which) || !length(which) || anyNA(which)) {
    stop("'which' must be a character vector of labels \"parameter:coefficient\".",
         call. = FALSE)
  }
  if (anyDuplicated(which)) {
    stop("'which' must not repeat a label: ",
         paste(unique(which[duplicated(which)]), collapse = ", "), ".",
         call. = FALSE)
  }
  spec <- fit@spec
  params <- spec@distrib@params
  design <- statmod_design(spec)
  if (length(attr(design, "structural"))) {
    stop("The profile of a model that carries a structural term (gas(), ",
         "regime(), a marginal\n  break-point) is not supported yet: the ",
         "contribution of the term is the state of a recursion\n  and not a ",
         "column times a coefficient, and the hold does not reach it.",
         call. = FALSE)
  }
  lab <- coef_labels(spec, design)
  bad <- setdiff(which, rownames(lab))
  if (length(bad)) {
    stop(sprintf("Not a coefficient of this fit: %s.\n  The labels are %s%s.",
                 paste(sprintf("'%s'", bad), collapse = ", "),
                 paste(utils::head(rownames(lab), 12L), collapse = ", "),
                 if (nrow(lab) > 12L) ", ..." else ""), call. = FALSE)
  }
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  param <- lab[which, "parameter"]
  coefname <- lab[which, "coefficient"]
  at <- vapply(seq_along(which), function(i) {
    offs[[match(param[i], params)]] +
      match(coefname[i], design[[param[i]]]$coef_names)
  }, numeric(1))
  isal <- rownames(lab)[at] %in% fit@aliased
  if (any(isal)) {
    stop(sprintf(paste0("%s aliased, so it has no estimate to profile:\n  ",
                        "its column is a combination of the others and the ",
                        "fit left it out."),
                 paste(sprintf("'%s' is", which[isal]), collapse = ", ")),
         call. = FALSE)
  }
  kc <- kinked_coords(spec, design)
  if (any(at %in% kc)) {
    stop(sprintf(paste0("%s under a kinked penalty (a lasso, a SCAD or an ",
                        "MCP). The profile of such a\n  coefficient is not ",
                        "supported yet: that block is fitted by a coordinate ",
                        "descent, which\n  cannot leave one coordinate at a ",
                        "held value."),
                 paste(sprintf("'%s' is", which[at %in% kc]), collapse = ", ")),
         call. = FALSE)
  }
  method <- fit@methods$smooth
  cfg <- inner_settings(method)
  outer <- fit@methods$outer
  blocks <- statmod_blocks(spec, design)
  hyper0 <- fit@hyper
  if (identical(hyper, "refit") && !is.null(fit@methods$sparse_criterion) &&
      nrow(path_rows(spec, blocks, hyper0))) {
    stop("hyper = \"refit\" is not supported yet in a model whose kinked ",
         "penalties are tuned by a path:\n  the path would have to run at ",
         "every point. Use hyper = \"held\".", call. = FALSE)
  }
  beta_fit <- unlist(fit@coefficients[params], use.names = FALSE)
  # which search each point runs. The coordinates the criterion estimates are
  # read with the profiled coefficients held, so a profiled coefficient is not
  # among them
  hold0 <- stats::setNames(as.list(beta_fit[at]), which)
  spec_at <- profile_hold(spec, param, coefname, beta_fit[at])
  gam <- if (is.null(outer)) NULL else if (outer_minimize(outer)) {
    marginal_coords(spec_at, design, NULL)
  } else {
    marginal_coords(spec_at, design, outer)
  }
  ng <- length(gam$where)
  nh <- if (is.null(outer) || identical(hyper, "held")) 0L else
    nrow(outer_hyper_index(spec_at, blocks))
  psi <- profile_psi_groups(spec, design, param, coefname)
  list(fit = fit, spec = spec, params = params, psi = psi, which = which, param = param,
       coefname = coefname, at = at, estimate = beta_fit[at],
       beta_fit = beta_fit, hyper = hyper, hyper0 = hyper0, inner = method,
       cfg = cfg, outer = outer,
       emethod = if (is.null(outer)) cfg$emethod else
         select_emethod(outer, NULL, cfg$emethod),
       use_outer = !is.null(outer) && (nh + ng > 0L), ng = ng, nh = nh)
}


#' The Positions of One Break-Point Term Among the Held Coefficients
#'
#' @description
#' Finds, for each continuous [modelterms7::seg()] term without subformulas
#' that has two or more of its positions among the labelled coefficients, the
#' indices into `which` of those positions in the order of the term's
#' break-points.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design.
#' @param param,coefname The equation and the name of each labelled
#'   coefficient.
#'
#' @return A list of integer vectors, possibly empty.
#'
#' @keywords internal
profile_psi_groups <- function(spec, design, param, coefname) {
  st <- attr(design, "state")
  out <- list()
  if (is.null(st)) return(out)
  for (p in unique(param)) {
    for (nm in names(st$terms[[p]])) {
      tm <- st$terms[[p]][[nm]]
      if (!S7::S7_inherits(tm, modelterms7::SegTerm) ||
          !identical(tm@kind, "seg") || length(tm@subformulas) ||
          tm@npsi < 2L) {
        next
      }
      cols <- design[[p]]$blocks[[nm]]
      cn <- design[[p]]$coef_names[cols]
      bp <- tm@build_info
      psi <- cn[vapply(paste0("psi", seq_len(tm@npsi)),
                       function(q) bp$index[[q]][1L], integer(1))]
      i <- which(param == p & coefname %in% psi)
      if (length(i) >= 2L) out[[length(out) + 1L]] <- i[order(match(coefname[i], psi))]
    }
  }
  out
}


#' Drop the Combinations With Unordered Break-Points
#'
#' @description
#' Removes from the grid the rows in which the positions of one break-point
#' term are not strictly increasing in the order of its break-points. The
#' positions of a model are ordered and distinct, and a refit with two of them
#' held at equal or crossed values has no meaning.
#'
#' @param ctx The context of the call.
#' @param grid The data frame of held values.
#'
#' @return `grid` without those rows, with its row names renumbered.
#'
#' @keywords internal
profile_prune <- function(ctx, grid) {
  keep <- rep(TRUE, nrow(grid))
  for (g in ctx$psi) {
    v <- as.matrix(grid[g])
    keep <- keep & apply(v, 1L, function(r) all(diff(r) > 0))
  }
  if (!any(keep)) {
    stop("No combination of the values keeps the positions of a break-point ",
         "term in increasing order.", call. = FALSE)
  }
  out <- grid[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' A Specification With Named Coefficients Held
#'
#' @description
#' Writes the values of the labelled coefficients into `held_coef`, replacing
#' a hold of the same coefficient and keeping the other holds.
#'
#' @param spec A [StatmodSpec()].
#' @param param,coefname The equation and the name of each coefficient.
#' @param value Their values.
#'
#' @return `spec`, with `held_coef` set.
#'
#' @keywords internal
profile_hold <- function(spec, param, coefname, value) {
  hc <- spec@held_coef
  for (p in unique(param)) {
    j <- param == p
    v <- stats::setNames(as.numeric(value[j]), coefname[j])
    old <- hc[[p]]
    if (length(old)) old <- old[setdiff(names(old), names(v))]
    hc[[p]] <- c(old, v)
  }
  S7::set_props(spec, held_coef = hc)
}


#' The Values of a Profile for Each Coefficient
#'
#' @description
#' Turns the argument `values` of [statmod_profile()] into one sorted numeric
#' vector for each label.
#'
#' @param ctx The context of the call.
#' @param values The argument `values`.
#'
#' @return A named list of numeric vectors.
#'
#' @keywords internal
profile_values <- function(ctx, values) {
  which <- ctx$which
  default <- NULL
  given <- list()
  if (S7::S7_inherits(values, ProfileRange)) {
    default <- values
  } else if (is.list(values)) {
    nms <- names(values)
    if (is.null(nms)) nms <- rep("", length(values))
    if (sum(nms == "") > 1L) {
      stop("'values' has more than one element without a name; at most one ",
           "sets the default.", call. = FALSE)
    }
    extra <- setdiff(nms[nms != ""], which)
    if (length(extra)) {
      stop("'values' names ", paste(sprintf("'%s'", extra), collapse = ", "),
           ", which is not in 'which'.", call. = FALSE)
    }
    if (anyDuplicated(nms[nms != ""])) {
      stop("'values' names a coefficient twice.", call. = FALSE)
    }
    if (any(nms == "")) default <- values[[which(nms == "")]]
    given <- values[nms != ""]
  } else {
    stop("'values' must be a profile_range() or a list, as described at ",
         "?statmod_profile; wrap a numeric vector in list().", call. = FALSE)
  }
  if (is.null(default)) default <- profile_range()
  se <- NULL
  out <- list()
  for (i in seq_along(which)) {
    v <- if (!is.null(given[[which[i]]])) given[[which[i]]] else default
    if (S7::S7_inherits(v, ProfileRange)) {
      if (is.null(se)) {
        se <- sqrt(diag(stats::vcov(ctx$fit)))
      }
      s <- se[which[i]]
      if (is.na(s) || !is.finite(s) || s <= 0) {
        stop(sprintf(paste0("The standard error of '%s' is not available, ",
                            "so profile_range() cannot place the grid.\n  ",
                            "Give the values of that coefficient as a ",
                            "numeric vector."), which[i]), call. = FALSE)
      }
      v <- ctx$estimate[i] + s * v@width * seq(-1, 1, length.out = v@n)
    }
    if (!is.numeric(v) || !length(v) || any(!is.finite(v))) {
      stop(sprintf("The values of '%s' must be finite numbers.", which[i]),
           call. = FALSE)
    }
    out[[which[i]]] <- sort(unique(as.numeric(v)))
  }
  out
}


#' Run Independent Units of a Profile
#'
#' @description
#' The one place through which a profile reaches the worker processes of
#' [worker_map()].
#'
#' @param ctx The context of the call.
#' @param cores The number of worker processes.
#' @param n How many units.
#' @param body The unit's body, a function of the unit index.
#'
#' @return The list of results, in unit order.
#'
#' @keywords internal
profile_map <- function(ctx, cores, n, body) {
  worker_map(S7::set_props(ctx$spec, workers = as.integer(cores)), n, body,
             what = "profile points")
}


#' One Refit of a Profile
#'
#' @description
#' Fits the model with the labelled coefficients held at `val`. Without an
#' outer search it is the restricted fit of [statmod_restrict()] with several
#' coefficients held. With one, the search runs over the coordinates that the
#' criterion estimates and, when `hyper = "refit"`, over the
#' hyperparameters, with the held coefficients excluded from it.
#'
#' @details
#' A break-point of a continuous [modelterms7::seg()] term is held as in
#' [pin_held_breakpoints()], with one refit for each break-point that the held
#' position can stand for, the one with the smallest objective being kept.
#'
#' @param ctx The context of the call.
#' @param val The values of the held coefficients, in the order of
#'   `ctx$which`.
#' @param start `NULL` to start at the coefficients of the fit, or the list
#'   returned by a previous call, whose `par` and `hyper` are the start.
#'
#' @return A list with `par`, `hyper`, `objective`, `loglik`, `criterion` and
#'   `converged`.
#'
#' @keywords internal
profile_point <- function(ctx, val, start = NULL) {
  spec <- profile_hold(ctx$spec, ctx$param, ctx$coefname, val)
  cfg <- ctx$cfg
  beta_fit <- ctx$beta_fit
  beta0 <- if (is.null(start)) beta_fit else start$par
  beta0[ctx$at] <- val
  hyper0 <- if (is.null(start)) ctx$hyper0 else start$hyper
  run <- function(variant) {
    design <- statmod_design(spec)
    blocks <- statmod_blocks(spec, design)
    pins <- pin_held_breakpoints(spec, design, beta0, from = beta_fit,
                                 variant = variant)
    res <- if (ctx$use_outer) {
      outer_fit(spec, design, blocks, hyper0, ctx$inner, ctx$outer, NULL,
                pins$beta, ctx$emethod, cfg$maxit, cfg$tol, verbosity(0),
                hold_hyper = identical(ctx$hyper, "held"))
    } else {
      statmod_alternate(spec, design, blocks, hyper0, ctx$inner, pins$beta,
                        cfg$expected, cfg$emethod, cfg$maxit, cfg$tol,
                        verbosity(0))
    }
    list(design = design, pins = pins, res = res)
  }
  best <- run(1L)
  for (v in seq_len(best$pins$variants)[-1L]) {
    other <- run(v)
    if (is.finite(other$res$value) && other$res$value < best$res$value) {
      best <- other
    }
  }
  res <- best$res
  cf <- res$obj$split(res$par)
  list(par = res$par,
       hyper = if (ctx$use_outer && !is.null(res$hyper)) res$hyper else hyper0,
       objective = res$value,
       loglik = statmod_loglik_at(spec, cf, best$design),
       criterion = if (ctx$use_outer) as.numeric(res$criterion) else NA_real_,
       converged = isTRUE(res$converged) || isTRUE(res$settled))
}


#' A Failed Point of a Profile
#'
#' @param ctx The context of the call.
#' @param msg The message of the error.
#'
#' @return A point whose numbers are all missing, with the message in `error`.
#'
#' @keywords internal
profile_failed <- function(ctx, msg) {
  list(par = NULL, hyper = NULL, objective = NA_real_, loglik = NA_real_,
       criterion = NA_real_, converged = FALSE, error = msg)
}


#' The Chains of Points That a Start Visits
#'
#' @description
#' Splits the rows of the grid into the sequences along which a refit starts
#' from the previous one (`"neighbor"`) or into single points (`"fit"`).
#'
#' @param ctx The context of the call.
#' @param grid The data frame of held values.
#' @param start One of `"neighbor"`, `"fit"` and `"both"`.
#'
#' @return A list of chains, each a list with the row indices `rows` and
#'   whether the start `carry`ies from one to the next.
#'
#' @keywords internal
profile_chains <- function(ctx, grid, start) {
  n <- nrow(grid)
  single <- lapply(seq_len(n), function(i) list(rows = i, carry = FALSE))
  if (identical(start, "fit")) return(single)
  # chains along the first coefficient, outward from the estimate on each side
  k <- ncol(grid)
  keys <- if (k > 1L) do.call(paste, c(grid[-1L], sep = "\r")) else
    rep("", n)
  out <- list()
  for (key in unique(keys)) {
    rows <- which(keys == key)
    x <- grid[[1L]][rows]
    up <- rows[x >= ctx$estimate[1L]]
    up <- up[order(grid[[1L]][up])]
    dn <- rows[x < ctx$estimate[1L]]
    dn <- dn[order(-grid[[1L]][dn])]
    if (length(up)) out[[length(out) + 1L]] <- list(rows = up, carry = TRUE)
    if (length(dn)) out[[length(out) + 1L]] <- list(rows = dn, carry = TRUE)
  }
  if (identical(start, "both")) out <- c(out, single)
  out
}


#' @rdname profile_compute
#' @keywords internal
S7::method(profile_compute, ProfileGrid) <- function(method, ctx, grid, cores) {
  chains <- profile_chains(ctx, grid, method@start)
  if (cores > 1L) ctx$spec <- S7::set_props(ctx$spec, threads = 1L)
  vals <- as.matrix(grid)
  body <- function(i) {
    ch <- chains[[i]]
    prev <- NULL
    out <- vector("list", length(ch$rows))
    for (m in seq_along(ch$rows)) {
      r <- ch$rows[m]
      pt <- tryCatch(profile_point(ctx, vals[r, ], prev),
                     error = function(e) profile_failed(ctx, conditionMessage(e)))
      out[[m]] <- pt
      prev <- if (ch$carry && is.null(pt$error)) pt else NULL
    }
    out
  }
  res <- profile_map(ctx, cores, length(chains), body)
  pts <- vector("list", nrow(grid))
  for (i in seq_along(chains)) {
    for (m in seq_along(chains[[i]]$rows)) {
      r <- chains[[i]]$rows[m]
      pt <- res[[i]][[m]]
      cur <- pts[[r]]
      if (is.null(cur) || (is.null(pt$error) && (!is.null(cur$error) ||
          pt$objective < cur$objective))) {
        pts[[r]] <- pt
      }
    }
  }
  pts
}


#' Assemble the Table of a Profile
#'
#' @description
#' Converts the points returned by [profile_compute()] into the data frame of
#' [statmod_profile()], with the row of the unrestricted fit as the attribute
#' `maximum`.
#'
#' @details
#' The objective of a point is the negative of the penalized log-likelihood,
#' \eqn{-\ell + \tfrac12 P_\lambda}, so that the penalty column is
#' \eqn{-(\mathrm{objective} + \ell)} and the penalized log-likelihood is the
#' negative of the objective.
#'
#' @param ctx The context of the call.
#' @param grid The data frame of held values.
#' @param pts The points.
#' @param method The profile method.
#'
#' @return The data frame.
#'
#' @keywords internal
profile_assemble <- function(ctx, grid, pts, method) {
  one <- function(ll, obj, crit, conv) {
    data.frame(profile_loglik = ll, profile_penalty = -(obj + ll),
               profile_pen_loglik = -obj, outer_criterion = crit,
               converged = conv)
  }
  num <- function(f) vapply(pts, function(p) as.numeric(f(p)), numeric(1))
  out <- cbind(grid, one(num(function(p) p$loglik), num(function(p) p$objective),
                         num(function(p) p$criterion),
                         vapply(pts, function(p) isTRUE(p$converged),
                                logical(1))))
  names(out) <- c(names(grid), "profile_loglik", "profile_penalty",
                  "profile_pen_loglik", "outer_criterion", "converged")
  failed <- Filter(Negate(is.null), lapply(pts, `[[`, "error"))
  if (length(failed)) {
    warning(sprintf(paste0("%d of %d refits signaled an error and are ",
                           "reported as missing. The first: %s"),
                    length(failed), length(pts), failed[[1L]]), call. = FALSE)
  }
  fit <- ctx$fit
  mx <- cbind(as.data.frame(as.list(stats::setNames(ctx$estimate, ctx$which)),
                            check.names = FALSE),
              one(fit@loglik, fit@objective,
                  if (ctx$use_outer && length(fit@criterion)) fit@criterion[[1L]]
                  else NA_real_, isTRUE(fit@converged)))
  names(mx) <- names(out)
  attr(out, "maximum") <- mx
  attr(out, "method") <- method
  attr(out, "hyper") <- ctx$hyper
  attr(out, "criterion") <- ctx$outer
  out
}
