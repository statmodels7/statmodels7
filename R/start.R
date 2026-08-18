#' @include statmod.R
NULL

#' @title S7 Class for a Starting-Value Strategy
#'
#' @description
#' Where a fit begins, as an object rather than as a vector of numbers.
#'
#' @details
#' The reason it is an object is that the interesting answers are not values
#' but PROCEDURES: draw around the intercept-only fit, or search the
#' likelihood for a basin. A strategy is asked once, before the alternation
#' between the coefficients and the hyperparameters begins, which is what
#' separates it from an optimizer: an optimizer runs at every step of the fit,
#' a strategy runs once at its start.
#'
#' @param label A short name, used when printing.
#'
#' @return An S7 object of class \code{start_strategy}.
#'
#' @seealso \code{\link{start_intercepts}}, \code{\link{start_origin}},
#'   \code{\link{start_random}}
#'
#' @examples
#' start_origin()
#'
#' @export
start_strategy <- S7::new_class(
  "start_strategy",
  properties = list(label = S7::class_character),
  abstract = TRUE
)

#' The start_strategy Class Object
#'
#' @description
#' Fetched rather than captured, so that a test of membership cannot be fooled
#' by the class being re-created --- which is what \pkg{covr} does.
#'
#' @return The \code{\link{start_strategy}} class object.
#'
#' @keywords internal
start_strategy_class <- function() start_strategy


#' @title Where a Fit Begins
#'
#' @description
#' The one generic a starting-value strategy implements: given the model, it
#' returns one starting vector per distribution parameter.
#'
#' @param strategy A \code{\link{start_strategy}}.
#' @param spec The specification.
#' @param design The design.
#' @param obj The objective.
#' @param ... Passed to methods.
#'
#' @return A named list, one numeric vector per distribution parameter, each
#'   as long as that parameter's design.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + 2 * dd$x + rnorm(40, sd = 0.3)
#' spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
#' design <- statmod_design(spec)
#' # a strategy that needs no objective is asked with NULL for it; the ones
#' # that search the likelihood are handed the objective by statmod()
#' start_at(start_origin(), spec, design, NULL)
#' start_at(start_intercepts(), spec, design, NULL)
#'
#' @export
start_at <- S7::new_generic("start_at", "strategy",
  function(strategy, spec, design, obj, ...) S7::S7_dispatch())


#' @title S7 Classes for the Shipped Strategies
#' @description The classes the constructors below instantiate.
#' @param fn The generator a random start draws from.
#' @param args Further arguments to it.
#' @param center Whether the draw is added to the intercept-only start.
#' @return An S7 object inheriting from \code{\link{start_strategy}}.
#' @name StartIntercepts-class
#' @aliases StartIntercepts StartOrigin StartRandom
#' @keywords internal
StartIntercepts <- S7::new_class("StartIntercepts", parent = start_strategy)

#' @rdname StartIntercepts-class
#' @keywords internal
StartOrigin <- S7::new_class("StartOrigin", parent = start_strategy)

#' @rdname StartIntercepts-class
#' @keywords internal
StartRandom <- S7::new_class("StartRandom", parent = start_strategy,
  properties = list(fn = S7::class_function, args = S7::class_list,
                    center = S7::class_logical))


#' @title Start at the Intercept-Only Fit
#'
#' @description
#' The default: each equation's intercept at the maximum likelihood estimate
#' of the model with every covariate removed, and every other coefficient at
#' zero.
#'
#' @details
#' The model with no covariates is the same model with every slope set to
#' zero, so it is exactly where the model with them should begin, and it costs
#' one small fit. The penalized blocks start at zero, where their penalty is
#' smallest. This is what \code{statmod()} does when \code{start} is
#' \code{NULL}; the constructor exists so that the default can be named,
#' compared against and passed on.
#'
#' @return A \code{\link{start_strategy}}.
#'
#' @examples
#' start_intercepts()
#'
#' @seealso \code{\link{start_origin}}, \code{\link{start_random}}
#' @export
start_intercepts <- function() {
  StartIntercepts(label = "intercept-only fit")
}

#' @title Start at the Origin
#'
#' @description
#' Every coefficient at zero on the unconstrained scale, which for each
#' distribution parameter is the value its link maps zero to.
#'
#' @details
#' The name is \code{start_origin} and not \code{start_zeros} because
#' \pkg{optimizers7} already exports the second for the starting POINT of an
#' optimizer, and the toolkit's members share no exported name --- a
#' collision would be reported by \code{\link{statmodels7_conflicts}} and, far
#' worse, would mean that which function a user got depended on the order the
#' packages were attached in.
#'
#' It is the plainest starting point and rarely the best: a location parameter
#' on the identity link then begins at zero whatever the response is, which on
#' a response centred at a thousand sends the run travelling. It is here
#' because it is the reference other strategies are judged against, and
#' because a model whose equations are all on a log or logit link is not
#' harmed by it.
#'
#' @return A \code{\link{start_strategy}}.
#'
#' @examples
#' start_origin()
#'
#' @seealso \code{\link{start_intercepts}}
#' @export
start_origin <- function() {
  StartOrigin(label = "zero")
}

#' @title Start From a Random Draw
#'
#' @description
#' Every coefficient drawn on the unconstrained scale, by default added to the
#' intercept-only fit rather than replacing it.
#'
#' @details
#' The draw is added to \code{\link{start_intercepts}}'s answer unless
#' \code{center = FALSE}, and that is not a detail: a coefficient drawn from a
#' standard normal is a sensible perturbation and a hopeless absolute value,
#' since the intercept of a location equation is on the scale of the response.
#' Centring keeps the scale and randomizes the direction, which is what a
#' caller wanting several starts is after.
#'
#' The stream is the caller's, so \code{set.seed()} governs the result and a
#' fit begun this way is reproducible only alongside its seed.
#'
#' @param fn A generator called as \code{fn(n, ...)}, returning \code{n}
#'   values. Defaults to \code{\link[stats]{rnorm}}.
#' @param ... Further arguments to \code{fn}, such as \code{sd} or
#'   \code{min} and \code{max}.
#' @param center Whether to add the draw to the intercept-only start rather
#'   than use it alone. Defaults to \code{TRUE}.
#'
#' @return A \code{\link{start_strategy}}.
#'
#' @examples
#' start_random()
#' start_random(stats::runif, min = -2, max = 2)
#'
#' @seealso \code{\link{start_intercepts}},
#'   \code{\link[optimizers7]{multistart}}
#' @export
start_random <- function(fn = stats::rnorm, ..., center = TRUE) {
  if (!is.function(fn)) {
    stop("'fn' must be a function called as fn(n, ...).", call. = FALSE)
  }
  if (length(center) != 1L || !is.logical(center) || is.na(center)) {
    stop("'center' must be TRUE or FALSE.", call. = FALSE)
  }
  StartRandom(label = if (center) "random around the intercept-only fit"
              else "random", fn = fn, args = list(...), center = center)
}


#' @rdname StartIntercepts-class
#' @keywords internal
StartSearch <- S7::new_class("StartSearch", parent = start_strategy,
  properties = list(optimizer = S7::class_any, over = S7::class_any,
                    hyper = S7::class_logical))


#' @title Search the Likelihood for a Starting Point
#'
#' @description
#' Runs a global search ONCE, before the fit begins, over the coordinates
#' where the likelihood is not convex.
#'
#' @details
#' \strong{Once, and not inside the fit.} A search belongs to the starting
#' value and not to the scoring step. Handed to \code{inner_optimizer} instead
#' --- as \code{chain(sa(), iwls())} --- it would rerun at every
#' hyperparameter the outer criterion tried, which on an ordinary fit is 46
#' times and inside \code{\link{cv}} is the folds times the path, each time
#' returning the same answer.
#'
#' \strong{On the likelihood alone.} The penalties are off. What a starting
#' value has to get right is the BASIN of the likelihood, which is the one
#' thing the fit will not correct by itself; the penalties enter afterwards,
#' when their hyperparameters are estimated, and at the probe values they
#' represent nobody's choice --- a lasso at \eqn{\lambda = 1} against an
#' unaveraged log-likelihood empties whole blocks and would search a different
#' model.
#'
#' \strong{Over the coordinates where the problem is not convex}, which the
#' toolkit already knows how to name. They are the parameters of a structural
#' term (\code{\link[modelterms7]{gas}}, \code{\link[modelterms7]{regime}}),
#' the coefficients of a term that recomputes its own block
#' (\code{\link[modelterms7]{nl}}, \code{\link[modelterms7]{seg}} and the
#' break-point terms), and each equation's intercept. Everything else keeps
#' the default: a smooth, a ridge or a random effect is a convex block whose
#' optimum the scoring step reaches from anywhere, and searching over a
#' thousand random-effect coefficients would spend the whole budget on the
#' one part of the model that does not need it. \code{over} overrides the
#' choice by name.
#'
#' \strong{The hyperparameters} are left where they were unless
#' \code{hyper = TRUE}, which extends the search to the smooth ones on their
#' log scale. It is off by default because each of those coordinates costs a
#' full refit at every proposal rather than one likelihood evaluation.
#' Kinked penalties are never searched: their hyperparameter has a known
#' upper end and is swept by a warm-started path, which a random jump would
#' both fail to improve on and destroy.
#'
#' @param optimizer The optimizer to search with. Defaults to
#'   \code{optimizers7::sa()}.
#' @param over Optional names of the coefficients to search over, overriding
#'   the choice described above.
#' @param hyper Whether to search the smooth hyperparameters too. Defaults to
#'   \code{FALSE}.
#'
#' @return A \code{\link{start_strategy}}.
#'
#' @examples
#' start_search()
#' start_search(optimizers7::sa(maxit = 20))
#'
#' @seealso \code{\link{start_intercepts}},
#'   \code{\link[optimizers7]{sa}}, \code{\link[optimizers7]{chain}}
#' @export
start_search <- function(optimizer = optimizers7::sa(), over = NULL,
                         hyper = FALSE) {
  if (!S7::S7_inherits(optimizer, optimizers7::optimizer)) {
    stop("'optimizer' must be an optimizers7 optimizer, e.g. sa().",
         call. = FALSE)
  }
  if (!is.null(over) && (!is.character(over) || !length(over))) {
    stop("'over' must be a character vector of coefficient names, or NULL.",
         call. = FALSE)
  }
  if (length(hyper) != 1L || !is.logical(hyper) || is.na(hyper)) {
    stop("'hyper' must be TRUE or FALSE.", call. = FALSE)
  }
  StartSearch(label = paste0("search with ", optimizer@name),
              optimizer = optimizer, over = over, hyper = hyper)
}


#' @title Print a Starting-Value Strategy
#' @name print.start_strategy
#' @param x A \code{\link{start_strategy}}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @examples
#' print(start_random())
#' @keywords internal
S7::method(print, start_strategy) <- function(x, ...) {
  cat("<start> ", x@label, "\n", sep = "")
  invisible(x)
}


#' @title Starting Values From the Intercept-Only Fit
#' @name start_at.StartIntercepts
#' @description The default strategy's answer.
#' @param strategy A \code{StartIntercepts} object.
#' @param spec,design,obj,... As in \code{\link{start_at}}.
#' @return A named list of numeric vectors.
#' @keywords internal
S7::method(start_at, StartIntercepts) <-
  function(strategy, spec, design, obj, ...) {
    params <- spec@distrib@params
    out <- stats::setNames(lapply(design, function(d) numeric(d$npar)),
                           params)
    eta0 <- statmod_intercepts(spec)
    for (p in params) {
      if (design[[p]]$npar == 0L) next
      v <- eta0[[p]]
      if (is.null(v) || !is.finite(v)) next
      # AN OFFSET IS PART OF THE PREDICTOR AND THE INTERCEPT-ONLY FIT DOES NOT
      # SEE IT. statmod_intercepts() fits the distribution to the response
      # alone, so it answers with the predictor the model should have on
      # average; the intercept has to carry that MINUS what the offset already
      # contributes, or the run starts wherever the offset happens to sit.
      # On a count model over person-years the offset averages 7.56, so an
      # uncorrected intercept put the starting mean at exp(7.47) = 1750 where
      # the data average 0.92 -- out by a factor of nineteen hundred, and the
      # scoring iteration spent its budget crawling back.
      off <- spec@offsets[[p]]
      if (!is.null(off)) {
        m <- mean(off)
        if (is.finite(m)) v <- v - m
      }
      # the intercept carries the whole of it when there is one
      ii <- parametric_intercept(spec, design, p)
      if (!is.na(ii)) out[[p]][ii] <- v
    }
    # and each term says where its own block begins. The base method of
    # term_coef_start() is zero everywhere, so an ordinary block is
    # unchanged; a term that recomputes its design from its coefficients
    # answers with the start it built itself, because zero is degenerate
    # there rather than neutral. A jump reads its break-point off
    # -g_k/delta_k, which at zero puts every break-point at the same
    # clamped position and leaves the block singular, and nothing in a
    # scoring step moves it off that point.
    for (p in params) {
      if (design[[p]]$npar == 0L) next
      # THE RESPONSE ON THE SCALE OF THE PREDICTOR, which is what a term with
      # parameters of its own needs to estimate them and the one thing it
      # cannot work out: the term knows its formula and its charts, the layer
      # knows the distribution, the link and the equation. It exists only
      # where the response reads the parameter directly -- a mean or a
      # location -- and a term in a scale's equation is handed nothing rather
      # than a quantity invented for it.
      tg <- predictor_target(spec, p)
      for (nm in names(spec@terms[[p]])) {
        idx <- design[[p]]$blocks[[nm]]
        if (is.null(idx) || !length(idx)) next
        v <- tryCatch(modelterms7::term_coef_start(spec@terms[[p]][[nm]],
                                                   target = tg),
                      error = function(e) NULL)
        if (is.null(v) || length(v) != length(idx) || !all(is.finite(v))) next
        # A term asking for zeros is asking for nothing, and writing them
        # would undo the intercept set above: the intercept is a column of
        # the parametric block, whose start is zero everywhere by the base
        # method. Only a term that wants something says so.
        if (all(v == 0)) next
        out[[p]][idx] <- as.numeric(v)
      }
    }
    out
  }

#' @title Where an Equation's Intercept Is
#' @name parametric_intercept
#'
#' @description
#' The position of the intercept of one equation's PARAMETRIC block, or
#' \code{NA} where it has none.
#'
#' @details
#' It is a column of the parametric block and not merely a coefficient whose
#' name ends in \code{(Intercept)}. \code{\link[modelterms7]{nl}} names the
#' intercept of each of its own parameters the same way, and those live on
#' those parameters' charts rather than on the predictor's, so a model
#' written \code{y ~ 0 + nl(...)} puts one of them first. Writing the
#' intercept-only fit there sets the parameter to \code{linkinv} of a value
#' that was never on its scale: measured on a logistic growth curve whose
#' asymptote rides a log link, \code{mean(y) = 23.9} became a starting
#' \eqn{\phi} of \code{2.5e10}, an objective of \code{7.0e20} and a gradient
#' of \code{1.4e21}, on data whose every scale is ordinary. The lasso path
#' built at those coefficients then spanned \code{2.8e15} to \code{2.8e19}
#' where the block empties at about 300, so all of its points were the same
#' empty fit and every subject deviation was estimated as exactly zero.
#'
#' @param spec The specification.
#' @param design The design.
#' @param p The distribution parameter naming the equation.
#'
#' @return An integer position into the equation's coefficient vector, or
#'   \code{NA_integer_}.
#'
#' @seealso \code{\link{start_at}}, \code{\link{statmod_intercepts}}
#' @keywords internal
parametric_intercept <- function(spec, design, p) {
  for (nm in names(spec@terms[[p]])) {
    tm <- spec@terms[[p]][[nm]]
    if (!S7::S7_inherits(tm, modelterms7::LinparTerm)) next
    idx <- design[[p]]$blocks[[nm]]
    if (is.null(idx) || !length(idx)) next
    cn <- tryCatch(modelterms7::term_coef_names(tm), error = function(e) NULL)
    if (is.null(cn) || length(cn) != length(idx)) next
    j <- which(cn == "(Intercept)" | grepl("\\.\\(Intercept\\)$", cn))
    if (length(j)) return(as.integer(idx[j[1L]]))
  }
  NA_integer_
}

#' @title The Response on the Scale of a Predictor
#' @name predictor_target
#'
#' @description
#' The response carried onto the scale of one equation's linear predictor,
#' which is what \code{\link[modelterms7]{term_coef_start}} needs to estimate
#' a term's own parameters from the data.
#'
#' @details
#' It exists only where the response reads the parameter directly, which
#' \code{params_interpretation} says: a mean or a location. For a scale or a
#' shape there is no per-observation reading of the parameter, and
#' \code{NULL} is returned rather than a quantity invented for the occasion.
#'
#' The scale matters and is not a detail. Measured on a Poisson whose
#' predictor is a logistic growth curve with \eqn{\phi = 4}, a term handed
#' the raw response estimated \eqn{\phi} between 52.7 and 54.8 over five
#' samples, and one handed \eqn{g(y)} estimated it between 4.03 and 4.07.
#' What does NOT matter, measured on the same shape beside another term, is
#' the other terms' contribution: subtracting it moved the estimate from
#' 39.68 to 39.87 against a truth of 40, so the response is passed as it
#' stands and nothing is residualized.
#'
#' A response on a bound of the parameter's own domain is moved half way to
#' the nearest admissible value, since a link is not defined at its bound --
#' for counts under a log link that is the classical half.
#'
#' @param spec The specification.
#' @param p The distribution parameter naming the equation.
#'
#' @return A numeric vector, one value per observation, or \code{NULL}.
#'
#' @seealso \code{\link{start_at}},
#'   \code{\link[modelterms7]{term_coef_start}}
#' @keywords internal
predictor_target <- function(spec, p) {
  d <- spec@distrib
  if (!("params_interpretation" %in% S7::prop_names(d))) return(NULL)
  ip <- d@params_interpretation[[p]]
  if (is.null(ip) || !ip %in% c("mean", "location")) return(NULL)
  y <- suppressWarnings(as.numeric(spec@response))
  if (!length(y) || length(y) != spec@n_obs) return(NULL)
  lk <- d@link_params[[p]]
  if (is.null(lk)) return(NULL)
  b <- lk@link_bounds
  inside <- is.finite(y) & y > b[1L] & y < b[2L]
  if (!any(inside)) return(NULL)
  if (is.finite(b[1L])) {
    lo <- b[1L] + (min(y[inside]) - b[1L]) / 2
    y[is.finite(y) & y <= b[1L]] <- lo
  }
  if (is.finite(b[2L])) {
    hi <- b[2L] - (b[2L] - max(y[inside])) / 2
    y[is.finite(y) & y >= b[2L]] <- hi
  }
  tg <- suppressWarnings(linkfunctions7::linkfun(lk, y))
  if (!any(is.finite(tg))) return(NULL)
  off <- spec@offsets[[p]]
  if (!is.null(off) && length(off) == length(tg)) tg <- tg - off
  tg
}

#' @title Starting Values at Zero
#' @name start_at.StartOrigin
#' @description Every coefficient at zero.
#' @param strategy A \code{StartOrigin} object.
#' @param spec,design,obj,... As in \code{\link{start_at}}.
#' @return A named list of numeric vectors.
#' @keywords internal
S7::method(start_at, StartOrigin) <-
  function(strategy, spec, design, obj, ...) {
    stats::setNames(lapply(design, function(d) numeric(d$npar)),
                    spec@distrib@params)
  }

#' @title Starting Values From a Random Draw
#' @name start_at.StartRandom
#' @description The draw, added to the intercept-only fit unless told not to.
#' @param strategy A \code{StartRandom} object.
#' @param spec,design,obj,... As in \code{\link{start_at}}.
#' @return A named list of numeric vectors.
#' @keywords internal
S7::method(start_at, StartRandom) <-
  function(strategy, spec, design, obj, ...) {
    base <- if (strategy@center) {
      start_at(start_intercepts(), spec, design, obj)
    } else {
      stats::setNames(lapply(design, function(d) numeric(d$npar)),
                      spec@distrib@params)
    }
    for (p in names(base)) {
      n <- length(base[[p]])
      if (!n) next
      draw <- do.call(strategy@fn, c(list(n), strategy@args))
      if (length(draw) != n || !is.numeric(draw)) {
        stop(sprintf("'fn' must return %d numeric values for '%s'.", n, p),
             call. = FALSE)
      }
      base[[p]] <- base[[p]] + as.numeric(draw)
    }
    base
  }


#' Which Coefficients a Search Should Cover
#'
#' @description
#' The positions, in the stacked coefficient vector, of the coordinates where
#' the likelihood is not convex: each equation's intercept, and the blocks of
#' the terms that recompute their own design.
#'
#' @details
#' The predicate for the second is \code{\link{refreshes_own_block}}, the same
#' one \code{\link{unfittable_reason}} uses, so a term written later is
#' covered without an edit here. A convex block is left out deliberately
#' rather than forgotten: the scoring step reaches its optimum from anywhere,
#' and a search over it would spend the budget where it buys nothing.
#'
#' @param spec The specification.
#' @param design The design.
#'
#' @return An integer vector of positions, possibly empty.
#'
#' @keywords internal
search_coords <- function(spec, design) {
  params <- spec@distrib@params
  npar <- vapply(design[params], function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  idx <- integer(0)
  for (a in seq_along(params)) {
    p <- params[a]
    if (design[[p]]$npar == 0L) next
    nms <- design[[p]]$coef_names
    # the intercept, wherever the equation puts it
    hit <- which(nms == "(Intercept)" | grepl("\\(Intercept\\)$", nms))
    for (nm in names(spec@terms[[p]])) {
      if (refreshes_own_block(spec@terms[[p]][[nm]])) {
        hit <- c(hit, design[[p]]$blocks[[nm]])
      }
    }
    idx <- c(idx, offs[a] + sort(unique(hit)))
  }
  sort(unique(idx))
}


#' @title Starting Values From a Global Search
#' @name start_at.StartSearch
#' @description
#' Runs the optimizer once on the likelihood, over the non-convex coordinates
#' and, where the model has one, over a structural term's own parameters.
#' @details
#' A structural term's parameters do not live in the coefficient vector; they
#' are held in the design's structural state, and the objective reads them
#' from there. The search therefore sets them into that state and evaluates
#' the likelihood, exactly as \code{statmod_fit_joint()} does, and leaves the
#' best it found in place --- which is how a strategy returning coefficients
#' can nonetheless start a filter somewhere better than its own
#' \code{term_start()}.
#' @param strategy A \code{StartSearch} object.
#' @param spec,design,obj,... As in \code{\link{start_at}}.
#' @return A named list of numeric vectors.
#' @keywords internal
S7::method(start_at, StartSearch) <-
  function(strategy, spec, design, obj, ...) {
    base <- start_at(start_intercepts(), spec, design, obj)
    beta0 <- obj$stack(base)

    idx <- if (is.null(strategy@over)) {
      search_coords(spec, design)
    } else {
      nms <- unlist(lapply(spec@distrib@params,
                           function(p) design[[p]]$coef_names),
                    use.names = FALSE)
      hit <- match(strategy@over, nms)
      if (anyNA(hit)) {
        stop(sprintf("'over' names %s, which is not a coefficient.",
                     paste(strategy@over[is.na(hit)], collapse = ", ")),
             call. = FALSE)
      }
      sort(hit)
    }

    # a structural term's own parameters, which the objective reads from the
    # design's state rather than from the coefficient vector
    sst <- statmod_structural_state(design)
    su <- attr(design, "structural")
    zkey <- NULL
    zfree <- character(0)
    if (length(su)) {
      zkey <- su[[1L]]$term
      zfree <- setdiff(names(sst$zeta[[zkey]]), sst$held[[zkey]])
    }
    setz <- function(z) {
      v <- sst$zeta[[zkey]]
      v[zfree] <- as.numeric(z)
      sst$zeta[[zkey]] <- v
      sst$key <- NULL
      sst$value <- NULL
    }

    nz <- length(zfree)
    if (!length(idx) && !nz) return(base)
    z0 <- if (nz) as.numeric(sst$zeta[[zkey]][zfree]) else numeric(0)
    u0 <- c(beta0[idx], z0)

    # the LIKELIHOOD alone: what a starting value has to get right is the
    # basin, and the penalties enter when their hyperparameters are estimated
    fn <- function(u) {
      v <- tryCatch({
        b <- beta0
        if (length(idx)) b[idx] <- u[seq_along(idx)]
        if (nz) setz(u[length(idx) + seq_len(nz)])
        -statmod_loglik_at(spec, obj$split(b), design)
      }, error = function(e) NA_real_)
      if (!is.finite(v)) Inf else v
    }

    if (!is.finite(fn(u0))) {
      # the search cannot start where the model cannot be evaluated, and
      # saying so here names the strategy rather than leaving an optimizer to
      # report an objective it was handed
      stop("the likelihood is not finite at the intercept-only start, so ",
           "there is nothing to search from.", call. = FALSE)
    }
    res <- optimizers7::minimize(strategy@optimizer, fn, u0)

    # the state keeps whatever the search left in it, so the best point found
    # is where the filter starts; the coefficients go back as a list
    if (nz) setz(res@par[length(idx) + seq_len(nz)])
    if (length(idx)) beta0[idx] <- res@par[seq_along(idx)]
    obj$split(beta0)
  }
