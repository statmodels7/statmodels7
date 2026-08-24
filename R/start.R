#' @include statmod.R
NULL

#' @title S7 Class for a Starting-Value Strategy
#'
#' @description
#' The abstract parent of the starting-value strategies. A strategy says where
#' a fit begins as a **procedure**, not as a vector of numbers: fit the
#' intercept-only model, draw around it, or search the likelihood for a basin.
#' `statmod(start =)` accepts one, or a plain named list of coefficients.
#'
#' @details
#' A strategy is asked exactly once, before the alternation between the
#' coefficients and the hyperparameters begins. An optimizer runs at every
#' step of that alternation, which is why a global search belongs here: given
#' to `inner_optimizer` it would rerun at every hyperparameter the outer
#' search tried and return the same answer each time.
#'
#' **The class is abstract and cannot be instantiated.** Calling
#' `start_strategy()` signals an error; what the class is for is to be
#' subclassed and to make a membership test possible. To write a strategy of
#' your own, subclass it and register a [start_at()] method on the subclass.
#' That method receives the specification, the design and the objective, and
#' returns one vector per distribution parameter.
#'
#' @param label A short name, printed by [print.start_strategy()]. A single
#'   string. Set by each subclass's constructor, since the class itself
#'   cannot be built.
#'
#' @return Nothing: the class is abstract and the constructor signals an
#'   error. Used as a class object, for `S7::S7_inherits()` and for
#'   registering methods, it has one property, `label`, which every subclass
#'   inherits.
#'
#' @seealso [start_intercepts()] (the default), [start_origin()],
#'   [start_random()] and [start_search()] for the four shipped strategies,
#'   [start_at()] for the generic they implement.
#'
#' @examples
#' # The four shipped strategies all inherit from this.
#' S7::S7_inherits(start_origin(), start_strategy)
#' S7::S7_inherits(start_search(), start_strategy)
#'
#' # Every one of them carries the label the class defines.
#' vapply(list(start_intercepts(), start_origin(), start_random(),
#'             start_search()),
#'        function(s) s@label, character(1))
#'
#' # The class itself is abstract.
#' try(start_strategy("my own"))
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
#' Fetches the [start_strategy()] class object from the namespace at the
#' moment it is asked for, so a membership test reads the class as it now
#' exists.
#'
#' @details
#' A class captured in a variable at build time makes such a test fail under
#' \pkg{covr}, which re-evaluates a package's code and so re-creates the
#' class. `S7::S7_inherits()` against the captured copy is then
#' `FALSE` for an object of the live class. `distributions7` and
#' `linkfunctions7` record the same trap.
#'
#' @return The `start_strategy` S7 class object.
#'
#' @seealso [start_strategy()] for the class itself.
#'
#' @keywords internal
start_strategy_class <- function() start_strategy


#' @title Where a Fit Begins
#'
#' @description
#' The one generic a starting-value strategy implements. Given the model, it
#' returns one starting vector per distribution parameter, on the coefficient
#' scale, ready for the fit to begin from.
#'
#' @details
#' Write a method on your own subclass of [start_strategy()] to add a
#' strategy. The four shipped methods show the range: [start_origin()] reads
#' only the design's widths, [start_intercepts()] fits a small model,
#' [start_random()] draws from the caller's generator, and [start_search()]
#' runs an optimizer over `obj`.
#'
#' A method must return a full-length vector for every distribution
#' parameter, including the ones it has nothing to say about. Zero is the
#' conventional filler, and it is the value a penalized block wants anyway,
#' its penalty being smallest there.
#'
#' @param strategy A [start_strategy()], which decides the method.
#' @param spec The [StatmodSpec()] being fitted.
#' @param design The design, as [statmod_design()] returns it.
#' @param obj The objective, as [statmod_objective()] returns it, or `NULL`.
#'   Only [start_search()] reads it; the other three accept `NULL`.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A named list, one numeric vector per distribution parameter in the
#'   family's order, each as long as that parameter's design is wide.
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
#' @description
#' The three classes [start_intercepts()], [start_origin()] and
#' [start_random()] instantiate. Each inherits from [start_strategy()] and
#' carries its `label`; only `StartRandom` adds properties of its own. Use the
#' constructors, which validate; these raw ones do not.
#' @param fn The generator a random start draws from, called as `fn(n, ...)`.
#'   `StartRandom` only.
#' @param args A list of further arguments to `fn`. `StartRandom` only.
#' @param center A single logical: whether the draw is added to the
#'   intercept-only start. `StartRandom` only.
#' @return An S7 object inheriting from [start_strategy()]:
#'   `StartIntercepts` and `StartOrigin` with the one inherited property
#'   `label`, `StartRandom` with `fn`, `args` and `center` besides.
#' @seealso [start_intercepts()], [start_origin()], [start_random()] for the
#'   constructors, [start_at()] for the methods registered on these classes.
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
#' The model with no covariates is the model with every slope set to zero, so
#' its maximum likelihood estimate is exactly where the model with covariates
#' should begin. It costs one small fit, on an intercept-only design.
#'
#' Every other coefficient starts at zero, which is where a penalized block's
#' penalty is smallest.
#'
#' Two adjustments make it usable where a naive intercept would not be:
#'
#' - An **offset** is subtracted before the value is written. Without that, a
#'   model with an offset averaging 6.74 begins at
#'   \eqn{\exp(0.897 + 6.74) = 2080} against a sample mean of 2.45. Measured
#'   on a negative binomial over person-years: 4.9 s and 9 scoring iterations
#'   with the correction, more than 25 minutes without.
#' - The value is written to a **parametric intercept** and to nothing else.
#'   A term such as [modelterms7::nl()] names the intercept of each of its own
#'   parameters `(Intercept)` too, and those live on that parameter's chart,
#'   not the predictor's. Writing a predictor-scale value onto a `phi` that
#'   rides a log link gave \eqn{\exp(23.9) = 2.5 \times 10^{10}}.
#'
#' This is what [statmod()] does when `start` is `NULL`. The constructor
#' exists so the default can be named, compared against and passed on.
#'
#' @return A `StartIntercepts` object, inheriting from [start_strategy()].
#'
#' @examples
#' start_intercepts()
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(60, 0, 10))
#' dd$y <- 500 + 20 * dd$x + rnorm(60, sd = 5)
#' spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
#'
#' # The intercept starts at the response's own scale, and the slope at zero.
#' start_at(start_intercepts(), spec, statmod_design(spec), NULL)
#'
#' @seealso [start_origin()], [start_random()], [start_search()]
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
#' # When it is a poor choice
#'
#' A location parameter on the identity link begins at zero whatever the
#' response is. On a response centered at a thousand the run then starts a
#' thousand units away and has to travel there, and on a badly scaled model
#' it may not arrive at all. [start_intercepts()], the default, exists for
#' that reason.
#'
#' It is harmless where every equation is on a log or a logit link, since
#' zero is then an ordinary interior value of the parameter, and it is the
#' reference the other strategies are judged against.
#'
#' # The name
#'
#' `start_origin`, not `start_zeros`. \pkg{optimizers7} already exports
#' `start_zeros()` for the starting point of an optimizer, and no two members
#' of this toolkit export the same name. A collision would be reported by
#' [statmodels7_conflicts()], and, worse, which function a caller reached
#' would depend on the order the packages were attached in.
#'
#' @return A `StartOrigin` object, inheriting from [start_strategy()].
#'
#' @examples
#' start_origin()
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 1 + 2 * dd$x + rnorm(40, sd = 0.3)
#' spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
#'
#' # Every coefficient of every equation is zero, intercepts included.
#' s <- start_at(start_origin(), spec, statmod_design(spec), NULL)
#' s
#' all(unlist(s) == 0)
#'
#' @seealso [start_intercepts()], the default and usually the better choice.
#' @export
start_origin <- function() {
  StartOrigin(label = "zero")
}

#' @title Start From a Random Draw
#'
#' @description
#' Draws every coefficient from a generator of the caller's choosing and, by
#' default, adds the draw to the intercept-only fit. Built for running a model
#' from several starts and comparing where they land, which is how a
#' multimodal likelihood is explored.
#'
#' @details
#' # Why the draw is centered by default
#'
#' A coefficient drawn from a standard normal is a sensible perturbation and a
#' hopeless absolute value: the intercept of a location equation is on the
#' scale of the response, which may be a thousand. Adding the draw to
#' [start_intercepts()]'s answer keeps that scale and randomizes the
#' direction, leaving the several starts something to differ in.
#'
#' `center = FALSE` uses the draw alone, which is right when the generator is
#' already scaled to the problem.
#'
#' # The random stream
#'
#' The draw comes from the session's generator, so [set.seed()] governs the
#' result and a fit begun this way repeats only alongside its seed. Nothing
#' is seeded internally.
#'
#' @param fn A generator called as `fn(n, ...)` and returning `n` values.
#'   Defaults to [stats::rnorm()]. Any function of that shape serves.
#' @param ... Further arguments to `fn`, captured at construction and stored
#'   on the object: `sd = 0.1`, or `min` and `max` for [stats::runif()].
#' @param center A single logical. `TRUE`, the default, adds the draw to the
#'   intercept-only start; `FALSE` uses the draw alone.
#'
#' @return A `StartRandom` object, inheriting from [start_strategy()], with
#'   properties `fn`, `args` and `center`.
#'
#' @examples
#' start_random()
#' start_random(stats::runif, min = -2, max = 2)
#' start_random(sd = 0.1)@args
#'
#' set.seed(1)
#' dd <- data.frame(x = runif(40))
#' dd$y <- 100 + 2 * dd$x + rnorm(40, sd = 0.3)
#' spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
#' design <- statmod_design(spec)
#'
#' # Centered, the intercept stays on the response's scale.
#' set.seed(2)
#' start_at(start_random(sd = 0.5), spec, design, NULL)$mu
#'
#' # Uncentered, it is the raw draw, which is nowhere near 100.
#' set.seed(2)
#' start_at(start_random(sd = 0.5, center = FALSE), spec, design, NULL)$mu
#'
#' @seealso [start_intercepts()] for what the draw is added to,
#'   [optimizers7::multistart()] for running several starts inside an
#'   optimizer.
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
  properties = list(optimizer = S7::class_any, over = S7::class_any))


#' @title Search the Likelihood for a Starting Point
#'
#' @description
#' Runs a global search once, before the fit begins, over the coordinates
#' where the likelihood is not convex. Built for a model whose objective has
#' several local optima: a break-point term, a nonlinear term, a score-driven
#' filter or a latent-state mixture.
#'
#' @details
#' # Once, before the fit
#'
#' A search belongs to the starting value, never to the scoring step. Given
#' to `inner_optimizer` instead, as `chain(sa(), iwls())`, it would rerun at
#' every hyperparameter the outer criterion tried, which on an ordinary fit
#' is 46 times and inside [cv()] is the folds times the path, returning the
#' same answer each time.
#'
#' # On the likelihood alone
#'
#' The penalties are off. What a starting value has to get right is the
#' basin of the likelihood, the one thing the fit will not correct by itself.
#' The penalties enter afterwards, when their hyperparameters are estimated,
#' and at the probe values they represent nobody's choice: a lasso at
#' \eqn{\lambda = 1} against an unaveraged log-likelihood empties whole
#' blocks, so a search there would explore a different model.
#'
#' # Which coordinates are searched
#'
#' The non-convex ones, which the toolkit can already name: the parameters of
#' a structural term ([modelterms7::gas()], [modelterms7::regime()]), the
#' coefficients of a term that recomputes its own block
#' ([modelterms7::nl()], [modelterms7::seg()] and the break-point terms), and
#' each equation's intercept.
#'
#' Everything else keeps its default. A smooth, a ridge or a random effect is
#' a convex block whose optimum the scoring step reaches from anywhere, and
#' searching over a thousand random-effect coefficients would spend the whole
#' budget on the part of the model that needs it least.
#'
#' A penalized coordinate is excluded **wherever it sits**, including inside
#' a non-convex term. A sub-formula develops a break-point or a nonlinear
#' parameter over groups, and those deviations are columns of the term's own
#' block. They are the case the rule exists for: on the likelihood alone
#' nothing identifies them, since the penalty is what does, so a search over
#' them fits each group's own points and moves away from the penalized mode.
#' Measured on `jseg(x, npsi = 2, by = ~random(~1|id))` over thirty groups,
#' 210 of the 219 coordinates are such deviations.
#'
#' `over` overrides the choice by name.
#'
#' # The hyperparameters are not searched
#'
#' They cannot be, from here. The objective is the likelihood with the
#' penalties off, in which a hyperparameter does not appear, so no proposal
#' could change one.
#'
#' A global search over them is a search over the **outer** criterion, where
#' each proposal costs a full refit instead of one likelihood evaluation, and
#' that is `statmod(outer_optimizer = optimizers7::sa())`. A kinked penalty
#' is outside even that: its hyperparameter has a known upper end and is
#' swept by a warm-started path, which a random jump would fail to improve on
#' and would destroy.
#'
#' @param optimizer The optimizer to search with, any \pkg{optimizers7}
#'   optimizer. Defaults to `optimizers7::sa()`, simulated annealing, which
#'   is the one built for a multimodal surface.
#' @param over Optional character vector naming the coefficients to search
#'   over, overriding the choice above. Names not in the model are an error.
#'
#' @return A `StartSearch` object, inheriting from [start_strategy()], with
#'   properties `optimizer` and `over`.
#'
#' @examples
#' start_search()
#' start_search(optimizers7::sa(maxit = 20))
#'
#' # It changes nothing on a convex fit: there are no non-convex coordinates
#' # to search but the intercepts, and the scoring step reaches those anyway.
#' set.seed(1)
#' dd <- data.frame(x = runif(60))
#' dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.3)
#' a <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' b <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd,
#'              start = start_search(optimizers7::sa(maxit = 30)))
#' max(abs(unlist(a@coefficients) - unlist(b@coefficients)))
#'
#' @seealso [start_intercepts()] for the default,
#'   [optimizers7::sa()] for the search, [statmod()] for the
#'   `outer_optimizer` argument that searches the hyperparameters instead.
#' @export
start_search <- function(optimizer = optimizers7::sa(), over = NULL) {
  if (!S7::S7_inherits(optimizer, optimizers7::optimizer)) {
    stop("'optimizer' must be an optimizers7 optimizer, e.g. sa().",
         call. = FALSE)
  }
  if (!is.null(over) && (!is.character(over) || !length(over))) {
    stop("'over' must be a character vector of coefficient names, or NULL.",
         call. = FALSE)
  }
  StartSearch(label = paste0("search with ", optimizer@name),
              optimizer = optimizer, over = over)
}


#' @title Print a Starting-Value Strategy
#' @name print.start_strategy
#' @description
#' Prints a one-line description of a starting-value strategy: its label, and
#' for [start_random()] whether the draw is centered, for [start_search()]
#' which optimizer it searches with.
#' @param x A [start_strategy()] or any object inheriting from it.
#' @param ... Unused.
#' @return `x`, invisibly, as a print method should.
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
#' @param strategy A `StartIntercepts` object.
#' @param spec,design,obj,... As in [start_at()].
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
#' The position of the intercept of one equation's parametric block, or
#' `NA` where it has none.
#'
#' @details
#' What is looked for is a column of the parametric block, which is a
#' stricter test than a coefficient whose name ends in `(Intercept)`.
#' [modelterms7::nl()] names the intercept of each of its own parameters the
#' same way, and those live on those parameters' own charts, not on the
#' predictor's, so a model written `y ~ 0 + nl(...)` puts one of them
#' first. Writing the
#' intercept-only fit there sets the parameter to `linkinv` of a value
#' that was never on its scale: measured on a logistic growth curve whose
#' asymptote rides a log link, `mean(y) = 23.9` became a starting
#' \eqn{\phi} of `2.5e10`, an objective of `7.0e20` and a gradient
#' of `1.4e21`, on data whose every scale is ordinary. The lasso path
#' built at those coefficients then spanned `2.8e15` to `2.8e19`
#' where the block empties at about 300, so all of its points were the same
#' empty fit and every subject deviation was estimated as exactly zero.
#'
#' @param spec The specification.
#' @param design The design.
#' @param p The distribution parameter naming the equation.
#'
#' @return An integer position into the equation's coefficient vector, or
#'   `NA_integer_`.
#'
#' @seealso [start_at()], [statmod_intercepts()]
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
#' The response carried onto the scale of one equation's linear predictor.
#' [modelterms7::term_coef_start()] takes it to estimate a term's own
#' parameters from the data, which is the one thing a term cannot work out
#' for itself.
#'
#' @details
#' It exists only where the response reads the parameter directly, which
#' `params_interpretation` says: a mean or a location. A scale or a shape has
#' no per-observation reading, so `NULL` is returned and nothing is invented
#' for the occasion.
#'
#' The scale matters and is not a detail. Measured on a Poisson whose
#' predictor is a logistic growth curve with \eqn{\phi = 4}, a term handed
#' the raw response estimated \eqn{\phi} between 52.7 and 54.8 over five
#' samples, and one handed \eqn{g(y)} estimated it between 4.03 and 4.07.
#' What does not matter, measured on the same shape beside another term, is
#' the other terms' contribution: subtracting it moved the estimate from
#' 39.68 to 39.87 against a truth of 40, so the response is passed as it
#' stands and nothing is residualized.
#'
#' A response sitting on a bound of the parameter's own domain is moved half
#' way to the nearest admissible value, a link not being defined at its
#' bound. For counts under a log link that is the classical half.
#'
#' @param spec The specification.
#' @param p The distribution parameter naming the equation.
#'
#' @return A numeric vector, one value per observation, or `NULL`.
#'
#' @seealso [start_at()],
#'   [modelterms7::term_coef_start()]
#' @keywords internal
predictor_target <- function(spec, p) {
  d <- spec@distrib
  if (!("params_interpretation" %in% S7::prop_names(d))) return(NULL)
  ip <- d@params_interpretation[[p]]
  if (is.null(ip) || !ip %in% c("mean", "location")) return(NULL)
  # A RESPONSE THAT IS NOT A VECTOR OF NUMBERS has no target: a censored one
  # is an S7 object carrying values and their statuses, and the values alone
  # are not observations of the parameter -- a right-censored one is a lower
  # bound. `as.numeric` RAISES on such an object rather than warning, so the
  # suppression here caught nothing and a censored fit stopped before its
  # first iteration.
  y <- tryCatch(suppressWarnings(as.numeric(spec@response)),
                error = function(e) numeric(0))
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
#' @param strategy A `StartOrigin` object.
#' @param spec,design,obj,... As in [start_at()].
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
#' @param strategy A `StartRandom` object.
#' @param spec,design,obj,... As in [start_at()].
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
#' The predicate for the second is [refreshes_own_block()], the same one
#' [unfittable_reason()] uses, so a term written later is covered with no
#' edit here. A convex block is left out deliberately: the scoring step
#' reaches its optimum from anywhere, and a search over it would spend the
#' budget where it buys nothing.
#'
#' A penalized coordinate is then removed wherever it sits, which is not the
#' same question as whether its term is convex: a penalty declared through a
#' sub-formula of a break-point or nonlinear term covers columns of that
#' term's block, and the loop takes the whole block.
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
  idx <- sort(unique(idx))
  # A PENALIZED COORDINATE IS LEFT OUT WHEREVER IT SITS, and asking whether
  # its TERM is convex answers for the wrong thing: a random effect declared
  # through a sub-formula of a break-point term is a column of that term's
  # block, so the loop above takes it. Measured on a jseg over thirty groups
  # whose seven own coefficients each carry `random(~1|id)`, that was 210 of
  # 219 coordinates, and with a differenced gradient one iteration cost 2.3 s
  # against the milliseconds of a scoring step. The budget is the smaller
  # half of it: the search runs on the LIKELIHOOD ALONE, where a deviation
  # is identified by nothing -- the penalty is what identifies it, exactly as
  # for a filter's deviations -- so the search fits each group's own points
  # and hands the fit a start further from the penalized mode than the one it
  # began with. A structural term's penalty covers positions among the term's
  # own parameters and indexes no column, so it removes nothing here.
  pen <- unlist(lapply(statmod_penalized(spec, design), function(u) {
    if (isTRUE(u$structural)) integer(0) else u$index
  }), use.names = FALSE)
  setdiff(idx, as.integer(pen))
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
#' the likelihood, exactly as `statmod_fit_joint()` does, and leaves the
#' best it found in place --- which is how a strategy returning coefficients
#' can nonetheless start a filter somewhere better than its own
#' `term_start()`.
#' @param strategy A `StartSearch` object.
#' @param spec,design,obj,... As in [start_at()].
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
