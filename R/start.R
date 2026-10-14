#' @include statmod.R
NULL

#' @title S7 Class for a Starting-Value Strategy
#'
#' @description
#' The abstract parent of the starting-value strategies. A strategy defines
#' where a fit begins as a procedure and not as a vector of numbers: fit the
#' intercept-only model, draw around it, search the likelihood for a basin,
#' or reuse the estimates of another fit. `statmod(start =)` accepts a
#' strategy, or a named list of coefficient vectors.
#'
#' @details
#' [statmod()] calls a strategy exactly once, before the alternation between
#' the coefficients and the hyperparameters begins. An optimizer runs at every
#' step of that alternation, which is why a global search belongs here: given
#' to `inner_optimizer` it would rerun at every hyperparameter value that the
#' outer search tries and reach the same point each time.
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
#'   [start_random()], [start_search()] and [start_from()] for the five
#'   shipped strategies,
#'   [start_at()] for the generic they implement.
#'
#' @examples
#' # Every shipped strategy inherits from this.
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
#' Returns the [start_strategy()] class object from the namespace at the time
#' of the call, so that a membership test reads the class as it exists then.
#'
#' @details
#' A class captured in a variable at build time makes such a test fail under
#' \pkg{covr}, which re-evaluates a package's code and so re-creates the
#' class. `S7::S7_inherits()` against the captured copy is then
#' `FALSE` for an object of the live class. `distributions7` and
#' `linkfunctions7` document the same behavior of `identical()` on a class
#' that has been re-created.
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
#' strategy. The shipped methods show the range: [start_origin()] reads
#' only the design's widths, [start_intercepts()] fits a small model,
#' [start_random()] draws from the caller's generator, and [start_search()]
#' runs an optimizer over `obj`.
#'
#' A method must return a full-length vector for every distribution
#' parameter, including the ones it does not set. Zero is the conventional
#' filler, and it is the value at which the penalty of a penalized block is
#' smallest.
#'
#' @param strategy A [start_strategy()], which decides the method.
#' @param spec The [StatmodSpec()] being fitted.
#' @param design The design, as [statmod_design()] returns it.
#' @param obj The objective, as [statmod_objective()] returns it, or `NULL`.
#'   Only [start_search()] reads it, and [start_search()] signals an error
#'   when `obj` is `NULL`. The other four shipped strategies run with `NULL`;
#'   [start_from()] passes `obj` on to its `rest` strategy.
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
#' # a strategy that needs no objective is called with NULL for it; the ones
#' # that search the likelihood receive the objective from statmod()
#' start_at(start_origin(), spec, design, NULL)
#' start_at(start_intercepts(), spec, design, NULL)
#'
#' @export
start_at <- S7::new_generic("start_at", "strategy",
  function(strategy, spec, design, obj, ...) S7::S7_dispatch())


#' @title S7 Classes for the Shipped Strategies
#' @description
#' The five classes that [start_intercepts()], [start_origin()],
#' [start_random()], [start_search()] and [start_from()] instantiate. Each
#' inherits from [start_strategy()] and carries its `label`; `StartRandom`
#' adds `fn`, `args` and `center`, `StartSearch` adds `optimizer` and `over`,
#' and `StartFrom` adds `fit` and `rest`. Use the constructors, which
#' validate; these raw ones do not.
#' @param label A single string, the name that the strategy prints.
#' @param fn The generator a random start draws from, called as `fn(n, ...)`.
#'   `StartRandom` only.
#' @param args A list of further arguments to `fn`. `StartRandom` only.
#' @param center A single logical: whether the draw is added to the
#'   intercept-only start. `StartRandom` only.
#' @param optimizer,over An \pkg{optimizers7} optimizer and an optional
#'   character vector of coefficient labels. `StartSearch` only.
#' @param fit,rest A [StatmodFit()] and a [start_strategy()]. `StartFrom`
#'   only.
#' @return An S7 object inheriting from [start_strategy()]:
#'   `StartIntercepts` and `StartOrigin` with the one inherited property
#'   `label`, `StartRandom` with `fn`, `args` and `center` besides,
#'   `StartSearch` with `optimizer` and `over`, and `StartFrom` with `fit` and
#'   `rest`.
#' @seealso [start_intercepts()], [start_origin()], [start_random()],
#'   [start_search()], [start_from()] for the constructors, [start_at()] for
#'   the methods registered on these classes.
#' @name StartIntercepts-class
#' @aliases StartIntercepts StartOrigin StartRandom StartSearch StartFrom
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
#' of the model with every covariate removed. Where the first parameter is a
#' location or a mean on the identity link, the slopes of its parametric
#' columns start at a least-squares fit, and a term with a starting value of
#' its own starts at the value that the term declares. The remaining
#' coefficients start at zero.
#'
#' @details
#' The model with no covariates is the model with every slope set to zero, so
#' its maximum likelihood estimate is a natural starting point for the
#' intercepts of the model with covariates. It requires one small fit, on an
#' intercept-only design.
#'
#' Three further rules apply:
#'
#' - A location or mean on the identity link whose equation has an intercept
#'   and at least one other parametric column starts its parametric
#'   coefficients at a weighted least-squares fit net of the offset, and the
#'   intercepts of the other parameters are fitted to the residuals. This
#'   applies when every term of that equation is a parametric, smooth, random
#'   or penalized term.
#' - A term whose [modelterms7::term_coef_start()] method returns values that
#'   are not all zeros, such as a break-point or a nonlinear term, starts at
#'   those values.
#' - All other coefficients start at zero, which is where the penalty of a
#'   penalized block is smallest.
#'
#' Two adjustments apply to the intercepts:
#'
#' - The mean of an **offset** is subtracted from the intercept before the
#'   value is written, because the intercept-only fit does not see the offset
#'   and the predictor would otherwise start displaced by it.
#' - The value is written to a **parametric intercept** and to nothing else.
#'   A term such as [modelterms7::nl()] also names the intercept of each of
#'   its own parameters `(Intercept)`, but those coefficients lie on the
#'   chart of that parameter and not on the scale of the predictor.
#'
#' This is what [statmod()] does when `start` is `NULL`. The constructor
#' exists so that the default can be named, compared with other strategies and
#' passed on.
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
#' # The intercept starts at the response's own scale, and the slope at its
#' # least-squares value.
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
#' # Use of the origin
#'
#' A location parameter on the identity link begins at zero whatever the
#' response is. On a response centered at a thousand the run then starts a
#' thousand units from the response and has to travel there, and on a badly
#' scaled model it may not arrive. [start_intercepts()], the default, avoids
#' this by starting the intercepts at the intercept-only fit.
#'
#' The origin is appropriate where every equation is on a log or a logit
#' link, since zero is then an ordinary interior value of the parameter.
#'
#' # The name
#'
#' The function is named `start_origin()` and not `start_zeros()`, because
#' \pkg{optimizers7} already exports `start_zeros()` for the starting point of
#' an optimizer, and no two members of this toolkit export the same name. A
#' collision would be reported by [statmodels7_conflicts()], and the function
#' reached by a caller would depend on the order in which the packages were
#' attached.
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
#' Draws every coefficient from a generator chosen by the caller and, by
#' default, adds the draw to the intercept-only fit. It is intended for
#' running a model from several starts and comparing where they land, which
#' is one way to explore a multimodal likelihood.
#'
#' @details
#' # Centering
#'
#' A coefficient drawn from a standard normal is a reasonable perturbation but
#' not a reasonable absolute value, because the intercept of a location
#' equation is on the scale of the response, which may be a thousand. Adding
#' the draw to the answer of [start_intercepts()] keeps that scale and
#' randomizes the position around it, so that the several starts differ from
#' one another.
#'
#' With `center = FALSE` the draw is used alone, which is appropriate when the
#' generator is already scaled to the problem.
#'
#' # The random stream
#'
#' The draw comes from the session's generator, so [set.seed()] governs the
#' result, and a fit begun this way repeats only with the same seed. The
#' intercept-only fit to which a centered draw is added sets a seed of its own
#' and restores the caller's stream afterwards, so it does not change the
#' draw.
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

#' @rdname StartIntercepts-class
#' @keywords internal
StartFrom <- S7::new_class("StartFrom", parent = start_strategy,
  properties = list(fit = S7::class_any, rest = S7::class_any))


#' @title Start From Another Fit's Estimates
#'
#' @description
#' Takes the coefficients that were estimated in a fitted model, matches them
#' to the model that is about to be fitted, and starts the fit there. The
#' coefficients that the two models do not share are set by a second
#' strategy.
#'
#' @details
#' Two models of the same response often share most of their coefficients: one
#' drops a covariate, one adds a term, one changes the family. The estimates
#' of the first model are then a closer starting point for the second than a
#' rule that uses only the response.
#'
#' # What is matched, and what is projected
#'
#' The two designs are compared block by block, a block being the columns of
#' one term, and a block is treated in one of three ways:
#' \describe{
#'   \item{the parametric block}{is matched column by column, by name. The
#'     names of a model matrix are variables, levels and interactions, so a
#'     covariate that is added or dropped leaves the other columns matched.}
#'   \item{a block written on the same basis}{is carried across as it stands,
#'     because its coefficient names are the same coordinates.}
#'   \item{a block written on another basis}{is estimated by least squares
#'     against the part of the reference's predictor that the blocks already
#'     carried across do not explain.}
#' }
#'
#' The third rule makes two bases comparable. The coefficients of a smooth are
#' coordinates in a basis that the fit rotates, so
#' `s(x, bspline_smooth(k = 6))` and `s(x, bspline_smooth(k = 10))` carry the
#' names `s(x).z1` to `s(x).z4` in common but give them different meanings,
#' and one name does not identify one function. What the two bases share is
#' the fitted function, so the coefficients are found by projecting it,
#' \deqn{\hat\beta_S = \arg\min_\beta \lVert X_S\beta - r\rVert^2, \qquad
#' r = \eta^{\mathrm{ref}} - \eta_{-S},}
#' over the pending columns \eqn{S} alone, where \eqn{\eta^{\mathrm{ref}}} is
#' the reference's predictor and \eqn{\eta_{-S}} is the new model's predictor
#' with the pending columns at zero. Both predictors are the actual predictors
#' of their models, including the value of a block that moves with its
#' coefficients (a break-point or a nonlinear term). The equality of the two
#' blocks, \eqn{B_{\mathrm{new}}\beta_{\mathrm{new}} =
#' B_{\mathrm{old}}\beta_{\mathrm{old}}}, has a solution only if the old span
#' lies inside the new one, and the projection gives the closest function that
#' the new basis can represent. Where the spans are nested the projection is
#' exact, and where the two bases coincide it returns the coefficients it was
#' given, so the second rule above is a shortcut and not a different answer.
#'
#' Two blocks are paired by the stem of their coefficient names and by the
#' class of the term, and not by the key of the block, which is the deparsed
#' call of the term and therefore differs as soon as `k` does.
#'
#' A distribution parameter that the reference does not have, such as `theta`
#' when a negative binomial is fitted from a Poisson, is set by `rest`, so the
#' families need not agree. So is any block for which the reference has no
#' counterpart. A block that matches by name is carried across whatever the
#' response is. The projection of a block written on another basis is skipped,
#' and the block is set by `rest`, when the two responses differ or when
#' either model carries a structural term. A pending block that itself moves
#' with its coefficients (for example `seg(x)` refitted as
#' `seg(x, npsi = 2)`) is not projected either, because its coefficients are
#' not those of a linear block, and it keeps the values given by `rest`.
#'
#' The starting values are most useful where the same model is refitted many
#' times, as with a coefficient held at a sequence of values in an interval by
#' inversion, or with a model built up one term at a time.
#'
#' @param fit A [StatmodFit()] to take the estimates from.
#' @param rest The strategy for whatever the two models do not share.
#'   [start_intercepts()] by default, which is what a fit given no strategy
#'   uses. It receives the objective, so [start_search()] can serve as `rest`.
#'
#' @return A `StartFrom` object, inheriting from [start_strategy()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(200), z = runif(200))
#' dd$y <- rpois(200, exp(0.4 + 0.6 * dd$x - 0.3 * dd$z))
#' full <- statmod(y ~ x + z, distributions7::poisson_distrib(), dd)
#'
#' # the same model with one covariate dropped, started where the other
#' # one ended
#' sub <- statmod(y ~ z, distributions7::poisson_distrib(), dd,
#'                start = start_from(full))
#' coef(sub)
#'
#' # what was taken, and from where
#' spec <- statmod_spec(y ~ z, distributions7::poisson_distrib(), dd)
#' attr(start_at(start_from(full), spec, statmod_design(spec), NULL), "taken")
#'
#' @seealso [start_intercepts()], the default; [start_at()], the generic.
#' @export
start_from <- function(fit, rest = start_intercepts()) {
  if (!S7::S7_inherits(fit, StatmodFit)) {
    stop("'fit' must be a statmod fit to take the estimates from.",
         call. = FALSE)
  }
  if (!S7::S7_inherits(rest, start_strategy_class())) {
    stop("'rest' must be a start strategy, e.g. start_intercepts().",
         call. = FALSE)
  }
  if (S7::S7_inherits(rest, StartFrom)) {
    # a chain of references would have to be resolved in some order and the
    # order would decide the answer; one reference and one fallback is the
    # whole contract
    stop("'rest' cannot itself be start_from(): name one reference fit.",
         call. = FALSE)
  }
  StartFrom(label = "another fit's estimates", fit = fit, rest = rest)
}


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
#' A search is a starting value and is not part of the scoring step. If it
#' were given to `inner_optimizer`, it would be rerun at every hyperparameter
#' value that the outer criterion tries and, inside [cv()], at every point of
#' the path for every fold, and each run would reach the same point.
#'
#' # On the likelihood alone
#'
#' The penalties are off. A starting value has to select the basin of the
#' likelihood, which is the one thing that the fit does not correct by itself.
#' The penalties enter afterwards, when their hyperparameters are estimated,
#' and at the probe values of the hyperparameters they are not a chosen
#' penalty, so a search with them would explore the likelihood of a different
#' model.
#'
#' # Which coordinates are searched
#'
#' The non-convex ones, which the toolkit can name: the parameters of a
#' structural term ([modelterms7::gas()], [modelterms7::regime()]), the
#' coefficients of a term that recomputes its own block
#' ([modelterms7::nl()], [modelterms7::seg()] and the break-point terms), and
#' each equation's intercept.
#'
#' All other coefficients keep their default values. A smooth, a ridge or a
#' random effect is a convex block whose optimum the scoring step reaches from
#' any start, and a search over a thousand random-effect coefficients would
#' spend the budget on the part of the model that needs a search least.
#'
#' A penalized coordinate is excluded **wherever it sits**, including inside
#' a non-convex term. A sub-formula can develop a break-point or a nonlinear
#' parameter over groups, and the group deviations are columns of the block of
#' that term. On the likelihood alone nothing identifies these deviations,
#' because the penalty identifies them, so a search over them would fit the
#' points of each group and move away from the penalized mode. In a `jseg()`
#' term with `by = ~random(~1|id)`, most coordinates are such deviations.
#'
#' The argument `over` replaces this choice by naming the coefficients.
#' A structural term's own parameters are searched in either case.
#'
#' # The hyperparameters
#'
#' The search does not cover the hyperparameters. Its objective is the
#' likelihood with the penalties off, in which no hyperparameter appears, so
#' no proposal can change one.
#'
#' A global search over the hyperparameters is a search over the **outer**
#' criterion, where each proposal requires a full refit and not one likelihood
#' evaluation. It is requested by `statmod(outer_optimizer =
#' optimizers7::sa())`. The hyperparameters of a kinked penalty are outside
#' that search as well: they are chosen along a warm-started path scored by
#' `sparse_criterion`, which a random jump would not improve.
#'
#' @param optimizer The optimizer to search with, any \pkg{optimizers7}
#'   optimizer. Defaults to `optimizers7::sa()`, simulated annealing, which
#'   is the one built for a multimodal surface.
#' @param over Optional character vector naming the coefficients to search
#'   over, overriding the choice above. A coefficient is named by its label,
#'   the parameter and the coefficient joined by a colon (`"sigma:(Intercept)"`);
#'   a bare name is accepted where only one equation carries it. Names not in
#'   the model, and a bare name that several equations carry, are an error.
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
#' @description
#' The answer of the default strategy: the intercept-only fit, the
#' least-squares start of a location where [location_ols_start()] applies, and
#' the start that each term declares through
#' [modelterms7::term_coef_start()]. See [start_intercepts()].
#' @param strategy A `StartIntercepts` object.
#' @param spec,design,obj,... As in [start_at()].
#' @return A named list of numeric vectors.
#' @keywords internal
S7::method(start_at, StartIntercepts) <-
  function(strategy, spec, design, obj, ...) {
    params <- spec@distrib@params
    out <- stats::setNames(lapply(design, function(d) numeric(d$npar)),
                           params)
    # A LOCATION ON THE IDENTITY LINK starts from a least-squares fit of its
    # parametric columns, and the other parameters from the residuals: the
    # intercept-only fit reads the response's marginal shape, which a strong
    # covariate distorts. Measured on 17 location families, 4 slopes and 3
    # seeds against a fit started at the truth: the skew normal in its direct
    # chart stopped 19.6 and 8.2 log-likelihood units below its maximum in two
    # of 12 regressions from the intercept-only start and in none from this
    # one, and the 204 fits took 85 s against 198 s.
    ols <- location_ols_start(spec, design)
    eta0 <- if (is.null(ols)) statmod_intercepts(spec) else ols$eta0
    for (p in params) {
      if (design[[p]]$npar == 0L) next
      v <- eta0[[p]]
      if (is.null(v) || !is.finite(v)) next
      if (!is.null(ols) && identical(p, ols$param)) {
        # the least-squares coefficients were fitted net of the offset
        out[[p]][ols$cols] <- ols$coef
        next
      }
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

#' @title A Least-Squares Start for a Location on the Identity Link
#' @name location_ols_start
#'
#' @description
#' Fits the parametric columns of the first parameter's equation to the
#' response by least squares (weighted by the prior weights, net of that
#' equation's offset), and fits the distribution without covariates to the
#' residuals shifted by the fitted intercept. Returns the coefficients of
#' those columns and the intercepts of every parameter, or `NULL` where the
#' start does not apply.
#'
#' @details
#' It applies when the family is univariate, its first parameter is a
#' location or a mean on the identity link, that parameter's equation has an
#' intercept and at least one other parametric column, and every term of the
#' equation is a [modelterms7::LinparTerm], a [modelterms7::SmoothTerm], a
#' [modelterms7::RandomTerm] or a [modelterms7::PenalizedTerm]. The penalized
#' blocks start at zero, as they do from the intercept-only start. If the
#' equation has a term of any other kind (a break-point, a nonlinear term, a
#' filter), the least-squares start does not apply and the equation keeps the
#' intercept-only start. An aliased column takes the coefficient zero.
#'
#' @param spec The specification.
#' @param design The design.
#'
#' @return `NULL`, or a list with `param` (the first parameter's name),
#'   `cols` (positions in its coefficient vector), `coef` (their starting
#'   values) and `eta0` (link-scale intercepts, as from
#'   [statmod_intercepts()]).
#'
#' @seealso [start_at()], [statmod_intercepts()]
#' @keywords internal
location_ols_start <- function(spec, design) {
  d <- spec@distrib
  if (S7::S7_inherits(d, distributions7::multivariate_distrib)) return(NULL)
  y <- spec@response
  if (!is.numeric(y) || !is.null(dim(y))) return(NULL)
  p <- d@params[[1L]]
  if (!identical(d@link_params[[p]]@link_name, "identity")) return(NULL)
  if (!tolower(d@params_interpretation[[1L]]) %in% c("location", "mean")) {
    return(NULL)
  }
  if (is.null(design[[p]]) || design[[p]]$npar == 0L) return(NULL)
  terms <- spec@terms[[p]]
  plain <- vapply(terms, function(tm)
    S7::S7_inherits(tm, modelterms7::LinparTerm) ||
      S7::S7_inherits(tm, modelterms7::SmoothTerm) ||
      S7::S7_inherits(tm, modelterms7::RandomTerm) ||
      S7::S7_inherits(tm, modelterms7::PenalizedTerm), logical(1))
  if (!length(terms) || !all(plain)) return(NULL)
  ii <- parametric_intercept(spec, design, p)
  if (is.na(ii)) return(NULL)
  cols <- sort(unique(unlist(lapply(names(terms), function(nm)
    if (S7::S7_inherits(terms[[nm]], modelterms7::LinparTerm))
      design[[p]]$blocks[[nm]]), use.names = FALSE)))
  if (length(cols) < 2L) return(NULL)
  X <- as.matrix(design[[p]]$X[, cols, drop = FALSE])
  off <- spec@offsets[[p]]
  ya <- if (is.null(off)) y else y - off
  w <- spec@weights
  ok <- is.finite(ya) & apply(is.finite(X), 1L, all)
  if (length(w) == length(y)) ok <- ok & is.finite(w) & w > 0
  if (sum(ok) <= length(cols)) return(NULL)
  fit <- tryCatch(if (length(w) == length(y))
    stats::lm.wfit(X[ok, , drop = FALSE], ya[ok], w[ok]) else
      stats::lm.fit(X[ok, , drop = FALSE], ya[ok]),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  b <- fit$coefficients
  b[!is.finite(b)] <- 0
  jj <- match(ii, cols)
  r <- as.numeric(ya - X %*% b)
  z <- (r + b[[jj]])[ok]
  spec_r <- tryCatch(S7::set_props(spec, response = z), error = function(e) NULL)
  if (is.null(spec_r)) return(NULL)
  eta0 <- statmod_intercepts(spec_r)
  if (is.null(eta0[[p]]) || !is.finite(eta0[[p]])) return(NULL)
  # the location of the residual distribution is the intercept
  b[[jj]] <- eta0[[p]]
  list(param = p, cols = cols, coef = as.numeric(b), eta0 = eta0)
}

#' @title Where an Equation's Intercept Is
#' @name parametric_intercept
#'
#' @description
#' The position of the intercept of one equation's parametric block, or
#' `NA` where it has none.
#'
#' @details
#' The position searched for is a column of the parametric block, which is a
#' stricter test than a coefficient whose name ends in `(Intercept)`.
#' [modelterms7::nl()] names the intercept of each of its own parameters the
#' same way, and those coefficients lie on the charts of the parameters of
#' the term and not on the scale of the predictor, so a model written
#' `y ~ 0 + nl(...)` has one of them first. Writing the intercept-only fit
#' there would set the parameter to the inverse link of a value that was never
#' on its scale, and for a parameter on a log link the starting value would be
#' far outside the range of the data.
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
#' A target exists only where the response reads the parameter directly,
#' which `params_interpretation` records as a mean or a location. A scale or a
#' shape has no per-observation reading, so `NULL` is returned.
#'
#' The response is passed on the scale of the predictor, \eqn{g(y)}, because a
#' term that estimates its parameters from the raw response would work on the
#' wrong scale. The contribution of the other terms of the equation is not
#' subtracted, and the offset of the equation is.
#'
#' A response that lies on a bound of the domain of the parameter is moved
#' half way from the bound to the smallest (or, for an upper bound, the
#' largest) response inside the domain, because a link is not defined at its
#' bound. For counts under a log link this gives one half when the smallest
#' positive count is 1.
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
#' @description The draw, added to the intercept-only fit unless `center` is
#'   `FALSE`.
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
#' The predicate for the second is [refreshes_own_block()], which tests
#' whether the term registers a [modelterms7::term_refresh()] method of its
#' own, so a term written later is covered with no edit here. A convex block
#' is left out deliberately: the scoring step reaches its optimum from any
#' start, and a search over it would use the budget where it gains nothing.
#'
#' A penalized coordinate is then removed wherever it sits, which differs from
#' testing whether its term is convex: a penalty declared through a
#' sub-formula of a break-point or nonlinear term covers columns of the block
#' of that term, and the loop takes the whole block.
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
#'
#' The objective `obj` is required, and an error is signaled when it is
#' `NULL`. A coefficient named in `over` is resolved by its label, with a bare
#' name accepted where one equation alone carries it.
#' @param strategy A `StartSearch` object.
#' @param spec,design,obj,... As in [start_at()].
#' @return A named list of numeric vectors.
#' @keywords internal
S7::method(start_at, StartSearch) <-
  function(strategy, spec, design, obj, ...) {
    if (is.null(obj)) {
      stop(paste0("start_search() evaluates the objective, which statmod()",
                  " passes as 'obj'."), call. = FALSE)
    }
    base <- start_at(start_intercepts(), spec, design, obj)
    beta0 <- obj$stack(base)

    idx <- if (is.null(strategy@over)) {
      search_coords(spec, design)
    } else {
      params <- spec@distrib@params
      nms <- unlist(lapply(params, function(p) design[[p]]$coef_names),
                    use.names = FALSE)
      eq <- rep(params, vapply(params, function(p)
        length(design[[p]]$coef_names), integer(1)))
      labs <- paste(eq, nms, sep = ":")
      # a coefficient is named by its label, parameter:name, as coef_labels()
      # writes it; a bare name is accepted where one equation alone carries it
      hit <- vapply(strategy@over, function(o) {
        k <- match(o, labs)
        if (!is.na(k)) return(k)
        k <- which(nms == o)
        if (length(k) > 1L) {
          stop(sprintf(paste0("'over' names '%s', which is a coefficient of",
                              " several equations.\n  Name one by its label:",
                              " %s."), o, paste(labs[k], collapse = ", ")),
               call. = FALSE)
        }
        if (!length(k)) NA_integer_ else k
      }, integer(1))
      if (anyNA(hit)) {
        stop(sprintf("'over' names %s, which is not a coefficient.",
                     paste(strategy@over[is.na(hit)], collapse = ", ")),
             call. = FALSE)
      }
      sort(unique(unname(hit)))
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


#' @title Starting Values From Another Fit
#' @name start_at.StartFrom
#' @description
#' The reference fit's estimates where the two models share a coefficient,
#' the fallback strategy's answer everywhere else.
#' @details
#' A block is matched by the name that the formula gave its term: column by
#' column for the parametric block, and as a whole for every other block where
#' the coefficient names agree. The remaining blocks, which are the same term
#' written on another basis and paired by [block_stem()] and
#' [same_term_kind()], go to [project_blocks()], which estimates them against
#' the partial residual of the reference's predictor. See [start_from()] for
#' the three rules.
#'
#' The result carries an attribute `"taken"`, a data frame that names every
#' coefficient taken from the reference and records whether it was `matched`
#' or `projected`, so that the reused coefficients can be read directly.
#' @param strategy A `StartFrom` object.
#' @param spec,design,obj,... As in [start_at()].
#' @return A named list of numeric vectors, with the attribute `"taken"`.
#' @keywords internal
S7::method(start_at, StartFrom) <-
  function(strategy, spec, design, obj, ...) {
    out <- start_at(strategy@rest, spec, design, obj)
    ref <- strategy@fit
    rdesign <- statmod_design(ref@spec)
    took <- list()
    pending <- list()

    for (p in spec@distrib@params) {
      # a parameter the reference does not carry is left to the fallback,
      # which is what lets a negative binomial start from a Poisson
      if (!p %in% ref@spec@distrib@params) next
      if (design[[p]]$npar == 0L) next
      rb <- rdesign[[p]]$blocks
      rn <- rdesign[[p]]$coef_names
      rstem <- vapply(names(rb), function(z) block_stem(rn[rb[[z]]]),
                      character(1))
      for (nm in names(design[[p]]$blocks)) {
        idx <- design[[p]]$blocks[[nm]]
        if (!length(idx)) next
        here <- design[[p]]$coef_names[idx]
        if (!nm %in% names(rb)) {
          # A BLOCK IS KEYED BY THE TERM'S OWN CALL, so `s(x, bspline_smooth(k = 6))` and
          # `s(x, bspline_smooth(k = 10))` are two keys and never meet by name. What they
          # share is the stem of their coefficient names and the term's
          # class, which is what pairs them for the projection.
          if (!is.na(block_stem(here))) {
            cand <- names(rb)[!is.na(rstem) & rstem == block_stem(here)]
            if (length(cand) == 1L &&
                same_term_kind(spec@terms[[p]][[nm]],
                               ref@spec@terms[[p]][[cand]])) {
              pending[[p]] <- c(pending[[p]], list(list(nm = nm, idx = idx)))
            }
          }
          next
        }
        ridx <- rb[[nm]]
        if (!length(ridx)) next
        there <- rn[ridx]
        take <- if (is_parametric_block(spec, p, nm)) {
          # BY COLUMN: adding or dropping a covariate is the case this
          # exists for, and a model matrix names variables and levels
          k <- match(here, there)
          data.frame(to = idx[!is.na(k)], from = ridx[k[!is.na(k)]],
                     name = here[!is.na(k)], stringsAsFactors = FALSE)
        } else if (identical(here, there)) {
          # THE SAME BASIS: the coefficients are the same coordinates, so
          # they are carried across as they stand, which is exact and costs
          # nothing. Re-estimating them would be work for no answer.
          data.frame(to = idx, from = ridx, name = here,
                     stringsAsFactors = FALSE)
        } else {
          # THE SAME TERM ON ANOTHER BASIS: a basis the fit rotates gives the
          # same names to different coordinates, so the coefficients are not
          # comparable and the FUNCTION is. The block is left to the
          # projection below rather than to the fallback.
          pending[[p]] <- c(pending[[p]], list(list(nm = nm, idx = idx)))
          NULL
        }
        if (is.null(take) || !nrow(take)) next
        out[[p]][take$to] <- ref@coefficients[[p]][take$from]
        took[[length(took) + 1L]] <-
          data.frame(parameter = p, term = nm, coefficient = take$name,
                     value = ref@coefficients[[p]][take$from],
                     how = "matched", stringsAsFactors = FALSE)
      }
    }
    pr <- project_blocks(spec, design, ref, rdesign, out, pending)
    out <- pr$start
    took <- c(took, pr$taken)

    taken <- if (length(took)) do.call(rbind, took) else
      data.frame(parameter = character(0), term = character(0),
                 coefficient = character(0), value = numeric(0),
                 how = character(0), stringsAsFactors = FALSE)
    rownames(taken) <- NULL
    structure(out, taken = taken)
  }


#' Project a Reference Fit's Predictor onto the Blocks That Did Not Match
#'
#' @description
#' Estimates the coefficients of the blocks that the `StartFrom` method of
#' [start_at()] left pending, which are terms that the two models share
#' written on a different basis. The estimate is by least squares against the
#' part of the reference's predictor that the blocks already carried across do
#' not explain.
#'
#' @details
#' Two bases of the same term span different subspaces, so
#' \eqn{B_{\mathrm{new}}\beta_{\mathrm{new}} = B_{\mathrm{old}}
#' \beta_{\mathrm{old}}} has no solution in general, and the projection
#' \deqn{\hat\beta_S = \arg\min_\beta \lVert X_S\beta - r\rVert^2, \qquad
#' r = \eta^{\mathrm{ref}} - \eta_{-S},}
#' gives the closest reproduction of the fitted function that the new basis
#' allows. Only the pending columns \eqn{S} are estimated. The vector \eqn{r}
#' is the partial residual of the reference's predictor
#' \eqn{\eta^{\mathrm{ref}}} after the blocks already matched by name are
#' removed at the values they were given (\eqn{\eta_{-S}} is the new model's
#' predictor with the pending columns at zero). Both predictors are read
#' through [statmod_design_at()], so a block that moves with its coefficients
#' (a break-point or a nonlinear term) contributes its value and not its
#' linearization \eqn{X\beta}. Where every block matched exactly there is
#' nothing to solve, and a refit whose bases all agree returns the
#' coefficients it was given.
#'
#' The projection is exact when the reference's function lies in the span of
#' the new basis, and it reduces to the identity when the two bases coincide,
#' so the name match above is a shortcut and not a different answer. Whatever
#' the reference explains and the new model has no column for, such as a
#' covariate that the formula dropped, stays in \eqn{r} and is absorbed by
#' the pending blocks.
#'
#' Three conditions restrict the projection. The two responses must be
#' identical, so that the two designs are read at the same rows: the
#' projection compares two predictors observation by observation and has no
#' meaning across different data. Neither model may carry a structural term,
#' whose contribution is the state of a recursion and not a function of the
#' coefficients. Where either of these two conditions fails, the coefficients
#' are returned unchanged. A pending block that itself moves with its
#' coefficients (for example `seg(x)` refitted as `seg(x, npsi = 2)`) is not
#' projected, because its coefficients are not those of a linear block, and it
#' keeps the values given by the fallback strategy. The pending columns of an
#' equation are also left to the fallback when their design has more than
#' \eqn{5 \times 10^7} entries.
#'
#' @param spec The [StatmodSpec()] being started.
#' @param design Its design.
#' @param ref The reference [StatmodFit()].
#' @param rdesign The reference's design.
#' @param out The starting coefficients as they stand, with the matched
#'   blocks already written in.
#' @param pending A list, one entry per distribution parameter, of
#'   `list(nm, idx)` pairs naming a block and the columns it occupies.
#'
#' @return A list with `start`, the coefficients with the pending blocks
#'   estimated, and `taken`, a list of one-row-per-coefficient data frames
#'   recording what was projected.
#'
#' @seealso [start_from()], whose method calls it.
#'
#' @keywords internal
project_blocks <- function(spec, design, ref, rdesign, out, pending) {
  taken <- list()
  if (!length(pending)) return(list(start = out, taken = taken))
  # the two predictors are compared observation by observation, so they have
  # to be read at the same rows; an identical response is what says so, and
  # is what start_from()'s own page asks the caller for
  if (!identical(spec@response, ref@spec@response)) {
    return(list(start = out, taken = taken))
  }
  if (length(attr(design, "structural")) ||
      length(attr(rdesign, "structural"))) {
    return(list(start = out, taken = taken))
  }
  # the predictor of one equation at given coefficients: X b, plus the
  # adjustment that statmod_design_at() records for a block that moves with
  # its coefficients
  eta_at <- function(sp, des, cf, p) {
    d <- statmod_design_at(sp, cf, des)[[p]]
    e <- as.numeric(d$X %*% cf[[p]])
    if (!is.null(d$adj)) e <- e + d$adj
    e
  }
  for (p in names(pending)) {
    # a block that moves with its own coefficients is not linear in them
    pending[[p]] <- Filter(function(b)
      !refreshes_own_block(spec@terms[[p]][[b$nm]]), pending[[p]])
    if (!length(pending[[p]])) next
    cols <- sort(unique(unlist(lapply(pending[[p]], function(b) b$idx))))
    A <- design[[p]]$X[, cols, drop = FALSE]
    # a smooth's block is dense by construction, so this densifies nothing
    # that was sparse for a reason; the guard is for a factor `by`, where the
    # block is one copy per level and the projection is not worth a
    # gigabyte
    if (prod(dim(A)) > 5e7) next
    held <- out
    held[[p]][cols] <- 0
    r <- eta_at(ref@spec, rdesign, ref@coefficients, p) -
      eta_at(spec, design, held, p)
    b <- tryCatch(qr.coef(qr(as.matrix(A)), r), error = function(e) NULL)
    if (is.null(b)) next
    # a column the pivot drops has no least-squares value and keeps whatever
    # the fallback gave it, which is the same rule vcov() follows for an
    # aliased coordinate
    keep <- is.finite(b)
    if (!any(keep)) next
    out[[p]][cols[keep]] <- b[keep]
    for (bl in pending[[p]]) {
      ok <- bl$idx %in% cols[keep]
      if (!any(ok)) next
      taken[[length(taken) + 1L]] <-
        data.frame(parameter = p, term = bl$nm,
                   coefficient = design[[p]]$coef_names[bl$idx[ok]],
                   value = out[[p]][bl$idx[ok]],
                   how = "projected", stringsAsFactors = FALSE)
    }
  }
  list(start = out, taken = taken)
}


#' The Stem Shared by a Block's Coefficient Names
#'
#' @description
#' Returns the part every coefficient name of a block has in common before
#' its own suffix, which is what identifies the term across two models
#' written on different bases.
#'
#' @details
#' A block is keyed in the design by the term's deparsed call, so
#' `s(x, bspline_smooth(k = 6))` and `s(x, bspline_smooth(k = 10))` are two different keys and never meet by
#' name. Their coefficients, on the other hand, are `s(x).lin`, `s(x).z1`,
#' ... on both sides, so the stem `s(x)` names the term without naming the
#' basis. The suffix is taken from the last dot and not the first, so a
#' covariate whose own name carries one -- `s(my.var)` -- keeps it.
#'
#' `NA` is returned where the names do not reduce to one common stem, as for
#' the parametric block of several columns, whose columns are variables and
#' levels and which is matched column by column instead. A block of a single
#' column whose name contains a dot, such as `my.var`, reduces to the stem
#' `my`.
#'
#' @param nms A block's coefficient names.
#'
#' @return A single string, or `NA_character_` where the block has no stem.
#'
#' @seealso [start_from()], whose method pairs blocks with it.
#'
#' @keywords internal
block_stem <- function(nms) {
  if (!length(nms)) return(NA_character_)
  st <- unique(sub("\\.[^.]*$", "", nms))
  if (length(st) != 1L || !nzchar(st[[1L]]) || identical(st[[1L]], nms[[1L]]))
    return(NA_character_)
  st[[1L]]
}


#' Whether Two Terms Are of the Same Kind
#'
#' @description
#' `TRUE` when two built terms have the same class, which is what keeps the
#' stem of [block_stem()] from pairing a smooth with something else that
#' happens to name its coefficients the same way.
#'
#' @param a,b Two built terms, either of which may be `NULL`.
#'
#' @return A single logical.
#'
#' @seealso [block_stem()], the other half of the pairing.
#'
#' @keywords internal
same_term_kind <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  identical(class(a), class(b))
}


#' Whether a Term Is the Equation's Parametric Block
#'
#' @description
#' `TRUE` for the block of unpenalized parametric columns, which is the one
#' [start_from()] matches column by column.
#'
#' @details
#' The test reads the class of the term and not the name of the term. A caller
#' may give a term any name, and the deciding property is that the block
#' carries ordinary model-matrix columns whose names are variables and levels,
#' so that a coefficient of one model has the same meaning as the coefficient
#' of the same name in another. Every other kind of block is a basis, a set of
#' levels or a standardized design, where a shared name does not imply a
#' shared meaning.
#'
#' @param spec A [StatmodSpec()].
#' @param param The distribution parameter.
#' @param nm The term's name in that parameter's equation.
#'
#' @return A single logical.
#'
#' @seealso [start_from()], which uses it to decide how a block is matched.
#'
#' @keywords internal
is_parametric_block <- function(spec, param, nm) {
  tm <- spec@terms[[param]][[nm]]
  if (is.null(tm)) return(FALSE)
  S7::S7_inherits(tm, modelterms7::LinparTerm)
}
