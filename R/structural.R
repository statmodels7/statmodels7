#' @include refresh.R
NULL

# A structural term rewrites the predictor instead of contributing columns to
# it. What the fitting layer has to do differently is three things: leave the
# term out of the design, add the level its filter produces to the equation it
# sits in, and correct the gradient of every equation for the fact that the
# level was driven by scores read at predictors those equations also enter.
# The last is the reverse recursion of modelterms7::term_adjoint(): the
# derivative of the objective in a coefficient is not the block times the
# score, and treating it as though it were would converge to the wrong point
# while reporting a gradient of zero.

#' The Structural Terms of a Specification
#'
#' @description
#' Returns the parameter, name and shape of every structural term in a
#' specification.
#'
#' @details
#' The shape is `"filter"` for a term that produces a predictor, which
#' the layer adds to the equation in which the term sits, and `"loglik"` for
#' a term whose contribution is a likelihood mixed over components (the
#' states of a regime chain, the intervals or nodes of a marginal break-point
#' term), which cannot be written as a predictor.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A list of entries with `param`, `term` and `kind`,
#'   possibly empty.
#'
#' @seealso [statmod_filter_at()]
#'
#' @keywords internal
statmod_structural <- function(spec) {
  out <- list()
  for (p in names(spec@terms)) {
    for (nm in names(spec@terms[[p]])) {
      tm <- spec@terms[[p]][[nm]]
      if (!S7::S7_inherits(tm, modelterms7::structural_term)) next
      out[[length(out) + 1L]] <- list(param = p, term = nm,
                                      kind = structural_kind(tm))
    }
  }
  out
}


#' Which Shape a Structural Term Implements
#'
#' @description
#' Returns `"filter"` when the term registers its own
#' [modelterms7::term_filter()] method, `"loglik"` when it registers its own
#' [modelterms7::term_loglik()] method, and `""` when it registers neither.
#'
#' @details
#' The shape is read from the methods that the term registers, not from a
#' list of class names, so a structural term written outside the package is
#' routed without an edit here. The class on which a method was registered is
#' `attr(m, "signature")[[1]]`, and a method registered on `structural_term`
#' itself does not count. The filter shape is tested first.
#'
#' @param term A term.
#'
#' @return A single string.
#'
#' @keywords internal
structural_kind <- function(term) {
  owns <- function(generic) {
    m <- tryCatch(S7::method(generic, S7::S7_class(term)),
                  error = function(e) NULL)
    if (is.null(m)) return(FALSE)
    cls <- attr(m, "signature")[[1L]]
    !identical(attr(cls, "name"), "structural_term")
  }
  if (owns(modelterms7::term_filter)) return("filter")
  if (owns(modelterms7::term_loglik)) return("loglik")
  ""
}


#' The State of the Structural Terms
#'
#' @description
#' Returns the environment in which a design carries the own parameters of
#' the structural terms, on the unconstrained scale of
#' [modelterms7::term_links()].
#'
#' @param design The design.
#'
#' @return An environment, or `NULL`.
#'
#' @keywords internal
statmod_structural_state <- function(design) attr(design, "structure")


#' Reuse a Structural Quantity Computed at the Same Point
#'
#' @description
#' A depth-one exact memo on the design's structural state: where the last
#' call's key is `identical()` to this one's, the stored value is
#' returned; otherwise `compute()` runs and replaces it.
#'
#' @details
#' The criterion, its gradient and the joint step's curvature all read the
#' same mode, so the same coefficients and the same term parameters recur
#' within a fit. The
#' cache returns the previously computed object itself, so a hit is
#' bit-identical to recomputing by construction, and the key is compared
#' with `identical()` on the full numeric inputs, so a collision
#' cannot happen.
#'
#' The memo is not used where the design carries refreshable terms: a
#' break-point block advances its rescaling schedule as the alternation
#' commits, so the same coefficients do not imply the same design there,
#' and a key that does not include the schedule cannot justify returning a
#' stored value.
#'
#' @param design The design whose structural state holds the slots.
#' @param slot A short name, one cache per slot.
#' @param key The exact inputs on which the value depends.
#' @param compute A function of no arguments.
#'
#' @return `compute()`'s value, possibly from the cache.
#'
#' @keywords internal
structural_memo <- function(design, slot, key, compute) {
  sst <- statmod_structural_state(design)
  if (is.null(sst) || !is.null(attr(design, "refresh"))) return(compute())
  k_slot <- paste0(slot, "_key")
  v_slot <- paste0(slot, "_value")
  k <- if (exists(k_slot, envir = sst, inherits = FALSE)) {
    get(k_slot, envir = sst)
  } else NULL
  if (!is.null(k) && identical(k, key)) return(get(v_slot, envir = sst))
  val <- compute()
  assign(k_slot, key, envir = sst)
  assign(v_slot, val, envir = sst)
  val
}


#' The Starting Parameters of a Structural Term
#'
#' @description
#' Returns what the term itself declares through
#' [modelterms7::term_start()]: the parameters as near the model without the
#' term as its charts allow, which only the term can determine. A score
#' loading on the log chart has no coordinate for zero, so it starts at a weak
#' response instead.
#'
#' @param term A built structural term.
#' @param target The response on the scale of the term's equation, as
#'   [predictor_target()] returns it, or `NULL`. A term whose start is read
#'   off the data uses it: the marginal break-point term (its exact profile)
#'   and the regime term (quantiles of the target) do so, and the other
#'   methods ignore it through the dots.
#'
#' @return A named numeric vector.
#'
#' @keywords internal
structural_zeta_start <- function(term, target = NULL) {
  modelterms7::term_start(term, target = target)
}


#' The Starting Points of a Structural Term
#'
#' @description
#' Returns the starts given by [modelterms7::term_starts()] for the
#' structural term of the model, when that function gives more than one and
#' the specification carries no values from a previous fit.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design, as [statmod_design()] returns it.
#'
#' @return `NULL` where there is nothing to try, otherwise a list with `term`,
#'   the unit's key, and `zeta`, the list of starts, whose first element is
#'   the one that [statmod_design()] already used.
#'
#' @keywords internal
structural_multistarts <- function(spec, design) {
  su <- attr(design, "structural")
  for (u in su) {
    if (!is.null(spec@structural[[u$term]])) next
    z <- modelterms7::term_starts(spec@terms[[u$param]][[u$term]],
                                  target = predictor_target(spec, u$param))
    if (length(z) > 1L) return(list(term = u$term, zeta = z))
  }
  NULL
}


#' From the Unconstrained Scale to the Term's Parameters
#'
#' @param term A built structural term.
#' @param zeta The unconstrained values.
#'
#' @return A named list on the parameter scale.
#'
#' @keywords internal
structural_psi <- function(term, zeta) {
  links <- modelterms7::term_links(term)
  nm <- modelterms7::term_params(term)
  stats::setNames(lapply(nm, function(j)
    linkfunctions7::linkinv(links[[j]], zeta[[j]])), nm)
}


#' The Score and Curvature That Drive a Filter
#'
#' @description
#' Returns the derivative of the log-density in one distribution parameter's
#' predictor, and its second derivative, as functions of that predictor at
#' one observation, while every other parameter stays where it is.
#'
#' @details
#' They are read one observation at a time because a score-driven filter
#' evaluates them at the predictor that it has just produced, and that
#' predictor is not known before the recursion reaches the observation. A
#' term whose callbacks do not depend on the state is given the whole index
#' vector instead. A regime chain is such a term: its levels shift a
#' predictor known in advance.
#'
#' @param spec A [StatmodSpec()].
#' @param theta The per-observation parameters, as
#'   [statmod_eta()] returns them.
#' @param p The distribution parameter in whose equation the term sits.
#' @param scaling The exponent \eqn{d} of [modelterms7::gas()]: where it is
#'   not zero, `score` returns \eqn{u = s\,\mathcal{I}^{-d}} and `curvature`
#'   its derivative in the predictor,
#'   \eqn{u' = s'\,\mathcal{I}^{-d} - d\,s\,\mathcal{I}^{-d-1}\mathcal{I}'},
#'   and the compiled context carries `scaling`, and its kernel composes the
#'   same quantities from the family's scalar entry points.
#'
#' @return A list with `score`, `curvature`, `logdens` and `fast`, where
#'   `fast` is the context for the compiled filter, or `NULL` where the
#'   registries of \pkg{distributions7} and \pkg{linkfunctions7} do not cover
#'   the family and the link.
#'
#' @keywords internal
structural_callbacks <- function(spec, theta, p, scaling = 0) {
  d <- spec@distrib
  y <- spec@response
  lk <- d@link_params[[p]]
  n <- spec@n_obs
  params <- d@params
  # A filter calls these once per TIME STEP, so anything that does not
  # depend on the step belongs out here: the recycling of each parameter to
  # the sample length, and the component key, which is a function of the
  # parameter names alone and was being rebuilt with match() and paste() on
  # every call. Measured, this is worth about 12 per cent of a gas fit and
  # no more -- the cost of a step is dominated by the two generic calls
  # below, not by their arguments. It is NOT a change of order: rep_len()
  # on a vector that is already the right length returns it untouched, so
  # the old form was linear in n like this one, which the timings confirm.
  theta_n <- lapply(params, function(q) rep_len(theta[[q]], n))
  names(theta_n) <- params

  # distributions7 resolves the family's methods and the link's once and
  # applies the chain rule for the ONE component wanted. Through the
  # generics each of these calls validated its arguments, aligned theta by
  # name, dispatched, and assembled every component of its order to keep
  # one of them -- which a filter pays once per observation per iteration,
  # and which a profile of a fitted model charged for the whole of its time.
  k <- distributions7::distrib_kernel(d, p)

  # The list of scalars the kernel takes was rebuilt by an lapply over every
  # parameter on every call, and the score and the curvature at one step each
  # built their own -- so a filter of n steps allocated 2np lists per
  # iteration. What varies with the step is only a parameter that varies BY
  # OBSERVATION, which is a property of the fit and not of the step, so the
  # list is built once and only those entries are written; a parameter carried
  # by an intercept alone is written never. The buffer is reused rather than
  # copied, which is what removes the allocation, so a caller must read it
  # before asking for the next step -- the kernel does, and nothing else holds
  # it.
  #
  # Measured on a panel of 25 groups and 750 observations, 18.08 s to 16.64 s
  # and the log-likelihood identical to ten decimals: 1.09x, which is worth
  # keeping and is not where the time is. A profile charges 41.6 per cent of
  # the total to this closure and that number is misleading -- almost all of
  # it is the two kernel calls underneath, which the allocation merely sits
  # in front of. What the fit really spends is filter runs, and the number of
  # those is decided by the fitting strategy rather than by this function.
  varying <- which(vapply(params, function(q) length(theta[[q]]) > 1L,
                          logical(1)))
  buf <- lapply(theta_n, function(v) v[1L])
  names(buf) <- params
  # A regime term is handed the WHOLE index vector at once, its callbacks not
  # depending on the state, and then the buffer is the wrong shape and the
  # comparison below is a vector: that route keeps the general form, which is
  # what it was already paying for k calls rather than 2nk.
  last <- -1L
  general <- function(i) lapply(theta_n, function(v) v[i])
  at <- if (!length(varying)) {
    function(i) if (length(i) == 1L) buf else general(i)
  } else {
    function(i) {
      if (length(i) != 1L) return(general(i))
      if (i != last) {
        for (q in varying) buf[[q]] <<- theta_n[[q]][i]
        last <<- i
      }
      buf
    }
  }

  # the fast context of piano_parallel.txt section 2a: where the C
  # registries of distributions7 and linkfunctions7 cover this family
  # and this link, the filter's kernel reads the score and the
  # curvature through their scalar entry points instead of calling back
  # into R -- the same composition, bit-identical, held to that by
  # modelterms7's twin test -- and an uncovered pair leaves the context
  # inert, the callbacks below running as before. distrib_scalar_route()
  # gives the name the registry knows the distribution by and the constants
  # its entries read after the parameters (the binomial's size), and NULL
  # where the registry does not cover it; link_scalar_route() does the same
  # for the link, with the link's own parameters
  route <- if (!is.matrix(y)) distributions7::distrib_scalar_route(d)
  lroute <- linkfunctions7::link_scalar_route(lk)
  fast <- if (!is.null(route) && !is.null(lroute)) {
    consts <- lapply(route$constants, function(v) {
      v <- as.numeric(v)
      if (length(v) == 1L) v else rep_len(v, n)
    })
    list(family = route$name,
         link = lroute$name,
         link_par = lroute$par,
         k = match(p, params),
         bounds = as.numeric(lk@link_bounds),
         y = as.numeric(y),
         theta = c(theta_n, consts))
  } else NULL

  if (scaling != 0) {
    check_filter_scaling(d, y, theta_n)
    if (!is.null(fast)) fast$scaling <- scaling
    return(list(
      score = function(e, i) {
        th <- at(i)
        as.numeric(k$score(y[i], th, e)) *
          as.numeric(k$information(y[i], th, e))^(-scaling)
      },
      curvature = function(e, i) {
        th <- at(i)
        s <- as.numeric(k$score(y[i], th, e))
        info <- as.numeric(k$information(y[i], th, e))
        as.numeric(k$curvature(y[i], th, e)) * info^(-scaling) -
          scaling * s * info^(-scaling - 1) *
          as.numeric(k$dinformation(y[i], th, e))
      },
      logdens = function(e, i) as.numeric(k$logdens(y[i], at(i), e)),
      fast = fast))
  }

  list(
    score = function(e, i) as.numeric(k$score(y[i], at(i), e)),
    curvature = function(e, i) as.numeric(k$curvature(y[i], at(i), e)),
    logdens = function(e, i) as.numeric(k$logdens(y[i], at(i), e)),
    fast = fast
  )
}


#' Run the Terms of the Likelihood Shape at the Current Parameters
#'
#' @description
#' Computes the smoothed posterior over the components of each structural
#' term of the likelihood shape (the states of a regime term, the intervals
#' or nodes of a marginal break-point term), the shifts that the components
#' add to the predictor, and the posterior-weighted predictor that these
#' give.
#'
#' @details
#' A term of this shape does not report a predictor: its contribution is a
#' likelihood mixed over components, so the log-likelihood of the model comes
#' from [modelterms7::term_loglik()] and not from the density at a single
#' predictor.
#'
#' The remaining quantities follow from Fisher's identity: the derivative of
#' that mixed likelihood in any predictor of the model is the
#' posterior-weighted derivative of the ordinary one. A caller therefore
#' differentiates its own log-density once per component, vectorized, and
#' weights.
#'
#' The predictor reported for the equation in which the term sits is the
#' posterior-weighted one, and a fitted value here is that predictor.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param eta_static The static predictors.
#' @param theta_static The parameters that they give.
#'
#' @return A named list, one entry per structural term of the likelihood
#'   shape (a regime term or a marginal break-point term), each with `param`,
#'   `term`, `tm`, `psi`, `gamma` (the posterior), `mu` (the shifts), `cb`
#'   (the callbacks), `eta_static` and `eta_bar` (the posterior-weighted
#'   predictor).
#'
#' @seealso [statmod_structural()]
#'
#' @keywords internal
statmod_regime_at <- function(spec, design, eta_static, theta_static) {
  su <- attr(design, "structural")
  st <- statmod_structural_state(design)
  out <- list()
  for (u in su) {
    if (!identical(u$kind, "loglik")) next
    tm <- spec@terms[[u$param]][[u$term]]
    psi <- structural_psi(tm, st$zeta[[u$term]])
    cb <- structural_callbacks(spec, theta_static, u$param)
    # the shifts each mixture component adds to the predictor are the term's
    # own answer, so a regime chain and a marginal break-point ride one
    # route; a component whose shift varies by observation -- a quadrature
    # node's hinge value -- comes back as a matrix, aligned with the
    # posterior's columns, and the callbacks let the term rebuild its nodes
    mu <- modelterms7::term_levels(tm, psi, eta = eta_static[[u$param]],
                                   y = spec@response, logdens = cb$logdens)
    gam <- modelterms7::term_posterior(tm, eta_static[[u$param]],
                                       spec@response, cb$logdens, psi)
    eta_bar <- if (is.matrix(mu)) {
      rowSums(gam * mu) + eta_static[[u$param]]
    } else {
      as.numeric(gam %*% mu) + eta_static[[u$param]]
    }
    out[[u$term]] <- list(param = u$param, term = u$term, tm = tm, psi = psi,
                          gamma = gam, mu = mu, cb = cb,
                          eta_static = eta_static[[u$param]],
                          eta_bar = eta_bar)
  }
  out
}


#' The Shift of One Mixture Component
#'
#' @description
#' Returns the level that a component adds to the predictor: one number for
#' a regime or a side pattern, one value per observation for a quadrature
#' node.
#'
#' @param mu What [modelterms7::term_levels()] returned.
#' @param k The component.
#'
#' @return A single number or a vector over the observations.
#'
#' @keywords internal
component_shift <- function(mu, k) {
  if (is.matrix(mu)) mu[, k] else mu[[k]]
}


#' How Many Mixture Components a Term Reports
#'
#' @param mu What [modelterms7::term_levels()] returned.
#'
#' @return A single integer.
#'
#' @keywords internal
component_count <- function(mu) {
  if (is.matrix(mu)) ncol(mu) else length(mu)
}


#' The Parameters Under One Regime
#'
#' @description
#' Returns the distribution's parameters at every observation, with the
#' equation of the term shifted by the level of one component.
#'
#' @param spec A [StatmodSpec()].
#' @param eta_static The static predictors.
#' @param param The equation in which the term sits.
#' @param shift The level to add.
#'
#' @return A named list of per-observation parameters.
#'
#' @keywords internal
statmod_theta_shifted <- function(spec, eta_static, param, shift) {
  links <- spec@distrib@link_params
  params <- spec@distrib@params
  stats::setNames(lapply(params, function(q) {
    e <- eta_static[[q]]
    if (identical(q, param)) e <- e + shift
    linkfunctions7::linkinv(links[[q]], e)
  }), params)
}


#' The Name of a Third-Derivative Component
#'
#' @description
#' Locates the \eqn{(a, b, c)} entry of a distribution's third-derivative
#' list, which is keyed by [distributions7::deriv_names()].
#'
#' @details
#' The name is built from the parameter names in the family's own order and
#' is not parsed out of an existing name, so that a parameter whose own name
#' contains the separator is handled as \pkg{distributions7} handles it.
#'
#' @param params The parameter names, in the family's order.
#' @param a,b,c Indices into `params`.
#'
#' @return A single string.
#'
#' @keywords internal
deriv3_key <- function(params, a, b, c) {
  paste(params[sort(c(a, b, c))], collapse = "_")
}


#' The Name of a Fourth-Derivative Component
#'
#' @description
#' Locates the \eqn{(a, b, c, d)} entry of a distribution's fourth-derivative
#' list, built the same way [deriv3_key()] builds its own.
#'
#' @details
#' The fourth order is needed where a filter's third derivative is: each
#' order of differentiating the predictor through the recursion requires one
#' more order of the family, because the score that drives the recursion is
#' read at the predictor that the recursion produces. The values come from
#' [distributions7::distrib_deriv4()].
#'
#' @param params The parameter names, in the family's order.
#' @param a,b,c,d Indices into `params`.
#'
#' @return A single string.
#'
#' @keywords internal
deriv4_key <- function(params, a, b, c, d) {
  paste(params[sort(c(a, b, c, d))], collapse = "_")
}


#' The Name of a Fifth-Derivative Component
#'
#' @description
#' Locates the \eqn{(a, b, c, d, e)} entry of a distribution's
#' fifth-derivative list, built the same way [deriv4_key()] builds its
#' own.
#'
#' @details
#' The fifth order is needed where a filter's fourth derivative is, by the
#' rule that the page of [deriv4_key()] states one order down: the criterion's
#' second derivative in a pair of hyperparameters reads
#' [modelterms7::term_fourth()], and the score that drives that recursion is
#' read at the predictor that the recursion produces.
#'
#' [distributions7::distrib_deriv5()] supplies the values in closed form for
#' the families that register their own fifth-order method, and as one
#' central difference of the fourth order for the others.
#'
#' @param params The parameter names, in the family's order.
#' @param a,b,c,d,e Indices into `params`.
#'
#' @return A single string.
#'
#' @keywords internal
deriv5_key <- function(params, a, b, c, d, e) {
  paste(params[sort(c(a, b, c, d, e))], collapse = "_")
}


#' The Observed Information Over the Coefficients and a Filter's Parameters
#'
#' @description
#' Computes the information of a model carrying a structural filter, over the
#' coefficients of every equation and the term's own parameters, with the
#' recursion's own curvature in it.
#'
#' @details
#' The gradient of such a model is exact already, through
#' [modelterms7::term_adjoint()]. The matrix that the scoring step inverts is
#' assembled as though the level were an offset; it serves as a scoring
#' matrix and differs from the observed information, which this function
#' computes. Writing
#' \eqn{u} for the coefficients followed by the term's parameters on the
#' unconstrained scale, \eqn{V_{q,t}} for the derivative of equation
#' \eqn{q}'s predictor in \eqn{u}, and \eqn{E_t} for the second derivative
#' of the predictor that the filter produces,
#'
#' \deqn{-\frac{\partial^2 \ell}{\partial u \partial u^\top}
#'   = -\sum_t w_t \sum_{q,r} \ell_{qr,t} V_{q,t}^\top V_{r,t}
#'     - \sum_t w_t \ell_{p,t} E_t.}
#'
#' Only the equation carrying the filter has a \eqn{V} that is not its own
#' design: there it is the forward Jacobian of the recursion, which
#' [modelterms7::term_curvature()] returns beside the contracted
#' \eqn{E}. The third derivatives that the second sum needs are those of
#' \pkg{distributions7}, from [distributions7::distrib_deriv3()].
#'
#' The term's parameters are the last columns, which is the convention that
#' `term_curvature()` shares with its caller.
#'
#' A term that mixes over latent states, in place of shifting the predictor,
#' is routed to [statmod_regime_information()], whose Hessian
#' comes from the same forward recursion that computes the likelihood; for a
#' model that carries neither kind of term, the ordinary information is used.
#'
#' @param spec A [StatmodSpec()].
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#'
#' @return A square matrix over the coefficients followed by the term's free
#'   parameters (a level held because an intercept carries it leaves the
#'   matrix), or the ordinary observed information where there is no
#'   structural term.
#'
#' @seealso [statmod_structural_score()],
#'   [statmod_regime_information()]
#'
#' @keywords internal
statmod_full_information <- function(spec, coef,
                                     design = statmod_design(spec)) {
  params <- spec@distrib@params
  ev <- statmod_eta(spec, design, coef)
  if (!length(ev$filters) && !length(ev$regimes)) {
    return(statmod_information_at(spec, coef, design, expected = FALSE))
  }
  # the criterion, its gradient and the joint step all read this matrix at
  # the mode they share, so the same point is asked for several times over
  sst <- statmod_structural_state(design)
  structural_memo(design, "info",
                  list(coef = coef, zeta = sst$zeta), function() {
    statmod_full_information_impl(spec, coef, design, params, ev)
  })
}

#' @rdname statmod_full_information
#' @param params,ev The parameter names and the evaluated predictors,
#'   already in hand at the caller.
#' @keywords internal
statmod_full_information_impl <- function(spec, coef, design, params, ev) {
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  nb <- sum(npar)
  n <- spec@n_obs
  w <- spec@weights

  if (length(ev$regimes)) {
    return(statmod_regime_information(spec, ev, design, npar, offs, nb, n, w))
  }

  f <- ev$filters[[1L]]
  ap <- match(f$param, params)
  np <- length(modelterms7::term_params(f$tm))
  m <- nb + np

  # one row per observation, each equation's design placed in its own
  # columns and zero elsewhere, the term's block included.
  #
  # as_dense() on the block, and the reason is not cosmetic: an equation
  # carrying a random effect has a SPARSE design, and writing a dgCMatrix
  # into a slice of a base matrix is a length error rather than a
  # conversion -- `y ~ x + random(~1|id) + gas(...)` stopped here with
  # "number of items to replace is not a multiple of replacement length".
  # The densification is LOCAL and buys back nothing sparsity was keeping:
  # this matrix is one column per unknown with a single block filled, it is
  # read one ROW at a time by the callback below (the access a
  # compressed-column matrix is worst at), and the fit, the penalties and
  # the solve all go on seeing the sparse design.
  Vs <- lapply(seq_along(params), function(a) {
    M <- matrix(0, n, m)
    if (npar[a] > 0L) {
      M[, offs[a] + seq_len(npar[a])] <- as_dense(design[[params[a]]]$X)
    }
    M
  })

  gl <- distributions7::distrib_gradient(spec@distrib, spec@response,
                                         ev$theta, scale = "link", threads = spec@threads)
  H <- distributions7::distrib_hessian(spec@distrib, spec@response, ev$theta,
                                       scale = "link", threads = spec@threads)
  D3 <- distributions7::distrib_deriv3(spec@distrib, spec@response, ev$theta,
                                       scale = "link", threads = spec@threads)

  # The derivatives are read at the FITTED predictor, which the recursion
  # inside term_curvature() reproduces exactly: it runs at the same
  # parameters over the same static predictor, so the two agree by
  # construction and the callback need not evaluate the family again.
  #
  if (isTRUE(f$joint)) {
    cv <- joint_filter_curvature(spec, design, f, Vs, gl, H, D3, ev$theta)
    for (k in seq_along(f$params)) {
      Vs[[match(f$params[k], params)]] <- cv$jacobian[[k]]
    }
  } else {
    cv <- filter_curvature(spec, design, f, ap, Vs, gl, H, D3, ev$theta)
    Vs[[ap]] <- cv$jacobian
  }

  out <- matrix(0, m, m)
  for (a in seq_along(params)) {
    for (b in seq_along(params)) {
      wv <- w * rep_len(H[[hess_key(params, a, b)]], n)
      out <- out + wcrossprod(Vs[[a]], wv, Vs[[b]], spec@threads)
    }
  }
  # the (a, b) and (b, a) terms are transposes of one another and are formed
  # by separate crossproducts, so the sum is symmetric to a rounding and not
  # to the bit; a caller about to factor it wants the symmetry exact
  out <- -(out + cv$curvature)
  out <- (out + t(out)) / 2
  # a level the equation's intercept carries is held rather than estimated,
  # so its row and column leave the information: keeping them would report a
  # matrix that is singular along exactly that direction
  held <- statmod_structural_state(design)$held[[f$term]]
  if (length(held)) {
    zn <- modelterms7::term_params(f$tm)
    drop <- nb + match(held, zn)
    out <- out[-drop, -drop, drop = FALSE]
  }
  out
}


#' A Filter's Second-Order Recursion at the Fitted Point
#'
#' @description
#' Runs [modelterms7::term_curvature()] at the parameters at which the filter
#' was fitted. The family's derivatives are looked up instead of evaluated
#' again, and the function returns the forward Jacobian of the filtered
#' predictor beside the contracted second derivative.
#'
#' @details
#' Three callers read it at the same point: [statmod_full_information()], for
#' the joint information; [filter_joint_movement()], for the rows along which
#' the exact outer gradient reads how the mode's movement reaches the
#' determinant; and [filter_joint_jacobian()], for the standard error of a
#' filtered predictor. All go through the one memo slot, so the recursion runs
#' once at a point.
#'
#' The recursion is run again here at parameters at which it has already been
#' run, so the predictor that it reaches is the one at which the derivatives
#' were read, and the callbacks read the stored values instead of evaluating
#' the family again. In the filter itself the score is evaluated at a
#' predictor that the recursion has just produced, which cannot be known in
#' advance. The seed is part of the memo's key, so a caller whose layout
#' differed would miss the cache and not take the wrong shape.
#'
#' The assembly of the callback is written once, in `.structural_blocks()`:
#' the exact gradient of a penalized filter needs the same pieces at one order
#' higher. The same pieces are also passed as data, with which an eligible
#' term runs its second-order recursion compiled, and the callback is kept for
#' the cases that the kernel does not handle.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param f The filter's entry of [statmod_eta()]'s `filters`.
#' @param ap Which distribution parameter the filter sits in.
#' @param Vs The static rows, one matrix per distribution parameter over the
#'   coefficients followed by the term's parameters.
#' @param gl,H,D3 The family's first three derivatives on the link scale at
#'   the fitted predictors.
#' @param theta The per-observation parameters at which they were read, from
#'   which the expected information is read where the filter scales its score.
#'
#' @return The list that [modelterms7::term_curvature()] returns.
#'
#' @keywords internal
filter_curvature <- function(spec, design, f, ap, Vs, gl, H, D3, theta) {
  params <- spec@distrib@params
  n <- spec@n_obs
  w <- spec@weights
  # the recursion reads the derivatives of what drives it, the contraction
  # weighs by the log-likelihood's own score; the two coincide at scaling 0
  dr <- filter_driving(spec, theta, ap, filter_scaling(f$tm), gl, H, D3)
  blocks <- .structural_blocks(params, ap, Vs, dr$H, dr$D3, NULL, n)(NULL)
  s_at <- rep_len(gl[[f$param]], n)
  sd_at <- rep_len(dr$gl[[f$param]], n)
  cd_at <- rep_len(dr$H[[hess_key(params, ap, ap)]], n)
  bd_data <- structural_blocks_data(params, ap, Vs, dr$H, dr$D3, n)
  structural_memo(design, "curv",
                  list(zeta = f$psi, eta = f$eta_static,
                       g = w * s_at, seed = Vs[[ap]]), function() {
    modelterms7::term_curvature(
      f$tm, f$eta_static, spec@response,
      function(e, i) sd_at[i], function(e, i) cd_at[i], f$psi,
      w * s_at, Vs[[ap]], blocks,
      score_values = sd_at, curvature_values = cd_at,
      blocks_data = bd_data, threads = spec@threads)
  })
}


#' The Observed Information of a Model Carrying a Mixture Term
#'
#' @description
#' Computes the observed information over the coefficients of every equation
#' together with the own parameters of a structural term of the likelihood
#' shape (a regime term or a marginal break-point term), from the exact
#' Hessian of the mixed log-likelihood.
#'
#' @details
#' The matrix that [statmod_information_at()] assembles for a mixture is
#' the complete-data information, the ordinary one averaged over the
#' smoothed states. It serves as a scoring matrix and differs from the
#' observed information: by Louis's missing-information principle the two
#' differ by the conditional variance of the complete-data score, which is
#' positive semidefinite, so the complete-data matrix is the larger and a
#' standard error read off it is too small.
#'
#' [modelterms7::term_hessian()] returns the exact Hessian by
#' propagating first and second derivatives through the same scaled forward
#' recursion that computes the likelihood. This function supplies the part of
#' the calculation that belongs to the model: how each equation's unknowns
#' reach each predictor, and the family's first and second derivatives at the
#' predictor to which each component shifts.
#'
#' Unlike the filter's, these derivatives cannot be looked up from a single
#' evaluation: a component shifts the predictor by a level of its own, so the
#' family is evaluated once per component, vectorized over observations.
#'
#' @param spec A [StatmodSpec()].
#' @param ev The evaluated predictors, from [statmod_eta()].
#' @param design The design.
#' @param npar,offs,nb The per-equation coefficient counts, their offsets and
#'   their total.
#' @param n,w The number of observations and their weights.
#'
#' @return A square matrix over the stacked coefficients followed by the
#'   term's own free parameters.
#'
#' @keywords internal
statmod_regime_information <- function(spec, ev, design, npar, offs, nb, n,
                                       w) {
  params <- spec@distrib@params
  links <- spec@distrib@link_params
  r <- ev$regimes[[1L]]
  ap <- match(r$param, params)
  zn <- modelterms7::term_params(r$tm)
  np <- length(zn)
  m <- nb + np

  Vs <- lapply(seq_along(params), function(a) {
    M <- matrix(0, n, m)
    # a random effect's block is a sparse Matrix beside a structural term
    # (densify_small() leaves it so), and it does not assign into a base
    # matrix: summary() of a marginal jump beside random(~1 | id) stopped here
    if (npar[a] > 0L) {
      M[, offs[a] + seq_len(npar[a])] <- as.matrix(design[[params[a]]]$X)
    }
    M
  })

  # the parameters at a predictor the term has shifted, which is where the
  # family has to be evaluated again rather than looked up
  theta_at <- function(e, i) {
    stats::setNames(lapply(seq_along(params), function(q) {
      x <- if (q == ap) e else rep_len(ev$eta_static[[params[q]]], n)[i]
      linkfunctions7::linkinv(links[[params[q]]], x)
    }), params)
  }
  yof <- function(i) {
    if (is.matrix(spec@response)) spec@response[i, , drop = FALSE]
    else spec@response[i]
  }

  ld <- function(e, i) {
    distributions7::distrib_pdf(spec@distrib, yof(i), theta_at(e, i),
                                log = TRUE, threads = spec@threads)
  }
  gr <- function(e, i) {
    g <- distributions7::distrib_gradient(spec@distrib, yof(i),
                                          theta_at(e, i), scale = "link", threads = spec@threads)
    vapply(params, function(q) rep_len(g[[q]], length(i)), numeric(length(i)))
  }
  he <- function(e, i) {
    H <- distributions7::distrib_hessian(spec@distrib, yof(i),
                                         theta_at(e, i), scale = "link", threads = spec@threads)
    out <- array(0, c(length(i), length(params), length(params)))
    for (a in seq_along(params)) {
      for (b in seq_along(params)) {
        out[, a, b] <- rep_len(H[[hess_key(params, a, b)]], length(i))
      }
    }
    out
  }

  oh <- modelterms7::term_hessian(r$tm, r$eta_static, spec@response,
                                  logdens = ld, grad = gr, hess = he,
                                  psi = r$psi, seed = Vs,
                                  cols = nb + seq_len(np), level = ap,
                                  weights = w)
  out <- -oh$hessian
  out <- (out + t(out)) / 2
  held <- statmod_structural_state(design)$held[[r$term]]
  if (length(held)) {
    drop <- nb + match(held, zn)
    out <- out[-drop, -drop, drop = FALSE]
  }
  out
}


#' The Structural Terms' Estimated Parameters
#'
#' @description
#' Returns one named vector per structural term, on the scale that the
#' term's own documentation uses for its parameters, together with the
#' unconstrained values that were estimated.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design the fit ran on.
#'
#' @return A named list, one entry per structural term, each a list with
#'   `parameter` (the named vector on the parameter scale), `unconstrained`,
#'   `equation` (the distribution parameter in which the term sits) and `held`
#'   (the names of the parameters that an intercept of the equation carries).
#'   The list is empty when the model has no such term.
#'
#' @seealso [statmod_fit_structural()]
#'
#' @keywords internal
statmod_structural_par <- function(spec, design) {
  su <- attr(design, "structural")
  sst <- statmod_structural_state(design)
  if (is.null(sst)) return(list())
  out <- list()
  for (u in su) {
    tm <- spec@terms[[u$param]][[u$term]]
    z <- sst$zeta[[u$term]]
    out[[u$term]] <- list(
      parameter = unlist(structural_psi(tm, z)),
      unconstrained = z,
      equation = u$param,
      # what the equation's own intercept carries instead
      held = sst$held[[u$term]])
  }
  out
}


#' Fit the Structural Terms' Own Parameters
#'
#' @description
#' Estimates each structural term's parameters on the unconstrained scale of
#' [modelterms7::term_links()], while the coefficients stay where they are,
#' by a general optimizer with the exact gradient.
#'
#' @details
#' The scoring step does not apply to these parameters, because they are not
#' coefficients of a design block: the predictor is a recursion and not a
#' product, so there is no \eqn{X'WX} to invert. A general optimizer
#' ([optimizers7::lbfgs()] by default) is used instead.
#'
#' The gradient is exact. [statmod_structural_score()] chains the filter's
#' own Jacobian, which is already the total derivative of the predictor in
#' those parameters, and the derivative of the term's own penalties is added.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param obj The objective, from [statmod_objective()].
#' @param beta The stacked coefficients, held fixed.
#' @param hyper The hyperparameters.
#' @param optimizer An \pkg{optimizers7} optimizer, or `NULL` for
#'   [optimizers7::lbfgs()].
#' @param verbose Whether to print a line per term.
#'
#' @return A list with `value`, `converged` and `iterations`.
#'
#' @seealso [statmod_structural()]
#'
#' @keywords internal
statmod_fit_structural <- function(spec, design, obj, beta, hyper, optimizer,
                                   verbose = FALSE) {
  sst <- statmod_structural_state(design)
  su <- attr(design, "structural")
  if (is.null(sst)) return(list(value = obj$fn(beta), converged = TRUE,
                                iterations = 0L))
  if (is.null(optimizer)) optimizer <- optimizers7::lbfgs()
  cf <- obj$split(beta)
  its <- 0L
  ok <- TRUE
  for (u in su) {
    if (!nzchar(u$kind)) next
    key <- u$term
    nm <- names(sst$zeta[[key]])
    # a level a linear intercept already carries is held, not estimated
    free <- setdiff(nm, sst$held[[key]])
    # A COORDINATE A KINKED PENALTY COVERS IS LEFT HERE, and that is measured
    # rather than assumed. [sparse_fit_structural()] fits it afterwards in
    # the same pass and has the last word, so holding it here changes nothing:
    # swept over `lambda` of 1, 5, 20 and 80 on the panel that route's page
    # describes, the log-likelihood, all ten parameters, the count of exact
    # zeros and the convergence flag are identical to the printed digit with
    # the hold and without it. What the hold would cost is one walk over the
    # model's penalties per call, on every fit carrying a structural term
    # whether or not any penalty there has a kink.
    if (!length(free)) next
    set <- function(z) {
      v <- sst$zeta[[key]]
      v[free] <- as.numeric(z)
      sst$zeta[[key]] <- v
      sst$key <- NULL
      sst$value <- NULL
    }
    z0 <- sst$zeta[[key]][free]
    raw <- function(z) {
      set(z)
      -statmod_loglik_at(spec, cf, design) +
        statmod_penalty_at(spec, cf, hyper, design, "value")
    }
    # The objective is evaluated once here, outside the guard below, so that
    # a term or a callback that cannot be evaluated at all raises where it
    # can be read. Inside the search the same failure is a statement about
    # the point: a filter whose loadings put it outside the region where its
    # recursion is bounded diverges, and the density then rejects a
    # non-finite parameter. That is a point the search must step back from,
    # not a run to abandon.
    v0 <- raw(as.numeric(z0))
    fn <- function(z) {
      v <- tryCatch(raw(z), error = function(e) NA_real_)
      if (!is.finite(v)) Inf else v
    }
    gr <- function(z) {
      v <- tryCatch({
        set(z)
        g <- -statmod_structural_score(spec, cf, design)[[key]][free]
        # the objective above adds the term's own penalties, so the gradient
        # adds their derivative: two functions that differ by a penalty are
        # not each other's gradient, and an optimizer handed the pair walks
        # until its budget runs out
        pg <- statmod_structural_penalty(spec, design, hyper, "gradient")
        if (!is.null(pg[[key]])) g <- g + pg[[key]][free]
        g
      }, error = function(e) NULL)
      if (is.null(v) || !all(is.finite(v))) numeric(length(z)) else v
    }
    res <- optimizers7::minimize(optimizer, fn, as.numeric(z0), gr = gr)
    set(res@par)
    its <- its + res@iterations
    ok <- ok && isTRUE(res@converged)
    if (verbose) {
      vb_say("%s: %d iterations, converged %s", short_keys(key),
             as.integer(res@iterations), res@converged)
    }
  }
  list(value = obj$fn(beta), converged = ok, iterations = its)
}

#' Fit a Kinked Penalty Over a Structural Term's Own Parameters
#'
#' @description
#' Runs a proximal gradient iteration over the parameters of a score-driven
#' filter that a penalty with a kink covers, while the coefficients and the
#' term's other parameters stay where the rest of the pass left them.
#'
#' @details
#' The two routes that [sparse_fit()] takes for a block of coefficients cannot
#' be taken here. A coordinate descent reads the block's columns and its running
#' residual, and a structural term contributes no column; the proximal branch
#' reads `beta[block$index]`, and such a penalty covers no entry of the
#' stacked vector. The parameters are in the design's structural state
#' instead, so the iteration is written on that subvector: the smooth part is
#' the negative log-likelihood from which this penalty is taken back out, its
#' gradient is [statmod_structural_score()]'s with the same correction, and the
#' operator is the penalty's own.
#'
#' The iteration stops when the gradient criterion `crit_grad(tol)` is met. A
#' coordinate at exactly zero is one whose gradient of the smooth part lies
#' inside the interval that the kink opens.
#'
#' The hyperparameter is held. A prediction-error or cross-validation
#' criterion ([aic()], [bic()] or [cv()]) that would select it is rejected by
#' [assert_criterion_reach()], because such a criterion works only over the
#' coefficients, where this penalty covers nothing.
#'
#' The convergence flag of a fit that uses this route is the flag of the
#' alternation. Beside a structural term or a kinked block the alternation
#' requires the point to be settled in the sense of [alternation_settled()]
#' (the Newton decrement of the smooth part and the decrement of the kinked
#' zeros both at most [mode_error_limit()]).
#'
#' @param obj The full objective.
#' @param beta The stacked coefficients, returned unchanged: this route moves
#'   the structural state and no coefficient.
#' @param block One entry of `statmod_blocks()$sparse`, whose `structural` is
#'   `TRUE`.
#' @param hyper The hyperparameters.
#' @param spec A [StatmodSpec()].
#' @param design The design, whose structural state this writes into.
#' @param maxit The iteration budget.
#' @param tol The stopping tolerance.
#' @param verbose Whether the optimizer prints its own trace.
#'
#' @return A list shaped like [sparse_fit()]'s.
#'
#' @seealso [sparse_fit()], [statmod_fit_structural()], which fits the term's
#'   other parameters in the same pass.
#'
#' @keywords internal
sparse_fit_structural <- function(obj, beta, block, hyper, spec, design,
                                  maxit = 500, tol = 1e-8, verbose = FALSE) {
  none <- list(par = beta, value = obj$fn(beta), converged = TRUE,
               method = "proximal gradient (structural)", iterations = 0L)
  if (is.null(spec) || is.null(design)) return(none)
  sst <- statmod_structural_state(design)
  key <- block$zterm
  if (is.null(sst) || is.null(key) || is.null(sst$zeta[[key]])) return(none)
  nm <- names(sst$zeta[[key]])
  free <- setdiff(nm, sst$held[[key]])
  # a coordinate a linear intercept already carries is held, and a penalty
  # over it has nothing to shrink
  pos <- match(nm[block$cols], free)
  pos <- pos[!is.na(pos)]
  if (!length(pos)) return(none)
  cf <- obj$split(beta)
  th <- as.list(hyper[[block$param]][[block$term]])
  pen <- block$penalty
  set <- function(z) {
    v <- sst$zeta[[key]]
    v[free] <- as.numeric(z)
    sst$zeta[[key]] <- v
    sst$key <- NULL
    sst$value <- NULL
  }
  z0 <- as.numeric(sst$zeta[[key]][free])
  # the objective minus THIS penalty, which the operator applies instead
  smooth_fn <- function(z) {
    v <- tryCatch({
      set(z)
      -statmod_loglik_at(spec, cf, design) +
        statmod_penalty_at(spec, cf, hyper, design, "value") -
        penalties7::penalty_value(pen, z[pos], th)
    }, error = function(e) NA_real_)
    if (!is.finite(v)) Inf else v
  }
  smooth_gr <- function(z) {
    v <- tryCatch({
      set(z)
      g <- -statmod_structural_score(spec, cf, design)[[key]][free]
      pg <- statmod_structural_penalty(spec, design, hyper, "gradient")
      if (!is.null(pg[[key]])) g <- g + pg[[key]][free]
      g[pos] <- g[pos] - penalties7::penalty_gradient(pen, z[pos], th)
      g
    }, error = function(e) NULL)
    if (is.null(v) || !all(is.finite(v))) numeric(length(z)) else v
  }
  prox <- function(v, step) {
    out <- v
    out[pos] <- penalties7::penalty_prox(pen, v[pos], step, th)
    out
  }
  gval <- function(z) penalties7::penalty_value(pen, z[pos], th)
  opt <- optimizers7::prox_grad(prox = prox, g = gval,
                                criterion = optimizers7::crit_grad(tol),
                                maxit = maxit, verbose = verbose)
  res <- tryCatch(optimizers7::minimize(opt, smooth_fn, z0, gr = smooth_gr),
                  error = function(e) NULL)
  if (is.null(res)) return(none)
  set(res@par)
  list(par = beta, value = obj$fn(beta), converged = isTRUE(res@converged),
       method = "proximal gradient (structural)",
       iterations = res@iterations)
}


#' The Penalty Over a Structural Term's Free Parameters, as a Block
#'
#' @description
#' Returns the Hessian of whatever penalty covers a structural term's own
#' parameters, over the free ones and in their order, which is the order in
#' which the tail of [statmod_full_information()] carries them.
#'
#' @details
#' It is written once because its three readers (the variance matrix of
#' [vcov.StatmodFit()], the marginal criterion and
#' [statmod_structural_table()]) would otherwise each place the block for
#' themselves. The result is `NULL` where the structural terms carry no
#' penalty, where the only such penalty belongs to a mixed covariance class
#' (placed by `joint_penalty_at()`), or where every term has a number of free
#' parameters other than `nfree`, and the caller then pads with zeros.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param hyper The hyperparameters.
#' @param nfree The number of free parameters that the caller expects, or
#'   `NULL` to accept whatever number the term has.
#'
#' @return A square matrix, or `NULL`.
#'
#' @seealso [statmod_structural_penalty()]
#'
#' @keywords internal
structural_penalty_block <- function(spec, design, hyper, nfree = NULL) {
  ph <- statmod_structural_penalty(spec, design, hyper, "hessian")
  if (!length(ph)) return(NULL)
  st <- statmod_structural_state(design)
  for (nm in names(ph)) {
    free <- setdiff(names(st$zeta[[nm]]), st$held[[nm]])
    if (!is.null(nfree) && length(free) != nfree) next
    return(ph[[nm]][free, free, drop = FALSE])
  }
  NULL
}


#' Fit the Coefficients and a Structural Term's Parameters in One System
#'
#' @description
#' Minimizes the objective over the stacked coefficients of every equation
#' and the free parameters of a structural term (of the filter shape or of the
#' likelihood shape) in one system, instead of alternating between the two
#' blocks.
#'
#' @details
#' The default optimizer is [optimizers7::newton()], with the exact gradient
#' of both blocks and the exact observed information over both together
#' ([statmod_full_information()], plus the penalty's Hessian), run to its own
#' stopping rule. An alternation would hand the term's parameters to an
#' optimizer of their own at every pass, and every iteration of that
#' optimizer would run the recursion and its adjoint again while the
#' coefficients are held at a point that is about to move.
#'
#' The unknowns are in the same order as in the information (the coefficients
#' first, then the term's free parameters), so the system needs no
#' permutation; a level that an intercept in the same equation carries is
#' held and leaves the system, exactly as it leaves the information.
#'
#' The objective is evaluated once outside the guard, so that a term that
#' cannot be evaluated at all raises where it can be read; inside the search
#' a non-finite value belongs to the point (a filter whose loadings put it
#' outside the region where its recursion is bounded), and the search steps
#' back from that point instead of abandoning the run.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param obj The objective, as [statmod_objective()] returns it.
#' @param beta The coefficients from which to start.
#' @param hyper The hyperparameters.
#' @param optimizer An `optimizers7` optimizer, or `NULL` for
#'   Newton's method with the exact joint information.
#' @param verbose Logical; report the run.
#'
#' @return A list with `par` (the stacked coefficients only; the filter's free
#'   parameters are written into the design's structural state), `value`,
#'   `converged` and `iterations`.
#'
#' @seealso [statmod_fit_structural()], which fits only the term's
#'   parameters, with the coefficients held.
#'
#' @keywords internal
statmod_fit_joint <- function(spec, design, obj, beta, hyper,
                              optimizer = NULL, verbose = FALSE) {
  jp <- statmod_joint_pieces(spec, design, obj, hyper,
                             kinds = c("filter", "loglik"))
  u0 <- c(beta, jp$zeta())
  jp$raw(u0)
  if (is.null(optimizer)) optimizer <- optimizers7::newton()
  # A FREE LEVEL IS SEARCHED ON THE LEVEL AROUND WHICH THE PREDICTOR
  # FLUCTUATES, L = omega / (1 - sum phi), and not on omega: with omega a
  # coordinate, raising the persistence at fixed omega moves that level, so
  # the likelihood runs along the ridge omega = L (1 - sum phi) and Newton
  # crawled on it. On the Nile, flow ~ 0 + gas(1, 1) stopped after 200
  # iterations at -642.67 (beta 0.22, a negative eigenvalue) from both
  # starts, where the intercept spelling of the same model, whose level IS
  # a coefficient, converges at -637.40.
  lv <- joint_level_chart(spec, jp)
  if (is.null(lv)) {
    res <- optimizers7::minimize(optimizer, jp$fn, u0, gr = jp$gr, he = jp$he)
    par <- res@par
  } else {
    res <- optimizers7::minimize(
      optimizer, function(v) jp$fn(lv$to(v)), lv$from(u0),
      gr = function(v) lv$grad(v, jp$gr(lv$to(v))),
      he = function(v) {
        u <- lv$to(v)
        lv$hess(v, jp$gr(u), jp$he(u))
      })
    par <- lv$to(res@par)
  }
  jp$setz(par[jp$ix])
  if (verbose) {
    vb_say("%d coefficients and %d parameters of %s: %d iterations, converged %s",
           jp$nb, length(jp$free), short_keys(jp$key),
           as.integer(res@iterations), res@converged)
  }
  list(par = par[seq_len(jp$nb)], value = res@value,
       converged = isTRUE(res@converged), iterations = res@iterations)
}


#' How Far Above Its Mode a Filter Fit Stopped
#'
#' @description
#' Computes the Newton decrement \eqn{\tfrac12 g^\top H^{-1} g} of the joint
#' objective over the coefficients and the free parameters of a score-driven
#' term, at the point that an inner fit returned.
#'
#' @details
#' [statmod()] reads it, beside the optimizer's flag, to determine whether
#' the fit from the term's own start failed: where the filter explodes the
#' objective is flat to rounding, the flag at that point depends on the
#' platform, and the decrement is many orders of magnitude above
#' [mode_error_limit()]. The coefficients that the specification holds (those
#' that the outer criterion estimates) are left out, because their score is
#' not a residual of the inner fit. A gradient or a curvature that is not
#' finite, or a curvature that cannot be solved, gives `Inf`.
#'
#' @param spec A [StatmodSpec()].
#' @param design Its design, whose structural state holds the filter's
#'   parameters at the point.
#' @param res The inner fit's result, a list with `par` (the stacked
#'   coefficients) and `obj` (the objective).
#' @param hyper The hyperparameters.
#'
#' @return A single non-negative number, or `Inf`.
#'
#' @seealso [statmod_joint_pieces()], [mode_error_limit()]
#'
#' @keywords internal
filter_decrement <- function(spec, design, res, hyper) {
  if (is.null(res$obj) || is.null(res$par)) return(Inf)
  v <- tryCatch({
    jp <- statmod_joint_pieces(spec, design, res$obj, hyper)
    u <- c(res$par, jp$zeta())
    keep <- setdiff(seq_along(u), held_stack(spec, design))
    g <- jp$gr(u)[keep]
    H <- as.matrix(jp$he(u))[keep, keep, drop = FALSE]
    if (!all(is.finite(g)) || !all(is.finite(H))) Inf else
      0.5 * sum(g * solve(H, g))
  }, error = function(e) Inf)
  if (length(v) == 1L && is.finite(v)) abs(v) else Inf
}


#' The Level Chart of a Filter's Joint Fit
#'
#' @description
#' The change of coordinates under which [statmod_fit_joint()] searches a
#' free level: the level \eqn{\omega} of a score-driven term is replaced by
#' \eqn{L = \omega / \prod_k (1 - r_k)}, the level around which the
#' predictor fluctuates, where \eqn{r_k} are the partial autocorrelations that
#' carry the persistence (by the Durbin-Levinson recursion,
#' \eqn{1 - \sum_k \xi_k = \prod_k (1 - r_k)}, where \eqn{\xi_k} are the
#' autoregressive coefficients that the term reports). The other coordinates
#' are unchanged.
#'
#' @details
#' The chart applies where the level is one free coordinate that is not
#' developed, the link is the identity, the term has at least one persistence,
#' and every persistence is a free scalar partial autocorrelation; in any
#' other case the result is `NULL` and the fit runs on the term's own
#' coordinates. With \eqn{\omega = L P(\rho)}, the gradient is
#' \eqn{J^\top g} and the Hessian
#' \eqn{J^\top H J + g_\omega \nabla^2 \omega}, \eqn{J} the Jacobian of the
#' map.
#'
#' @param spec A [StatmodSpec()].
#' @param jp The joint pieces, as [statmod_joint_pieces()] returns them.
#'
#' @return `NULL`, or a list of the functions `to` (chart to the term's
#'   coordinates), `from`, `grad` and `hess`.
#'
#' @keywords internal
joint_level_chart <- function(spec, jp) {
  tm <- spec@terms[[jp$param]][[jp$key]]
  lvl <- tryCatch(modelterms7::term_level_param(tm), error = function(e) NULL)
  if (length(lvl) != 1L || !lvl %in% jp$free) return(NULL)
  prm <- modelterms7::term_params(tm)
  # a level developed by a subformula (`omega ~ 1 + random(...)`) is left on
  # its own coordinates: its constant column is one of several, and with the
  # chart the outer gradient of a covariance class over it lost its O(h^2)
  if (any(startsWith(prm, paste0(lvl, "."))) || grepl(".", lvl, fixed = TRUE)) {
    return(NULL)
  }
  pacf <- grep("^pacf", prm, value = TRUE)
  if (!length(pacf) || !all(grepl("^pacf[0-9]+$", pacf)) ||
      !all(pacf %in% jp$free)) {
    return(NULL)
  }
  lks <- modelterms7::term_links(tm)
  x <- c(-1.5, 0, 2.5)
  if (!isTRUE(all.equal(linkfunctions7::linkinv(lks[[lvl]], x), x)) ||
      !isTRUE(all.equal(linkfunctions7::dlinkinv(lks[[lvl]], x), rep(1, 3)))) {
    return(NULL)
  }
  il <- jp$nb + match(lvl, jp$free)
  ip <- jp$nb + match(pacf, jp$free)
  # P = prod (1 - r_k) and its first two derivatives in the unconstrained
  # partial autocorrelations; each written as a product over the other
  # factors, which stays finite where one factor is zero
  pr <- function(v) {
    rho <- v[ip]
    k <- length(rho)
    r <- p1 <- p2 <- numeric(k)
    for (j in seq_len(k)) {
      r[j] <- linkfunctions7::linkinv(lks[[pacf[j]]], rho[j])
      p1[j] <- -linkfunctions7::dlinkinv(lks[[pacf[j]]], rho[j])
      p2[j] <- -linkfunctions7::d2linkinv(lks[[pacf[j]]], rho[j])
    }
    p <- 1 - r
    rest <- function(drop) prod(p[-drop])
    d1 <- vapply(seq_len(k), function(j) p1[j] * rest(j), numeric(1))
    d2 <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in seq_len(k)) {
      d2[i, j] <- if (i == j) p2[i] * rest(i) else p1[i] * p1[j] * rest(c(i, j))
    }
    list(P = prod(p), d1 = d1, d2 = d2)
  }
  list(
    to = function(v) {
      v[il] <- v[il] * pr(v)$P
      v
    },
    from = function(u) {
      u[il] <- u[il] / pr(u)$P
      u
    },
    grad = function(v, g) {
      q <- pr(v)
      out <- g
      out[il] <- g[il] * q$P
      out[ip] <- g[ip] + g[il] * v[il] * q$d1
      out
    },
    hess = function(v, g, H) {
      q <- pr(v)
      J <- diag(length(v))
      J[il, il] <- q$P
      J[il, ip] <- v[il] * q$d1
      out <- crossprod(J, as.matrix(H) %*% J)
      K <- matrix(0, length(v), length(v))
      K[il, ip] <- q$d1
      K[ip, il] <- q$d1
      K[ip, ip] <- v[il] * q$d2
      out + g[il] * K
    }
  )
}


#' The Joint Objective Over the Coefficients and a Filter's Parameters
#'
#' @description
#' Returns, as closures, the objective, the gradient and the exact
#' information with which [statmod_fit_joint()] minimizes over the
#' coefficients followed by the filter's free parameters, so that
#' [statmod_certificate()] reads the same ones at a fitted point.
#'
#' @param spec,design,obj,hyper As [statmod_fit_joint()] takes them.
#' @param kinds The kinds of structural unit from which to take the term:
#'   `"filter"` for the joint step, and `"loglik"` as well where a mixture
#'   over latent states is read at a fitted point.
#'
#' @return A list with `fn`, `gr`, `he`, `raw`, `setz` (writes the filter's
#'   free parameters into the design's structural state), `zeta` (reads
#'   them), `ix`, `nb`, `free`, `key` and `param` (the distribution parameter
#'   in which the term sits).
#'
#' @keywords internal
statmod_joint_pieces <- function(spec, design, obj, hyper, kinds = "filter") {
  sst <- statmod_structural_state(design)
  su <- Filter(function(u) u$kind %in% kinds, attr(design, "structural"))
  key <- su[[1L]]$term
  nm <- names(sst$zeta[[key]])
  free <- setdiff(nm, sst$held[[key]])
  nb <- sum(obj$npar)
  ix <- nb + seq_along(free)

  setz <- function(z) {
    v <- sst$zeta[[key]]
    v[free] <- as.numeric(z)
    sst$zeta[[key]] <- v
    sst$key <- NULL
    sst$value <- NULL
  }
  njoint <- nb + length(free)
  # A COVARIANCE CLASS SPLIT BETWEEN THE TWO HALVES is placed by one function
  # rather than by the two blocks below: its Hessian has a cross block, which
  # is the only place the correlation between a coefficient and a filter's own
  # parameter enters, and neither half's assembly has anywhere to put it.
  mixed <- any(vapply(statmod_penalized(spec, design),
                      function(z) isTRUE(z$mixed), TRUE))
  raw <- function(u) {
    setz(u[ix])
    # the VALUE of a mixed class arrives through obj$fn, which reads
    # statmod_penalty_at(); only its derivatives are placed here
    obj$fn(u[seq_len(nb)])
  }
  fn <- function(u) {
    v <- tryCatch(raw(u), error = function(e) NA_real_)
    if (!is.finite(v)) Inf else v
  }
  gr <- function(u) {
    v <- tryCatch({
      b <- u[seq_len(nb)]
      setz(u[ix])
      g <- -statmod_structural_score(spec, obj$split(b), design)[[key]][free]
      pg <- statmod_structural_penalty(spec, design, hyper, "gradient")
      if (!is.null(pg[[key]])) g <- g + pg[[key]][free]
      out <- c(obj$gr(b), g)
      if (mixed) {
        # obj$gr() already carries the coefficient part of a mixed class
        jg <- joint_penalty_at(spec, design, obj$split(b), hyper,
                               "gradient", njoint)
        jg[seq_len(nb)] <- 0
        out <- out + jg
      }
      out
    }, error = function(e) NULL)
    if (is.null(v) || !all(is.finite(v))) numeric(length(u)) else v
  }
  he <- function(u) {
    b <- u[seq_len(nb)]
    setz(u[ix])
    cf <- obj$split(b)
    H <- statmod_full_information(spec, cf, design)
    P <- matrix(0, nrow(H), ncol(H))
    # as_dense() for the reason statmod_full_information() records: this
    # matrix spans the coefficients AND the term's own parameters, so it is
    # dense by construction, and writing a sparse block into a slice of a
    # dense one signals that the number of items to replace is not a multiple
    # of the replacement length -- from inside the objective, naming nothing
    # a caller wrote.
    P[seq_len(nb), seq_len(nb)] <-
      as_dense(statmod_penalty_at(spec, cf, hyper, design, "hessian"))
    ph <- statmod_structural_penalty(spec, design, hyper, "hessian")
    if (!is.null(ph[[key]])) P[ix, ix] <- ph[[key]][free, free, drop = FALSE]
    if (mixed) {
      P <- P + joint_penalty_at(spec, design, cf, hyper, "hessian", njoint)
    }
    H + P
  }

  list(fn = fn, gr = gr, he = he, raw = raw, setz = setz,
       zeta = function() as.numeric(sst$zeta[[key]][free]),
       ix = ix, nb = nb, free = free, key = key, param = su[[1L]]$param)
}





#' Run the Structural Terms at the Current Parameters
#'
#' @description
#' Returns the level that each filter adds to the equation in which it sits,
#' together with the derivative of that level in the term's own parameters.
#'
#' @details
#' The function keeps no memo of its own. The caller, [statmod_eta()], keeps
#' one keyed on the coefficients and the term's parameters, because the
#' objective, its gradient and its curvature are evaluated at the same point
#' in turn and a filter is the expensive part of each.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param eta_static The static predictors, one per distribution parameter.
#' @param theta_static The parameters that they give.
#'
#' @return A named list, one entry per structural term, carrying its
#'   predictor, the derivative of that predictor in the term's own
#'   parameters, and the curvature that the recursion read at each step.
#'
#' @seealso [statmod_structural()]
#'
#' @keywords internal
statmod_filter_at <- function(spec, design, eta_static, theta_static) {
  su <- attr(design, "structural")
  st <- statmod_structural_state(design)
  out <- list()
  for (u in su) {
    if (!identical(u$kind, "filter")) next
    tm <- spec@terms[[u$param]][[u$term]]
    if (is_joint_filter(tm)) {
      out[[u$term]] <- joint_filter_entry(
        spec, u, tm, structural_psi(tm, st$zeta[[u$term]]), eta_static)
      next
    }
    cb <- structural_callbacks(spec, theta_static, u$param,
                               scaling = filter_scaling(tm))
    psi <- structural_psi(tm, st$zeta[[u$term]])
    f <- modelterms7::term_filter(tm, eta_static[[u$param]], spec@response,
                                  cb$score, cb$curvature, psi,
                                  fast = cb$fast, threads = spec@threads)
    out[[u$term]] <- list(param = u$param, term = u$term, tm = tm,
                          eta = f$eta, jacobian = f$jacobian,
                          curv = f$curv,
                          psi = psi, cb = cb,
                          eta_static = eta_static[[u$param]])
  }
  out
}


#' A Structural Term's Parameters, With Standard Errors
#'
#' @description
#' Returns the estimated parameters of every structural term, on their own
#' scale, together with standard errors and intervals that are taken from
#' the joint observed information.
#'
#' @details
#' A structural term contributes no design block, so the coefficient tables
#' do not report it. [statmod_full_information()] spans the coefficients and
#' the term's own parameters, and the tail block of its inverse (with the
#' penalty Hessians added) is the variance of those parameters on the
#' unconstrained scale. A covariance class split between the coefficients and
#' the term's own parameters is placed by `joint_penalty_at()`, as in
#' [vcov.StatmodFit()].
#'
#' The quantities reported are those of [modelterms7::term_readable()], which
#' are not always the coordinates on which the term was estimated (a
#' persistence of a score-driven term is an autoregressive coefficient, a
#' function of the partial autocorrelations). The standard error of a quantity
#' is the delta method over the joint variance of the term's parameters,
#' \eqn{\sqrt{J V J^\top}}, with \eqn{J} the Jacobian that
#' `term_readable()` returns (a persistence of a term with \eqn{q > 1} reads
#' several coordinates). The interval is built on the scale that the term
#' declares for the quantity (for a score-driven term, the log scale for a
#' positive loading, the scale of the partial autocorrelation's link for the
#' persistence when \eqn{q = 1}, and the identity scale for a persistence when
#' \eqn{q > 1}): the standard error is carried to that scale by the derivative
#' of its link and the end points are mapped back. The interval of a
#' persistence is therefore confined to \eqn{(-1, 1)} only when
#' \eqn{q = 1}.
#'
#' A level held because an intercept in the same equation carries it is not
#' estimated, and its standard error is reported as missing, as is that of any
#' quantity whose Jacobian row reads a held parameter.
#'
#' A quantity that depends on a coordinate past the edge of its chart (see
#' [structural_edge_rows()]) is reported without a standard error and an
#' interval, for the reason that a correlation at the boundary is: the
#' chart's derivative tends to zero there, so the delta method reports a
#' standard error that collapses exactly where the quantity stops being
#' identified.
#'
#' @param fit A [StatmodFit()].
#' @param level The interval's level, `0.95` by default.
#'
#' @return A data frame with one row per quantity that each structural term
#'   reports and the columns `parameter`, `term`, `name`, `estimate`, `se`,
#'   `lower`, `upper`, `held` (the quantity is a held parameter), `edge`
#'   (it depends on a coordinate past the edge of its chart), `component` and
#'   `position` (the one of the term's own parameters to which the quantity
#'   belongs, and where in the parameter vector it sits; both are read off the
#'   Jacobian's support). The value is `NULL` where the model carries no
#'   structural term.
#'
#' @seealso [statmod_full_information()] for the joint matrix inverted,
#'   [coef.StatmodFit()] for the same parameters without their uncertainty
#'
#' @keywords internal
statmod_structural_table <- function(fit, level = 0.95) {
  spec <- fit@spec
  design <- statmod_design(spec)
  su <- attr(design, "structural")
  if (is.null(su) || !length(su)) return(NULL)

  nb <- sum(vapply(design, function(d) d$npar, integer(1)))
  V <- tryCatch({
    I <- statmod_full_information(spec, fit@coefficients, design)
    # with a penalty over the term's own parameters the matrix to invert is
    # the penalized one: the deviations of a panel are identified by that
    # penalty and by nothing else, so the unpenalized information is singular
    # along the direction that trades a population value against them
    ps <- structural_penalty_block(spec, design, fit@hyper, nrow(I) - nb)
    if (!is.null(ps)) {
      ix <- nb + seq_len(nrow(I) - nb)
      I[ix, ix] <- I[ix, ix] + ps
    }
    # AND THE DESIGN'S OWN PENALTIES, for the same reason. A random effect is
    # identified by its penalty and by nothing else, so the unpenalized
    # information is singular along the direction that trades its population
    # value against its deviations; with a penalized block beside a
    # structural term that direction is in this matrix, and inverting it left
    # every structural standard error missing. Measured on a break-point term
    # with a random development beside a score-driven level: the smallest
    # eigenvalue is 2.7e-13 and its direction is the six deviations at +0.378
    # against the intercept at -0.378, exactly the flat one.
    Sd <- tryCatch(statmod_penalty_at(spec, fit@coefficients, fit@hyper,
                                      design, "hessian"),
                   error = function(e) NULL)
    if (!is.null(Sd) && nrow(Sd) == nb) {
      I[seq_len(nb), seq_len(nb)] <- I[seq_len(nb), seq_len(nb)] +
        as.matrix(Sd)
    }
    # and a covariance class split between the coefficients and the term's
    # own parameters, placed whole by the function that vcov() reads: its
    # cross block is what identifies the two halves together
    I <- I + as.matrix(joint_penalty_at(spec, design, fit@coefficients,
                                        fit@hyper, "hessian", nrow(I)))
    # the slot of a held break-point is a column of zeros by design, as
    # vcov() records, and is left out of the inverse rather than making the
    # whole matrix singular: beside a jump() the filter's loading and
    # persistence lost their standard errors to that one zero row (pig2 with
    # jump(tt) + gas(1, 1) | alpha ~ s(tt))
    keep <- setdiff(seq_len(nrow(I)), term_held_stack(spec, design))
    Vk <- matrix(0, nrow(I), nrow(I))
    Vk[keep, keep] <- solve(I[keep, keep, drop = FALSE])
    Vk
  }, error = function(e) NULL)
  z <- stats::qnorm(1 - (1 - level) / 2)

  rows <- list()
  for (u in su) {
    tm <- spec@terms[[u$param]][[u$term]]
    st <- statmod_structural_state(design)
    zeta <- st$zeta[[u$term]]
    held <- st$held[[u$term]]
    nm <- modelterms7::term_params(tm)
    free <- setdiff(nm, held)
    # What is reported is what the TERM says it reports, which is not always
    # the coordinate it was estimated on: a score-driven persistence rides a
    # partial autocorrelation, and the literature's xi_j is the
    # autoregressive coefficient, a function of the whole chart. The Jacobian
    # comes with it, so the standard error is the delta method over the JOINT
    # variance rather than one entry of its diagonal -- above q = 1 a
    # coefficient reads several coordinates and their covariance matters.
    rd <- modelterms7::term_readable(tm, zeta)
    # the joint matrix drops a held level, so its tail is the free ones in
    # the term's own order
    Vz <- matrix(0, length(nm), length(nm), dimnames = list(nm, nm))
    ok <- !is.null(V) && nrow(V) == nb + length(free)
    if (ok) {
      Vz[free, free] <- V[nb + seq_along(free), nb + seq_along(free),
                          drop = FALSE]
    }
    J <- rd$jacobian
    vq <- if (ok) diag(J %*% Vz %*% t(J)) else rep(NA_real_, length(nm))
    # WHICH of the term's own parameters a reported quantity belongs to, and
    # WHERE in the parameter vector it sits, read off the Jacobian's support
    # rather than off the quantity's name. A quantity reading coordinates of
    # several parameters belongs to none of them and stays with the term: an
    # autoregressive coefficient is a function of the whole chart, which is
    # why it is not filed under one partial autocorrelation.
    cps <- tryCatch(modelterms7::term_components(tm),
                    error = function(e) list())
    owner <- function(i) {
      k <- which(J[i, ] != 0)
      if (!length(k)) return(list("", NA_integer_))
      hit <- vapply(cps, function(z) all(k %in% z$index), logical(1))
      if (sum(hit) != 1L) return(list("", NA_integer_))
      list(names(cps)[which(hit)],
           if (length(k) == 1L) as.integer(k) else NA_integer_)
    }
    # a held parameter is not estimated, and a quantity that reads one is not
    # either: it would be reported with the variance of the rest alone
    touches_held <- apply(J[, held, drop = FALSE] != 0, 1L, any)
    # WHICH ROWS ARE A HELD PARAMETER, read off the Jacobian's support and not
    # by position: a term may report more quantities than it has parameters
    # (a regime chain, k^2 probabilities over k(k-1) log-ratios), so row i is
    # not parameter i
    is_held <- apply(J != 0, 1L, function(r) any(r) && all(nm[r] %in% held))
    edge_q <- structural_edge_rows(tm, zeta, nm)
    for (i in seq_along(rd$name)) {
      s <- if (ok && !touches_held[[i]] && !edge_q[[i]]) sqrt(max(vq[[i]], 0))
           else NA_real_
      val <- rd$value[[i]]
      lo <- hi <- NA_real_
      if (is.finite(s)) {
        # the interval is built on the scale that keeps the quantity in its
        # own set and mapped back, as every interval in the toolkit is
        g <- rd$scale[[i]]
        t0 <- linkfunctions7::linkfun(g, val)
        st <- s * abs(linkfunctions7::dlinkfun(g, val))
        ends <- sort(c(linkfunctions7::linkinv(g, t0 - z * st),
                       linkfunctions7::linkinv(g, t0 + z * st)))
        lo <- ends[1L]
        hi <- ends[2L]
      }
      ow <- owner(i)
      ds <- structural_display(tm, u$param, rd$name[[i]])
      nm_i <- if (is_joint_filter(tm)) ds$name else rd$name[[i]]
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = ds$param, term = u$term, name = nm_i,
        estimate = as.numeric(val), se = s, lower = lo, upper = hi,
        held = is_held[[i]], edge = edge_q[[i]], component = ow[[1L]],
        position = ow[[2L]], stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}


#' The Quantities of a Structural Term That Read a Coordinate at an Edge
#'
#' @description
#' Marks each quantity that [modelterms7::term_readable()] reports for a
#' structural term according to whether it depends on a free coordinate that
#' has run past the edge of a chart mapping onto a bounded set.
#'
#' @details
#' The coordinates checked are those that [modelterms7::term_charted()]
#' names, and a coordinate is at the edge when its free value exceeds `edge`
#' in absolute value, the rule that [edge_violations()] and
#' [statmod_certificate()] use, on a side where its chart maps onto a finite
#' bound: a log link has an edge towards zero and none towards infinity. The
#' dependence is read off the Jacobian evaluated after those coordinates are
#' moved back to \eqn{\pm}`edge`: at the fitted point the derivative of a
#' saturated chart can underflow to zero, which would hide the dependence
#' that it is meant to reveal.
#'
#' @param tm The structural term.
#' @param zeta Its parameters on the unconstrained scale.
#' @param nm Its parameter names, as [modelterms7::term_params()] gives them.
#' @param edge The free value past which a coordinate is at the edge.
#'
#' @return A logical vector with one element per reported quantity.
#'
#' @seealso [statmod_structural_table()], [readable_joint()]
#'
#' @keywords internal
structural_edge_rows <- function(tm, zeta, nm, edge = 8) {
  z <- unlist(zeta)[nm]
  ch <- intersect(tryCatch(modelterms7::term_charted(tm),
                           error = function(e) character(0)), nm)
  lk <- tryCatch(modelterms7::term_links(tm), error = function(e) list())
  # ON THE SIDE WHERE THE CHART HAS A BOUND. A log link saturates towards
  # zero and not towards infinity, so a loading at a free value of 8.88 (a
  # gas() on the Nile flows, kappa1 = 7188 in units of a variance) is an
  # ordinary estimate. A coordinate whose link is the identity and which is
  # still charted is one whose chart lies elsewhere (the log-ratios of a
  # regime chain), and is bounded on both sides.
  bounded_side <- function(j) {
    g <- lk[[j]]
    if (is.null(g) || identical(g@link_name, "identity")) return(TRUE)
    far_v <- linkfunctions7::linkinv(g, sign(z[[j]]) * 50)
    b <- g@link_bounds
    b <- b[is.finite(b)]
    length(b) > 0L && is.finite(far_v) &&
      min(abs(far_v - b)) <= 1e-8 * max(1, abs(b))
  }
  far <- ch[is.finite(z[ch]) & abs(z[ch]) > edge]
  far <- far[vapply(far, bounded_side, logical(1))]
  rd0 <- modelterms7::term_readable(tm, zeta)
  if (!length(far)) return(rep(FALSE, length(rd0$name)))
  zc <- z
  zc[far] <- sign(z[far]) * edge
  J <- modelterms7::term_readable(tm, as.list(zc))$jacobian
  apply(J[, far, drop = FALSE] != 0, 1L, any)
}


#' The Latent Variables of a Fitted Structural Term
#'
#' @description
#' Returns the posterior summary of the latent variable over which a
#' structural term of the likelihood shape integrates: the posterior mean and
#' standard deviation of each group's break-points for a marginal break-point
#' term
#' ([modelterms7::jump()], [modelterms7::seg()] or [modelterms7::jseg()]
#' with `marginal = TRUE`), and the smoothed probability of each regime at
#' each observation for [modelterms7::regime()].
#'
#' @details
#' The quantities come from the same decomposition on which the marginal
#' likelihood is computed. For a break-point term it is the posterior over a
#' group's intervals or quadrature nodes, and the within-interval moments are
#' those of the fitted prior truncated to the interval. For a regime term it
#' is the forward and backward recursions, which give \eqn{P(S_t = j \mid y_1,
#' \dots, y_n)}. The computation is [modelterms7::term_latent()]'s; this
#' function supplies the quantities that the term does not hold, the fitted
#' predictors and the model's log-density.
#'
#' @param fit A [StatmodFit()] whose model carries a structural
#'   term of the likelihood shape.
#'
#' @return For a break-point term, a data frame with one row per group and
#'   break-point: `group`, `psi` (the index of the break-point within the
#'   term), `mean` and `sd` (the posterior mean and standard deviation of its
#'   position). For a regime term, a data
#'   frame with one row per observation, in the order of the data, and one
#'   column per regime, `state1`, `state2`, and so on.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(id = rep(1:4, each = 8), x = rep(1:8, 4))
#' dd$psi <- 4.5 + rep(rnorm(4, 0, 0.4), each = 8)
#' dd$y <- 1 + 2 * (dd$x >= dd$psi) + rnorm(32, 0, 0.3)
#' fit <- statmod(y ~ jump(x, psi ~ random(~1 | id), marginal = TRUE),
#'                distributions7::gaussian1_distrib(), dd)
#' statmod_latent(fit)
#'
#' # The smoothed probability of the second regime is read at a change of level.
#' dr <- data.frame(t = 1:60, y = c(rnorm(30), rnorm(30, 3)))
#' fr <- statmod(y ~ regime(k = 2, time = t),
#'               distributions7::gaussian1_distrib(), dr)
#' round(statmod_latent(fr)$state2[c(1, 29, 30, 31, 32, 60)], 3)
#'
#' @seealso [modelterms7::term_latent()],
#'   [statmod()]
#'
#' @export
statmod_latent <- function(fit) {
  spec <- fit@spec
  design <- statmod_design(spec)
  ev <- statmod_eta(spec, design, fit@coefficients)
  if (!length(ev$regimes)) {
    stop(paste("the model carries no structural term of the likelihood",
               "shape, so there is no latent variable to report."),
         call. = FALSE)
  }
  r <- ev$regimes[[1L]]
  modelterms7::term_latent(r$tm, r$eta_static, spec@response, r$cb$logdens,
                           r$psi)
}
