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

#' Which Terms Rewrite the Likelihood
#'
#' @description
#' The parameter, name and shape of every structural term in a
#' specification.
#'
#' @details
#' The shape is \code{"filter"} for a term that produces a predictor, which
#' the layer adds to the equation the term sits in, and \code{"loglik"} for
#' one whose contribution is a likelihood mixed over states and cannot be
#' written as a predictor at all.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#'
#' @return A list of entries with \code{param}, \code{term} and \code{kind},
#'   possibly empty.
#'
#' @seealso \code{\link{statmod_filter_at}}
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


#' Which Shape of the Structural Contract a Term Implements
#'
#' @description
#' \code{"filter"} when the term reports a predictor,
#' \code{"loglik"} when it reports a likelihood, and \code{""} when it
#' implements neither.
#'
#' @details
#' The question is asked of the methods the term registers, not of a list of
#' class names, so a structural term written outside the package is routed
#' without an edit here. The class a method was registered on is
#' \code{attr(m, "signature")[[1]]}, compared by name and package.
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
#' The environment a design carries the structural terms' own parameters in,
#' on the unconstrained scale of \code{\link[modelterms7]{term_links}}.
#'
#' @param design The design.
#'
#' @return An environment, or \code{NULL}.
#'
#' @keywords internal
statmod_structural_state <- function(design) attr(design, "structure")


#' Reuse a Structural Quantity Computed at the Same Point
#'
#' @description
#' A depth-one exact memo on the design's structural state: where the last
#' call's key is \code{identical()} to this one's, the stored value is
#' returned; otherwise \code{compute()} runs and replaces it.
#'
#' @details
#' Measured on the gas panel at 60 groups, 62 of the 154 curvature
#' recursions of one fit recompute a point already visited -- the same
#' coefficients and the same term parameters, up to five times each,
#' because the criterion, its gradient and the joint step's curvature all
#' read the same mode -- and the recursion is 35 per cent of the fit. The
#' cache returns the previously computed object itself, so a hit is
#' bit-identical to recomputing by construction, and the key is compared
#' with \code{identical()} on the full numeric inputs, so a collision
#' cannot happen.
#'
#' It stands aside where the design carries REFRESHABLE terms: a
#' break-point block advances its rescaling schedule as the alternation
#' commits, so the same coefficients do not imply the same design there,
#' and a key that cannot see the schedule must not answer.
#'
#' @param design The design whose structural state holds the slots.
#' @param slot A short name, one cache per slot.
#' @param key The exact inputs the value depends on.
#' @param compute A function of no arguments.
#'
#' @return \code{compute()}'s value, possibly from the cache.
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


#' The Parameters a Structural Term Starts From
#'
#' @description
#' What the term itself declares through
#' \code{\link[modelterms7]{term_start}}: as near the model without the
#' term as its charts allow, which only the term can say -- a score
#' loading on the log chart has no coordinate for zero, and starts at a
#' weak response instead.
#'
#' @param term A built structural term.
#' @param target The response on the scale of the term's equation, from
#'   \code{\link{predictor_target}}, or \code{NULL}. A term whose start is
#'   read off the data -- the marginal break-point term's exact profile --
#'   consumes it; every other method ignores it through the dots.
#'
#' @return A named numeric vector.
#'
#' @keywords internal
structural_zeta_start <- function(term, target = NULL) {
  modelterms7::term_start(term, target = target)
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


#' The Score and Curvature a Filter Is Driven By
#'
#' @description
#' The derivative of the log-density in one distribution parameter's
#' predictor, and its second derivative, as functions of that predictor at
#' one observation, with every other parameter held where it is.
#'
#' @details
#' They are read one observation at a time because a score-driven filter
#' evaluates them at the predictor it has just produced, which is not known
#' before the recursion reaches that observation. A term whose callbacks do
#' not depend on the state -- a regime chain, whose levels shift a predictor
#' known in advance -- is given the whole index vector instead.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param theta The per-observation parameters, as
#'   \code{\link{statmod_eta}} returns them.
#' @param p The distribution parameter the term sits in.
#'
#' @return A list with \code{score}, \code{curvature} and \code{logdens}.
#'
#' @keywords internal
structural_callbacks <- function(spec, theta, p) {
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

  list(
    score = function(e, i) as.numeric(k$score(y[i], at(i), e)),
    curvature = function(e, i) as.numeric(k$curvature(y[i], at(i), e)),
    logdens = function(e, i) as.numeric(k$logdens(y[i], at(i), e)),
    # the fast context of piano_parallel.txt section 2a: where the C
    # registries of distributions7 and linkfunctions7 cover this family
    # and this link, the filter's kernel reads the score and the
    # curvature through their scalar entry points instead of calling back
    # into R -- the same composition, bit-identical, held to that by
    # modelterms7's twin test -- and an uncovered pair leaves the context
    # inert, the callbacks above running as before
    fast = if (!S7::S7_inherits(d, distributions7::multivariate_distrib) &&
               !is.matrix(y)) {
      list(family = attr(S7::S7_class(d), "name"),
           link = lk@link_name,
           k = match(p, params),
           bounds = as.numeric(lk@link_bounds),
           y = as.numeric(y),
           theta = theta_n)
    } else NULL
  )
}


#' Run the Structural Terms at the Current Parameters
#'
#' @description
#' The level each filter adds to the equation it sits in, together with the
#' derivative of that level in the term's own parameters.
#'
#' @details
#' The result is memoized on the coefficients and the term's parameters,
#' since the objective, its gradient and its curvature are asked for at the
#' same point in turn and a filter is the expensive part of each.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param eta_static The static predictors, one per distribution parameter.
#' @param theta_static The parameters they give.
#'
#' @return A named list, one entry per structural term, each with
#'   \code{param}, \code{level}, \code{jacobian} and the callbacks used.
#'
#' @seealso \code{\link{statmod_structural}}
#'
#' @keywords internal
#' Run the Regime Terms at the Current Parameters
#'
#' @description
#' The smoothed state probabilities of each latent Markov term, the levels
#' its regimes shift the predictor by, and the per-state predictors those
#' give.
#'
#' @details
#' A term of this shape does not report a predictor: its contribution is a
#' likelihood mixed over states, so the model's log-likelihood comes from
#' \code{\link[modelterms7]{term_loglik}} and not from the density at a
#' single predictor. What makes everything else follow is Fisher's identity
#' -- the derivative of that mixed likelihood in ANY predictor the model
#' carries is the posterior-weighted derivative of the ordinary one, so a
#' caller differentiates its own log-density once per state, vectorized,
#' and weights.
#'
#' The predictor reported for the equation the term sits in is the
#' posterior-weighted one, which is what a fitted value means here.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param eta_static The static predictors.
#' @param theta_static The parameters they give.
#'
#' @return A list, one entry per regime term.
#'
#' @seealso \code{\link{statmod_structural}}
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
#' The level a component adds to the predictor: one number for a regime or
#' a side pattern, one value per observation for a quadrature node.
#'
#' @param mu What \code{\link[modelterms7]{term_levels}} returned.
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
#' @param mu What \code{\link[modelterms7]{term_levels}} returned.
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
#' The distribution's parameters at every observation with the regime term's
#' equation shifted by one state's level.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param eta_static The static predictors.
#' @param param The equation the term sits in.
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
#' list, which is keyed by \code{\link[distributions7]{deriv_names}}.
#'
#' @details
#' The name is BUILT from the parameter names in the family's own order, not
#' parsed out of one, which is the discipline \pkg{distributions7} records
#' for a parameter whose own name contains the separator.
#'
#' @param params The parameter names, in the family's order.
#' @param a,b,c Indices into \code{params}.
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
#' list, built the same way \code{\link{deriv3_key}} builds its own.
#'
#' @details
#' The fourth order is wanted where a filter's third derivative is: each
#' order of differentiating the predictor through the recursion pulls in one
#' more order of the family, the score the recursion is driven by being read
#' at the predictor it produces. Every family of \pkg{distributions7} carries
#' it in closed form.
#'
#' @param params The parameter names, in the family's order.
#' @param a,b,c,d Indices into \code{params}.
#'
#' @return A single string.
#'
#' @keywords internal
deriv4_key <- function(params, a, b, c, d) {
  paste(params[sort(c(a, b, c, d))], collapse = "_")
}


#' The Observed Information Over the Coefficients and a Filter's Parameters
#'
#' @description
#' The information of a model carrying a structural filter, over the
#' coefficients of every equation and the term's own parameters, with the
#' recursion's own curvature in it.
#'
#' @details
#' The gradient of such a model is exact already, through
#' \code{\link[modelterms7]{term_adjoint}}. Its curvature is not: the matrix
#' the scoring step inverts is assembled as though the level were an offset,
#' which is a legitimate scoring matrix and is not the information. Writing
#' \eqn{u} for the coefficients followed by the term's parameters on the
#' unconstrained scale, \eqn{V_{q,t}} for the derivative of equation
#' \eqn{q}'s predictor in \eqn{u}, and \eqn{E_t} for the second derivative
#' of the predictor the filter produces,
#'
#' \deqn{-\frac{\partial^2 \ell}{\partial u \partial u^\top}
#'   = -\sum_t w_t \sum_{q,r} \ell_{qr,t} V_{q,t}^\top V_{r,t}
#'     - \sum_t w_t \ell_{p,t} E_t.}
#'
#' Only the equation carrying the filter has a \eqn{V} that is not its own
#' design: there it is the forward Jacobian of the recursion, which
#' \code{\link[modelterms7]{term_curvature}} returns beside the contracted
#' \eqn{E}. The third derivatives the second sum needs are
#' \pkg{distributions7}'s, in closed form for every family.
#'
#' The term's parameters are the LAST columns, which is the convention
#' \code{term_curvature()} shares with its caller.
#'
#' A term that mixes over latent states rather than shifting the predictor
#' is routed to \code{\link{statmod_regime_information}}, whose Hessian
#' comes from the same forward recursion the likelihood does; a model
#' carrying neither gets the ordinary information.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param coef A named list of coefficient vectors.
#' @param design The design.
#'
#' @return A square matrix over the coefficients and the term's parameters,
#'   or the ordinary information where there is no structural term.
#'
#' @seealso \code{\link{statmod_structural_score}},
#'   \code{\link{statmod_regime_information}}
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
  # The assembly itself is written once, in .structural_blocks(): the exact
  # gradient needs the same cross and M at one order higher, and two copies
  # would drift the moment either was touched. Given no direction it returns
  # exactly these two pieces and never looks at a fourth derivative.
  blocks <- .structural_blocks(params, ap, Vs, H, D3, NULL, n)(NULL)

  # The recursion is re-run here at parameters it has already been run at,
  # so the predictor it reaches is the one the derivatives above were read
  # at and the callbacks can LOOK THEM UP rather than ask the family again.
  # That is the difference between this and the filter itself, where the
  # score is evaluated at a predictor the recursion has just produced and
  # cannot be known in advance: measured, it is where the time goes.
  s_at <- rep_len(gl[[f$param]], n)
  c_at <- rep_len(H[[hess_key(params, ap, ap)]], n)
  # the same recursion at the same point is asked for by the exact
  # gradient's shared parts too, so the memo sits on THIS call and both
  # read one computation; the seed is part of the key, so a caller whose
  # layout differed would miss the cache rather than take the wrong shape
  # the same pieces as DATA beside the callback: with them (and the score
  # and curvature as lookups) an eligible term runs its second-order
  # recursion compiled, with the callback kept for the cases the kernel
  # declines (a developed autoregressive chart, the third order)
  bd_data <- structural_blocks_data(params, ap, Vs, H, D3, n)
  cv <- structural_memo(design, "curv",
                        list(zeta = f$psi, eta = f$eta_static,
                             g = w * s_at, seed = Vs[[ap]]), function() {
    modelterms7::term_curvature(
      f$tm, f$eta_static, spec@response,
      function(e, i) s_at[i], function(e, i) c_at[i], f$psi,
      w * s_at, Vs[[ap]], blocks,
      score_values = s_at, curvature_values = c_at,
      blocks_data = bd_data, threads = spec@threads)
  })
  Vs[[ap]] <- cv$jacobian

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


#' The Observed Information of a Model Carrying a Regime Term
#'
#' @description
#' The observed information over the coefficients of every equation together
#' with the regime term's own parameters, from the exact Hessian of the
#' mixed log-likelihood.
#'
#' @details
#' The matrix \code{\link{statmod_information_at}} assembles for a mixture is
#' the complete-data information, the ordinary one averaged over the
#' smoothed states. That is a legitimate scoring matrix and it is not the
#' observed information: by Louis's missing-information principle the two
#' differ by the conditional variance of the complete-data score, which is
#' positive semidefinite, so the complete-data matrix is the larger and a
#' standard error read off it is too small.
#'
#' \code{\link[modelterms7]{term_hessian}} returns the exact Hessian by
#' propagating first and second derivatives through the same scaled forward
#' recursion that computes the likelihood. What this function supplies is
#' the model's side of that contract: how each equation's unknowns reach
#' each predictor, and the family's first and second derivatives at the
#' predictor each regime shifts to.
#'
#' Unlike the filter's, these derivatives cannot be looked up from a single
#' evaluation: a regime shifts the predictor by a level of its own, so the
#' family is evaluated once per regime, vectorized over observations. That
#' is the same property that made the forward recursion compilable.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param ev The evaluated predictors, from \code{\link{statmod_eta}}.
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
    if (npar[a] > 0L) M[, offs[a] + seq_len(npar[a])] <- design[[params[a]]]$X
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
#' One named vector per structural term, on the scale the term's own
#' documentation describes them on, together with the unconstrained values
#' they were estimated as.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design the fit ran on.
#'
#' @return A named list, one entry per structural term, each a list with
#'   \code{parameter} and \code{unconstrained}. Empty when there is none.
#'
#' @seealso \code{\link{statmod_fit_structural}}
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
#' \code{\link[modelterms7]{term_links}}, the coefficients held where they
#' are, by a general optimizer with the exact gradient.
#'
#' @details
#' It is a general optimizer and not the scoring step because these are not
#' coefficients of a design block: there is no \eqn{X'WX} to invert, the
#' predictor being a recursion rather than a product. The gradient is exact
#' -- \code{\link{statmod_structural_score}} chains the filter's own
#' jacobian, which is the total derivative of the predictor in those
#' parameters -- so a quasi-Newton method has everything it needs.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param obj The objective, from \code{\link{statmod_objective}}.
#' @param beta The stacked coefficients, held fixed.
#' @param hyper The hyperparameters.
#' @param optimizer An \pkg{optimizers7} optimizer, or \code{NULL} for
#'   \code{\link[optimizers7]{lbfgs}}.
#' @param verbose Whether to print a line per term.
#'
#' @return A list with \code{value}, \code{converged} and \code{iterations}.
#'
#' @seealso \code{\link{statmod_structural}}
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

#' The Penalty Over a Structural Term's Free Parameters, as a Block
#'
#' @description
#' The Hessian of whatever penalty covers a structural term's own parameters,
#' over the FREE ones and in their order, which is the order the tail of
#' \code{\link{statmod_full_information}} carries them in.
#'
#' @details
#' Written once because three readers need it and a block each of them placed
#' for itself would agree only by accident: the joint fit, the variance
#' matrix, and the marginal criterion. Where no structural term carries a
#' penalty the answer is \code{NULL}, and a caller pads with zeros as before.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param hyper The hyperparameters.
#' @param nfree The number of free parameters the caller expects, or
#'   \code{NULL} to take whatever the term has.
#'
#' @return A square matrix, or \code{NULL}.
#'
#' @seealso \code{\link{statmod_structural_penalty}}
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


#' Fit the Coefficients and a Filter's Parameters in One System
#'
#' @description
#' One Newton step over the stacked coefficients of every equation and the
#' free parameters of a structural term of the filter shape, instead of
#' alternating between the two blocks.
#'
#' @details
#' The alternation was not a statement about the model: the exact gradient of
#' both blocks and the exact observed information over both together are
#' already available, the second as \code{\link{statmod_full_information}},
#' which was built for \code{\link{vcov.StatmodFit}} and discarded for the
#' fit. What the alternation cost is filter runs -- each sweep handed the
#' term's parameters to an optimizer of their own, and every one of ITS
#' iterations ran the recursion and its adjoint again, with the coefficients
#' held at a point that was about to move.
#'
#' The unknowns are ordered as the information orders them, the coefficients
#' first and the term's free parameters after, so no permutation is needed;
#' a level an intercept in the same equation carries is held and leaves the
#' system, exactly as it leaves the information.
#'
#' The objective is evaluated once outside the guard, so that a term that
#' cannot be evaluated at all raises where it can be read; inside the search
#' a non-finite value is a statement about the point, a filter whose loadings
#' put it outside the region where its recursion is bounded, and the search
#' must step back from it rather than abandon the run.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param obj The objective, as \code{\link{statmod_objective}} returns it.
#' @param beta The coefficients to start from.
#' @param hyper The hyperparameters.
#' @param optimizer An \code{optimizers7} optimizer, or \code{NULL} for
#'   Newton's method with the exact joint information.
#' @param verbose Logical; report the run.
#'
#' @return A list with \code{par}, \code{value}, \code{converged} and
#'   \code{iterations}.
#'
#' @seealso \code{\link{statmod_fit_structural}}, which fits the term's
#'   parameters alone and is what a term of the likelihood shape still uses.
#'
#' @keywords internal
statmod_fit_joint <- function(spec, design, obj, beta, hyper,
                              optimizer = NULL, verbose = FALSE) {
  sst <- statmod_structural_state(design)
  su <- Filter(function(u) identical(u$kind, "filter"),
               attr(design, "structural"))
  key <- su[[1L]]$term
  nm <- names(sst$zeta[[key]])
  free <- setdiff(nm, sst$held[[key]])
  nb <- length(beta)
  ix <- nb + seq_along(free)

  setz <- function(z) {
    v <- sst$zeta[[key]]
    v[free] <- as.numeric(z)
    sst$zeta[[key]] <- v
    sst$key <- NULL
    sst$value <- NULL
  }
  raw <- function(u) {
    setz(u[ix])
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
      c(obj$gr(b), g)
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
    H + P
  }

  u0 <- c(beta, as.numeric(sst$zeta[[key]][free]))
  raw(u0)
  if (is.null(optimizer)) optimizer <- optimizers7::newton()
  res <- optimizers7::minimize(optimizer, fn, u0, gr = gr, he = he)
  setz(res@par[ix])
  if (verbose) {
    vb_say("%d coefficients and %d parameters of %s: %d iterations, converged %s",
           nb, length(free), short_keys(key), as.integer(res@iterations),
           res@converged)
  }
  list(par = res@par[seq_len(nb)], value = res@value,
       converged = isTRUE(res@converged), iterations = res@iterations)
}





#' Run the Structural Terms at the Current Parameters
#'
#' @description
#' The level each filter adds to the equation it sits in, together with the
#' derivative of that level in the term's own parameters.
#'
#' @details
#' The result is memoized on the coefficients and the term's parameters,
#' since the objective, its gradient and its curvature are asked for at the
#' same point in turn and a filter is the expensive part of each.
#'
#' @param spec A \code{\link{StatmodSpec}}.
#' @param design The design.
#' @param eta_static The static predictors, one per distribution parameter.
#' @param theta_static The parameters they give.
#'
#' @return A named list, one entry per structural term.
#'
#' @seealso \code{\link{statmod_structural}}
#'
#' @keywords internal
statmod_filter_at <- function(spec, design, eta_static, theta_static) {
  su <- attr(design, "structural")
  st <- statmod_structural_state(design)
  out <- list()
  for (u in su) {
    if (!identical(u$kind, "filter")) next
    tm <- spec@terms[[u$param]][[u$term]]
    cb <- structural_callbacks(spec, theta_static, u$param)
    psi <- structural_psi(tm, st$zeta[[u$term]])
    f <- modelterms7::term_filter(tm, eta_static[[u$param]], spec@response,
                                  cb$score, cb$curvature, psi,
                                  fast = cb$fast, threads = spec@threads)
    out[[u$term]] <- list(param = u$param, term = u$term, tm = tm,
                          eta = f$eta, jacobian = f$jacobian,
                          psi = psi, cb = cb,
                          eta_static = eta_static[[u$param]])
  }
  out
}


#' A Structural Term's Parameters, With Standard Errors
#'
#' @description
#' The estimated parameters of every structural term, on their own scale,
#' with standard errors and intervals taken from the joint observed
#' information.
#'
#' @details
#' A structural term contributes no design block, so nothing in the
#' coefficient tables can report it and the values were reachable only
#' through \code{fit@structural}. The information for them exists:
#' \code{\link{statmod_full_information}} spans the coefficients AND the
#' term's own parameters, and the tail block of its inverse is their
#' variance on the unconstrained scale.
#'
#' The interval is built on that scale and mapped back through the link,
#' which is what keeps a persistence inside \eqn{(-1, 1)} and a positive
#' loading positive; the standard error on the parameter scale is the delta
#' method, \eqn{|h'(\zeta)|\,\mathrm{se}(\zeta)}, and is reported for
#' reading rather than for building the interval from. A level held because
#' an intercept in the same equation carries it is not estimated and is
#' reported with a missing standard error rather than a zero one.
#'
#' @param fit A \code{\link{StatmodFit}}.
#' @param level The interval level.
#'
#' @return A data frame with one row per parameter of per structural term,
#'   or \code{NULL} where the model carries none.
#'
#' @seealso \code{\link{statmod_full_information}}
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
    solve(I)
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
    # partial autocorrelation, and the literature's beta_j is the
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
    # a held parameter is not estimated, and a quantity that reads one is not
    # either: it would be reported with the variance of the rest alone
    touches_held <- apply(J[, held, drop = FALSE] != 0, 1L, any)
    for (i in seq_along(rd$name)) {
      s <- if (ok && !touches_held[[i]]) sqrt(max(vq[[i]], 0)) else NA_real_
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
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = u$param, term = u$term, name = rd$name[[i]],
        estimate = as.numeric(val), se = s, lower = lo, upper = hi,
        held = nm[[i]] %in% held, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}


#' The Posterior Break-Points of a Marginal Term
#'
#' @description
#' The posterior mean and standard deviation of each group's latent
#' break-points in a fitted model carrying a marginal break-point term
#' (\code{\link[modelterms7]{jump}}, \code{\link[modelterms7]{seg}} or
#' \code{\link[modelterms7]{jseg}} with \code{marginal = TRUE}).
#'
#' @details
#' The quantities come from the same decomposition the marginal likelihood
#' is computed on: the posterior over a group's intervals or quadrature
#' nodes, with the within-interval moments those of the fitted prior
#' truncated to it. The computation is
#' \code{\link[modelterms7]{term_latent}}'s; this function supplies what
#' the term cannot see, the fitted predictors and the model's log-density.
#'
#' @param fit A \code{\link{StatmodFit}} whose model carries a structural
#'   term of the likelihood shape.
#'
#' @return A data frame with one row per group and break-point:
#'   \code{group}, \code{psi}, \code{mean} and \code{sd}.
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
#' @seealso \code{\link[modelterms7]{term_latent}},
#'   \code{\link{statmod}}
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
