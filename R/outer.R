#' @include inference.R
NULL

# The hyperparameters, estimated by a marginal criterion.
#
# The penalized objective is F(b; t) = -l(b) + sum_j rho_j(b; t), and every
# penalties7 penalty keeps its normalizing constant, so rho IS minus a log
# prior density and nothing has to be added by hand. The marginal likelihood
# is then
#
#   L(t) = integral exp(-F(b; t)) db
#
# and its Laplace approximation at the penalized mode b(t) is
#
#   log L(t) = l(b) - rho(b; t) + (q/2) log(2 pi) - (1/2) log|A'(H + S)A|
#
# with H the information of the likelihood, S the penalty's Hessian, and A an
# orthonormal basis of whichever subspace is integrated over. Writing it out
# with the constants in place reproduces Wood's (2011) REML criterion term for
# term, the -(M_p/2) log(2 pi) of that paper appearing here as the difference
# between the rank of the penalty and the dimension of the coefficient vector.

#' How the Hyperparameters Are Estimated
#'
#' @description
#' The class every outer criterion is an instance of: [reml()], [ml()],
#' [aic()], [bic()] and [cv()] all return one. It carries which subspace of
#' the coefficients is integrated over, which information matrix enters the
#' determinant, and the settings a path or a cross-validation needs.
#'
#' One class serves all five, so a property a given criterion does not read
#' is present and unused. Use the five constructors, which validate; this raw
#' one does not.
#'
#' @param kind `"reml"`, `"ml"`, `"aic"`, `"bic"` or `"cv"`. Decides which
#'   criterion is evaluated and, through [outer_minimize()], whether it is
#'   made small or large.
#' @param hessian `"expected"` or `"observed"`. Which information enters the
#'   Laplace determinant, and whether the exact gradient is available.
#' @param k The price of one effective degree of freedom, for a
#'   prediction-error criterion. `NA_real_` for [bic()], which resolves it to
#'   \eqn{\log n} at fit time, and for the criteria that do not charge per
#'   degree of freedom.
#' @param nfolds How many folds [cv()] uses. Read by `"cv"` alone.
#' @param rule `"min"` or `"1se"`. Read by `"cv"` alone.
#' @param folds A fold number per observation, or `integer(0)` for folds
#'   drawn at fit time. Read by `"cv"` alone.
#'
#' @return An object of class `OuterMethod` with one property per argument
#'   above.
#'
#' @seealso [reml()] and [ml()] for the marginal criteria, [aic()] and
#'   [bic()] for the prediction-error ones, [cv()] for cross-validation,
#'   [statmod()] for where one is passed.
#'
#' @examples
#' reml()
#' ml(hessian = "observed")
#' aic()
#' cv(nfolds = 5)
#'
#' @name OuterMethod-class
#' @aliases OuterMethod
#' @keywords internal
#' @export
OuterMethod <- S7::new_class("OuterMethod",
  properties = list(
    kind = S7::class_character,
    hessian = S7::class_character,
    k = S7::class_numeric,
    nfolds = S7::class_numeric,
    rule = S7::class_character,
    folds = S7::class_numeric
  ),
  validator = function(self) {
    if (!identical(length(self@kind), 1L) ||
        !self@kind %in% c("reml", "ml", "aic", "bic", "cv")) {
      return(paste0("Property 'kind' must be \"reml\", \"ml\", \"aic\", ",
                    "\"bic\" or \"cv\"."))
    }
    if (!identical(length(self@hessian), 1L) ||
        !self@hessian %in% c("expected", "observed")) {
      return("Property 'hessian' must be \"expected\" or \"observed\".")
    }
    if (!identical(length(self@rule), 1L) ||
        !self@rule %in% c("min", "1se")) {
      return("Property 'rule' must be \"min\" or \"1se\".")
    }
    if (length(self@nfolds) != 1L || self@nfolds < 2) {
      return("Property 'nfolds' must be a single number, at least 2.")
    }
    NULL
  }
)


#' The Properties Every Criterion Carries
#'
#' @description
#' The defaults for the three properties only [cv()] reads, so that one class
#' can serve every criterion and each constructor names only what it needs.
#'
#' @details
#' What a **path** does is not among them. How many values it visits, how far
#' down it reaches, and whether a term's own hyperparameters are combined or
#' swept one at a time all belong to the term. The same criterion object is
#' put to the model's smooth hyperparameters as well, and those are read at
#' the mode instead of being swept, so a path setting on the criterion would
#' be meaningless for most of what it is asked about.
#'
#' @return A named list with `nfolds`, `rule` and `folds`, to be spliced into
#'   an [OuterMethod()] call.
#'
#' @seealso [path_fallbacks()] for the settings a term that declares a kinked
#'   penalty without offering a grid argument falls back to,
#'   [modelterms7::term_search()] for where a path's shape is set.
#'
#' @keywords internal
outer_path_defaults <- function() {
  list(nfolds = 10, rule = "min", folds = numeric(0))
}


#' Estimate the Hyperparameters by a Marginal Likelihood
#'
#' @description
#' Build the object that tells [statmod()] to choose a model's smooth
#' hyperparameters by a marginal likelihood: the coefficients are integrated
#' out by a Laplace approximation at the penalized mode, and what is left is
#' maximized in the hyperparameters. Pass the result as `outer_criterion`.
#'
#' `reml()` integrates **every** coefficient out. `ml()` integrates only the
#' penalized directions and profiles the rest. `reml()` is [statmod()]'s
#' default.
#'
#' @details
#' # The criterion
#'
#' At the penalized mode \eqn{\hat\beta(\theta)},
#' \deqn{\log L(\theta) = \ell(\hat\beta) - \rho(\hat\beta;\theta)
#'   + \frac{q}{2}\log 2\pi - \frac12\log|A'(H+S)A|,}
#' with \eqn{H} the information of the log-likelihood, \eqn{S} the penalty's
#' second derivative in the coefficients, and \eqn{A} an orthonormal basis of
#' the subspace integrated over. Nothing is added to \eqn{\rho} to make this
#' work: a \pkg{penalties7} penalty keeps its normalizing constant, so it is
#' exactly minus a log prior density, and for a quadratic penalty that constant
#' carries the \eqn{-\frac{r}{2}\log\lambda} and the log pseudo-determinant
#' that a marginal criterion needs. Written out, the expression reproduces
#' Wood's (2011) REML criterion term for term.
#'
#' # What each one integrates
#'
#' `reml()` takes \eqn{A = I}: every coefficient is integrated, the
#' unpenalized ones under the flat prior their absence of a penalty amounts
#' to.
#'
#' `ml()` takes \eqn{A} spanning the range space of the penalty, so an
#' unpenalized coefficient is profiled. An ordinary covariate is one; so is
#' the linear component of a Demmler-Reinsch smooth, which that smooth's
#' penalty leaves alone.
#'
#' This is the distinction between REML and ML for a variance component in a
#' mixed model, and `reml()` is the default for the same reason: profiling a
#' fixed effect leaves the variance estimate biased downwards.
#'
#' # Which hyperparameters
#'
#' Those of the terms fitted in one system, meaning those whose penalty is
#' twice differentiable. A lasso, a SCAD or an MCP has a kink, its
#' coefficients are estimated by a method of their own, and a Laplace
#' approximation at a point where the second derivative does not exist would
#' be arithmetic with no meaning. Those hyperparameters stay where their term
#' left them, and [bic()] is what [statmod()] puts to them instead.
#'
#' # The exact gradient
#'
#' The criterion has one where the information is the observed one and every
#' penalty under estimation has a Hessian linear in its hyperparameters,
#' which covers [modelterms7::s()], [modelterms7::te()] and any
#' [penalties7::quadratic_penalty()]. It is then supplied to the search and
#' [optimizers7::lbfgs()] becomes the default optimizer; otherwise the search
#' compares values.
#'
#' Measured in evaluations of the criterion, each a whole inner fit, against
#' [optimizers7::nelder_mead()]: 40 against 32 with one smoothing parameter,
#' 40 against 135 with two, 41 against 269 with three, and 12 against 283
#' with three and a modeled scale. It does not pay in one dimension and pays
#' from two on, a simplex needing a vertex per dimension and a quasi-Newton
#' method not.
#'
#' # ML needs a null basis
#'
#' For every penalty that has one, since that is what says which directions
#' are profiled. [penalties7::is_proper()] answers for a penalty with no null
#' space at all, and [penalties7::penalty_null_basis()] for the quadratic and
#' structured branches. A penalty offering neither is rejected by name,
#' instead of being integrated over a subspace guessed at.
#'
#' @param hessian Which information enters the determinant: `"observed"` (the
#'   default) or `"expected"`. Matched with [match.arg()]. The observed one
#'   is what the exact gradient needs.
#'
#' @return An [OuterMethod()] object of kind `"reml"` or `"ml"`, with
#'   `hessian` as supplied and the path settings unused.
#'
#' @references
#' Wood, S. N. (2011). Fast stable restricted maximum likelihood and marginal
#' likelihood estimation of semiparametric generalized linear models.
#' *Journal of the Royal Statistical Society, Series B*, 73(1), 3--36.
#'
#' @seealso [aic()], [bic()] and [cv()] for the prediction-error criteria,
#'   [statmod()] for where these are passed, [hyper()] for reading back what
#'   they chose.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(200, -2, 2))
#' dd$y <- sin(1.4 * dd$x) + rnorm(200, sd = 0.3)
#'
#' fit <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dd,
#'                outer_criterion = reml())
#'
#' # The smoothing parameter was estimated, and hyper() says by what.
#' hyper(fit)
#'
#' # ML profiles the unpenalized directions instead of integrating them, so
#' # it shrinks a little less. The gap is small here because only two of the
#' # ten coefficients are unpenalized; it widens with the fixed effects.
#' fml <- statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dd,
#'                outer_criterion = ml())
#' c(reml = hyper(fit)$estimate, ml = hyper(fml)$estimate)
#' c(reml = sum(fit@edf$edf), ml = sum(fml@edf$edf))
#'
#' # The marginal log-likelihood is what these maximize, and it is only
#' # available where one of them ran.
#' logLik(fit, type = "marginal")
#'
#' @export
reml <- function(hessian = c("observed", "expected")) {
  do.call(OuterMethod, c(list(kind = "reml", hessian = match.arg(hessian),
                             k = NA_real_), outer_path_defaults()))
}

#' @rdname reml
#' @export
ml <- function(hessian = c("observed", "expected")) {
  do.call(OuterMethod, c(list(kind = "ml", hessian = match.arg(hessian),
                             k = NA_real_), outer_path_defaults()))
}


#' @title Print an Outer Method
#' @name print.OuterMethod
#' @description
#' Prints one line naming the criterion and the information it is built on,
#' as `<REML>  observed information`, and for [cv()] the number of folds and
#' the selection rule.
#' @param x An [OuterMethod()].
#' @param ... Unused.
#' @return `x`, invisibly, as a print method should.
#' @seealso [reml()], [aic()] and [cv()] for the objects printed.
#' @keywords internal
print.OuterMethod <- function(x, ...) {
  if (identical(x@kind, "cv")) {
    cat(sprintf("<CV>  %s folds, rule \"%s\"\n",
                if (length(x@folds)) "given" else format(x@nfolds), x@rule))
    return(invisible(x))
  }
  cat(sprintf("<%s>  %s information\n", toupper(x@kind), x@hessian))
  invisible(x)
}
S7::method(print, OuterMethod) <- print.OuterMethod


#' The Hyperparameters an Outer Method Estimates
#'
#' @description
#' One row per hyperparameter of a penalized term that is fitted in the joint
#' system, with the link that carries it onto the whole line.
#'
#' @details
#' A term whose penalty has a kink is left out. Its coefficients are
#' estimated by a coordinate descent with everything else held fixed, and the
#' criterion is a Laplace approximation, which asks for a second derivative
#' that does not exist at a kink. Those hyperparameters go to a path instead.
#'
#' A hyperparameter its own term holds is left out too, whatever kind of
#' penalty it belongs to: [statmod_held()] is what says so.
#'
#' The result is the index the whole outer search runs on. Its row count is
#' the dimension of that search, and zero rows is what makes the criterion
#' not run at all.
#'
#' @param spec A [StatmodSpec()].
#' @param blocks The block split, as [statmod_blocks()] returns it.
#'
#' @return A data frame with one row per estimated hyperparameter and three
#'   character columns, `parameter`, `term` and `name`. The links that carry
#'   each onto the whole line ride on the `"link"` attribute, a list the same
#'   length as the frame has rows, since a list column would not survive the
#'   subsetting the search does.
#'
#' @seealso [hyper_to_eta()] for the map those links define,
#'   [statmod_held()] for what is excluded.
#'
#' @keywords internal
outer_hyper_index <- function(spec, blocks) {
  kinked <- vapply(blocks$sparse, function(b)
    paste(b$param, b$term, sep = "\r"), character(1))
  rows <- list()
  links <- list()
  for (u in statmod_penalty_keys(spec)) {
    if (paste(u$param, u$key, sep = "\r") %in% kinked) next
    for (h in u$penalty@params) {
      # HELD by the term, so there is nothing here to estimate. The term is
      # where the penalty is named and where the caller says so.
      if (h %in% names(u$fixed)) next
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = u$param, term = u$key, name = h,
        stringsAsFactors = FALSE)
      links[[length(links) + 1L]] <- u$penalty@link_params[[h]]
    }
  }
  if (!length(rows)) {
    return(structure(data.frame(parameter = character(0), term = character(0),
                                name = character(0),
                                stringsAsFactors = FALSE),
                     links = list()))
  }
  structure(do.call(rbind, rows), links = links)
}


#' Move Between the Hyperparameters and the Free Vector
#'
#' @description
#' `hyper_to_eta()` carries the estimated hyperparameters onto the whole
#' line through their links; `eta_to_hyper()` puts a free vector back.
#'
#' @details
#' The outer search runs on the free scale for the reason every other search
#' in the toolkit does: a smoothing parameter is positive, and an optimizer
#' that does not know it will step outside the domain. Each hyperparameter's
#' own link comes from the `"link"` attribute [outer_hyper_index()] carries.
#'
#' The two are inverses over the estimated hyperparameters:
#' `eta_to_hyper(hyper_to_eta(h, idx), idx, h)` returns `h`. The
#' hyperparameters `idx` does not name are untouched by both.
#'
#' @param hyper The hyperparameter structure, keyed by distribution parameter
#'   and then by term.
#' @param idx The index, as [outer_hyper_index()] returns it, whose rows say
#'   which hyperparameters are on the search and whose `"link"` attribute
#'   says how each is carried.
#' @param eta A free vector, as long as `idx` has rows.
#'
#' @return `hyper_to_eta()` gives an unnamed numeric vector with one entry
#'   per row of `idx`. `eta_to_hyper()` gives the hyperparameter structure,
#'   with the estimated entries replaced and every other one as it was.
#'
#' @seealso [outer_hyper_index()] for the index and its links,
#'   [outer_fit()] for the search that runs on this scale.
#'
#' @keywords internal
hyper_to_eta <- function(hyper, idx) {
  links <- attr(idx, "links")
  vapply(seq_len(nrow(idx)), function(i) {
    v <- hyper[[idx$parameter[i]]][[idx$term[i]]][[idx$name[i]]]
    linkfunctions7::linkfun(links[[i]], v)
  }, numeric(1))
}

#' @rdname hyper_to_eta
#' @keywords internal
eta_to_hyper <- function(eta, idx, hyper) {
  links <- attr(idx, "links")
  for (i in seq_len(nrow(idx))) {
    hyper[[idx$parameter[i]]][[idx$term[i]]][[idx$name[i]]] <-
      linkfunctions7::linkinv(links[[i]], eta[[i]])
  }
  hyper
}


#' The Subspace a Marginal Criterion Integrates Over
#'
#' @description
#' `NULL` for REML, which integrates everything, and an orthonormal basis
#' of the penalty's range space for ML.
#'
#' @details
#' The total penalty is block diagonal over the terms, each owning its own
#' columns, so its null space is the direct sum of the unpenalized columns
#' and each penalty's own null space. The basis is assembled that way.
#'
#' It is never read off the eigenvalues of the assembled matrix. Such a count
#' measures the arithmetic and not the family: with two smoothing parameters
#' ten orders of magnitude apart the small contributions sink below any
#' tolerance and the count falls, while the null space itself does not move.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param kind `"reml"` or `"ml"`.
#'
#' @return `NULL` for `"reml"`, which integrates every direction and needs no
#'   basis. For `"ml"`, a `p x r` matrix with orthonormal columns spanning
#'   the penalty's range space, `p` the total coefficient count and `r` its
#'   rank.
#'
#' @seealso [penalty_range_basis()] for one penalty's contribution,
#'   [reml()] and [ml()] for the distinction this implements.
#'
#' @keywords internal
integrated_basis <- function(spec, design, kind) {
  if (identical(kind, "reml")) return(NULL)
  params <- spec@distrib@params
  npar <- vapply(design, function(d) d$npar, integer(1))
  offs <- cumsum(npar) - npar
  total <- sum(npar)
  cols <- list()
  for (u in statmod_penalized(spec, design)) {
    # a penalty over a structural term's own parameters owns no column of the
    # design; its directions are appended by statmod_marginal_full()
    if (isTRUE(u$structural)) next
    pen <- u$penalty
    k <- length(u$cols)
    R <- penalty_range_basis(pen, k, u$param, u$key)
    if (!ncol(R)) next
    M <- matrix(0, total, ncol(R))
    M[u$index, ] <- R
    cols[[length(cols) + 1L]] <- M
  }
  if (!length(cols) && structural_penalized(spec, design)) {
    return(matrix(0, total, 0))
  }
  if (!length(cols)) {
    stop(paste0("ml() found no penalized direction to integrate over.\n",
                "  Every penalty here is entirely null, which is not a model",
                " ml() can\n  read; reml() integrates the whole coefficient",
                " vector instead."), call. = FALSE)
  }
  do.call(cbind, cols)
}


#' An Orthonormal Basis of a Penalty's Range Space
#'
#' @description
#' Returns an orthonormal basis of the directions a penalty constrains, which
#' are the ones [ml()] integrates over and [reml()] would integrate anyway.
#'
#' @details
#' A proper penalty constrains every direction, so the basis is the identity
#' and [penalties7::is_proper()] is all that has to be asked. A penalty with
#' a null space supplies it through
#' [penalties7::penalty_null_basis()], and the range basis is the orthogonal
#' complement.
#'
#' A penalty answering neither question is **rejected by name**. Integrating
#' over a subspace guessed at would give a criterion that looks like a
#' number and is not one.
#'
#' @param pen A \pkg{penalties7} penalty object.
#' @param k The number of coefficients in the term's block.
#' @param p The distribution parameter, named in the error message.
#' @param nm The term's name, named in the error message.
#'
#' @return A `k x r` matrix with orthonormal columns, `r` being the penalty's
#'   rank. The `k x k` identity for a proper penalty.
#'
#' @seealso [integrated_basis()], which assembles these blockwise.
#'
#' @keywords internal
penalty_range_basis <- function(pen, k, p, nm) {
  if (isTRUE(penalties7::is_proper(pen))) return(diag(1, k, k))
  N <- tryCatch(penalties7::penalty_null_basis(pen), error = function(e) NULL)
  if (is.null(N)) {
    stop(sprintf(paste0("ml() cannot read the null space of the penalty on",
                        " '%s' in '%s'.\n  Its class does not expose one, so",
                        " which directions are profiled is\n  not something",
                        " to guess at. Use reml(), which integrates every",
                        "\n  coefficient and needs no such basis."), nm, p),
         call. = FALSE)
  }
  N <- as.matrix(N)
  if (!ncol(N)) return(diag(1, k, k))
  # the range is the orthogonal complement of the null space inside the block
  Q <- qr.Q(qr(N), complete = TRUE)
  Q[, seq.int(ncol(N) + 1L, k), drop = FALSE]
}


#' Does a Structural Term Carry a Penalty of Its Own?
#'
#' @description
#' Reports whether any penalty of the model covers the own parameters of a
#' structural term, as opposed to a block of design columns. The deviations
#' of a score-driven filter over a panel are the case.
#'
#' @details
#' The answer decides which matrix the Laplace determinant is taken over.
#' `TRUE` sends the criterion to [statmod_marginal_full()], which spans the
#' coefficients **and** the term's free parameters; `FALSE` leaves it over
#' the coefficients alone.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design, whose structural state holds the terms.
#'
#' @return A single logical. `FALSE` for a model with no structural term, and
#'   for one whose structural term declares no penalty.
#'
#' @seealso [statmod_marginal_full()] for the matrix it selects,
#'   [statmod_structural_penalty()] for the penalties it looks for.
#'
#' @keywords internal
structural_penalized <- function(spec, design) {
  for (u in statmod_penalized(spec, design)) {
    if (isTRUE(u$structural)) return(TRUE)
  }
  FALSE
}


#' The Penalized Curvature Over the Coefficients and a Filter's Parameters
#'
#' @description
#' The matrix whose determinant the Laplace approximation reads, where a
#' penalty covers a structural term's own parameters instead of a block of
#' design columns.
#'
#' @details
#' # Why the determinant has to span them
#'
#' A marginal criterion integrates the quantities its penalty shrinks. Where
#' the penalty is a structural term's own, as with the deviations of a panel,
#' those quantities are not coefficients.
#'
#' Taken over the coefficients alone the determinant then does not depend on
#' the hyperparameter at all, and the criterion reduces to the penalized
#' likelihood, whose maximum in a shrinkage parameter is at no shrinkage. It
#' does not merely lose accuracy: it answers a different question.
#'
#' # Both pieces already exist
#'
#' [statmod_full_information()] spans the coefficients followed by the term's
#' free parameters, which is the order the joint fit uses, and
#' [statmod_structural_penalty()] gives the penalty's Hessian in those same
#' parameters.
#'
#' The information here is the **observed** one whatever the method asks for.
#' The expected information over a filter's own parameters is not a quantity
#' the toolkit carries, the recursion's state entering the expectation.
#'
#' # Under ML
#'
#' The tail is integrated over the penalized coordinates alone, through the
#' penalty's own range basis, so a parameter the penalty does not cover, such
#' as a population value, is profiled instead of integrated. That is what
#' happens to an unpenalized coefficient too.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients, at the penalized mode.
#' @param hyper The hyperparameters.
#' @param basis The integrated subspace over the coefficients, as
#'   [integrated_basis()] returns it, or `NULL` for [reml()].
#'
#' @return A square symmetric matrix spanning the integrated coefficients
#'   followed by the integrated structural parameters. `NULL` where the
#'   structural term could not be run at the given point, which a caller
#'   reads as an unusable hyperparameter.
#'
#' @seealso [statmod_marginal()], the caller,
#'   [structural_penalized()] for the predicate that selects this route.
#'
#' @keywords internal
statmod_marginal_full <- function(spec, design, coef, hyper, basis = NULL) {
  sst <- statmod_structural_state(design)
  su <- Filter(function(u) identical(u$kind, "filter"),
               attr(design, "structural"))
  if (is.null(sst) || !length(su)) return(NULL)
  key <- su[[1L]]$term
  free <- setdiff(names(sst$zeta[[key]]), sst$held[[key]])
  nb <- sum(vapply(design, function(d) d$npar, integer(1)))
  ix <- nb + seq_along(free)

  K <- tryCatch(statmod_full_information(spec, coef, design),
                error = function(e) NULL)
  if (is.null(K) || nrow(K) != nb + length(free)) return(NULL)
  S <- matrix(0, nrow(K), ncol(K))
  # as_dense() for the same reason the joint objective's own Hessian gives:
  # this matrix spans the coefficients AND the structural term's parameters
  # and is dense by construction, so a sparse block cannot be written into a
  # slice of it. This site is the twin of the one in statmod_fit_joint(), and
  # only that one failed when the penalty's accumulator followed the design --
  # a formula carrying both a filter and a random effect is what reaches it.
  S[seq_len(nb), seq_len(nb)] <-
    as_dense(statmod_penalty_at(spec, coef, hyper, design, "hessian"))
  ps <- structural_penalty_block(spec, design, hyper, length(free))
  if (!is.null(ps)) S[ix, ix] <- ps
  S <- zap_nonfinite(S)
  M <- K + S
  if (is.null(basis)) return(M)
  # the same subspace the gradient projects onto, composed once: two callers
  # building it separately would agree only by accident
  A <- structural_joint_basis(spec, design, key, free, nb, basis)
  crossprod(A, M %*% A)
}


#' Which of a Structural Term's Free Parameters a Penalty Covers
#'
#' @description
#' Returns the positions, among a structural term's free parameters, that
#' some penalty of the model shrinks. Those are the directions [ml()]
#' integrates over; the rest are profiled.
#'
#' @details
#' A term declares its penalties through `modelterms7::term_penalties()`,
#' each entry naming a subset of the term's own parameters. This collects
#' those subsets and maps them onto positions in `free`.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param key The term's key in the specification, its call as written.
#' @param free The term's free parameter names, in the order the term holds
#'   them.
#'
#' @return An integer vector of positions into `free`, sorted and without
#'   repeats. Empty where no penalty covers the term.
#'
#' @seealso [statmod_marginal_full()], the caller.
#'
#' @keywords internal
structural_range_cols <- function(spec, design, key, free) {
  sst <- statmod_structural_state(design)
  nm <- names(sst$zeta[[key]])
  out <- integer(0)
  for (u in statmod_penalized(spec, design)) {
    if (!isTRUE(u$structural) || !identical(u$term, key)) next
    out <- c(out, match(nm[u$cols], free))
  }
  sort(unique(out[!is.na(out)]))
}


#' The Marginal Criterion at Given Coefficients and Hyperparameters
#'
#' @description
#' The Laplace approximation to the log marginal likelihood, evaluated at the
#' penalized mode.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#' @param coef The coefficients, at the penalized mode.
#' @param hyper The hyperparameters.
#' @param method An [OuterMethod()].
#' @param approx The approximation for the expected information.
#' @param basis The integrated subspace, from [integrated_basis()].
#'
#' @return A list with `value`, `loglik`, `penalty`,
#'   `logdet` and `q`, or `NULL` where the determinant does not
#'   exist.
#'
#' @keywords internal
statmod_marginal <- function(spec, design, coef, hyper, method,
                             approx = "bartlett", basis = NULL, ctx = NULL) {
  expected <- identical(method@hessian, "expected")
  ll <- statmod_loglik_at(spec, coef, design)
  rho <- statmod_penalty_at(spec, coef, hyper, design, "value")
  if (structural_penalized(spec, design)) {
    M <- statmod_marginal_full(spec, design, coef, hyper, basis)
    if (is.null(M)) return(NULL)
  } else {
    # ⚠️ Sharing the factorization with ctx_penalized() was MEASURED AND NOT
    # TAKEN. Where the criterion's matrix is the one that context holds -- the
    # observed information, nothing projected away -- the same matrix is
    # factorized here and there, which at p = 503 is 12.4 ms spent twice, and
    # reading the determinant from the context removes one of them. End to
    # end it is worth nothing: 1.01x at p = 503 and 0.92x to 1.04x over the
    # four shapes, because the criterion is evaluated at many points the
    # gradient never reaches (every trial point of a line search) and there
    # is nothing to share at those. The cheap route below is sparse where the
    # matrix is, which serves the criterion-only points as well.
    H <- ctx_information(ctx, spec, design, coef, hyper, expected, approx)
    S <- ctx_penalty(ctx, spec, design, coef, hyper)
    M <- H + S
    if (!is.null(basis)) M <- crossprod(basis, M %*% basis)
  }
  # a determinant that does not exist is a hyperparameter the search should
  # step away from, not an error: at a far-out value the penalized information
  # can lose definiteness while the fit itself is sound.
  #
  # It must be a statement about the MATRIX, though, and reading whether
  # chol() raised is not one. Measured on a hierarchical score-driven panel,
  # K+S reaches a condition number of 8.0e15 while still being positive
  # definite, and there whether the factorization succeeds is decided by
  # rounding: the outer search saw a dozen consecutive points as unavailable
  # while backtracking towards one that had been available a moment before,
  # which is a hole in the domain that is not there. pd_logdet() keeps the
  # cheap route where it is safe and pays for the eigendecomposition only
  # where the condition estimate says the cheap one cannot be trusted.
  # A COORDINATE AT A BOUNDARY IS NOT ONE THE LAPLACE APPROXIMATION INTEGRATES
  # OVER. Its curvature there is NaN, which denies the whole matrix a
  # factorization, and it is pinned to the identity instead: the determinant
  # gains log(1) = 0 and the dimension q drops by one, so the coordinate
  # contributes neither a determinant term nor a 2*pi. Measured on a Student
  # t whose nu had reached double.xmax, the criterion went from unavailable
  # -- which stopped the whole fit at its first evaluation -- to a finite
  # value at a converged mode. See pin_boundary().
  M <- pin_boundary(M)
  held <- attr(M, "held")
  ld <- pd_logdet(M)
  if (!isTRUE(ld$ok)) return(NULL)
  q <- nrow(M) - (if (is.null(held)) 0L else held)
  logdet <- ld$logdet
  list(value = ll - rho + q / 2 * log(2 * pi) - logdet / 2,
       loglik = ll, penalty = rho, logdet = logdet, q = q)
}


#' Estimate the Hyperparameters
#'
#' @description
#' Runs `outer_optimizer` on the marginal criterion, refitting the
#' coefficients at every hyperparameter it tries.
#'
#' @details
#' Each evaluation is a whole inner fit, warm-started from the previous one,
#' which is what makes the search affordable: after the first few
#' hyperparameters the coefficients move very little and the inner loop
#' converges in two or three iterations.
#'
#' **The optimizer is chosen by whether the gradient exists.** Where
#' [statmod_marginal_grad()] applies -- the observed information, and
#' penalties whose Hessian is linear in their hyperparameters -- the criterion
#' is handed its exact derivative and [optimizers7::lbfgs()] is the
#' default; otherwise the search compares values and
#' [optimizers7::nelder_mead()] is. An optimizer given explicitly is
#' used as given, and one that needs a gradient it cannot be given will say so
#' itself.
#'
#' @param spec The specification.
#' @param design The design.
#' @param blocks The block split.
#' @param hyper The starting hyperparameters.
#' @param inner_optimizer How the smooth block is fitted.
#' @param method An [OuterMethod()].
#' @param optimizer An \pkg{optimizers7} optimizer, or `NULL` to let the
#'   availability of the exact gradient decide.
#' @param beta The starting coefficients, stacked.
#' @param approx The approximation for the expected information.
#' @param maxit,tol The alternation's budget and tolerance.
#' @param vb The resolved verbosity.
#'
#' @return A list with `par`, `hyper`, `value`, `criterion`,
#'   `converged`, `history` and the inner results.
#'
#' @seealso [reml()], [statmod()]
#'
#' @keywords internal
outer_fit <- function(spec, design, blocks, hyper, inner_optimizer, method,
                      optimizer, beta, approx, maxit, tol, vb) {
  idx <- outer_hyper_index(spec, blocks)
  if (!nrow(idx)) {
    # statmod() does not reach this: a criterion applies to the smooth
    # penalties and comes into play only where there is one, so a model with
    # none leaves it unused rather than failing. What this guards is the
    # internal contract, for a caller reaching the outer loop directly.
    stop(paste0("there is no hyperparameter for a marginal criterion to\n",
                "  estimate: no term here carries a penalty that is twice\n",
                "  differentiable. A lasso, a scad or an mcp is chosen by a",
                "\n  path over its own values instead -- cv() or bic() --",
                " its penalty\n  having a kink."),
         call. = FALSE)
  }
  expected <- identical(method@hessian, "expected")
  basis <- integrated_basis(spec, design, method@kind)
  labels <- paste(idx$parameter, idx$term, idx$name, sep = "/")
  # what the TRACE prints. The full label names the column of hist_outer,
  # where it has to stay unique and reconstructible; a reader of a running
  # fit needs the term and not its whole specification.
  shown <- paste(idx$parameter, short_keys(idx$term), idx$name, sep = "/")
  pe <- outer_minimize(method)
  basis <- if (pe) NULL else basis
  exact <- outer_gradient_ok(spec, design, idx, method, 1L)
  exact2 <- exact && outer_gradient_ok(spec, design, idx, method, 2L)
  # whether the STOPPING RULE is this package's business too. An optimizer
  # given by name comes with its own rule and keeps it; one chosen here is
  # chosen whole, the rule included, because only this package knows what the
  # criterion it is being pointed at can resolve.
  chose_optimizer <- is.null(optimizer)
  if (chose_optimizer) optimizer <- outer_default_optimizer(exact, exact2)
  # Whether the TRACE prints a gradient. It reports what the search is using,
  # so where the search uses none there is none to report, and asking would
  # compute the quantity the laziness below exists to avoid.
  trace_grad <- exact &&
    "gradient" %in% optimizers7::optimizer_provides(optimizer)

  # the search minimizes; a marginal likelihood is maximized and goes in with
  # its sign turned, a prediction-error criterion as it is. The method says
  # which it is rather than the search assuming.
  sgn <- if (pe) 1 else -1

  state <- new.env(parent = emptyenv())
  state$beta <- beta
  state$inner <- NULL
  state$rows <- list()
  state$evals <- 0L
  state$key <- NULL
  # the worst value seen ON THE SEARCH'S OWN SCALE, which is what an
  # unavailable point is made worse than
  state$worst <- NA_real_
  # What the criterion could resolve at the best-located mode seen so far. It
  # is a RUNNING MINIMUM rather than the latest reading, and the asymmetry of
  # the two failures is what chooses that: a resolution smaller than the truth
  # leaves the search where it was, while one larger stops a healthy search
  # short. A minimum can only shrink, so it can only become more conservative.
  state$resolution <- NA_real_

  # WHICH criterion this search is running, written once. A prediction-error
  # method and a marginal one are read by different functions and only one of
  # them takes the integrated basis, so a second reader that picked for itself
  # would answer for the wrong quantity -- which is what
  # criterion_resolution() did on an aic() fit before this existed.
  criterion_at <- function(cf, hy, par, ctx = NULL) {
    if (pe) {
      statmod_pe(spec, design, cf, hy, method, approx,
                 statmod_active(spec, blocks, par, hy))
    } else {
      statmod_marginal(spec, design, cf, hy, method, approx, basis, ctx)
    }
  }

  # one inner fit serves the value and the gradient: an optimizer asks for
  # both at the same point, and refitting for the second would double the cost
  # of every step
  evaluate <- function(eta) {
    key <- paste(format(eta, digits = 17), collapse = ",")
    if (identical(state$key, key)) return(state$last)
    hy <- eta_to_hyper(eta, idx, hyper)
    # a hyperparameter far enough out takes the model somewhere its own
    # machinery cannot go -- a variance beyond what a density can represent,
    # a curvature that is no longer definite. That is a step the search should
    # take back, not an error. At the STARTING point it is neither: nothing
    # has moved yet, so a failure there is the caller's and is raised.
    m <- NULL
    res <- NULL
    cf <- NULL
    ctx <- NULL
    err <- NULL
    # A step the search takes back must leave NO TRACE, and until now only
    # half of it did. The coefficients are protected -- `state$beta` is
    # written below, after the point is known to be usable -- but a
    # structural term's own parameters live in the design's structural state,
    # which is an ENVIRONMENT the inner fit writes into as it goes, and
    # `statmod_fit_structural()` stores its optimizer's last point whether
    # that optimizer converged or not. So an unavailable point moved the
    # filter's parameters and the next evaluation started from wherever the
    # failure had left them, which ratchets: measured on a panel of twenty
    # groups, the search took one long first step, and then every one of the
    # thirty backtracked points was unavailable, including points beside the
    # start that had converged at the first evaluation.
    sst <- statmod_structural_state(design)
    z_save <- if (is.null(sst)) NULL else sst$zeta
    body <- function() {
      res <<- statmod_alternate(spec, design, blocks, hy, inner_optimizer,
                                state$beta, expected, approx, maxit, tol,
                                vb_inner(vb), hold_refresh = TRUE)
      cf <<- res$obj$split(res$par)
      # Both criteria are read AT THE MODE -- a Laplace approximation there,
      # or a log-likelihood and an edf there -- so where the inner fit did not
      # reach one there is no criterion to compare, and the point is
      # unavailable rather than merely worse. This used to be enforced by
      # accident: a wild hyperparameter made the inner fit RAISE, the catch
      # below turned that into an unavailable point, and the search stepped
      # back. Once the inner fit stopped raising and started returning its
      # last usable point, the criterion at that point looked ordinary and the
      # search walked to a smoothing parameter of 1.8e308.
      # one context per point: the criterion, the gradient and the Hessian
      # read the same information, the same penalty and the same factorization
      # instead of each assembling its own. It is built BEFORE the point is
      # judged, because what judges it reads that same factorization.
      ctx <<- outer_context(spec, design, cf, hy, approx)
      # ⚠️ AND THE FLAG IS NOT WHAT DECIDES, because it answers another
      # question. The flag says whether the inner stopping rule fired;
      # availability asks whether the criterion, a Laplace expansion AT THE
      # MODE, is valid here. Measured on `y ~ s(x) | sigma ~ s(z)`, of 38
      # inner fits in one search 38 are at their mode by the second reading
      # -- 1e-09 to 3e-09 against a limit of 1e-03 -- and FOUR by the first.
      # The other 34 stopped on the objective-stall guard with the objective
      # already fixed to twelve significant digits and the score oscillating
      # between 2.5e-06 and 3.3e-06, just above its absolute tolerance of
      # 1e-06. Read as unavailable, they made this line search backtrack
      # eleven times an iteration and accept a step of 0.0026 where the
      # Newton step is 1.4, so the search moved 0.005 in eta over 38
      # evaluations and stopped 4.0 criterion units short of the optimum its
      # own gradient was correctly pointing at.
      #
      # The rule ADDS points and never removes one -- a converged run stays
      # usable whatever the mode error reads -- so no model that fitted
      # before can stop fitting. Making the flag STRICTER is the other
      # direction, and piano_stabilita.txt 13d measured it and withdrew it:
      # it cost a false negative on a good fit.
      if (!isTRUE(res$converged)) {
        q <- inner_mode_error(ctx, spec, design, cf, hy,
                              tryCatch(res$obj$gr(res$par),
                                       error = function(e) NULL), expected)
        if (!is.finite(q) || q > mode_error_limit()) {
          m <<- NULL
          return(invisible(NULL))
        }
      }
      m <<- criterion_at(cf, hy, res$par, ctx)
      if (is.null(m)) return(invisible(NULL))
      # The derivatives are NOT computed here. Which of them the search reads
      # is the OPTIMIZER's business and it varies: nelder_mead reads neither,
      # lbfgs and bfgs read the gradient and never the Hessian, newton reads
      # both. Computing them because the criterion COULD supply them paid for
      # the dearest quantity in the fit -- the Hessian is 49 to 78 per cent of
      # it -- and threw it away. They are produced on demand instead, which
      # needs no predicate about the optimizer and follows what it actually
      # does rather than what its class declares.
      invisible(NULL)
    }
    if (state$evals == 0L) {
      body()
    } else {
      err <- tryCatch({
        body()
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(err)) m <- NULL
    }
    state$evals <- state$evals + 1L
    if (is.null(m)) {
      # the half that was missing: put the structural parameters back where
      # the last USABLE point left them, as the coefficients already are
      if (!is.null(sst) && !is.null(z_save)) {
        sst$zeta <- z_save
        sst$key <- NULL
        sst$value <- NULL
      }
      if (vb$outer) {
        # the PARAMETER scale, as the accepted line uses: printing eta here
        # made a prior scale look negative
        hv <- tryCatch(paste(signif(vapply(seq_along(labels), function(j)
          hyper_value(hy, idx, j), numeric(1)), 5), collapse = ", "),
          error = function(e) paste(signif(eta, 5), collapse = ", "))
        cat(sprintf("  %6d %16s %12s   %s%s\n", state$evals, "unavailable",
                    "--", hv,
                    if (is.null(err)) "" else paste0("  [", err, "]")))
      }
      # The intent stated above, completed for the branch that does not
      # raise: at the STARTING point an unavailable criterion is the
      # caller's problem and is reported with its cause. Handing the
      # optimizer a non-finite value instead has it stop with "the
      # objective is not finite at the starting value", which names the
      # point and not the reason -- measured on a gamma model whose inner
      # fit CONVERGED to a degenerate parameter, so the penalized
      # information had no Cholesky factor at any hyperparameter and the
      # message pointed nowhere.
      if (state$evals == 1L) {
        cause <- if (!is.null(err)) {
          paste0(": ", err)
        } else if (!is.null(res) && !isTRUE(res$converged)) {
          ": the inner fit did not converge there"
        } else {
          paste0(": the inner fit converged to a point where the ",
                 "penalized information has no Cholesky factor, which is ",
                 "what a degenerated parameter (a coordinate run to a ",
                 "boundary, a direction the data do not identify) leaves ",
                 "behind")
        }
        stop(paste0("the ", toupper(method@kind), " criterion is ",
                    "unavailable at the starting hyperparameters", cause,
                    ". Fit with outer_criterion = NULL and inspect the ",
                    "coefficients."), call. = FALSE)
      }
      nh <- nrow(idx)
      # A value the search moves away from, and it has to be FINITE. An
      # infinite one is worse than useless to a method that differences its
      # own gradient: the probe lands in the unavailable region, the
      # difference comes back non-finite, the direction is meaningless and the
      # line search stops -- measured on a panel whose outer bfgs, having no
      # analytic gradient over a structural penalty, died at its 73rd
      # evaluation with the search still improving. A barrier strictly worse
      # than everything seen keeps the objective's own scale, so a difference
      # taken across it points back into the feasible region and the step is
      # simply rejected.
      w <- if (is.finite(state$worst)) state$worst else 1
      bar <- w + abs(w) + 1
      out <- list(value = sgn * bar, grad = rep(0, nh),
                  hess = if (exact2) diag(if (pe) 1 else -1, nh, nh) else NULL,
                  ok = FALSE)
      state$key <- key
      state$last <- out
      return(out)
    }
    state$beta <- res$par
    state$inner <- res
    state$hyper <- hy
    vfn <- sgn * m$value
    if (is.finite(vfn)) {
      state$worst <- if (is.finite(state$worst)) max(state$worst, vfn) else vfn
    }
    row <- data.frame(evaluation = state$evals, criterion = m$value,
                      loglik = m$loglik, penalty = m$penalty)
    for (j in seq_along(labels)) {
      row[[labels[j]]] <- hyper_value(hy, idx, j)
    }
    state$rows[[length(state$rows) + 1L]] <- row
    if (vb$outer) {
      # the same columns the inner optimizer prints, and the hyperparameters
      # on ONE scale: the PARAMETER scale, which is what a reader of a fit
      # wants and what the returned object carries. The unavailable line used
      # to print the free scale instead, so a prior scale read as negative --
      # it was log(sigma) under a heading that looked like sigma.
      # the trace reports the gradient the search is using, so it asks for one
      # only where the search has one to use: printing it under a
      # derivative-free method would compute the very quantity nobody wants.
      gv <- if (exact && trace_grad) derivs(1L) else NULL
      gmax <- if (!is.null(gv) && all(is.finite(gv)))
        sprintf("%12.4g", max(abs(gv))) else sprintf("%12s", "--")
      cat(sprintf("  %6d %16.6f %s   %s\n", state$evals, m$value, gmax,
                  paste(signif(vapply(seq_along(labels), function(j)
                    hyper_value(hy, idx, j), numeric(1)), 5), collapse = ", ")))
    }
    # `par`, `split` and `score` are carried for criterion_resolution(), which
    # needs the mode the criterion was read at and the score the inner fit
    # stopped short by. The score costs one gradient of the inner objective
    # against a whole inner fit, so it is not worth deferring.
    out <- list(value = m$value, ok = TRUE, ctx = ctx, cf = cf, hy = hy,
                par = res$par, split = res$obj$split,
                score = res$obj$gr(res$par))
    state$key <- key
    state$last <- out
    # RECOMPUTED AT EVERY USABLE POINT rather than once at the start. The mode
    # is least well located at the first evaluation, which is a cold start, so
    # a reading taken there is the worst of the run and governs every step
    # after it. One criterion assembly per evaluation against a whole inner
    # fit is what it costs.
    if (chose_optimizer) {
      r <- criterion_resolution(out, spec, design, method, criterion_at)
      if (is.finite(r) && r > 0) {
        state$resolution <- if (is.finite(state$resolution))
          min(state$resolution, r) else r
      } else if (!is.null(attr(r, "mode_error"))) {
        # REFUSED because the inner fit is not at a mode, and that is worth
        # saying: the search then runs with no resolution at all, which is
        # what it did before one existed, and the reason is a property of the
        # fit rather than of the search. Said ONCE -- it is the same fit
        # every evaluation and a line per evaluation would be noise.
        if (vb$outer && !isTRUE(state$said_mode_error)) {
          state$said_mode_error <- TRUE
          vb_say(paste0("no resolution: the inner fit sits %.3g above its own",
                        " minimum, so the displacement it implies is not a",
                        " correction"), attr(r, "mode_error"), indent = 5L)
        }
      }
    }
    out
  }

  # The derivatives, produced when the search asks and cached beside the point
  # they belong to. `order` is 1 for the gradient and 2 for the Hessian.
  derivs <- function(order) {
    st <- state$last
    nh <- nrow(idx)
    slot <- if (order == 1L) "grad" else "hess"
    if (!is.null(st[[slot]])) return(st[[slot]])
    if (!isTRUE(st$ok)) {
      # an unavailable point keeps the barrier's own answers
      return(if (order == 1L) rep(0, nh) else
        if (exact2) diag(if (pe) 1 else -1, nh, nh) else NULL)
    }
    fallback <- function(o) if (o == 1L) rep(0, nh) else
      diag(if (pe) 1 else -1, nh, nh)
    if (pe) {
      # the prediction-error route computes the gradient on the way to the
      # Hessian, so both are kept rather than the second order being asked for
      # twice
      d <- statmod_pe_derivs(spec, design, st$cf, st$hy, method, idx, order)
      if (is.null(st$grad)) {
        st$grad <- if (is.null(d$grad)) fallback(1L) else d$grad
      }
      if (order >= 2L) st$hess <- if (is.null(d$hess)) fallback(2L) else d$hess
    } else if (order == 1L) {
      v <- statmod_marginal_grad(spec, design, st$cf, st$hy, method, idx,
                                 basis, ctx = st$ctx)
      st$grad <- if (is.null(v)) fallback(1L) else v
    } else {
      v <- statmod_marginal_hess(spec, design, st$cf, st$hy, method, idx,
                                 basis, ctx = st$ctx)
      st$hess <- if (is.null(v)) fallback(2L) else v
    }
    state$last <- st
    st[[slot]]
  }

  fn <- function(eta) sgn * evaluate(eta)$value
  gr <- function(eta) { evaluate(eta); sgn * derivs(1L) }
  he <- function(eta) { evaluate(eta); sgn * derivs(2L) }

  # optimizers7 checks a caller-supplied gradient against a directional
  # difference of the objective, and on a criterion this curved the step it
  # uses is too long: measured on a variance component, it read a rate of
  # 0.996 where the gradient is 0.4786 and numDeriv agrees with the gradient
  # to 1.3e-4. The guard is for a gradient that may not match the objective;
  # here both come from this package and the agreement is what the tests
  # establish, so it is turned off for this search and for nothing else.
  if (exact) {
    old_check <- getOption("optimizers7.check_gradient")
    options(optimizers7.check_gradient = FALSE)
    on.exit(options(optimizers7.check_gradient = old_check), add = TRUE)
  }

  # The legend, printed ONCE. The names are what identify a hyperparameter and
  # they are long -- a term's key is the call that produced it -- so repeating
  # them on every evaluation is what made a trace of 130 evaluations
  # unreadable. Here they are said once and the lines below carry the values.
  if (vb$outer) {
    vb_rule(sprintf("outer search: %s over %d hyperparameter%s",
                    toupper(method@kind), length(shown),
                    if (length(shown) == 1L) "" else "s"),
            vb_name(optimizer), indent = 2L)
    for (j in seq_along(shown)) {
      vb_say("h%-3d %s", j, shown[j], indent = 5L)
    }
    # the column header the evaluation lines below fill, in the shape the
    # inner optimizer prints its own. `criterion` is what the search
    # maximizes and `|grad|max` its exact gradient where there is one; the
    # hyperparameters are on the PARAMETER scale throughout.
    cat(sprintf("  %6s %16s %12s   %s\n", "eval", toupper(method@kind),
                "|grad|max",
                if (length(shown) == 1L) "h" else
                  paste0("h1..h", length(shown))))
  }

  eta0 <- hyper_to_eta(hyper, idx)

  # THE CRITERION HAS A RESOLUTION AND THE STOPPING RULE HAS TO KNOW IT.
  # Every evaluation refits the coefficients from the RUNNING warm start, so
  # the value at one hyperparameter depends on the path taken to it, and a
  # line search cannot verify a decrease smaller than that. Both halves of an
  # optimizer's default rule are then unreachable: `crit_grad()` asks 1e-6 of
  # a gradient that bottoms out where the decrease stops being verifiable, and
  # `crit_rel_obj()` asks 1e-12 of a value of order 1e4, which is an absolute
  # 1e-8. The run reports failure at a point it does not leave -- measured on
  # a t-prior random effect, three consecutive evaluations identical to seven
  # digits with the flag FALSE, while the inner fit at those hyperparameters
  # converges and the log-likelihood agrees to nine digits.
  #
  # The resolution is COMPUTED at this fit rather than declared from the inner
  # tolerance, and the difference is not one of accuracy. Declared as the
  # criterion's scale times the tolerance it cannot be both safe and useful:
  # measured, the resolution divided by that scale ranges over 140 times across
  # shapes at ONE tolerance, so a bound covering a random intercept over 500
  # levels is four orders too large for a gaussian smooth, where it stops the
  # search with the outer gradient still 1.8e-2. `criterion_resolution()`
  # displaces the mode by the error the inner fit's own score implies and reads
  # how far the criterion moves, which costs one assembly and no refit; see
  # there for what it was measured against.
  #
  # It is ADDED to the optimizer's own rule rather than replacing it, so the
  # run can only stop earlier and never later, and it goes to the LINE SEARCH
  # as well, the two covering different exits. The rule is read between
  # iterations; a search that exhausts its backtracking budget never gets
  # there, the loop asking the rule with no previous iterate, where anything
  # reading a CHANGE in the objective returns FALSE by construction. Told what
  # the criterion can resolve, the search returns at the first step whose
  # predicted improvement is below it instead of paying for thirty trials,
  # each of them a whole inner fit.
  if (chose_optimizer) {
    fn(eta0)
    resolution <- state$resolution
    if (is.finite(resolution) && resolution > 0) {
      # the RULE keeps the reading from the starting point, being a property of
      # the optimizer object and settled before the run; only the line search,
      # which is asked again at every iteration, follows the running minimum
      optimizer <- S7::set_props(
        optimizer,
        criterion = optimizers7::crit_any(
          optimizer@criterion, optimizers7::crit_abs_obj(resolution)))
    }
    # ⚠️ THE LINE SEARCH IS GIVEN A TENTH OF IT, and the two reasons both point
    # the same way. What is computed is the spread between evaluations reached
    # from DIFFERENT warm starts, while a line search compares trials taken one
    # after another from nearly the same state, where the criterion is far more
    # reproducible -- repeated at one hyperparameter from a running warm start
    # its spread is exactly zero. And the reading is the more aggressive
    # consumer, ending a whole iteration rather than one comparison. Measured
    # at the full number on a gaussian smooth: the search stops after 3
    # evaluations against 29 and gives up 3.2e-06 of criterion against a
    # resolution of 9.7e-07, so it stops just before it has to; at a tenth it
    # reaches the same optimum as a search that was told nothing.
    #
    # A derivative-free search has no line search, so the property is asked for
    # rather than assumed.
    #
    # ⚠️ AND IT IS GIVEN ITS CLOSURE WHETHER OR NOT THE FIRST READING WAS
    # USABLE, which is the whole difference between the two consumers.
    # The criterion's rule is a number settled before the run, so it can only
    # be built from the reading at the starting point -- a cold start, the
    # worst-located mode of the whole fit. The line search reads a CLOSURE over
    # the running minimum, so it needs no reading now and picks up whatever the
    # search refines later.
    #
    # Nesting it inside the test above threw that away: measured on
    # `seg(x, psi ~ random(~1|id))`, whose first reading is refused because the
    # cold-started mode sits 0.046 above its own minimum while the RUNNING
    # MINIMUM is 3.1e-11, the fit went from 5 evaluations to 18 and from
    # converged to not -- at an identical answer, cor 0.9932 and rmse 0.0674
    # either way. Outside the test the closure serves that fit from the good
    # readings that follow.
    #
    # Where no reading is ever usable the closure returns NA, which
    # optimizers7 reads as no resolution at all -- which is the answer for a
    # fit whose mode is never located.
    if ("line_search" %in% S7::prop_names(optimizer)) {
      optimizer <- S7::set_props(
        optimizer,
        line_search = S7::set_props(
          optimizer@line_search,
          resolution = function() state$resolution / 10))
    }
    # THE BACKTRACKING BUDGET, and it is a COST decision rather than a
    # numerical one, which is why it sits outside the resolution's block and
    # is set whether or not one could be computed.
    #
    # optimizers7 defaults to 30 backtracks, which is right for an objective
    # costing microseconds. Here a trial is a penalized refit, and a line
    # search that is going to fail spends the whole budget discovering it --
    # the traces that started this work show two blocks of 31 identical
    # evaluations, which is that budget burned twice.
    #
    # ⚠️ WHAT IT BUYS IS SMALLER THAN THE EVALUATION COUNT SUGGESTS, and the
    # reason is in evaluate() above: `state$beta` is written only after a
    # point is known usable, so every trial of a line search warm-starts from
    # the last ACCEPTED point, and as the step shrinks the trial begins at
    # very nearly its own answer. Measured, removing 22 of 38 evaluations
    # removed 2.8 s of 30.8 -- 0.13 s each against an average evaluation's
    # 0.81. An evaluation count is not a cost here.
    #
    # ⚠️ AND IT MUST BE MEASURED ON THIS PATH AND NOT ANOTHER. Swept with the
    # optimizer NAMED -- which turns the resolution above off -- a short budget
    # flips the convergence flag from TRUE to FALSE on healthy shapes, and that
    # would have been a bad trade. It does not happen here: measured at 30, 12
    # and 8 backtracks on a smooth, two smooths with a random effect, and a
    # random intercept, all three are unchanged in evaluations, criterion to
    # six decimals, effective degrees of freedom and flag, the resolution rule
    # stopping the search before the budget is ever reached. What moves is the
    # expensive shape: a hierarchical break-point model goes from 31
    # evaluations and 25.6 s to 13 and 20.3, with the criterion 1.3e-04 BETTER
    # and the edf identical.
    #
    # 12 rather than 8 because 12 takes 5.3 s of the 6.6 available and 8 takes
    # the rest at a criterion the named-optimizer sweep shows starting to
    # degrade: against a budget of 30 the criterion given up there is 4.8e-05
    # at 12, 7.8e-04 at 8 and 3.3e-03 at 6, and the last is past the 1.6e-03
    # that criterion can resolve.
    if ("line_search" %in% S7::prop_names(optimizer)) {
      optimizer <- S7::set_props(
        optimizer,
        line_search = S7::set_props(optimizer@line_search,
                                    max_step = outer_backtracks()))
    }
  }

  res <- if (exact2) {
    optimizers7::minimize(optimizer, fn, eta0, gr = gr, he = he)
  } else if (exact) {
    optimizers7::minimize(optimizer, fn, eta0, gr = gr)
  } else {
    optimizers7::minimize(optimizer, fn, eta0)
  }

  # the last evaluation is not necessarily the optimum, so the fit is taken at
  # the reported point rather than at whatever was tried last
  hy <- eta_to_hyper(res@par, idx, hyper)
  inner <- statmod_alternate(spec, design, blocks, hy, inner_optimizer,
                             state$beta, expected, approx, maxit, tol, vb,
                             hold_refresh = TRUE)
  cff <- inner$obj$split(inner$par)
  m <- if (pe) statmod_pe(spec, design, cff, hy, method, approx,
                          statmod_active(spec, blocks, inner$par, hy)) else
    statmod_marginal(spec, design, cff, hy, method, approx, basis)
  list(par = inner$par, hyper = hy, value = inner$value,
       criterion = if (is.null(m)) NA_real_ else m$value,
       converged = res@converged && inner$converged,
       obj = inner$obj, hist_blocks = inner$hist_blocks,
       hist_inner = inner$hist_inner,
       hist_outer = if (length(state$rows)) do.call(rbind, state$rows) else
         NULL,
       iterations = res@iterations, evaluations = state$evals,
       exact_gradient = exact, exact_hessian = exact2,
       # the optimizer that actually searched, which is the caller's where
       # one was given and otherwise the one chosen from what the criterion
       # can supply: a fit says what fitted it rather than leaving a reader
       # to reconstruct the default
       optimizer = optimizer)
}


#' One Hyperparameter of an Index
#'
#' @description
#' The value the `j`-th row of an outer index points at.
#'
#' @param hyper The hyperparameter structure.
#' @param idx The index.
#' @param j A row.
#'
#' @return A single number.
#'
#' @keywords internal
hyper_value <- function(hyper, idx, j) {
  as.numeric(hyper[[idx$parameter[j]]][[idx$term[j]]][[idx$name[j]]])[1L]
}


#' The Verbosity of an Inner Fit Inside the Outer Search
#'
#' @description
#' The same switches with the block trace off, since one line per pass per
#' hyperparameter tried is not a trace anybody reads.
#'
#' @param vb The resolved verbosity.
#'
#' @return A list of switches.
#'
#' @keywords internal
vb_inner <- function(vb) {
  vb$blocks <- vb$blocks && vb$inner
  vb
}


#' The Outer Line Search's Backtracking Budget
#'
#' @description
#' How many backtracks a line search inside [outer_fit()] may take
#' before it gives up, where the optimizer is the one this package chose.
#'
#' @details
#' \pkg{optimizers7} defaults to 30, which suits an objective costing
#' microseconds. Every trial here is a penalized refit, and a line search that
#' is going to fail spends its whole budget finding that out.
#'
#' The saving is smaller than the evaluation count implies. Every trial
#' warm-starts from the last accepted point, so as the step shrinks the refit
#' begins at nearly its own answer: measured, removing 22 evaluations of 38
#' removed 2.8 seconds of 30.8, which is 0.13 seconds each against an average
#' evaluation's 0.81.
#'
#' What it costs is nothing on the shapes where nothing was wrong. Measured at
#' 30, 12 and 8 backtracks on a smooth, on two smooths with a random effect and
#' on a random intercept, all three are unchanged in evaluations, in criterion
#' to six decimals, in effective degrees of freedom and in the convergence
#' flag: the resolution [criterion_resolution()] supplies stops the
#' search before the budget is ever reached. A hierarchical break-point model
#' goes from 31 evaluations and 25.6 seconds to 13 and 20.3, with the criterion
#' 1.3e-04 better and the same degrees of freedom.
#'
#' The value is 12 rather than 8 because 12 takes 5.3 seconds of the 6.6 there
#' are to take, and because the criterion 8 gives up is larger: swept with the
#' optimizer named, so that the resolution does not mask it, the criterion lost
#' against a budget of 30 is 4.8e-05 at 12, 7.8e-04 at 8 and 3.3e-03 at 6, the
#' last being past the 1.6e-03 that criterion can resolve.
#'
#' An optimizer the caller named keeps its own budget, as it keeps its own
#' stopping rule.
#'
#' @return A single number.
#'
#' @seealso [outer_fit()], [criterion_resolution()]
#'
#' @keywords internal
outer_backtracks <- function() 12


#' Which Optimizer the Outer Search Uses When the Caller Names None
#'
#' @description
#' The choice is made from what the criterion can supply: its exact Hessian,
#' its exact gradient, or neither.
#'
#' @param exact Whether the criterion has an exact gradient.
#' @param exact2 Whether it has an exact Hessian as well.
#'
#' @return An \pkg{optimizers7} optimizer.
#'
#' @seealso [outer_fit()], [outer_gradient_ok()]
#'
#' @keywords internal
outer_default_optimizer <- function(exact, exact2) {
  if (exact2) return(optimizers7::newton())
  if (exact) return(optimizers7::lbfgs())
  optimizers7::nelder_mead()
}
