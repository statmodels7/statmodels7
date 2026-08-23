#' @include formula.R
NULL

#' Bind a Model's Term Blocks Side by Side
#'
#' @description
#' The assembled design, sparse when any term's block is.
#'
#' @details
#' A grouping indicator is sparse by construction -- a row belongs to one
#' group, so a random effect over \eqn{m} of them has density \eqn{1/m} --
#' and \pkg{modelterms7} builds it that way. Binding it beside a dense block
#' with `cbind()` does not work: base dispatch reads the sparse block as
#' a vector and reports that the number of items to replace is not a multiple
#' of the replacement length, three frames from anything a caller wrote.
#'
#' The result is sparse whenever ONE block is, which is the right rule
#' because sparsity is a property of the assembled matrix rather than of its
#' pieces: a dense fixed block beside a large indicator leaves a matrix that
#' is still overwhelmingly zero, and its factorization stays sparse under a
#' fill-reducing ordering.
#'
#' @param mats A list of blocks.
#' @param n The number of observations, for the empty case.
#'
#' @return A matrix, or a \pkg{Matrix} sparse matrix.
#'
#' @keywords internal
bind_blocks <- function(mats, n) {
  if (!length(mats)) return(matrix(0, n, 0L))
  if (!any(vapply(mats, isS4, logical(1)))) return(do.call(cbind, mats))
  Reduce(function(a, b) Matrix::cbind2(a, b), mats)
}


#' Is a Design Sparse, and the Zero Matrix to Accumulate It Into
#'
#' @description
#' `design_sparse()` reports whether any equation's block is sparse, and
#' `zero_information()` gives the square zero matrix of the right kind to
#' accumulate the information into.
#'
#' @details
#' The information is assembled one `crossprod` per parameter pair and
#' placed into a square accumulator. With a sparse design each product is
#' sparse, and placing it into a dense accumulator signals that the number of
#' items to replace is not a multiple of the replacement length -- from
#' inside the assembly, naming nothing a caller wrote. The accumulator
#' follows the design instead.
#'
#' @param design The design.
#' @param total The number of stacked coefficients.
#'
#' @return A logical, or a square matrix of zeros.
#'
#' @keywords internal
design_sparse <- function(design) {
  any(vapply(design, function(d) isS4(d$X), logical(1)))
}

#' @rdname design_sparse
#' @keywords internal
as_dense <- function(A) if (isS4(A)) as.matrix(A) else A

#' @rdname design_sparse
#' @keywords internal
as_sparse <- function(A) {
  if (isS4(A)) A else methods::as(methods::as(A, "denseMatrix"), "CsparseMatrix")
}

#' Take the Offsets Out of an Equation
#'
#' @description
#' Splits a one-sided formula into the `offset()` terms it names and
#' whatever is left, which is what the interpreter is given.
#'
#' @details
#' An offset is a column of the linear predictor whose coefficient is known to
#' be one, and `y ~ x + offset(log_n)` is how R has always written it.
#' Without this the term reached `model.matrix` through
#' [modelterms7::linpar()], where `terms()` marks it in the
#' `"offset"` attribute and the design EXCLUDES it: the term contributed
#' no column, no offset and no message, and the model fitted was the one
#' without it. On a count model over person-years that moved the intercept
#' from -7.5 to -0.6, which is the difference between a log rate and a log
#' count.
#'
#' Only a top-level additive term is taken, which is where R recognizes one
#' too; `stats::offset(x)` is recognized beside `offset(x)`. Several
#' are summed, as [stats::glm()] sums them.
#'
#' @param eq A one-sided formula.
#'
#' @return A list with `formula`, the equation with the offsets removed,
#'   and `offsets`, a list of the expressions taken out.
#'
#' @seealso [statmod_spec()], [eval_offsets()]
#'
#' @keywords internal
split_offsets <- function(eq) {
  is_offset <- function(e) {
    if (!is.call(e) || length(e) != 2L) return(FALSE)
    f <- e[[1L]]
    if (is.name(f)) return(identical(as.character(f), "offset"))
    # stats::offset(x) and stats:::offset(x)
    is.call(f) && length(f) == 3L &&
      as.character(f[[1L]]) %in% c("::", ":::") &&
      identical(as.character(f[[3L]]), "offset")
  }
  found <- list()
  strip <- function(e) {
    if (is.call(e) && length(e) == 3L && identical(e[[1L]], quote(`+`))) {
      a <- strip(e[[2L]])
      b <- strip(e[[3L]])
      if (is.null(a)) return(b)
      if (is.null(b)) return(a)
      return(call("+", a, b))
    }
    if (is_offset(e)) {
      found[[length(found) + 1L]] <<- e[[2L]]
      return(NULL)
    }
    e
  }
  kept <- strip(eq[[length(eq)]])
  # an equation that was nothing but an offset keeps its intercept, which is
  # what `y ~ offset(x)` means in R
  if (is.null(kept)) kept <- quote(1)
  out <- eq
  out[[length(out)]] <- kept
  list(formula = out, offsets = found)
}


#' Reject an Offset Buried Inside a Term
#'
#' @description
#' Stops where `offset()` appears anywhere other than as a term of the
#' equation itself.
#'
#' @details
#' [split_offsets()] takes the top-level additive terms, which is
#' where R recognizes an offset. One written inside another term's formula --
#' `ridge(~ z + offset(o))`, `random(~ 1 + offset(o) | g)`,
#' `nl(a ~ 0 + ridge(~ g + offset(o)))` -- reaches that term's own
#' `model.matrix`, where it is dropped exactly as it used to be dropped
#' at the equation level: measured, the fit ran, the block had the columns of
#' the model WITHOUT it, and the intercept came back at 1.33 against the
#' -5.01 the offset gives, a factor of 566.
#'
#' It is REFUSED rather than routed up, and the reason is that the meaning
#' differs by where it sits. In a penalized term's formula an offset would be
#' a contribution to the equation's predictor, which is what writing it at the
#' equation level already says. In a SUBFORMULA it would be a contribution to
#' that parameter's own chart, which is a different quantity. Rather than pick
#' one reading and apply it everywhere, the construct is rejected with the
#' place it belongs named, which costs a user one edit and cannot silently
#' fit the wrong model.
#'
#' @param eq A one-sided formula, with the equation's own offsets already
#'   taken out.
#' @param param The parameter the equation belongs to, for the message.
#'
#' @return `NULL`, invisibly; an error where one is found.
#'
#' @seealso [split_offsets()]
#'
#' @keywords internal
reject_nested_offsets <- function(eq, param) {
  found <- NULL
  walk <- function(e) {
    if (!is.call(e) || !is.null(found)) return(invisible(NULL))
    f <- e[[1L]]
    nm <- if (is.name(f)) as.character(f) else
      if (is.call(f) && length(f) == 3L &&
          as.character(f[[1L]]) %in% c("::", ":::")) as.character(f[[3L]]) else ""
    if (identical(nm, "offset") && length(e) == 2L) {
      found <<- e
      return(invisible(NULL))
    }
    for (i in seq_along(e)[-1L]) {
      ei <- e[[i]]
      # an omitted argument -- the blank in `x[, 1]` -- is an empty symbol,
      # and recursing into it errors. missing() does not answer this: it is
      # about a formal argument of a function, not about a call's slot.
      if (is.symbol(ei) && !nzchar(as.character(ei))) next
      walk(ei)
    }
    invisible(NULL)
  }
  walk(eq[[length(eq)]])
  if (!is.null(found)) {
    stop(sprintf(paste0("'%s' in the equation for '%s' is inside another",
                        " term, where it\n  is not an offset: a term builds",
                        " its own block and drops it, so the\n  model would",
                        " be fitted without it.\n  Write it as a term of the",
                        " equation instead -- '... + %s' -- which is\n  where",
                        " R recognizes an offset and where it enters the",
                        " linear\n  predictor."),
                 deparse(found), param, deparse(found)), call. = FALSE)
  }
  invisible(NULL)
}


#' Evaluate the Offsets a Formula Names
#'
#' @description
#' The offset per distribution parameter, evaluated in the data.
#'
#' @details
#' The expressions are re-evaluated rather than carried as numbers, which is
#' what lets an offset survive prediction: [statmod_respec()] calls
#' this against the new data, where a vector supplied through the
#' `offsets` argument at fitting time has the wrong length and cannot be
#' reused.
#'
#' @param formula The model formula, before the offsets are stripped.
#' @param params The distribution's parameter names.
#' @param data A data frame.
#' @param env The environment the formula carried.
#' @param n The number of observations.
#'
#' @return A named list, one entry per parameter, `NULL` where the
#'   equation names no offset.
#'
#' @seealso [split_offsets()]
#'
#' @keywords internal
eval_offsets <- function(formula, params, data, env, n) {
  out <- stats::setNames(vector("list", length(params)), params)
  eqs <- statmod_equations(formula, params)$equations
  for (p in names(eqs)) {
    ex <- split_offsets(eqs[[p]])$offsets
    if (!length(ex)) next
    total <- numeric(n)
    for (e in ex) {
      v <- tryCatch(eval(e, data, env), error = function(err)
        stop(sprintf("The offset '%s' in the equation for '%s' could not be\n  evaluated: %s",
                     deparse(e), p, conditionMessage(err)), call. = FALSE))
      if (!is.numeric(v) || !length(v) %in% c(1L, n)) {
        stop(sprintf(paste0("The offset '%s' in the equation for '%s' must be",
                            " numeric of\n  length 1 or %d, and it is %s of",
                            " length %d."),
                     deparse(e), p, n, class(v)[1L], length(v)),
             call. = FALSE)
      }
      if (anyNA(v) || !all(is.finite(v))) {
        stop(sprintf(paste0("The offset '%s' in the equation for '%s' is not",
                            " finite everywhere.\n  An offset enters the",
                            " linear predictor as it stands, so a missing or",
                            "\n  infinite value has no fitted answer."),
                     deparse(e), p), call. = FALSE)
      }
      total <- total + rep_len(as.numeric(v), n)
    }
    out[[p]] <- total
  }
  out
}


#' Add Two Sets of Offsets
#'
#' @description
#' Combines the offsets a formula names with those the `offsets` argument
#' supplies, per parameter.
#'
#' @details
#' They are SUMMED where both are given, which is what
#' [stats::glm()] does with a formula offset and an `offset`
#' argument together.
#'
#' @param a,b Two named lists of offsets, either entry possibly `NULL`.
#'
#' @return A named list.
#'
#' @keywords internal
add_offsets <- function(a, b) {
  for (p in names(a)) {
    if (is.null(b[[p]])) next
    a[[p]] <- if (is.null(a[[p]])) b[[p]] else a[[p]] + b[[p]]
  }
  a
}


#' Zero the Non-Finite Entries of a Penalty's Hessian
#'
#' @description
#' A penalty at a hyperparameter far enough out returns non-finite entries,
#' and every consumer of [statmod_penalty_at()] zeroes them before
#' using the matrix.
#'
#' @details
#' Seven places wrote `S[!is.finite(S)] <- 0`, which is correct on a base
#' matrix and a trap on a sparse one: the logical index is a DENSE \eqn{p
#' \times p} matrix, so the storage the accumulator exists to keep is thrown
#' away at the first consumer. On a sparse matrix only the STORED values can
#' be non-finite -- a structural zero is finite by construction -- so the same
#' answer comes from the value slot alone, in \eqn{O(\mathrm{nnz})}.
#'
#' Written once here rather than seven times, which is the discipline this
#' package records for a shape of mistake rather than a line of it.
#'
#' @param S A penalty's Hessian, sparse or dense.
#'
#' @return `S` with its non-finite entries replaced by zero.
#'
#' @seealso [statmod_penalty_at()]
#'
#' @keywords internal
zap_nonfinite <- function(S) {
  if (isS4(S)) {
    bad <- !is.finite(S@x)
    if (any(bad)) S@x[bad] <- 0
    return(S)
  }
  S[!is.finite(S)] <- 0
  S
}

#' @rdname design_sparse
#' @keywords internal
zero_information <- function(design, total) {
  if (design_sparse(design)) {
    Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                         dims = c(total, total))
  } else {
    matrix(0, total, total)
  }
}

#' The Specification of a Model, Before It Is Fitted
#'
#' @description
#' Everything the formula, the data and the distribution produce: one
#' equation per distribution parameter, the terms each equation names built
#' against the data, the response, and the prior weights and offsets.
#'
#' @param formula The model formula, as given.
#' @param distrib The \pkg{distributions7} object.
#' @param equations A named list of one-sided formulas, one per parameter, in
#'   the family's order.
#' @param terms A named list, one entry per parameter, each a named list of
#'   built \pkg{modelterms7} terms.
#' @param response The evaluated left-hand side.
#' @param n_obs The number of observations.
#' @param weights Prior weights, one per observation.
#' @param offsets A named list of offsets, one per parameter or `NULL`.
#' @param intercepts A named logical, whether each equation carried one.
#'
#' @return An object of class `StatmodSpec`.
#'
#' @seealso [statmod_spec()]
#'
#' @examples
#' dd <- data.frame(y = rnorm(10), x = runif(10))
#' S7::S7_inherits(statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd),
#'                 StatmodSpec)
#'
#' @name StatmodSpec-class
#' @aliases StatmodSpec
#' @keywords internal
#' @export
StatmodSpec <- S7::new_class("StatmodSpec",
  properties = list(
    formula = S7::class_any,
    distrib = S7::class_any,
    equations = S7::class_list,
    terms = S7::class_list,
    response = S7::class_any,
    n_obs = S7::class_integer,
    weights = S7::class_numeric,
    offsets = S7::class_list,
    intercepts = S7::class_logical,
    newdata = S7::class_any,
    structural = S7::class_list,
    # what the IMPLICIT linpar was built with. It is kept on the
    # specification because a rebuild has to reproduce it: a fold of cv()
    # that built a dense design where the fit built a sparse one would be
    # fitting a different model's storage, and paying for it.
    linpar = S7::new_property(S7::class_list, default = quote(list())),
    # the thread count statmod(threads =) was given, carried here because
    # the specification is what every consumer of the assembly already
    # holds. The default keeps every other constructor -- a respec, a fold
    # of cv() -- on the sequential path.
    threads = S7::new_property(S7::class_integer, default = 1L),
    # the worker-process count for the independent fits of a
    # cross-validation's folds. A fold's own spec is built fresh and takes
    # the defaults, which is what keeps the two levels from nesting: a fit
    # inside a worker is sequential by construction.
    workers = S7::new_property(S7::class_integer, default = 1L)
  )
)


#' Options for the Unpenalized Parametric Block
#'
#' @description
#' How the design of the linear predictor's parametric part is built: the
#' storage, and the contrasts for its factors.
#'
#' @details
#' It governs the IMPLICIT [modelterms7::linpar()] term -- the one
#' the bare covariates of a formula collapse into, which a caller never
#' writes -- so this is the only place its arguments can be given. A
#' `linpar()` written out takes them directly.
#'
#' **Sparse storage.** `sparse = TRUE` builds the block through
#' [Matrix::sparse.model.matrix()], which BUILDS it sparse rather
#' than building a dense matrix and compressing it -- the second would cost
#' the memory the choice exists to avoid. Measured at 20000 rows and a factor
#' of 1000 levels, 0.002 s and 1.8 MB against `stats::model.matrix`'s
#' 0.100 s and 161.5 MB, the numbers identical; and a design that would be
#' 32 GB dense builds in 0.02 s and 19 MB, which is what says there is no
#' dense intermediate. It pays where the formula carries a factor of many
#' levels and costs more than it saves on numeric covariates, whose block is
#' dense whatever is asked for.
#'
#' **There is no rescaling here, and that is measured rather than
#' omitted.** Scaling the columns and carrying the coefficients back is the
#' remedy for a conditioning that squares, which is what forming \eqn{X'X}
#' does; [iwls()] fits through a QR of the design and never forms
#' it. On columns spanning fifteen decades the raw fit and the scaled one
#' converge in the same number of iterations, and both agree with
#' [stats::lm()] to \eqn{10^{-14}}. What does move is the SCORE the
#' fit reports, 1.5e+02 against 9.2e-05, and that is a reading rather than an
#' answer: the final verdict is already arbitrated on a dimensionless scale.
#'
#' @param sparse Whether the block is a `dgCMatrix`. `NULL`, the
#'   default, leaves it to [modelterms7::linpar()], which settles it
#'   at build from the size of the design.
#' @param contrasts The contrasts for the block's factors, as a named list of
#'   the kind [stats::model.matrix()]'s `contrasts.arg` takes,
#'   or `NULL` for the session's `options("contrasts")`.
#'
#' @return A named list, for [statmod()]'s `linpar_control`.
#'   The argument and this function are named differently on purpose: with
#'   one name for both, the argument's default would resolve to its own
#'   promise. [stats::glm()] and [stats::glm.control()]
#'   keep them apart for the same reason.
#'
#' @seealso [statmod()], [modelterms7::linpar()]
#'
#' @examples
#' linpar_options(sparse = TRUE)
#'
#' @export
linpar_options <- function(sparse = NULL, contrasts = NULL) {
  if (!is.null(sparse) &&
      (!is.logical(sparse) || length(sparse) != 1L || is.na(sparse))) {
    stop(paste0("'sparse' must be TRUE, FALSE, or NULL to settle it from",
                " the design."), call. = FALSE)
  }
  if (!is.null(contrasts) && !is.list(contrasts)) {
    stop(paste0("'contrasts' must be a named list, one entry per factor,",
                " or NULL."), call. = FALSE)
  }
  out <- list(sparse = sparse)
  if (!is.null(contrasts)) out$contrasts <- contrasts
  out
}


#' Build a Model Specification
#'
#' @description
#' Splits the formula into one equation per distribution parameter, interprets
#' each with \pkg{modelterms7} and builds its terms against the data.
#'
#' @details
#' The equations are interpreted in an environment where \pkg{modelterms7}'s
#' term constructors shadow whatever the user has attached, so that `s()`
#' means ours even with \pkg{mgcv} on the search path. See
#' [statmod_equations()] for the split itself, which is not the
#' obvious one.
#'
#' Prior weights enter the log-likelihood as \eqn{\sum_i w_i \ell_i} and are
#' taken as given. They are deliberately NOT normalized: dividing by their sum
#' would turn the log-likelihood into a mean, shrinking every standard error
#' by \eqn{\sqrt{n}} and making the information criteria incomparable with an
#' unweighted fit of the same model.
#'
#' @param formula The model formula.
#' @param distrib A \pkg{distributions7} distribution object.
#' @param data A data frame.
#' @param weights Optional prior weights, one per observation.
#' @param offsets Optional named list of offsets, one per parameter.
#' @param need_response Whether the left-hand side must evaluate. A likelihood
#'   needs it; a prediction does not, and new data routinely has no response
#'   column.
#' @param linpar How the IMPLICIT parametric block is built, as
#'   [linpar_options()] returns it. It is kept on the
#'   specification, so a rebuild -- a fold of [cv()] -- reproduces
#'   the storage rather than quietly densifying.
#'
#' @return An object of class [StatmodSpec()].
#'
#' @seealso [statmod_equations()], [statmod()]
#'
#' @examples
#' dd <- data.frame(y = rnorm(20), x = runif(20), z = runif(20))
#' spec <- statmod_spec(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd)
#' names(spec@terms)
#'
#' @export
statmod_spec <- function(formula, distrib, data, weights = NULL,
                         offsets = NULL, need_response = TRUE,
                         linpar = list()) {
  if (!is.data.frame(data)) {
    # rstatmod() returns the truth beside the data rather than attached to
    # it, so a caller who passed the whole result is one field away
    if (inherits(data, "StatmodSim")) {
      stop("'data' is an rstatmod() result, which holds the truth beside ",
           "the data.
  Pass its 'data' field.", call. = FALSE)
    }
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!S7::S7_inherits(distrib, distributions7::distrib)) {
    stop("'distrib' must be a distributions7 distribution object.",
         call. = FALSE)
  }
  params <- distrib@params
  split <- statmod_equations(formula, params)

  # the response is evaluated once, in the data, exactly as
  # interpret_formula() would evaluate it for a single-parameter model
  env <- environment(formula)
  if (is.null(env)) env <- baseenv()
  # prediction needs the design and not the response, and new data routinely
  # has no response column; a likelihood needs it and says so
  response <- tryCatch(eval(split$response, data, env),
                       error = function(e) if (need_response) stop(e) else NULL)
  if (is.null(response)) {
    response <- rep(NA_real_, nrow(data))
  }
  # A CENSORED RESPONSE IS REFUSED HERE, where it can be named. The pieces
  # of a censored likelihood exist -- `cens()` marks the statuses and
  # distributions7 carries the derivatives of a distribution function --
  # and nothing in this package assembles them, so the log-density was
  # being asked for at an object it cannot read and the run died inside
  # `dnorm` three frames from the cause.
  if (S7::S7_inherits(response, modelterms7::censored_response)) {
    stop("A censored response is not fitted: `cens()` marks the statuses ",
         "and this package has no censored likelihood to give them to. ",
         "Fit the observed values, or use distributions7's cdf ",
         "derivatives directly.", call. = FALSE)
  }
  n <- if (is.matrix(response)) nrow(response) else length(response)
  if (n == 0L) stop("The response is empty.", call. = FALSE)

  # the offsets come out of the equations before the interpreter sees them:
  # left in, model.matrix() drops them and the model silently loses them
  stripped <- lapply(split$equations, function(eq) split_offsets(eq)$formula)
  # what is left cannot name an offset: it would be inside another term, and
  # that term's own model.matrix would drop it without a word
  for (p in names(stripped)) reject_nested_offsets(stripped[[p]], p)
  from_formula <- eval_offsets(formula, params, data, env, n)

  built <- statmod_terms(stripped, data, env, response, linpar)
  terms_by_param <- built$terms
  intercepts <- built$intercepts

  weights <- check_weights(weights, n)
  offsets <- add_offsets(from_formula, check_offsets(offsets, params, n))

  StatmodSpec(
    formula = formula, distrib = distrib,
    equations = stripped, terms = terms_by_param,
    response = response, n_obs = as.integer(n),
    weights = weights, offsets = offsets, intercepts = intercepts,
    newdata = NULL, linpar = linpar
  )
}


#' The Same Model Read on Other Rows
#'
#' @description
#' A specification carrying the fitted terms and a new data frame, so that
#' every block is reapplied to those rows rather than rebuilt from them.
#'
#' @details
#' A term records how its block was made -- a factor's levels and contrasts, a
#' spline's knots, a basis reparametrization -- and
#' [modelterms7::term_predict()] reapplies that record. Rebuilding
#' instead gives a block of the same shape, multiplying the same coefficients,
#' that means something else: measured on `y ~ s(x, k = 10)` at 200
#' observations, predicting on 40 of the rows the model was fitted to differed
#' from the fitted values there by 0.237, and on the 51 rows with
#' \eqn{|x| < 0.5}, where the rebuilt knots move furthest, by 1.19. The whole
#' data handed back agrees exactly, which is why nothing noticed.
#'
#' @param spec The fitted [StatmodSpec()].
#' @param data The rows to read it on.
#' @param need_response Whether the response has to be there.
#'
#' @return A [StatmodSpec()] whose `newdata` is set.
#'
#' @seealso [statmod_design()]
#'
#' @keywords internal
statmod_respec <- function(spec, data, need_response = TRUE) {
  if (!is.data.frame(data)) {
    # rstatmod() returns the truth beside the data rather than attached to
    # it, so a caller who passed the whole result is one field away
    if (inherits(data, "StatmodSim")) {
      stop("'data' is an rstatmod() result, which holds the truth beside ",
           "the data.
  Pass its 'data' field.", call. = FALSE)
    }
    stop("'data' must be a data frame.", call. = FALSE)
  }
  env <- environment(spec@formula)
  if (is.null(env)) env <- baseenv()
  split <- statmod_equations(spec@formula, spec@distrib@params)
  response <- tryCatch(eval(split$response, data, env),
                       error = function(e) if (need_response) stop(e) else NULL)
  if (is.null(response)) response <- rep(NA_real_, nrow(data))
  # A CENSORED RESPONSE IS REFUSED HERE, where it can be named. The pieces
  # of a censored likelihood exist -- `cens()` marks the statuses and
  # distributions7 carries the derivatives of a distribution function --
  # and nothing in this package assembles them, so the log-density was
  # being asked for at an object it cannot read and the run died inside
  # `dnorm` three frames from the cause.
  if (S7::S7_inherits(response, modelterms7::censored_response)) {
    stop("A censored response is not fitted: `cens()` marks the statuses ",
         "and this package has no censored likelihood to give them to. ",
         "Fit the observed values, or use distributions7's cdf ",
         "derivatives directly.", call. = FALSE)
  }
  n <- if (is.matrix(response)) nrow(response) else length(response)
  if (n == 0L) stop("The response is empty.", call. = FALSE)
  # An offset the FORMULA names is re-evaluated here, which is the whole
  # reason the expressions are not carried as numbers: a vector supplied
  # through the `offsets` argument at fitting time has the length of the
  # fitting data and cannot be reused, so prediction used to drop the offset
  # and return the predictor of a model without one.
  S7::set_props(spec, response = response, n_obs = as.integer(n),
                weights = rep(1, n),
                offsets = eval_offsets(spec@formula, spec@distrib@params,
                                       data, env, n),
                newdata = data)
}


#' Interpret and Build Each Parameter's Terms
#'
#' @description
#' Runs \pkg{modelterms7}'s interpreter on every equation and builds the terms
#' it names against the data.
#'
#' @details
#' The equations are interpreted with \pkg{modelterms7}'s constructors in front
#' of the search path, so that `s()` means ours whatever the user has
#' attached. A factor covariate needs no special handling: the interpreter
#' collects bare covariates into one `linpar()`, whose block comes from
#' `model.matrix` and therefore carries the contrasts.
#'
#' A break-point term whose starting positions the caller did not name has
#' them chosen on a grid rather than left at the interior quantiles of the
#' covariate; see [seg_grid_start()].
#'
#' @param equations A named list of one-sided formulas.
#' @param data A data frame.
#' @param env The environment the original formula carried.
#' @param response The evaluated left-hand side, or `NULL`.
#'
#' @return A list with `terms` (a named list per parameter) and
#'   `intercepts` (a named logical).
#'
#' @keywords internal
statmod_terms <- function(equations, data, env, response = NULL,
                          linpar = list()) {
  shim <- terms_first(env)
  params <- names(equations)
  out_terms <- stats::setNames(vector("list", length(params)), params)
  intercepts <- stats::setNames(logical(length(params)), params)
  for (p in params) {
    eq <- equations[[p]]
    environment(eq) <- shim
    out <- modelterms7::interpret_formula(eq, data, linpar)
    intercepts[[p]] <- out$intercept
    out_terms[[p]] <- lapply(out$terms, function(tm)
      modelterms7::term_build(seg_grid_start(tm, data, response), data))
    names(out_terms[[p]]) <- names(out$terms)
  }
  reject_unfittable(out_terms)
  list(terms = out_terms, intercepts = intercepts)
}


#' Choose a Break-Point Term's Starting Positions on a Grid
#'
#' @description
#' Runs [modelterms7::seg_start()] on a
#' [modelterms7::seg()], [modelterms7::jump()] or
#' [modelterms7::jseg()] term whose starting positions the caller
#' did not name, and returns the specification unchanged for anything else.
#'
#' @details
#' The objective has local optima in the break-points and the iteration
#' converges from within a basin around where it starts, so where a run
#' begins decides what it finds. Measured over eight samples and four
#' starting positions on a joint jump and change of slope, the fraction of
#' runs recovering the break-point is 0 to 0.5 from a single conventional
#' start and 1 from the grid. The term's own default is a conventional
#' start: the interior quantiles of the covariate, which look at the
#' covariate and not at the response.
#'
#' The rule costs `k` linear fits and is exact for a gaussian
#' response, so it places a starting value and does not fit. Two things
#' are therefore not asked of it. It is applied whatever equation the term
#' sits in, the response being what there is to score against even where
#' the term develops a scale; and it is skipped where the response is not
#' plain numbers -- a censored one, or a matrix -- rather than being given
#' a reading of its own.
#'
#' A caller who names `psi` has said where to begin and is left
#' alone, which is also how the grid is turned off.
#'
#' @param tm One term specification.
#' @param data The data frame the term is built against.
#' @param response The evaluated left-hand side, or `NULL`.
#'
#' @return The specification, with `psi` set where the rule applies.
#'
#' @seealso [modelterms7::seg_start()],
#'   [statmod_terms()]
#'
#' @keywords internal
seg_grid_start <- function(tm, data, response) {
  if (!S7::S7_inherits(tm, modelterms7::SegTerm)) return(tm)
  if (!is.null(tm@spec$psi)) return(tm)
  if (is.null(response) || !is.numeric(response) || is.matrix(response) ||
      anyNA(response) || length(response) != nrow(data)) {
    return(tm)
  }
  # a starting rule is allowed to fail: the term then keeps the quantiles
  # it would have used, and the fit reports what it finds from there
  tryCatch(modelterms7::seg_start(tm, data, response),
           error = function(e) tm)
}


#' Reject a Term the Fitting Scheme Does Not Cover
#'
#' @description
#' Signals an error naming any term whose block is not a fixed design, which
#' is what the alternation of [statmod()] assembles.
#'
#' @details
#' One shape is outside that assembly, and it is read off the term rather than
#' from a list of class names, so a term written later is covered without an
#' edit here.
#'
#' A **structural** term rewrites the likelihood instead of contributing a
#' predictor, so it has no design block at all and answers neither
#' `term_matrix()` nor `term_npar()`. Reaching it through the design
#' produced an error naming one of those generics, which says nothing about the
#' cause. [statmod_structural()] routes those, and what remains here
#' is the term class that is structural and implements neither shape of the
#' contract.
#'
#' A term whose block depends on its own coefficients was rejected here too
#' until the alternation learned to refresh one; it is fitted now, by
#' [statmod_design_at()].
#'
#' Every equation is examined before the error is raised, so a model carrying
#' one such term in the mean and another in the scale reports both rather than
#' the first.
#'
#' @param terms The built terms, a named list of named lists, one per
#'   distribution parameter.
#'
#' @return `NULL`, invisibly; called for the error.
#'
#' @seealso [statmod_terms()], [statmod()]
#'
#' @keywords internal
reject_unfittable <- function(terms) {
  found <- character(0)
  for (p in names(terms)) {
    for (nm in names(terms[[p]])) {
      why <- unfittable_reason(terms[[p]][[nm]])
      if (nzchar(why)) {
        found <- c(found, sprintf("'%s' in '%s': %s", nm, p, why))
      }
    }
  }
  reject_incompatible(terms)
  if (length(found) == 0L) return(invisible(NULL))
  stop(if (length(found) == 1L) {
    paste0("statmod() cannot fit ", found)
  } else {
    paste0("statmod() cannot fit these terms:\n",
           paste0("  ", found, collapse = "\n"))
  }, call. = FALSE)
}


#' Which Structural Levels a Linear Intercept Already Carries
#'
#' @description
#' The parameters of the structural terms that must be held rather than
#' estimated, because the equation they sit in already spans the constant
#' they would shift it by.
#'
#' @details
#' A score-driven level and a regime's first level both add a constant to
#' their equation's predictor. With an intercept there too the two are
#' EXACTLY confounded: shifting the intercept by \eqn{c} and the level by
#' \eqn{-c(1 - \sum_j b_j)} leaves every predictor unchanged, since the
#' recursion is affine in the level given the score path, the score path
#' depends on the predictor alone, and the starting level moves by the same
#' \eqn{c}. The likelihood is flat along that direction, and a fit reaches
#' the ridge without failing -- the score is small because the surface is
#' flat, not because it is a maximum.
#'
#' **The linear intercept wins.** Where both are present the term's
#' level is held at zero and the coefficient carries it, which is what makes
#' `y ~ x + gas(...)` an ordinary thing to write. Nothing about the
#' model is lost: what a constant cannot express is the dynamics, or the
#' difference between one regime and another, and those are the parameters
#' that remain free.
#'
#' The question is asked of the **span** of the equation's design and
#' not of a column named `"(Intercept)"`: a factor coded without one,
#' or any set of columns summing to a constant, spans it just as well.
#' Which parameter is the level is the term's own answer, through
#' [modelterms7::term_level_param()].
#'
#' **A developed level asks the same question of a subspace.** With
#' `omega ~ Z gamma` the confounding is no longer with one constant
#' but with whatever `span(Z)` shares with the span of the equation's
#' design. The constant coordinates are the term's own answer, held as
#' above. For the rest, an unpenalized coordinate whose column lies in the
#' equation's span is flagged with a warning rather than held: holding it
#' would change the model where the confounding is not exact (a
#' time-varying shared column is exactly flat only when its lags stay in
#' the development's span, which depends on \eqn{q} and on the column),
#' while a penalized coordinate is identified by its penalty, exactly as a
#' deviation is. Where the direction really is flat, the variance matrix
#' names it.
#'
#' @param spec A [StatmodSpec()].
#' @param design The design.
#'
#' @return A named list, one character vector per structural term.
#'
#' @seealso [reject_incompatible()]
#'
#' @keywords internal
statmod_held_levels <- function(spec, design) {
  su <- attr(design, "structural")
  out <- list()
  if (is.null(su) || !length(su)) return(out)
  n <- spec@n_obs
  for (u in su) {
    out[[u$term]] <- character(0)
    tm <- spec@terms[[u$param]][[u$term]]
    lvl <- modelterms7::term_level_param(tm)
    X <- design[[u$param]]$X
    if (is.null(X) || !ncol(X)) next
    qrX <- qr(X)
    # the constant lies in the span exactly when regressing it on the design
    # leaves nothing behind
    if (length(lvl) && max(abs(qr.resid(qrX, rep(1, n)))) < 1e-8) {
      out[[u$term]] <- lvl
    }
    # the subspace half of the question, for a developed level: an
    # unpenalized coordinate whose column the equation already spans
    Zl <- modelterms7::term_level_design(tm)
    if (is.null(Zl)) next
    pen_idx <- unlist(lapply(modelterms7::term_penalties(tm),
                             function(e) e$index))
    nm_all <- modelterms7::term_params(tm)
    for (k in seq_len(ncol(Zl))) {
      cnm <- colnames(Zl)[k]
      if (cnm %in% out[[u$term]]) next
      if (match(cnm, nm_all) %in% pen_idx) next
      col <- Zl[, k]
      if (max(abs(qr.resid(qrX, col))) < 1e-8 * max(abs(col), 1)) {
        warning(sprintf(
          paste("the development of the level of '%s' carries '%s', whose",
                "column the design of '%s' already spans; the two may be",
                "confounded, and the coordinate carries no penalty to",
                "identify it."),
          u$term, cnm, u$param), call. = FALSE)
      }
    }
  }
  out
}


#' Combinations of Terms That Are Not a Model
#'
#' @description
#' Rejects a formula carrying more than one structural term, whatever
#' equations they sit in.
#'
#' @details
#' Two of them are not a model that the layer could fit and then report. A
#' filter is driven by the score of the log-likelihood at the predictor it
#' has just produced, so two filters in one equation are two levels adding
#' up, with nothing to tell one from the other; in two equations, each is
#' driven by a score that depends on the other's state, and the pair is one
#' recursion written as two, which neither term implements. A term whose
#' contribution is a likelihood mixed over latent states does not report a
#' predictor at all, so it cannot be combined with anything that expects
#' one.
#'
#' The count is over the whole formula rather than per equation because that
#' is the honest boundary: what makes one admissible is that everything else
#' in the model is a predictor it can be driven by.
#'
#' @param terms The built terms, a named list of named lists.
#'
#' @return `NULL`, invisibly; called for the error.
#'
#' @seealso [reject_unfittable()], [statmod_structural()]
#'
#' @keywords internal
reject_incompatible <- function(terms) {
  found <- character(0)
  for (p in names(terms)) {
    for (nm in names(terms[[p]])) {
      if (S7::S7_inherits(terms[[p]][[nm]], modelterms7::structural_term)) {
        found <- c(found, sprintf("'%s' in '%s'", nm, p))
      }
    }
  }
  if (length(found) <= 1L) return(invisible(NULL))
  stop(sprintf(paste0("statmod() fits at most one structural term and this",
                      " model has %d:\n  %s.\n",
                      "  Each rewrites the predictor or the likelihood the",
                      " others would be read at,\n  so the pair is one",
                      " recursion written as two and neither term",
                      " implements it."),
               length(found), paste(found, collapse = ", ")), call. = FALSE)
}


#' Why a Term Is Outside the Fitting Scheme
#'
#' @description
#' Returns the reason a term cannot be assembled as a fixed design block, or
#' the empty string when it can.
#'
#' @param term One built term.
#'
#' @return A single string.
#'
#' @seealso [reject_unfittable()]
#'
#' @keywords internal
unfittable_reason <- function(term) {
  if (S7::S7_inherits(term, modelterms7::structural_term)) {
    kind <- structural_kind(term)
    if (identical(kind, "filter") || identical(kind, "loglik")) return("")
    return(paste("it is a structural term implementing neither shape of the",
                 "contract: a term that rewrites the likelihood has to",
                 "answer either term_filter(), reporting a predictor, or",
                 "term_loglik(), reporting a likelihood."))
  }
  ""
}


#' Does a Term Recompute Its Own Block?
#'
#' @description
#' `TRUE` when the term registers a
#' [modelterms7::term_refresh()] method of its own rather than
#' inheriting the identity registered on `model_term`.
#'
#' @details
#' The owning class of a method is `attr(m, "signature")[[1]]`, and it is
#' compared by name and package rather than by `identical()`: an S7 class
#' re-created from the same definition is not identical to the original, which
#' is what happens whenever a package's code is re-evaluated rather than
#' loaded.
#'
#' @param term One built term.
#'
#' @return A single logical.
#'
#' @seealso [unfittable_reason()]
#'
#' @keywords internal
refreshes_own_block <- function(term) {
  m <- tryCatch(S7::method(modelterms7::term_refresh, S7::S7_class(term)),
                error = function(e) NULL)
  if (is.null(m)) return(FALSE)
  owner <- attr(m, "signature")[[1]]
  base <- modelterms7::model_term
  !(identical(attr(owner, "name"), attr(base, "name")) &&
    identical(attr(owner, "package"), attr(base, "package")))
}


#' Validate Prior Weights
#'
#' @description
#' Returns a vector of prior weights of the right length, defaulting to one
#' per observation.
#'
#' @details
#' They are not normalized. Making them sum to one would turn
#' \eqn{\sum_i w_i \ell_i} into a mean, which is the averaging trap the
#' objective's own scale already has to avoid: every standard error would
#' shrink by \eqn{\sqrt{n}} and the information criteria would stop being
#' comparable.
#'
#' @param weights The weights, or `NULL`.
#' @param n The number of observations.
#'
#' @return A numeric vector of length `n`.
#'
#' @keywords internal
check_weights <- function(weights, n) {
  if (is.null(weights)) return(rep(1, n))
  if (!is.numeric(weights)) {
    stop("'weights' must be numeric.", call. = FALSE)
  }
  if (length(weights) != n) {
    stop(sprintf("'weights' has length %d but there are %d observations.",
                 length(weights), n), call. = FALSE)
  }
  if (anyNA(weights) || any(weights < 0)) {
    stop("'weights' must be non-negative and complete.", call. = FALSE)
  }
  if (all(weights == 0)) {
    stop("'weights' are all zero: there is nothing to fit.", call. = FALSE)
  }
  as.numeric(weights)
}


#' Validate Offsets
#'
#' @description
#' Returns a named list with one offset per parameter, `NULL` where none
#' was given.
#'
#' @param offsets A named list, or `NULL`.
#' @param params The parameter names.
#' @param n The number of observations.
#'
#' @return A named list of length `length(params)`.
#'
#' @keywords internal
check_offsets <- function(offsets, params, n) {
  out <- stats::setNames(vector("list", length(params)), params)
  if (is.null(offsets)) return(out)
  if (!is.list(offsets) || is.null(names(offsets))) {
    stop("'offsets' must be a named list, one entry per parameter.",
         call. = FALSE)
  }
  bad <- setdiff(names(offsets), params)
  if (length(bad)) {
    stop(sprintf(paste0("'%s' is not a parameter of this distribution.\n",
                        "  Its parameters are: %s."),
                 bad[1L], paste(params, collapse = ", ")), call. = FALSE)
  }
  for (p in names(offsets)) {
    o <- offsets[[p]]
    if (!is.numeric(o) || (length(o) != 1L && length(o) != n)) {
      stop(sprintf("The offset for '%s' must be numeric of length 1 or %d.",
                   p, n), call. = FALSE)
    }
    out[[p]] <- rep_len(as.numeric(o), n)
  }
  out
}


#' The Design of a Specification
#'
#' @description
#' Returns, per parameter, the terms' blocks side by side and the names of the
#' coefficients they carry.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A named list with one entry per parameter, each a list with
#'   `X`, `coef_names`, `npar` and `blocks` (the column
#'   range each term occupies).
#'
#' @seealso [statmod_spec()]
#'
#' @examples
#' dd <- data.frame(y = rnorm(20), x = runif(20))
#' d <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
#' vapply(statmod_design(d), function(z) z$npar, integer(1))
#'
#' @export
statmod_design <- function(spec) {
  out <- statmod_design_blocks(spec)
  su <- statmod_structural(spec)
  if (length(su)) {
    sst <- new.env(parent = emptyenv())
    # a specification a fit has been through carries the parameters it
    # arrived at; without them a prediction would rebuild the filter at its
    # starting values, which is the model without the term
    sst$zeta <- stats::setNames(
      lapply(su, function(u) {
        z <- spec@structural[[u$term]]
        # a fresh start may read the data: the marginal break-point term's is
        # the exact per-group profile of the target, without which the fit
        # begins on a surface that is flat in the prior's own parameters
        if (is.null(z)) {
          structural_zeta_start(spec@terms[[u$param]][[u$term]],
                                target = predictor_target(spec, u$param))
        } else z
      }),
      vapply(su, function(u) u$term, character(1)))
    sst$key <- NULL
    sst$value <- NULL
    attr(out, "structural") <- su
    attr(out, "structure") <- sst
    # asked once the blocks exist, since what it compares is the equation's
    # column space against the constant
    sst$held <- statmod_held_levels(spec, out)
    # An UNHELD level starts at the equation's own data-based intercept
    # rather than at zero -- the mirror of "the intercept wins": held, the
    # intercept carries the response's scale and zero is the right start;
    # unheld, which is the volatility spelling `sigma ~ gas(...) - 1`, the
    # level itself must absorb it, and a filter started at zero on a
    # response of another scale reads a score of y^2/sigma^2 at its first
    # step and leaves the representable range before the search's guard
    # can step back (measured: sigma = NaN at y scaled by 100, where the
    # intercept spelling of the same model fits). Only a FRESH start is
    # touched: a specification a fit has been through keeps the parameters
    # it arrived at.
    fresh <- vapply(su, function(u) is.null(spec@structural[[u$term]]),
                    logical(1))
    lvl_all <- lapply(su, function(u)
      modelterms7::term_level_param(spec@terms[[u$param]][[u$term]]))
    need <- vapply(seq_along(su), function(i) {
      fresh[i] && length(lvl_all[[i]]) &&
        !all(lvl_all[[i]] %in% sst$held[[su[[i]]$term]])
    }, logical(1))
    if (any(need)) {
      eta0 <- statmod_intercepts(spec)
      for (i in which(need)) {
        u <- su[[i]]
        v <- eta0[[u$param]]
        if (is.null(v) || !is.finite(v)) next
        tm <- spec@terms[[u$param]][[u$term]]
        lvl <- setdiff(lvl_all[[i]], sst$held[[u$term]])[1L]
        lk <- modelterms7::term_links(tm)[[lvl]]
        z <- tryCatch(linkfunctions7::linkfun(lk, v),
                      error = function(e) NA_real_)
        if (is.finite(z)) sst$zeta[[u$term]][[lvl]] <- z
      }
    }
  }
  rf <- statmod_refreshable(spec)
  if (length(rf)) {
    # whether a term's block is the Jacobian of its contribution decides how
    # the objective may read it: a Jacobian block is recomputed at every
    # trial point and a scoring step on it is Gauss-Newton, while a frozen
    # working linearization (jump, jseg) is read as the fit last committed
    # it and moves only through statmod_commit_refresh() -- the fixed-point
    # iteration it belongs to is not a descent method on the objective, so
    # dragging its read-off through a line search stalls both
    rf <- lapply(rf, function(r) {
      r$frozen <- !isTRUE(modelterms7::term_jacobian_block(
        spec@terms[[r$param]][[r$term]]))
      r
    })
    # the terms a refresh chains from, and the cache of the last point they
    # were refreshed at; an environment because a design is copied by value
    # and the state has to be the same one wherever it is read
    st <- new.env(parent = emptyenv())
    st$terms <- spec@terms
    st$working <- FALSE
    st$key <- NULL
    st$value <- NULL
    attr(out, "refresh") <- rf
    attr(out, "state") <- st
  }
  out
}


#' The Blocks of a Design
#'
#' @description
#' The per-parameter blocks, before any refreshable term is recomputed at
#' coefficients.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A named list, one entry per distribution parameter.
#'
#' @keywords internal
statmod_design_blocks <- function(spec) {
  lapply(spec@terms, function(tms) {
    # a structural term contributes no columns at all
    tms <- tms[!vapply(tms, S7::S7_inherits, logical(1),
                       modelterms7::structural_term)]
    if (!length(tms)) {
      return(list(X = matrix(0, spec@n_obs, 0L), coef_names = character(0),
                  npar = 0L, blocks = list()))
    }
    # on other rows every block is REAPPLIED, never rebuilt: a term records
    # how its block was made and term_predict() replays that record, where
    # building again would give a block of the same shape multiplying the same
    # coefficients and meaning something else
    mats <- if (is.null(spec@newdata)) lapply(tms, modelterms7::term_matrix)
      else lapply(tms, modelterms7::term_predict, newdata = spec@newdata)
    nms <- lapply(tms, modelterms7::term_coef_names)
    widths <- vapply(mats, ncol, integer(1))
    ends <- cumsum(widths)
    starts <- ends - widths + 1L
    list(
      X = bind_blocks(mats, spec@n_obs),
      coef_names = unlist(nms, use.names = FALSE),
      npar = as.integer(sum(widths)),
      blocks = stats::setNames(
        lapply(seq_along(tms), function(j)
          if (widths[j] == 0L) integer(0) else seq.int(starts[j], ends[j])),
        names(tms))
    )
  })
}
