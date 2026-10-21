#' @include formula.R
NULL

#' Bind a Model's Term Blocks Side by Side
#'
#' @description
#' Binds one equation's term blocks side by side into its design, keeping the
#' result sparse whenever any block is.
#'
#' @details
#' A grouping indicator is sparse by construction. A row belongs to one group,
#' so a random effect over \eqn{m} groups has density \eqn{1/m}, and
#' \pkg{modelterms7} builds it that way.
#'
#' Blocks are bound with `cbind()` when all of them are base matrices. When
#' any block is sparse, they are bound with `Matrix::cbind2()`, so that the
#' result is a \pkg{Matrix} object whatever the order of the blocks.
#'
#' The result is sparse when one block is. Sparsity is a property of the
#' assembled matrix and not of its pieces: a dense fixed block beside a large
#' indicator leaves a matrix that is still overwhelmingly zero, and its
#' factorization stays sparse under a fill-reducing ordering.
#'
#' @param mats A list of design blocks, all with the same number of rows,
#'   each dense or sparse. May be empty.
#' @param n The number of observations, used only when `mats` is empty, where
#'   there is no block from which to read a row count.
#'
#' @return An `n x p` matrix with `p` the total width of the blocks: a
#'   \pkg{Matrix} object when any input is sparse, a base matrix otherwise.
#'   For an empty `mats` the result is an `n x 0` base matrix, which is the
#'   design of an equation that carries only a structural term.
#'
#' @seealso [design_sparse()] for the same test applied to a whole design,
#'   [statmod_design()] for the assembly that this serves.
#'
#' @keywords internal
bind_blocks <- function(mats, n) {
  if (!length(mats)) return(matrix(0, n, 0L))
  if (!any(vapply(mats, isS4, logical(1)))) return(do.call(cbind, mats))
  Reduce(function(a, b) Matrix::cbind2(a, b), mats)
}


#' Is a Design Sparse, and the Zero Matrix into Which to Accumulate It
#'
#' @description
#' `design_sparse()` reports whether any equation's block is sparse, and
#' `zero_information()` gives the square zero matrix of the right kind, into
#' which the information is accumulated. `as_dense()` and `as_sparse()` convert
#' a matrix to the other storage.
#'
#' @details
#' The information is assembled one `crossprod` per parameter pair and each
#' product is written into a square accumulator. With a sparse design each
#' product is sparse, and writing a sparse block into a dense accumulator
#' signals that the number of items to replace is not a multiple of the
#' replacement length, from inside the assembly. The accumulator therefore
#' follows the design.
#'
#' @param design The design, a list with one entry per distribution
#'   parameter, and each entry holds its block as `X`.
#' @param total The number of stacked coefficients across the equations,
#'   which is the side of the accumulator.
#' @param A A matrix, dense or sparse.
#'
#' @return `design_sparse()` gives a single logical. `zero_information()`
#'   gives a `total x total` matrix of zeros, a `dgCMatrix` when the design
#'   is sparse and a base matrix otherwise. `as_dense()` gives `A` as a base
#'   matrix, and `as_sparse()` gives `A` as a \pkg{Matrix} object, converting
#'   a base matrix to compressed column storage.
#'
#' @seealso [statmod_information_at()], which accumulates into it,
#'   [zap_nonfinite()] for the other place where the two storages meet.
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
#' Splits a one-sided formula into the `offset()` terms that it names and
#' whatever is left, and passes the second part to the interpreter.
#'
#' @details
#' An offset is a column of the linear predictor whose coefficient is known
#' to be one, and it is written in R as in `y ~ x + offset(log_n)`.
#'
#' The offset has to be taken out before the interpreter reads the equation.
#' Left in, the term would reach `model.matrix` through
#' [modelterms7::linpar()], where `terms()` marks it in the `"offset"`
#' attribute and the design excludes it: the term would contribute no column,
#' no offset and no message, and the model fitted would be the one without it.
#'
#' Only a top-level additive term is taken, which is where R recognizes one
#' too, and `stats::offset(x)` is recognized beside `offset(x)`. Several in
#' one equation are summed, as [stats::glm()] sums them.
#'
#' @param eq A one-sided formula, one equation of the model.
#'
#' @return A list of two:
#'   \describe{
#'     \item{`formula`}{`eq`, from which the offset terms are removed. The
#'       formula is `~ 1` where the equation held nothing else.}
#'     \item{`offsets`}{a list of the expressions taken out, unevaluated.
#'       The list is empty when there were none.}
#'   }
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
#' [split_offsets()] takes the top-level additive terms, which is where R
#' recognizes an offset. One written inside another term's formula, as in
#' `ridge(~ z + offset(o))`, `random(~ 1 + offset(o) | g)` or
#' `nl(a ~ 0 + ridge(~ g + offset(o)))`, reaches that term's own
#' `model.matrix` and is dropped there, as it would be dropped at the
#' equation level if it were left in. The block then carries the columns of
#' the model without the offset.
#'
#' Such a term is rejected and is not moved to the equation, because the
#' meaning of an offset depends on where it is written. In a penalized term's
#' formula an offset would be a contribution to the equation's predictor, and
#' writing it at the equation level already expresses that. In a subformula it
#' would be a contribution to that parameter's own chart, a different
#' quantity. Picking one reading would fit the wrong model in some cases, so
#' the error names the place where the offset belongs.
#'
#' The whole expression tree is walked, so a term written later is covered.
#'
#' @param eq A one-sided formula, with the equation's own offsets already
#'   taken out by [split_offsets()].
#' @param param The parameter to which the equation belongs, named in the
#'   message.
#'
#' @return `NULL`, invisibly, when the walk finds no offset. Signals an error
#'   naming `param` and the offending call otherwise.
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


#' Reject Missing Values in the Data That a Fit Reads
#'
#' @description
#' Signals an error naming the response, or the first variable of the data
#' that an equation names, where it holds a missing value, with the rows. A
#' fit does not drop rows: in a time series or a panel a dropped row changes
#' what the neighboring rows mean, so the choice of which rows to remove is
#' the caller's.
#'
#' @param response The response, as evaluated.
#' @param equations The equations, one formula per parameter.
#' @param data The data frame.
#'
#' @return `NULL`, invisibly, when every value is present.
#'
#' @keywords internal
reject_missing <- function(response, equations, data) {
  rows <- function(v) {
    w <- if (is.matrix(v)) which(rowSums(is.na(v)) > 0) else which(is.na(v))
    paste0(if (length(w) == 1L) "row " else "rows ",
           paste(utils::head(w, 6L), collapse = ", "),
           if (length(w) > 6L) sprintf(", ... (%d in all)", length(w)) else "")
  }
  tail_msg <- paste0("\n  statmod() does not drop rows. Remove them from the ",
                     "data, or replace the missing values.")
  if (anyNA(response)) {
    stop("The response is missing at ", rows(response), ".", tail_msg,
         call. = FALSE)
  }
  for (p in names(equations)) {
    vars <- intersect(all.vars(equations[[p]]), names(data))
    for (v in vars) {
      if (anyNA(data[[v]])) {
        stop(sprintf("The variable '%s' in the equation for '%s' is missing at %s.",
                     v, p, rows(data[[v]])), tail_msg, call. = FALSE)
      }
    }
  }
  invisible(NULL)
}


#' Evaluate the Offsets That a Formula Names
#'
#' @description
#' Returns the offset of each distribution parameter, evaluated in the data.
#'
#' @details
#' The expressions are kept and re-evaluated, not carried as numbers, so that
#' an offset written in the formula applies in prediction.
#' [statmod_respec()] calls this against the new data; a vector supplied
#' through [statmod()]'s `offsets` argument at fitting time has the length of
#' the fitting data and cannot be reused on other rows.
#'
#' Each expression is evaluated in `data` with `env` behind it, so a symbol
#' resolves as a column first and as a variable of the caller's environment
#' second. A result of length one is recycled to `n`, so a single number is a
#' constant offset, and a result of any other length except `n` is an error.
#' A value that is missing or not finite is an error. Several offsets in one
#' equation are summed.
#'
#' @param formula The model formula, before [split_offsets()] has stripped
#'   anything.
#' @param params The distribution's parameter names, in the family's order.
#' @param data A data frame to evaluate in.
#' @param env The environment that the formula carried, the enclosure of the
#'   evaluation.
#' @param n The number of observations.
#'
#' @return A named list with one entry per element of `params`, each a
#'   numeric vector of length `n` or `NULL` where that equation names no
#'   offset.
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
#' Combines the offsets that a formula names with those that the `offsets`
#' argument supplies, per parameter.
#'
#' @details
#' Where both name an offset for the same parameter the two are summed,
#' as [stats::glm()] does with a formula offset and an `offset` argument
#' together. Where only one does, that one is taken.
#'
#' @param a,b Two named lists of offsets keyed by distribution parameter,
#'   either list's entries possibly `NULL`. Either list may itself be `NULL`.
#'
#' @return `a`, in which each entry is summed with the entry of `b` for the
#'   same parameter where both are present, and is replaced by `b`'s where only
#'   `b` names one. A parameter named only in `b` is not added, and the result
#'   is `NULL` when `a` is.
#'
#' @seealso [eval_offsets()] for the formula's own,
#'   [statmod()] for the `offsets` argument.
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
#' The assignment `S[!is.finite(S)] <- 0` is correct on a base matrix. On a
#' sparse matrix the logical index is a dense \eqn{p \times p} matrix, so the
#' sparse storage that the accumulator is meant to preserve would be lost.
#'
#' On a sparse matrix only the stored values can be non-finite, because a
#' structural zero is finite by construction, so the same result comes
#' from the value slot alone in \eqn{O(\mathrm{nnz})} operations.
#'
#' @param S A penalty's Hessian, a square matrix, sparse or dense.
#'
#' @return `S`, in which the non-finite entries are replaced by zero, and
#'   which keeps the storage in which it arrived.
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


#' A Symmetric Matrix From Its Upper Blocks
#'
#' @description
#' Binds a symmetric matrix out of the blocks on and above its block
#' diagonal; the blocks below are the transposes of those above.
#'
#' @details
#' It replaces writing each block into a zero accumulator. With a sparse
#' design that accumulator is a `dgCMatrix`, and a sub-assignment into it
#' rebuilds the compressed columns at every assignment. Binding the blocks
#' once avoids that and does no arithmetic on them, so the entries are the
#' ones that the products gave.
#'
#' @param blocks A list of lists, `blocks[[a]][[b]]` the block of rows `a`
#'   and columns `b` for `b >= a`. Entries below the diagonal are not read.
#' @param sparse Whether the result is a `dgCMatrix`, as
#'   [design_sparse()] determines for the accumulator that it replaces.
#'
#' @return The assembled matrix, without dimnames: a `dgCMatrix` without
#'   stored zeros when `sparse` is `TRUE`, otherwise a base matrix.
#'
#' @seealso [zero_information()], the accumulator that this replaces where
#'   every block is written exactly once.
#'
#' @keywords internal
assemble_blocks <- function(blocks, sparse) {
  K <- length(blocks)
  get <- function(a, b) {
    m <- if (b >= a) blocks[[a]][[b]] else t(blocks[[b]][[a]])
    if (sparse) methods::as(methods::as(m, "CsparseMatrix"), "generalMatrix")
    else as.matrix(m)
  }
  rows <- lapply(seq_len(K), function(a) {
    r <- lapply(seq_len(K), function(b) get(a, b))
    if (sparse) Reduce(Matrix::cbind2, r) else do.call(cbind, r)
  })
  # a zero block written into the accumulator it replaces stored nothing,
  # and the pattern decides CHOLMOD's ordering, so stored zeros are dropped
  out <- if (sparse) Matrix::drop0(Reduce(Matrix::rbind2, rows)) else
    do.call(rbind, rows)
  # the accumulator it replaces carried no names, and the products do
  dimnames(out) <- NULL
  out
}

#' The Specification of a Model, Before It Is Fitted
#'
#' @description
#' The specification holds everything that the formula, the data and the
#' distribution produce before any fitting: one equation per distribution
#' parameter, the terms that each equation names built against the data, the
#' response, and the prior weights and offsets. A fit keeps one of these in
#' its `spec` property, so it is also what `summary()`, `predict()` and the
#' accessors read.
#'
#' Build one with [statmod_spec()], which interprets the formula and
#' validates. This raw constructor takes the pieces already made.
#'
#' @param formula The model formula as given, including the bars.
#' @param distrib The \pkg{distributions7} distribution object.
#' @param equations A named list of one-sided formulas, one per parameter, in
#'   the family's order.
#' @param terms A named list with one entry per parameter, each a named list
#'   of built \pkg{modelterms7} terms keyed by the term's call as written.
#' @param response The evaluated left-hand side.
#' @param n_obs The number of observations, a single integer.
#' @param weights Prior weights, a numeric vector of length `n_obs`.
#' @param offsets A named list of offsets, one entry per parameter, `NULL`
#'   where an equation has none.
#' @param intercepts A named logical, one per parameter: whether that
#'   equation's parametric block carried an intercept.
#' @param newdata The rows that a respecification reads, or `NULL` for the
#'   fitting rows.
#' @param structural A named list of the parameters of the structural term,
#'   keyed by the term; empty until a fit sets them.
#' @param linpar A list, the options of the implicit parametric block as
#'   [linpar_options()] returns them.
#' @param held_coef A named list, one entry per distribution parameter, of
#'   coefficients held at a value, each a named numeric vector keyed by
#'   coefficient name. The list is empty by default.
#' @param threads The number of threads, a single integer.
#' @param workers The number of worker processes for the independent fits of
#'   a cross-validation, a single integer.
#'
#' @return An object of class `StatmodSpec` with one property per argument
#'   above and the properties described next. `linpar` holds the options of
#'   the implicit parametric block and is set by [statmod_spec()]. `newdata`
#'   holds the rows
#'   that a respecification reads, `NULL` until [statmod_respec()] sets it.
#'   `structural` holds the parameters of a structural term where a fit has
#'   set them. `held_coef` holds coefficients held at a value, by name, empty
#'   by default. `threads` and `workers` hold the thread and worker counts
#'   that [statmod()] sets.
#'
#' @seealso [statmod_spec()], the constructor to use,
#'   [statmod_design()] for what is assembled from one.
#'
#' @examples
#' dd <- data.frame(y = rnorm(10), x = runif(10))
#' spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
#' S7::S7_inherits(spec, StatmodSpec)
#'
#' # The specification has one equation per parameter of the family, whatever
#' # the formula wrote.
#' names(spec@equations)
#' spec@n_obs
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
    # COEFFICIENTS HELD AT A VALUE: one entry per distribution parameter, a
    # named numeric vector keyed by coefficient name. A held coordinate is
    # written into the starting values and then dropped from every scoring
    # step's solve, so the rest is fitted CONDITIONALLY on it and the score
    # in the held direction is what the constrained optimum leaves. That
    # score is the whole content of the Rao and Terrell statistics, and the
    # restricted maximum is the whole content of the likelihood ratio; a
    # hold is how all three are computed without rewriting the formula,
    # which would rebuild every term and move a basis's knots.
    held_coef = S7::new_property(S7::class_list, default = quote(list())),
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
#' Specifies how the parametric part of each linear predictor is built: in
#' what storage, and with which contrasts for its factors. Pass the result as
#' [statmod()]'s `linpar_control`.
#'
#' It is meant for a formula naming a factor of many levels, where
#' `sparse = TRUE` keeps the design in sparse storage instead of a dense
#' matrix with one column per level.
#'
#' @details
#' # Which term it reaches
#'
#' This function configures the implicit [modelterms7::linpar()] term, which
#' is the term into which the bare covariates of a formula collapse. A caller
#' never writes this term and so has no other way to configure it. A
#' `linpar()` written out in the formula takes these arguments directly and
#' ignores this.
#'
#' # Sparse storage
#'
#' `sparse = TRUE` builds the block through
#' [Matrix::sparse.model.matrix()], which builds it sparse instead of
#' building a dense matrix and compressing it, so that it forms no dense
#' intermediate. Sparse storage pays where the formula carries a factor of many
#' levels. A numeric covariate has no zeros to exploit, and its sparse block
#' costs more than it saves.
#'
#' # Scaling of the columns
#'
#' The block is not rescaled. Scaling the columns and carrying the
#' coefficients back is a remedy for a conditioning that squares, as forming
#' \eqn{X'X} does. [iwls()] with its default `decomposition = "qr"` fits
#' through a QR factorization of the augmented design and does not form
#' \eqn{X'X}, and its final convergence verdict is read on a dimensionless
#' scale, so the scale of a column does not affect that verdict.
#'
#' # The name of the argument
#'
#' The argument of [statmod()] and this function are named differently on
#' purpose. With one name for both, the argument's default would resolve to
#' its own promise and R would report *promise already under evaluation*.
#' [stats::glm()] and [stats::glm.control()] are kept apart for the same
#' reason.
#'
#' @param sparse A single logical, or `NULL`. `TRUE` makes the block a
#'   `dgCMatrix`, `FALSE` a base matrix. `NULL`, the default, leaves it to
#'   [modelterms7::linpar()], which settles it at build time from the size of
#'   the design.
#' @param contrasts The contrasts for the block's factors, a named list of
#'   the kind that [stats::model.matrix()]'s `contrasts.arg` takes, or `NULL`
#'   for the session's `options("contrasts")`. The contrasts are carried on
#'   the specification, so that a fold of [cv()] reproduces them instead of
#'   re-reading the option.
#'
#' @return A named list with the element `sparse` (possibly `NULL`) and, when
#'   `contrasts` is not `NULL`, the element `contrasts`, to be passed as
#'   [statmod()]'s `linpar_control`.
#'
#' @seealso [statmod()] for where it is passed,
#'   [modelterms7::linpar()] for the term that it configures.
#'
#' @examples
#' linpar_options(sparse = TRUE)
#'
#' # Sparse storage matters for a factor of many levels.
#' set.seed(1)
#' n <- 2000
#' dd <- data.frame(y = rnorm(n), g = factor(sample(200, n, replace = TRUE)))
#' dense  <- statmod_spec(y ~ g, distributions7::gaussian1_distrib(), dd,
#'                        linpar = linpar_options(sparse = FALSE))
#' sparse <- statmod_spec(y ~ g, distributions7::gaussian1_distrib(), dd,
#'                        linpar = linpar_options(sparse = TRUE))
#' Xd <- statmod_design(dense)$mu$X
#' Xs <- statmod_design(sparse)$mu$X
#'
#' # The two blocks hold the same numbers, and the sparse one uses much less
#' # memory.
#' c(dense = class(Xd)[1], sparse = class(Xs)[1])
#' all.equal(as.matrix(Xs), Xd, check.attributes = FALSE)
#' c(dense = object.size(Xd), sparse = object.size(Xs))
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
#' Builds everything from which [statmod()] fits: splits the formula into one
#' equation per distribution parameter, interprets each with
#' \pkg{modelterms7}, and builds the terms that it names against the data.
#' [statmod()] calls this first, and calling it directly inspects a model
#' without fitting it.
#'
#' @details
#' # The formula
#'
#' The equations are separated by `|`, and the first of them carries the
#' response. A parameter with no equation is given an intercept.
#' [statmod_equations()] does the split, which is not the obvious one: R's
#' precedence reduces the whole right-hand side of a three-equation formula to
#' the last term.
#'
#' Each equation is interpreted in an environment where \pkg{modelterms7}'s
#' term constructors sit in front of the search path, so `s()` means this
#' toolkit's even with \pkg{mgcv} attached.
#'
#' # Rejected before any fit
#'
#' The function signals an error for a multivariate family (which
#' [distributions7::fit_distrib()] fits with constant parameters), a censored
#' response written with `cens()`, a family method whose signature has no
#' `...`, a missing value in the response or in a variable that an equation
#' names, a two-column response on a family without `size`, an offset written
#' inside another term, and a structural term that cannot be fitted or a
#' formula with more than one structural term.
#'
#' # Prior weights are not normalized
#'
#' They enter the log-likelihood as \eqn{\sum_i w_i \ell_i} and are taken as
#' given, so they act as frequency weights: a constant weight \eqn{c} on every
#' row is \eqn{c} copies of the data. Dividing the weights by their sum would
#' turn the log-likelihood into a mean, multiplying every standard error by
#' \eqn{\sqrt{n}} and making the information criteria incomparable with an
#' unweighted fit of the same model.
#'
#' The restricted (REML) criterion counts \eqn{c n} observations against the
#' \eqn{p} coefficients when the weights are the constant \eqn{c}. The REML
#' estimate of a dispersion parameter therefore changes under a constant
#' rescaling of the weights, and the ML estimate does not.
#'
#' @param formula The model formula, with the parameters' equations separated
#'   by `|`.
#' @param distrib A \pkg{distributions7} distribution object, which determines
#'   how many equations there are and what they are called.
#' @param data A data frame holding the response, the covariates and any
#'   matrix columns that the terms name. A missing value in the response, or in
#'   a variable that an equation names, is an error, because a fit does not drop
#'   rows.
#' @param weights Optional prior weights, a numeric vector of length
#'   `nrow(data)` (not the name of a column), or `NULL` for all ones. The
#'   weights must be non-negative, finite and complete, and not all zero. They
#'   act as frequency weights and are not normalized.
#' @param offsets Optional named list of offsets, one entry per parameter,
#'   each a finite numeric vector of length 1 or `nrow(data)`, summed with any
#'   offset the formula names for the same parameter.
#' @param need_response `TRUE`, the default, requires the left-hand side to
#'   evaluate. Prediction uses `FALSE`, because new data routinely has no
#'   response column.
#' @param linpar How the implicit parametric block is built, as
#'   [linpar_options()] returns it. The options are kept on the specification,
#'   so that a rebuild such as a fold of [cv()] reproduces the storage instead
#'   of quietly densifying.
#'
#' @return A [StatmodSpec()] object.
#'
#' @seealso [statmod_equations()] for the split,
#'   [statmod_design()] for the assembly, [statmod()] to fit it.
#'
#' @examples
#' dd <- data.frame(y = rnorm(20), x = runif(20), z = runif(20))
#' spec <- statmod_spec(y ~ x | sigma ~ z,
#'                      distributions7::gaussian1_distrib(), dd)
#'
#' # The terms have one entry per parameter of the family, in the family's
#' # order.
#' names(spec@terms)
#' spec@equations
#'
#' # Unweighted, so the weights are ones.
#' c(n = spec@n_obs, total_weight = sum(spec@weights))
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
  # A MULTIVARIATE FAMILY IS REJECTED HERE, BY NAME. Its parameters are one
  # value for the whole sample, so a model that lets them vary by observation
  # is a branch this layer does not have; without this check the fit reached
  # distributions7's own parameter check three layers down, whose message
  # tells the caller to use a model -- which is what they were doing.
  if (S7::S7_inherits(distrib, distributions7::multivariate_distrib)) {
    stop(sprintf(paste0(
      "statmod() does not fit a multivariate response: '%s' describes a ",
      "vector of\n  responses, and this layer models one response per ",
      "observation. The family\n  can be fitted with constant parameters by ",
      "distributions7::fit_distrib()."),
      distrib@distrib_name), call. = FALSE)
  }
  assert_method_signatures(distrib)
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
         "and this package has no censored likelihood that accepts them. ",
         "Fit the observed values, or use distributions7's cdf ",
         "derivatives directly.", call. = FALSE)
  }
  # cbind(successes, failures) writes the trials onto the family, and a
  # two-level factor is how a binary outcome is ordinarily recorded: every
  # distribution wants a number
  diviso <- split_binomial_matrix(response, distrib)
  response <- coerce_response(diviso$response)
  distrib <- diviso$distrib
  n <- if (is.matrix(response)) nrow(response) else length(response)
  if (n == 0L) stop("The response is empty.", call. = FALSE)

  # the offsets come out of the equations before the interpreter sees them:
  # left in, model.matrix() drops them and the model silently loses them
  stripped <- lapply(split$equations, function(eq) split_offsets(eq)$formula)
  # what is left cannot name an offset: it would be inside another term, and
  # that term's own model.matrix would drop it without a word
  for (p in names(stripped)) reject_nested_offsets(stripped[[p]], p)
  if (need_response) reject_missing(response, stripped, data)
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
#' Builds a specification that carries the fitted terms and a new data frame,
#' so that every block is reapplied to those rows rather than rebuilt from
#' them.
#'
#' @details
#' A term records how its block was made: a factor's levels and contrasts, a
#' spline's knots, a basis reparametrization. [modelterms7::term_predict()]
#' reapplies that record to new rows.
#'
#' Rebuilding the block on the new rows instead would give a block of the
#' same shape, multiplying the same coefficients, with a different meaning:
#' the knots of a spline, for example, would move with the rows. A structural
#' term is the exception. Where the rows carry a complete response, it is
#' rebuilt on them with [modelterms7::term_build()], so that its recursion
#' runs over their times and groups from its own start; where they do not, it
#' continues the fitted series.
#'
#' The offsets named in the formula are re-evaluated against `data` rather
#' than carried across, since a vector of the fitting data's length does not
#' apply to other rows.
#'
#' # The rows are the data's
#'
#' `n_obs` is `nrow(data)` and never the length of whatever the response
#' expression evaluates to. `eval(expr, data, env)` falls through to the
#' formula's environment for anything that `data` does not carry, and a
#' prediction grid carries the covariates and not the response. A response
#' expression would then resolve to a variable of that name in scope, such as
#' the fit's own response, and the row count would be that variable's.
#'
#' So a response whose variables are not columns of `data` is absent
#' where `need_response` is `FALSE`, and one of a length other than
#' `nrow(data)` is an error where it is `TRUE`, rather than a vector
#' recycled into a different model.
#'
#' @param spec The fitted [StatmodSpec()].
#' @param data The rows to read the model on.
#' @param need_response `TRUE` where the response must be present, as for a
#'   log-likelihood; `FALSE` for a prediction.
#'
#' @return A [StatmodSpec()] that carries the fitted terms. Its `newdata` is
#'   `data`, its `n_obs` is the row count of `data`, and the response and the
#'   offsets named in the formula are evaluated on `data`. The weights are
#'   reset to ones, the family carries the numbers of trials of these rows
#'   (missing where the rows have no response and the family holds one per
#'   observation), and a structural term is rebuilt on them where they carry a
#'   complete response. Every other property is `spec`'s.
#'
#' @seealso [statmod_design()], which reapplies the terms,
#'   [predict.StatmodFit()], the caller.
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
  n_rows <- nrow(data)
  # THE ROW COUNT IS THE DATA'S AND NEVER THE RESPONSE'S. `eval(expr, data,
  # env)` falls through to the formula's environment for anything `data`
  # does not carry, and a prediction grid carries the covariates and NOT the
  # response -- which is what the page of predict() promises. So the
  # response resolved to whatever variable of that name was in scope, which
  # in an ordinary session is the fit's own response, and `n_obs` came back
  # as the FITTING count: predict(fit, "mu", newdata) returned the fitted
  # values, of the wrong length, in silence. Measured on 37 new rows against
  # 400 fitted: 37 with no such variable in scope, 400 with one.
  wanted <- all.vars(split$response)
  from_here <- !length(wanted) || all(wanted %in% names(data))
  response <- if (!from_here && !need_response) {
    NULL
  } else {
    tryCatch(eval(split$response, data, env),
             error = function(e) if (need_response) stop(e) else NULL)
  }
  # A response of another length is not these rows' either, whatever it
  # resolved from, and a caller who needs one is told rather than fitted to
  # a recycled vector.
  n_resp <- function(r) if (is.matrix(r)) nrow(r) else length(r)
  if (!is.null(response) && n_resp(response) != n_rows) {
    if (need_response) {
      stop(sprintf(paste0(
        "the response has %d value%s and 'data' has %d row%s. A response",
        " that\n  is not a column of 'data' is read from the formula's",
        " environment, which is the source\n  of a length like this."),
        n_resp(response), if (n_resp(response) == 1L) "" else "s",
        n_rows, if (n_rows == 1L) "" else "s"), call. = FALSE)
    }
    response <- NULL
  }
  absent <- is.null(response)
  if (absent) response <- rep(NA_real_, n_rows)
  # A CENSORED RESPONSE IS REFUSED HERE, where it can be named. The pieces
  # of a censored likelihood exist -- `cens()` marks the statuses and
  # distributions7 carries the derivatives of a distribution function --
  # and nothing in this package assembles them, so the log-density was
  # being asked for at an object it cannot read and the run died inside
  # `dnorm` three frames from the cause.
  if (S7::S7_inherits(response, modelterms7::censored_response)) {
    stop("A censored response is not fitted: `cens()` marks the statuses ",
         "and this package has no censored likelihood that accepts them. ",
         "Fit the observed values, or use distributions7's cdf ",
         "derivatives directly.", call. = FALSE)
  }
  # The trials are recomputed from THESE rows. That is the whole reason they
  # are read off the response rather than taken from the family: a `size`
  # fixed at fitting time cannot follow a fold or a prediction, and recycling
  # it against a different row count returns a number rather than an error.
  # The size the fitted family carries is the row sums of the FITTING rows,
  # written there by the same function. A response matrix at these rows
  # brings its own row sums, so it is compared with nothing: measured, a
  # prediction at rows whose trials differ from the fitting rows' stopped on
  # "the number of trials is given twice and the two disagree".
  famiglia <- spec@distrib
  if (is.matrix(response) && "size" %in% S7::prop_names(famiglia) &&
      !identical(as.character(famiglia@dimension), "multivariate")) {
    famiglia <- S7::set_props(famiglia, size = 1)
  }
  diviso <- split_binomial_matrix(response, famiglia)
  # The SAME mapping as the fit, not one read off these rows: a fold or a
  # prediction can hold only one of the two values, and re-reading it there
  # would code the response the other way round.
  response <- coerce_response(diviso$response,
                              levels = attr(spec@response, "response_levels"))
  n <- if (is.matrix(response)) nrow(response) else length(response)
  if (n == 0L) stop("The response is empty.", call. = FALSE)
  distrib <- diviso$distrib
  # A prediction at rows that carry no response has no trials either. The
  # parameters do not depend on them -- the probability of a binomial is a
  # function of the predictor alone -- so the numbers of trials the fit
  # carried are replaced by missing values, which a quantity that does need
  # them (a moment, a prediction interval) refuses through
  # check_trials_known(), and nothing is recycled against the wrong rows.
  if (absent && !need_response && "size" %in% S7::prop_names(distrib) &&
      length(distrib@size) > 1L) {
    distrib <- S7::set_props(distrib, size = rep(NA_real_, n))
  } else {
    distrib <- check_trials(distrib, n)
  }
  # An offset the FORMULA names is re-evaluated here, which is the whole
  # reason the expressions are not carried as numbers: a vector supplied
  # through the `offsets` argument at fitting time has the length of the
  # fitting data and cannot be reused, so prediction used to drop the offset
  # and return the predictor of a model without one.
  out <- S7::set_props(spec, response = response, n_obs = as.integer(n),
                       distrib = distrib,
                       weights = rep(1, n),
                       offsets = eval_offsets(spec@formula, spec@distrib@params,
                                              data, env, n),
                       newdata = data)
  # ROWS THAT CARRY THE RESPONSE ARE A SERIES OF THEIR OWN: a structural term
  # is rebuilt on them, with the fitted parameters, so the recursion runs
  # over their times and groups from its own start, and every route then
  # reads them as it reads the fitting rows. Left with the fitting rows'
  # blueprint, predict(newdata = <the fitting rows>) on a regime() model was
  # rejected as a prediction past the series, and on a gas() model the
  # standard error read the intercept's row alone (46.5 at every year of the
  # Nile, where the fitted values carry 46.5, 34.7, 33.1, ...). Rows without
  # it continue the fitted series (statmod_eta_continued()).
  su <- statmod_structural(spec)
  known <- if (is.matrix(response)) stats::complete.cases(response) else
    !is.na(response)
  if (!absent && length(su) && all(known)) {
    tms <- out@terms
    for (u in su) {
      tms[[u$param]][[u$term]] <-
        modelterms7::term_build(tms[[u$param]][[u$term]], data)
    }
    out@terms <- tms
  }
  out
}


#' Interpret and Build Each Parameter's Terms
#'
#' @description
#' Runs \pkg{modelterms7}'s interpreter on every equation and builds the terms
#' that it names against the data.
#'
#' @details
#' The equations are interpreted with \pkg{modelterms7}'s constructors in
#' front of the search path, so `s()` means this toolkit's whatever the
#' caller has attached. See [terms_first()] for the shim and the collision that
#' it removes.
#'
#' A factor covariate needs no handling of its own: the interpreter collects
#' the bare covariates of an equation into one [modelterms7::linpar()], whose
#' block comes from `model.matrix`, or from `Matrix::sparse.model.matrix()`
#' when the block is sparse, and so carries the contrasts.
#'
#' The starting positions of a break-point term are chosen on a grid, through
#' [seg_grid_start()], when the caller did not name them, in place of the
#' interior quantiles of the covariate, which are the term's default.
#'
#' @param equations A named list of one-sided formulas, one per distribution
#'   parameter.
#' @param data The data frame against which the terms are built.
#' @param env The environment that the original formula carried, which becomes
#'   the parent of the interpreting environment.
#' @param response The evaluated left-hand side, or `NULL` where there is
#'   none. Only [seg_grid_start()] reads it.
#' @param linpar The options of the implicit parametric block, as
#'   [linpar_options()] returns them.
#'
#' @return A list of two:
#'   \describe{
#'     \item{`terms`}{a named list with one entry per parameter, each a named
#'       list of built terms keyed by the term's call as written.}
#'     \item{`intercepts`}{a named logical, one per parameter: whether that
#'       equation's parametric block carried an intercept.}
#'   }
#'
#' @seealso [statmod_spec()], the caller, [terms_first()] for the shim.
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
  reject_duplicate_names(out_terms)
  # the covariance classes are assembled HERE as well as where the penalties
  # are enumerated, so that a label used on two groupings, or a joint prior
  # named twice or of the wrong dimension, is reported where the formula was
  # read rather than at the first evaluation of a criterion
  statmod_classes(out_terms)
  # several gas() terms in different equations are one recursion, and the
  # rest of the layer sees them as one structural term
  out_terms <- join_gas_terms(out_terms)
  list(terms = out_terms, intercepts = intercepts)
}

#' Reject Two Terms of One Equation That Name Their Coefficients Alike
#'
#' @description
#' The coefficients of an equation are addressed by name in `coef()`,
#' `vcov()`, `confint()` and the summary, so two terms of one equation must
#' not give their coefficients the same names. A term prefixes its names
#' with its label, which for `s()` and `te()` is built from the covariates:
#' two smooths of the same covariate with different constructions, such as
#' `s(x, bspline_smooth(k = 10)) + s(x, pspline_smooth())`, would share it.
#' Such a formula is rejected here, naming the terms and the argument
#' `label` that separates them.
#'
#' @param terms The built terms, a list with one element per distribution
#'   parameter, each a named list of terms as [statmod_terms()] builds them.
#'
#' @return `NULL`, invisibly. The function is called for its error.
#'
#' @keywords internal
reject_duplicate_names <- function(terms) {
  for (p in names(terms)) {
    tms <- terms[[p]]
    tms <- tms[!vapply(tms, S7::S7_inherits, logical(1),
                       modelterms7::structural_term)]
    if (length(tms) < 2L) next
    nms <- lapply(tms, modelterms7::term_coef_names)
    owner <- rep(names(tms), lengths(nms))
    all_nm <- unlist(nms, use.names = FALSE)
    dup <- unique(all_nm[duplicated(all_nm)])
    if (!length(dup)) next
    who <- unique(owner[all_nm %in% dup])
    stop(sprintf(paste0(
      "the terms %s of '%s' give their coefficients the same names",
      " (%s).\n  Give one of them its own prefix with the argument",
      " 'label', for example label = \"%s2\"."),
      paste(sprintf("'%s'", who), collapse = " and "), p,
      paste(dup[seq_len(min(3L, length(dup)))], collapse = ", "),
      sub("\\..*$", "", dup[[1L]])), call. = FALSE)
  }
  invisible(NULL)
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
#' The objective has local optima in the break-points, and the iteration
#' converges from within a basin around its starting position, so the start
#' determines which optimum is found. The term's own default is a conventional
#' start, the interior quantiles of the covariate, which depend on the
#' covariate and not on the response. Scoring a grid of candidate positions on
#' the response gives a start in the basin of the best fit of the linear part.
#'
#' [modelterms7::seg_start()] performs one linear fit for each increasing
#' combination of `npsi` points of a grid of ten, \code{choose(10, npsi)} fits
#' in all. The score is exact for a gaussian response, so the rule places a
#' starting value and does not fit. The rule is applied in whichever equation
#' the term sits, because the response is the quantity against which the rule
#' scores, even where the term develops a scale. It is skipped where the
#' response is not plain numbers (a censored one, a matrix or one with
#' missing values), because those have no reading here.
#'
#' A term whose `psi` the caller names keeps it, and naming `psi` also turns
#' off the grid.
#'
#' @param tm One term specification, as the formula interpreter produced it.
#' @param data The data frame against which the term is built.
#' @param response The evaluated left-hand side, or `NULL`.
#'
#' @return `tm`, in which `psi` is set to the grid's choice where the rule
#'   applies, and `tm` unchanged where it does not: a term of another kind, a
#'   term whose `psi` the caller named, or a response against which the rule
#'   cannot score.
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


#' Reject a Term That the Fitting Scheme Does Not Cover
#'
#' @description
#' Signals an error naming any term that the fitting scheme of [statmod()]
#' does not cover, and any formula that carries more than one structural term.
#'
#' @details
#' Terms are examined by what they implement and not by a list of class
#' names, so a term written later is covered without an edit here.
#'
#' A structural term rewrites the likelihood or the filtered predictor
#' instead of contributing a fixed design block, so it has no design block and
#' does not implement `term_matrix()` or `term_npar()`. A term of the filter
#' shape, such as `gas()`, contributes a predictor through `term_filter()`. A
#' term of the likelihood shape, such as `regime()`, contributes a likelihood
#' through `term_loglik()`. [statmod_structural()] routes both. The function
#' rejects a structural term that implements neither `term_filter()` nor
#' `term_loglik()`, and a structural term of the likelihood shape that
#' carries a covariance label ([unfittable_reason()]).
#' [reject_incompatible()] is called before these reasons are raised, so a
#' formula with more than one structural term is rejected with that
#' function's message.
#'
#' A term whose block depends on its own coefficients is fitted by
#' [statmod_design_at()] and is not rejected here.
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
#' Returns the parameters of the structural terms that must be held rather
#' than estimated, because the equation in which they sit already spans the
#' constant by which they would shift it.
#'
#' @details
#' A score-driven level and a regime's first level both add a constant to
#' their equation's predictor. With an intercept there too the two are
#' exactly confounded: shifting the intercept by \eqn{c} and the level by
#' \eqn{-c(1 - \sum_j \xi_j)} leaves every predictor unchanged, since the
#' recursion is affine in the level given the score path, the score path
#' depends on the predictor alone, and the starting level moves by the same
#' \eqn{c}. The likelihood is flat along that direction, and a fit can reach
#' the ridge without failing, because the score is small wherever the surface
#' is flat.
#'
#' The linear intercept takes precedence. Where both are present the term's
#' level is held at zero and the coefficient carries it, which makes
#' `y ~ x + gas(...)` an ordinary formula. This leaves the model able to
#' express the same fits, since a constant cannot express the dynamics or the
#' difference between one regime and another, and the parameters that express
#' them remain free.
#'
#' The function tests the span of the equation's design and not a column
#' named `"(Intercept)"`: a factor coded without one, or any set of
#' columns summing to a constant, spans the constant just as well.
#' [modelterms7::term_level_param()] gives which parameter of the term is
#' the level.
#'
#' For a developed level, the same test applies to a subspace. With
#' `omega ~ Z gamma` the confounding is no longer with one constant but with
#' whatever `span(Z)` shares with the span of the equation's design. The
#' constant coordinates are held as above. For the rest, an unpenalized
#' coordinate whose column lies in the equation's span is flagged with a
#' warning and is not held: holding it would change the model where the
#' confounding is not exact (a time-varying shared column is exactly flat
#' only when its lags stay in the development's span, which depends on
#' \eqn{q} and on the column). A penalized coordinate is identified by its
#' penalty. Where the direction really is flat, the variance matrix shows
#' it.
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
    if (is_joint_filter(tm)) {
      # each level adds a constant to its own equation, so each is held when
      # that equation's design spans the constant
      lv <- modelterms7::term_level_param(tm)
      for (k in seq_along(tm@params)) {
        Xk <- design[[tm@params[k]]]$X
        if (is.null(Xk) || !ncol(Xk)) next
        if (max(abs(qr.resid(qr(Xk), rep(1, n)))) < 1e-8) {
          out[[u$term]] <- c(out[[u$term]], lv[k])
        }
      }
      next
    }
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
    # the development is built on the rows the term was fitted on, so the
    # question is one about the fitting data; at new rows (a forecast past
    # the series) the two designs do not share their rows and there is
    # nothing to compare. Asked there, qr.resid() stopped every forecast of a
    # panel whose level is developed: 308 rows of development against six of
    # newdata on nlme::Ovary.
    if (is.null(Zl) || nrow(Zl) != nrow(X)) next
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
#' Rejects a formula carrying more than one structural term, in whatever
#' equations they sit.
#'
#' @details
#' A filter is driven by the score of the log-likelihood at the predictor
#' that it has just produced, so two filters in one equation are two levels
#' adding up, with nothing to tell one from the other, and are rejected.
#' Score-driven terms in different equations are accepted: each is driven by
#' a score that depends on the others' state, so together they are one
#' recursion, which [join_gas_terms()] assembles. Any other structural term
#' beside a second one is rejected.
#' A term whose contribution is a likelihood mixed over latent states does not
#' report a predictor at all, so it cannot be combined with anything that
#' expects one.
#'
#' The count is over the whole formula and not per equation, because a
#' structural term is admissible only where everything else in the model is a
#' predictor that can drive it.
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
  gas_eq <- character(0)
  other <- 0L
  for (p in names(terms)) {
    for (nm in names(terms[[p]])) {
      tm <- terms[[p]][[nm]]
      if (S7::S7_inherits(tm, modelterms7::structural_term)) {
        found <- c(found, sprintf("'%s' in '%s'", nm, p))
        if (S7::S7_inherits(tm, modelterms7::GasTerm)) {
          gas_eq <- c(gas_eq, p)
        } else {
          other <- other + 1L
        }
      }
    }
  }
  if (length(found) <= 1L) return(invisible(NULL))
  # score-driven terms in DIFFERENT equations, with nothing structural beside
  # them, are filtered as one recursion by join_gas_terms()
  if (other == 0L && !anyDuplicated(gas_eq)) return(invisible(NULL))
  if (other == 0L) {
    p <- gas_eq[anyDuplicated(gas_eq)]
    stop(sprintf(paste0("the equation of '%s' carries two gas() terms:\n",
                        "  %s.\n",
                        "  Two levels in one predictor add up, with nothing to",
                        " tell one from the other.\n  One gas() term in each",
                        " of several equations is filtered jointly."),
                 p, paste(found[gas_eq == p], collapse = ", ")),
         call. = FALSE)
  }
  stop(sprintf(paste0("statmod() fits one structural term, or one gas() term",
                      " in each of several\n  equations, and this model has",
                      " %d:\n  %s.\n",
                      "  Each rewrites the predictor or the likelihood at",
                      " which the others would be read,\n  and a structural",
                      " term other than gas() beside another one is not",
                      " implemented."),
               length(found), paste(found, collapse = ", ")), call. = FALSE)
}


#' The Covariance Labels That a Term Carries, Including Its Sub-Terms
#'
#' @description
#' Returns every covariance label reachable from a term: its own, and those
#' of the sub-terms that develop its own parameters, found by walking to any
#' depth.
#'
#' @details
#' A label written inside a subformula is not visible where the equation's
#' terms are enumerated: `seg(x, psi ~ random(~ 1 | u | id))` is one
#' `SegTerm`, whose own [modelterms7::term_tag()] is `NA`, and the labeled
#' term is a sub-term of its break-point. The sub-terms are reached through
#' [modelterms7::term_components()], whose `subs` field is the list of them,
#' and the walk recurses because a sub-term is an ordinary term and may
#' develop parameters of its own.
#'
#' The names of the result record where each label was found, so a message can
#' name the sub-term rather than only the equation's term.
#'
#' @param term One built term.
#' @param what How to name this term in the result's names. Defaults to the
#'   term itself, and a recursive call passes the path by which the term was
#'   reached.
#'
#' @return A named character vector, empty when the term has no reachable
#'   label.
#'
#' @seealso [unfittable_reason()], its caller;
#'   [modelterms7::term_tag()] for one term's answer.
#'
#' @keywords internal
term_tags_deep <- function(term, what = "it") {
  out <- character(0)
  tag <- tryCatch(modelterms7::term_tag(term), error = function(e) NA_character_)
  if (length(tag) == 1L && !is.na(tag)) {
    out <- stats::setNames(as.character(tag), what)
  }
  comp <- tryCatch(modelterms7::term_components(term), error = function(e) list())
  for (nm in names(comp)) {
    subs <- comp[[nm]]$subs
    if (!length(subs)) next
    for (s in subs) {
      out <- c(out, term_tags_deep(
        s, sprintf("the term developing '%s' inside it", nm)))
    }
  }
  out
}


#' Why a Term Is Outside the Fitting Scheme
#'
#' @description
#' Returns the reason why a term cannot be assembled as a fixed design block,
#' or the empty string when it can.
#'
#' @details
#' There are two reasons, and both are read from the term rather than from
#' its class.
#'
#' A covariance label ([modelterms7::term_tag()]) under a structural term of
#' the likelihood shape would make a latent that the likelihood integrates
#' out share a prior with coefficients that are estimated. The two are
#' integrated by different routes and there is no one prior to share, so the
#' term is rejected. Fitting it as though the label were absent would report
#' a different model under the name of the one requested. A label under a
#' filter is fitted: its parameters are numbers estimated beside the
#' coefficients, and the block is read among them.
#'
#' Whether a class is admissible depends on all its members, and
#' [class_space()] determines it where every member is visible.
#'
#' A structural term that implements neither `term_filter()` nor
#' `term_loglik()` is rejected for the reason that its message gives.
#'
#' @param term One built term.
#'
#' @return A single string.
#'
#' @seealso [reject_unfittable()]
#'
#' @keywords internal
unfittable_reason <- function(term) {
  # A label under a structural term of the LIKELIHOOD shape is the one that
  # cannot be reached, and the obstruction is the model's rather than the
  # addressing's: what such a term carries is a latent the likelihood
  # integrates out, and a covariance block is a prior over coefficients, so
  # correlating the two would mean one joint prior over quantities integrated
  # two different ways. Under a FILTER the parameters are estimated beside the
  # coefficients and the block is read among them, and under an ADDITIVE parent
  # there is nothing extra to do at all -- measured, a labelled effect inside
  # seg(x, psi ~ ...) occupies columns of that term's block and so sits in the
  # same vector as one written in an equation.
  if (S7::S7_inherits(term, modelterms7::structural_term) &&
      !identical(structural_kind(term), "filter")) {
    lab <- term_tags_deep(term)
    lab <- lab[names(lab) != "it"]
    if (length(lab)) {
      return(sprintf(paste("%s carries the covariance label '%s', and this",
                           "term's contribution is a likelihood mixed over a",
                           "latent state rather than a predictor. A covariance",
                           "block is a prior over coefficients, and the latent",
                           "is integrated out by a different route, so there is",
                           "no one prior to share. Drop the middle bar, or move",
                           "the labeled effect to an equation of its own."),
                     names(lab)[1L], lab[[1L]]))
    }
  }
  if (S7::S7_inherits(term, modelterms7::structural_term)) {
    kind <- structural_kind(term)
    if (identical(kind, "filter") || identical(kind, "loglik")) return("")
    return(paste("it is a structural term with a method for neither",
                 "term_filter() nor term_loglik(): a term that rewrites the",
                 "likelihood implements term_filter(), which returns a",
                 "predictor, or term_loglik(), which returns a likelihood."))
  }
  ""
}


#' Does a Term Recompute Its Own Block?
#'
#' @description
#' Returns `TRUE` when the term registers a
#' [modelterms7::term_refresh()] method of its own rather than
#' inheriting the identity registered on `model_term`.
#'
#' @details
#' The owning class of a method is `attr(m, "signature")[[1]]`, and it is
#' compared by name and package rather than by `identical()`: an S7 class
#' re-created from the same definition is not identical to the original, as
#' happens whenever a package's code is re-evaluated rather than loaded.
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


#' Read a Two-Column Response as Successes and Trials
#'
#' @description
#' Splits a `cbind(successes, failures)` response into the successes and the
#' number of trials, returning the response and the distribution onto which
#' the trials have been written.
#'
#' @details
#' R writes aggregated binomial data as `cbind(y, n - y) ~ x`.
#' The two columns are read as the successes and the failures: the response
#' becomes the successes and the family's `size` becomes the row sums.
#'
#' The number of trials is taken from the response and not from the family.
#' A size given to the family is a vector fixed at fitting time, so it cannot
#' follow the rows anywhere else, whereas the row sums of a response
#' expression are recomputed wherever that expression is evaluated, so a fold
#' of [cv()] or a prediction at new data uses the right trials. A size of
#' another length than the rows is recycled by the density without any signal
#' (a size of length 200 against 20 observations returns
#' 200 log-densities). See [check_trials()], which rejects that case where the
#' row sums cannot be had.
#'
#' A size given to the family as well must agree with the row sums, and the
#' function signals an error when the two statements disagree; the default
#' size of one is not a statement and is replaced.
#'
#' A matrix response is left alone for a multivariate family, where it is the
#' ordinary response and each row is one observation.
#'
#' @param response The evaluated left-hand side.
#' @param distrib The family.
#'
#' @return A list with `response`, the successes, and `distrib`, carrying the
#'   row sums as its `size`. Both are unchanged where the response is not a
#'   matrix or the family is multivariate. A matrix response on a univariate
#'   family signals an error where it has a number of columns other than two,
#'   is not numeric, holds a negative count, or where the family has no
#'   `size`.
#'
#' @seealso [coerce_response()], [statmod_spec()], [statmod_respec()]
#'
#' @keywords internal
split_binomial_matrix <- function(response, distrib) {
  invariato <- list(response = response, distrib = distrib)
  if (!is.matrix(response)) return(invariato)
  # a multivariate family reads a matrix response as one observation per row,
  # which is what it is for
  if (identical(as.character(distrib@dimension), "multivariate")) return(invariato)

  nome <- distrib@distrib_name
  if (ncol(response) != 2L) {
    stop(sprintf(
      "The response is a matrix with %d columns and '%s' takes one number per
  observation. Two columns are read as successes and failures, as glm()
  reads them; any other number of columns is not read here.",
      ncol(response), nome), call. = FALSE)
  }
  if (!("size" %in% S7::prop_names(distrib))) {
    stop(sprintf(
      "A two-column response counts successes and failures, and '%s' has no
  number of trials onto which to write them. Use binomial_distrib() or a
  beta-binomial, or pass the response as a single column.",
      nome), call. = FALSE)
  }
  if (!is.numeric(response)) {
    stop("A two-column response must be numeric counts.", call. = FALSE)
  }

  successi <- response[, 1L]
  fallimenti <- response[, 2L]
  if (any(successi < 0 | fallimenti < 0, na.rm = TRUE)) {
    stop("A two-column response counts successes and failures, and one of
  them is negative.", call. = FALSE)
  }
  prove <- successi + fallimenti

  # A size the caller ALSO gave is not overwritten in silence: it either
  # agrees with the row sums, and then there is nothing to choose, or the two
  # statements disagree and only the modeller knows which was meant. The
  # default of one is not a statement and is replaced.
  vecchia <- distrib@size
  predefinita <- length(vecchia) == 1L && isTRUE(vecchia == 1)
  if (!predefinita &&
      !isTRUE(all.equal(as.numeric(rep_len(vecchia, length(prove))),
                        as.numeric(prove)))) {
    stop("The number of trials is given twice and the two disagree: the
  response's row sums give one value and the family's 'size' another.
  Give it once -- either cbind(successes, failures), or a single-column
  response with the numbers of trials as the family's 'size'.",
         call. = FALSE)
  }

  list(response = successi, distrib = S7::set_props(distrib, size = prove))
}


#' Reject a Number of Trials That Cannot Follow These Rows
#'
#' @description
#' Signals an error where the family carries one number of trials per
#' observation and the rows to which it is being applied are a different
#' number.
#'
#' @details
#' A `size` given to `binomial_distrib()` is a vector fixed when the fit was
#' written, so at other rows it is recycled and the density returns a value
#' without any signal: a size of length 200 against 20 observations returns
#' 200 log-densities, so the family and not the data determines the number of
#' terms.
#'
#' The function signals an error because it cannot return a correct number:
#' the fit does not record to which rows the trials belonged. The caller can
#' write the response as `cbind(successes, failures)` instead, whose row sums
#' are recomputed wherever the expression is evaluated, and the error message
#' gives that remedy.
#'
#' @param distrib The family.
#' @param n The number of rows to which it is about to be applied.
#'
#' @return `distrib`, unchanged.
#'
#' @seealso [split_binomial_matrix()], [statmod_respec()]
#'
#' @keywords internal
check_trials <- function(distrib, n) {
  if (!("size" %in% S7::prop_names(distrib))) return(distrib)
  taglia <- distrib@size
  if (length(taglia) <= 1L || length(taglia) == n) return(distrib)
  stop(sprintf(
    "The family carries %d numbers of trials and these are %d rows. A 'size'
  given to the family is fixed when the fit is written and cannot follow
  other rows: recycled against them it returns a number rather than an
  error. Write the response as cbind(successes, failures) instead -- its
  row sums are recomputed wherever the expression is.",
    length(taglia), n), call. = FALSE)
}


#' Reject a Quantity That Needs Trials That the Rows Lack
#'
#' @description
#' Signals an error where the family's numbers of trials are missing, as
#' [statmod_respec()] leaves them when a prediction is made at rows with
#' no response.
#'
#' @details
#' A parameter of a binomial family, its probability, is a function of the
#' linear predictor alone and is predicted at any rows. A moment of the
#' response or a prediction interval for a new observation depends on the
#' number of trials, and at rows that carry no `cbind(successes, failures)`
#' there is none to read.
#'
#' @param distrib The family, as [statmod_respec()] returns it.
#' @param what The quantity requested, in words, for the message.
#'
#' @return `distrib`, unchanged.
#'
#' @seealso [check_trials()], [statmod_respec()], [predict.StatmodFit()]
#'
#' @keywords internal
check_trials_known <- function(distrib, what) {
  if (!("size" %in% S7::prop_names(distrib))) return(distrib)
  if (!anyNA(distrib@size)) return(distrib)
  stop(sprintf(
    "%s needs the number of trials of each row, and 'newdata' has no
  response from which to read it. Add to 'newdata' the variables with which
  the response is written, cbind(successes, failures) in the formula, with
  the real numbers of trials in them; the parameters themselves need no
  response.", what),
    call. = FALSE)
}


#' Carry a Categorical Response onto the Numeric Scale
#'
#' @description
#' Returns the response as a number, converting a two-level `factor`,
#' `character` or `logical` to 0 and 1 and leaving anything else alone.
#'
#' @details
#' A binary outcome is ordinarily recorded as a factor, and `glm()` accepts
#' one: its binomial family reads the first level as the failure. The
#' conversion makes the response a number before it reaches the density,
#' which would otherwise receive a factor.
#'
#' The first level is the failure, which is `glm()`'s rule, so a model moved
#' from `glm()` keeps its signs. The conversion applies to every family and
#' not only to the binary ones, which differs from base R, where
#' `glm(f ~ x, gaussian)` on a factor `f` signals *"NA/NaN/Inf in 'y'"* after
#' three warnings from `Ops.factor`. A two-level factor has one numeric
#' reading and this is it; on a gaussian family the result is the linear
#' probability model for the response `as.numeric(f == "b")`.
#'
#' The levels used are returned as the attribute `"response_levels"`, which
#' lets [statmod_respec()] reproduce the same mapping on other rows. A factor
#' carries its levels through a subset, but a character vector does not, so a
#' cross-validation fold holding only one of the two values would otherwise be
#' coded the other way round, and the deviance would be wrong with no
#' signal.
#'
#' @param response The evaluated left-hand side.
#' @param levels The mapping to reuse, from a previous call's attribute, or
#'   `NULL` to read it off this response.
#'
#' @return `response` unchanged, or a numeric vector of 0 and 1 carrying the
#'   attribute `"response_levels"`.
#'
#' @seealso [statmod_spec()], [statmod_respec()]
#'
#' @keywords internal
coerce_response <- function(response, levels = NULL) {
  if (is.logical(response)) {
    out <- as.integer(response)
    attr(out, "response_levels") <- c("FALSE", "TRUE")
    return(out)
  }
  if (!is.factor(response) && !is.character(response)) return(response)

  lev <- if (!is.null(levels)) levels
         else if (is.factor(response)) base::levels(response)
         else sort(unique(response[!is.na(response)]))
  chr <- as.character(response)

  if (length(lev) != 2L) {
    stop(sprintf(
      "The response is a %s with %d level%s (%s), and a distribution needs a
  number. Two levels are read as 0 and 1, as glm() reads them; any other
  number of levels has no numeric reading defined here. Convert it
  yourself -- as.integer(y == \"<level>\") -- and say which level is the
  one being modeled.",
      if (is.factor(response)) "factor" else "character vector",
      length(lev), if (length(lev) == 1L) "" else "s",
      paste(utils::head(lev, 5), collapse = ", ")), call. = FALSE)
  }
  fuori <- setdiff(unique(chr[!is.na(chr)]), lev)
  if (length(fuori)) {
    stop(sprintf(
      "The response holds %s, which is not among the levels being coded
  (%s). This happens where a fold or a prediction carries a value that was
  not in the data used to fit.",
      paste(sprintf("'%s'", utils::head(fuori, 3)), collapse = ", "),
      paste(lev, collapse = ", ")), call. = FALSE)
  }

  out <- as.numeric(chr == lev[2L])
  out[is.na(chr)] <- NA_real_
  attr(out, "response_levels") <- lev
  out
}


#' Validate Prior Weights
#'
#' @description
#' Returns a vector of prior weights of the right length, defaulting to one
#' per observation.
#'
#' @details
#' They are not normalized. Making them sum to one would turn
#' \eqn{\sum_i w_i \ell_i} into a mean: every standard error would be
#' multiplied by \eqn{\sqrt{n}} and the information criteria would stop being
#' comparable with those of an unweighted fit. The weights act as frequency
#' weights.
#'
#' The weights must be numeric, of length `n`, non-negative, finite and
#' complete, and not all zero; anything else is an error.
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
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("'weights' must be non-negative, finite and complete.",
         call. = FALSE)
  }
  if (all(weights == 0)) {
    stop("'weights' are all zero: there is nothing to fit.", call. = FALSE)
  }
  as.numeric(weights)
}


#' Validate Offsets
#'
#' @description
#' Returns a named list with one offset per parameter, `NULL` for a
#' parameter whose offset was not given.
#'
#' @details
#' An entry must be numeric, of length 1 or `n`, and finite everywhere, as
#' an offset written in the formula must be; the names must be parameters of
#' the distribution.
#'
#' @param offsets A named list, or `NULL`.
#' @param params The parameter names.
#' @param n The number of observations.
#'
#' @return A named list of length `length(params)`, and each entry is a
#'   numeric vector of length `n` or `NULL`.
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
    if (anyNA(o) || !all(is.finite(o))) {
      stop(sprintf(paste0("The offset for '%s' is not finite everywhere.\n",
                          "  An offset enters the linear predictor as it",
                          " stands, so a missing or\n  infinite value has no",
                          " fitted answer."), p), call. = FALSE)
    }
    out[[p]] <- rep_len(as.numeric(o), n)
  }
  out
}


#' The Design of a Specification
#'
#' @description
#' Returns, per parameter, the terms' blocks side by side and the names of the
#' coefficients that they carry.
#'
#' @param spec A [StatmodSpec()].
#' @param unseen A named list, one entry per distribution parameter, of the
#'   keys of the [modelterms7::random()] terms whose rows of a level that was
#'   absent from the fitting data are given zeros instead of an error; `NULL`
#'   (the default) for none. [predict.StatmodFit()] uses it to read a new
#'   group under `random = "zero"` or `"marginal"`.
#' @param skip The distribution parameters whose equations are not rebuilt
#'   from `spec@newdata`: each of their blocks is a matrix of zeros of the
#'   width that the fit gave it. [predict.StatmodFit()] passes the parameters
#'   that a prediction does not read and whose covariates `newdata` lacks;
#'   `character(0)` (the default) for none.
#'
#' @return A named list with one entry per parameter, each a list with `X`,
#'   `asked_sparse` (whether a block is sparse because the caller requested
#'   it, absent for an equation that holds only structural terms),
#'   `coef_names`, `npar` and `blocks` (the column range that each term
#'   occupies).
#'
#'   The list carries attributes where the model needs them. `structural` and
#'   `structure` are present when a structural term is; `structural` lists
#'   the structural units and `structure` is an environment with the term's
#'   current parameters and held levels. `refresh` and `state` are present when
#'   a term recomputes its own block; `refresh` lists those terms and `state`
#'   is an environment with the terms from which they chain. When every block
#'   is independent of the coefficients and the design has no structural term,
#'   the attribute `eta_memo` is an environment holding the last predictors
#'   that [statmod_eta()] computed on this design, so a second call at the
#'   same coefficients returns them without recomputing.
#'
#' @seealso [statmod_spec()]
#'
#' @examples
#' dd <- data.frame(y = rnorm(20), x = runif(20))
#' d <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
#' vapply(statmod_design(d), function(z) z$npar, integer(1))
#'
#' @export
statmod_design <- function(spec, unseen = NULL, skip = character(0)) {
  out <- statmod_design_blocks(spec, unseen, skip)
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
    fresh <- vapply(su, function(u) is.null(spec@structural[[u$term]]),
                    logical(1))
    structural_start_fixups(spec, sst, su, fresh)
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
  # the last point statmod_eta() was asked for, where nothing in the design
  # moves with the coefficients; an environment for the reason above
  if (!length(su) && !length(rf)) {
    out <- densify_small(out)
    mm <- new.env(parent = emptyenv())
    mm$key <- NULL
    mm$value <- NULL
    attr(out, "eta_memo") <- mm
  }
  out
}


#' Store a Small Design as Base Matrices
#'
#' @description
#' Converts the sparse blocks of a design to base matrices when the
#' equations together carry fewer than `min_dim` coefficients and the caller
#' has not requested sparse storage for any block.
#'
#' @details
#' A random effect builds its block sparse whatever its size. Below about a
#' hundred coefficients the sparse route costs more than it saves, because its
#' fixed cost is the coercions and the S4 dispatch around each product, which
#' do not shrink with the matrix; [worth_sparse()] records the same crossover
#' for the factorization. The conversion is applied only where every block is
#' independent of the coefficients and the design has no structural term,
#' since other designs are rebuilt during the fit and a rebuilt block would
#' come back sparse. A block that the caller requested sparse
#' (`sparse = TRUE` on a term, `linpar_options(sparse = TRUE)`, a sparse
#' matrix as input) is never converted, and the rest of the design is not
#' converted either: the design keeps the storage that it was given.
#'
#' @param design The design, a list with one entry per distribution
#'   parameter, and each entry holds its block as `X`.
#' @param min_dim The number of coefficients below which the design is
#'   stored dense.
#'
#' @return The design, with every `X` a base matrix where the total is below
#'   `min_dim`, and unchanged otherwise.
#'
#' @keywords internal
densify_small <- function(design, min_dim = 100L) {
  total <- sum(vapply(design, function(d) as.integer(d$npar), integer(1)))
  if (total >= min_dim) return(design)
  if (any(vapply(design, function(d) isTRUE(d$asked_sparse), logical(1)))) {
    return(design)
  }
  for (p in names(design)) {
    if (isS4(design[[p]]$X)) design[[p]]$X <- as.matrix(design[[p]]$X)
  }
  design
}


#' Settle a Fresh Start of a Structural Term Against Its Equation
#'
#' @description
#' Adjusts the start of a structural term's own parameters once the design
#' is known: a level held by an intercept in the same equation starts at zero,
#' and a level that is not held and that the term's own start left at zero
#' starts at the equation's data-based intercept. A filter also
#' receives a second start, kept beside the first, in which each score loading
#' equals \eqn{0.1\,\mathcal{I}^{d-1}}, where \eqn{d} is the `scaling` of the
#' filter (so that the loading is \eqn{0.1/\mathcal{I}} for \eqn{d = 0}) and
#' \eqn{\mathcal{I}} is the expected information of its equation's predictor
#' at the intercept-only fit ([predictor_information()]); [statmod()] fits
#' from it where the fit from the first start does not converge.
#'
#' @details
#' [statmod_design()] applies these adjustments to the start that
#' [modelterms7::term_start()] gives, and [statmod()] applies them again to
#' every further start that [modelterms7::term_starts()] gives, so each start
#' of a multistart fit is treated as the first one is.
#'
#' @param spec A [StatmodSpec()].
#' @param sst The design's structural state, an environment carrying `zeta`
#'   and `held`.
#' @param su The structural units, as `attr(design, "structural")`.
#' @param fresh A logical vector, one per unit: `TRUE` where the unit's
#'   values are a fresh start rather than values at which a fit arrived.
#'
#' @return `NULL`, invisibly; `sst$zeta` is modified in place.
#'
#' @keywords internal
structural_start_fixups <- function(spec, sst, su, fresh) {
  # a HELD coordinate is held at zero, the value the intercept beside it
  # takes over from: a term's start may read the data and put the level
  # elsewhere, and a held level left there would shift every regime by a
  # constant the fit can never remove
  for (i in seq_along(su)) {
    u <- su[[i]]
    if (!fresh[i]) next
    h <- sst$held[[u$term]]
    if (length(h)) sst$zeta[[u$term]][h] <- 0
  }
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
      # a term whose own start already read the data keeps it: a regime
      # term places its first level at a low quantile of the response,
      # which the intercept-only fit would move to the mean
      if (!isTRUE(sst$zeta[[u$term]][[lvl]] == 0)) next
      lk <-modelterms7::term_links(tm)[[lvl]]
      z <- tryCatch(linkfunctions7::linkfun(lk, v),
                    error = function(e) NA_real_)
      if (is.finite(z)) sst$zeta[[u$term]][[lvl]] <- z
    }
  }
  # A SECOND START FOR A SCORE LOADING, ON THE SCALE OF THE SCORE IT
  # MULTIPLIES. The score is used unscaled, so its size is the family's:
  # (y - mu)/sigma^2 for a gaussian mean, of order phi for a beta on the
  # logit. The term's own start, a loading of 0.1, is a weak response only
  # where the score is of order one; on a beta share with phi near 130 it
  # made the filter explode at the first step and the fit stopped,
  # unconverged, at -3326.66 where the static model reaches 413.58. Read on
  # the scale of a score scaled by the inverse information, 0.1 is the same
  # weak response for every family: 0.1 / I, I the expected information of
  # the equation's predictor at the intercept-only fit, from which the beta
  # converges at 450.48. It is kept in sst$info_start, and statmod() fits from
  # it where the fit from the term's own start does not converge (Giovanni,
  # 2026-10-03). Used as THE start it converged on 34 of 40 simulated
  # gaussian panels with a random level where 0.1 converged on 37, and on the
  # battery it took gas-panel-omega to a lower REML optimum (-609.528 against
  # -609.490), so it is the second start and not the first.
  # With gas(scaling = d) the score is multiplied by I^-d, so the same weak
  # response is 0.1 * I^(d - 1): 0.1 / I at d = 0, and the term's own 0.1 at
  # d = 1, where the score already is on that scale.
  filt <- which(vapply(seq_along(su), function(i)
    fresh[i] && identical(su[[i]]$kind, "filter"), logical(1)))
  if (length(filt)) {
    info <- tryCatch(predictor_information(spec), error = function(e) NULL)
    for (i in filt) {
      u <- su[[i]]
      I <- info[[u$param]]
      if (is.null(I) || !is.finite(I) || I <= 0) next
      tm <- spec@terms[[u$param]][[u$term]]
      lks <- modelterms7::term_links(tm)
      z <- sst$zeta[[u$term]]
      # 0.1 / I written as it was at d = 0: 0.1 * I^-1 can differ from it in
      # the last bit, and a fit started from either can end elsewhere
      d_sc <- filter_scaling(tm)
      a_start <- if (d_sc == 0) 0.1 / I else 0.1 * I^(d_sc - 1)
      for (a in grep("^kappa[0-9]+$", modelterms7::term_params(tm),
                     value = TRUE)) {
        lk <- lks[[a]]
        if (is.null(lk)) lk <- linkfunctions7::log_link()
        new <- tryCatch(linkfunctions7::linkfun(lk, a_start),
                        error = function(e) NA_real_)
        old <- tryCatch(linkfunctions7::linkfun(lk, 0.1),
                        error = function(e) NA_real_)
        if (!is.finite(new)) next
        if (a %in% names(z)) {
          z[[a]] <- new
        } else if (is.finite(old) && old != 0) {
          # a developed loading starts at the constant projected onto its
          # design, which is linear in the constant
          dev <- startsWith(names(z), paste0(a, "."))
          z[dev] <- z[dev] * (new / old)
        }
      }
      if (is.null(sst$info_start)) sst$info_start <- list()
      sst$info_start[[u$term]] <- z
    }
  }
  invisible(NULL)
}


#' The Expected Information of Each Predictor at the Intercept-Only Fit
#'
#' @description
#' Computes, for each distribution parameter, the mean over the observations
#' of the expected information of its linear predictor,
#' \eqn{-\mathrm{E}[\partial^2 \ell / \partial \eta^2]}, at the parameters of
#' the intercept-only fit.
#'
#' @details
#' It sets the scale of a score-driven term's second start
#' ([structural_start_fixups()]): the score of a predictor has variance equal
#' to this information, so a loading of \eqn{0.1/\mathcal{I}} on the unscaled
#' score moves the predictor by a tenth of the score's natural unit whatever
#' the family. For a filter with `scaling` \eqn{d} the loading is
#' \eqn{0.1\,\mathcal{I}^{d-1}}, which gives the same step.
#'
#' @param spec A [StatmodSpec()].
#'
#' @return A named numeric vector, one entry per distribution parameter, or
#'   `NULL` where the intercept-only fit is not available.
#'
#' @keywords internal
predictor_information <- function(spec) {
  d <- spec@distrib
  y <- spec@response
  eta0 <- statmod_intercepts(spec)
  if (is.null(eta0)) return(NULL)
  n <- spec@n_obs
  th <- lapply(d@params, function(q) {
    v <- eta0[[q]]
    if (is.null(v) || !is.finite(v)) v <- 0
    rep_len(linkfunctions7::linkinv(d@link_params[[q]], v), n)
  })
  names(th) <- d@params
  eh <- distributions7::distrib_expected_hessian(d, y, th, scale = "link")
  vapply(d@params, function(q) {
    v <- eh[[paste(q, q, sep = "_")]]
    if (is.null(v)) NA_real_ else -mean(as.numeric(v))
  }, numeric(1))
}


#' The Blocks of a Design
#'
#' @description
#' Returns the per-parameter blocks, before any refreshable term is
#' recomputed at coefficients.
#'
#' @param spec A [StatmodSpec()].
#' @param unseen,skip As in [statmod_design()].
#'
#' @return A named list, one entry per distribution parameter.
#'
#' @keywords internal
statmod_design_blocks <- function(spec, unseen = NULL, skip = character(0)) {
  stats::setNames(lapply(names(spec@terms), function(p) {
    tms <- spec@terms[[p]]
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
      else if (p %in% skip) lapply(tms, function(tm)
        matrix(0, spec@n_obs, length(modelterms7::term_coef_names(tm))))
      else stats::setNames(lapply(names(tms), function(k) {
        if (k %in% unseen[[p]]) {
          modelterms7::term_predict(tms[[k]], newdata = spec@newdata,
                                    unseen = "zero")
        } else modelterms7::term_predict(tms[[k]], newdata = spec@newdata)
      }), names(tms))
    nms <- lapply(tms, modelterms7::term_coef_names)
    widths <- vapply(mats, ncol, integer(1))
    ends <- cumsum(widths)
    starts <- ends - widths + 1L
    list(
      X = bind_blocks(mats, spec@n_obs),
      # a block sparse because the caller asked for it, as distinct from the
      # block random() builds sparse whatever its size; densify_small() never
      # touches the first
      asked_sparse = any(vapply(names(tms), function(k)
        isS4(mats[[k]]) && !S7::S7_inherits(tms[[k]], modelterms7::RandomTerm),
        logical(1))),
      coef_names = unlist(nms, use.names = FALSE),
      npar = as.integer(sum(widths)),
      blocks = stats::setNames(
        lapply(seq_along(tms), function(j)
          if (widths[j] == 0L) integer(0) else seq.int(starts[j], ends[j])),
        names(tms))
    )
  }), names(spec@terms))
}


#' Check the Signatures of a Family's Own Methods
#'
#' @description
#' Signals an error naming the method and the signature that it needs where a
#' family registers a method of its own, for the density or for a
#' derivative, whose formals have no `...`.
#'
#' @details
#' [statmod()] passes further arguments to these methods, including the
#' number of threads, and a method without `...` fails from inside the fit
#' with "unused argument (threads = 1)". The derivative methods also take
#' `scale`, which the generic resolves and passes on. The methods of the
#' families shipped with \pkg{distributions7} all carry `...`, so the error
#' arises for a family written by a user.
#'
#' @param distrib A univariate family.
#'
#' @return `distrib`, invisibly.
#'
#' @keywords internal
assert_method_signatures <- function(distrib) {
  need <- list(
    distrib_pdf = "function(distrib, y, theta, log = FALSE, ...)",
    distrib_gradient = "function(distrib, y, theta, scale = c(\"parameter\", \"link\"), ...)",
    distrib_hessian = "function(distrib, y, theta, scale = c(\"parameter\", \"link\"), ...)",
    distrib_deriv3 = "function(distrib, y, theta, scale = c(\"parameter\", \"link\"), ...)",
    distrib_deriv4 = "function(distrib, y, theta, scale = c(\"parameter\", \"link\"), ...)",
    distrib_expected_hessian = "function(distrib, y, theta, scale = c(\"parameter\", \"link\"), ...)")
  cls <- S7::S7_class(distrib)
  for (g in names(need)) {
    m <- tryCatch(S7::method(getExportedValue("distributions7", g), cls),
                  error = function(e) NULL)
    if (is.null(m) || "..." %in% names(formals(m))) next
    stop(sprintf(paste0(
      "the %s() method of '%s' has no `...` in its signature, and statmod() ",
      "passes\n  further arguments to it (the number of threads). Write it ",
      "as\n  %s."), g, distrib@distrib_name, need[[g]]), call. = FALSE)
  }
  invisible(distrib)
}
