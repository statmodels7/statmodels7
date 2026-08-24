#' Split a Multi-Parameter Formula Into One Equation Per Parameter
#'
#' @description
#' Takes the single formula [statmod()] accepts, in which the equations of the
#' distribution's parameters are separated by `|`, and returns the response
#' together with one one-sided formula per parameter of the family. Every
#' parameter gets an equation: one a caller did not write is filled in as
#' `~ 1`, an intercept.
#'
#' @details
#' # The syntax
#'
#' \preformatted{    y ~ <terms>  |  p2 ~ <terms>  |  p3 ~ <terms>}
#'
#' The first equation carries the response and models the family's first
#' parameter, and each `|` introduces the next. The parameters may be named
#' in any order and any of them may be omitted.
#'
#' # Why the recovery is a walk
#'
#' R's own precedence decides this. `~` binds looser than `|` and associates
#' to the left, so `y ~ a | p2 ~ b | p3 ~ c` parses as
#'
#' \preformatted{    ((y ~ (a | p2)) ~ (b | p3)) ~ c}
#'
#' The right-hand side of the whole formula is `c` alone. Splitting it on `|`
#' returns one piece and drops two equations without an error.
#'
#' What the associativity calls for is a walk down the left spine of the
#' nested `~` calls, collecting each level's right-hand side and then
#' reversing. The innermost left operand is the response. Each collected
#' piece is either `<terms> | <name of the next parameter>`, or, for the
#' last, `<terms>` alone.
#'
#' # Bars inside a call survive
#'
#' The walk descends only through `~` and through a `|` at the top of a
#' collected piece. A `|` anywhere inside a call is left as it stands, so
#' `random(1 | id)` and `gas(by = ~ random(1 | id))` reach the term
#' constructor whole.
#'
#' # What is refused
#'
#' Five conditions, each with an error naming the offending name:
#'
#' - `formula` is not a formula.
#' - `formula` has no left-hand side, so there is no response.
#' - A `|` is not followed by an equation, as in `y ~ x | sigma`.
#' - A name is not a parameter of this distribution. The message lists the
#'   ones that are.
#' - A parameter is given two equations.
#'
#' @param formula A two-sided formula, with the parameters' equations
#'   separated by `|` as above. Its environment is carried onto every
#'   equation returned, so a term's symbols resolve where the caller wrote
#'   them; a formula with a `NULL` environment gets [baseenv()].
#' @param params The distribution's parameter names in the family's own
#'   order, as `distrib@params` gives them. A character vector. The first
#'   element is the parameter the response's own equation belongs to.
#'
#' @return A list of three:
#'   \describe{
#'     \item{`response`}{the unevaluated left-hand side, a language object,
#'       usually a symbol but any expression the caller wrote.}
#'     \item{`equations`}{a named list of one-sided formulas, one per element
#'       of `params` and in that order, whatever order the caller wrote them
#'       in. A parameter with no equation of its own holds `~ 1`.}
#'     \item{`given`}{the parameter names the formula supplied, in the order
#'       it supplied them. A subset of `params`, possibly of length one.}
#'   }
#'
#' @seealso [statmod()], which calls this on the formula it is given,
#'   [statmod_spec()] for what is built from the result.
#'
#' @examples
#' e <- statmod_equations(y ~ x1 + x2 | sigma ~ z, c("mu", "sigma"))
#' e$response
#' e$equations
#'
#' # Every parameter gets an equation; the ones not written get an intercept.
#' f <- statmod_equations(y ~ x, c("mu", "sigma", "nu"))
#' vapply(f$equations, function(q) deparse(q[[2]]), character(1))
#' f$given
#'
#' # Three equations survive, where splitting the right-hand side on "|"
#' # would keep only the last.
#' g <- y ~ a | sigma ~ b | nu ~ cc
#' deparse(g[[3]])                       # the whole right-hand side: just "cc"
#' statmod_equations(g, c("mu", "sigma", "nu"))$given
#'
#' # A bar inside a call is not a separator.
#' h <- statmod_equations(y ~ random(1 | id) | sigma ~ z, c("mu", "sigma"))
#' deparse(h$equations$mu[[2]])
#'
#' # A name that is not a parameter is refused, and the message says which
#' # names are.
#' try(statmod_equations(y ~ x | zz ~ z, c("mu", "sigma")))
#'
#' @export
statmod_equations <- function(formula, params) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (length(formula) != 3L) {
    stop("'formula' must have a left-hand side: the response.", call. = FALSE)
  }
  env <- environment(formula)
  if (is.null(env)) env <- baseenv()

  # walk the left spine, collecting each level's right-hand side
  rhs <- list()
  e <- formula
  while (is.call(e) && identical(e[[1L]], quote(`~`)) && length(e) == 3L) {
    rhs[[length(rhs) + 1L]] <- e[[3L]]
    e <- e[[2L]]
  }
  response <- e
  rhs <- rev(rhs)

  # each piece but the last names the parameter of the piece after it
  names_ <- params[1L]
  terms_ <- vector("list", length(rhs))
  for (i in seq_along(rhs)) {
    p <- rhs[[i]]
    if (is.call(p) && identical(p[[1L]], quote(`|`)) && length(p) == 3L) {
      terms_[[i]] <- p[[2L]]
      names_ <- c(names_, deparse1(p[[3L]]))
    } else {
      terms_[[i]] <- p
    }
  }
  if (length(names_) > length(terms_)) {
    stop(sprintf(paste0("The equation for '%s' is missing: a '|' introduces a\n",
                        "  parameter and must be followed by '%s ~ <terms>'."),
                 names_[length(names_)], names_[length(names_)]),
         call. = FALSE)
  }
  given <- names_[seq_along(terms_)]

  bad <- setdiff(given, params)
  if (length(bad)) {
    stop(sprintf(paste0("'%s' is not a parameter of this distribution.\n",
                        "  Its parameters are: %s."),
                 bad[1L], paste(params, collapse = ", ")), call. = FALSE)
  }
  dup <- given[duplicated(given)]
  if (length(dup)) {
    stop(sprintf("Parameter '%s' has more than one equation.", dup[1L]),
         call. = FALSE)
  }

  # a parameter with no equation gets an intercept
  eqs <- stats::setNames(vector("list", length(params)), params)
  for (i in seq_along(given)) {
    eqs[[given[i]]] <- one_sided(terms_[[i]], env)
  }
  for (p in params) {
    if (is.null(eqs[[p]])) eqs[[p]] <- one_sided(quote(1), env)
  }

  list(response = response, equations = eqs, given = given)
}


#' Build a One-Sided Formula From an Expression
#'
#' @description
#' Wraps a term expression as `~ expr` and attaches `env` to it, so that the
#' symbols in a term resolve where the caller wrote them. Used by
#' [statmod_equations()] on each piece it collects, and on `quote(1)` for a
#' parameter with no equation.
#'
#' @details
#' The formula is built with `eval(bquote(~ .(expr)), envir = env)`, and the
#' environment is then assigned again. The second step makes the result
#' reliable: `eval()` gives the formula whatever environment it was evaluated
#' in, and the assignment states the intended one, so the two cannot come
#' apart if the construction changes.
#'
#' @param expr A language object, the right-hand side to wrap: a symbol, a
#'   call, or the literal `1`. Substituted unevaluated, so a term call is not
#'   run here.
#' @param env The environment to attach. Passed through with no check, so a
#'   caller supplying `NULL` gets a formula with a `NULL` environment rather
#'   than an error.
#'
#' @return A one-sided formula of length 2, whose `[[2]]` is `expr` and whose
#'   `environment()` is `env`.
#'
#' @seealso [statmod_equations()], its only caller.
#'
#' @keywords internal
one_sided <- function(expr, env) {
  f <- eval(bquote(~ .(expr)), envir = env)
  environment(f) <- env
  f
}


#' Evaluate a Formula's Terms With modelterms7 in Front
#'
#' @description
#' Returns a fresh environment holding every function \pkg{modelterms7}
#' exports, whose parent is `env`. A term call evaluated there reaches
#' \pkg{modelterms7}'s function whatever the caller has attached, and every
#' other name the caller wrote resolves as usual, one step further out.
#'
#' @details
#' The masking this removes is real and its symptom names nothing.
#' \pkg{mgcv} exports `s()` and `te()`, and \pkg{segmented} exports `seg()`.
#' With either package attached, a caller writing a `statmod()` formula gets
#' that package's function, whose value is not a model term; the failure then
#' surfaces inside `model.matrix()` as `invalid type (list) for variable
#' 's(x)'`, which names neither the call nor the mask.
#'
#' Putting \pkg{modelterms7} in front removes the ambiguity at the point
#' where the formula is read. A caller who wants the other package's function
#' writes `mgcv::s(x)`, which reaches past the shim; \pkg{modelterms7}'s
#' formula interpreter then rejects it, naming the class returned and the
#' package that supplied it.
#'
#' Only functions are copied. \pkg{modelterms7} exports objects that are not
#' functions, S7 classes among them, and shadowing those would change what a
#' name means without any of the benefit.
#'
#' @param env The environment the formula carried, which becomes the parent
#'   of the result. Everything visible from `env` stays visible.
#'
#' @return A new environment holding one binding per exported function of
#'   \pkg{modelterms7}, currently 82 of them, with `env` as its parent.
#'
#' @seealso [statmod_equations()] for the formula this is used to interpret,
#'   [modelterms7::interpret_formula()] for the reading it feeds.
#'
#' @keywords internal
terms_first <- function(env) {
  shim <- new.env(parent = env)
  ns <- asNamespace("modelterms7")
  for (nm in getNamespaceExports("modelterms7")) {
    obj <- get(nm, envir = ns)
    if (is.function(obj)) assign(nm, obj, envir = shim)
  }
  shim
}
