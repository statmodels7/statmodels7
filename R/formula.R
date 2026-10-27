#' Split a Multi-Parameter Formula Into One Equation Per Parameter
#'
#' @description
#' Takes the single formula that [statmod()] accepts, in which the equations
#' of the distribution's parameters are separated by `|`, and returns the
#' response together with one one-sided formula per parameter of the family.
#' Every parameter receives an equation, and a parameter for which the caller
#' did not write an equation receives `~ 1`, an intercept.
#'
#' @details
#' # The syntax
#'
#' \preformatted{    y ~ <terms>  |  p2 ~ <terms>  |  p3 ~ <terms>}
#'
#' The first equation carries the response and models the family's first
#' parameter, and each `|` introduces the next. The parameters may be named
#' in any order, and any of them may be omitted.
#'
#' # Splitting the formula
#'
#' The split follows R's own precedence. `~` binds looser than `|` and
#' associates to the left, so `y ~ a | p2 ~ b | p3 ~ c` parses as
#'
#' \preformatted{    ((y ~ (a | p2)) ~ (b | p3)) ~ c}
#'
#' The right-hand side of the whole formula is `c` alone. Splitting it on `|`
#' returns one piece and drops two equations without an error.
#'
#' The function therefore walks down the left spine of the nested `~` calls,
#' collects the right-hand side of each level and reverses the collection.
#' The innermost left operand is the response. Each collected piece is
#' either `<terms> | <name of the next parameter>` or, for the last piece,
#' `<terms>` alone.
#'
#' # Bars inside a call
#'
#' The walk descends only through `~` and through a `|` at the top of a
#' collected piece. A `|` anywhere inside a call is left as it stands, so
#' `random(~ 1 | id)` and `gas(omega ~ random(~ 1 | id), by = id)` reach
#' the term constructor whole.
#'
#' # Invalid formulas
#'
#' Five conditions signal an error:
#'
#' - `formula` is not a formula.
#' - `formula` has no left-hand side, so there is no response.
#' - A `|` is not followed by an equation, as in `y ~ x | sigma`. The
#'   message names the parameter whose equation is missing.
#' - A name is not a parameter of this distribution. The message names it
#'   and lists the parameters of the distribution.
#' - A parameter is given two equations. The message names the parameter.
#'
#' @param formula A two-sided formula in which the equations of the
#'   parameters are separated by `|` as above. Its environment is carried
#'   onto every equation returned, so a term's symbols resolve where the
#'   caller wrote them; if the formula has a `NULL` environment, [baseenv()]
#'   is used.
#' @param params The distribution's parameter names in the family's own
#'   order, as `distrib@params` gives them. It is a character vector. The
#'   first element is the parameter to which the response's own equation
#'   belongs.
#'
#' @return A list of three:
#'   \describe{
#'     \item{`response`}{the unevaluated left-hand side, a language object,
#'       usually a symbol but any expression that the caller wrote.}
#'     \item{`equations`}{a named list of one-sided formulas, one per element
#'       of `params` and in that order, regardless of the order in which the
#'       caller wrote them. A parameter with no equation of its own holds
#'       `~ 1`.}
#'     \item{`given`}{the parameter names that the formula supplied, in the
#'       order in which it supplied them. It is a subset of `params`, possibly
#'       of length one.}
#'   }
#'
#' @seealso [statmod()], which calls this on the formula that it is given,
#'   [statmod_spec()] for what is built from the result.
#'
#' @examples
#' e <- statmod_equations(y ~ x1 + x2 | sigma ~ z, c("mu", "sigma"))
#' e$response
#' e$equations
#'
#' # Every parameter receives an equation; a parameter without one receives an
#' # intercept.
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
#' h <- statmod_equations(y ~ random(~ 1 | id) | sigma ~ z, c("mu", "sigma"))
#' deparse(h$equations$mu[[2]])
#'
#' # A name that is not a parameter is an error, and the message lists the
#' # parameters of the distribution.
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
#' symbols in a term resolve where the caller wrote them.
#' [statmod_equations()] uses it on each piece that it collects, and on
#' `quote(1)` for a parameter with no equation.
#'
#' @details
#' The formula is built with `eval(bquote(~ .(expr)), envir = env)`, and the
#' environment is then assigned again. `eval()` gives the formula the
#' environment in which it was evaluated, and the assignment states the
#' intended environment explicitly, so the result carries `env` whatever the
#' construction does.
#'
#' @param expr A language object, the right-hand side to wrap: a symbol, a
#'   call, or the literal `1`. It is substituted unevaluated, so a term call
#'   is not run here.
#' @param env The environment to attach. It is not checked, and `NULL` gives a
#'   formula with a `NULL` environment.
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
#' Returns a fresh environment holding every function that \pkg{modelterms7}
#' exports, whose parent is `env`. A term call evaluated there reaches
#' \pkg{modelterms7}'s function whatever the caller has attached, and every
#' other name that the caller wrote resolves as usual, one step further out.
#'
#' @details
#' \pkg{mgcv} exports `s()` and `te()`, and \pkg{segmented} exports `seg()`.
#' If either package is attached, a `statmod()` formula that uses one of these
#' names would call that package's function, whose value is not a model term.
#' The failure then occurs inside `model.matrix()` as `invalid type (list) for
#' variable 's(x)'`, a message that does not mention the masking package.
#'
#' Putting \pkg{modelterms7} in front removes the ambiguity at the point
#' where the formula is read. A caller who wants the other package's function
#' writes `mgcv::s(x)`, which reaches past the shim; the formula interpreter
#' of \pkg{modelterms7} then rejects the value, naming the class returned
#' and the call as written.
#'
#' The function binds each export of \pkg{modelterms7} for which
#' `is.function()` is `TRUE`. The S7 class generators of that package are
#' functions and are bound with the others. Each name is an active binding
#' made by [export_binding()], which reads the object from the namespace when
#' the name is evaluated. A term keeps this environment as the environment of
#' its formula, and copies of the objects would be serialized with every fit:
#' the S7 generics and classes of \pkg{modelterms7} weigh about 16 MB, against
#' about 0.1 MB for the bindings.
#'
#' @param env The environment that the formula carried, which becomes the
#'   parent of the result. Everything visible from `env` stays visible.
#'
#' @return A new environment holding one binding per exported function of
#'   \pkg{modelterms7}, with `env` as its parent.
#'
#' @seealso [statmod_equations()] for the formula that is interpreted with
#'   this function, [modelterms7::interpret_formula()] for the reading that
#'   it feeds.
#'
#' @keywords internal
terms_first <- function(env) {
  shim <- new.env(parent = env)
  ns <- asNamespace("modelterms7")
  for (nm in getNamespaceExports("modelterms7")) {
    if (is.function(get(nm, envir = ns))) {
      makeActiveBinding(nm, export_binding(nm), shim)
    }
  }
  shim
}

#' The Active Binding of One Export of modelterms7
#'
#' @description
#' Returns a function of no arguments that gets the object called `name` from
#' the namespace of \pkg{modelterms7}. It is the function that
#' [makeActiveBinding()] takes in [terms_first()].
#'
#' @details
#' The function is built from an expression, with the environment of the base
#' package, so that it holds nothing but the name. A closure written in the
#' body of [terms_first()] would carry the environment that the caller passed,
#' and a closure written here would carry its source reference: a package
#' installed with the source kept (`R_KEEP_PKG_SOURCE=yes`) serializes the
#' source of the whole package with each such reference.
#'
#' @param name A single string, the name of an export.
#'
#' @return A function of no arguments.
#'
#' @seealso [terms_first()]
#'
#' @keywords internal
export_binding <- function(name) {
  as.function(list(bquote(get(.(name), envir = asNamespace("modelterms7")))),
              envir = baseenv())
}
