#' Split a Multi-Parameter Formula Into One Equation Per Parameter
#'
#' @description
#' Takes the single formula \code{\link{statmod}} accepts, in which the
#' equations of the distribution's parameters are separated by \code{|}, and
#' returns the response together with one one-sided formula per parameter.
#'
#' @details
#' The syntax is
#' \preformatted{    y ~ <terms>  |  p2 ~ <terms>  |  p3 ~ <terms>}
#' with the first equation carrying the response and modelling the family's
#' first parameter, and each \code{|} introducing the next.
#'
#' The recovery is not the obvious one, and the reason is R's own precedence.
#' \code{~} binds looser than \code{|} and associates to the left, so
#' \code{y ~ a | p2 ~ b | p3 ~ c} is the tree
#' \code{((y ~ (a | p2)) ~ (b | p3)) ~ c}: the right-hand side of the whole
#' formula is \code{c} alone, and splitting it on \code{|} returns one piece
#' and silently drops two equations. What the associativity calls for instead
#' is a walk down the left spine of the nested \code{~} calls, collecting each
#' level's right-hand side; the innermost left operand is the response, and
#' each collected piece is either \code{<terms> | <name of the next
#' parameter>} or, for the last, \code{<terms>} alone.
#'
#' A \code{|} inside a call is untouched, the walk descending only through
#' \code{~} and the top-level \code{|}, so \code{random(1 | id)} and
#' \code{gas(by = ~ random(1 | id))} survive intact.
#'
#' @param formula The model formula.
#' @param params The distribution's parameter names, in the family's order.
#'
#' @return A list with \code{response} (the unevaluated left-hand side),
#'   \code{equations} (a named list of one-sided formulas, one per element of
#'   \code{params}, in that order) and \code{given} (the names the formula
#'   actually supplied).
#'
#' @seealso \code{\link{statmod}}
#'
#' @examples
#' statmod_equations(y ~ x1 + x2 | sigma ~ z, c("mu", "sigma"))
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
#' Wraps a term expression as \code{~ expr} carrying the environment the
#' original formula had, so that a term's symbols resolve where the user
#' wrote them.
#'
#' @param expr A language object.
#' @param env The environment to attach.
#'
#' @return A one-sided formula.
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
#' Returns an environment in which a term call resolves to modelterms7's
#' function whatever the user has attached, with everything else visible
#' behind it.
#'
#' @details
#' \pkg{mgcv} exports \code{s()} and \code{te()} and \pkg{segmented} exports
#' \code{seg()}. With either attached, a user writing a statmod formula gets
#' the other package's function, and the failure surfaces inside
#' \code{model.matrix} naming neither the call nor the mask. Interpreting the
#' formula in an environment whose parent chain reaches modelterms7 first
#' removes the ambiguity rather than reporting it; a user who wants the other
#' package's term writes \code{mgcv::s(x)}, which the interpreter then
#' rejects with the message it already has.
#'
#' @param env The environment the formula carried.
#'
#' @return An environment.
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
