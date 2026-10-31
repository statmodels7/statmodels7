#' @rdname statmod_variant_search
#' @format `statmod_variants` is an environment holding, during a search,
#'   whether the search is running (`active`) and which terms replace which
#'   in the formula (`override`, a list by parameter of terms keyed by the
#'   call as written, which [statmod_terms()] reads); it is empty outside a
#'   search.
#' @keywords internal
statmod_variants <- new.env(parent = emptyenv())


#' The Terms That Carry a Choice for a Search
#'
#' @description
#' Returns, for every term of the specification whose
#' [modelterms7::term_pilot()] is not `NULL`, the parameter, the call as
#' written, the built term and its pilot.
#'
#' @param spec A [StatmodSpec].
#'
#' @return A list, empty where no term needs a search.
#'
#' @keywords internal
variant_terms <- function(spec) {
  out <- list()
  for (p in names(spec@terms)) {
    tms <- spec@terms[[p]]
    for (nm in names(tms)) {
      tm <- tms[[nm]]
      if (!S7::S7_inherits(tm, modelterms7::model_term)) next
      pl <- modelterms7::term_pilot(tm)
      if (is.null(pl)) next
      out[[length(out) + 1L]] <- list(param = p, name = nm, key = paste(p, nm, sep = ":"),
                                      term = tm, pilot = pl)
    }
  }
  out
}


#' The Score That Ranks the Versions of a Term
#'
#' @description
#' Returns a number that is smaller for a better fit: the outer criterion
#' with the sign that makes it a minimum, or the penalized objective where no
#' outer criterion ran.
#'
#' @param fit A [StatmodFit].
#'
#' @return A single number.
#'
#' @keywords internal
variant_score <- function(fit) {
  m <- fit@methods$outer
  if (!is.null(m) && length(fit@criterion) && is.finite(fit@criterion)) {
    return(if (outer_minimize(m)) fit@criterion else -fit@criterion)
  }
  fit@objective
}


#' A Fit Searched over the Versions of a Term
#'
#' @description
#' Fits the model once with each term that carries a choice replaced by its
#' pilot ([modelterms7::term_pilot()]), then, one such term at a time, once
#' for each of its versions ([modelterms7::term_variants()]), and returns the
#' fit whose score ([variant_score()]) is the smallest.
#'
#' @details
#' Every fit is the call of [statmod()] as written, evaluated where it was
#' made, with the replacement recorded in [statmod_variants]; the search is
#' not entered again from inside. Where several terms carry a choice, each is
#' settled in turn with the earlier ones at their chosen version and the later
#' ones at their pilot. The versions tried, with their scores, are in
#' `fit@history$variants`, and the fit's `elapsed` covers the whole search.
#'
#' @param cl The call of [statmod()].
#' @param caller The environment in which it was made.
#' @param auto The terms that carry a choice, as [variant_terms()] returns
#'   them.
#'
#' @return A [StatmodFit].
#'
#' @keywords internal
statmod_variant_search <- function(cl, caller, auto) {
  t0 <- proc.time()[["elapsed"]]
  old_active <- statmod_variants$active
  old_override <- statmod_variants$override
  statmod_variants$active <- TRUE
  on.exit({
    statmod_variants$active <- old_active
    statmod_variants$override <- old_override
  }, add = TRUE)
  chosen <- list()
  set_override <- function(current = NULL, current_term = NULL) {
    ov <- list()
    for (a in auto) {
      tm <- if (identical(a$key, current)) current_term
            else if (!is.null(chosen[[a$key]])) chosen[[a$key]] else a$pilot
      ov[[a$param]][[a$name]] <- tm
    }
    statmod_variants$override <- ov
  }
  set_override()
  pilot_fit <- eval(cl, caller)
  rows <- list()
  best_fit <- NULL
  for (a in auto) {
    pt <- pilot_fit@spec@terms[[a$param]][[a$name]]
    cf <- stats::coef(pilot_fit)[[a$param]][modelterms7::term_coef_names(pt)]
    pv <- tryCatch(modelterms7::term_value(pt, coef = cf), error = function(e) NULL)
    vars <- modelterms7::term_variants(a$term, pilot = pv)
    if (is.null(names(vars))) names(vars) <- sprintf("version %d", seq_along(vars))
    best <- NULL
    for (vn in names(vars)) {
      set_override(a$key, vars[[vn]])
      f <- tryCatch(eval(cl, caller), error = function(e) e)
      bad <- inherits(f, "error")
      sc <- if (bad) NA_real_ else variant_score(f)
      rows[[length(rows) + 1L]] <- data.frame(
        param = a$param, term = a$name, version = vn, score = sc,
        chosen = FALSE, error = if (bad) conditionMessage(f) else "",
        stringsAsFactors = FALSE)
      if (is.finite(sc) && (is.null(best) || sc < best$score)) {
        best <- list(score = sc, fit = f, term = vars[[vn]], name = vn,
                     row = length(rows))
      }
    }
    if (is.null(best)) {
      stop(sprintf("no version of the term '%s' of %s could be fitted.",
                   a$name, a$param), call. = FALSE)
    }
    rows[[best$row]]$chosen <- TRUE
    chosen[[a$key]] <- best$term
    best_fit <- best$fit
  }
  best_fit@history$variants <- do.call(rbind, rows)
  best_fit@elapsed <- proc.time()[["elapsed"]] - t0
  best_fit
}
