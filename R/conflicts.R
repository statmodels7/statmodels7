#' Exports of the Toolkit That Mask One Another
#'
#' @description
#' Reports the names exported by more than one attached member package,
#' together with the packages exporting them, ordered so that the package that
#' a bare name reaches comes first. Attaching several packages puts several
#' environments on the search path, and a name exported by two of them
#' resolves to whichever was attached last; this function shows which
#' function that is. The members currently export no name in common, so the
#' result is empty.
#'
#' @details
#' # Definition of a conflict
#'
#' A conflict is a name exported by two members. A name that one member
#' exports and on which another member only registers an S7 or S3 method is
#' not a conflict and does not appear, since methods dispatch on the class of
#' the argument and do not mask.
#'
#' Only members that are attached are examined. A package that is merely
#' installed puts nothing on the search path and masks nothing, so it cannot
#' contribute. If fewer than two members are attached, the result is empty by
#' construction.
#'
#' # Scope of the comparison
#'
#' The comparison is between the members. The exports of `statmodels7`
#' itself are not compared with theirs, so a name that this package and a
#' member both export would not be reported here. At present every name that
#' this package exports differs from the names that the members export, and
#' the example below performs the check.
#'
#' # Namespaces and attached environments
#'
#' The names are read with [getNamespaceExports()] and not from the attached
#' environment. The two agree for an installed package. They differ under
#' \pkg{pkgload}, which attaches the internal objects of a package along with
#' its exports and adds shims of its own, `system.file` and
#' `library.dynam.unload` among them. Reading the attached environment there
#' would report those shims as names that every member exports, a conflict
#' between packages that export no such thing.
#'
#' @return A named list, one entry per masked name, sorted by name. Each entry
#'   is a character vector of the packages exporting that name, the most
#'   recently attached first, which is the one that a bare call reaches. If no
#'   name is masked, the result is an empty named list, so `length()` gives
#'   the number of conflicts and `names()` gives the masked names.
#'
#' @seealso [statmodels7_packages()] for the members examined,
#'   [format_conflicts()] for the one-line-per-name rendering that the attach
#'   message uses.
#'
#' @examples
#' # Empty: the members export disjoint sets of names.
#' statmodels7_conflicts()
#' length(statmodels7_conflicts())
#'
#' # The following check, which this function does not perform, counts for
#' # each member the names shared with the exports of statmodels7 itself.
#' own <- getNamespaceExports("statmodels7")
#' vapply(statmodels7_packages(),
#'        function(p) length(intersect(own, getNamespaceExports(p))),
#'        integer(1))
#'
#' @export
statmodels7_conflicts <- function() {
  pkgs <- statmodels7_packages()
  attached <- paste0("package:", pkgs)
  # The search path is ordered from most to least recently attached, which is
  # exactly the order a name resolves in.
  on_path <- intersect(search(), attached)
  if (length(on_path) < 2L) return(structure(list(), names = character()))

  nms <- sub("^package:", "", on_path)
  exports <- lapply(nms, getNamespaceExports)
  names(exports) <- nms

  all_names <- unlist(exports, use.names = FALSE)
  dup <- unique(all_names[duplicated(all_names)])
  if (!length(dup)) return(structure(list(), names = character()))

  out <- lapply(dup, function(nm) {
    names(exports)[vapply(exports, function(e) nm %in% e, logical(1))]
  })
  names(out) <- dup
  out[order(names(out))]
}


#' Render a Conflict Report
#'
#' @description
#' Formats the result of [statmodels7_conflicts()] as one line per masked
#' name, in the shape that the attach message prints:
#' `x pkgA::name masks pkgB::name`. The winner is named on the left of
#' `masks`, and every package that it hides is named on the right, comma
#' separated.
#'
#' @param conflicts A named list as [statmodels7_conflicts()] returns: one
#'   entry per masked name, each a character vector of at least two package
#'   names with the winner first. An empty list gives an empty result. The
#'   argument is not validated, and a vector of length one would render as
#'   masking nothing.
#'
#' @return A character vector with one element per entry of `conflicts`, in
#'   the order in which they arrived, with no names. The result is
#'   `character(0)` for an empty input, so a caller can `c()` the result into
#'   a message unconditionally.
#'
#' @seealso [statmodels7_conflicts()] for the input,
#'   [statmodels7_attach_message()] for the other half of the same message.
#'
#' @keywords internal
format_conflicts <- function(conflicts) {
  vapply(names(conflicts), function(nm) {
    who <- conflicts[[nm]]
    sprintf("x %s::%s masks %s", who[1L], nm,
            paste(sprintf("%s::%s", who[-1L], nm), collapse = ", "))
  }, character(1), USE.NAMES = FALSE)
}
