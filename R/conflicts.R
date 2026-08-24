#' Exports of the Toolkit That Mask One Another
#'
#' @description
#' Reports the names exported by more than one attached member package,
#' together with the packages exporting them, ordered so that the package a
#' bare name reaches comes first. Attaching several packages puts several
#' environments on the search path, and a name exported by two of them
#' resolves to whichever was attached last; this is how a caller sees which
#' function that is. As of `statmodels7` 0.88.1 the eight members export no
#' name in common, so the result is empty.
#'
#' @details
#' # What counts as a conflict
#'
#' Two members exporting the same name. A name that one member exports and
#' another only registers an S7 or S3 method on is not a conflict and does
#' not appear: methods dispatch on the class of the argument, they do not
#' mask.
#'
#' Only members that are attached are examined. A package that is merely
#' installed puts nothing on the search path and masks nothing, so it cannot
#' contribute, and with fewer than two members attached the result is empty
#' by construction.
#'
#' # What it does not cover
#'
#' The comparison is **between the members**. `statmodels7`'s own exports are
#' not compared against theirs, so a name this package and a member both
#' export would go unreported here. There is no such name today, and the
#' example below is the one-line check.
#'
#' # Namespaces, not attached environments
#'
#' The names are read with [getNamespaceExports()], and not off the attached
#' environment. The two agree for an installed package. They diverge under
#' \pkg{pkgload}, which attaches a package's internal objects along with its
#' exports and adds shims of its own, `system.file` and
#' `library.dynam.unload` among them. Reading the attached environment there
#' reports those shims as names every member exports, a conflict between
#' packages that export no such thing.
#'
#' @return A named list, one entry per masked name, sorted by name. Each entry
#'   is a character vector of the packages exporting that name, the most
#'   recently attached first, which is the one a bare call reaches. An empty
#'   named list when nothing is masked, so `length()` is the count of
#'   conflicts and `names()` are the masked names.
#'
#' @seealso [statmodels7_packages()] for the members examined,
#'   [format_conflicts()] for the one-line-per-name rendering the attach
#'   message uses.
#'
#' @examples
#' # Empty: the eight members are disjoint in what they export.
#' statmodels7_conflicts()
#' length(statmodels7_conflicts())
#'
#' # The check this function does not do, spelled out.
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
#' name, in the shape the attach message prints:
#' `x pkgA::name masks pkgB::name`. The winner is named on the left of
#' `masks` and every package it hides on the right, comma separated.
#'
#' @param conflicts A named list as [statmodels7_conflicts()] returns: one
#'   entry per masked name, each a character vector of at least two package
#'   names with the winner first. An empty list gives an empty result. Not
#'   validated; a vector of length one would render as masking nothing.
#'
#' @return A character vector with one element per entry of `conflicts`, in
#'   the order they arrived, with no names. `character(0)` for an empty
#'   input, so a caller can `c()` the result into a message
#'   unconditionally.
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
