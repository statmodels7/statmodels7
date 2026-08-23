#' Exports of the Toolkit That Mask One Another
#'
#' @description
#' The names exported by more than one member package, with the packages that
#' export them, in the order the search path resolves them.
#'
#' @details
#' Attaching several packages puts several environments on the search path,
#' and a name exported by two of them resolves to whichever was attached last.
#' The report names the winner first, so that a caller can see which function
#' a bare name reaches. A name that one member exports and another only
#' registers a method on is not a conflict and does not appear: methods
#' dispatch, they do not mask.
#'
#' Only members that are attached are examined, since a package that is merely
#' installed masks nothing.
#'
#' The names come from each member's namespace rather than from its attached
#' environment. The two agree for an installed package and do not under
#' \pkg{pkgload}, which attaches a package's internal objects along with its
#' exports and adds shims of its own (`system.file`,
#' `library.dynam.unload`); reading the attached environment reports those
#' shims as a name every member exports, which is a conflict between packages
#' that export no such thing.
#'
#' @return A named list, one entry per masked name, each a character vector of
#'   the packages exporting it, most recently attached first. Empty when there
#'   is nothing to report.
#'
#' @seealso [statmodels7_packages()]
#'
#' @examples
#' statmodels7_conflicts()
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
#' Formats the result of [statmodels7_conflicts()] as one line per
#' masked name.
#'
#' @param conflicts A list as returned by [statmodels7_conflicts()].
#'
#' @return A character vector, one element per masked name.
#'
#' @keywords internal
format_conflicts <- function(conflicts) {
  vapply(names(conflicts), function(nm) {
    who <- conflicts[[nm]]
    sprintf("x %s::%s masks %s", who[1L], nm,
            paste(sprintf("%s::%s", who[-1L], nm), collapse = ", "))
  }, character(1), USE.NAMES = FALSE)
}
