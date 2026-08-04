#' The Packages of the Toolkit
#'
#' @description
#' The names of the packages this one installs and attaches.
#'
#' @details
#' The list is read from the \code{Imports} field of this package's own
#' \code{DESCRIPTION}, keeping the entries whose names end in \code{7}, which
#' is the toolkit's naming convention. Reading it rather than writing it out
#' keeps one enumeration: a member added to \code{Imports} is a member here,
#' and the two cannot disagree. A base or third-party dependency does not end
#' in \code{7} and is therefore not reported as a member.
#'
#' @return A character vector of package names, sorted.
#'
#' @seealso \code{\link{statmodels7_versions}},
#'   \code{\link{statmodels7_conflicts}}
#'
#' @examples
#' statmodels7_packages()
#'
#' @export
statmodels7_packages <- function() {
  path <- system.file("DESCRIPTION", package = "statmodels7")
  imports <- read.dcf(path, fields = "Imports")[[1L]]
  if (is.na(imports)) return(character())
  pkgs <- trimws(strsplit(imports, ",")[[1L]])
  # a version constraint travels with the name and is not part of it
  pkgs <- trimws(sub("\\(.*", "", pkgs))
  sort(pkgs[nzchar(pkgs) & grepl("7$", pkgs)])
}


#' Installed Versions of the Toolkit
#'
#' @description
#' The version of each member package as currently installed.
#'
#' @return A data frame with columns \code{package} and \code{version}, the
#'   latter \code{NA} for a member that is not installed.
#'
#' @seealso \code{\link{statmodels7_packages}}, \code{\link{statmodels7_update}}
#'
#' @examples
#' statmodels7_versions()
#'
#' @export
statmodels7_versions <- function() {
  pkgs <- statmodels7_packages()
  version <- vapply(pkgs, function(p) {
    v <- tryCatch(as.character(utils::packageVersion(p)),
                  error = function(e) NA_character_)
    v
  }, character(1), USE.NAMES = FALSE)
  data.frame(package = pkgs, version = version, stringsAsFactors = FALSE)
}


#' The Body of the Attach Message
#'
#' @description
#' The member packages laid out in two columns with their versions, padded so
#' that the second column and the versions within it start at a fixed offset.
#'
#' @details
#' Padding is applied inside a line and removed at its end, since a trailing
#' run of spaces is invisible where it is written and visible wherever the
#' message is pasted. A member that is not installed is marked rather than
#' shown blank: the two are different states and only one of them is a
#' complete toolkit. Plain ASCII throughout, a non-ASCII character in a
#' package's own output being a portability question not worth answering for
#' a decoration.
#'
#' @param versions A data frame as \code{\link{statmodels7_versions}} returns.
#'
#' @return A character vector, one element per row of the layout, empty when
#'   there are no members.
#'
#' @seealso \code{\link{statmodels7_versions}}
#'
#' @keywords internal
statmodels7_attach_message <- function(versions) {
  n <- nrow(versions)
  if (n == 0L) return(character())
  w_pkg <- max(nchar(versions$package))
  w_ver <- max(nchar(ifelse(is.na(versions$version), "not installed",
                            versions$version)))
  entry <- sprintf("v %-*s %-*s", w_pkg, versions$package, w_ver,
                   ifelse(is.na(versions$version), "not installed",
                          versions$version))
  entry[is.na(versions$version)] <- sub("^v", "x", entry[is.na(versions$version)])

  rows <- ceiling(n / 2)
  left <- entry[seq_len(rows)]
  right <- c(entry[-seq_len(rows)], rep("", rows))[seq_len(rows)]
  sub("\\s+$", "", paste0(left, "   ", right))
}


#' Attach the Member Packages
#'
#' @description
#' Attaches every installed member that is not on the search path already.
#'
#' @details
#' Each member's own startup messages are suppressed, since the caller asked
#' for the toolkit and gets one message about it rather than one per package,
#' and so are the masking warnings, which
#' \code{\link{statmodels7_conflicts}} reports in one place instead.
#'
#' @return The names of the packages attached, invisibly.
#'
#' @seealso \code{\link{statmodels7_packages}}
#'
#' @keywords internal
statmodels7_attach <- function() {
  v <- statmodels7_versions()
  to_attach <- v$package[!is.na(v$version) & !paste0("package:", v$package) %in% search()]
  for (p in to_attach) {
    suppressPackageStartupMessages(
      library(p, character.only = TRUE, warn.conflicts = FALSE)
    )
  }
  invisible(to_attach)
}


.onAttach <- function(libname, pkgname) {
  v <- statmodels7_versions()
  statmodels7_attach()

  header <- sprintf("-- Attaching the statmodels7 toolkit %s",
                    utils::packageVersion("statmodels7"))
  msg <- c(header, statmodels7_attach_message(v))

  conf <- statmodels7_conflicts()
  if (length(conf)) {
    msg <- c(msg, "-- Conflicts", format_conflicts(conf))
  }
  packageStartupMessage(paste(msg, collapse = "\n"))
}
