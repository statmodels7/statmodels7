#' The Packages of the Toolkit
#'
#' @description
#' Names the eight member packages that `statmodels7` installs and attaches:
#' `basis7`, `distributions7`, `linkfunctions7`, `modelterms7`,
#' `numericals7`, `optimizers7`, `parameters7` and `penalties7`. Installing
#' `statmodels7` installs all eight, and attaching it attaches all eight;
#' this function is how the rest of the package, and a caller, learn which
#' those are.
#'
#' @details
#' # Where the list comes from
#'
#' The names are read out of the `Imports` field of this package's own
#' installed `DESCRIPTION`, keeping the entries that end in `7`. That is the
#' toolkit's naming convention: every member package is named for what it
#' supplies with a `7` appended, for S7. Version constraints such as
#' `distributions7 (>= 0.26.0)` are stripped, since a constraint travels
#' with a name and is not part of it.
#'
#' Reading the field keeps one enumeration. A member added to `Imports` is a
#' member here on the next install, with no second list to keep in step.
#'
#' # Why `S7` itself is excluded
#'
#' `S7` ends in a `7` and is not a member. Its digit counts R's object
#' systems, and it is the system every member is built on. It has to be
#' declared in `Imports` because the code says `S7::` throughout, so the
#' name reaches this field and is removed by name.
#'
#' @return A character vector of package names, sorted alphabetically. Length
#'   eight for a complete installation. The result does not say whether a
#'   member is installed; [statmodels7_versions()] answers that.
#'
#' @seealso [statmodels7_versions()] for the installed version of each,
#'   [statmodels7_conflicts()] for the exports they mask between them,
#'   [statmodels7_update()] to install the ones that are out of date.
#'
#' @examples
#' pkgs <- statmodels7_packages()
#' pkgs
#'
#' # Every name ends in 7, and S7 is not among them.
#' all(grepl("7$", pkgs))
#' "S7" %in% pkgs
#'
#' @export
statmodels7_packages <- function() {
  path <- system.file("DESCRIPTION", package = "statmodels7")
  imports <- read.dcf(path, fields = "Imports")[[1L]]
  if (is.na(imports)) return(character())
  pkgs <- trimws(strsplit(imports, ",")[[1L]])
  # a version constraint travels with the name and is not part of it
  pkgs <- trimws(sub("\\(.*", "", pkgs))
  sort(pkgs[nzchar(pkgs) & grepl("7$", pkgs) & pkgs != "S7"])
}


#' Installed Versions of the Toolkit
#'
#' @description
#' Reports the version of each member package as it is installed on this
#' machine, one row per member. A member that is not installed reports `NA`
#' instead of being dropped, so the result always has one row per name
#' [statmodels7_packages()] returns and a missing member is visible.
#'
#' @details
#' The version is read with [utils::packageVersion()] and coerced to a
#' character string, so `"0.38.0"` sorts and prints as written and no
#' comparison is implied between two members' numbers. The members version
#' independently.
#'
#' This is the source of the table the attach message prints, and of the
#' comparison [statmodels7_update()] makes against GitHub.
#'
#' @return A data frame with one row per member and two character columns:
#'   \describe{
#'     \item{`package`}{the member's name, sorted alphabetically.}
#'     \item{`version`}{its installed version, or `NA` when the package
#'       cannot be found in the library path.}
#'   }
#'   Row names are the default integers, and `stringsAsFactors` is `FALSE`,
#'   so both columns are character.
#'
#' @seealso [statmodels7_packages()] for the names alone,
#'   [statmodels7_update()] to install the members that are behind.
#'
#' @examples
#' v <- statmodels7_versions()
#' v
#'
#' # One row per member, and every column is character.
#' nrow(v) == length(statmodels7_packages())
#' vapply(v, class, character(1))
#'
#' # Which members, if any, are missing from this library.
#' v$package[is.na(v$version)]
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
#' Lays the member packages out in two columns with their versions, padded so
#' that the second column and the versions within it start at a fixed offset.
#' Each entry is prefixed `v` when the member is installed and `x` when it is
#' not, and the two columns are filled down the left one first.
#'
#' @details
#' Padding is applied inside a line and stripped at its end, since a trailing
#' run of spaces is invisible where it is written and visible wherever the
#' message is pasted. A member that is not installed is marked `x` and named
#' with the text `not installed`; leaving it blank would hide the difference
#' between an absent member and a present one.
#'
#' Plain ASCII throughout. A non-ASCII character in a package's startup
#' message is a portability risk taken for a decoration.
#'
#' @param versions A data frame as [statmodels7_versions()] returns, with
#'   character columns `package` and `version` and `NA` in `version` for a
#'   member that is not installed. Not validated: this is called from
#'   `.onAttach()` with that function's own result.
#'
#' @return A character vector with `ceiling(nrow(versions) / 2)` elements, one
#'   per row of the two-column layout, each already trimmed of trailing
#'   spaces. `character(0)` when `versions` has no rows.
#'
#' @seealso [statmodels7_versions()] for the input,
#'   [format_conflicts()] for the other half of the same message.
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
#' Attaches every installed member that is not on the search path already,
#' with [library()]. A member already attached is left alone, so calling this
#' twice attaches nothing the second time, and a member that is not installed
#' is skipped without an error.
#'
#' @details
#' Each member's own startup messages are suppressed: the caller asked for
#' the toolkit and gets one message about it, not eight. The masking warnings
#' are suppressed too, and [statmodels7_conflicts()] reports the same
#' information once, at the end, in one block.
#'
#' Attachment order follows [statmodels7_packages()], which is alphabetical.
#' The members export no name in common, so nothing about the result depends
#' on that order; [statmodels7_conflicts()] is what would report it if one
#' day they did.
#'
#' @return The names of the packages this call attached, a character vector,
#'   invisibly. Empty when every installed member was already on the search
#'   path.
#'
#' @seealso [statmodels7_packages()] for the members,
#'   [statmodels7_conflicts()] for what the suppressed warnings would have
#'   said.
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
