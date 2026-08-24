#' Install or Update the Toolkit
#'
#' @description
#' Installs the eight member packages of the toolkit from their GitHub
#' repositories under `github.com/statmodels7`, or prints what is installed
#' and stops. The default is to print: installing packages is a side effect,
#' and this function asks for it to be requested by name.
#'
#' @details
#' # Where the packages come from
#'
#' None of the toolkit is on CRAN, so the install path is GitHub, at
#' `statmodels7/<package>` for each member [statmodels7_packages()] names.
#' The work is handed to \pkg{pak}, which resolves the dependencies among the
#' members itself and installs only what is missing or behind. The order the
#' names are passed in therefore does not matter.
#'
#' \pkg{pak} is not a dependency of this package. When `action = "install"`
#' is asked for and \pkg{pak} is not installed, the call stops with a message
#' giving both the `install.packages("pak")` line and the equivalent
#' [remotes::install_github()] call spelled out with all eight repositories,
#' so a caller who prefers \pkg{remotes} can paste it.
#'
#' # Installing over a loaded package
#'
#' A member that is attached in the current session holds its own DLL open on
#' Windows, and \pkg{pak} cannot replace a file that is in use. Restart R
#' before installing, or the install of a compiled member fails partway and
#' leaves that member in a state where `packageVersion()` reports it absent.
#'
#' @param action `"report"` (the default) to print the installed versions and
#'   the calls that would update them, or `"install"` to install and update.
#'   Matched with [match.arg()], so an unambiguous prefix such as `"i"` is
#'   accepted and anything else is an error.
#' @param quiet Accepted for symmetry with the installers and currently not
#'   read: \pkg{pak} is called with its own defaults. A single logical.
#'
#' @return A data frame as [statmodels7_versions()] returns, with columns
#'   `package` and `version`. Returned invisibly under both actions: the
#'   report prints to the console with [cat()] and is not the return value.
#'   Under `"install"` the versions are read **after** the install, so the
#'   result is the state the call left behind.
#'
#' @seealso [statmodels7_versions()] for the versions alone,
#'   [statmodels7_packages()] for the members installed.
#'
#' @examples
#' # Prints the table and the two calls that would update it. Installs
#' # nothing.
#' v <- statmodels7_update()
#' identical(v, statmodels7_versions())
#'
#' \dontrun{
#' # Installs from GitHub, so it needs a network connection and writes to
#' # the library. Restart R first if a member is attached.
#' statmodels7_update("install")
#' }
#'
#' @export
statmodels7_update <- function(action = c("report", "install"),
                               quiet = FALSE) {
  action <- match.arg(action)
  pkgs <- statmodels7_packages()
  spec <- paste0("statmodels7/", pkgs)

  if (action == "report") {
    v <- statmodels7_versions()
    missing <- v$package[is.na(v$version)]
    cat("statmodels7 toolkit, installed versions:\n")
    for (i in seq_len(nrow(v))) {
      cat(sprintf("  %-16s %s\n", v$package[i],
                  if (is.na(v$version[i])) "not installed" else v$version[i]))
    }
    if (length(missing)) {
      cat(sprintf("\n%d not installed.\n", length(missing)))
    }
    cat("\nTo install or update every member:\n")
    cat("  statmodels7_update(\"install\")\n")
    cat(sprintf("  pak::pak(c(%s))\n",
                paste0("\"", spec, "\"", collapse = ", ")))
    return(invisible(v))
  }

  if (!requireNamespace("pak", quietly = TRUE)) {
    stop(paste0(
      "'pak' is needed to install from GitHub and is not available.\n",
      "  install.packages(\"pak\"), or run the equivalent yourself:\n",
      "  remotes::install_github(c(",
      paste0("\"", spec, "\"", collapse = ", "), "))"
    ), call. = FALSE)
  }
  pak::pak(spec, ask = FALSE)
  invisible(statmodels7_versions())
}
