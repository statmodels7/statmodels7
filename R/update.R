#' Install or Update the Toolkit
#'
#' @description
#' Installs the member packages of the toolkit from their GitHub
#' repositories under `github.com/statmodels7`, or prints what is installed
#' and stops. The default is to print: installing packages is a side effect,
#' so the installation runs only when `action = "install"` is given.
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
#' is given and \pkg{pak} is not installed, the call stops with a message
#' that contains the `install.packages("pak")` line and the equivalent
#' [remotes::install_github()] call written out for every member repository,
#' so that it can be pasted where \pkg{remotes} is used instead.
#'
#' # Installing over a loaded package
#'
#' A member that is attached in the current session holds its own DLL open on
#' Windows, and \pkg{pak} cannot replace a file that is in use. R should be
#' restarted before the installation. Otherwise the installation of a compiled
#' member can fail partway and leave that member in a state in which
#' `packageVersion()` reports it as absent.
#'
#' # Printing
#'
#' Both actions print, and no argument suppresses the output. The return value
#' of the report is the result of [statmodels7_versions()], so a caller who
#' wants the versions without the console output calls that function instead.
#' The installation is carried out by \pkg{pak}, which prints its own progress
#' and has no `quiet` argument.
#'
#' @param action `"report"` (the default) to print the installed versions and
#'   the calls that would update them, or `"install"` to install and update.
#'   Matched with [match.arg()], so an unambiguous prefix such as `"i"` is
#'   accepted and anything else is an error.
#'
#' @return A data frame as [statmodels7_versions()] returns, with columns
#'   `package` and `version`. Returned invisibly under both actions: the
#'   report prints to the console with [cat()] and is not the return value.
#'   Under `"install"` the versions are read after the installation, so the
#'   result is the state that the call left behind.
#'
#' @seealso [statmodels7_versions()] for the versions alone,
#'   [statmodels7_packages()] for the names of the members.
#'
#' @examples
#' # Prints the table and the two calls that would update it. Installs
#' # nothing.
#' v <- statmodels7_update()
#'
#' # The return value is the versions, so that is the call to make when the
#' # table is wanted and the printing is not.
#' identical(v, statmodels7_versions())
#'
#' \dontrun{
#' # Installs from GitHub, so it needs a network connection and writes to
#' # the library. Restart R first if a member is attached.
#' statmodels7_update("install")
#' }
#'
#' @export
statmodels7_update <- function(action = c("report", "install")) {
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
