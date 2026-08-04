#' Install or Update the Toolkit
#'
#' @description
#' Installs every member package from its GitHub repository, or reports what
#' would be installed.
#'
#' @details
#' The toolkit is not on CRAN, so the install path is GitHub and the work is
#' delegated to \pkg{pak}, which resolves the dependencies among the members
#' itself and installs only what is missing or out of date. The order of the
#' names therefore does not matter.
#'
#' With \code{action = "report"}, the default, nothing is installed and the
#' function returns the versions currently installed together with the call
#' that would update them. Installing packages is a side effect worth asking
#' for explicitly.
#'
#' @param action \code{"report"} to describe what is installed, or
#'   \code{"install"} to install and update.
#' @param quiet Passed to \pkg{pak}; suppresses its output when \code{TRUE}.
#'
#' @return For \code{"report"}, a data frame as
#'   \code{\link{statmodels7_versions}} returns, invisibly for
#'   \code{"install"}. The report is printed as a side effect.
#'
#' @seealso \code{\link{statmodels7_versions}}
#'
#' @examples
#' statmodels7_update()
#'
#' \dontrun{
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
