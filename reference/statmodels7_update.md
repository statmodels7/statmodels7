# Install or Update the Toolkit

Installs the eight member packages of the toolkit from their GitHub
repositories under `github.com/statmodels7`, or prints what is installed
and stops. The default is to print: installing packages is a side
effect, and this function asks for it to be requested by name.

## Usage

``` r
statmodels7_update(action = c("report", "install"))
```

## Arguments

- action:

  `"report"` (the default) to print the installed versions and the calls
  that would update them, or `"install"` to install and update. Matched
  with [`match.arg()`](https://rdrr.io/r/base/match.arg.html), so an
  unambiguous prefix such as `"i"` is accepted and anything else is an
  error.

## Value

A data frame as
[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
returns, with columns `package` and `version`. Returned invisibly under
both actions: the report prints to the console with
[`cat()`](https://rdrr.io/r/base/cat.html) and is not the return value.
Under `"install"` the versions are read **after** the install, so the
result is the state the call left behind.

## Where the packages come from

None of the toolkit is on CRAN, so the install path is GitHub, at
`statmodels7/<package>` for each member
[`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
names. The work is handed to pak, which resolves the dependencies among
the members itself and installs only what is missing or behind. The
order the names are passed in therefore does not matter.

pak is not a dependency of this package. When `action = "install"` is
asked for and pak is not installed, the call stops with a message giving
both the `install.packages("pak")` line and the equivalent
`remotes::install_github()` call spelled out with all eight
repositories, so a caller who prefers remotes can paste it.

## Installing over a loaded package

A member that is attached in the current session holds its own DLL open
on Windows, and pak cannot replace a file that is in use. Restart R
before installing, or the install of a compiled member fails partway and
leaves that member in a state where
[`packageVersion()`](https://rdrr.io/r/utils/packageDescription.html)
reports it absent.

## Printing

Both actions print, and neither takes an argument to suppress it. The
report exists to be read, and its return value is
[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md),
so a caller who wants the versions without the console output calls that
instead. The install hands the work to pak, which prints its own
progress and takes no `quiet` argument.

## See also

[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
for the versions alone,
[`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
for the members installed.

## Examples

``` r
# Prints the table and the two calls that would update it. Installs
# nothing.
v <- statmodels7_update()
#> statmodels7 toolkit, installed versions:
#>   basis7           0.6.0
#>   distributions7   0.43.0
#>   linkfunctions7   0.3.0
#>   modelterms7      0.67.0
#>   numericals7      0.12.0
#>   optimizers7      0.8.0
#>   parameters7      0.20.0
#>   penalties7       0.20.0
#> 
#> To install or update every member:
#>   statmodels7_update("install")
#>   pak::pak(c("statmodels7/basis7", "statmodels7/distributions7", "statmodels7/linkfunctions7", "statmodels7/modelterms7", "statmodels7/numericals7", "statmodels7/optimizers7", "statmodels7/parameters7", "statmodels7/penalties7"))

# The return value is the versions, so that is the call to make when the
# table is wanted and the printing is not.
identical(v, statmodels7_versions())
#> [1] TRUE

if (FALSE) { # \dontrun{
# Installs from GitHub, so it needs a network connection and writes to
# the library. Restart R first if a member is attached.
statmodels7_update("install")
} # }
```
