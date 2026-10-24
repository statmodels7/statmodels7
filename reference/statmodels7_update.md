# Install or Update the Toolkit

Installs the member packages of the toolkit from their GitHub
repositories under `github.com/statmodels7`, or prints what is installed
and stops. The default is to print: installing packages is a side
effect, so the installation runs only when `action = "install"` is
given.

## Usage

``` r
statmodels7_update(action = c("report", "install"))
```

## Arguments

- action:

  `"report"` (the default) to print the installed versions and the calls
  that would update them, or `"install"` to install and update. The
  argument is matched with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html), so an
  unambiguous prefix such as `"i"` is accepted and anything else is an
  error.

## Value

A data frame as
[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
returns, with columns `package` and `version`. The data frame is
returned invisibly under both actions: the report prints to the console
with [`cat()`](https://rdrr.io/r/base/cat.html) and is not the return
value. Under `"install"` the versions are read after the installation,
so the result is the state that the call left behind.

## Where the packages come from

The members of the toolkit are not on CRAN, so the install path is
GitHub, at `statmodels7/<package>` for each member that
[`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
names. The work is handed to pak, which resolves the dependencies among
the members itself and installs only what is missing or behind. The
order in which the names are passed therefore does not matter.

pak is not a dependency of this package. When `action = "install"` is
given and pak is not installed, the call stops with a message that
contains the `install.packages("pak")` line and the equivalent
`remotes::install_github()` call written out for every member
repository, so that it can be pasted where remotes is used instead.

## Installing over a loaded package

A member that is attached in the current session holds its own DLL open
on Windows, and pak cannot replace a file that is in use. R should be
restarted before the installation. Otherwise the installation of a
compiled member can fail partway and leave that member in a state in
which
[`packageVersion()`](https://rdrr.io/r/utils/packageDescription.html)
reports it as absent.

## Printing

Both actions print, and an argument cannot suppress the output. The
return value of the report is the result of
[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md),
so a caller who wants the versions without the console output calls that
function instead. The installation is carried out by pak, which prints
its own progress and has no `quiet` argument.

## See also

[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
for the versions alone,
[`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
for the names of the members.

## Examples

``` r
# Prints the table and the two calls that would update it. Installs
# nothing.
v <- statmodels7_update()
#> statmodels7 toolkit, installed versions:
#>   basis7           0.15.2
#>   distributions7   0.87.0
#>   linkfunctions7   0.5.6
#>   modelterms7      0.107.0
#>   numericals7      0.23.5
#>   optimizers7      0.15.2
#>   parameters7      0.26.1
#>   penalties7       0.33.1
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
