# Install or Update the Toolkit

Installs every member package from its GitHub repository, or reports
what would be installed.

## Usage

``` r
statmodels7_update(action = c("report", "install"), quiet = FALSE)
```

## Arguments

- action:

  `"report"` to describe what is installed, or `"install"` to install
  and update.

- quiet:

  Passed to pak; suppresses its output when `TRUE`.

## Value

For `"report"`, a data frame as
[`statmodels7_versions`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
returns, invisibly for `"install"`. The report is printed as a side
effect.

## Details

The toolkit is not on CRAN, so the install path is GitHub and the work
is delegated to pak, which resolves the dependencies among the members
itself and installs only what is missing or out of date. The order of
the names therefore does not matter.

With `action = "report"`, the default, nothing is installed and the
function returns the versions currently installed together with the call
that would update them. Installing packages is a side effect worth
asking for explicitly.

## See also

[`statmodels7_versions`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)

## Examples

``` r
statmodels7_update()
#> statmodels7 toolkit, installed versions:
#>   basis7           0.3.1
#>   distributions7   0.0.0.9000
#>   linkfunctions7   0.0.0.9000
#>   optimizers7      0.0.0.9000
#>   parameters7      0.3.0
#> 
#> To install or update every member:
#>   statmodels7_update("install")
#>   pak::pak(c("statmodels7/basis7", "statmodels7/distributions7", "statmodels7/linkfunctions7", "statmodels7/optimizers7", "statmodels7/parameters7"))

if (FALSE) { # \dontrun{
statmodels7_update("install")
} # }
```
