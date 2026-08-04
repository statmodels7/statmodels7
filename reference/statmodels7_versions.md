# Installed Versions of the Toolkit

The version of each member package as currently installed.

## Usage

``` r
statmodels7_versions()
```

## Value

A data frame with columns `package` and `version`, the latter `NA` for a
member that is not installed.

## See also

[`statmodels7_packages`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md),
[`statmodels7_update`](https://statmodels7.github.io/statmodels7/reference/statmodels7_update.md)

## Examples

``` r
statmodels7_versions()
#>          package    version
#> 1         basis7      0.3.1
#> 2 distributions7 0.0.0.9000
#> 3 linkfunctions7 0.0.0.9000
#> 4    optimizers7 0.0.0.9000
#> 5    parameters7      0.3.0
```
