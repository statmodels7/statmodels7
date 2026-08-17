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
#>          package version
#> 1         basis7   0.4.1
#> 2 distributions7  0.24.0
#> 3 linkfunctions7   0.1.0
#> 4    modelterms7  0.49.0
#> 5    numericals7   0.7.0
#> 6    optimizers7   0.4.0
#> 7    parameters7  0.11.0
#> 8     penalties7  0.15.0
```
