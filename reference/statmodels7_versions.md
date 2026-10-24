# Installed Versions of the Toolkit

Reports the version of each member package as it is installed on this
machine, one row per member. A member that is not installed has `NA` as
its version and is kept, so the result always has one row per name that
[`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
returns and a missing member is visible.

## Usage

``` r
statmodels7_versions()
```

## Value

A data frame with one row per member and two character columns:

- `package`:

  the member's name, sorted alphabetically.

- `version`:

  its installed version, or `NA` when the package cannot be found in the
  library path.

Row names are the default integers, and `stringsAsFactors` is `FALSE`,
so both columns are character.

## Details

The version is read with
[`utils::packageVersion()`](https://rdrr.io/r/utils/packageDescription.html)
and coerced to a character string, so a version such as `"0.38.0"` is
kept as written. The members are versioned independently, and the
version numbers of two members are not comparable.

This is the source of the table that the attach message prints and of
the installed versions that
[`statmodels7_update()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_update.md)
reports.

## See also

[`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
for the names alone,
[`statmodels7_update()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_update.md)
to install or update the members.

## Examples

``` r
v <- statmodels7_versions()
v
#>          package version
#> 1         basis7  0.15.2
#> 2 distributions7  0.87.0
#> 3 linkfunctions7   0.5.6
#> 4    modelterms7 0.107.0
#> 5    numericals7  0.23.5
#> 6    optimizers7  0.15.2
#> 7    parameters7  0.26.1
#> 8     penalties7  0.33.1

# One row per member, and every column is character.
nrow(v) == length(statmodels7_packages())
#> [1] TRUE
vapply(v, class, character(1))
#>     package     version 
#> "character" "character" 

# Which members, if any, are missing from this library.
v$package[is.na(v$version)]
#> character(0)
```
