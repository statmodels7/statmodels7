# The Packages of the Toolkit

The names of the packages this one installs and attaches.

## Usage

``` r
statmodels7_packages()
```

## Value

A character vector of package names, sorted.

## Details

The list is read from the `Imports` field of this package's own
`DESCRIPTION`, keeping the entries whose names end in `7`, which is the
toolkit's naming convention. Reading it rather than writing it out keeps
one enumeration: a member added to `Imports` is a member here, and the
two cannot disagree. A base or third-party dependency does not end in
`7` and is therefore not reported as a member.

## See also

[`statmodels7_versions`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md),
[`statmodels7_conflicts`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)

## Examples

``` r
statmodels7_packages()
#> [1] "basis7"         "distributions7" "linkfunctions7" "numericals7"   
#> [5] "optimizers7"    "parameters7"    "penalties7"    
```
