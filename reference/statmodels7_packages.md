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
two cannot disagree.

`S7` is excluded, and it is the one name the convention cannot tell
apart on its own: its `7` counts R's object systems and not this
toolkit's, and it is what every member is built on rather than a member.
Declaring it in `Imports`, which `R CMD check` requires because the code
says `S7::` throughout, therefore made it report as a ninth package
until this line existed.

## See also

[`statmodels7_versions`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md),
[`statmodels7_conflicts`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)

## Examples

``` r
statmodels7_packages()
#> [1] "basis7"         "distributions7" "linkfunctions7" "modelterms7"   
#> [5] "numericals7"    "optimizers7"    "parameters7"    "penalties7"    
```
