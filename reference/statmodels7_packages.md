# The Packages of the Toolkit

Names the eight member packages that `statmodels7` installs and
attaches: `basis7`, `distributions7`, `linkfunctions7`, `modelterms7`,
`numericals7`, `optimizers7`, `parameters7` and `penalties7`. Installing
`statmodels7` installs all eight, and attaching it attaches all eight;
this function is how the rest of the package, and a caller, learn which
those are.

## Usage

``` r
statmodels7_packages()
```

## Value

A character vector of package names, sorted alphabetically. Length eight
for a complete installation. The result does not say whether a member is
installed;
[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
answers that.

## Where the list comes from

The names are read out of the `Imports` field of this package's own
installed `DESCRIPTION`, keeping the entries that end in `7`. That is
the toolkit's naming convention: every member package is named for what
it supplies with a `7` appended, for S7. Version constraints such as
`distributions7 (>= 0.26.0)` are stripped, since a constraint travels
with a name and is not part of it.

Reading the field keeps one enumeration. A member added to `Imports` is
a member here on the next install, with no second list to keep in step.

## Why `S7` itself is excluded

`S7` ends in a `7` and is not a member. Its digit counts R's object
systems, and it is the system every member is built on. It has to be
declared in `Imports` because the code says `S7::` throughout, so the
name reaches this field and is removed by name.

## See also

[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
for the installed version of each,
[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
for the exports they mask between them,
[`statmodels7_update()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_update.md)
to install the ones that are out of date.

## Examples

``` r
pkgs <- statmodels7_packages()
pkgs
#> [1] "basis7"         "distributions7" "linkfunctions7" "modelterms7"   
#> [5] "numericals7"    "optimizers7"    "parameters7"    "penalties7"    

# Every name ends in 7, and S7 is not among them.
all(grepl("7$", pkgs))
#> [1] TRUE
"S7" %in% pkgs
#> [1] FALSE
```
