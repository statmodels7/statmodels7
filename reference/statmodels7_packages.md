# The Packages of the Toolkit

Names the member packages that `statmodels7` installs and attaches:
`basis7`, `distributions7`, `linkfunctions7`, `modelterms7`,
`numericals7`, `optimizers7`, `parameters7` and `penalties7`. Installing
`statmodels7` installs all of them, and attaching it attaches all of
them. The rest of the package and the caller use this function to obtain
the list.

## Usage

``` r
statmodels7_packages()
```

## Value

A character vector of package names, sorted alphabetically, with one
element per member. The result does not say whether a member is
installed;
[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
reports that.

## Source of the list

The names are read from the `Imports` field of the installed
`DESCRIPTION` of this package, keeping the entries that end in `7`. This
follows the naming convention of the toolkit: every member name ends in
`7`, for S7. Version constraints such as `distributions7 (>= 0.26.0)`
are stripped, since a constraint is attached to a name and is not part
of it.

Reading the field keeps one enumeration. A member added to `Imports` is
a member here after the next install, and no second list has to be kept
in step.

## The exclusion of `S7`

`S7` ends in a `7` but is not a member: it is the object system on which
every member is built. It has to be declared in `Imports` because the
code calls `S7::` throughout, so the name reaches this field and is
removed by name.

## See also

[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
for the installed version of each,
[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
for the exports that more than one member makes,
[`statmodels7_update()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_update.md)
to install or update the members.

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
