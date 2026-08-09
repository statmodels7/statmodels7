# Exports of the Toolkit That Mask One Another

The names exported by more than one member package, with the packages
that export them, in the order the search path resolves them.

## Usage

``` r
statmodels7_conflicts()
```

## Value

A named list, one entry per masked name, each a character vector of the
packages exporting it, most recently attached first. Empty when there is
nothing to report.

## Details

Attaching several packages puts several environments on the search path,
and a name exported by two of them resolves to whichever was attached
last. The report names the winner first, so that a caller can see which
function a bare name reaches. A name that one member exports and another
only registers a method on is not a conflict and does not appear:
methods dispatch, they do not mask.

Only members that are attached are examined, since a package that is
merely installed masks nothing.

The names come from each member's namespace rather than from its
attached environment. The two agree for an installed package and do not
under pkgload, which attaches a package's internal objects along with
its exports and adds shims of its own (`system.file`,
`library.dynam.unload`); reading the attached environment reports those
shims as a name every member exports, which is a conflict between
packages that export no such thing.

## See also

[`statmodels7_packages`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)

## Examples

``` r
statmodels7_conflicts()
#> named list()
```
