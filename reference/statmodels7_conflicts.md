# Exports of the Toolkit That Mask One Another

Reports the names exported by more than one attached member package,
together with the packages exporting them, ordered so that the package a
bare name reaches comes first. Attaching several packages puts several
environments on the search path, and a name exported by two of them
resolves to whichever was attached last; this is how a caller sees which
function that is. As of `statmodels7` 0.88.1 the eight members export no
name in common, so the result is empty.

## Usage

``` r
statmodels7_conflicts()
```

## Value

A named list, one entry per masked name, sorted by name. Each entry is a
character vector of the packages exporting that name, the most recently
attached first, which is the one a bare call reaches. An empty named
list when nothing is masked, so
[`length()`](https://rdrr.io/r/base/length.html) is the count of
conflicts and [`names()`](https://rdrr.io/r/base/names.html) are the
masked names.

## What counts as a conflict

Two members exporting the same name. A name that one member exports and
another only registers an S7 or S3 method on is not a conflict and does
not appear: methods dispatch on the class of the argument, they do not
mask.

Only members that are attached are examined. A package that is merely
installed puts nothing on the search path and masks nothing, so it
cannot contribute, and with fewer than two members attached the result
is empty by construction.

## What it does not cover

The comparison is **between the members**. `statmodels7`'s own exports
are not compared against theirs, so a name this package and a member
both export would go unreported here. There is no such name today, and
the example below is the one-line check.

## Namespaces, not attached environments

The names are read with
[`getNamespaceExports()`](https://rdrr.io/r/base/ns-reflect.html), and
not off the attached environment. The two agree for an installed
package. They diverge under pkgload, which attaches a package's internal
objects along with its exports and adds shims of its own, `system.file`
and `library.dynam.unload` among them. Reading the attached environment
there reports those shims as names every member exports, a conflict
between packages that export no such thing.

## See also

[`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
for the members examined,
[`format_conflicts()`](https://statmodels7.github.io/statmodels7/reference/format_conflicts.md)
for the one-line-per-name rendering the attach message uses.

## Examples

``` r
# Empty: the eight members are disjoint in what they export.
statmodels7_conflicts()
#> named list()
length(statmodels7_conflicts())
#> [1] 0

# The check this function does not do, spelled out.
own <- getNamespaceExports("statmodels7")
vapply(statmodels7_packages(),
       function(p) length(intersect(own, getNamespaceExports(p))),
       integer(1))
#>         basis7 distributions7 linkfunctions7    modelterms7    numericals7 
#>              0              0              0              0              0 
#>    optimizers7    parameters7     penalties7 
#>              0              0              0 
```
