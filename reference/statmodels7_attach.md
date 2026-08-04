# Attach the Member Packages

Attaches every installed member that is not on the search path already.

## Usage

``` r
statmodels7_attach()
```

## Value

The names of the packages attached, invisibly.

## Details

Each member's own startup messages are suppressed, since the caller
asked for the toolkit and gets one message about it rather than one per
package, and so are the masking warnings, which
[`statmodels7_conflicts`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
reports in one place instead.

## See also

[`statmodels7_packages`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
