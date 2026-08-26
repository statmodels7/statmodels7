# Attach the Member Packages

Attaches every installed member that is not on the search path already,
with [`library()`](https://rdrr.io/r/base/library.html). A member
already attached is left alone, so calling this twice attaches nothing
the second time, and a member that is not installed is skipped without
an error.

## Usage

``` r
statmodels7_attach()
```

## Value

The names of the packages this call attached, a character vector,
invisibly. Empty when every installed member was already on the search
path.

## Details

Each member's own startup messages are suppressed: the caller asked for
the toolkit and gets one message about it, not eight. The masking
warnings are suppressed too, and
[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
reports the same information once, at the end, in one block.

Attachment order follows
[`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md),
which is alphabetical. The members export no name in common, so nothing
about the result depends on that order;
[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
is what would report it if one day they did.

## See also

[`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
for the members,
[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
for what the suppressed warnings would have said.
