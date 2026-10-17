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

The names of the packages that this call attached, a character vector,
invisibly. The vector is empty when every installed member was already
on the search path.

## Details

Each member's own startup messages are suppressed, so that attaching the
toolkit prints one message instead of one per member. The masking
warnings are suppressed too, and
[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
reports the same information once, at the end, in one block.

Attachment order follows
[`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md),
which is alphabetical. The members export no name in common, so the
result does not depend on that order. If a name were exported by two
members,
[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
would report it.

## See also

[`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
for the members,
[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
for the conflicts that the suppressed warnings would have reported.
