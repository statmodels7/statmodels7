# The Body of the Attach Message

Lays the member packages out in two columns with their versions, padded
so that the second column and the versions within it start at a fixed
offset. Each entry is prefixed `v` when the member is installed and `x`
when it is not, and the two columns are filled down the left one first.

## Usage

``` r
statmodels7_attach_message(versions)
```

## Arguments

- versions:

  A data frame as
  [`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
  returns, with character columns `package` and `version` and `NA` in
  `version` for a member that is not installed. The argument is not
  validated, because `.onAttach()` calls this function with the result
  of that function.

## Value

A character vector with `ceiling(nrow(versions) / 2)` elements, one per
row of the two-column layout, in which trailing spaces are already
trimmed from each element. The result is `character(0)` when `versions`
has no rows.

## Details

Padding is applied inside a line and stripped at its end, because a
trailing run of spaces is invisible where it is written and becomes
visible when the message is pasted elsewhere. A member that is not
installed is marked `x` and shown with the text `not installed`, so that
an absent member differs visibly from a present one.

The message uses plain ASCII only. A non-ASCII character in a startup
message can print incorrectly in a session with another encoding, and
the layout needs none.

## See also

[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
for the input,
[`format_conflicts()`](https://statmodels7.github.io/statmodels7/reference/format_conflicts.md)
for the other half of the same message.
