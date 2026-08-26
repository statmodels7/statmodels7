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
  `version` for a member that is not installed. Not validated: this is
  called from `.onAttach()` with that function's own result.

## Value

A character vector with `ceiling(nrow(versions) / 2)` elements, one per
row of the two-column layout, each already trimmed of trailing spaces.
`character(0)` when `versions` has no rows.

## Details

Padding is applied inside a line and stripped at its end, since a
trailing run of spaces is invisible where it is written and visible
wherever the message is pasted. A member that is not installed is marked
`x` and named with the text `not installed`; leaving it blank would hide
the difference between an absent member and a present one.

Plain ASCII throughout. A non-ASCII character in a package's startup
message is a portability risk taken for a decoration.

## See also

[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
for the input,
[`format_conflicts()`](https://statmodels7.github.io/statmodels7/reference/format_conflicts.md)
for the other half of the same message.
