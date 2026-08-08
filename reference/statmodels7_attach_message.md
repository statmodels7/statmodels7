# The Body of the Attach Message

The member packages laid out in two columns with their versions, padded
so that the second column and the versions within it start at a fixed
offset.

## Usage

``` r
statmodels7_attach_message(versions)
```

## Arguments

- versions:

  A data frame as
  [`statmodels7_versions`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
  returns.

## Value

A character vector, one element per row of the layout, empty when there
are no members.

## Details

Padding is applied inside a line and removed at its end, since a
trailing run of spaces is invisible where it is written and visible
wherever the message is pasted. A member that is not installed is marked
rather than shown blank: the two are different states and only one of
them is a complete toolkit. Plain ASCII throughout: a non-ASCII
character in a package's own output is a portability risk taken for a
decoration.

## See also

[`statmodels7_versions`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
