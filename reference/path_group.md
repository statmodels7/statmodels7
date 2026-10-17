# Collapse the Shared Hyperparameters of a Path

Turns one line per estimated hyperparameter into one line per axis:
those carrying the same label become a single axis, swept once, whose
value is written into each of them.

## Usage

``` r
path_group(rows)
```

## Arguments

- rows:

  One line per estimated hyperparameter, with an `id` column.

## Value

One line per axis, with a `members` list column and no `id`.

## Details

It plays for a path the role that
[`index_group()`](https://statmodels7.github.io/statmodels7/reference/index_group.md)
plays for the outer search. The axis is identified by its first member,
so the code that reads `parameter`, `term` and `name` (to find the block
to which the axis belongs, to find whether it scales the kink, or to
read the history) continues to read the same fields, and the `members`
list column lists the hyperparameters that the axis represents.

If a line shares its label with no other line, it is its own axis and
carries itself as its single member.
