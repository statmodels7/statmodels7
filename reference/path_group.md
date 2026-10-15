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
so everything that reads `parameter`/`term`/`name` – the block it
belongs to, whether it scales the kink, the history a reader reads –
goes on reading the same fields, and the `members` list column says
which hyperparameters the axis stands for.

Where nothing is shared, each line is its own axis and carries itself as
its single member.
