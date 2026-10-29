# Drop the Combinations With Unordered Break-Points

Removes from the grid the rows in which the positions of one break-point
term are not strictly increasing in the order of its break-points. The
positions of a model are ordered and distinct, and a refit with two of
them held at equal or crossed values has no meaning.

## Usage

``` r
profile_prune(ctx, grid)
```

## Arguments

- ctx:

  The context of the call.

- grid:

  The data frame of held values.

## Value

`grid` without those rows, with its row names renumbered.
