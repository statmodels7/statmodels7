# The Chains of Points That a Start Visits

Splits the rows of the grid into the sequences along which a refit
starts from the previous one (`"neighbor"`) or into single points
(`"fit"`).

## Usage

``` r
profile_chains(ctx, grid, start)
```

## Arguments

- ctx:

  The context of the call.

- grid:

  The data frame of held values.

- start:

  One of `"neighbor"`, `"fit"` and `"both"`.

## Value

A list of chains, each a list with the row indices `rows` and whether
the start `carry`ies from one to the next.
