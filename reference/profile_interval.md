# The Interval of a Position from Its Profile

Returns the smallest and largest grid positions of the run, connected to
`psi`, over which the increase of the criterion stays at or under
`cutoff`.

## Usage

``` r
profile_interval(x, delta, psi, cutoff)
```

## Arguments

- x:

  The grid of positions, increasing.

- delta:

  The increase of the criterion at each position.

- psi:

  The estimated position.

- cutoff:

  The cutoff of the increase.

## Value

A numeric vector of length two; `NA` where the profile is empty.
