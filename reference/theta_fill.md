# A Value Standing in for a Row That Cannot Be Read

A Value Standing in for a Row That Cannot Be Read

## Usage

``` r
theta_fill(v)
```

## Arguments

- v:

  A vector of parameter values.

## Value

The first finite value of `v`, or 1 where no value of `v` is finite, so
that a row with no covariance is carried through the vectorized
arithmetic and blanked afterwards.
