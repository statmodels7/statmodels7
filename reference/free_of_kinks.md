# A Score Over the Coordinates That a Kink Leaves Free

Zeroes the entries of a score at the coordinates that
[`zero_kinked()`](https://statmodels7.github.io/statmodels7/reference/zero_kinked.md)
returns.

## Usage

``` r
free_of_kinks(score, spec, design, coef)
```

## Arguments

- score:

  A score over the stacked coefficients.

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

- coef:

  The coefficients, a named list.

## Value

`score`, with those entries set to zero.

## Details

At a coordinate that a kinked penalty holds at zero, the objective has
no derivative, and the entry that the smooth part reports is the
log-likelihood's score, which lies in the subgradient interval of the
kink. Read as a residual, this entry would inflate the mode error of a
fit sitting at its mode, so it is set to zero. Whether such a coordinate
should leave zero is the condition of the kinked block, which
[`alternation_readings()`](https://statmodels7.github.io/statmodels7/reference/alternation_readings.md)
reads.
