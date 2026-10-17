# How Far Above Its Mode a Restricted Fit Stopped

Computes \\\tfrac12 g'K^{-1}g\\ over the free coordinates of a
restricted fit, which is the decrease of the penalized objective that a
Newton step from the point returned by
[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md)
would predict, in log-likelihood units.

## Usage

``` r
restricted_mode_error(score, information, at)
```

## Arguments

- score:

  The restricted fit's score, all coordinates.

- information:

  Its penalized information, all coordinates.

- at:

  The position of the held coordinate.

## Value

A single number, or `NA` where the score on the free coordinates is not
finite, where there is no free coordinate, or where the information on
the free coordinates cannot be inverted.

## Details

[`inner_mode_error()`](https://statmodels7.github.io/statmodels7/reference/inner_mode_error.md)
computes the same quantity inside a search; here it is applied to a
restricted fit. The inner optimizer's flag records whether a stopping
rule fired, which is a boolean about a threshold on a score whose size
depends on the model, whereas whether a point is at its mode is a matter
of distance and has a natural scale. A stopping rule reads a change in
the objective and this quantity reads the point, so the two need not
agree, and the print method of
[`StatmodTest()`](https://statmodels7.github.io/statmodels7/reference/StatmodTest-class.md)
reads the mode error and not the flag.

The held coordinate is excluded because it is not free: its score is the
quantity that Rao's statistic reads, and it does not vanish under the
restriction.

## See also

[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
the limit against which it is read, and
[`inner_mode_error()`](https://statmodels7.github.io/statmodels7/reference/inner_mode_error.md),
the same reading inside a search.
