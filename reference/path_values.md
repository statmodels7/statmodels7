# The Values a Path Visits

A geometric grid of kink sizes from the one that empties the block down
to `min_ratio` of it, carried back onto the hyperparameter.

## Usage

``` r
path_values(pen, theta, name, s_max, n_values = 40L, min_ratio = 1e-04)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters in force.

- name:

  Which one the path varies.

- s_max:

  The size of the kink at the top of the path.

- n_values:

  How many points.

- min_ratio:

  The smallest kink size, as a fraction of `s_max`.

## Value

A numeric vector of values for `name`, from the emptiest fit to the
fullest.

## Details

The grid is geometric in the size of the kink rather than in the
hyperparameter, so that a penalty whose kink narrows as its
hyperparameter grows is swept in the same order as one whose kink
widens: from the empty model towards the full one. Values the penalty
cannot reach are dropped.

The exponent relating the two is read ONCE and every target inverted
through it, rather than each target bracketed on its own. Measured, a
bracketing solve costs 4.18 ms against a fit's 62.5 ms, so a path of
twenty-five points spent 6.7 per cent of itself locating the values it
would visit; through the exponent the whole grid costs four evaluations
of the size.
[`kink_by_power`](https://statmodels7.github.io/statmodels7/reference/kink_by_power.md)
checks the relation before the values are used and returns `NA` where it
does not hold, and those fall back to
[`kink_solve`](https://statmodels7.github.io/statmodels7/reference/kink_solve.md)
one at a time.

## See also

[`kink_power`](https://statmodels7.github.io/statmodels7/reference/kink_power.md),
[`path_forced`](https://statmodels7.github.io/statmodels7/reference/path_forced.md)
