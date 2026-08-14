# Invert the Size of the Kink Through a Power Law

The value giving each target size, checked against the size it produces.

## Usage

``` r
kink_by_power(pen, theta, name, target, pw)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters in force.

- name:

  Which one to solve for.

- target:

  The sizes the kink should have.

- pw:

  What
  [`kink_power`](https://statmodels7.github.io/statmodels7/reference/kink_power.md)
  returned, or `NULL`.

## Value

A numeric vector as long as `target`, `NA` where the power law did not
answer.

## Details

The check is what licenses the closed route: the exponent came from two
points and the relation is asserted at a third before the answer is
used. A target that misses by more than a rounding, or that lands
outside the hyperparameter's own interval, comes back `NA` and the
caller falls back to bracketing.

## See also

[`kink_power`](https://statmodels7.github.io/statmodels7/reference/kink_power.md),
[`kink_solve`](https://statmodels7.github.io/statmodels7/reference/kink_solve.md)
