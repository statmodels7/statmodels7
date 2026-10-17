# Invert the Size of the Kink Through a Power Law

Returns, for each target size, the value that gives it, checked against
the size that the value produces.

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

  The sizes that the kink should have.

- pw:

  What
  [`kink_power()`](https://statmodels7.github.io/statmodels7/reference/kink_power.md)
  returned, or `NULL`.

## Value

A numeric vector as long as `target`, `NA` where the power law gives no
value.

## Details

The exponent comes from two points, and the relation is checked at the
first and the last admissible target before any answer is used. If
either check misses by more than a rounding, every target comes back
`NA`; a target whose value lands outside the hyperparameter's own
interval comes back `NA` on its own. The caller then falls back to
bracketing.

## See also

[`kink_power()`](https://statmodels7.github.io/statmodels7/reference/kink_power.md),
[`kink_solve()`](https://statmodels7.github.io/statmodels7/reference/kink_solve.md)
