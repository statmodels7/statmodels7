# The Hyperparameter That Gives the Kink a Chosen Size

Solves `kink_scale(pen, theta) == target` in one named hyperparameter.

## Usage

``` r
kink_solve(pen, theta, name, target)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters in force; the value of `name` is where the search
  starts.

- name:

  Which one to solve for.

- target:

  The size that the kink should have.

## Value

A single value, or `NA` where the target is out of reach.

## Details

The size of the kink is monotone in such a hyperparameter but not
necessarily increasing: a Laplace prior written by its scale has a kink
of \\1/\sigma\\, which narrows as the hyperparameter grows. The
direction of the walk is therefore measured before the root is
bracketed, by comparing the size at the current value with the size at
twice it. The root is then bracketed by doubling or halving and found
with [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html).

Where the size is a power of the hyperparameter, as it is for the
shipped kinked penalties because the hyperparameter enters a separable
penalty as a scale, the inversion is closed and the search is not run.
[`kink_power()`](https://statmodels7.github.io/statmodels7/reference/kink_power.md)
measures the exponent and
[`kink_by_power()`](https://statmodels7.github.io/statmodels7/reference/kink_by_power.md)
checks the answer against the size before it is returned, so that a
penalty that does not obey a power law falls back to the search.

## See also

[`kink_power()`](https://statmodels7.github.io/statmodels7/reference/kink_power.md),
[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
