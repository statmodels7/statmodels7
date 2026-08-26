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

  The other hyperparameters.

- name:

  Which one to solve for.

- target:

  The size the kink should have.

## Value

A single value, or `NA` where the target is out of reach.

## Details

The size of the kink is monotone in such a hyperparameter but not
necessarily increasing: a Laplace prior written by its scale has a kink
of \\1/\sigma\\, which narrows as the hyperparameter grows. Which way to
walk is therefore measured before walking, by comparing the size at the
current value with the size at twice it, and only then is the root
bracketed by doubling and found with
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html). A version
that assumed the size increases returned `NA` for every target on the
Laplace, having walked away from the answer.

Where the size is a power of the hyperparameter, which every kinked
penalty here turns out to be, the hyperparameter entering a separable
penalty as a scale is, the inversion is closed and the search is not run
at all.
[`kink_power()`](https://statmodels7.github.io/statmodels7/reference/kink_power.md)
measures the exponent and the answer is checked against the size before
it is returned, so a penalty that does not obey a power law falls back
to the search and is never assumed into one.

## See also

[`kink_power()`](https://statmodels7.github.io/statmodels7/reference/kink_power.md),
[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
