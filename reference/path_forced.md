# The Values a Caller Wrote Out

The grid the term named, ordered from the emptiest fit towards the
fullest so that the warm starts run the way every other path here runs.

## Usage

``` r
path_forced(pen, theta, name, values)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters in force.

- name:

  Which one the path varies.

- values:

  What the term wrote out.

## Value

A numeric vector, from the emptiest fit to the fullest.

## Details

Which end is the sparse one is a property of the penalty and not of the
numbers: the kink of a lasso widens with \\\lambda\\ and that of a
Laplace prior written by its scale narrows with \\\sigma\\, so the order
is settled by asking the penalty which way its kink moves rather than by
sorting downwards. Nothing else is applied – the value that empties the
block does not cap the grid and `min_ratio` does not extend it, both of
those being ways to build one.

## See also

[`term_values`](https://statmodels7.github.io/modelterms7/reference/term_values.html),
[`path_values`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
