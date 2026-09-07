# Read a Two-Column Response as Successes and Trials

Splits a `cbind(successes, failures)` response into the successes and
the number of trials, returning the response and the distribution the
trials have been written onto.

## Usage

``` r
split_binomial_matrix(response, distrib)
```

## Arguments

- response:

  The evaluated left-hand side.

- distrib:

  The family.

## Value

A list with `response`, the successes, and `distrib`, carrying the row
sums as its `size`. Both unchanged where the response is not a
two-column matrix on a univariate family.

## Details

`cbind(y, n - y) ~ x` is how aggregated binomial data are written in R,
and it was not read here: the matrix went down as a vector of twice the
length and the run died on *"Parameter dimension mismatch. All
parameters should have length 1 or 400"*, a message about a length
nobody had asked for.

THE TRIALS COME FROM THE RESPONSE AND NOT FROM THE FAMILY, and that is
the half of this worth having. A size handed to the family is a vector
fixed at fitting time, so it cannot follow the rows anywhere else; the
row sums of a response expression are recomputed wherever that
expression is evaluated, which is what makes a fold of
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) or a
prediction at new data come out right. That the alternative is not
merely awkward but silently wrong is measured on the density itself: a
size of length 200 against 20 observations returns 200 log-densities,
summing to -1125.63, and signals nothing. See
[`check_trials()`](https://statmodels7.github.io/statmodels7/reference/check_trials.md),
which refuses that case where the row sums cannot be had.

A matrix response is left alone for a multivariate family, where it is
the ordinary thing and each row is one observation.

## See also

[`coerce_response()`](https://statmodels7.github.io/statmodels7/reference/coerce_response.md),
[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md),
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
