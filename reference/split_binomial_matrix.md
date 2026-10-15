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
sums as its `size`. Both are unchanged where the response is not a
matrix or the family is multivariate. A matrix response on a univariate
family signals an error where it has a number of columns other than two,
is not numeric, holds a negative count, or where the family has no
`size`.

## Details

`cbind(y, n - y) ~ x` is how aggregated binomial data are written in R.
The two columns are read as the successes and the failures: the response
becomes the successes and the family's `size` becomes the row sums.

The number of trials is taken from the response and not from the family.
A size given to the family is a vector fixed at fitting time, so it
cannot follow the rows anywhere else, whereas the row sums of a response
expression are recomputed wherever that expression is evaluated, which
is what makes a fold of
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) or a
prediction at new data use the right trials. A size of another length
than the rows is recycled by the density without any signal (a size of
length 200 against 20 observations returns 200 log-densities). See
[`check_trials()`](https://statmodels7.github.io/statmodels7/reference/check_trials.md),
which rejects that case where the row sums cannot be had.

A size given to the family as well must agree with the row sums, and the
two statements disagreeing is an error; the default size of one is not a
statement and is replaced.

A matrix response is left alone for a multivariate family, where it is
the ordinary response and each row is one observation.

## See also

[`coerce_response()`](https://statmodels7.github.io/statmodels7/reference/coerce_response.md),
[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md),
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
