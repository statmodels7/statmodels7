# Check the `marginal` Argument of a Marginal Criterion

Checks the shape of `marginal` before any family is known: `NULL`, or a
character vector with no missing and no empty entry. Whether the names
belong to the family is checked by
[`marginal_params()`](https://statmodels7.github.io/statmodels7/reference/marginal_params.md),
at fit time.

## Usage

``` r
check_marginal_arg(marginal)
```

## Arguments

- marginal:

  The argument as given.

## Value

`marginal`, unchanged, or an error naming the argument.
