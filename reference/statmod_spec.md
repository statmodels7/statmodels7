# Build a Model Specification

Splits the formula into one equation per distribution parameter,
interprets each with modelterms7 and builds its terms against the data.

## Usage

``` r
statmod_spec(
  formula,
  distrib,
  data,
  weights = NULL,
  offsets = NULL,
  need_response = TRUE
)
```

## Arguments

- formula:

  The model formula.

- distrib:

  A distributions7 distribution object.

- data:

  A data frame.

- weights:

  Optional prior weights, one per observation.

- offsets:

  Optional named list of offsets, one per parameter.

- need_response:

  Whether the left-hand side must evaluate. A likelihood needs it; a
  prediction does not, and new data routinely has no response column.

## Value

An object of class
[`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Details

The equations are interpreted in an environment where modelterms7's term
constructors shadow whatever the user has attached, so that `s()` means
ours even with mgcv on the search path. See
[`statmod_equations`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
for the split itself, which is not the obvious one.

Prior weights enter the log-likelihood as \\\sum_i w_i \ell_i\\ and are
taken as given. They are deliberately NOT normalized: dividing by their
sum would turn the log-likelihood into a mean, shrinking every standard
error by \\\sqrt{n}\\ and making the information criteria incomparable
with an unweighted fit of the same model.

## See also

[`statmod_equations`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md),
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
dd <- data.frame(y = rnorm(20), x = runif(20), z = runif(20))
spec <- statmod_spec(y ~ x | sigma ~ z, distributions7::gaussian1_distrib(), dd)
names(spec@terms)
#> [1] "mu"    "sigma"
```
