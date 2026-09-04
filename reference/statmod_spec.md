# Build a Model Specification

Builds everything
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
fits from: splits the formula into one equation per distribution
parameter, interprets each with modelterms7, and builds the terms it
names against the data.
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
calls this first; calling it directly is how a model is inspected
without being fitted.

## Usage

``` r
statmod_spec(
  formula,
  distrib,
  data,
  weights = NULL,
  offsets = NULL,
  need_response = TRUE,
  linpar = list()
)
```

## Arguments

- formula:

  The model formula, with the parameters' equations separated by `|`.

- distrib:

  A distributions7 distribution object, which decides how many equations
  there are and what they are called.

- data:

  A data frame holding the response, the covariates and any matrix
  columns the terms name.

- weights:

  Optional prior weights, a numeric vector of length `nrow(data)`, or
  `NULL` for all ones.

- offsets:

  Optional named list of offsets, one entry per parameter, summed with
  any the formula names.

- need_response:

  `TRUE`, the default, requires the left-hand side to evaluate. `FALSE`
  is what prediction uses: new data routinely has no response column.

- linpar:

  How the implicit parametric block is built, as
  [`linpar_options()`](https://statmodels7.github.io/statmodels7/reference/linpar_options.md)
  returns it. Kept on the specification, so a rebuild such as a fold of
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  reproduces the storage instead of quietly densifying.

## Value

A
[`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
object.

## The formula

The equations are separated by `|`, the first carrying the response, and
a parameter with no equation gets an intercept.
[`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
does the split, which is not the obvious one: R's precedence makes the
whole right-hand side of a three-equation formula the last term alone.

Each equation is interpreted in an environment where modelterms7's term
constructors sit in front of the search path, so `s()` means this
toolkit's even with mgcv attached.

## Prior weights are not normalized

They enter the log-likelihood as \\\sum_i w_i \ell_i\\ and are taken as
given. Dividing by their sum would turn the log-likelihood into a mean,
shrinking every standard error by \\\sqrt{n}\\ and making the
information criteria incomparable with an unweighted fit of the same
model.

## See also

[`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
for the split,
[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
for the assembly,
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
to fit it.

## Examples

``` r
dd <- data.frame(y = rnorm(20), x = runif(20), z = runif(20))
spec <- statmod_spec(y ~ x | sigma ~ z,
                     distributions7::gaussian1_distrib(), dd)

# One entry per parameter of the family, in the family's order.
names(spec@terms)
#> [1] "mu"    "sigma"
spec@equations
#> $mu
#> ~x
#> <environment: 0x561f064afec0>
#> 
#> $sigma
#> ~z
#> <environment: 0x561f064afec0>
#> 

# Unweighted, so the weights are ones.
c(n = spec@n_obs, total_weight = sum(spec@weights))
#>            n total_weight 
#>           20           20 
```
