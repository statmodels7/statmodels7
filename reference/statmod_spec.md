# Build a Model Specification

Builds everything from which
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
fits: splits the formula into one equation per distribution parameter,
interprets each with modelterms7, and builds the terms that it names
against the data.
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
calls this first, and calling it directly inspects a model without
fitting it.

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

  A distributions7 distribution object, which determines how many
  equations there are and what they are called.

- data:

  A data frame holding the response, the covariates and any matrix
  columns that the terms name. A missing value in the response, or in a
  variable that an equation names, is an error, because a fit does not
  drop rows.

- weights:

  Optional prior weights, a numeric vector of length `nrow(data)` (not
  the name of a column), or `NULL` for all ones. The weights must be
  non-negative, finite and complete, and not all zero. They act as
  frequency weights and are not normalized.

- offsets:

  Optional named list of offsets, one entry per parameter, each a finite
  numeric vector of length 1 or `nrow(data)`, summed with any offset the
  formula names for the same parameter.

- need_response:

  `TRUE`, the default, requires the left-hand side to evaluate.
  Prediction uses `FALSE`, because new data routinely has no response
  column.

- linpar:

  How the implicit parametric block is built, as
  [`linpar_options()`](https://statmodels7.github.io/statmodels7/reference/linpar_options.md)
  returns it. The options are kept on the specification, so that a
  rebuild such as a fold of
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  reproduces the storage instead of quietly densifying.

## Value

A
[`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
object.

## The formula

The equations are separated by `|`, and the first of them carries the
response. A parameter with no equation is given an intercept.
[`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
does the split, which is not the obvious one: R's precedence reduces the
whole right-hand side of a three-equation formula to the last term.

Each equation is interpreted in an environment where modelterms7's term
constructors sit in front of the search path, so `s()` means this
toolkit's even with mgcv attached.

## Rejected before any fit

The function signals an error for a multivariate family (which
[`distributions7::fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.html)
fits with constant parameters), a censored response written with
`cens()`, a family method whose signature has no `...`, a missing value
in the response or in a variable that an equation names, a two-column
response on a family without `size`, an offset written inside another
term, and a structural term that cannot be fitted or a formula with more
than one structural term.

## Prior weights are not normalized

They enter the log-likelihood as \\\sum_i w_i \ell_i\\ and are taken as
given, so they act as frequency weights: a constant weight \\c\\ on
every row is \\c\\ copies of the data. Dividing the weights by their sum
would turn the log-likelihood into a mean, multiplying every standard
error by \\\sqrt{n}\\ and making the information criteria incomparable
with an unweighted fit of the same model.

The restricted (REML) criterion counts \\c n\\ observations against the
\\p\\ coefficients when the weights are the constant \\c\\. The REML
estimate of a dispersion parameter therefore changes under a constant
rescaling of the weights, and the ML estimate does not.

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

# The terms have one entry per parameter of the family, in the family's
# order.
names(spec@terms)
#> [1] "mu"    "sigma"
spec@equations
#> $mu
#> ~x
#> <environment: 0x5588b13f3bc8>
#> 
#> $sigma
#> ~z
#> <environment: 0x5588b13f3bc8>
#> 

# Unweighted, so the weights are ones.
c(n = spec@n_obs, total_weight = sum(spec@weights))
#>            n total_weight 
#>           20           20 
```
