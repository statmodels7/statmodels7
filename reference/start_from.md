# Start From Another Fit's Estimates

Takes the coefficients that were estimated in a fitted model, matches
them to the model that is about to be fitted, and starts the fit there.
The coefficients that the two models do not share are set by a second
strategy.

## Usage

``` r
start_from(fit, rest = start_intercepts())
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
  to take the estimates from.

- rest:

  The strategy for whatever the two models do not share.
  [`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md)
  by default, which is what a fit given no strategy uses. It receives
  the objective, so
  [`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md)
  can serve as `rest`.

## Value

A `StartFrom` object, inheriting from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).

## Details

Two models of the same response often share most of their coefficients:
one drops a covariate, one adds a term, one changes the family. The
estimates of the first model are then a closer starting point for the
second than a rule that uses only the response.

## What is matched, and what is projected

The two designs are compared block by block, a block being the columns
of one term, and a block is treated in one of three ways:

- the parametric block:

  is matched column by column, by name. The names of a model matrix are
  variables, levels and interactions, so a covariate that is added or
  dropped leaves the other columns matched.

- a block written on the same basis:

  is carried across as it stands, because its coefficient names are the
  same coordinates.

- a block written on another basis:

  is estimated by least squares against the part of the reference's
  predictor that the blocks already carried across do not explain.

The third rule makes two bases comparable. The coefficients of a smooth
are coordinates in a basis that the fit rotates, so
`s(x, bspline_smooth(k = 6))` and `s(x, bspline_smooth(k = 10))` carry
the names `s(x).z1` to `s(x).z4` in common but give them different
meanings, and one name does not identify one function. What the two
bases share is the fitted function, so the coefficients are found by
projecting it, \$\$\hat\beta_S = \arg\min\_\beta \lVert X_S\beta -
r\rVert^2, \qquad r = \eta^{\mathrm{ref}} - \eta\_{-S},\$\$ over the
pending columns \\S\\ alone, where \\\eta^{\mathrm{ref}}\\ is the
reference's predictor and \\\eta\_{-S}\\ is the new model's predictor
with the pending columns at zero. Both predictors are the actual
predictors of their models, including the value of a block that moves
with its coefficients (a break-point or a nonlinear term). The equality
of the two blocks, \\B\_{\mathrm{new}}\beta\_{\mathrm{new}} =
B\_{\mathrm{old}}\beta\_{\mathrm{old}}\\, has a solution only if the old
span lies inside the new one, and the projection gives the closest
function that the new basis can represent. Where the spans are nested
the projection is exact, and where the two bases coincide it returns the
coefficients it was given, so the second rule above is a shortcut and
not a different answer.

Two blocks are paired by the stem of their coefficient names and by the
class of the term, and not by the key of the block, which is the
deparsed call of the term and therefore differs as soon as `k` does.

A distribution parameter that the reference does not have, such as
`theta` when a negative binomial is fitted from a Poisson, is set by
`rest`, so the families need not agree. So is any block for which the
reference has no counterpart. A block that matches by name is carried
across whatever the response is. The projection of a block written on
another basis is skipped, and the block is set by `rest`, when the two
responses differ or when either model carries a structural term. A
pending block that itself moves with its coefficients (for example
`seg(x)` refitted as `seg(x, npsi = 2)`) is not projected either,
because its coefficients are not those of a linear block, and it keeps
the values given by `rest`.

The starting values are most useful where the same model is refitted
many times, as with a coefficient held at a sequence of values in an
interval by inversion, or with a model built up one term at a time.

## See also

[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
the default;
[`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md),
the generic.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(200), z = runif(200))
dd$y <- rpois(200, exp(0.4 + 0.6 * dd$x - 0.3 * dd$z))
full <- statmod(y ~ x + z, distributions7::poisson_distrib(), dd)

# the same model with one covariate dropped, started where the other
# one ended
sub <- statmod(y ~ z, distributions7::poisson_distrib(), dd,
               start = start_from(full))
coef(sub)
#> $mu
#> (Intercept)           z 
#>   0.7674242  -0.2566149 
#> 

# what was taken, and from where
spec <- statmod_spec(y ~ z, distributions7::poisson_distrib(), dd)
attr(start_at(start_from(full), spec, statmod_design(spec), NULL), "taken")
#>   parameter   term coefficient      value     how
#> 1        mu linpar (Intercept)  0.4752490 matched
#> 2        mu linpar           z -0.2913813 matched
```
