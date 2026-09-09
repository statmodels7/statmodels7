# Start From Another Fit's Estimates

Takes the coefficients a model already fitted found, matches them to the
model about to be fitted, and starts there. What the two models do not
share is left to a second strategy.

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
  by default, which is what a fit given no strategy at all uses.

## Value

A `StartFrom` object, inheriting from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).

## Details

Two models of the same response often share most of their coefficients:
one drops a covariate, one adds a term, one changes the family. The
estimates of the first are then a far better starting point for the
second than any rule that looks only at the response, and there is no
reason to find them twice.

## What is matched, and what is projected

The two designs are compared BLOCK BY BLOCK, a block being one term's
columns, and three things can happen to a block:

- the parametric block:

  is matched COLUMN BY COLUMN, by name. Adding or dropping a covariate
  is the case this exists for, and the names of a model matrix are
  variables, levels and interactions, which mean what they say.

- a block written on the same basis:

  is carried across as it stands, its coefficient names being the same
  coordinates. This is exact and costs nothing.

- a block written on ANOTHER basis:

  is estimated by least squares against what the reference's predictor
  leaves once the blocks already carried across are removed.

The third rule is what makes two bases comparable at all. A smooth's
coefficients are coordinates in a basis the fit rotates, so
`s(x, bspline_smooth(k = 6))` and `s(x, bspline_smooth(k = 10))` carry
the names `s(x).z1` to `s(x).z4` in common and mean something different
by each: measured on one data set, `s(x).z1` is 0.0267 at `k = 6` and
-0.0277 at `k = 10`, opposite in sign. What the two do share is the
FUNCTION, so the coefficients are found by projecting it,
\$\$\hat\beta_S = \arg\min\_\beta \lVert X_S\beta - r\rVert^2, \qquad r
= X^{\mathrm{ref}}\beta^{\mathrm{ref}} - X\_{-S}\beta\_{-S},\$\$ over
the pending columns \\S\\ alone. The equality
\\B\_{\mathrm{new}}\beta\_{\mathrm{new}} =
B\_{\mathrm{old}}\beta\_{\mathrm{old}}\\ has no solution unless the old
span sits inside the new one, and the projection is what remains: the
closest the new basis can come to the function it is started from. Where
the spans do nest it is exact, and where the two bases coincide it
returns the coefficients it was given, so the second rule above is a
fast path rather than a different answer.

Two blocks are paired by the stem of their coefficient names and by the
term's class, not by the block's own key, which is the term's deparsed
call and therefore differs the moment `k` does.

A distribution parameter the reference does not have – fitting a
negative binomial from a Poisson, where `theta` is new – is left to
`rest`, so the families need not agree. So is any block the reference
has no counterpart for, and any block at all when the two responses
differ, when either model carries a structural term, or when the
equation carries a block that moves with its coefficients.

## What it buys

On the blocks that MATCH, between 1.1 and 2.1 times, and the honest
reading is that it is not more because nothing was starting from zero:
[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
the default, already fits the intercept-only model. On a Poisson at \\n
= 4000\\ refitted with one covariate dropped it is 1.12 to 1.57 times,
on a two-equation gaussian 1.89, and on a negative binomial with 40
columns at \\n = 20000\\ it is 2.05 – 0.81 s against 0.40 s. Fitting a
negative binomial from a Poisson's estimates is 1.29.

On a block that is PROJECTED the starting predictor is a different order
of magnitude closer. Measured at \\n = 4000\\, widening
`s(x, bspline_smooth(k = 6))` to `s(x, bspline_smooth(k = 12))`, the
root mean square gap between the starting predictor and the reference's
is \\3\times 10^{-15}\\ where leaving the block to the fallback gives
1.046; a quadratic basis carried onto a cubic one at the same `k` gives
0.0101 against 1.049. What that is worth in time depends on the family:
a gaussian on the identity link solves its inner problem in one step
whatever the start, so it is 1.09 to 1.26 times, while a Poisson goes
from 18 criterion evaluations to 6 and 1.94 times, and a Poisson that
also drops a covariate 1.71.

The gain is worth having where the same model is refitted many times
over: a coefficient held at a sequence of values, which is what an
interval by inversion walks, or a model built up one term at a time.

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
