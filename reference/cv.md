# Choose the Hyperparameters by Cross-Validation

Builds the object with which
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
chooses the hyperparameters of the penalties with a kink by
cross-validation: the points of a path of values are scored by the
log-likelihood that the fit assigns to observations on which it was not
fitted, and the best is kept. The result is passed as
`sparse_criterion`. It is not an outer criterion:
`outer_criterion = cv()` is an error, because the smoothing parameters
are estimated by
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md).

## Usage

``` r
cv(nfolds = 10, folds = NULL, rule = c("min", "1se"))
```

## Arguments

- nfolds:

  How many folds, a single number of at least 2. The default is `10`. It
  is ignored when `folds` is given.

- folds:

  A fold number per observation, an integer vector as long as the data,
  for a user-defined partition or to compare two criteria on the same
  one. `NULL`, the default, draws a partition with a fixed seed.

- rule:

  `"min"` for the best criterion, `"1se"` for the sparsest fit within
  one standard error of it. It is matched with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

## Value

An
[`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
object of kind `"cv"`, carrying `nfolds`, `rule` and `folds`. Its
`hessian` is `"observed"` and its `k` is `NA_real_`; this criterion
reads neither.

## Cross-validation for a penalty with a kink

A marginal criterion,
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
approximates an integral by a Laplace expansion at the penalized mode,
which requires the second derivative of the penalty there. The mode of a
lasso, a SCAD or an MCP sits at the kink for every coefficient that it
sets to zero, and the derivative does not exist at the kink, so the
criterion is undefined at the point at which it would be read.

Cross-validation scores prediction and not a posterior, and it needs
nothing from the penalty beyond a fit.

Passed as `sparse_criterion` to a model with no hyperparameter of a
kinked penalty left to estimate, `cv()` has no effect.

## The criterion

The criterion is the mean over folds of \\-2\ell_f/w_f\\, where
\\\ell_f\\ is the weighted log-likelihood of the fold left out at the
coefficients fitted on the other folds, and \\w_f\\ is the sum of the
weights of its rows. Each training fit reapplies to its rows the terms
built on all the rows (knots, factor levels and scales of
standardization), as
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
does for new data, and carries the weights and offsets of the fit.

`rule = "1se"` takes the largest kink whose criterion is within one
standard error of the smallest, which is the sparsest fit that is not
measurably worse.

## The folds

Unless `folds` is given, the observations are assigned to folds by a
permutation drawn with a fixed seed, and the caller's random-number
state is restored afterwards, so that the partition does not depend on
the caller's seed.

## The path

The penalized mode is a piecewise smooth function of the hyperparameter:
differentiable while the active set holds, turning a corner whenever a
coefficient joins it or leaves. A criterion read there inherits the
corners, so a gradient search would read a slope about to change.

The hyperparameter is therefore swept over a grid, geometric in the size
of the kink, from the value that leaves the block empty down to
`min_ratio` of it. Every fit starts from the previous one's
coefficients.

## Which hyperparameters

The path covers every hyperparameter of the term's penalty that the
constructor did not hold at a number.

## A product within a term, an alternation between terms

For a term with several kinked hyperparameters, every combination is
visited under `search = "grid"` and one hyperparameter is swept at a
time under `"cyclic"`.

The choice belongs to the term, through `enet(X, search =)` and
[`modelterms7::term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.html),
and not to this criterion. The criterion is also applied to the model's
smooth hyperparameters, which are read at the mode and not swept, so
that an argument of this kind would have no meaning for them.

Between two terms the sweep alternates whichever search each one names,
so that `y ~ lasso(X) + enet(R)` costs the two blocks added and not
multiplied, and a term that sets a product does not make the other term
pay for it.

A term that names neither is swept over the product, because the cyclic
sweep traverses a cross through the point in hand and can stop where
each coordinate is separately best without being jointly so. Its cost is
the product of the term's grids, whereas the cost of the cyclic sweep is
their sum. At two hyperparameters this is `n_lambda * n_alpha` fits
against `n_lambda + n_alpha` per pass. If three or more hyperparameters
are estimated, the cost of the product grows exponentially, and
`"cyclic"` is there for that case.

## How long the path is, and how far down it reaches

The length of the path and its depth are settings of the term, written
on its own signature: `lasso(x, n_lambda = 25, min_ratio = 1e-4)` and
`enet(x, n_lambda = 25, n_alpha = 5)`.

The two numbers differ because the axes do. \\\lambda\\ descends the
size of the kink over four decades, while \\\alpha\\ spans one bounded
interval, so that the default product is 25 by 5.

## The shape of the grid

The grid is not a rectangle for the elastic net and is a rectangle for
SCAD and MCP. The elastic net's kink is \\\lambda\alpha\\, so the value
emptying the block is \\\lambda\_{\max} = \kappa/\alpha\\ and every
\\\alpha\\ carries its own \\\lambda\\ axis, descending from its own
top. The shapes of SCAD and MCP do not change the kink, so that
\\\lambda\_{\max}\\ is one number whatever the shape and the two axes
form a rectangle. Each axis is built at the settings of the axes outside
it, which covers both cases without a rule about families.

## The cost

Each point of the path costs `nfolds` fits. The number of points is the
term's `n_lambda` for one kinked hyperparameter, and the product of its
axes for several. The fits along a path are warm-started within each
fold.

The folds run over separate R processes when
`statmod(threads = n_threads(workers = 4))` requests them. One seed is
drawn per fold in the parent, so that the answer is identical at any
worker count.

## References

Breiman, L., Friedman, J. H., Olshen, R. A. and Stone, C. J. (1984).
*Classification and Regression Trees*. Wadsworth.

Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization paths
for generalized linear models via coordinate descent. *Journal of
Statistical Software* 33(1), 1–22.

## See also

[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
for the other prediction-error criteria,
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
for the marginal ones,
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for where this is passed,
[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
for reading back what it chose.

## Examples

``` r
# Two of eight columns carry signal; the rest are noise for the path to
# shrink away.
set.seed(1)
X <- matrix(rnorm(120 * 8), 120, 8)
dd <- data.frame(y = 1 + X %*% c(2, -1.5, rep(0, 6)) + rnorm(120, sd = 0.5))
dd$x <- X

fit <- statmod(y ~ lasso(x, n_lambda = 12),
               distributions7::gaussian1_distrib(), dd,
               sparse_criterion = cv(nfolds = 5))
hyper(fit)
#>   parameter                    term   name estimate  held source   id
#> 1        mu lasso(x, n_lambda = 12) lambda 14.03983 FALSE     cv <NA>

# The two columns that carry signal are recovered.
round(coef(fit)$mu[2:3], 3)
#> lasso.1 lasso.2 
#>   1.995  -1.512 

# "1se" takes the largest kink within one standard error of the best, so
# it never keeps more coefficients than "min". On this data the two agree,
# because the only sparser point of the path has a criterion more than one
# standard error above the minimum.
one_se <- statmod(y ~ lasso(x, n_lambda = 12),
                  distributions7::gaussian1_distrib(), dd,
                  sparse_criterion = cv(nfolds = 5, rule = "1se"))
c(min = hyper(fit)$estimate, one_se = hyper(one_se)$estimate)
#>      min   one_se 
#> 14.03983 14.03983 
c(min = sum(coef(fit)$mu[-1] != 0), one_se = sum(coef(one_se)$mu[-1] != 0))
#>    min one_se 
#>      7      7 
```
