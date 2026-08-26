# Choose the Hyperparameters by Cross-Validation

Builds the object that tells
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
to choose a model's hyperparameters by cross-validation: a path of
values is scored by the log-likelihood the fit assigns to observations
it was not fitted on, and the best is kept. Pass the result as
`sparse_criterion`, or as `outer_criterion`, or as both.

## Usage

``` r
cv(nfolds = 10, folds = NULL, rule = c("min", "1se"))
```

## Arguments

- nfolds:

  How many folds, a single number of at least 2. `10` by default.
  Ignored when `folds` is given.

- folds:

  A fold number per observation, an integer vector as long as the data,
  for a partition of your own or to compare two criteria on the same
  one. `integer(0)`, the default, draws them at fit time.

- rule:

  `"min"` for the best criterion, `"1se"` for the sparsest fit within
  one standard error of it. Matched with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

## Value

An
[`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
object of kind `"cv"`, carrying `nfolds`, `rule` and `folds`. Its
`hessian` is `"observed"` and its `k` is `NA_real_`, neither being read
by this criterion.

## Why a penalty with a kink needs this

A marginal criterion,
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
approximates an integral by a Laplace expansion at the penalized mode,
which asks for the second derivative of the penalty there. The mode of a
lasso, a SCAD or an MCP sits at the kink for every coefficient it sets
to zero, and that is where the derivative does not exist, so the
criterion is undefined at the point it would be read at.

Cross-validation asks a different question, about prediction and not
about a posterior, and asking it needs nothing from the penalty beyond a
fit.

## The criterion

The mean over folds of \\-2\ell/n_f\\ on the fold left out, each
training fit rebuilding the design on its own rows.

`rule = "1se"` takes the largest kink whose criterion is within one
standard error of the smallest, which is the sparsest fit that is not
measurably worse.

## A path, not a search

The penalized mode is a piecewise smooth function of the hyperparameter:
differentiable while the active set holds, turning a corner whenever a
coefficient joins it or leaves. A criterion read there inherits the
corners, so a gradient search would read a slope about to change.

The hyperparameter is therefore swept over a grid, geometric in the size
of the kink, from the value that leaves the block empty down to
`min_ratio` of it. Every fit starts from the previous one's
coefficients.

## Which hyperparameters

The term's answer: every one of its penalty's that the constructor did
not hold at a number. What the criterion decides is how they are
covered.

## A product within a term, an alternation between terms

A term carrying several kinked hyperparameters has every combination
visited under `search = "grid"` and one swept at a time under
`"cyclic"`.

The choice belongs to the term, through `enet(X, search =)` and
[`modelterms7::term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.html),
never to this criterion. The same criterion is put to the model's smooth
hyperparameters as well, and those are read at the mode instead of being
swept, so most of what it is asked about could not use such an argument.

Between two terms the sweep alternates whichever each one named, so
`y ~ lasso(X) + enet(R)` costs the two blocks added, never multiplied,
and one term asking for a product does not make the other pay for it.

A term that names neither gets the product, because the cyclic sweep
traverses a cross through the point in hand and can stop where each
coordinate is separately best without being jointly so. Its cost is the
product of the term's grids where the cyclic sweep's is their sum, which
at two hyperparameters is `n_lambda * n_alpha` fits against
`n_lambda + n_alpha` per pass; with three or more estimated it grows
exponentially, and `"cyclic"` is there for that.

## How long the path is, and how far down it reaches

The term's answers again, written on its own signature where a reader
can see them: `lasso(x, n_lambda = 25, min_ratio = 1e-4)`,
`enet(x, n_lambda = 25, n_alpha = 5)`.

The two numbers differ because the axes do. \\\lambda\\ descends the
size of the kink over four decades, while \\\alpha\\ spans one bounded
interval, so the shipped product is 25 by 5.

## The grid is not a rectangle

For two different reasons. The elastic net's kink is \\\lambda\alpha\\,
so the value emptying the block is \\\lambda\_{\max} = \kappa/\alpha\\
and every \\\alpha\\ carries its own \\\lambda\\ axis, descending from
its own top. The shapes of SCAD and MCP leave the kink alone, so there
\\\lambda\_{\max}\\ is one number whatever the shape and the two axes
really are a rectangle. Each axis is built at the settings of the axes
outside it, which covers both cases with no rule about families written
down anywhere.

## The cost

`nfolds` fits per point of the path, and how many points there are is
the term's `n_lambda` for one kinked hyperparameter, or the product of
its axes for several.

The warm starts are worth 1.8 times, and building each fold's design
once instead of once per point another 4 per cent. What remains is the
proximal iteration: measured at 200 observations and 20 columns, 0.88
seconds a fit, against `cv.glmnet`'s 0.03 seconds for its whole path of
100 values on five folds. That distance is why `n_lambda` is 25 here and
100 there.

The folds themselves run over separate R processes when
`statmod(threads = n_threads(workers = 4))` asks for it, with one seed
drawn per fold in the parent, so the answer is identical at any worker
count. Measured, a ten-fold lasso cross-validation goes from 31.45 s to
14.18 s at four workers.

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
               outer_criterion = cv(nfolds = 5))
hyper(fit)
#>   parameter                    term   name estimate  held source
#> 1        mu lasso(x, n_lambda = 12) lambda 14.03983 FALSE    bic

# The two columns that carry signal are recovered.
round(coef(fit)$mu[2:3], 3)
#> lasso.1 lasso.2 
#>   1.995  -1.512 

# "1se" takes the largest kink within one standard error of the best, so
# it never keeps more coefficients than "min". On this data the two agree,
# the criterion being flat enough that the best point is already within
# one standard error of itself.
one_se <- statmod(y ~ lasso(x, n_lambda = 12),
                  distributions7::gaussian1_distrib(), dd,
                  outer_criterion = cv(nfolds = 5, rule = "1se"))
c(min = hyper(fit)$estimate, one_se = hyper(one_se)$estimate)
#>      min   one_se 
#> 14.03983 14.03983 
c(min = sum(coef(fit)$mu[-1] != 0), one_se = sum(coef(one_se)$mu[-1] != 0))
#>    min one_se 
#>      7      7 
```
