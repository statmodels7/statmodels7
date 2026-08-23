# Choose the Hyperparameters by Cross-Validation

`cv()` scores a path of hyperparameter values by the log-likelihood the
fit assigns to observations it was not fitted on.

## Usage

``` r
cv(nfolds = 10, folds = NULL, rule = c("min", "1se"))
```

## Arguments

- nfolds:

  How many folds. Ignored when `folds` is given.

- folds:

  A fold number per observation, for a partition of your own or to
  compare two criteria on the same one.

- rule:

  `"min"` for the best criterion, `"1se"` for the sparsest fit within
  one standard error of it.

## Value

An
[`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

## Details

**Why a penalty with a kink needs this.** A marginal criterion –
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
– approximates an integral by a Laplace expansion at the penalized mode,
which asks for the second derivative of the penalty there. The mode of a
lasso, a SCAD or an MCP sits at the kink for every coefficient it sets
to zero, which is where that derivative does not exist, so the criterion
is not defined at the point it would be read at. Cross-validation asks a
different question, about prediction rather than about a posterior, and
asking it needs nothing from the penalty beyond a fit.

**The criterion** is the mean over folds of \\-2\ell/n_f\\ on the fold
left out, each training fit rebuilding the design on its own rows.
`rule = "1se"` takes the largest kink whose criterion is within one
standard error of the smallest, which is the sparsest fit that is not
measurably worse.

**The path, not a search.** The penalized mode is a piecewise smooth
function of the hyperparameter: differentiable while the active set
holds, turning a corner whenever a coefficient joins it or leaves. A
criterion read there inherits the corners, so the hyperparameter is
swept over a grid rather than searched by slope. The grid is geometric
in the size of the kink, from the value that leaves the block empty down
to `min_ratio` of it, and every fit starts from the previous one's
coefficients.

**Which hyperparameters** is the TERM's answer: every one of its
penalty's that the constructor did not hold at a number. What the
criterion decides is how they are covered.

**A product within a term, an alternation between terms.** A term
carrying several kinked hyperparameters has every combination of them
visited under `search = "grid"` and one swept at a time under
`"cyclic"`. The choice is the TERM's – `enet(X, search =)`,
[`term_search`](https://statmodels7.github.io/modelterms7/reference/term_search.html)
– and not this criterion's, since the same criterion is put to the
smooth hyperparameters of the model, which are read at the mode rather
than swept, so most of what it is asked about could not use such an
argument. Between two terms the sweep alternates whichever each one
named, so `y ~ lasso(X) + enet(R)` costs the two blocks added and not
multiplied, and one term asking for a product does not make the other
pay for it.

A term that names neither gets the product, because the cyclic sweep
traverses a cross through the point in hand and can stop where each
coordinate is separately best without being jointly so. Its cost is the
product of the term's grids where the cyclic sweep's is their sum, which
at two hyperparameters is `n_lambda * n_alpha` fits against
`n_lambda + n_alpha` per pass; with three or more estimated it grows
exponentially, and `"cyclic"` is there for that.

**How long the path is, and how far down it reaches**, are the term's
too, and are on its own signature where a reader can see them:
`lasso(x, n_lambda = 25, min_ratio = 1e-4)`,
`enet(x, n_lambda = 25, n_alpha = 5)`. The two numbers differ because
the axes do – \\\lambda\\ descends the size of the kink over four
decades and \\\alpha\\ spans one bounded interval – so the shipped
product is 25 by 5.

**The grid is not a rectangle**, and for two different reasons. The
elastic net's kink is \\\lambda\alpha\\, so the value emptying the block
is \\\lambda\_{\max} = \kappa/\alpha\\ and every \\\alpha\\ carries its
own \\\lambda\\ axis, descending from its own top. The shapes of SCAD
and MCP leave the kink alone, so there \\\lambda\_{\max}\\ is one number
whatever the shape and the two axes really are a rectangle. Each axis is
built at the settings of the axes outside it, which covers both without
a rule about families.

**The cost** is `nfolds` fits per point of the path, and how many points
there are is the term's `n_lambda` for one kinked hyperparameter and the
product of its axes for several. The warm starts are worth 1.8 times,
and building each fold's design once rather than once per point another
4 per cent, but what remains is the proximal iteration: measured at 200
observations and 20 columns, 0.88 seconds a fit, against `cv.glmnet`'s
0.03 seconds for its whole path of 100 values on five folds. That
distance is the reason `n_lambda` is 25 here and 100 there, and closing
it needs the compiled coordinate descent that a separable penalty on a
linear predictor admits.

## References

Breiman, L., Friedman, J. H., Olshen, R. A. and Stone, C. J. (1984).
*Classification and Regression Trees*. Wadsworth.

Friedman, J., Hastie, T. and Tibshirani, R. (2010). Regularization paths
for generalized linear models via coordinate descent. *Journal of
Statistical Software* 33(1), 1–22.

## See also

[`reml`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(y = rnorm(60))
dd$x <- matrix(rnorm(60 * 5), 60, 5)
statmod(y ~ lasso(x, n_lambda = 6), distributions7::gaussian1_distrib(), dd,
        outer_criterion = cv(nfolds = 3))
#> Warning: The path for 'lasso(x, n_lambda = 6)' in 'mu' stopped at its sparse end (lambda = 15.48).
#>   The criterion was still falling there, so widen the path with min_ratio
#>   or set the value yourself.
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ lasso(x, n_lambda = 6), distrib = distributions7::gaussian1_distrib(), 
#>             data = dd, outer_criterion = cv(nfolds = 3))
#> 
#> Distribution: gaussian1
#> Observations: 60
#> 
#>   mu         ~ lasso(x, n_lambda = 6)
#>                linpar                   1 coef
#>                lasso(x, n_lambda = 6)   5 coef, edf 0.00
#>   sigma      ~ 1
#>                linpar                   1 coef
#> 
#> fitted in 306 ms   search: converged
#> 1 pass(es) over 2 block(s)
```
