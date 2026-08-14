# Choose the Hyperparameters by Cross-Validation

`cv()` scores a path of hyperparameter values by the log-likelihood the
fit assigns to observations it was not fitted on.

## Usage

``` r
cv(
  nfolds = 10,
  folds = NULL,
  rule = c("min", "1se"),
  n_values = 25,
  min_ratio = 1e-04,
  over = NULL
)
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

- n_values:

  How many points the path visits.

- min_ratio:

  The smallest kink the path reaches, as a fraction of the one that
  empties the block.

- over:

  Which hyperparameters to sweep. Defaults to the ones that set the size
  of the kink.

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

**Which hyperparameters.** Those whose value sets the size of the kink,
read from the penalty by probing the subdifferential rather than taken
from a list of families: \\\lambda\\ for the lasso, SCAD and MCP, and
the elastic net's \\\lambda\\ with \\\alpha\\ held, as glmnet holds it
and ncvreg holds \\\gamma\\. Name others in `over` to sweep them too.

**The cost** is `nfolds * n_values` fits per hyperparameter. The warm
starts are worth 1.8 times, and building each fold's design once rather
than once per point another 4 per cent, but what remains is the proximal
iteration: measured at 200 observations and 20 columns, 0.88 seconds a
fit, against `cv.glmnet`'s 0.03 seconds for its whole path of 100 values
on five folds. That distance is the reason `n_values` is 25 here and 100
there, and closing it needs the compiled coordinate descent that a
separable penalty on a linear predictor admits.

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
statmod(y ~ lasso(x), distributions7::gaussian1_distrib(), dd,
        outer_criterion = cv(nfolds = 3, n_values = 6))
#> Warning: The path for 'lasso(x)' in 'mu' stopped at its sparse end (lambda = 15.48).
#>   The criterion was still falling there, so widen the path with min_ratio
#>   or set the value yourself.
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ lasso(x), distrib = distributions7::gaussian1_distrib(), 
#>             data = dd, outer_criterion = cv(nfolds = 3, n_values = 6))
#> 
#> Distribution: gaussian1
#> Observations: 60
#> 
#>   mu         ~ lasso(x)
#>                linpar           1 coef
#>                lasso(x)         5 coef, edf 0.00
#>   sigma      ~ 1
#>                linpar           1 coef
#> 
#> log-likelihood -75.244717    objective 65.014193
#> fitted in 919 ms, converged
#> 1 sweeps over 2 block(s)
```
