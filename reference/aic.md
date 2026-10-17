# Prediction-Error Criteria for the Hyperparameters

Build the object that specifies, for
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
that a model's hyperparameters are chosen by an estimate of prediction
error. Pass the result as `outer_criterion` for the smooth
hyperparameters, as `sparse_criterion` for a lasso's or a SCAD's, or as
both.

The criterion is \$\$C(\theta) = -2\ell(\hat\beta(\theta)) +
\kappa\\\tau(\theta),\$\$ with \\\tau = \mathrm{tr}\[(H+S)^{-1}H\]\\ the
effective degrees of freedom and \\\kappa\\ the price of one of them.
`aic()` charges 2 by default and `bic()` charges \\\log n\\, with \\n\\
the number of observations (the total of the weights, which are
frequencies). Both are minimized.

## Usage

``` r
aic(k = 2, hessian = c("observed", "expected"))

bic(hessian = c("observed", "expected"))
```

## Arguments

- k:

  The price of one degree of freedom, a single non-negative finite
  number. Defaults to 2, which is Akaike's. `k = log(n)` by hand is
  `bic()`. A value that is not a single non-negative finite number is an
  error. `bic()` does not take it: its \\\kappa\\ is `NA` on the object
  and is resolved to \\\log n\\ when the model is fitted, because the
  sample size is unknown until then.

- hessian:

  The information on which the criterion is built, `"observed"` (the
  default) or `"expected"`. The argument is matched with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html). The exact
  gradient is available on either wherever the family writes its
  expected information out; the Hessian is exact on the observed one and
  is not available on the expected one, where the search uses an
  optimizer that needs only the gradient. A family that does not write
  it out would read the outer product of its scores at the data under
  the default `iwls(approx = "opg")`, and
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  rejects `"expected"` there: see
  [`assert_criterion_information()`](https://statmodels7.github.io/statmodels7/reference/assert_criterion_information.md).

## Value

An
[`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
object with properties `kind` (`"aic"` or `"bic"`), `hessian` as
supplied, `k` (the number given, or `NA_real_` for `bic()`), `marginal`
(`NULL`), and the path settings `nfolds`, `rule` and `folds`, which only
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
reads.

## Marginal and prediction-error criteria

[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
integrate the coefficients out and maximize what is left, which needs a
penalty that is twice differentiable at the mode. A kinked penalty puts
the mode on the kink for every coefficient that it sets to zero, where
the penalized information does not exist. A marginal criterion therefore
leaves the coordinates that such a penalty covers out of its determinant
and does not select the hyperparameter of the kinked penalty. These two
criteria select it, and `bic()` is
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)'s
default `sparse_criterion`.

For a smooth penalty both families apply. A marginal criterion selects
the hyperparameter that makes the data most likely after the
coefficients are integrated out, and these select the one that predicts
best.

## Exact derivatives

The derivatives follow from the implicit function theorem, using the
same pieces as the marginal criteria. The envelope theorem does not
apply here: \\\ell\\ alone is not stationary at the penalized mode, so
its derivative carries \\d\hat\beta/d\theta\\ from the first order. This
is computable because of the mode's own condition,
\\\partial\ell/\partial\beta = \partial\rho/\partial\beta\\.

The gradient is exact with either information. The Hessian is exact with
`hessian = "observed"`, which is the default, and is not available with
`"expected"`.

## Which effective degrees of freedom

The trace runs over the whole coefficient vector, so a term's
contribution reads the full penalized information, including the term's
couplings with the other blocks.
[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
reads each term's share from the diagonal block of the same full
smoother matrix \\(H+S)^{-1}H\\, so the per-term numbers add up to
\\\tau\\.

Where a penalty has a kink the trace is restricted to the active
coordinates, so for a lasso \\\tau\\ is the number of surviving
coefficients, the unbiased count of Zou, Hastie and Tibshirani (2007).

## GCV is not offered

Classical GCV divides a residual sum of squares by \\(n-\tau)^2\\, which
estimates an unknown scale. Here every distribution parameter carries
its own equation and the fit profiles no parameter, so there is no
unknown scale for that ratio to estimate and the criterion degenerates
to `aic()`. That is the substitution that Wood (2008) makes in the other
direction when the scale is known.

A GCV on the squared error of the fitted mean is a different and
well-defined object. It would need the derivative of that mean in the
parameters, which is not one of distributions7's generics.

## Kinked penalties

The penalized mode is only piecewise smooth in such a hyperparameter,
turning a corner whenever a coefficient joins the active set or leaves
it, so the criterion inherits the corners and is swept over a grid. How
the path covers a term carrying more than one kinked hyperparameter
belongs to the term, through
[`modelterms7::term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.html).
The same criterion object also selects the model's smooth
hyperparameters, which are not swept.

## References

Zou, H., Hastie, T. and Tibshirani, R. (2007). On the degrees of freedom
of the lasso. *The Annals of Statistics* **35**(5), 2173–2192.

Wood, S. N. (2008). Fast stable direct fitting and smoothness selection
for generalized additive models. *Journal of the Royal Statistical
Society, Series B* **70**(3), 495–518.

## See also

[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
for the marginal criteria,
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) for
cross-validation,
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for where these are passed.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(200, -2, 2))
dd$y <- sin(1.4 * dd$x) + rnorm(200, sd = 0.3)

fa <- statmod(y ~ s(x, bspline_smooth(k = 10)), distributions7::gaussian1_distrib(), dd,
              outer_criterion = aic())
fb <- statmod(y ~ s(x, bspline_smooth(k = 10)), distributions7::gaussian1_distrib(), dd,
              outer_criterion = bic())

# BIC charges log(n) = 5.3 per degree of freedom against AIC's 2, so it
# buys a smoother fit: a larger smoothing parameter and fewer edf.
c(aic = unlist(fa@hyper), bic = unlist(fb@hyper))
#> aic.mu.s(x, bspline_smooth(k = 10)).lambda 
#>                                  0.5696677 
#> bic.mu.s(x, bspline_smooth(k = 10)).lambda 
#>                                 13.4428649 
c(aic = sum(fa@edf$edf), bic = sum(fb@edf$edf))
#>      aic      bic 
#> 9.914632 7.258561 

# The object is a specification and carries no data.
aic()
#> <AIC>  observed information
aic(k = 4)@k
#> [1] 4
```
