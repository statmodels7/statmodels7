# Prediction-Error Criteria for the Hyperparameters

`aic()` and `bic()` choose the hyperparameters by an estimate of
prediction error rather than by a marginal likelihood.

## Usage

``` r
aic(k = 2, hessian = c("observed", "expected"))

bic(hessian = c("observed", "expected"))
```

## Arguments

- k:

  The price of one degree of freedom. Defaults to 2; `bic()` uses \\\log
  n\\, resolved when the model is fitted.

- hessian:

  Which information is used, `"expected"` or `"observed"`. The exact
  derivatives need the observed one.

## Value

An
[`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

## Details

**The criterion** is \$\$C(\theta) = -2\ell(\hat\beta(\theta)) +
\kappa\\\tau(\theta),\$\$ with \\\tau = \mathrm{tr}\[(H+S)^{-1}H\]\\ the
effective degrees of freedom and \\\kappa\\ the price of one of them:
\\2\\ for `aic()` and \\\log n\\ for `bic()`, which is `aic()` with that
\\\kappa\\ and is offered separately because it is what anybody would
look for.

**Both derivatives are exact**, by the implicit function theorem, and
come from the same pieces the marginal criterion uses. One thing differs
and it is what separates the two families: the envelope theorem does not
apply here. \\\ell\\ alone is not stationary at the penalized mode, so
its derivative carries \\d\hat\beta/d\theta\\ from the first order, and
what makes it computable is that at the mode
\\\partial\ell/\partial\beta\\ is \\\partial\rho/\partial\beta\\.

**Which \\\tau\\.** The trace is taken over the whole coefficient
vector, so a term's contribution reads the full penalized information
and not only its own block.
[`statmod_edf`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)'s
per-term numbers invert the block instead, which is what a per-term
reading has to do, and the two agree when the blocks are orthogonal and
differ when they are not.

**GCV is not among these**, and the reason is the framework rather than
the work: classical GCV divides a residual sum of squares by
\\(n-\tau)^2\\, which estimates an unknown scale. Here every
distribution parameter has its own equation and nothing is profiled, so
there is no unknown scale for that ratio to estimate, and the criterion
it degenerates to is `aic()` – the substitution Wood (2008) makes in the
other direction when the scale is known. A GCV on the squared error of
the fitted mean is a different and well-defined object; it needs the
derivative of that mean in the parameters, which is not one of
distributions7's generics.

**Over a penalty with a kink** these criteria sweep a path rather than
differentiate. How that path covers a term carrying more than one such
hyperparameter is the TERM's own,
[`term_search`](https://statmodels7.github.io/modelterms7/reference/term_search.html),
since the same criterion is asked of the smooth hyperparameters of the
model as well and those are not swept at all.

## References

Wood, S. N. (2008). Fast stable direct fitting and smoothness selection
for generalized additive models. *Journal of the Royal Statistical
Society, Series B*, 70(3), 495–518.

## See also

[`reml`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(200, -2, 2))
dd$y <- sin(1.4 * dd$x) + rnorm(200, sd = 0.3)
statmod(y ~ s(x, k = 10), distributions7::gaussian1_distrib(), dd,
        outer_criterion = aic())
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ s(x, k = 10), distrib = distributions7::gaussian1_distrib(), 
#>             data = dd, outer_criterion = aic())
#> 
#> Distribution: gaussian1
#> Observations: 200
#> 
#>   mu         ~ s(x, k = 10)
#>                linpar           1 coef
#>                s(x, k = 10)     9 coef, edf 7.92
#>   sigma      ~ 1
#>                linpar           1 coef
#> 
#> log-likelihood -33.778973    objective 44.594796
#> AIC 87.386042 over 5 hyperparameter evaluation(s)
#> fitted in 550 ms, converged
```
