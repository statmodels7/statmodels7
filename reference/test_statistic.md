# Statistics for Testing a Coefficient

Four statistics test the hypothesis that one coefficient of a fit equals
a value: Wald's, the likelihood ratio, Rao's score and Terrell's
gradient. Each is an object built by
[`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md),
[`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md),
[`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md)
or
[`gradient_test()`](https://statmodels7.github.io/statmodels7/reference/gradient_test.md),
and the object is passed as the `test` argument of
[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
and
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md).
The default is
[`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md).

`test_statistic` is the abstract parent of the four classes. It cannot
be instantiated, and serves to test membership with
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html).

## Usage

``` r
test_statistic(name = character(0), label = character(0))
```

## Arguments

- name:

  A single string, the short name of the statistic (`"wald"`, `"lr"`,
  `"score"` or `"gradient"`), which also gives the name of the
  constructor
  ([`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md)
  and so on).

- label:

  A single string, the name under which a result prints the statistic
  (for example `"Likelihood-ratio"`).

## Value

Nothing: the class is abstract and its constructor signals an error.
Every object that inherits from it has the properties `name` and
`label`.

## Details

Write \\\ell_p\\ for the penalized log-likelihood (the log-likelihood
where the model carries no penalty), \\\hat\beta\\ for the estimates of
the fit, \\\tilde\beta\\ for the estimates in which the coefficient
\\\beta_j\\ is held at the tested value \\b\\, \\U = \partial\ell_p\\
for the penalized score and \\K\\ for the penalized information. The
four statistics are \$\$W = \frac{(\hat\beta_j -
b)^2}{\widehat{\mathrm{Var}}(\hat\beta_j)}, \qquad LR =
2\\\ell_p(\hat\beta) - \ell_p(\tilde\beta)\\,\$\$ \$\$S =
U_j(\tilde\beta)^2\\\[K(\tilde\beta)^{-1}\]\_{jj}, \qquad T =
U_j(\tilde\beta)\\(\hat\beta_j - b),\$\$ and each is compared with a
\\\chi^2_1\\ distribution. In
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
the value of a table is the signed root of the statistic, with the sign
of the estimate, and the p-value is that of the \\\chi^2_1\\ value.

|  |  |  |  |
|----|----|----|----|
| Constructor | Statistic | What it needs | Cost |
| [`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md) | \\W\\ | the unrestricted fit and its variance matrix | no refit |
| [`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md) | \\LR\\ | the restricted refit and the objective at the unrestricted joint mode | one refit per value \\b\\ |
| [`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md) | \\S\\ | the restricted refit and its penalized information | one refit and one Hessian per value \\b\\ |
| [`gradient_test()`](https://statmodels7.github.io/statmodels7/reference/gradient_test.md) | \\T\\ | the restricted refit and its score | one refit per value \\b\\, no matrix inversion |

The three statistics that need the restricted refit are the restricted
tests, and
[`test_restricted()`](https://statmodels7.github.io/statmodels7/reference/test_restricted.md)
returns `TRUE` for them. The refit maximizes \\\ell_p\\ over every
coefficient but \\\beta_j\\, at the hyperparameters of the fit, with the
inner method that the fit used, and it starts from the unrestricted
estimates; the outer search over the hyperparameters is not repeated.
The statistics are therefore conditional on the smoothing that the data
chose once. Their cost in each consumer is as follows.

- [`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md)
  runs one refit.

- [`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
  runs one refit per coefficient row, against zero.

- [`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
  inverts the statistic, which takes a few dozen refits per interval,
  and it requires `readable = FALSE` for a restricted test, because a
  readable quantity is a function of several coefficients at once and
  one coordinate cannot be held in its place.

The restricted refit is defined for a coefficient that is not under a
kinked penalty (a lasso, a SCAD or an MCP), that is not aliased and that
belongs to a model without a structural term. For any other coefficient
the restricted statistics are `NA`.

The Wald statistic depends on how the coefficient is written: a
nonlinear reparametrization changes it, whereas the likelihood ratio, a
difference of objective values at the same two points, is not changed.
It also loses power against a distant alternative (the Hauck-Donner
effect), because the curvature is read at \\\hat\beta\\ and not under
the null. The interval that inverts the Wald statistic is symmetric
about the estimate. The intervals that invert the other three follow the
shape of the likelihood and are asymmetric where it is. The variance
that the Wald statistic reads is set by the argument `vcov_type`,
described at
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md).

## See also

[`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md),
[`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md),
[`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md)
and
[`gradient_test()`](https://statmodels7.github.io/statmodels7/reference/gradient_test.md)
for the four statistics,
[`test_restricted()`](https://statmodels7.github.io/statmodels7/reference/test_restricted.md)
for the distinction between them,
[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
and
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
for the functions that take one, and
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
for the variance that the Wald statistic reads.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(80))
dd$y <- 1 + 2 * dd$x + rnorm(80, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)

# the four statistics of one hypothesis, that the slope is 1.5
tests <- list(wald_test(), lr_test(), score_test(), gradient_test())
vapply(tests, function(tt) {
  statmod_test(fit, param = "mu", coefname = "x", value = 1.5,
               test = tt)@statistic
}, numeric(1))
#> [1] 11.67136 11.15541 10.41256 10.41256

# which of them refit the model under the restriction
vapply(tests, test_restricted, logical(1))
#> [1] FALSE  TRUE  TRUE  TRUE

# the class itself is abstract
try(test_statistic(name = "mine", label = "Mine"))
#> Error in S7::new_object(S7::S7_object(), name = name, label = label) : 
#>   Can't construct an object from abstract class <test_statistic>
```
