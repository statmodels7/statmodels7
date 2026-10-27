# Sandwich Variances with a Leverage Correction

Select the meat of the frequentist variance \\V = V_b B V_b\\ with the
scores of the observations rescaled by their leverage: `hc2_sandwich()`
and `hc3_sandwich()`, which with a `cluster` are the cluster-robust
estimators CR2 and CR3, and `hc4_sandwich()`. Pass the object as the
`sandwich` argument of
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md).

## Usage

``` r
hc2_sandwich(cluster = NULL, data = NULL, reference = c("t", "normal"))

hc3_sandwich(cluster = NULL, data = NULL)

hc4_sandwich()
```

## Arguments

- cluster:

  `NULL` (the default) or a one-sided formula with a single variable of
  the data, such as `~ id`.

- data:

  `NULL` (the default) or the data frame that holds the cluster
  variable, with one row for each observation of the fit.

- reference:

  For `hc2_sandwich()` with a `cluster`: `"t"` (the default) or
  `"normal"`, the reference distribution of the Wald intervals and
  tests.

## Value

A `SandwichMethod` object with the properties `type` (`"HC2"`, `"HC3"`
or `"HC4"`), `cluster`, `data` and `reference`.

## Details

Let \\s_i\\ be the vector of the derivatives of the log-likelihood of
observation \\i\\ in the \\d\\ predictors of the model, \\Z_i\\ the \\d
\times P\\ matrix that places the design rows of the observation in the
stacked coefficients, \\W_i\\ the \\d \times d\\ information of the
observation, and \\H\_{ii}\\ its leverage block, as defined at
[`statmod_leverage()`](https://statmodels7.github.io/statmodels7/reference/statmod_leverage.md).
The score of the observation is \\u_i = Z_i' s_i\\. Each estimator
replaces \\s_i\\ by \$\$s_i^\* = W_i^{1/2}\\ (I - H\_{ii})^{-p}\\
W_i^{-1/2}\\ s_i,\$\$ and then forms \\B = \sum_i u_i^\*
u_i^{\*\prime}\\ with \\u_i^\* = Z_i' s_i^\*\\. The estimators differ in
the exponent \\p\\:

- HC2 (`hc2_sandwich()`): \\p = 1/2\\, as in MacKinnon and White (1985);

- HC3 (`hc3_sandwich()`): \\p = 1\\, as in MacKinnon and White (1985);

- HC4 (`hc4_sandwich()`): \\p = \delta_i/2\\ with \\\delta_i = \min\\4,
  \mathrm{tr}(H\_{ii})/\bar h\\\\, where \\\bar h\\ is the mean of
  \\\mathrm{tr}(H\_{ii})\\ over the observations, as in Cribari-Neto
  (2004).

The matrix power is taken on the eigenvalues of the symmetric matrix
\\H\_{ii}\\. When a model has one equation, \\W_i\\ is a scalar and
\\s_i^\* = s_i (1 - h_i)^{-p}\\, the textbook estimators; for a Gaussian
linear model with a constant scale they equal the estimators built from
[`stats::hatvalues()`](https://rdrr.io/r/stats/influence.measures.html).
An eigenvalue of \\H\_{ii}\\ equal to one (a direction that the fit
reproduces exactly) is given the factor zero, which is the
pseudo-inverse, and a warning reports it.

With `cluster = ~ id`, the observations of a cluster \\g\\ are corrected
together: \\s_g^\* = W_g^{1/2}\\(I - H\_{gg})^{-p}\\W_g^{-1/2}\\s_g\\,
where \\s_g\\ stacks the vectors \\s_i\\ of its observations, \\W_g\\ is
block diagonal and \\H\_{gg}\\ is the diagonal block of the hat matrix
for the cluster, of dimension equal to the number of observations in the
cluster times \\d\\. The exponent \\1/2\\ is CR2 (Bell and McCaffrey,
2002; Pustejovsky and Tipton, 2018) and the exponent \\1\\ is CR3. A
cluster with one observation each reduces them to HC2 and HC3.
`hc4_sandwich()` has no `cluster` argument, since no version of HC4 for
clusters is established. The estimators apply no degrees-of-freedom
factor, so that `hc2_sandwich()` and `hc3_sandwich()` coincide with
[`hc0_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md)
when the leverage is zero.

The correction is the exact bias removal for the linearized fit when the
working model is the model itself (the observations independent, with
variance \\W_i^{-1}\\ for their scores), and it uses the information
selected by the `information` argument of
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
when it is given, for the roots \\W_i^{1/2}\\ as well as for the inverse
\\K^{-1}\\. When that argument is `NULL`, the leverage is computed on
the expected information, whatever information the bread \\V_b\\
resolves for the family. With the observed information, an observation
whose matrix \\W_i\\ is indefinite contributes the positive part, with
the pseudo-inverse root
([`statmod_leverage()`](https://statmodels7.github.io/statmodels7/reference/statmod_leverage.md)
returns how many observations that concerns). With penalties the
leverage is that of the penalized fit, so that the sum of the traces is
the effective degrees of freedom of the model. Prior weights, random
effects, aliased coordinates and the models that are rejected are
treated as at
[`hc0_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md).
Prior weights enter the leverage through \\W_i\\; without a cluster the
leverage of one of the \\w_i\\ copies of a row is \\1/w_i\\ times that
of the row, which is what the estimator uses, and with a cluster the
copies of a row are corrected within the cluster.

**Few clusters.** With CR2 and `reference = "t"` (the default),
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
and the Wald rows of
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
and
[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md)
use a \\t\\ distribution whose degrees of freedom are computed for each
coefficient by the Satterthwaite approximation of Pustejovsky and Tipton
(2018). For a coefficient \\c'\beta\\, with \\t_g = (I - H\_{gg})^{-1/2}
W_g^{1/2} Z_g V_b c\\, the variance estimate is \\\sum_g (t_g' \hat
r_g)^2\\ in the standardized scores \\\hat r = W^{-1/2}\hat s\\. If the
true standardized scores \\r\\ had identity covariance, then \\\hat r =
(I - H) r\\ and the estimate is a quadratic form of \\r\\ whose degrees
of freedom are \\(\sum_g q_g'q_g)^2 / \sum\_{g,h} (q_g'q_h)^2\\, where
\\q_g\\ is the vector \\(I - H)\\ t_g\\ with \\t_g\\ placed in the rows
of cluster \\g\\. `reference = "normal"` keeps the normal distribution.
The other sandwich types always use the normal distribution.

## References

MacKinnon, J. G. and White, H. (1985). Some
heteroskedasticity-consistent covariance matrix estimators with improved
finite sample properties. *Journal of Econometrics*, 29(3), 305–325.

Bell, R. M. and McCaffrey, D. F. (2002). Bias reduction in standard
errors for linear regression with multi-stage samples. *Survey
Methodology*, 28(2), 169–181.

Cribari-Neto, F. (2004). Asymptotic inference under heteroskedasticity
of unknown form. *Computational Statistics & Data Analysis*, 45(2),
215–233.

Pustejovsky, J. E. and Tipton, E. (2018). Small-sample methods for
cluster-robust variance estimation and hypothesis testing in fixed
effects models. *Journal of Business & Economic Statistics*, 36(4),
672–683.

## See also

[`hc0_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md),
[`statmod_leverage()`](https://statmodels7.github.io/statmodels7/reference/statmod_leverage.md),
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150), id = rep(1:30, each = 5))
dd$y <- 1 + 2 * dd$x + rnorm(150, sd = 0.2 + dd$x)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
sqrt(diag(vcov(fit, vcov_type = frequentist_vcov(
  sandwich = hc3_sandwich()), param = "mu")))
#> mu:(Intercept)           mu:x 
#>     0.09543413     0.24453148 
vcov(fit, vcov_type = frequentist_vcov(
  sandwich = hc2_sandwich(cluster = ~ id)), param = "mu")
#>                mu:(Intercept)        mu:x
#> mu:(Intercept)     0.01074648 -0.02488978
#> mu:x              -0.02488978  0.07485765
```
