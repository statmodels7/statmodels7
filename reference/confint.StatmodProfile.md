# Profile Likelihood Intervals

Inverts the profile of each label. Writing \\\hat\ell\\ for the maximum
of the criterion and \\\ell_p(\theta)\\ for the profile at the value
\\\theta\\, the interval is the set of values at which the likelihood
ratio statistic \\2\\\hat\ell - \ell_p(\theta)\\\\ does not exceed the
quantile of a chi-squared distribution.

## Usage

``` r
# S3 method for class 'StatmodProfile'
confint(
  object,
  parm = NULL,
  level = 0.95,
  type = c("marginal", "simultaneous"),
  ...
)
```

## Arguments

- object:

  A
  [StatmodProfile](https://statmodels7.github.io/statmodels7/reference/StatmodProfile.md)
  object.

- parm:

  The labels, a subset of `object@profiled_labels`; by default all of
  them.

- level:

  The coverage.

- type:

  `"marginal"` or `"simultaneous"`.

- ...:

  Unused.

## Value

A data frame with one row for each label, named by it, and the columns
`estimate`, `lower` and `upper`, on the scale of the profile.

## Details

Between two points of the profile the limit is found by a cubic spline
through the square root of the statistic, taken with the sign of
\\\theta - \hat\theta\\, which is close to linear in \\\theta\\ where
the profile is close to quadratic. A limit outside the values of the
profile is `NA`, and the function warns.

With `type = "marginal"` the quantile has one degree of freedom, and
each interval covers its label with probability `level`. With
`type = "simultaneous"` it has \\k\\ degrees of freedom, \\k\\ being the
number of labels of the profile. Each interval is then the projection on
its label of the joint region in which the statistic of the \\k\\ labels
does not exceed \\\chi^2_k(\mathrm{level})\\, because the profile of one
label is the maximum of the criterion over the other labels, and the
\\k\\ intervals cover the labels jointly with probability at least
`level`.

The profile of the position of a sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
is constant between two consecutive values of the covariate. Its
interval runs from the value of the covariate below the first interval
of values in which the statistic does not exceed the quantile to the
value above the last one.

The one-dimensional profile of a label is exact when the labels are
profiled one at a time (`joint = FALSE` in
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md)).
In a joint profile it is the maximum over the grid of the other labels,
so its accuracy is limited by the spacing of that grid.

The function signals an error when the profile does not give a
likelihood ratio statistic (the property `is_likelihood_ratio`): for the
hyperparameter of a kinked penalty that a path chooses, and for a
coefficient without a penalty held under
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
where
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
gives the likelihood ratio.

## See also

[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md),
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
prof <- statmod_profile(fit, which = "mu:x",
                        values = profile_range(width = 3, n = 13))
confint(prof)
#>      estimate    lower    upper
#> mu:x 1.956143 1.622843 2.289442
```
