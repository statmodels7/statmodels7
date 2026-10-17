# The Hyperparameters of a Fitted Model

Reports every hyperparameter of every penalty that the model carries:
one row each, with the value, whether it was held, and what put it
there. A smoothing parameter, a prior scale, a lasso's \\\lambda\\ and
an elastic net's \\\alpha\\ all appear here.

## Usage

``` r
hyper(fit, scale = c("parameter", "link"))
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- scale:

  `"parameter"` (the default) or `"link"`. The argument is matched with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

## Value

A data frame with one row per hyperparameter and the columns:

- `parameter`:

  the distribution parameter whose equation contains the penalized term.
  For a covariance class that spans several equations, the names of all
  of them, separated by commas.

- `term`:

  the term's key, its call as written.

- `name`:

  the hyperparameter's own name, as the penalty names it.

- `estimate`:

  its value, on the requested scale.

- `held`:

  a logical: whether the term fixed it.

- `source`:

  what produced the value, as above.

- `id`:

  the sharing label, `NA` where the hyperparameter is not shared.

The result is a data frame of no rows, with those columns, where the
model carries no penalty at all.

## A hyperparameter is not a coefficient

It governs the coefficients under it instead of sitting beside them, and
the two are estimated by different routes and reported with different
qualifications.
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
reports the coefficients, and this function reports the hyperparameters.

## The two scales

`scale = "parameter"`, the default, is the scale on which the penalty is
written: a smoothing parameter is a positive number, and the `sigma` of
a gaussian prior is a scale.

`scale = "link"` is the free scale on which the outer search runs,
obtained through each hyperparameter's own link, and it is the scale on
which the searches of two fits are compared. Where a hyperparameter
carries no link, the two scales coincide.

## The columns `held` and `source`

`held` is a logical that records whether a term fixed the value,
directly or through the label that the hyperparameter shares. `source`
records what produced the value:

- `"fixed"` for a hyperparameter that a term holds, as in
  `s(x, bspline_smooth(), hyper = c(lambda = 2))`.

- the name of the outer criterion (`"reml"`, `"ml"`, `"aic"` or `"bic"`)
  for a hyperparameter that the outer search chose.

- the name of the sparse criterion (`"aic"`, `"bic"` or `"cv"`) for a
  hyperparameter of a kinked penalty chosen along a path over its own
  values.

- `""` where no criterion ran, so that the value stayed at the starting
  value of its term.

The distinction matters for inference. A value chosen by
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md),
whether by the outer search or along a path, is the argument of a
minimum of a prediction-error criterion and not the root of the
derivative of a likelihood, so the value does not yield a standard
error. A value that
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
estimated carries a standard error, computed from the curvature of that
criterion, and
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
prints it.

## Why two rows may carry the same number

Terms whose
[`modelterms7::term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.html)
give the same label estimate one value, and that value is written under
each of their own keys. A shared hyperparameter therefore appears once
per term, with the same estimate and the same `source`, and the `id`
column marks the agreement as part of the model and not as a
coincidence. The penalties are not merged, so each term still has a row
of its own and the effective degrees of freedom are still counted term
by term.

## See also

[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
for the coefficients,
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
for the two printed together with standard errors,
[`statmod_held()`](https://statmodels7.github.io/statmodels7/reference/statmod_held.md)
for which are held.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(80, 0, 1))
d$y <- sin(3 * d$x) + rnorm(80, 0, 0.3)

# The hyperparameter is estimated by REML, which the source column records.
fit <- statmod(y ~ s(x, bspline_smooth(k = 6)),
               distributions7::gaussian1_distrib(), d)
hyper(fit)
#>   parameter                        term   name estimate  held source   id
#> 1        mu s(x, bspline_smooth(k = 6)) lambda 29.31279 FALSE   reml <NA>

# The same value is shown on the scale on which the outer search ran.
hyper(fit, scale = "link")
#>   parameter                        term   name estimate  held source   id
#> 1        mu s(x, bspline_smooth(k = 6)) lambda 3.378024 FALSE   reml <NA>

# Here the term holds the value instead, and it is reported as fixed.
held <- statmod(y ~ s(x, bspline_smooth(k = 6), hyper = c(lambda = 2)),
                distributions7::gaussian1_distrib(), d)
hyper(held)[, c("name", "estimate", "held", "source")]
#>     name estimate held source
#> 1 lambda        2 TRUE  fixed
```
