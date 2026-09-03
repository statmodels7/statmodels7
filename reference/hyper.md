# The Hyperparameters of a Fitted Model

Reports every hyperparameter of every penalty the model carries: one row
each, with the value, whether it was held, and what put it there. A
smoothing parameter, a prior scale, a lasso's \\\lambda\\ and an elastic
net's \\\alpha\\ all appear here.

## Usage

``` r
hyper(fit, scale = c("parameter", "link"))
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- scale:

  `"parameter"` (the default) or `"link"`. Matched with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

## Value

A data frame with one row per hyperparameter and six columns:

- `parameter`:

  the distribution parameter whose equation the penalized term sits in.

- `term`:

  the term's key, its call as written.

- `name`:

  the hyperparameter's own name, as the penalty names it.

- `estimate`:

  its value, on the scale asked for.

- `held`:

  a logical: whether the term fixed it.

- `source`:

  what put the value there, as above.

- `id`:

  the sharing label, `NA` where the hyperparameter is not shared.

A data frame of no rows, with those columns, where the model carries no
penalty at all.

## A hyperparameter is not a coefficient

It governs the coefficients under it instead of sitting beside them, and
the two are estimated by different routes and reported with different
qualifications.
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
holds the coefficients and this holds the hyperparameters; neither holds
the other.

## The two scales

`scale = "parameter"`, the default, is the scale the penalty is written
on: a smoothing parameter is a positive number, and a gaussian prior's
`sigma` is a scale. That is what a reader wants.

`scale = "link"` is the free scale the outer search runs on, through
each hyperparameter's own link. That is what a caller comparing two
fits' searches wants. Where a hyperparameter carries no link the two
coincide.

## What `source` says that `held` cannot

`held` is a logical and says only whether the value moved. `source` says
what put it there:

- `"fixed"` for one the term itself held, as `s(x, lambda = 2)`.

- the criterion's name, `"reml"` or `"ml"`, for one a marginal criterion
  maximized.

- the criterion that scored the path, `"bic"` and so on, for one chosen
  along a path over its own values.

The distinction has a consequence a reader needs. A value chosen along a
path is the argument of a minimum over a grid, not the root of a
derivative, so no standard error follows from it. One a marginal
criterion reached carries one, from the curvature of that criterion, and
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
prints it.

## Why two rows may carry the same number

Terms whose
[`modelterms7::term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.html)
give the same label estimate one value, and that value is written under
each of their own keys. So a shared hyperparameter appears once per
term, with the same estimate and the same `source`, and the `id` column
is what says the agreement is the model rather than a coincidence. The
penalties are not merged, which is why each still has a row of its own
and why the effective degrees of freedom are still counted term by term.

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

# Estimated by REML, which is what source says.
fit <- statmod(y ~ s(x, k = 6), distributions7::gaussian1_distrib(), d)
hyper(fit)
#>   parameter        term   name estimate  held source   id
#> 1        mu s(x, k = 6) lambda 29.67929 FALSE   reml <NA>

# The same value on the scale the outer search ran on.
hyper(fit, scale = "link")
#>   parameter        term   name estimate  held source   id
#> 1        mu s(x, k = 6) lambda  3.39045 FALSE   reml <NA>

# Held by the term instead, and reported as fixed.
held <- statmod(y ~ s(x, k = 6, lambda = 2),
                distributions7::gaussian1_distrib(), d)
hyper(held)[, c("name", "estimate", "held", "source")]
#>     name estimate held source
#> 1 lambda        2 TRUE  fixed
```
