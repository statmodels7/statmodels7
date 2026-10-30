# The Profile Likelihood Returned by statmod_profile()

The object that
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md)
returns.

## Usage

``` r
StatmodProfile(
  points = NULL,
  profiled_labels = character(0),
  profiled_quantities = list(),
  profile_column = character(0),
  unrestricted_point = NULL,
  profile_method = profile_method(),
  hyper_handling = character(0),
  outer_method = NULL,
  is_likelihood_ratio = logical(0),
  joint = logical(0)
)
```

## Arguments

- points:

  The data frame of the points.

- profiled_labels:

  The labels of `which`.

- profiled_quantities:

  A list of
  [profiled_quantity](https://statmodels7.github.io/statmodels7/reference/profiled_quantity.md)
  objects, one for each label.

- profile_column:

  The name of the column of `points` that is the profile of each label,
  `"outer_criterion"` or `"profile_pen_loglik"`, a character vector
  named by the labels.

- unrestricted_point:

  The data frame of the unrestricted fit.

- profile_method:

  The
  [profile_method](https://statmodels7.github.io/statmodels7/reference/profile_method.md)
  object.

- hyper_handling:

  `"held"` or `"refit"`, the argument `hyper` of
  [`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md).

- outer_method:

  The outer criterion object of the fit, or `NULL`.

- is_likelihood_ratio:

  Whether the profile of each label gives a likelihood ratio statistic,
  so that confint() inverts it, a logical vector named by the labels.

- joint:

  Whether the labels were held jointly.

## Value

A `StatmodProfile` object.

## Details

With `joint = TRUE` the table `points` has one column for each label,
named by the label, and one row for each combination of values held.
With `joint = FALSE` each label is profiled by itself, with the other
labels free, and `points` has the columns `profiled_label`, the label
held in the row, and `value`, its value. The columns that follow are, in
both cases, `profile_loglik`, `profile_penalty`, `profile_pen_loglik`,
`outer_criterion` and `converged`, described at
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md).
`unrestricted_point` has the same columns at the fit itself, with one
row for each label when `joint = FALSE`.

## See also

[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md),
[`confint.StatmodProfile()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodProfile.md),
[`autoplot.StatmodProfile()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodProfile.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
prof <- statmod_profile(fit, which = "mu:x",
                        values = profile_range(width = 2, n = 9))
prof@profile_column
#>                 mu:x 
#> "profile_pen_loglik" 
head(as.data.frame(prof), n = 3)
#>       mu:x profile_loglik profile_penalty profile_pen_loglik outer_criterion
#> 1 1.615742      -23.56354               0          -23.56354       -25.70965
#> 2 1.700842      -22.70458               0          -22.70458       -24.86501
#> 3 1.785942      -22.07563               0          -22.07563       -24.24654
#>   converged
#> 1      TRUE
#> 2      TRUE
#> 3      TRUE
```
