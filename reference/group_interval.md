# The Interval of a New Group's Parameter

Selects the interval of the parameter of a new group, that is, the
parameter at a level of a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term that is not in the fitting data. It adds the variation of the
effects of a new group to the uncertainty of the estimates. It requires
`random = "zero"` or `"marginal"` in
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
and implies `se = TRUE`.

## Usage

``` r
group_interval()
```

## Value

An S7 object of class `GroupInterval`, inheriting from
[`interval_type()`](https://statmodels7.github.io/statmodels7/reference/interval_type.md),
with the one property `label`.

## Details

The ends are quantiles of the new group's predictor carried through the
inverse link; the standard error is the standard deviation of the new
group's parameter. See
[predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md)
for the formulas and the restrictions.

## See also

[`confidence_interval()`](https://statmodels7.github.io/statmodels7/reference/confidence_interval.md),
[`prediction_interval()`](https://statmodels7.github.io/statmodels7/reference/prediction_interval.md),
[predict_intervals](https://statmodels7.github.io/statmodels7/reference/predict_intervals.md),
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(g = factor(rep(1:8, each = 5)), x = runif(40))
dd$y <- 1 + dd$x + rnorm(8, sd = 0.5)[dd$g] + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
               dd)
predict(fit, what = "mu", newdata = data.frame(x = 0.5, g = "new"),
        random = "zero", interval = group_interval())
#>        fit        se     lower   upper
#> 1 1.470069 0.5548627 0.3825577 2.55758
```
