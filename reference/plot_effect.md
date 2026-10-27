# Graph of the Effect of a Term

Selects the graph of the contribution of a term, or of a variable, to
the predictor of its equation, with a pointwise confidence band and the
partial residuals. For a numeric variable the contribution is drawn as a
curve over the range of the variable; for a factor, as one estimate with
its interval per level.

## Usage

``` r
plot_effect(
  residuals = TRUE,
  binned = NULL,
  scale = "link",
  at = NULL,
  rug = TRUE,
  level = 0.95,
  vcov_type = bayesian_vcov(),
  information = distributions7::expected_info(),
  n_grid = 200
)
```

## Arguments

- residuals:

  Whether to draw the partial residuals.

- binned:

  `NULL`, or the number of classes of the variable within which the
  means of the partial residuals are drawn.

- scale:

  `"link"` (the default), the contribution to the predictor, or
  `"response"`, the parameter.

- at:

  A named list of values for the other covariates on the response scale,
  or `NULL`.

- rug:

  Whether to draw the observed values of the variable on the horizontal
  axis.

- level:

  The confidence level of the band.

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object, the variance of the band.

- information:

  An object from
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  (the default) or
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html),
  the weight of the working response.

- n_grid:

  The number of points of the grid of a numeric variable.

## Value

A `PlotEffect` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

On the scale of the link the contribution of a term is \\f(x) = X b\\,
evaluated on a grid of the variable, and the band is \\f(x) \pm
z\_{(1+\text{level})/2}\\\sqrt{x' V x}\\, where \\V\\ is the block of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the coefficients of the term. A smooth is drawn as estimated, under
the constraint that identifies it. The partial residuals are \\f(x_i) +
(z_i - \eta_i)\\, where \\z_i = \eta_i + s_i / w_i\\ is the working
response of the equation, \\s_i\\ the score with respect to the
predictor and \\w_i\\ the information selected by `information`. For a
gaussian mean with the identity link they are the ordinary partial
residuals. For other equations they are often skewed (for a scale
parameter, the cloud lies below the curve), and `binned` adds the means
of the partial residuals within classes of the variable.

With `scale = "response"` the graph shows the parameter itself,
\\g^{-1}(\eta)\\, from
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
with the other covariates at `at` (or at their median and first level)
and the random effects at zero; on this scale the effect depends on the
other covariates, and the partial residuals are not drawn.

Other variables that the term reads are held at their median or first
level, except for a factor with at most twelve levels, which is crossed
with the grid and drawn as one curve per level (`s(x, by = g)`).

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_derivative()`](https://statmodels7.github.io/statmodels7/reference/plot_derivative.md),
[`plot_surface()`](https://statmodels7.github.io/statmodels7/reference/plot_surface.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(200), z = rnorm(200))
dd$y <- sin(2 * pi * dd$x) + 0.5 * dd$z + rnorm(200, sd = exp(-1 + dd$x))
fit <- statmod(y ~ z + s(x) | sigma ~ x, distributions7::gaussian1_distrib(),
               dd)
ggplot2::autoplot(fit, term = "s(x)", type = plot_effect(), data = dd)

ggplot2::autoplot(fit, term = "x", param = "sigma",
                  type = plot_effect(binned = 10), data = dd)

ggplot2::autoplot(fit, term = "x", param = "sigma",
                  type = plot_effect(scale = "response"), data = dd)
```
