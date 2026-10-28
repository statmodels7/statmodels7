# Graph of the Derivative of a Smooth

Selects the graph of the first or second derivative of the effect of a
smooth of one numeric variable, with a confidence band.

## Usage

``` r
plot_derivative(
  order = 1,
  level = 0.95,
  simultaneous = FALSE,
  n_sim = 10000,
  seed = NULL,
  vcov_type = bayesian_vcov(),
  n_grid = 200
)
```

## Arguments

- order:

  The order of the derivative, 1 or 2.

- level:

  The confidence level of the band.

- simultaneous:

  Whether the band is simultaneous over the grid.

- n_sim:

  The number of draws that give the simultaneous critical value.

- seed:

  `NULL`, or a single number that seeds the draws; the state of the
  random number generator is restored afterwards.

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object, the variance of the band.

- n_grid:

  The number of points of the grid of the variable.

## Value

A `PlotDerivative` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

The effect of the smooth is \\f(x) = X(x) b\\, where \\X(x)\\ is the row
of the block of the term at \\x\\ and \\b\\ its coefficients, so that
\\f^{(m)}(x) = X^{(m)}(x) b\\ is linear in \\b\\. The rows
\\X^{(m)}(x)\\ are obtained from the design of the term by one
finite-difference stencil of order \\m\\ and accuracy four, with the
step and the weights of numericals7
([`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html),
[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html));
near the ends of the range of the variable the nodes are shifted inside
the range. The standard error of the derivative is the exact
delta-method value \\\sqrt{X^{(m)}(x)\\ V\\ X^{(m)}(x)'}\\, where \\V\\
is the block of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
for the coefficients of the term. The graph has the derivative with
respect to the variable on the axis, on the scale of the link.

With `simultaneous = FALSE` the band is the derivative \\\pm
z\_{(1+\text{level})/2}\\ standard errors. With `simultaneous = TRUE`
the factor \\z\\ is replaced by the `level` quantile of \\\max_x
\|X^{(m)}(x)\\ e\| / \mathrm{se}(x)\\, where \\e\\ is drawn `n_sim`
times from \\N(0, V)\\, so that the band covers the whole derivative
over the grid with probability `level` (the simultaneous bands of
Ruppert, Wand and Carroll). A smooth with a `by` factor is drawn as one
curve for each level, each with its own critical value.

## See also

[`plot_effect()`](https://statmodels7.github.io/statmodels7/reference/plot_effect.md),
[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150))
dd$y <- sin(2 * pi * dd$x) + rnorm(150, sd = 0.3)
fit <- statmod(y ~ s(x), distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, term = "s(x)", plot_type = plot_derivative(), data = dd)

ggplot2::autoplot(fit, term = "s(x)",
                  plot_type = plot_derivative(order = 2, simultaneous = TRUE,
                                         n_sim = 2000, seed = 1),
                  data = dd)
```
