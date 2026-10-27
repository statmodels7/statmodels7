# The Layers of a Graph

Turns the list returned by
[`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md)
into a list of ggplot2 components (layers, scales, guides), which
[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md)
adds to an empty
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).
It dispatches on the type.

## Usage

``` r
plot_layers(type, data, ...)
```

## Arguments

- type:

  An object inheriting from
  [plot_type](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

- data:

  The list returned by
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

- ...:

  Unused.

## Value

A list of ggplot2 components.

## See also

[`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md),
[plot_type](https://statmodels7.github.io/statmodels7/reference/plot_type.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(100))
dd$y <- sin(2 * pi * dd$x) + rnorm(100, sd = 0.3)
fit <- statmod(y ~ s(x), distributions7::gaussian1_distrib(), dd)
d <- plot_data(plot_effect(), fit@spec@terms$mu[["s(x)"]], fit,
               ctx = list(param = "mu", name = "s(x)", vars = "x",
                          cols = NULL, data = dd))
ggplot2::ggplot() + plot_layers(plot_effect(), d)
```
