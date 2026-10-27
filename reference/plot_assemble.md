# Assemble a Graph from Its Data and Its Layers

Calls
[`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md)
and
[`plot_layers()`](https://statmodels7.github.io/statmodels7/reference/plot_layers.md)
and adds the theme, the layers and the labels recorded on the data to an
empty
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).
The theme comes first, so that a list of layers may also carry
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
elements.

## Usage

``` r
plot_assemble(type, x, fit, ctx)
```

## Arguments

- type:

  An object inheriting from
  [plot_type](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

- x:

  The term or family on which
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md)
  dispatches.

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ctx:

  The context described at
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

## Value

A ggplot2 object.
