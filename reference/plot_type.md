# Graphical Types for a Fitted Model

`plot_type` is the abstract S7 class of the objects that select a graph
of a
[`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).
A graph is drawn by
[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md)
(which returns it) or
[`plot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md)
(which prints it), and the object passed as `type` decides which graph
is drawn and with which options. Three abstract subclasses divide the
types by what they read:

- `term_plot`: one term of one equation, named by the argument `term`
  ([`plot_effect()`](https://statmodels7.github.io/statmodels7/reference/plot_effect.md),
  [`plot_derivative()`](https://statmodels7.github.io/statmodels7/reference/plot_derivative.md),
  [`plot_shrinkage()`](https://statmodels7.github.io/statmodels7/reference/plot_shrinkage.md)
  and others);

- `model_plot`: the fitted distribution as a function of the covariates,
  without a term
  ([`plot_param()`](https://statmodels7.github.io/statmodels7/reference/plot_param.md),
  [`plot_quantiles()`](https://statmodels7.github.io/statmodels7/reference/plot_quantiles.md),
  [`plot_density()`](https://statmodels7.github.io/statmodels7/reference/plot_density.md),
  [`plot_moments()`](https://statmodels7.github.io/statmodels7/reference/plot_moments.md));

- `diagnostic_plot`: the fit against the data, through the quantile
  residuals
  ([`plot_qq()`](https://statmodels7.github.io/statmodels7/reference/plot_qq.md),
  [`plot_worm()`](https://statmodels7.github.io/statmodels7/reference/plot_worm.md),
  [`plot_calibration()`](https://statmodels7.github.io/statmodels7/reference/plot_calibration.md)
  and others).

## Usage

``` r
plot_type(label = character(0))

term_plot(label = character(0))

model_plot(label = character(0))

diagnostic_plot(label = character(0))
```

## Arguments

- label:

  A single string, the name of the graph.

## Value

Nothing: the classes are abstract and their constructors signal an
error. Every object that inherits from `plot_type` has the property
`label`.

## Details

A graph is produced in two steps, each an S7 generic.
[`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md)
computes the quantities to draw and returns them as a list of data
frames; it dispatches on the type and on a second object, the term for a
`term_plot` and the family for the other two.
[`plot_layers()`](https://statmodels7.github.io/statmodels7/reference/plot_layers.md)
turns that list into ggplot2 layers and dispatches on the type alone. A
new type is a subclass of one of the three classes, with a method for
each generic.

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md),
[`plot_layers()`](https://statmodels7.github.io/statmodels7/reference/plot_layers.md),
[`default_plot()`](https://statmodels7.github.io/statmodels7/reference/default_plot.md)

## Examples

``` r
S7::S7_inherits(plot_effect(), plot_type)
#> [1] TRUE
S7::S7_inherits(plot_effect(), term_plot)
#> [1] TRUE
try(plot_type(label = "mine"))
#> Error in S7::new_object(S7::S7_object(), label = label) : 
#>   Can't construct an object from abstract class <plot_type>
```
