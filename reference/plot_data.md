# The Data of a Graph

Computes the quantities that a graph draws and returns them as a named
list of data frames, one for each layer. It dispatches on the type and
on `x`: the term for a
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md),
the family of the fit for a
[model_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md)
or a
[diagnostic_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Usage

``` r
plot_data(type, x, fit, ctx, ...)
```

## Arguments

- type:

  An object inheriting from
  [plot_type](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

- x:

  The term or the family on which the method dispatches.

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ctx:

  The context described in the details.

- ...:

  Unused.

## Value

A named list of data frames.

## Details

A method has the signature `function(type, x, fit, ctx, ...)`. The
context `ctx` is a list with the elements `param` (the equation, or
`NULL`), `name` (the name of the term in that equation, or `NULL`),
`vars` (the variables whose effect is drawn), `cols` (the positions
within the term's coefficients that belong to those variables, or `NULL`
for all of them), `data` (the data frame passed to
[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
or `NULL`) and `title` (what `term` named: the term, or the variable).
The list that a method returns may carry the attribute `"labels"`, a
list of arguments of
[`ggplot2::labs()`](https://ggplot2.tidyverse.org/reference/labs.html)
(`title`, `x`, `y` and others), which
[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md)
applies to the graph.

## See also

[`plot_layers()`](https://statmodels7.github.io/statmodels7/reference/plot_layers.md),
[plot_type](https://statmodels7.github.io/statmodels7/reference/plot_type.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(100))
dd$y <- sin(2 * pi * dd$x) + rnorm(100, sd = 0.3)
fit <- statmod(y ~ s(x), distributions7::gaussian1_distrib(), dd)
str(plot_data(plot_effect(), fit@spec@terms$mu[["s(x)"]], fit,
              ctx = list(param = "mu", name = "s(x)", vars = "x",
                         cols = NULL, data = dd)), max.level = 1)
#> List of 4
#>  $ curve :'data.frame':  200 obs. of  4 variables:
#>  $ points:'data.frame':  100 obs. of  2 variables:
#>  $ bins  :'data.frame':  0 obs. of  2 variables:
#>  $ rug   :'data.frame':  100 obs. of  1 variable:
#>  - attr(*, "kind")= chr "numeric"
#>  - attr(*, "labels")=List of 5
```
