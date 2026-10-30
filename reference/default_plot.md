# The Default Graph of a Term

Returns the
[plot_type](https://statmodels7.github.io/statmodels7/reference/plot_type.md)
object that
[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md)
uses for a term when `plot_type` is not given. The method for
[modelterms7::model_term](https://statmodels7.github.io/modelterms7/reference/model_term.html)
returns
[`plot_effect()`](https://statmodels7.github.io/statmodels7/reference/plot_effect.md),
and a class of term with a more informative graph has a method of its
own.

## Usage

``` r
default_plot(term, ...)
```

## Arguments

- term:

  A model term, an object inheriting from
  [modelterms7::model_term](https://statmodels7.github.io/modelterms7/reference/model_term.html).

- ...:

  Unused.

## Value

An object inheriting from
[plot_type](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## See also

[plot_type](https://statmodels7.github.io/statmodels7/reference/plot_type.md),
[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md)

## Examples

``` r
default_plot(modelterms7::s(x))
#> <statmodels7::PlotEffect>
#>  @ label      : chr "plot_effect()"
#>  @ legend     : logi TRUE
#>  @ residuals  : logi TRUE
#>  @ binned     : NULL
#>  @ scale      : chr "link"
#>  @ at         : NULL
#>  @ rug        : logi TRUE
#>  @ level      : num 0.95
#>  @ vcov_type  : <statmodels7::BayesianVcov>
#>  .. @ name       : chr "bayesian"
#>  .. @ information: NULL
#>  @ information: <distributions7::ExpectedInfo>
#>  .. @ method: NULL
#>  @ n_grid     : num 200
```
