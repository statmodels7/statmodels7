# The Variance the Priors That Are Not Gaussian Add to Each Predictor

\\z^\top \mathrm{Var}(b)\\ z\\ for every Student t and every other prior
[`prior_parts()`](https://statmodels7.github.io/statmodels7/reference/prior_parts.md)
lists, in each parameter's equation. A Student t with scale \\\Sigma\\
has variance \\\Sigma\nu/(\nu - 2)\\, infinite for \\\nu \le 2\\; a
univariate prior of another family has the variance its family reports;
a prior over several coordinates that is not a Student t has none read
here.

## Usage

``` r
heavy_variance(spec, parts)
```

## Arguments

- spec:

  The specification at the rows predicted.

- parts:

  What
  [`prior_parts()`](https://statmodels7.github.io/statmodels7/reference/prior_parts.md)
  returns.

## Value

A named list of numeric vectors, one value a row, `NA` where the
variance is infinite or not read.
