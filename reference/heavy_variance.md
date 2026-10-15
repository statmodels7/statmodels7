# The Variance Added to Each Predictor by Priors Other Than Gaussian

\\z^\top \mathrm{Var}(b)\\ z\\ for every Student t and every other prior
that
[`prior_parts()`](https://statmodels7.github.io/statmodels7/reference/prior_parts.md)
lists, in each parameter's equation. A Student t with scale \\\Sigma\\
has variance \\\Sigma\nu/(\nu - 2)\\, infinite for \\\nu \le 2\\; a
univariate prior of another family has the variance that its family
reports; the variance of a prior over several coordinates that is not a
Student t is not read and is returned as `NA`.

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
