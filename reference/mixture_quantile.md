# A Quantile of a Mixture of Gaussians and Point Masses

Returns the value at which \\\sum_c w_c \Phi((e - m_c)/s_c)\\ reaches
`p`, one value per row; a component with \\s_c = 0\\ is a point mass at
\\m_c\\. The value is found by bisection, which needs no derivative and
is exact to the tolerance for a step function as well. A row in which
some mean or standard deviation is not finite (a missing covariate) is
`NA`, and the bisection runs over the other rows.

## Usage

``` r
mixture_quantile(M, S, W, p)
```

## Arguments

- M, S:

  Matrices with a row per observation and a column per component.

- W:

  The components' weights, rescaled to sum to one.

- p:

  The probability.

## Value

One value a row, `NA` for a row that has no finite mixture.
