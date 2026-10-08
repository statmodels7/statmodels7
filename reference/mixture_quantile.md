# A Quantile of a Mixture of Gaussians and Point Masses

The value at which \\\sum_c w_c \Phi((e - m_c)/s_c)\\ reaches `p`, one a
row, a component with \\s_c = 0\\ being a point mass at \\m_c\\. Found
by bisection, which needs no derivative and is exact to the tolerance
for a step function as well.

## Usage

``` r
mixture_quantile(M, S, W, p)
```

## Arguments

- M, S:

  Matrices with a row per observation and a column per component.

- W:

  The components' weights, summing to one.

- p:

  The probability.

## Value

One value a row.
