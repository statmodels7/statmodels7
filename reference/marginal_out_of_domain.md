# Rows Where a New Group's Predictor Leaves the Domain of Its Link

A predictor whose domain is not the whole line (the square root and
inverse links) leaves it with positive probability under a Gaussian
effect, and there the average over new groups does not exist:
\\E\[1/(\eta_0 + u)\]\\ is infinite for a Gaussian \\u\\. Where that
probability is below 1e-8 the average over the nodes is the meaningful
number; above it the row is reported as `NA`, with a warning.

## Usage

``` r
marginal_out_of_domain(link, eta0, Zs, nodes, idx)
```

## Arguments

- link:

  The parameter's link.

- eta0:

  The predictor at the fixed part, one value a row.

- Zs:

  The within-group rows of the marginal terms.

- nodes:

  What
  [`random_nodes()`](https://statmodels7.github.io/statmodels7/reference/random_nodes.md)
  returns.

- idx:

  The marginal terms in this parameter's equation.

## Value

A logical vector, `TRUE` at a row to report as `NA`.

## Details

The variance of the effects' contribution at each row is read off the
nodes, \\\sum_k w_k (z^\top b_k)^2\\, which a Gauss-Hermite grid gives
exactly, the effects having mean zero.
