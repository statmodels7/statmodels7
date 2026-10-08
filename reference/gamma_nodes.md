# Nodes and Weights for a Gamma(a, a) Mixing Variable

A trapezoidal rule in \\v = \log w\\ for the \\\mathrm{Gamma}(a, a)\\
density, over the range between its quantiles \\10^{-12}\\ and \\1 -
10^{-12}\\, the weights normalized to sum to one.

## Usage

``` r
gamma_nodes(a, h = 0.6, k_min = 40L)
```

## Arguments

- a:

  The shape, \\\nu/2\\.

- h:

  The spacing of the nodes in \\\log w\\.

- k_min:

  The smallest number of nodes.

## Value

A list with `w`, the nodes, and `a`, the weights, which sum to one.

## Details

A prediction's distribution function, averaged over \\w\\, depends on
\\w\\ through \\\sqrt{w/(w + c)}\\, which is not smooth at zero, and its
tails come from small \\w\\. A Gauss-Laguerre rule in \\w\\ converges
slowly there; in \\\log w\\ the integrand decays exponentially at both
ends and the trapezoidal rule converges exponentially. Measured on the
probability that a Gaussian measurement around a Student t effect falls
below the 2.5 per cent quantile, against adaptive quadrature: at \\\nu =
1\\ the error is 1e-04 with 40 nodes and 2e-08 with 80, against 1.2e-02
for Gauss-Laguerre with 96; at \\\nu = 2.54\\ it is 1e-07 and 4e-13,
against 5e-05.

The range widens as \\\nu\\ falls (25 units of \\\log w\\ at \\\nu =
2.54\\, 62 at 0.95, 190 at 0.3), so the rule fixes the spacing of the
nodes rather than their number. At a spacing of 0.6 the error is at most
2e-08 from \\\nu = 0.3\\ to 2.54, against 7e-05 with 64 nodes at \\\nu =
0.5\\; at least 40 nodes are used, which covers a large \\\nu\\, where
the range is short.

## Examples

``` r
q <- statmodels7:::gamma_nodes(1.5)
c(length(q$w), sum(q$a), sum(q$a * q$w))
#> [1] 40.000000  1.000000  1.000009
```
