# Nodes and Weights for a Gamma(a, a) Mixing Variable

Returns a trapezoidal rule in \\v = \log w\\ for the \\\mathrm{Gamma}(a,
a)\\ density, over the range between its quantiles \\10^{-12}\\ and
\\1 - 10^{-12}\\. The weights are normalized to sum to one.

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
ends and the trapezoidal rule converges exponentially.

The range of \\\log w\\ widens as \\\nu\\ falls, so the rule fixes the
spacing `h` of the nodes and not their number. At least `k_min` nodes
are used, so that a large \\\nu\\, where the range is short, is still
covered.

## Examples

``` r
q <- statmodels7:::gamma_nodes(1.5)
c(length(q$w), sum(q$a), sum(q$a * q$w))
#> [1] 40.000000  1.000000  1.000009
```
