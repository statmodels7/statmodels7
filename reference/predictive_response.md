# The Predictive Distribution of the Response

Returns the median, the ends of an interval and the standard deviation
of the response at each row, under a mixture of the family over the
predictors.

## Usage

``` r
predictive_response(spec, comps, level)
```

## Arguments

- spec:

  The specification at the rows predicted.

- comps:

  A list of components, each a list with `eta`, a named list of the
  predictors' means, `C`, their covariance as an array `P x P x n` or
  `NULL` for none, and optionally `weight`, equal weights otherwise. If
  the list has the attribute `infinite_variance` equal to `TRUE`, the
  standard deviation is reported as `NA`.

- level:

  The interval's level.

## Value

A data frame with `fit` (the median), `se`, `lower` and `upper`.

## Details

Each component \\c\\ gives, at row \\i\\, the predictors' mean
\\\eta\_{ci}\\ and their covariance \\C\_{ci}\\, and the predictors are
taken jointly Gaussian within it. The distribution of \\Y\\ is the
mixture of the family over the components, with weights \\w_c\\ that sum
to one (equal where a component carries none), and, within each
component, over the predictors: \$\$G_i(y) = \sum\_{c} w_c E\_{\eta \sim
\mathrm{N}(\eta\_{ci}, C\_{ci})}\[F(y \mid h^{-1}(\eta))\].\$\$ Each
expectation is evaluated at the single node \\\eta\_{ci}\\ where
\\C\_{ci}\\ is zero at every row. Otherwise \\C\_{ci} = LL^\top\\ is
split by its eigenvectors, and \\r\\ directions carry variance. The
directions after the first are a Gauss-Hermite product grid with
\\k^{r-1}\\ nodes, where \\k = \max(2, \min(8, \lfloor 64^{1/(r-1)}
\rfloor))\\, and the first, the one of largest variance, is a 20-node
Gauss-Hermite rule. Where that variance is large against the spread of
the conditional law, the law's distribution function turns from 0 to 1
between two adjacent nodes, which the rule integrates badly. A direction
cell whose distribution function differs by more than 0.3 between two
adjacent nodes is therefore integrated by
[`numericals7::quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.html)
instead, at its default tolerances. A direction cell is left to
Gauss-Hermite whatever its values when it is among the lightest, whose
weights add up to at most 1e-8. One component carries the default
interval, and two sources of variation are folded into its covariance;
the parametric bootstrap gives one component per replica. A quantile of
\\G_i\\ lies between the smallest and the largest quantile of the
conditional laws at the nodes, because a mixture's distribution function
is an average of theirs. For a continuous family it is found by Newton's
method on \\G_i(y) = p\\ with the mixture's density as the derivative,
started at the weighted median of those quantiles and kept inside a
bracket that every evaluation shrinks; for a discrete family it is found
by bisection on the integers, as the smallest value at which \\G_i\\
reaches the level. The quantiles read the family's distribution function
and quantile function, and the density for a continuous family, so the
family does not need a location parameter. The standard deviation reads
the family's mean and variance, and is
\\\sqrt{E\[\mathrm{Var}(Y\mid\eta)\] + \mathrm{Var}(E\[Y\mid\eta\])}\\,
`NA` where the family's mean or variance does not exist.
