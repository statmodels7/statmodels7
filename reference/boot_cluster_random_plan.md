# The Draw of a Cluster Replicate Beside Random Effects

Returns the function of
[`bootstrap_plan()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_plan.md)
for a cluster bootstrap of a fit with
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
terms. A draw is the vector of the clusters drawn, and the attribute
`build` of the result turns a draw into the specification of the
replicate.

## Usage

``` r
boot_cluster_random_plan(object, resample, spec0, design0)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- resample:

  A
  [ClusterResample](https://statmodels7.github.io/statmodels7/reference/ParametricResample.md).

- spec0, design0:

  Its specification and design.

## Value

A function of no arguments returning a list with `pick`, with attributes
`kind` (`"regrouped"`), `targets`, `random` and `build`.

## Details

The clusters are drawn with replacement, \\G\\ of \\G\\. A cluster drawn
\\k\\ times enters the replicate \\k\\ times, once for each draw, and
each copy is a new cluster. The grouping variable of every random-effect
term is relabeled in the replicate so that the groups of two copies are
different groups, each with an effect of its own, as in Stata's
`bootstrap, cluster() idcluster()`. A group is relabeled when every one
of its observations lies in a single cluster; a group that spans several
clusters keeps its label, and its effect is shared by the copies.

The design of the replicate is rebuilt on the rows drawn. The terms
other than the random-effect terms are reapplied with
[`modelterms7::term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.html),
so a basis keeps the knots of the fit and a factor its contrasts. A
random effect term is built again with
[`modelterms7::term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.html)
on the relabeled groups; its hyperparameters are those of the fit. The
weights, the offsets and the response are those of the rows drawn.

The coordinates that a replicate has in common with the fit are those of
the terms that are not random effects. `build()` returns their positions
in the replicate (`map`, with `NA` for a coordinate of the fit that has
no counterpart), the starting values of the replicate, and `expand()`,
which returns the coefficients of a refit in the stacked positions of
the fit with `NA` at the random effects. A draw that does not reproduce
the design of the fit when every cluster is drawn once is an error: it
means that a variable of the model is not a column of the data.
