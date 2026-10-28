# The Scores of the Increments of a Regime Likelihood

Returns the derivative of the increment \\\log p(y_t \mid y_1, \dots,
y\_{t-1})\\ of the scaled forward recursion of a
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term, in the unknowns of
[`joint_scores()`](https://statmodels7.github.io/statmodels7/reference/joint_scores.md).

## Usage

``` r
regime_scores(spec, coef, design, ev)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design.

- ev:

  The evaluated predictors, from
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md).

## Value

A matrix as
[`joint_scores()`](https://statmodels7.github.io/statmodels7/reference/joint_scores.md)
returns it.

## Details

The recursion is that of
[`modelterms7::term_scores()`](https://statmodels7.github.io/modelterms7/reference/term_scores.html),
which receives the log-density and its score in the predictors of the
distribution, the derivative of each predictor in the unknowns, and the
columns of the term's own parameters. The groups run in the order of the
term's `build_info`. The sum of the rows, weighted, is the `gradient` of
[`modelterms7::term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.html).
The columns of a level that the term holds are dropped.
