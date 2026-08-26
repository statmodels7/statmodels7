# Fit the Structural Terms' Own Parameters

Estimates each structural term's parameters on the unconstrained scale
of
[`modelterms7::term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.html),
the coefficients held where they are, by a general optimizer with the
exact gradient.

## Usage

``` r
statmod_fit_structural(
  spec,
  design,
  obj,
  beta,
  hyper,
  optimizer,
  verbose = FALSE
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- obj:

  The objective, from
  [`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md).

- beta:

  The stacked coefficients, held fixed.

- hyper:

  The hyperparameters.

- optimizer:

  An optimizers7 optimizer, or `NULL` for
  [`optimizers7::lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.html).

- verbose:

  Whether to print a line per term.

## Value

A list with `value`, `converged` and `iterations`.

## Details

A general optimizer serves here and the scoring step does not, because
these are not coefficients of a design block: there is no \\X'WX\\ to
invert, the predictor being a recursion and not a product.

The gradient is exact.
[`statmod_structural_score()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_score.md)
chains the filter's own Jacobian, which is already the total derivative
of the predictor in those parameters, so a quasi-Newton method has
everything it needs.

## See also

[`statmod_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural.md)
