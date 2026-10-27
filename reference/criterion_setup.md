# The Hyperparameter of a Term for the Criterion Graph

Finds the row of the outer index that belongs to the hyperparameter
drawn, and gathers what the evaluation of the criterion needs.

## Usage

``` r
criterion_setup(fit, ctx, hyper)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ctx:

  The context described at
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

- hyper:

  `NULL`, or the name of the hyperparameter.

## Value

A list with the design, the blocks, the outer index, the free
coordinates `eta0` of every hyperparameter, the selected `row`, its
`name` and its fitted `value`.

## Details

An error is signaled where the term has no hyperparameter that the outer
criterion estimates (a held value, a lasso, scad or mcp penalty, which
is chosen along a path), where the fit ran no outer criterion, where the
hyperparameter is not positive, and where the model has a structural
term.
