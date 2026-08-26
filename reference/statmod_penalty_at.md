# The Penalty of a Specification at Given Coefficients

The sum of every penalized term's `penalty_value()`, with its gradient
and Hessian placed in the columns that term owns.

## Usage

``` r
statmod_penalty_at(
  spec,
  coef,
  hyper,
  design = statmod_design(spec),
  what = c("value", "gradient", "hessian")
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- hyper:

  A named list with one entry per distribution parameter, each a named
  list of hyperparameter vectors keyed by term.

- design:

  The design.

- what:

  Which quantity: `"value"`, `"gradient"` or `"hessian"`.

## Value

Depends on `what`:

- `"value"`:

  a single number, the penalties summed. `0` for an unpenalized model.

- `"gradient"`:

  a named list of numeric vectors, one per distribution parameter, zero
  in every unpenalized coordinate.

- `"hessian"`:

  a symmetric `p x p` base matrix over the stacked coefficients, zero
  outside the penalized blocks.

## Details

Each penalized term is asked for the quantity at its own coefficients
and its own hyperparameters, and the answer is placed in the columns
that term owns. A term declaring several penalties, each over a subset
of its coefficients, contributes each one at its own coordinates.

The accumulator is a base matrix even where a penalty answers with a
sparse Hessian, and the coercion is written once at that boundary.
Making it the design's own kind was measured at 0.8x end to end against
a blast radius of twenty-two consumers.

## See also

[`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
for the other half of the objective,
[`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
for penalties over a structural term's own parameters, which this does
not cover.
