# The Penalty of a Specification at Given Coefficients

Computes the sum of every penalized term's `penalty_value()`. The
gradient and Hessian of that sum are placed in the columns that each
term owns.

## Usage

``` r
statmod_penalty_at(
  spec,
  coef,
  hyper,
  design = statmod_design(spec),
  what = c("value", "gradient", "hessian"),
  jeffreys = TRUE
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

- jeffreys:

  Whether the Jeffreys penalties of
  [`jfp_penalty()`](https://statmodels7.github.io/statmodels7/reference/jfp_penalty.md)
  are included (the default).
  [`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  passes `FALSE` and places their curvature itself, split into the part
  with the sparsity of the information and a product of low rank.

## Value

Depends on `what`:

- `"value"`:

  a single number, the sum of the penalties. The value is `0` for an
  unpenalized model, and `Inf` where the information of a Jeffreys block
  is not positive definite.

- `"gradient"`:

  a named list of numeric vectors, one per distribution parameter, zero
  in every unpenalized coordinate.

- `"hessian"`:

  a symmetric `p x p` matrix over the stacked coefficients, zero outside
  the penalized blocks: a Matrix object when any equation's design is
  sparse, a base matrix otherwise.

## Details

For each penalized term the function evaluates the quantity at the
coefficients and the hyperparameters of that term and places the result
in the columns that the term owns. A term declaring several penalties,
each over a subset of its coefficients, contributes each one at its own
coordinates. A covariance class that spans several equations is gathered
from the stacked vector, and its gradient is scattered back to the
equations that own the coordinates.

The accumulator for the Hessian follows the design, by the rule of
[`zero_information()`](https://statmodels7.github.io/statmodels7/reference/design_sparse.md):
it is a sparse matrix where any equation's design is sparse and a base
matrix otherwise. The Hessian of a penalty is added to a sparse
accumulator as it is, and it is converted to a base matrix when the
accumulator is a base matrix.

A penalty over the parameters of a structural term enters the value
only, and its derivatives are returned by
[`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md).
A covariance class whose coordinates are partly coefficients and partly
parameters of a structural term enters the value only, and its
derivatives are placed in the joint vector by
[`joint_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/joint_penalty_at.md).

## See also

[`statmod_loglik_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_loglik_at.md)
for the other half of the objective,
[`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
for the derivatives of penalties over a structural term's own
parameters, whose value this function includes.
