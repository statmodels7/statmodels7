# The Joint Variance of a Convention, in the Order of the Joint Information

Returns the matrix that
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
reports for `vcov_type` with `readable = FALSE`, put back in the order
of the joint information (coefficients, then the free parameters of the
structural terms).

## Usage

``` r
structural_convention_block(fit, spec, design, vcov_type, n, held)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec, design:

  Its specification and design.

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object.

- n:

  The dimension of the joint information.

- held:

  The held break-point slots, as
  [`term_held_stack()`](https://statmodels7.github.io/statmodels7/reference/term_held_stack.md)
  returns them.

## Value

A square matrix of dimension `n`.

## Details

[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
orders its output by equation; the permutation is undone here. A
coordinate without a variance is missing in its row and column, except
the slot of a held break-point, which is zero as in the inverse
information.
