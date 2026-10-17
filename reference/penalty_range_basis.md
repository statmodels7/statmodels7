# An Orthonormal Basis of a Penalty's Range Space

Returns an orthonormal basis of the directions that a penalty
constrains.
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
integrates over these directions, and
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
integrates them together with the others.

## Usage

``` r
penalty_range_basis(pen, k, p, nm)
```

## Arguments

- pen:

  A penalties7 penalty object.

- k:

  The number of coefficients in the term's block.

- p:

  The distribution parameter, named in the error message.

- nm:

  The term's name, named in the error message.

## Value

A `k x r` matrix with orthonormal columns, where `r` is the penalty's
rank. It is the `k x k` identity for a proper penalty.

## Details

A proper penalty constrains every direction, so the basis is the
identity and
[`penalties7::is_proper()`](https://statmodels7.github.io/penalties7/reference/is_proper.html)
is the only function that is called. A penalty with a null space
supplies it through
[`penalties7::penalty_null_basis()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.html),
and the range basis is the orthogonal complement.

A penalty that supplies neither is rejected by name, because integrating
over a guessed subspace would not give the marginal likelihood.

## See also

[`integrated_basis()`](https://statmodels7.github.io/statmodels7/reference/integrated_basis.md),
which assembles these blockwise.
