# The Jeffreys Penalty Over a Filter's Own Parameters

Returns the part of the gradient or of the Hessian of
[`jfp_penalty()`](https://statmodels7.github.io/statmodels7/reference/jfp_penalty.md)
that the objective over the coefficients does not carry: the entries in
the filter's free parameters, and for the Hessian the cross block, with
the block over the coefficients set to zero.

## Usage

``` r
jfp_joint_extra(spec, design, coef, what = c("gradient", "hessian"), n)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  A named list of coefficient vectors.

- what:

  `"gradient"` or `"hessian"`.

- n:

  The length of the joint vector.

## Value

A numeric vector of length `n`, or an `n x n` matrix.

## Details

A joint assembly that places
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
on the coefficients adds this to complete the Jeffreys penalty over the
joint vector. It is zero, of the requested size, where the model carries
no Jeffreys term or no filter.
