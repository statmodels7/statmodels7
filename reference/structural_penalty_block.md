# The Penalty Over a Structural Term's Free Parameters, as a Block

The Hessian of whatever penalty covers a structural term's own
parameters, over the free ones and in their order, which is the order
the tail of
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
carries them in.

## Usage

``` r
structural_penalty_block(spec, design, hyper, nfree = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- hyper:

  The hyperparameters.

- nfree:

  The number of free parameters the caller expects, or `NULL` to take
  whatever the term has.

## Value

A square matrix, or `NULL`.

## Details

Written once because three readers need it and a block each of them
placed for itself would agree only by accident: the joint fit, the
variance matrix, and the marginal criterion. Where no structural term
carries a penalty the answer is `NULL`, and a caller pads with zeros as
before.

## See also

[`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
