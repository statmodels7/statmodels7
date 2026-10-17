# The Penalty Over a Structural Term's Free Parameters, as a Block

Returns the Hessian of whatever penalty covers a structural term's own
parameters, over the free ones and in their order, which is the order in
which the tail of
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
carries them.

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

  The number of free parameters that the caller expects, or `NULL` to
  accept whatever number the term has.

## Value

A square matrix, or `NULL`.

## Details

It is written once because its three readers (the variance matrix of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
the marginal criterion and
[`statmod_structural_table()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_table.md))
would otherwise each place the block for themselves. The result is
`NULL` where the structural terms carry no penalty, where the only such
penalty belongs to a mixed covariance class (placed by
[`joint_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/joint_penalty_at.md)),
or where every term has a number of free parameters other than `nfree`,
and the caller then pads with zeros.

## See also

[`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
