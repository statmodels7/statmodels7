# Reject the Expected Information Beside a Term of the Likelihood Shape

Signals an error where the model carries a structural term of the
likelihood shape (a
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term or a marginal break-point term) and returns `NULL` invisibly
otherwise.

## Usage

``` r
structural_expected_check(design, what)
```

## Arguments

- design:

  The design.

- what:

  The argument that named the expected information, as it is to be
  printed.

## Value

`NULL`, invisibly.

## Details

The likelihood of such a model is a mixture over latent components: the
states of a regime term, the positions of a marginal break-point term.
The expectation of its second derivatives over the response has no
closed form, so neither the inner fit, nor an outer criterion, nor a
variance can be computed on the expected information there.

## See also

[`structural_outer_check()`](https://statmodels7.github.io/statmodels7/reference/structural_outer_check.md),
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
