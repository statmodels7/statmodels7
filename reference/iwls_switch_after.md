# When a Run on the Expected Information Continues on the Observed One

The iteration after which an
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
left at `hessian = "auto"`, and settled on the expected information,
continues on the observed information.

## Usage

``` r
iwls_switch_after()
```

## Value

A single number.

## Details

Fisher scoring is robust far from the mode, the expected information
being positive definite wherever the family is defined, and it converges
only linearly near the mode wherever the expected and the observed
information differ there – which they do for every model whose
dispersion has an equation of its own. Newton's step converges
quadratically near the mode. A run that has not converged after this
many scoring steps is near enough for the second to pay, and the
expected information still takes the step wherever the observed
penalized information is not positive definite or its step finds no
acceptable point.

## See also

[`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md),
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md).
