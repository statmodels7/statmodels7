# The Family's Second-Derivative Components at a Point

Asks the distribution for its Hessian, or for its expected information,
on the link scale at one value of the parameters.

## Usage

``` r
statmod_family_hessian(spec, theta, expected = TRUE, approx = "opg")
```

## Arguments

- spec:

  The specification, read for the distribution, the response and the
  thread count.

- theta:

  The parameters, one value per observation.

- expected:

  `TRUE` for the expected information, `FALSE` for the observed Hessian.

- approx:

  The strategy the expected information uses when the family has no
  closed form, passed through to
  [`distributions7::distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.html).

## Value

The named list of second-derivative components on the link scale, as
[`distributions7::distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.html)
returns it.

## Details

Written once because two routes read it at the same point.
[`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
builds a square-root decomposition out of the per-observation blocks and
falls back to the assembled information when that decomposition cannot
be formed, and both of those used to ask the family themselves. Where
the expected information IS the outer product of scores the
per-observation block is rank one, so the square-root route is abandoned
at every point and the work done for it was thrown away; passing the
components down instead of recomputing them is what stops that.

## See also

[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
and
[`info_blocks()`](https://statmodels7.github.io/statmodels7/reference/info_blocks.md),
its two readers.
