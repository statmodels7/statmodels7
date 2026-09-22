# Whether a Model Carries a Sharp Break-Point Term

`TRUE` where some equation carries
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
without `smoothed`, whose contribution has a kink in the break-point at
the observations.

## Usage

``` r
has_sharp_breakpoint(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A single logical.

## Details

[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
reads it to keep
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
from escalating the Levenberg damping on a rejected step: at the minimum
of an objective with a kink every Gauss-Newton step is rejected whatever
the damping, so the escalation spends eight attempts at every inner fit
and ends where a stop would have.

## See also

[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md),
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)
