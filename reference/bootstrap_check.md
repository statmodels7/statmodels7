# Reject a Fit That the Bootstrap Does Not Handle

Signals an error for a model that the bootstrap of
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md)
does not cover.

## Usage

``` r
bootstrap_check(object, resample, refit, spec0, design0)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- resample:

  The
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md).

- refit:

  `"coefficients"` or `"full"`.

- spec0, design0:

  Its specification and design.

## Value

`NULL`, invisibly.

## Details

A response that is not one number per observation is rejected. So are
`refit = "full"` with a
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
criterion, which needs the folds of the data, the case scheme beside a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term, and the cluster scheme beside a random-effect term that
[`boot_cluster_random_check()`](https://statmodels7.github.io/statmodels7/reference/boot_cluster_random_check.md)
cannot regroup. A model with a structural term
([`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html),
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html))
takes
[`parametric_resample()`](https://statmodels7.github.io/statmodels7/reference/parametric_resample.md)
only; a model with a block that moves with its coefficients
([`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html),
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html),
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html))
takes any scheme and is refitted by
[`boot_refit_breakpoint()`](https://statmodels7.github.io/statmodels7/reference/boot_refit_breakpoint.md).
