# The Joint Objective Over the Coefficients and a Filter's Parameters

The objective, gradient and exact information
[`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)
minimizes, over the coefficients followed by the filter's free
parameters, returned as closures so
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
reads the same ones at a fitted point.

## Usage

``` r
statmod_joint_pieces(spec, design, obj, hyper)
```

## Arguments

- spec, design, obj, hyper:

  As
  [`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)
  takes them.

## Value

A list with `fn`, `gr`, `he`, `raw`, `setz` (writes the filter's free
parameters into the design's structural state), `zeta` (reads them),
`ix`, `nb`, `free` and `key`.
