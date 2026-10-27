# The Joint Objective Over the Coefficients and a Filter's Parameters

Returns, as closures, the objective, the gradient and the exact
information with which
[`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)
minimizes over the coefficients followed by the filter's free
parameters, so that
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
reads the same ones at a fitted point.

## Usage

``` r
statmod_joint_pieces(
  spec,
  design,
  obj,
  hyper,
  kinds = "filter",
  information = distributions7::observed_info()
)
```

## Arguments

- spec, design, obj, hyper:

  As
  [`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)
  takes them.

- kinds:

  The kinds of structural unit from which to take the term: `"filter"`
  for the joint step, and `"loglik"` as well where a mixture over latent
  states is read at a fitted point.

- information:

  The information that `he` returns:
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
  (the default) or
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html),
  for the expected information conditional on the past (see
  [`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)).

## Value

A list with `fn`, `gr`, `he`, `raw`, `setz` (writes the filter's free
parameters into the design's structural state), `zeta` (reads them),
`ix`, `nb`, `free`, `key` and `param` (the distribution parameter in
which the term sits).
