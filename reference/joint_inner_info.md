# The Information of the Joint Step Beside a Filter

Returns the information that
[`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)
uses for a model with a score-driven filter: the
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
object bound to `method` where the caller named it as
`inner_information`, and
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
otherwise.

## Usage

``` r
joint_inner_info(method)
```

## Arguments

- method:

  The inner method, bound by
  [`inner_bind()`](https://statmodels7.github.io/statmodels7/reference/inner_bind.md).

## Value

An object that inherits from
[distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html).

## Details

The joint step runs over the coefficients and the filter's own
parameters, where the expected information is the one conditional on the
past (see
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)).
Fisher scoring on that matrix converges linearly near the mode, so an
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
that the family settled on the expected information leaves the joint
step on the observed one, and only an
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
that the caller named is used. The settlement is read from the method:
an
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
settled from
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
on the expected information switches to the observed one, and an
optimizers7 optimizer settled from it carries the attribute
`"inner_auto"`.

## See also

[`inner_info()`](https://statmodels7.github.io/statmodels7/reference/inner_info.md),
[`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)
