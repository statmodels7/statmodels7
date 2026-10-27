# Settle the Inner Information and Bind It to the Method

Settles `information` against `distrib` and stores the result in the
inner method, so that every function that receives the method reads the
same matrix. For
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
the result is held in the properties `hessian`, `backup`, `switching`
and `repair_floor` (see
[`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)).
For an optimizers7 optimizer it is held in the attribute
`"inner_information"`.

## Usage

``` r
inner_bind(method, distrib, information)
```

## Arguments

- method:

  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

- distrib:

  The distribution to which the model is fitted.

- information:

  An object that inherits from
  [distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html).

## Value

`method`, with the settled information bound to it.

## Details

A named
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
or
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
is stored as it is. An
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
is settled by the rule of
[`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md):
the expected information where the family writes it exactly and at an
acceptable cost, and the observed information otherwise. An optimizer
has no fallback and no switch, since both belong to the loop of
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
so the settled object is the whole of the choice there.

## See also

[`inner_info()`](https://statmodels7.github.io/statmodels7/reference/inner_info.md),
[`inner_settings()`](https://statmodels7.github.io/statmodels7/reference/inner_settings.md),
[`iwls_resolve()`](https://statmodels7.github.io/statmodels7/reference/iwls_resolve.md)
