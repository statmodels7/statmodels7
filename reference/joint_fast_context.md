# The Compiled Route's Context of a Joint Filter

Returns the list that
[`modelterms7::term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.html)
reads to run a joint filter compiled, or `NULL` where the response is a
matrix or the scalar registries do not know the family or one of the
links.

## Usage

``` r
joint_fast_context(spec, G, d, E)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- G:

  The positions of the filtered parameters.

- d:

  The scaling.

- E:

  The static predictors of every parameter, one column each.

## Value

A list, or `NULL`.

## Details

[`distributions7::distrib_scalar_route()`](https://statmodels7.github.io/distributions7/reference/distrib_scalar_route.html)
gives the name under which the registry knows the family and its
constants, and
[`linkfunctions7::link_scalar_route()`](https://statmodels7.github.io/linkfunctions7/reference/link_scalar_route.html)
the name and the parameters of each link. Whether the family's joint
entries are written out is read by the compiled code, which leaves the
context inert where they are not.
