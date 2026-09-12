# What One Coordinate Alone Would Buy

\\g_j^2/(2A\_{jj})\\, the decrement of the one-coordinate problem.

## Usage

``` r
coord_decrement(g, A)
```

## Arguments

- g:

  The outer criterion's gradient.

- A:

  The negated symmetric curvature.

## Value

A numeric vector as long as `g`.

## Details

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
reads it for the boundary label, where a verdict over the whole vector
would say nothing about which coordinate has stopped moving. It is the
same maximum
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md)
takes, restricted to one coordinate, so it is bounded by the joint
reading and a coordinate it calls settled had already passed the test
the verdict applies – which is what keeps `edge` a label and never a
verdict.

As a verdict of its own it is the weaker reading, and that is why it is
not one: over the same 1350 fits it misses 10 of the 552 that the joint
decrement finds.

`Inf` where the diagonal is not positive, so that such a coordinate is
never called settled.

## See also

[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md),
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
