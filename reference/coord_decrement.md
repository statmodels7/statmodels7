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
reads it for the boundary label, because a verdict over the whole vector
does not show which coordinate has stopped moving. It is the same
maximum that
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md)
takes, restricted to one coordinate, so it is bounded by the joint
reading, and a coordinate that it calls settled has already passed the
test that the verdict applies. This keeps `flat` a label and not a
verdict. It is a weaker reading than the joint decrement and is not used
as a verdict.

`Inf` where the diagonal or the gradient entry is not finite or the
diagonal is not positive, so that such a coordinate is never called
settled.

## See also

[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md),
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
