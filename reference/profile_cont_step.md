# One Step of a Continuation

Predicts the coefficients at `val` from the point `prev` and corrects
them. A failed corrector halves the step through a point that is not
returned.

## Usage

``` r
profile_cont_step(ctx, method, prev, val, nested, se, depth = 0L)
```

## Arguments

- ctx:

  The context of the call.

- method:

  The
  [`profile_continuation()`](https://statmodels7.github.io/statmodels7/reference/profile_continuation.md)
  object.

- prev:

  The previous point, with the held values in `val`.

- val:

  The held values to reach.

- nested:

  Whether the profile follows the nested problem.

- se:

  The scale of each label of `which`, for the nested problem.

- depth:

  The number of halvings so far.

## Value

The point at `val`, with `val` in it.
