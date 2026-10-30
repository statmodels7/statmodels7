# The Specification With the Labelled Quantities Held

The Specification With the Labelled Quantities Held

## Usage

``` r
profile_spec(ctx, val, spec = ctx$spec)
```

## Arguments

- ctx:

  The context of the call.

- val:

  The held values, one for each label of `which`.

- spec:

  The specification to hold them on.

## Value

`spec`, with the coefficients in `held_coef`, the hyperparameters in
`held_hyper` and the sharp break-point terms moved and held.
