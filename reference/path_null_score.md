# The Largest Score a Kinked Block Has to Beat

\\\max_j \|\partial(-\ell)/\partial\beta_j\|\\ over the block's own
coefficients, with the block held at the kink.

## Usage

``` r
path_null_score(obj, beta, block, hyper)
```

## Arguments

- obj:

  The stacked objective.

- beta:

  The current coefficients.

- block:

  One entry of `statmod_blocks()$sparse`.

- hyper:

  The hyperparameters.

## Value

A single number.

## Details

A coefficient leaves the kink when the unpenalized score there exceeds
the half-width of the subdifferential, so a kink at least this wide
leaves the whole block at zero. That is where a path starts: at the
smallest hyperparameter for which the term contributes nothing, which is
glmnet's `lambda.max` written for any separable penalty.

The other coefficients are held where the caller left them rather than
refitted, so the number is a starting point and not a boundary. The path
checks it: a top whose fit is not empty is doubled until it is.
