# Evaluate a Formula's Terms With modelterms7 in Front

Returns an environment in which a term call resolves to modelterms7's
function whatever the user has attached, with everything else visible
behind it.

## Usage

``` r
terms_first(env)
```

## Arguments

- env:

  The environment the formula carried.

## Value

An environment.

## Details

mgcv exports `s()` and `te()` and segmented exports `seg()`. With either
attached, a user writing a statmod formula gets the other package's
function, and the failure surfaces inside `model.matrix` naming neither
the call nor the mask. Interpreting the formula in an environment whose
parent chain reaches modelterms7 first removes the ambiguity rather than
reporting it; a user who wants the other package's term writes
`mgcv::s(x)`, which the interpreter then rejects with the message it
already has.
