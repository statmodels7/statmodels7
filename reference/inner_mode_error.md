# How Far Above Its Mode the Inner Fit Stopped

\\\tfrac12 g'K^{-1}g\\, the decrease the penalized likelihood's own
Newton correction predicts at the point the inner fit returned: how much
log-likelihood is still on the table there.

## Usage

``` r
inner_mode_error(ctx, spec, design, coef, hyper, score, expected = FALSE)
```

## Arguments

- ctx:

  The evaluation context, so the penalized factorization is the one the
  criterion will read rather than a second copy.

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients the inner fit returned.

- hyper:

  The hyperparameters.

- score:

  The inner objective's gradient at those coefficients.

- expected:

  Whether the penalized information is the expected one.

## Value

A single number, or `NA` where the penalized system could not be read
there – which is itself a reason to call the point unavailable.

## Details

It answers the question **availability** asks, which is a different
question from the one the inner optimizer's flag answers. The flag says
whether a stopping rule fired. Availability asks whether the criterion,
a Laplace expansion at the mode, is valid at this point.

The second is a matter of distance and has a natural scale,
log-likelihood units. The first is a boolean about a threshold on a
score whose size depends on the model.

The two come apart on nearly every point. Measured on
`y ~ s(x) | sigma ~ s(z)`, of 38 inner fits during one search, 38 are at
their mode by this reading, between 1e-09 and 3e-09 against a limit of
1e-03, and **four** report convergence. The other 34 stopped on the
objective-stall guard with the objective already fixed to twelve
significant digits and a score oscillating between 2.5e-06 and 3.3e-06,
just above the absolute tolerance of 1e-06. Read as unavailable, they
made the outer line search backtrack eleven times per iteration and
accept a step of 0.0026 where the Newton step is 1.4, so the search
moved 0.005 in eta over 38 evaluations and stopped 4.0 criterion units
below the optimum its own gradient was correctly pointing at.

It is used to add points and never to remove one: a run whose flag says
converged stays usable whatever this reads, so no model that fitted
before can stop fitting. That is also why it is not folded into the flag
itself, which `piano_stabilita.txt` section 13 measured and withdrew –
there the flag was made stricter, and it cost a false negative on a good
fit.

## See also

[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
