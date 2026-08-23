# The Predictor and Parameters of a Simulation

The linear predictors, the distribution's parameters, and – where a term
has state – the response it drew and the latent quantity behind it.

## Usage

``` r
rstatmod_eta(spec, design, coef, structural = NULL)
```

## Arguments

- spec:

  The specification.

- design:

  Its design.

- coef:

  The coefficients.

- structural:

  The structural term's own parameters, or `NULL`.

## Value

A list with `eta`, `theta`, `y`, `latent` and `structural`.

## Details

With no structural term this is
[`statmod_eta`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
exactly, which is the point: the simulated data comes from the assembly
a fit reads, so a term whose block moves with its coefficients
contributes what it contributes rather than a linearization.

With one, the static part is assembled the same way and the term is
asked to finish through
[`term_simulate`](https://statmodels7.github.io/modelterms7/reference/term_simulate.html).
A term that draws the response as it goes returns it; one that does not
returns `NULL` there and the caller draws at the predictor.

## See also

[`rstatmod`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md),
[`statmod_eta`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
