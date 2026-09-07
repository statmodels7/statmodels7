# How Far Above Its Mode a Restricted Fit Stopped

\\\tfrac12 g'K^{-1}g\\ over the FREE coordinates of a restricted fit:
how much of the penalized objective is still on the table at the point
[`statmod_restrict()`](https://statmodels7.github.io/statmodels7/reference/statmod_restrict.md)
returned, in log-likelihood units.

## Usage

``` r
restricted_mode_error(score, information, at)
```

## Arguments

- score:

  The restricted fit's score, all coordinates.

- information:

  Its penalized information, all coordinates.

- at:

  The position of the held coordinate.

## Value

A single number, or `NA` where the free block could not be inverted
there – which is itself a reason to doubt the point.

## Details

It is
[`inner_mode_error()`](https://statmodels7.github.io/statmodels7/reference/inner_mode_error.md)'s
question asked of a restricted fit, and it is asked for the reason
statmodels7 0.81.0 established: the inner optimizer's flag says whether
a stopping rule fired, which is a boolean about a threshold on a score
whose size depends on the model, while whether a point is at its mode is
a matter of distance and has a natural scale.

The two come apart here as they do there, and they come apart BACKWARDS:
the flag rejects the points that are located best. Measured on `y ~ x`
with a Poisson response at \\n = 200\\, holding the slope at nineteen
values from 0.6 to 1.5, all nineteen are at their mode – the worst at
1.04e-11 against a limit of 1e-03 – and fifteen report the flag. The
four it rejects are **the four best-located of the nineteen**, every one
of them near 1e-23 where the fifteen it accepts run out to 1.04e-11. On
a gaussian smooth holding `s(x).lin` at thirteen values it is the two
best-located of the thirteen, at 1.31e-12 against a worst of 6.67e-09.

The mechanism is that a stopping rule reads a CHANGE and a distance
reads a POINT: where the refit lands on its mode in one step the
objective does not move between iterations, the objective-stall guard
fires, and the run is labelled stopped rather than converged. The better
the point, the likelier the flag denies it, which is why a warning
printed off the flag would fire loudest where there is least to warn
about.

The held coordinate is excluded because it is not free: its score is the
quantity Rao's statistic reads and does not vanish under the
restriction.

## See also

[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
the limit it is read against, and
[`inner_mode_error()`](https://statmodels7.github.io/statmodels7/reference/inner_mode_error.md),
the same reading inside a search.
