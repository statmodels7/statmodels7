# The Step and the Tolerance of the Outer Hessian's Stencil

The step
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md)
differences at, and how far its two readings may disagree before it
refuses.

## Usage

``` r
hess_stencil_step()

hess_stencil_tol()
```

## Value

A single number.

## Details

Both are measured rather than taken from a library rule.
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html)
would give \\\epsilon^{1/3}\\, about 6e-6, which is right for
differencing a function evaluated to machine precision and wrong here:
the gradient is computed by refitting a mode, so what bounds the step
below is that reproducibility and not the rounding of a double.

Swept on a penalized filter against a second difference of the
criterion, the result is 8.0e-05 out at \\h = 0.1\\, where the
truncation still shows, and then FLAT at 1.06e-05 – the reference's own
error – for every \\h\\ from 1e-2 down to 3e-5. `1e-3` is two decades
below where truncation matters and two above where anything else begins.

The tolerance separates two regimes that are SIX ORDERS apart, with
nothing between them: 5.5e-08, 9.3e-07 and 1.1e-06 where the mode is
well located, against 6.0e-01 and 5.7e-01 where it is not. `1e-3` sits
three orders above the worst resolved reading and three below the best
unresolved one, and a relative error of that size in the curvature is
5e-4 in a standard error, under the fourth significant figure a summary
prints.

## See also

[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md)
