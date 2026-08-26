# Have the Refreshable Terms Settled?

Asks every term that recomputes its own block whether its own iteration
has anything further to say, through
[`modelterms7::term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.html),
and returns `TRUE` when none of them does. This is the verdict for the
refreshable half of a fit.

## Usage

``` r
statmod_refresh_settled(spec, design, which = "all")
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, whose refresh state holds the terms asked.

- which:

  Which entries to ask: `"all"` (the default), `"jacobian"` or
  `"frozen"`.

## Value

A single logical. `TRUE` when every term asked reports it has settled,
and `TRUE` when there is nothing to ask, so a model with no refreshable
term never blocks a fit's verdict on this.

## Why the score cannot answer this

Where a term's block is the Jacobian of its contribution, the gradient
of the model's objective is the model's own gradient and its vanishing
is the test. Where the block is a working linearization with a frozen
weight, as in a discontinuous break-point term, the gradient belongs to
the working model alone, and the profile objective there is a step
function in the break-point with no gradient to vanish at all.

Measured on `y ~ jump(x)` at \\n = 400\\, a Gaussian response with a
step of 2 at \\x = 6\\: the fit recovers the position at 6.004 and
reports `converged = TRUE`, while the score of the working model at that
point is \\7.8 \times 10^{7}\\. The size is the annealed rescaling
factor, whose auxiliary columns grow as the schedule tightens. A rule
reading that score would never stop.

## See also

[`modelterms7::term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.html)
for what each construction answers,
[`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md)
for the refresh,
[`statmod_commit_refresh()`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
for the state it reads.
