# The Refreshable Terms' Own Verdict

Reads
[`modelterms7::term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.html)
for every term that recomputes its own block and returns `TRUE` when all
of them report that their own iteration has settled. This is the verdict
for the refreshable half of a fit.

## Usage

``` r
statmod_refresh_settled(spec, design, which = "all")
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, whose refresh state holds the terms to be checked.

- which:

  Which entries to check: `"all"` (the default), `"jacobian"` or
  `"frozen"`.

## Value

A single logical. `TRUE` when every term checked reports that it has
settled, and `TRUE` when there is nothing to check, so a model with no
refreshable term never blocks a fit's verdict on this.

## The score of the working model

Where a term's block is the Jacobian of its contribution, the gradient
of the model's objective is the model's own gradient and its vanishing
is the test. Where the block is a working linearization with a frozen
weight, as in a discontinuous break-point term, the gradient belongs to
the working model alone, and the profile objective there is a step
function in the break-point with no gradient to vanish at all.

A fit of a discontinuous break-point term can report convergence at a
position where the score of the working model is very large, because the
auxiliary columns of the working block grow as the rescaling factor is
annealed, so that a stopping rule reading that score would not stop.

## See also

[`modelterms7::term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.html)
for what each construction answers,
[`statmod_design_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md)
for the refresh,
[`statmod_commit_refresh()`](https://statmodels7.github.io/statmodels7/reference/statmod_commit_refresh.md)
for the state that it reads.
