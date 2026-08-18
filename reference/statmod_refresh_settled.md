# Have the Refreshable Terms Settled?

`TRUE` when every term that recomputes its own block reports that its
own iteration has nothing further to say.

## Usage

``` r
statmod_refresh_settled(spec, design, which = "all")
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- which:

  Which refresh entries to ask: `"all"`, `"jacobian"` or `"frozen"`.

## Value

A single logical; `TRUE` when there is nothing to ask.

## Details

The question cannot always be answered by the score. Where a term's
block is the Jacobian of its contribution, the gradient of the model's
objective is the model's and its vanishing is the test; where the block
is a working linearization with a frozen weight, as in a discontinuous
break-point term, it is not, and the profile objective there is a step
function in the break-point with no gradient to vanish. Measured on
`jump()`: the fit reaches the break-point and the jump size to three
figures, the objective stops moving at the twelfth digit, and the score
of the working model stays at 0.176 forever.
[`term_converged`](https://statmodels7.github.io/modelterms7/reference/term_converged.html)
is what each construction answers instead.

## See also

[`statmod_design_at`](https://statmodels7.github.io/statmodels7/reference/statmod_design_at.md)
