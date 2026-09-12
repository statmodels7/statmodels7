# The Rise a Criterion Would Still Buy at a Point

The Newton decrement \\\tfrac{1}{2}g^\top A^{-1} g\\, in the criterion's
own units.

## Usage

``` r
joint_decrement(g, A)
```

## Arguments

- g:

  The outer criterion's gradient, over the coordinates under test.

- A:

  The negated symmetric curvature over the same coordinates.

## Value

A single number, or `NA_real_`.

## Details

Under a quadratic model of the criterion at the reported point, taking
the step \\A^{-1}g\\ raises it by exactly this, so the decrement answers
"how much is left here" in the units a reader compares criteria in. That
is what a gradient cannot do: measured over 1350 fits of eight shapes at
five sample sizes, against an independently located optimum, the 838 of
them whose gap exceeds `1e-6` give \$\$\log\_{10}(\mathrm{decrement}) =
0.0161 + 1.0018 \log\_{10}(\mathrm{gap}),\qquad R^2 = 0.9980\$\$ over a
gap running from `1.3e-06` to 137 – **eight** orders of magnitude –
where the gradient at those same points spans 3.9e-04 to 0.81 on healthy
fits of one shape alone as \\n\\ runs from 300 to 30000.

As a reading it separates the two classes completely where the gradient
does not. On 663 fits within `1e-3` of their optimum and 552 stopped
more than `1e-2` short of it:

|                       |                           |                       |
|-----------------------|---------------------------|-----------------------|
| reading               | flagged among the healthy | found among the short |
| \`max                 | g                         | \> 1e-2\`             |
| \`max                 | g                         | \> 1e-2 \* sqrt(n)\`  |
| \`max                 | g                         | \> 1e-2 \* n\`        |
| **decrement \> 1e-2** | **0 of 663**              | **552 of 552**        |

⚠️ The premise that the gradient grows with \\n\\, which the two middle
rows were written for, does not reproduce: at a fixed model structure
its slope in \\\log\_{10} n\\ is -0.42, -0.46 and +0.06 on the three
shapes that have one, and the control is the pair of shapes with the
same formula whose group count is fixed at 40 and proportional to \\n\\
– the first falls at -0.46 where the second rises at 1.24. What makes
the gradient grow is the number of penalized coefficients.

⚠️ And the perfect separation rests partly on the margin between the two
classes, healthy at or under `1e-3` against short by more than `1e-2`,
where the decrement's own spread within a class is about 1.5 orders.
What is solid is the calibration above.

`NA` where \\A\\ is not positive definite, the point then not being a
maximum in those directions and no rise following from it.

## See also

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
[`outer_curvature()`](https://statmodels7.github.io/statmodels7/reference/outer_curvature.md),
[`coord_decrement()`](https://statmodels7.github.io/statmodels7/reference/coord_decrement.md)
