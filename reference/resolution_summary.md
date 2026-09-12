# The Resolution the Search Is Given, Out of the Readings Taken So Far

The **largest** of the usable readings
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
has returned during this search.

## Usage

``` r
resolution_summary(x)
```

## Arguments

- x:

  Every usable reading so far.

## Value

A single positive number, or `NA_real_` where there is none.

## Details

A reading is taken at every usable evaluation, and each one asks how far
the criterion moves when the mode is displaced by the error **that**
evaluation's inner fit happened to stop at. The displacement varies
across evaluations by far more than the quantity it stands for does, so
which reading is kept decides the answer.

**The readings are already filtered to acceptable modes**, which is what
makes the largest of them the right one to keep:
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
returns `NA` where the mode error exceeds
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
so every reading that reaches here comes from a mode the layer would
accept, and what a resolution has to bound is how far the criterion
**can** move between two such modes – not how little it happened to move
at the luckiest evaluation.

**What the smallest cost**, measured directly. At one hyperparameter
reached from six different warm starts on a `pig1` smooth, the
criterion's own spread is `2.3e-06` under `iwls(hessian = "expected")`
and `4.7e-07` under `"observed"` – a factor of 5 between the two
branches. The readings taken during those searches run `1.8e-07` to
`3.5e-06` and `1.1e-10` to `1.2e-07`, so the largest of them is within a
factor of 4 of the spread in both (1.5 times high and 3.9 times low),
where the smallest is 13 times low in the first and **4300** times low
in the second. The inner score at which a Newton step stops varies by
three orders across evaluations, and a minimum over such a sequence
reads the luckiest one.

A resolution far below the truth does not stop a search early, it stops
it from ever returning: the line search goes on backtracking for
improvements the criterion cannot resolve, and the run ends by
exhausting
[`outer_backtracks()`](https://statmodels7.github.io/statmodels7/reference/outer_backtracks.md)
instead of on a rule. Measured over 60 fits of ten shapes at three seeds
each, the largest reading against the smallest: **368 criterion
evaluations against 577**, the search's own flag met in 58 against 43,
no flag lost anywhere, and
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
returning **the same verdict on every one of the 60** – the resolution
decides what a search costs and not where it arrives. Under
`hessian = "observed"` alone it is 173 evaluations and 30 flags of 30
against 334 and 18. What is given up is at most `1.1e-05` of criterion
on a criterion of order \\10^2\\ to \\10^3\\, with a median of `8.7e-07`
over the fits where anything is given up at all.

## See also

[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md),
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
