# Which Vector a Covariance Class Is Addressed In

Which of the three vectors a class's coordinates are positions in.

## Usage

``` r
class_space(cl)
```

## Arguments

- cl:

  One class, as
  [`statmod_classes()`](https://statmodels7.github.io/statmodels7/reference/statmod_classes.md)
  assembles it.

## Value

One of `"coef"`, `"zeta"` or `"mixed"`.

## Details

A member written in an equation, or in the subformula of an additive
term, has columns in the stacked coefficient vector. A member inside a
structural term has none: its coordinates are that term's own
parameters, which the design carries in its structural state.

A class whose members are all of one kind is read in that vector. A
class **split between them** is read in the JOINT vector the fit already
builds, \\\[\beta; \zeta\_{\mathrm{free}}\]\\: the inner step, the
marginal criterion and the variance all assemble it, with the
coefficients first and a filter's free parameters after them. What such
a class adds to those three is the CROSS block, which the penalty
produces on its own – it is one Hessian over the class's stacked vector
and knows nothing of the split.

## See also

[`statmod_classes()`](https://statmodels7.github.io/statmodels7/reference/statmod_classes.md),
its only caller;
[`class_joint_pieces()`](https://statmodels7.github.io/statmodels7/reference/class_joint_pieces.md)
for the positions a mixed class is read at.
