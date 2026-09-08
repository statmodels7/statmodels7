# Where a Mixed Covariance Class's Members Sit in the Joint Vector

One entry per member of a class split between the design's coefficients
and a filter's own parameters, each with its positions in the joint
vector \\\[\beta; \zeta\_{\mathrm{free}}\]\\ that the inner step, the
marginal criterion and the variance all assemble.

## Usage

``` r
class_joint_pieces(cl, design, params, offs)
```

## Arguments

- cl:

  One class, from
  [`statmod_classes()`](https://statmodels7.github.io/statmodels7/reference/statmod_classes.md),
  whose `space` is `"mixed"`.

- design:

  The design.

- params:

  The distribution's parameters, in order.

- offs:

  Where each parameter's coefficients start in the stacked vector.

## Value

A list of lists, one per member, each the piece with `joint` added and,
for a coefficient member, `cols` and `index`.

## Details

The joint order is the one those three already write and is not invented
here: the stacked coefficients first, then the free parameters of the
one structural term of the filter shape, in the order the term holds
them less whichever a linear intercept already carries.

A coefficient member's positions are its parameter's offset plus its
columns, exactly as
[`class_pieces()`](https://statmodels7.github.io/statmodels7/reference/class_pieces.md)
gives them. A structural member's are `nb` plus its place among the FREE
parameters, which is not its place among all of them: a held one is not
in the joint vector at all.

## A held coordinate cannot be in a class

The prior's dimension is fixed when the class is assembled, from the
members' widths. If one of the coordinates it collects is then held –
the level of a filter is, wherever a linear intercept in the same
equation already carries the constant – the class has fewer coordinates
than its prior describes and there is no honest matrix to estimate. It
is rejected here, where the held set is visible, naming the coordinate
and the two ways out.

## See also

[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md),
its caller;
[`joint_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/joint_penalty_at.md),
which reads the positions.
