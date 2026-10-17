# Where a Covariance Class's Members Sit in the Stacked Vector

`class_pieces()` gives one entry per member with its parameter, its
columns in that parameter's coefficients and their positions in the
stacked vector; `class_index()` interleaves those positions group by
group, in the order in which the class's penalty reads them.

## Usage

``` r
class_pieces(cl, design, params, offs)

class_interleave(cl, pos)

class_index(cl, design, params, offs)

class_zeta_cols(cl)
```

## Arguments

- cl:

  One class, from
  [`statmod_classes()`](https://statmodels7.github.io/statmodels7/reference/statmod_classes.md).

- design:

  The design.

- params:

  The distribution's parameters, in order.

- offs:

  Where each parameter's coefficients start in the stacked vector.

- pos:

  One integer vector per member, each group-major, giving that member's
  positions in the vector in which the class is addressed.

## Value

`class_pieces()` a list with one element per member: the member's piece
(see
[`label_pieces()`](https://statmodels7.github.io/statmodels7/reference/label_pieces.md)),
to which `cols` and `index` are added; `class_index()` and
`class_zeta_cols()` an integer vector of `m * dim` positions.

## Details

A blockwise penalty reads its argument in consecutive chunks of the
prior's dimension, one per group: penalties7 reshapes the vector by row.
The class's prior is over the \\d\\ columns that one group carries
across every member, so the index must list, for each group in turn,
that group's columns from each member.

Each member's own block is already group-major, so member \\k\\'s
columns for group \\i\\ are its positions \\(i-1)d_k + 1\\ to \\id_k\\.
The union across members is not contiguous (one member's columns all
precede the next member's in the stacked vector, and the two may be in
different equations), so the index is a permutation and not a range.
Reading and writing a matrix at `[index, index]` is correct for any
index, provided the penalty's own output is in the same order, and this
ordering arranges that.

A class whose members are all inside a structural term is addressed in
that term's own parameters instead. That vector is the one at which its
penalty is read, and the design has no column in it. `class_zeta_cols()`
interleaves those positions by the same rule, and `class_pieces()` is
not consulted at all, because there is no design block in which to look
up the columns. A class split between coefficients and a structural
term's parameters is addressed in the joint vector by
[`class_joint_pieces()`](https://statmodels7.github.io/statmodels7/reference/class_joint_pieces.md).

## See also

[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md),
their caller.
