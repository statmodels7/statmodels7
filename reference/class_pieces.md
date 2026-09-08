# Where a Covariance Class's Members Sit in the Stacked Vector

`class_pieces()` gives one entry per member with its parameter, its
columns in that parameter's coefficients and their positions in the
stacked vector; `class_index()` interleaves those positions group by
group, which is the order the class's penalty reads.

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
  positions in whichever vector the class is addressed in.

## Value

`class_pieces()` a list of lists with `param`, `term`, `cols` and
`index`; `class_index()` and `class_zeta_cols()` an integer vector of
`m * dim` positions.

## Details

A blockwise penalty reads its argument in consecutive chunks of the
prior's dimension, one per group: penalties7 reshapes the vector by row.
The class's prior is over the \\d\\ columns one group carries across
every member, so the index must list, for each group in turn, that
group's columns from each member.

Each member's own block is already group-major, so member \\k\\'s
columns for group \\i\\ are its positions \\(i-1)d_k + 1\\ to \\id_k\\.
The union across members is **not** contiguous – one member's columns
all precede the next member's in the stacked vector, and the two may be
in different equations – so what comes out is a permutation rather than
a range. Nothing downstream minds: reading and writing a matrix at
`[index, index]` is correct for any index, provided the penalty's own
output is in the same order, which is what this ordering arranges.

A class whose members are all inside a structural term is addressed in
that term's own parameters instead, which is the vector its penalty is
read at and where the design has no column. `class_zeta_cols()`
interleaves those positions by the same rule, and `class_pieces()` is
not consulted at all, there being no design block to look the columns up
in.

## See also

[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md),
their caller.
