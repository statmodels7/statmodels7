# The Hyperparameters That an Outer Method Estimates

Returns one row per hyperparameter of a penalized term that is fitted in
the joint system, with the link that carries it onto the whole line.

## Usage

``` r
outer_hyper_index(spec, blocks)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- blocks:

  The block split, as
  [`statmod_blocks()`](https://statmodels7.github.io/statmodels7/reference/statmod_blocks.md)
  returns it.

## Value

A data frame with one row per estimated hyperparameter, or per group of
shared ones, and three character columns, `parameter`, `term` and
`name`. The links that carry each onto the whole line are in the
`"links"` attribute, a list with one element per row of the frame (a
list column would not survive the subsetting that the search does). The
`"members"` attribute is a data frame of the same three columns plus
`row`, which lists the hyperparameters that each row represents.

## Details

A term whose penalty has a kink is left out. Its coefficients are
estimated by a coordinate descent while everything else is held fixed,
and the criterion is a Laplace approximation, which needs a second
derivative that does not exist at a kink. Those hyperparameters are
chosen by a path instead.

A hyperparameter that its own term holds is left out too, whatever the
kind of penalty to which it belongs, and so is every member of a shared
label when one member holds it.

The result is the index of the hyperparameters on the outer search. The
search also moves the coefficients that the criterion estimates (see
[`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)),
so its dimension is the number of rows plus the number of those
coefficients. If there are no rows and no such coefficients, the
criterion does not run.

## Shared hyperparameters

Hyperparameters carrying the same label through
[`modelterms7::term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.html)
are estimated at one value, so they occupy one row: the search has one
coordinate for the group. The row is identified by its first member, so
that every lookup that reads `parameter`, `term` and `name` keeps
working, and the `"members"` attribute lists the hyperparameters that
the row represents.

Sharing does not merge the penalties, which remain two objects estimated
at one value, so the derivative of the row is the sum of the derivatives
of its members and the value is written back under the key of each
member. If no hyperparameter is shared, the member table is the index
itself, one line per row.

Two members must carry the same link, since one free value has to lie in
both of their domains, and a label may not span a smooth penalty and a
kinked one, which are estimated by two different procedures and would be
given two different values. Both cases are rejected here, where both
sides are visible.

## See also

[`hyper_to_eta()`](https://statmodels7.github.io/statmodels7/reference/hyper_to_eta.md)
for the map that those links define,
[`statmod_held_ids()`](https://statmodels7.github.io/statmodels7/reference/statmod_held_ids.md)
for the shared labels that a term holds,
[`index_members()`](https://statmodels7.github.io/statmodels7/reference/index_members.md)
for the member table.
