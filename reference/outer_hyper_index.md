# The Hyperparameters an Outer Method Estimates

One row per hyperparameter of a penalized term that is fitted in the
joint system, with the link that carries it onto the whole line.

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
`name`. The links that carry each onto the whole line ride on the
`"links"` attribute, a list the same length as the frame has rows, since
a list column would not survive the subsetting the search does; the
`"members"` attribute is a data frame of the same three columns plus
`row`, saying which hyperparameters each row stands for.

## Details

A term whose penalty has a kink is left out. Its coefficients are
estimated by a coordinate descent with everything else held fixed, and
the criterion is a Laplace approximation, which asks for a second
derivative that does not exist at a kink. Those hyperparameters go to a
path instead.

A hyperparameter its own term holds is left out too, whatever kind of
penalty it belongs to:
[`statmod_held()`](https://statmodels7.github.io/statmodels7/reference/statmod_held.md)
is what says so.

The result is the index the whole outer search runs on. Its row count is
the dimension of that search, and at zero rows the criterion does not
run at all.

## Shared hyperparameters

Hyperparameters carrying the same label through
[`modelterms7::term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.html)
are estimated at one value, so they occupy ONE row: the search has one
coordinate for the group. The row is identified by its first member,
which is what keeps every lookup that reads `parameter`/`term`/`name`
working, and the `"members"` attribute says which hyperparameters that
row stands for.

Sharing does not merge the penalties, which stay two objects estimated
at one value, so the row's derivative is the SUM of its members' and the
value is written back under each member's own key. Where nothing is
shared the member table is the index itself, one line per row, and every
loop written over members is the loop that was there before.

Two members must carry the same link, since one free value has to land
in both their domains, and a label may not span a smooth penalty and a
kinked one, which are estimated by two different machines and would be
given two different values. Both are rejected here, where both sides are
visible.

## See also

[`hyper_to_eta()`](https://statmodels7.github.io/statmodels7/reference/hyper_to_eta.md)
for the map those links define,
[`statmod_held()`](https://statmodels7.github.io/statmodels7/reference/statmod_held.md)
for what is excluded,
[`index_members()`](https://statmodels7.github.io/statmodels7/reference/index_members.md)
for the member table.
