# Collapse the Shared Hyperparameters of an Index

Turns the member table into the index the search runs on: one row per
group, identified by the group's first member, with the member table and
the links carried alongside.

## Usage

``` r
index_group(mem, links, labs, kink_labs = character(0))
```

## Arguments

- mem:

  A data frame of `parameter`, `term` and `name`, one line per estimated
  hyperparameter.

- links:

  The links, one per line of `mem`.

- labs:

  The sharing labels, one per line of `mem`, `NA` where the
  hyperparameter is not shared.

- kink_labs:

  The labels carried by hyperparameters a path estimates.

## Value

The index, with its `"links"` and `"members"` attributes.

## Details

A label ties its members into one coordinate. Two checks run here
because this is where both sides are visible: every member of a group
must carry the same link, one free value having to land in each of their
domains, and no label may span a smooth penalty and a kinked one, which
are estimated by two different machines and would be given two values
under one name.

Where nothing is shared each row is its own group and the member table
is the index, so the result is what this function's caller produced
before sharing existed.
