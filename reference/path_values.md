# The Values a Path Visits

A geometric grid of kink sizes from the one that empties the block down
to `min_ratio` of it, carried back onto the hyperparameter.

## Usage

``` r
path_values(pen, theta, name, s_max, n_values = 40L, min_ratio = 0.001)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters in force.

- name:

  Which one the path varies.

- s_max:

  The size of the kink at the top of the path.

- n_values:

  How many points.

- min_ratio:

  The smallest kink size, as a fraction of `s_max`.

## Value

A numeric vector of values for `name`, from the emptiest fit to the
fullest.

## Details

The grid is geometric in the size of the kink rather than in the
hyperparameter, so that a penalty whose kink narrows as its
hyperparameter grows is swept in the same order as one whose kink
widens: from the empty model towards the full one. Values the penalty
cannot reach are dropped.
