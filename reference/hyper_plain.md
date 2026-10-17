# The Hyperparameters as Plain Named Vectors

The Hyperparameters as Plain Named Vectors

## Usage

``` r
hyper_plain(hyper)
```

## Arguments

- hyper:

  The nested list of hyperparameters, one list per distribution
  parameter and one entry per penalized term.

## Value

The same list, in which each entry is converted to a named numeric
vector and the distribution parameters that carry no penalized term are
dropped.
