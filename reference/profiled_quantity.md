# The Quantity That a Profile Holds

The class of one label of a profile, which records what the label is and
selects the graph that
[`autoplot.StatmodProfile()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodProfile.md)
draws for it. The five subclasses are

- `ProfiledCoefficient`:

  a coefficient without a penalty, of a linear predictor or of the
  subformula of a term;

- `ProfiledPenalizedCoefficient`:

  a coefficient under a quadratic penalty, of
  [`modelterms7::s()`](https://statmodels7.github.io/modelterms7/reference/s.html),
  [`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
  or
  [`modelterms7::ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.html);

- `ProfiledBreakpoint`:

  the position of a break-point of
  [`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
  without a subformula, with the values of its covariate in
  `covariate_values`;

- `ProfiledHyperparameter`:

  a hyperparameter that
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  or
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  estimates, on the scale of its link `link`;

- `ProfiledPathHyperparameter`:

  the hyperparameter of a kinked penalty that a path chooses, on the
  scale of its link `link`.

## Usage

``` r
profiled_quantity(
  label = character(0),
  estimate = integer(0),
  standard_error = integer(0)
)

ProfiledCoefficient(
  label = character(0),
  estimate = integer(0),
  standard_error = integer(0)
)

ProfiledPenalizedCoefficient(
  label = character(0),
  estimate = integer(0),
  standard_error = integer(0)
)

ProfiledBreakpoint(
  label = character(0),
  estimate = integer(0),
  standard_error = integer(0),
  covariate_values = integer(0),
  piecewise_constant = logical(0)
)

ProfiledHyperparameter(
  label = character(0),
  estimate = integer(0),
  standard_error = integer(0),
  link = NULL
)

ProfiledPathHyperparameter(
  label = character(0),
  estimate = integer(0),
  standard_error = integer(0),
  link = NULL,
  criterion_name = character(0),
  path_record = NULL
)
```

## Arguments

- label:

  The label in `which`.

- estimate:

  The estimate on the scale of the profile.

- standard_error:

  The standard error of the estimate on that scale, `NA` where it is not
  available.

- covariate_values:

  The values of the covariate of a break-point.

- piecewise_constant:

  Whether the profile of a break-point is constant between two
  consecutive values of the covariate, which is the case of the position
  of a sharp
  [`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
  or
  [`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html).

- link:

  The link of a hyperparameter, a
  [linkfunctions7::link](https://statmodels7.github.io/linkfunctions7/reference/link.html)
  object.

- criterion_name:

  The name of the criterion that chose the hyperparameter of a path
  (`"BIC"`, `"AIC"` or `"CV"`).

- path_record:

  The points that the path scored, a data frame with the columns `value`
  (on the scale of the link), `criterion` and `se` (the standard error
  of a cross-validated criterion, `NA` otherwise).

## Value

An object of the class.
