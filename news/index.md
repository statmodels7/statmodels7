# Changelog

## statmodels7 0.1.0

- First release. The package installs and attaches the five members of
  the toolkit – `linkfunctions7`, `distributions7`, `optimizers7`,
  `basis7` and `parameters7` – so that installing it installs all of
  them and
  [`library(statmodels7)`](https://statmodels7.github.io/statmodels7/)
  attaches all of them.

- [`statmodels7_packages()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_packages.md)
  reads the members from this package’s own `Imports`, keeping the names
  that end in `7`, so a member added there is a member everywhere and
  the list cannot disagree with itself.

- [`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
  reports the installed version of each member,
  [`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
  the exports that mask one another, and
  [`statmodels7_update()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_update.md)
  describes what is installed or installs it. The toolkit is not on
  CRAN, so the install path is GitHub through `pak`, which resolves the
  dependencies among the members itself.
