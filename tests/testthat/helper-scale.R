# The distance between two floating-point routes to one quantity, read on the
# scale of the expected values (floor 1e-4), as in distributions7's
# helper-distributions.R. The attributes, the pattern of non-finite entries
# and their values must match exactly; Inf is returned for any mismatch of
# structure.
scale_gap <- function(object, expected, floor = 1e-4) {
  if (length(object) != length(expected) ||
      !identical(attributes(object), attributes(expected))) return(Inf)
  fe <- is.finite(expected)
  if (!identical(is.finite(object), fe) ||
      !identical(object[!fe], expected[!fe])) return(Inf)
  if (!any(fe)) return(0)
  scale <- max(abs(expected[fe]), floor)
  max(abs(object[fe] - expected[fe])) / scale
}

expect_agrees_on_scale <- function(object, expected, tolerance = 1e-12,
                                   label = NULL) {
  lab <- if (is.null(label)) deparse(substitute(object)) else label
  expect_lte(scale_gap(object, expected), tolerance,
             label = paste(lab, "(gap on the scale of the expected values)"))
}
