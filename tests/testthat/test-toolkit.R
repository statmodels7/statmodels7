# What a meta-package promises: that its list of members is the list its
# DESCRIPTION declares, that attaching it attaches them, and that it says so.

test_that("the members are the ones DESCRIPTION imports", {
  pkgs <- statmodels7_packages()
  expect_identical(pkgs, c("basis7", "distributions7", "linkfunctions7",
                           "optimizers7", "parameters7"))

  # The list is read rather than written out, so it must agree with the field
  # it is read from. Parsed here a second way, so that a mistake in the
  # parsing is not confirmed by the same parsing.
  path <- system.file("DESCRIPTION", package = "statmodels7")
  imports <- read.dcf(path, fields = "Imports")[[1L]]
  declared <- trimws(unlist(strsplit(imports, "[,\n]")))
  declared <- declared[nzchar(declared)]
  expect_true(all(pkgs %in% declared))

  # a base or third-party dependency is not a member
  expect_false("utils" %in% pkgs)
  expect_true("utils" %in% declared)
})

test_that("every member is installed and reports a version", {
  v <- statmodels7_versions()
  expect_identical(v$package, statmodels7_packages())
  expect_false(anyNA(v$version),
    label = paste("not installed:",
                  paste(v$package[is.na(v$version)], collapse = ", ")))
})

test_that("attaching the package attaches every member", {
  # test_check() has already attached statmodels7, so the members are on the
  # search path by the time this runs.
  for (p in statmodels7_packages()) {
    expect_true(paste0("package:", p) %in% search(), label = p)
  }
})

test_that("the attach message lines up and marks what is missing", {
  v <- data.frame(
    package = c("aaa7", "bbbbbbbb7", "cc7"),
    version = c("1.0", "0.0.0.9000", NA_character_),
    stringsAsFactors = FALSE
  )
  msg <- statmodels7:::statmodels7_attach_message(v)

  # two columns, so three entries make two rows
  expect_length(msg, 2L)

  # What the padding is for: the second column starts at the same offset on
  # every line that has one. Asserting equal line WIDTHS instead would be
  # asserting the opposite, since a line whose right-hand entry is missing
  # must not be padded out to meet it.
  second <- regexpr("   [vx] ", msg)
  has_second <- second > 0
  expect_identical(length(unique(second[has_second])), 1L)
  expect_false(any(grepl("\\s$", msg)))

  # a member that is not installed is marked, not silently shown as blank
  expect_match(paste(msg, collapse = "\n"), "x cc7\\s+not installed")
  expect_match(paste(msg, collapse = "\n"), "v aaa7\\s+1\\.0")

  expect_identical(statmodels7:::statmodels7_attach_message(v[0, ]), character())
})

test_that("conflicts are reported in the order the search path resolves them", {
  conf <- statmodels7_conflicts()
  expect_true(is.list(conf))

  # Whatever the toolkit's own conflicts are, each reported name must really
  # be exported by each package named, and the first must be the one a bare
  # name reaches.
  for (nm in names(conf)) {
    who <- conf[[nm]]
    expect_gt(length(who), 1L)
    for (p in who) {
      expect_true(nm %in% getNamespaceExports(p),
                  label = sprintf("%s does not export %s", p, nm))
    }
    winner <- environmentName(environment(get(nm)))
    expect_identical(winner, who[1L], label = nm)
  }

  # and the rendering names both sides
  if (length(conf)) {
    line <- statmodels7:::format_conflicts(conf)[1L]
    expect_match(line, sprintf("%s::%s masks", conf[[1L]][1L], names(conf)[1L]))
  }
})

test_that("the report describes what is installed without installing it", {
  out <- utils::capture.output(res <- statmodels7_update())
  expect_identical(res$package, statmodels7_packages())
  txt <- paste(out, collapse = "\n")
  for (p in statmodels7_packages()) expect_match(txt, p, fixed = TRUE)
  expect_match(txt, "pak::pak", fixed = TRUE)

  expect_error(statmodels7_update("nonsense"))
})
