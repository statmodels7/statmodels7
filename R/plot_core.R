#' Graphical Types for a Fitted Model
#'
#' @description
#' `plot_type` is the abstract S7 class of the objects that select a graph of
#' a [StatmodFit()]. A graph is drawn by [autoplot.StatmodFit()] (which returns
#' it) or [plot.StatmodFit()] (which prints it), and the object passed as
#' `plot_type` decides which graph is drawn and with which options. Three
#' abstract subclasses divide the types by what they read:
#'
#' * `term_plot`: one term of one equation, named by the argument `term`
#'   ([plot_effect()], [plot_derivative()], [plot_shrinkage()] and others);
#' * `model_plot`: the fitted distribution as a function of the covariates,
#'   without a term ([plot_param()], [plot_quantiles()], [plot_density()],
#'   [plot_moments()]);
#' * `diagnostic_plot`: the fit against the data, through the quantile
#'   residuals ([plot_qq()], [plot_worm()], [plot_calibration()] and others).
#'
#' @details
#' A graph is produced in two steps, each an S7 generic. [plot_data()]
#' computes the quantities to draw and returns them as a list of data frames;
#' it dispatches on the type and on a second object, the term for a
#' `term_plot` and the family for the other two. [plot_layers()] turns that
#' list into ggplot2 layers and dispatches on the type alone. A new type is a
#' subclass of one of the three classes, with a method for each generic.
#'
#' @param label A single string, the name of the graph.
#' @param legend Whether the graph draws the legend that names its elements.
#'
#' @return Nothing: the classes are abstract and their constructors signal an
#'   error. Every object that inherits from `plot_type` has the properties
#'   `label` and `legend`.
#'
#' @seealso [autoplot.StatmodFit()], [plot_data()], [plot_layers()],
#'   [default_plot()]
#'
#' @examples
#' S7::S7_inherits(plot_effect(), plot_type)
#' S7::S7_inherits(plot_effect(), term_plot)
#' try(plot_type(label = "mine"))
#' @name plot_type
#' @aliases term_plot model_plot diagnostic_plot
#' @importFrom ggplot2 .data
#' @export
plot_type <- S7::new_class(
  "plot_type",
  properties = list(label = S7::class_character,
                    legend = S7::new_property(S7::class_logical,
                                              default = TRUE)),
  abstract = TRUE
)

#' @rdname plot_type
#' @export
term_plot <- S7::new_class("term_plot", parent = plot_type, abstract = TRUE)

#' @rdname plot_type
#' @export
model_plot <- S7::new_class("model_plot", parent = plot_type, abstract = TRUE)

#' @rdname plot_type
#' @export
diagnostic_plot <- S7::new_class("diagnostic_plot", parent = plot_type,
                                 abstract = TRUE)

#' The Class of the Graph Types
#'
#' @description
#' Returns the S7 class `plot_type`. The argument `plot_type` of
#' [autoplot.StatmodFit()] hides the class inside that function, which calls
#' this function to reach it.
#'
#' @return The S7 class `plot_type`.
#'
#' @keywords internal
plot_type_class <- function() plot_type


#' Check the Legend Argument of a Graph
#'
#' @description
#' Signals an error unless `legend` is `TRUE` or `FALSE`.
#'
#' @param legend The value to check.
#'
#' @return `legend`.
#'
#' @keywords internal
check_legend <- function(legend) {
  if (!is.logical(legend) || length(legend) != 1L || is.na(legend)) {
    stop("'legend' must be TRUE or FALSE.", call. = FALSE)
  }
  legend
}


#' The Default Graph of a Term
#'
#' @description
#' Returns the [plot_type] object that [autoplot.StatmodFit()] uses for a term
#' when `plot_type` is not given. The method for [modelterms7::model_term]
#' returns
#' [plot_effect()], and a class of term with a more informative graph has a
#' method of its own.
#'
#' @param term A model term, an object inheriting from
#'   [modelterms7::model_term].
#' @param ... Unused.
#'
#' @return An object inheriting from [plot_type].
#'
#' @seealso [plot_type], [autoplot.StatmodFit()]
#'
#' @examples
#' default_plot(modelterms7::s(x))
#' @export
default_plot <- S7::new_generic("default_plot", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(default_plot, modelterms7::model_term) <- function(term, ...) {
  plot_effect()
}


#' The Data of a Graph
#'
#' @description
#' Computes the quantities that a graph draws and returns them as a named list
#' of data frames, one for each layer. It dispatches on the type and on `x`:
#' the term for a [term_plot], the family of the fit for a [model_plot] or a
#' [diagnostic_plot].
#'
#' @details
#' A method has the signature `function(type, x, fit, ctx, ...)`. The context
#' `ctx` is a list with the elements `param` (the equation, or `NULL`),
#' `name` (the name of the term in that equation, or `NULL`), `vars` (the
#' variables whose effect is drawn), `cols` (the positions within the term's
#' coefficients that belong to those variables, or `NULL` for all of them),
#' `data` (the data frame passed to [autoplot.StatmodFit()], or `NULL`) and
#' `title` (what `term` named: the term, or the variable). The list that a
#' method returns may carry the attribute `"labels"`, a list of arguments of
#' `ggplot2::labs()` (`title`, `x`, `y` and others), which
#' [autoplot.StatmodFit()] applies to the graph.
#'
#' @param type An object inheriting from [plot_type].
#' @param x The term or the family on which the method dispatches.
#' @param fit A [StatmodFit()].
#' @param ctx The context described in the details.
#' @param ... Unused.
#'
#' @return A named list of data frames.
#'
#' @seealso [plot_layers()], [plot_type]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(100))
#' dd$y <- sin(2 * pi * dd$x) + rnorm(100, sd = 0.3)
#' fit <- statmod(y ~ s(x), distributions7::gaussian1_distrib(), dd)
#' str(plot_data(plot_effect(), fit@spec@terms$mu[["s(x)"]], fit,
#'               ctx = list(param = "mu", name = "s(x)", vars = "x",
#'                          cols = NULL, data = dd)), max.level = 1)
#' @export
plot_data <- S7::new_generic("plot_data", c("type", "x"),
  function(type, x, fit, ctx, ...) S7::S7_dispatch())

S7::method(plot_data, list(plot_type, S7::class_any)) <-
  function(type, x, fit, ctx, ...) {
    what <- if (S7::S7_inherits(x, modelterms7::model_term)) {
      sprintf("the term '%s'", ctx$name)
    } else "this family"
    stop(sprintf("the graph %s is not available for %s.", type@label, what),
         call. = FALSE)
  }


#' The Layers of a Graph
#'
#' @description
#' Turns the list returned by [plot_data()] into a list of ggplot2
#' components (layers, scales, guides), which [autoplot.StatmodFit()] adds to
#' an empty `ggplot2::ggplot()`. It dispatches on the type.
#'
#' @param type An object inheriting from [plot_type].
#' @param data The list returned by [plot_data()].
#' @param ... Unused.
#'
#' @return A list of ggplot2 components.
#'
#' @seealso [plot_data()], [plot_type]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(100))
#' dd$y <- sin(2 * pi * dd$x) + rnorm(100, sd = 0.3)
#' fit <- statmod(y ~ s(x), distributions7::gaussian1_distrib(), dd)
#' d <- plot_data(plot_effect(), fit@spec@terms$mu[["s(x)"]], fit,
#'                ctx = list(param = "mu", name = "s(x)", vars = "x",
#'                           cols = NULL, data = dd))
#' ggplot2::ggplot() + plot_layers(plot_effect(), d)
#' @export
plot_layers <- S7::new_generic("plot_layers", "type",
  function(type, data, ...) S7::S7_dispatch())


#' Graphs of a Fitted Model
#'
#' @description
#' `autoplot()` returns a ggplot2 object that shows one aspect of a
#' [StatmodFit()], and `plot()` prints that object and returns it invisibly.
#' The graph is chosen by `plot_type`, an object inheriting from [plot_type]:
#'
#' * with `term`, a graph of that term ([plot_effect()] when `plot_type` is not
#'   given, or the term's [default_plot()]);
#' * with a [model_plot] or a [diagnostic_plot] type, a graph of the whole
#'   fit, for which `term` and `param` are not given;
#' * with neither, the four diagnostic graphs of [diagnostics()] composed in
#'   one figure.
#'
#' @details
#' `term` is the name of a term as [summary.StatmodFit()] prints it (spaces
#' are ignored, so `"random(~1|id)"` names `random(~1 | id)`), or the name of
#' a variable. A variable selects every term of the equation in which it
#' appears and, within the linear part, the columns that involve it, so that
#' with `y ~ x + s(x)` the variable `x` draws the sum of the linear and the
#' smooth contribution. `param` is the equation; it may be omitted when the
#' term belongs to one equation only.
#'
#' The fit does not keep its data, so a graph that reads the covariates
#' (partial residuals, a grid over a covariate, residuals against a
#' covariate) requires the fitting data as `data`, and signals an error
#' without it.
#'
#' The data of the graph are in the returned object (`p$data` of each layer,
#' or `ggplot2::layer_data()`), and [plot_data()] returns them without the
#' graph.
#'
#' @param object,x A [StatmodFit()].
#' @param term The name of a term or of a variable, a single string, or
#'   `NULL`.
#' @param param The name of an equation (a parameter of the family), or
#'   `NULL`.
#' @param plot_type An object inheriting from [plot_type], or `NULL`.
#' @param data The data frame of the fit, or `NULL`.
#' @param ... Unused.
#'
#' @return `autoplot()` returns a ggplot2 object (a patchwork object for the
#'   composed diagnostics); `plot()` returns the same object invisibly, after
#'   printing it.
#'
#' @seealso [plot_type] for the list of graphs, [diagnostics()],
#'   [plot_data()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150), z = rnorm(150))
#' dd$y <- sin(2 * pi * dd$x) + 0.5 * dd$z + rnorm(150, sd = 0.3)
#' fit <- statmod(y ~ z + s(x), distributions7::gaussian1_distrib(), dd)
#' p <- ggplot2::autoplot(fit, term = "s(x)", data = dd)
#' p + ggplot2::labs(title = "The smooth of x")
#' ggplot2::autoplot(fit, term = "z", plot_type = plot_effect(residuals = FALSE),
#'                   data = dd)
#' @name autoplot.StatmodFit
#' @export
autoplot.StatmodFit <- function(object, term = NULL, param = NULL,
                                plot_type = NULL, data = NULL, ...) {
  check_dots("autoplot", ...)
  if (!is.null(data) && !is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!is.null(plot_type) && !S7::S7_inherits(plot_type, plot_type_class())) {
    stop("'plot_type' must be a plot_type object, such as plot_effect().",
         call. = FALSE)
  }
  if (is.null(term) && is.null(plot_type)) {
    if (!is.null(param)) {
      stop("'param' selects the equation of a term; give 'term' as well.",
           call. = FALSE)
    }
    return(diagnostics(object, data = data, compose = TRUE))
  }
  if (!is.null(plot_type) && !S7::S7_inherits(plot_type, term_plot)) {
    if (!is.null(term) || !is.null(param)) {
      stop(sprintf(paste0("the graph %s describes the whole fit and takes ",
                          "neither 'term' nor 'param'."), plot_type@label),
           call. = FALSE)
    }
    ctx <- list(param = NULL, name = NULL, vars = NULL, cols = NULL,
                data = data)
    return(plot_assemble(plot_type, object@spec@distrib, object, ctx))
  }
  if (is.null(term)) {
    stop(sprintf("the graph %s needs 'term'.", plot_type@label), call. = FALSE)
  }
  unit <- resolve_plot_term(object, term, param, data)
  if (is.null(plot_type)) {
    plot_type <- if (length(unit$parts) == 1L) {
      default_plot(unit$parts[[1L]]$term)
    } else plot_effect()
  }
  if (length(unit$parts) == 1L) {
    pt <- unit$parts[[1L]]
    ctx <- list(param = unit$param, name = pt$name, vars = unit$vars,
                cols = pt$cols, data = data, title = unit$title)
    plot_assemble(plot_type, pt$term, object, ctx)
  } else {
    ctx <- list(param = unit$param, name = term, vars = unit$vars,
                cols = NULL, data = data, title = unit$title)
    plot_assemble(plot_type, TermSum(parts = unit$parts), object, ctx)
  }
}

#' @rdname autoplot.StatmodFit
#' @export
plot.StatmodFit <- function(x, term = NULL, param = NULL, plot_type = NULL,
                            data = NULL, ...) {
  p <- autoplot.StatmodFit(x, term = term, param = param, plot_type = plot_type,
                           data = data, ...)
  print(p)
  invisible(p)
}

#' @importFrom ggplot2 autoplot
NULL

S7::method(autoplot, StatmodFit) <- autoplot.StatmodFit
S7::method(plot, StatmodFit) <- plot.StatmodFit


#' Assemble a Graph from Its Data and Its Layers
#'
#' @description
#' Calls [plot_data()] and [plot_layers()] and adds the theme, the layers
#' and the labels recorded on the data to an empty `ggplot2::ggplot()`. The
#' theme comes first, so that a list of layers may also carry
#' `ggplot2::theme()` elements. The legends are placed below the panel, and
#' removed when the property `legend` of the type is `FALSE`.
#'
#' @param type An object inheriting from [plot_type].
#' @param x The term or family on which [plot_data()] dispatches.
#' @param fit A [StatmodFit()].
#' @param ctx The context described at [plot_data()].
#'
#' @return A ggplot2 object.
#'
#' @keywords internal
plot_assemble <- function(type, x, fit, ctx) {
  d <- plot_data(type, x, fit, ctx = ctx)
  lab <- attr(d, "labels")
  # the theme first, so that a type's layers may adjust its elements
  p <- ggplot2::ggplot() + ggplot2::theme_bw() + plot_key_theme() +
    plot_layers(type, d)
  if (length(lab)) p <- p + do.call(ggplot2::labs, lab)
  if (!type@legend) p <- p + ggplot2::theme(legend.position = "none")
  p
}


#' The Legend of the Elements of a Graph
#'
#' @description
#' Returns the scales that name the elements of a graph in a legend below the
#' panel. Each layer maps a constant string, the description of the element,
#' to `colour` or `fill`, and `colour` and `fill` are named vectors whose
#' names are those strings and whose values are the colours drawn. A string
#' that no layer maps does not appear in the legend; the others appear in the
#' order given.
#'
#' @param colour,fill Named character vectors, or `NULL`.
#'
#' @return A list of ggplot2 components.
#'
#' @keywords internal
plot_key <- function(colour = NULL, fill = NULL) {
  out <- list(plot_key_theme())
  if (length(colour)) {
    out <- c(out, list(ggplot2::scale_colour_manual(
      name = NULL, values = colour, breaks = names(colour))))
  }
  if (length(fill)) {
    out <- c(out, list(ggplot2::scale_fill_manual(
      name = NULL, values = fill, breaks = names(fill))))
  }
  out
}

#' The Placement of the Legends of a Graph
#'
#' @description
#' Places the legends below the panel, one above the other, so that each
#' has the width of the graph.
#'
#' @return A `ggplot2::theme()` object.
#'
#' @keywords internal
plot_key_theme <- function() {
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                 legend.spacing.y = ggplot2::unit(2, "pt"),
                 legend.margin = ggplot2::margin(0, 0, 0, 0))
}

#' The Legend Key of an Estimate with Its Band
#'
#' @description
#' Draws the key of a line together with a light band behind it, in the
#' colour of the line, so that one legend entry names an estimate and its
#' interval. It is passed as `key_glyph` to the layer of the line, and the
#' layer of the band draws no key.
#'
#' @param data,params,size As for the `draw_key_*()` functions of ggplot2.
#'
#' @return A grob.
#'
#' @keywords internal
draw_key_band <- function(data, params, size) {
  col <- if (is.null(data$colour) || is.na(data$colour)) "#1f5fa8" else
    data$colour
  grid::grobTree(
    grid::rectGrob(gp = grid::gpar(col = NA, fill = col, alpha = 0.18)),
    ggplot2::draw_key_path(data, params, size))
}


#' The Name of an Estimate and Its Interval in a Legend
#'
#' @description
#' Writes the name of an estimate followed by its interval on a second line,
#' as in "estimated effect\\n(95% CI)", for a legend entry that draws the
#' estimate and its band together.
#'
#' @param what The name of the estimate.
#' @param level A number in (0, 1).
#' @param interval The name of the interval.
#'
#' @return A single string.
#'
#' @keywords internal
est_label <- function(what, level, interval = "CI") {
  sprintf("%s\n(%s%% %s)", what, format(100 * level), interval)
}

#' The Headings of the Panels of a Graph
#'
#' @description
#' Returns a labeller that heads each panel with the name of the variable
#' that divides the panels and the value of the panel, as in "x: (0.29,0.5]".
#'
#' @param name The name of the variable, or `NULL` for the value alone.
#'
#' @return A labeller for `ggplot2::facet_wrap()`.
#'
#' @keywords internal
strip_labeller <- function(name) {
  if (is.null(name)) return(ggplot2::label_value)
  force(name)
  ggplot2::as_labeller(function(v) paste0(name, ": ", v))
}

#' The Name of a Band in a Legend
#'
#' @description
#' Writes the level of a band as a percentage before its description, as in
#' "95% pointwise confidence band".
#'
#' @param level A number in (0, 1).
#' @param what The description of the band.
#'
#' @return A single string.
#'
#' @keywords internal
band_label <- function(level, what = "pointwise confidence band") {
  sprintf("%s%% %s", format(100 * level), what)
}


#' A Sum of Terms Drawn as One Effect
#'
#' @description
#' The class of the object on which [plot_data()] dispatches when a variable
#' names more than one term of an equation, as `x` does in `y ~ x + s(x)`.
#' Each element of `parts` is a list with the elements `name`, `term` and
#' `cols`, as in the context of [plot_data()].
#'
#' @param parts A list of the parts.
#'
#' @return A `TermSum` object.
#'
#' @keywords internal
TermSum <- S7::new_class("TermSum", properties = list(parts = S7::class_list))


#' Resolve the Term of a Graph
#'
#' @description
#' Finds the term, or the terms, that `term` names in the equation `param`,
#' and the columns of the linear part that belong to a variable.
#'
#' @param fit A [StatmodFit()].
#' @param term A single string, the name of a term or of a variable.
#' @param param The equation, or `NULL`.
#' @param data The data frame of the fit, or `NULL`; the columns of the
#'   linear part that belong to a variable are read from it.
#'
#' @return A list with the elements `param`, `vars`, `title` and `parts`; each
#'   element of `parts` is a list with `name`, `term` and `cols`.
#'
#' @keywords internal
resolve_plot_term <- function(fit, term, param, data) {
  if (!is.character(term) || length(term) != 1L || is.na(term)) {
    stop("'term' must be a single string.", call. = FALSE)
  }
  spec <- fit@spec
  params <- names(spec@terms)
  if (!is.null(param)) {
    if (!is.character(param) || length(param) != 1L || !param %in% params) {
      stop(sprintf("'param' must be one of %s.",
                   paste0("\"", params, "\"", collapse = ", ")), call. = FALSE)
    }
    params <- param
  }
  squash <- function(s) gsub("[[:space:]]", "", s)
  # a term named exactly: one part, whatever its variables
  hits <- Filter(function(p) squash(term) %in% squash(names(spec@terms[[p]])),
                 params)
  if (length(hits) > 1L) {
    stop(sprintf(paste0("the term '%s' appears in the equations %s; ",
                        "give 'param'."), term, paste(hits, collapse = ", ")),
         call. = FALSE)
  }
  if (length(hits) == 1L) {
    p <- hits
    k <- names(spec@terms[[p]])[squash(names(spec@terms[[p]])) == squash(term)]
    tm <- spec@terms[[p]][[k]]
    return(list(param = p, vars = term_variables(tm, k), title = k,
                parts = list(list(name = k, term = tm, cols = NULL))))
  }
  # a variable: every term of the equation that involves it
  found <- lapply(params, function(p) {
    tms <- spec@terms[[p]]
    parts <- list()
    for (k in names(tms)) {
      tm <- tms[[k]]
      if (S7::S7_inherits(tm, modelterms7::LinparTerm)) {
        if (!term %in% all.vars(tm@formula)) next
        if (is.null(data)) {
          stop(sprintf(paste0("the variable '%s' is in the linear part of ",
                              "'%s', whose columns are read from the data; ",
                              "give 'data'."), term, p), call. = FALSE)
        }
        cols <- linpar_columns(tm, term, data)
        if (length(cols)) parts[[k]] <- list(name = k, term = tm, cols = cols)
      } else if (term %in% term_variables(tm, k)) {
        parts[[k]] <- list(name = k, term = tm, cols = NULL)
      }
    }
    parts
  })
  names(found) <- params
  found <- Filter(length, found)
  if (!length(found)) {
    avail <- unlist(lapply(params, function(p) names(spec@terms[[p]])))
    stop(sprintf(paste0("'%s' is neither a term nor a variable of %s. ",
                        "The terms are: %s."), term,
                 if (is.null(param)) "the model" else
                   sprintf("the equation '%s'", param),
                 paste(unique(avail), collapse = ", ")), call. = FALSE)
  }
  if (length(found) > 1L) {
    stop(sprintf(paste0("the variable '%s' appears in the equations %s; ",
                        "give 'param'."), term,
                 paste(names(found), collapse = ", ")), call. = FALSE)
  }
  list(param = names(found), vars = term, title = term,
       parts = unname(found[[1L]]))
}


#' The Variables of a Term
#'
#' @description
#' Returns the names of the variables that a term reads: those of the formula
#' for the linear part, and those of the call that names the term otherwise;
#' for a nonlinear term, only those of the call that the build recorded as
#' data, since the call also names the term's parameters.
#'
#' @param tm A model term.
#' @param name The name of the term in its equation.
#'
#' @return A character vector.
#'
#' @keywords internal
term_variables <- function(tm, name) {
  if (S7::S7_inherits(tm, modelterms7::LinparTerm)) {
    return(all.vars(tm@formula))
  }
  out <- tryCatch(all.vars(str2lang(name)), error = function(e) character(0))
  # the call of nl() also names the term's own parameters; the variables of
  # the data are those among the ones the build recorded
  if (S7::S7_inherits(tm, modelterms7::NlTerm)) {
    out <- intersect(out, names(tm@build_info$data_vars))
  }
  setdiff(out, c("TRUE", "FALSE"))
}


#' The Variables That Share a Linear Term with a Variable
#'
#' @description
#' Returns the variables of the terms of a linear formula that involve
#' `var`: `var` alone for a main effect, and the other variables of each
#' interaction in which it appears.
#'
#' @param tm A [modelterms7::LinparTerm].
#' @param var A single string.
#'
#' @return A character vector.
#'
#' @keywords internal
linpar_partners <- function(tm, var) {
  fac <- attr(stats::terms(tm@formula), "factors")
  if (!length(fac) || !var %in% rownames(fac)) return(character(0))
  uses <- fac[var, ] > 0
  rownames(fac)[rowSums(fac[, uses, drop = FALSE] > 0) > 0]
}


#' The Columns of the Linear Part That Involve a Variable
#'
#' @description
#' Returns the positions, within the coefficients of a linear term, of the
#' columns whose model-matrix term involves `var`, read from the `assign`
#' attribute of the model matrix of `data`.
#'
#' @param tm A [modelterms7::LinparTerm].
#' @param var A single string.
#' @param data The data frame of the fit.
#'
#' @return An integer vector.
#'
#' @keywords internal
linpar_columns <- function(tm, var, data) {
  tt <- stats::terms(tm@formula)
  mf <- stats::model.frame(stats::delete.response(tt), data = data,
                           na.action = stats::na.pass)
  mm <- stats::model.matrix(stats::delete.response(tt), mf,
                            contrasts.arg = if (length(tm@contrasts))
                              tm@contrasts else NULL)
  fac <- attr(tt, "factors")
  labs <- attr(tt, "term.labels")
  uses <- vapply(labs, function(l) var %in% rownames(fac)[fac[, l] > 0],
                 logical(1))
  cols <- colnames(mm)[attr(mm, "assign") %in% which(uses)]
  match(cols, modelterms7::term_coef_names(tm))
}


#' The Coefficients and Variance of Part of a Term
#'
#' @description
#' Returns the coefficients of a term in its equation, or of the columns
#' `cols` within it, with their names in [vcov.StatmodFit()] and the block of
#' the variance matrix that they span.
#'
#' @param fit A [StatmodFit()].
#' @param param The equation.
#' @param name The name of the term.
#' @param cols Positions within the term's coefficients, or `NULL`.
#' @param vcov_type A [vcov_method] object.
#'
#' @return A list with `coef` (a numeric vector) and `vcov` (a matrix).
#'
#' @keywords internal
term_coef_block <- function(fit, param, name, cols, vcov_type) {
  des <- statmod_design(fit@spec)
  idx <- des[[param]]$blocks[[name]]
  if (!is.null(cols)) idx <- idx[cols]
  b <- fit@coefficients[[param]][idx]
  nms <- paste0(param, ":", des[[param]]$coef_names[idx])
  V <- vcov(fit, vcov_type = vcov_type, readable = FALSE)
  list(coef = b, vcov = V[nms, nms, drop = FALSE])
}


#' The Contribution of Part of a Term and Its Standard Error
#'
#' @description
#' Evaluates \eqn{X b} for the block of a term on `newdata` (or on the
#' fitting rows when `newdata` is `NULL`), restricted to the columns `cols`,
#' with the standard error \eqn{\sqrt{x_i' V x_i}} of each row.
#'
#' @param fit A [StatmodFit()].
#' @param param,name,cols As in [term_coef_block()].
#' @param tm The term.
#' @param newdata A data frame, or `NULL`.
#' @param vcov_type A [vcov_method] object.
#'
#' @return A list with `fit` and `se`, numeric vectors.
#'
#' @keywords internal
term_contribution <- function(fit, param, name, tm, cols, newdata, vcov_type) {
  X <- if (is.null(newdata)) modelterms7::term_matrix(tm)
    else modelterms7::term_predict(tm, newdata = newdata)
  X <- as.matrix(X)
  if (!is.null(cols)) X <- X[, cols, drop = FALSE]
  cb <- term_coef_block(fit, param, name, cols, vcov_type)
  list(fit = drop(X %*% cb$coef),
       se = sqrt(pmax(rowSums((X %*% cb$vcov) * X), 0)))
}


#' The Working Residuals of an Equation
#'
#' @description
#' Returns \eqn{z_i - \eta_i = s_i / w_i}, the difference between the working
#' response of an equation and its predictor, where \eqn{s_i} is the score of
#' the log-likelihood with respect to the predictor and \eqn{w_i} the
#' information, expected or observed, on the same scale. For a gaussian mean
#' with the identity link it is the ordinary residual.
#'
#' @param fit A [StatmodFit()].
#' @param param The equation.
#' @param information An object from [distributions7::expected_info()] or
#'   [distributions7::observed_info()].
#'
#' @return A numeric vector with one value per observation; `NA` where the
#'   information is not positive.
#'
#' @keywords internal
working_residuals <- function(fit, param, information) {
  spec <- fit@spec
  th <- fit@fitted
  g <- distributions7::distrib_gradient(spec@distrib, spec@response, th,
                                        scale = "link")[[param]]
  key <- paste(param, param, sep = "_")
  h <- if (S7::S7_inherits(information, distributions7::ExpectedInfo)) {
    distributions7::distrib_expected_hessian(spec@distrib, spec@response, th,
                                             scale = "link",
                                             method = information@method)
  } else if (S7::S7_inherits(information, distributions7::ObservedInfo)) {
    distributions7::distrib_hessian(spec@distrib, spec@response, th,
                                    scale = "link")
  } else {
    stop(paste0("'information' must be expected_info() or observed_info()."),
         call. = FALSE)
  }
  w <- -h[[key]]
  out <- g / w
  out[!is.finite(out) | w <= 0] <- NA_real_
  out
}


#' A Grid over the Variables of an Effect
#'
#' @description
#' Builds new data for an effect: `n_grid` points over the range of a
#' numeric variable, or its levels for a factor, crossed with the levels of a
#' factor among `cross` (when it has at most twelve), and every
#' other variable of `vars_all` held at its median or its first level.
#'
#' @param data The data frame of the fit.
#' @param var The variable on the horizontal axis.
#' @param vars_all Every variable the term reads.
#' @param n_grid The number of points for a numeric variable.
#' @param cross The variables that may be crossed with the grid: those that
#'   enter the drawn contribution together with `var`.
#'
#' @return A data frame with the attribute `"group"`, the name of the
#'   crossed factor or `NULL`.
#'
#' @keywords internal
effect_grid <- function(data, var, vars_all, n_grid, cross = vars_all) {
  v <- data[[var]]
  xs <- if (is.numeric(v)) {
    seq(min(v, na.rm = TRUE), max(v, na.rm = TRUE), length.out = n_grid)
  } else {
    f <- as.factor(v)
    factor(levels(f), levels = levels(f))
  }
  others <- setdiff(vars_all, var)
  group <- NULL
  for (o in intersect(others, cross)) {
    if (!is.numeric(data[[o]]) && is.null(group) &&
        nlevels(as.factor(data[[o]])) <= 12L) group <- o
  }
  out <- if (is.null(group)) data.frame(xs) else {
    f <- as.factor(data[[group]])
    expand.grid(xs, factor(levels(f), levels = levels(f)))
  }
  names(out) <- c(var, group)
  for (o in setdiff(others, group)) {
    vo <- data[[o]]
    out[[o]] <- if (is.matrix(vo)) {
      matrix(apply(vo, 2, stats::median, na.rm = TRUE), nrow(out), ncol(vo),
             byrow = TRUE)
    } else if (is.numeric(vo)) stats::median(vo, na.rm = TRUE)
    else {
      f <- as.factor(vo)
      factor(levels(f)[1L], levels = levels(f))
    }
  }
  attr(out, "group") <- group
  out
}


#' Check That a Graph Has Its Data
#'
#' @description
#' Signals an error naming the graph when `data` is `NULL`, and when it lacks
#' one of `vars` or has a number of rows different from the fit.
#'
#' @param ctx The context described at [plot_data()].
#' @param fit A [StatmodFit()].
#' @param type The [plot_type] object.
#' @param vars The variables the graph reads.
#'
#' @return `ctx$data`, invisibly.
#'
#' @keywords internal
plot_needs_data <- function(ctx, fit, type, vars = ctx$vars) {
  d <- ctx$data
  if (is.null(d)) {
    stop(sprintf(paste0("the graph %s reads the covariates, and a statmod ",
                        "fit does not keep its data; give the fitting data ",
                        "as 'data'."), type@label), call. = FALSE)
  }
  miss <- setdiff(vars, names(d))
  if (length(miss)) {
    stop(sprintf("'data' has no variable %s.",
                 paste0("'", miss, "'", collapse = ", ")), call. = FALSE)
  }
  if (nrow(d) != fit@spec@n_obs) {
    stop(sprintf(paste0("'data' has %d rows and the fit %d observations; ",
                        "give the data the model was fitted to."),
                 nrow(d), fit@spec@n_obs), call. = FALSE)
  }
  invisible(d)
}


#' Check a Level and a Number of Grid Points
#'
#' @description
#' Signals an error unless `level` is a single number in (0, 1) and `n_grid`
#' a single whole number of at least 2.
#'
#' @param level,n_grid The values to check.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
check_level_grid <- function(level, n_grid) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(n_grid) || length(n_grid) != 1L || !is.finite(n_grid) ||
      n_grid < 2 || n_grid != round(n_grid)) {
    stop("'n_grid' must be a single whole number of at least 2.",
         call. = FALSE)
  }
  invisible(NULL)
}


#' Diagnostic Graphs of a Fitted Model
#'
#' @description
#' Draws a set of [diagnostic_plot] graphs of a fit and returns them as a
#' named list of ggplot2 objects, or composed in one figure. The default set
#' is the observed density of the response with the densities of simulated
#' responses ([plot_ppc()]), the worm plot ([plot_worm()]), the calibration of
#' the predicted quantiles ([plot_calibration()]) and the quantile residuals
#' against the predicted median ([plot_residuals()] with
#' `against = "median"`). None of the four reads a particular parameter, so
#' the set applies to every family, including those without a location
#' parameter.
#'
#' @param fit A [StatmodFit()].
#' @param plot_types A list of objects inheriting from [diagnostic_plot], or
#'   `NULL` for the default set.
#' @param data The data frame of the fit, or `NULL`; it is required by the
#'   graphs that read a covariate.
#' @param compose Whether to compose the graphs in one figure with
#'   \pkg{patchwork}.
#'
#' @return A named list of ggplot2 objects, or a patchwork object when
#'   `compose = TRUE`.
#'
#' @seealso [autoplot.StatmodFit()], [diagnostic_plot]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150))
#' dd$y <- 1 + dd$x + rnorm(150, sd = 0.4)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' d <- diagnostics(fit)
#' names(d)
#' diagnostics(fit, plot_types = list(plot_qq(), plot_pit()), compose = TRUE)
#' @export
diagnostics <- function(fit, plot_types = NULL, data = NULL, compose = FALSE) {
  if (!S7::S7_inherits(fit, StatmodFit)) {
    stop("'fit' must be a statmod fit.", call. = FALSE)
  }
  if (is.null(plot_types)) {
    plot_types <- list(ppc = plot_ppc(), worm = plot_worm(),
                  calibration = plot_calibration(),
                  residuals = plot_residuals(against = "median"))
  }
  if (!is.list(plot_types) || !length(plot_types) ||
      !all(vapply(plot_types, S7::S7_inherits, logical(1), diagnostic_plot))) {
    stop("'plot_types' must be a list of diagnostic_plot objects, such as plot_qq().",
         call. = FALSE)
  }
  if (is.null(names(plot_types))) {
    labs <- vapply(plot_types, function(t) t@label, character(1))
    names(plot_types) <- sub("^plot_", "", sub("[(].*$", "", labs))
  }
  ctx <- list(param = NULL, name = NULL, vars = NULL, cols = NULL,
              data = data)
  out <- lapply(plot_types, function(t) plot_assemble(t, fit@spec@distrib, fit, ctx))
  if (!compose) return(out)
  # in a composed figure each panel is narrow, so the entries of a legend
  # are stacked one per row
  stacked <- ggplot2::guide_legend(ncol = 1)
  out <- lapply(out, function(p) p +
    ggplot2::guides(colour = stacked, fill = stacked) +
    ggplot2::theme(legend.box = "vertical",
                   legend.spacing.y = ggplot2::unit(0, "pt"),
                   legend.key.spacing.y = ggplot2::unit(0, "pt")))
  patchwork::wrap_plots(out)
}
