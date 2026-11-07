#' The Compiled Alternation of a Model with a Kinked Penalty
#'
#' @description
#' `engine_alternate_fit()` runs [statmod_alternate()] in compiled code
#' (`src/engine_kink.h`) for a model whose penalty has a kink: at every pass
#' the curvature of a scaled SCAD or MCP, the smooth block fitted on a subset
#' as [fit_smooth()] fits it, every kinked block fitted by the coordinate
#' descent of [coord_fit()], and at the end of the passes the reading of
#' [alternation_settled()]. It returns the list that [statmod_alternate()]
#' returns, or `NULL`, in which case [statmod_alternate()] runs the passes in
#' R from the same coefficients.
#'
#' @details
#' A model is covered when [engine_route()] covers its smooth block with the
#' closed-form expected information and `emethod` is `NULL`; nothing is
#' printed (`vb$blocks` and `vb$inner` are `FALSE`); the model has one
#' thread, a dense design and no structural, moving or frozen term; and every
#' kinked block is a lasso, an elastic net, a SCAD or an MCP over a dense
#' block, with a diagonal map or none, whose coordinates no other penalized
#' unit covers and none of which is held. Where the compiled alternation
#' meets a branch that it does not reproduce (the proximal fallback of
#' [sparse_fit()]) or a point where the R route signals an error, it
#' returns `NULL`.
#'
#' `engine_alternate_blocks()` describes the kinked blocks to the compiled
#' alternation: for each, its penalized unit, its equation, its columns and
#' the other columns of the equation, whether it is solved with the
#' equation's intercept profiled out ([coord_centers()]), whether
#' [statmod_curv()] writes its curvature and with which centering, the
#' previous point's kink on a path, and the columns themselves.
#'
#' @param spec,design,blocks,hyper,inner_optimizer,beta,expected,emethod,maxit,tol,vb
#'   As for [statmod_alternate()].
#' @param obj The objective at the start, [statmod_objective()].
#'
#' @return `engine_alternate_fit()`: the list of [statmod_alternate()], or
#'   `NULL`. `engine_alternate_blocks()`: a list of block descriptions, or
#'   `NULL` where a block is not covered.
#'
#' @seealso [statmod_alternate()], [engine_inner_fit()]
#' @keywords internal
engine_alternate_fit <- function(spec, design, blocks, hyper, inner_optimizer,
                                 beta, expected, emethod, maxit, tol, vb, obj) {
  if (!identical(spec@engine, "auto") || !is.null(emethod)) return(NULL)
  if (isTRUE(vb$blocks) || isTRUE(vb$inner)) return(NULL)
  if (!S7::S7_inherits(inner_optimizer, Iwls)) return(NULL)
  if (!identical(as.integer(spec@threads), 1L)) return(NULL)
  if (!length(blocks$sparse) || design_sparse(design)) return(NULL)
  method <- iwls_resolve(inner_optimizer, spec@distrib)
  rt <- engine_route(spec, design, obj, beta, hyper, method)
  if (is.null(rt) || !identical(rt$method$emethod, 0L)) return(NULL)
  if (any(vapply(rt$pen$units, function(u) u$kind %in% c(9L, 10L), logical(1)))) {
    return(NULL)
  }
  units <- statmod_penalized(spec, design)
  hf <- held_positions(spec, design, obj, beta)
  kb <- engine_alternate_blocks(spec, design, blocks, units, obj, beta, hf)
  if (is.null(kb)) return(NULL)
  lab <- tryCatch(coef_labels(spec, design), error = function(e) NULL)
  if (is.null(lab) || nrow(lab) != length(beta)) return(NULL)
  idx <- blocks$smooth
  gfull <- obj$split(seq_along(beta))
  groups <- Filter(length, lapply(gfull, function(ix) match(intersect(ix, idx), idx)))
  frozen <- match(intersect(hf$where, idx), idx)
  rt$method$damp_on_reject <- !has_sharp_breakpoint(spec)
  rt$method$guard <- FALSE
  res <- engine_alternate(rt$model, rt$pen, rt$method, as.numeric(beta), kb,
                          as.integer(idx), lapply(groups, as.integer),
                          as.integer(frozen), as.integer(hf$where),
                          as.numeric(hf$value), as.logical(lab$kinked),
                          as.integer(pinned_coords(spec, design)),
                          isTRUE(expected), as.integer(maxit), as.numeric(tol))
  if (!identical(res$status, 0L)) return(NULL)

  # the specification with the curvatures the last pass wrote, and the
  # objective on it, as statmod_alternate() leaves them
  for (b in seq_along(kb)) {
    if (!isTRUE(kb[[b]]$curv_todo)) next
    u <- units[[kb[[b]]$unit]]
    pen <- spec@terms[[u$param]][[u$term]]@penalty
    cv <- as.numeric(res$curvature[[b]])
    if (!identical(as.numeric(pen@curvature), cv)) {
      spec@terms[[u$param]][[u$term]]@penalty <- S7::set_props(pen, curvature = cv)
    }
  }
  if (any(vapply(kb, function(k) isTRUE(k$curv_todo), logical(1)))) {
    obj <- statmod_objective(spec, hyper, design, expected, emethod)
  }
  note <- NULL
  if (isTRUE(res$note$any)) {
    nk <- res$note
    note <- switch(as.character(nk$kind),
      "1" = paste0("the score is not finite at iteration ", nk$it,
                   ", so the run stopped at the last usable point"),
      "2" = paste0("the line search found no acceptable step at iteration ",
                   nk$it, if (isTRUE(rt$method$damp_on_reject)) ", even with damping"),
      "3" = paste0("the objective stopped moving at iteration ", nk$it),
      "4" = paste0("converged by the dimensionless rule (score ",
                   format(nk$dimless, digits = 3), " against ",
                   format(nk$score, digits = 3),
                   " absolute; the response's scale is the gap)"),
      NULL)
  }
  hb <- res$blocks
  hist_blocks <- if (length(hb$pass)) {
    list(data.frame(pass = hb$pass, block = hb$block, objective = hb$objective,
                    change = hb$change, iterations = hb$iterations,
                    converged = hb$converged))
  } else list()
  hi <- res$inner
  hist_inner <- if (length(hi$pass)) {
    list(data.frame(iteration = hi$iteration, objective = hi$objective,
                    score = hi$score, step = hi$step, rank = hi$rank,
                    route = hi$route, held = hi$held, damp = hi$damp,
                    pass = hi$pass, block = "smooth"))
  } else list()
  list(par = res$par, value = res$value, converged = isTRUE(res$converged),
       obj = obj, note = note, aliased = sort(as.integer(res$aliased)),
       hist_blocks = hist_blocks, hist_inner = hist_inner, spec = spec)
}

#' @rdname engine_alternate_fit
#' @param units The penalized units, [statmod_penalized()].
#' @param hf The held positions, [held_positions()].
#' @keywords internal
engine_alternate_blocks <- function(spec, design, blocks, units, obj, beta, hf) {
  keys <- vapply(units, function(u) paste(u$param, u$key, sep = "\r"), character(1))
  kinked_pos <- integer(0)
  out <- vector("list", length(blocks$sparse))
  for (b in seq_along(blocks$sparse)) {
    bl <- blocks$sparse[[b]]
    if (isTRUE(bl$structural)) return(NULL)
    k <- match(paste(bl$param, bl$term, sep = "\r"), keys)
    if (is.na(k)) return(NULL)
    u <- units[[k]]
    pen <- u$penalty
    if (!(S7::S7_inherits(pen, penalties7::ScadPenalty) ||
          S7::S7_inherits(pen, penalties7::McpPenalty) ||
          S7::S7_inherits(pen, penalties7::DistribPenalty))) return(NULL)
    if (!identical(as.integer(u$index), as.integer(bl$index))) return(NULL)
    if (length(intersect(hf$where, bl$index))) return(NULL)
    p <- bl$param
    d <- design[[p]]
    if (!is.matrix(d$X) || !is.double(d$X)) return(NULL)
    X <- d$X[, bl$cols, drop = FALSE]
    other <- setdiff(seq_len(d$npar), bl$cols)
    ctr <- coord_centers(spec, design, p, obj, beta)
    i0 <- if (ctr) obj$split(seq_along(beta))[[p]][parametric_intercept(spec, design, p)]
          else NA_integer_
    scad_mcp <- S7::S7_inherits(pen, penalties7::ScadPenalty) ||
      S7::S7_inherits(pen, penalties7::McpPenalty)
    todo <- scad_mcp && is.null(u$class) && !isTRUE(u$structural) &&
      !is.null(u$index) && identical(u$key, u$term)
    has_curv <- "curvature" %in% S7::prop_names(pen) && length(pen@curvature) > 0L
    kinked_pos <- c(kinked_pos, bl$index)
    out[[b]] <- list(unit = k, param = match(p, spec@distrib@params),
                     cols = as.integer(bl$cols), other = as.integer(other),
                     ctr = isTRUE(ctr), i0 = as.integer(i0), curv_todo = todo,
                     curv_centered = isTRUE(todo && curv_centered(spec, design, p)),
                     has_curv = has_curv,
                     prev_kink = if (is.null(bl$prev_kink)) NULL else as.numeric(bl$prev_kink),
                     label = paste(bl$param, bl$term, sep = "/"),
                     X = unname(X),
                     Xo = if (length(other)) unname(d$X[, other, drop = FALSE]) else NULL)
  }
  # a kinked unit's coordinates are covered by no other unit
  for (u in units) {
    if (any(vapply(blocks$sparse, function(bl) identical(u$param, bl$param) &&
                     identical(u$key, bl$term), logical(1)))) next
    if (length(intersect(u$index, kinked_pos))) return(NULL)
  }
  out
}
