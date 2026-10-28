#' @include bootstrap_engine.R
NULL

#' Reject a Random-Effect Term That the Cluster Bootstrap Cannot Regroup
#'
#' @description
#' Signals an error for a fit whose [modelterms7::random()] terms the cluster
#' bootstrap of [cluster_resample()] cannot rebuild on a resampled data set.
#'
#' @details
#' A replicate rebuilds every random-effect term on its own groups, with the
#' `build_info` of the fit for the within-group design. A term that shares a
#' covariance through a label, a smooth random effect
#' ([modelterms7::RandomSmoothTerm()]) and a random effect inside the
#' subformula of another term are not rebuilt in this version, and the fit is
#' rejected.
#'
#' @param spec0 The specification of the fit.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
boot_cluster_random_check <- function(spec0) {
  for (p in names(spec0@terms)) {
    for (k in names(spec0@terms[[p]])) {
      tm <- spec0@terms[[p]][[k]]
      if (!S7::S7_inherits(tm, modelterms7::RandomTerm)) next
      if (!is.na(modelterms7::term_tag(tm))) {
        stop(sprintf(paste0("cluster_resample() is not available beside '%s': ",
                            "the term shares a covariance with other terms ",
                            "through a label, and a replicate would have to ",
                            "rebuild the shared block."), k), call. = FALSE)
      }
      if (S7::S7_inherits(tm, modelterms7::RandomSmoothTerm)) {
        stop(sprintf(paste0("cluster_resample() is not available beside '%s': ",
                            "the smooth within-group design of the term is ",
                            "not rebuilt on the regrouped data."), k),
             call. = FALSE)
      }
    }
  }
  if (nrow(attr(random_modes(spec0, "zero"), "nested")) > 0L) {
    stop(paste0("cluster_resample() is not available beside a random() term ",
                "written inside the subformula of another term."),
         call. = FALSE)
  }
  invisible(NULL)
}


#' The Draw of a Cluster Replicate Beside Random Effects
#'
#' @description
#' Returns the function of [bootstrap_plan()] for a cluster bootstrap of a fit
#' with [modelterms7::random()] terms. A draw is the vector of the clusters
#' drawn, and the attribute `build` of the result turns a draw into the
#' specification of the replicate.
#'
#' @details
#' The clusters are drawn with replacement, \eqn{G} of \eqn{G}. A cluster drawn
#' \eqn{k} times enters the replicate \eqn{k} times, once for each draw, and
#' each copy is a new cluster. The grouping variable of every random-effect
#' term is relabeled in the replicate so that the groups of two copies are
#' different groups, each with an effect of its own, as in Stata's
#' `bootstrap, cluster() idcluster()`. A group is relabeled when every one of
#' its observations lies in a single cluster; a group that spans several
#' clusters keeps its label, and its effect is shared by the copies.
#'
#' The design of the replicate is rebuilt on the rows drawn. The terms other
#' than the random-effect terms are reapplied with [modelterms7::term_predict()],
#' so a basis keeps the knots of the fit and a factor its contrasts. A random
#' effect term is built again with [modelterms7::term_build()] on the
#' relabeled groups; its hyperparameters are those of the fit. The weights,
#' the offsets and the response are those of the rows drawn.
#'
#' The coordinates that a replicate has in common with the fit are those of the
#' terms that are not random effects. `build()` returns their positions in the
#' replicate (`map`, with `NA` for a coordinate of the fit that has no
#' counterpart), the starting values of the replicate, and `expand()`, which
#' returns the coefficients of a refit in the stacked positions of the fit with
#' `NA` at the random effects. A draw that does not reproduce the design of the
#' fit when every cluster is drawn once is an error: it means that a variable
#' of the model is not a column of the data.
#'
#' @param object The [StatmodFit()].
#' @param resample A [ClusterResample].
#' @param spec0,design0 Its specification and design.
#'
#' @return A function of no arguments returning a list with `pick`, with
#'   attributes `kind` (`"regrouped"`), `targets`, `random` and `build`.
#'
#' @keywords internal
boot_cluster_random_plan <- function(object, resample, spec0, design0) {
  boot_cluster_random_check(spec0)
  data <- boot_cluster_data(object, resample, spec0)
  id <- boot_cluster_ids(object, resample, spec0)
  n <- spec0@n_obs
  gi <- match(id, unique(id))
  ng <- max(gi)
  cl_rows <- split(seq_len(n), factor(gi, levels = seq_len(ng)))
  params <- spec0@distrib@params
  y0 <- spec0@response
  w0 <- spec0@weights

  rt <- list()
  for (p in names(spec0@terms)) {
    for (k in names(spec0@terms[[p]])) {
      tm <- spec0@terms[[p]][[k]]
      if (!S7::S7_inherits(tm, modelterms7::RandomTerm)) next
      gv <- eval(tm@build_info$gexpr, data, baseenv())
      per <- tapply(gi, as.character(gv), function(z) length(unique(z)))
      rt[[length(rt) + 1L]] <- list(param = p, key = k, gv = gv,
                                    nested = all(per == 1L))
    }
  }
  is_rand <- lapply(spec0@terms, function(tms)
    vapply(tms, S7::S7_inherits, logical(1), modelterms7::RandomTerm))
  widths <- function(tms) vapply(tms, function(tm)
    length(modelterms7::term_coef_names(tm)), integer(1))
  idx0 <- list()
  off <- 0L
  for (p in names(spec0@terms)) {
    w <- widths(spec0@terms[[p]])
    idx0[[p]] <- Map(function(a, b) off + a + seq_len(b),
                     c(0L, cumsum(w)[-length(w)]), w)
    off <- off + sum(w)
  }
  n0 <- off
  beta0 <- unlist(object@coefficients[params], use.names = FALSE)
  slice <- function(y, rows) {
    out <- if (is.matrix(y)) y[rows, , drop = FALSE] else y[rows]
    for (a in setdiff(names(attributes(y)), c("dim", "dimnames", "names"))) {
      attr(out, a) <- attr(y, a)
    }
    out
  }

  build <- function(pl) {
    pick <- pl$pick
    rows <- unlist(cl_rows[pick], use.names = FALSE)
    inst <- rep(seq_along(pick), lengths(cl_rows[pick]))
    db <- data[rows, , drop = FALSE]
    terms_b <- spec0@terms
    for (j in seq_along(rt)) {
      r <- rt[[j]]
      col <- paste0(".boot_group", j)
      db[[col]] <- if (r$nested) paste(r$gv[rows], inst, sep = "#") else
        as.character(r$gv[rows])
      tm <- terms_b[[r$param]][[r$key]]
      within <- tm@formula[[2L]][[2L]]
      tm@formula <- stats::as.formula(
        call("~", call("|", within, as.name(col))),
        env = environment(tm@formula))
      terms_b[[r$param]][[r$key]] <- modelterms7::term_build(tm, db)
    }
    dist <- spec0@distrib
    if ("size" %in% S7::prop_names(dist) && length(dist@size) > 1L) {
      dist <- S7::set_props(dist, size = dist@size[rows])
    }
    spec_b <- S7::set_props(
      spec0, response = slice(y0, rows), n_obs = length(rows),
      weights = w0[rows], distrib = dist,
      offsets = lapply(spec0@offsets, function(o)
        if (length(o) == n) o[rows] else o),
      newdata = db, terms = terms_b)
    map <- rep(NA_integer_, n0)
    offb <- 0L
    for (p in names(terms_b)) {
      wb <- widths(terms_b[[p]])
      for (q in seq_along(wb)) {
        if (!is_rand[[p]][[q]]) map[idx0[[p]][[q]]] <- offb + seq_len(wb[q])
        offb <- offb + wb[q]
      }
    }
    ok <- !is.na(map)
    start <- numeric(offb)
    start[map[ok]] <- beta0[ok]
    expand <- function(cf) {
      v <- unlist(cf[params], use.names = FALSE)
      out <- rep(NA_real_, n0)
      out[ok] <- v[map[ok]]
      out
    }
    list(spec = spec_b, map = map, common = map[ok], start = start,
         expand = expand, rows = rows)
  }

  chk <- build(list(pick = seq_len(ng)))
  db0 <- statmod_design(spec0)
  dbb <- statmod_design(chk$spec)
  same <- isTRUE(all.equal(as.numeric(slice(y0, chk$rows)),
                           as.numeric(chk$spec@response), tolerance = 0))
  for (p in params) {
    for (k in names(spec0@terms[[p]])) {
      if (is_rand[[p]][[k]]) next
      a <- as.matrix(db0[[p]]$X[chk$rows, db0[[p]]$blocks[[k]], drop = FALSE])
      b <- as.matrix(dbb[[p]]$X[, dbb[[p]]$blocks[[k]], drop = FALSE])
      same <- same && identical(dim(a), dim(b)) &&
        max(abs(a - b), 0) <= 1e-12 * (1 + max(abs(a), 0))
    }
  }
  if (!same) {
    stop(paste0("The design of the fit could not be rebuilt from the data of ",
                "the call: a variable of the model is not a column of the ",
                "data. Pass the data as cluster_resample(cluster = , data = )."),
         call. = FALSE)
  }

  out <- function() list(pick = sample.int(ng, ng, TRUE))
  attr(out, "kind") <- "regrouped"
  attr(out, "targets") <- list()
  attr(out, "random") <- "regrouped"
  attr(out, "build") <- build
  out
}


#' The Status of One Refitted Cluster Replicate With Regrouped Effects
#'
#' @description
#' Returns 0 for a replicate that is kept, 2 when its mode error exceeds
#' [mode_error_limit()] or cannot be read, and 3 when it leaves unidentified a
#' coordinate that the fit identifies.
#'
#' @details
#' It is [boot_validate()] with the stacked positions translated: the pinned
#' coordinates and the unidentified ones of the fit are mapped to the positions
#' of the replicate, and only the coordinates that the replicate shares with
#' the fit are compared. The effects of the regrouped terms have no
#' counterpart in the fit and are not compared.
#'
#' @param object The [StatmodFit()].
#' @param spec_b The specification of the replicate.
#' @param r The refit, as [bootstrap_refit()] returns it.
#' @param cfg The inner settings of the fit.
#' @param ref The reading of the fit itself, from [boot_read_mode()].
#' @param cb The result of the `build` function of [boot_cluster_random_plan()].
#'
#' @return A single integer.
#'
#' @keywords internal
boot_cluster_validate <- function(object, spec_b, r, cfg, ref, cb) {
  pin <- cb$map[as.integer(object@methods$pinned)]
  pin <- pin[!is.na(pin)]
  rd <- boot_read_mode(object, spec_b, r$design, r$coefficients, r$hyper, cfg,
                       pinned = pin)
  if (isTRUE(rd$singular)) return(if (isTRUE(ref$singular)) 0L else 3L)
  def_ref <- cb$map[ref$def]
  def_ref <- def_ref[!is.na(def_ref)]
  if (length(setdiff(intersect(rd$def, cb$common), def_ref))) return(3L)
  if (!is.finite(rd$me) || rd$me > mode_error_limit()) return(2L)
  0L
}
