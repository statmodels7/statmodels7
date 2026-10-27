#' @include predict.R
#' @title S7 Class for the Meat of a Sandwich Variance
#'
#' @description
#' The class that [hc0_sandwich()] and [hc1_sandwich()] instantiate. It
#' carries the type of the finite-sample correction, the optional cluster
#' variable and the optional data frame from which that variable is read.
#' Use the constructors, which set the properties.
#'
#' @param type A single string, `"HC0"` or `"HC1"`.
#' @param cluster `NULL` or a one-sided formula with one variable, such as
#'   `~ id`.
#' @param data `NULL` or the data frame that holds the cluster variable.
#'
#' @return An S7 object with the properties `type`, `cluster` and `data`.
#'
#' @seealso [hc0_sandwich()], [frequentist_vcov()]
#'
#' @examples
#' class(hc1_sandwich(cluster = ~ id))
#' hc1_sandwich(cluster = ~ id)@type
#'
#' @name SandwichMethod
#' @keywords internal
SandwichMethod <- S7::new_class(
  "SandwichMethod",
  properties = list(
    type = S7::class_character,
    cluster = S7::new_union(NULL, S7::class_formula),
    data = S7::new_union(NULL, S7::class_data.frame),
    reference = S7::new_property(S7::class_character, default = "normal")
  ),
  validator = function(self) {
    if (length(self@type) != 1L ||
        !self@type %in% c("HC0", "HC1", "HC2", "HC3", "HC4")) {
      return("@type must be \"HC0\", \"HC1\", \"HC2\", \"HC3\" or \"HC4\".")
    }
    if (identical(self@type, "HC4") && !is.null(self@cluster)) {
      return("HC4 has no version for clusters.")
    }
    if (length(self@reference) != 1L || !self@reference %in% c("t", "normal")) {
      return("@reference must be \"t\" or \"normal\".")
    }
    cl <- self@cluster
    if (!is.null(cl) && (length(cl) != 2L || !is.name(cl[[2L]]))) {
      return(paste0("must be a one-sided formula naming one variable, such ",
                    "as ~ id."))
    }
    NULL
  }
)

#' @title Sandwich Variance of the Coefficients
#'
#' @description
#' Select the meat of the frequentist variance \eqn{V = V_b B V_b}, with
#' \eqn{V_b = (H + S)^{-1}} the inverse of the penalized information and
#' \eqn{B} built from the scores of the observations. Pass the object as the
#' `sandwich` argument of [frequentist_vcov()].
#'
#' @details
#' The score \eqn{u_i} of observation \eqn{i} is the derivative of its
#' log-likelihood, with no penalty in it, in the whole vector of
#' coefficients: the fixed coefficients and the smooth coefficients of every
#' distribution parameter, and the random effects, which are treated as
#' penalized coefficients. The matrix \eqn{B} is
#' \eqn{\sum_i u_i u_i'}, or \eqn{\sum_g U_g U_g'} with
#' \eqn{U_g = \sum_{i \in g} u_i} over the clusters \eqn{g} when `cluster`
#' is given. The estimator is therefore the robust counterpart of the
#' frequentist variance of [frequentist_vcov()], conditional on the random
#' effects and on the hyperparameters. It does not contain the smoothing
#' bias of a penalized term, as the frequentist variance does not.
#'
#' `hc0_sandwich()` applies no correction (HC0, and CR0 with a cluster).
#' `hc1_sandwich()` multiplies by \eqn{n/(n - p)} (HC1), and by
#' \eqn{G/(G - 1) \cdot (n - 1)/(n - p)} with a cluster (CR1), where \eqn{G}
#' is the number of clusters. A cluster for each observation reduces CR0 to
#' HC0 and CR1 to HC1. Here \eqn{p} is the effective degrees of freedom of
#' the whole model, random effects included, which is not the number of
#' columns of the design where a term is penalized.
#'
#' The cluster variable is read from `data` when that argument is given and
#' from the data of the call to [statmod()] otherwise, which is evaluated
#' again in the environment of the model formula. A fit does not keep its
#' data, so the second route fails once the data frame is gone or has
#' changed.
#'
#' The meat always uses the score of the likelihood. With [ml()] the scores
#' sum to zero at the fit, apart from the penalty. Under [reml()] a
#' coordinate that the marginal criterion estimates, such as the scale of a
#' Gaussian, has a likelihood score that is not zero at the fit, and the
#' sandwich includes that offset.
#'
#' A cluster must contain every group of a [modelterms7::random()] term
#' entirely. With no `cluster` the observations are treated as independent
#' given the random effects.
#'
#' Models with a structural term ([modelterms7::gas()], [modelterms7::regime()]),
#' with a term whose block moves with its coefficients (the break-point terms
#' and [modelterms7::nl()]) are rejected. [hac_sandwich()] covers the
#' structural terms.
#'
#' The prior weights \eqn{w_i} of [statmod()] are frequency weights: row
#' \eqn{i} stands for \eqn{w_i} identical observations. Without a cluster
#' the copies are separate observations and
#' \eqn{B = \sum_i w_i g_i g_i'}, with \eqn{g_i} the score of one copy, and
#' \eqn{n = \sum_i w_i} in the HC1 factor. With a cluster the copies stay in
#' the cluster of the row, so \eqn{U_g = \sum_{i \in g} w_i g_i}, \eqn{G} is
#' the number of clusters with positive weight and \eqn{n = \sum_i w_i}. A
#' row with weight zero contributes nothing and is not counted. The result
#' equals the unweighted estimator on the data expanded to \eqn{w_i} copies
#' of each row, for integer weights.
#'
#' @param cluster `NULL` (the default) or a one-sided formula with a single
#'   variable of the data, such as `~ id`.
#' @param data `NULL` (the default) or the data frame that holds the cluster
#'   variable, with one row for each observation of the fit.
#'
#' @return A `SandwichMethod` object with the properties `type` (`"HC0"` or
#'   `"HC1"`), `cluster` and `data`.
#'
#' @references
#' White, H. (1980). A heteroskedasticity-consistent covariance matrix
#' estimator and a direct test for heteroskedasticity. *Econometrica*, 48(4),
#' 817--838.
#'
#' MacKinnon, J. G. and White, H. (1985). Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' *Journal of Econometrics*, 29(3), 305--325.
#'
#' Liang, K.-Y. and Zeger, S. L. (1986). Longitudinal data analysis using
#' generalized linear models. *Biometrika*, 73(1), 13--22.
#'
#' @seealso [frequentist_vcov()], which takes the object.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150))
#' dd$y <- 1 + 2 * dd$x + rnorm(150, sd = 0.2 + dd$x)
#' fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' sqrt(diag(vcov(fit, vcov_type = frequentist_vcov())))
#' sqrt(diag(vcov(fit, vcov_type = frequentist_vcov(
#'   sandwich = hc0_sandwich()))))
#' sqrt(diag(vcov(fit, vcov_type = frequentist_vcov(
#'   sandwich = hc1_sandwich()))))
#'
#' # observations in clusters
#' dd$id <- rep(1:30, each = 5)
#' fit2 <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
#' vcov(fit2, vcov_type = frequentist_vcov(
#'   sandwich = hc1_sandwich(cluster = ~ id)), param = "mu")
#' @name hc0_sandwich
#' @export
hc0_sandwich <- function(cluster = NULL, data = NULL) {
  SandwichMethod(type = "HC0", cluster = cluster, data = data)
}

#' @rdname hc0_sandwich
#' @export
hc1_sandwich <- function(cluster = NULL, data = NULL) {
  SandwichMethod(type = "HC1", cluster = cluster, data = data)
}

S7::method(print, SandwichMethod) <- function(x, ...) {
  cat("<sandwich> ", x@type,
      if (!is.null(x@cluster)) paste0(" clustered by ", all.vars(x@cluster)),
      "\n", sep = "")
  invisible(x)
}


#' The Cluster Index of a Sandwich Variance
#'
#' @description
#' Reads the cluster variable of a [SandwichMethod] object from its data or
#' from the data of the call that made the fit, checks it, and checks that
#' every group of a random term lies in one cluster.
#'
#' @param fit The [StatmodFit()].
#' @param design The design of the fit.
#' @param sw The [SandwichMethod] object.
#'
#' @return An integer vector of cluster numbers, one for each observation.
#'
#' @keywords internal
sandwich_cluster <- function(fit, design, sw) {
  spec <- fit@spec
  n <- spec@n_obs
  data <- sw@data
  if (is.null(data)) {
    data <- tryCatch({
      env <- environment(spec@formula)
      if (is.null(env)) env <- globalenv()
      eval(fit@call$data, env)
    }, error = function(e) NULL)
  }
  if (!is.data.frame(data) || nrow(data) != n) {
    stop(paste0("The data frame with the cluster variable could not be ",
                "recovered from the call of the fit. Pass it as ",
                "hc0_sandwich(cluster = , data = ) or hc1_sandwich(cluster ",
                "= , data = ), with one row for each of the ", n,
                " observations."), call. = FALSE)
  }
  v <- all.vars(sw@cluster)
  if (!v %in% names(data)) {
    stop(sprintf("The cluster variable '%s' is not a column of the data.", v),
         call. = FALSE)
  }
  cl <- data[[v]]
  if (anyNA(cl)) {
    stop(sprintf("The cluster variable '%s' has missing values.", v),
         call. = FALSE)
  }
  id <- as.integer(factor(cl))
  if (max(id) < 2L) {
    stop("The cluster variable has a single value; at least two clusters are needed.",
         call. = FALSE)
  }
  # every group of a random term in one cluster: a column of the block is the
  # indicator of a group (or a covariate times it), so the rows where the
  # column is nonzero are the group's rows
  for (p in spec@distrib@params) {
    d <- design[[p]]
    if (d$npar == 0L) next
    for (k in names(spec@terms[[p]])) {
      tm <- spec@terms[[p]][[k]]
      if (!S7::S7_inherits(tm, modelterms7::RandomTerm) &&
          !S7::S7_inherits(tm, modelterms7::RandomSmoothTerm)) next
      cols <- d$blocks[[k]]
      if (!length(cols)) next
      Xs <- methods::as(methods::as(d$X[, cols, drop = FALSE], "generalMatrix"),
                        "CsparseMatrix")
      sm <- Matrix::summary(Xs)
      sm <- sm[sm$x != 0, , drop = FALSE]
      pairs <- unique(data.frame(j = sm$j, g = id[sm$i]))
      if (anyDuplicated(pairs$j)) {
        stop(sprintf(paste0("A group of the random term '%s' (equation %s) ",
                            "falls in more than one cluster. Every group of ",
                            "a random term must lie in one cluster."),
                     k, p), call. = FALSE)
      }
    }
  }
  id
}


#' Reject a Fit That the Sandwich Variance Does Not Cover
#'
#' @description
#' Signals an error for a fit with a structural term, a term whose block moves
#' with its coefficients.
#'
#' @details
#' Prior weights are accepted; see [hc0_sandwich()] for their meaning.
#'
#' @param fit The [StatmodFit()].
#' @param design The design of the fit.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
sandwich_check <- function(fit, design) {
  spec <- fit@spec
  if (length(attr(design, "structural")) > 0L) {
    stop(paste0("A sandwich variance is not available for a model with a ",
                "structural term (gas() or regime()). It is rejected."),
         call. = FALSE)
  }
  rf <- attr(design, "refresh")
  if (length(rf) > 0L) {
    stop(paste0("A sandwich variance is not available for a model with a ",
                "term whose block moves with its coefficients (seg(), ",
                "jump(), jseg() or nl()). It is rejected."), call. = FALSE)
  }
  invisible(NULL)
}


#' The Sandwich Variance of the Coefficients
#'
#' @description
#' Computes \eqn{c\,V_b B V_b} for the coordinates in `cols`, with \eqn{V_b}
#' the inverse penalized information that [vcov.StatmodFit()] has built,
#' \eqn{B} the sum over observations (or over clusters) of the outer product
#' of the scores, and \eqn{c} the correction of the type.
#'
#' @param fit The [StatmodFit()].
#' @param design The design of the fit.
#' @param Vb The inverse penalized information over the coordinates `cols`.
#' @param cols The positions, among the coefficients stacked by equation, of
#'   the coordinates that `Vb` covers.
#' @param sw The [SandwichMethod] object.
#'
#' @return A symmetric matrix of the dimension of `Vb`.
#'
#' @seealso [statmod_score_obs()] for the scores,
#'   [vcov.StatmodFit()], which calls it.
#'
#' @keywords internal
sandwich_vcov <- function(fit, design, Vb, cols, sw, expected = FALSE,
                          emethod = NULL) {
  spec <- fit@spec
  n <- spec@n_obs
  id <- if (is.null(sw@cluster)) NULL else sandwich_cluster(fit, design, sw)
  edf <- fit@edf
  p <- if (is.data.frame(edf)) sum(edf$edf) else NA_real_
  if (!is.finite(p)) {
    stop("The effective degrees of freedom of the fit are not available.",
         call. = FALSE)
  }
  w <- spec@weights
  g <- statmod_score_obs(spec, fit@coefficients, design)
  df_fn <- NULL
  if (sw@type %in% c("HC2", "HC3", "HC4")) {
    adj <- sandwich_adjusted(fit, design, Vb, cols, g, id, sw, expected,
                             emethod)
    g <- adj$g
    df_fn <- adj$df_fn
  }
  # g already carries w_i. Without a cluster the meat is sum w_i g_i g_i',
  # the meat of w_i separate copies of row i, so each row enters as
  # g_i sqrt(w_i) = (w_i g_i) / sqrt(w_i). With a cluster the copies stay in
  # their cluster, and U_g sums w_i g_i as it stands.
  if (is.null(id)) {
    sw_row <- ifelse(w > 0, 1 / sqrt(w), 0)
    g <- lapply(g, function(v) v * sw_row)
  }
  n <- sum(w)
  blocks <- list()
  for (q in spec@distrib@params) {
    d <- design[[q]]
    if (d$npar == 0L) next
    blocks[[length(blocks) + 1L]] <- if (isS4(d$X)) {
      Matrix::Diagonal(x = g[[q]]) %*% d$X
    } else d$X * g[[q]]
  }
  U <- if (any(vapply(blocks, isS4, TRUE))) do.call(Matrix::cbind2, blocks) else
    do.call(cbind, blocks)
  U <- U[, cols, drop = FALSE]
  if (!is.null(id)) {
    G <- length(unique(id[w > 0]))
    A <- Matrix::sparseMatrix(i = id, j = seq_along(id), x = 1,
                              dims = c(max(id), length(id)))
    U <- A %*% U
  }
  B <- as.matrix(Matrix::crossprod(U))
  scale <- if (sw@type %in% c("HC0", "HC2", "HC3", "HC4")) 1 else if (is.null(id)) {
    if (n - p <= 0) stop("HC1 needs fewer effective parameters than observations.",
                         call. = FALSE)
    n / (n - p)
  } else {
    if (n - p <= 0) stop("CR1 needs fewer effective parameters than observations.",
                         call. = FALSE)
    G / (G - 1) * (n - 1) / (n - p)
  }
  V <- scale * (Vb %*% B %*% Vb)
  V <- (V + t(V)) / 2
  if (!is.null(df_fn)) attr(V, "df_fn") <- df_fn
  V
}
