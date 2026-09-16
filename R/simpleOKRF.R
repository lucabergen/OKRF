##' A simple implementation of Output Kernel Random Forests using forest-wise
##' feature approximations based on the R package \code{simpleRF}.
##' Uses reference classes and only plain \code{R}.
##' Not optimized for computation speed.
##'
##' Unordered factor variables can be handled in different ways.
##' Use "ignore" to treat them as ordered in the order of the factor levels.
##' With "partition" all 2-partitions of the factor levels are considered for splitting.
##'
##' @title simpleOKRF
##' @param K Gram matrix of targets. If possible, use Phi instead of K.
##' @param Phi Feature matrix of targets (n \times d, n > d).
##' @param X Covariate data of class \code{data.frame}, with one row per target observation.
##' @param tol Error tolerance of approximation. Default 0.001.
##' @param num_trees Number of trees.
##' @param mtry Number of variables to possibly split at in each node.
##' @param min_node_size Minimal node size. Default 5.
##' @param replace Sample with replacement. Default TRUE.
##' @param unordered_factors How to handle unordered factor variables. Either "ignore" or "partition" with default "ignore".
##' @param num_threads Number of threads used for mclapply, set to 1 for debugging.
##' @examples
##' \donttest{
##' # TODO: Include new example
##' }
##' @author Luca Bergen
##' @import stats
##' @export
simpleOKRF <- function(K = NULL, Phi = NULL, X,
                       tol = 1e-3, num_trees = 200, mtry = NULL,
                       min_node_size = NULL, replace = TRUE,
                       unordered_factors = "ignore",
                       num_threads = 1) {

  ## Check parameters

  # Covariate df
  if (!is.data.frame(X)) {
    stop("X must be a data.frame.")
  }
  if (nrow(X) < 2L) {
    stop("X must contain at least two observations.")
  }
  if (ncol(X) < 1L) {
    stop("X must contain at least one covariate.")
  }
  if (anyNA(X)) {
    stop("Missing values in split covariates are not supported.")
  }

  # Feature representation
  if (is.null(K) && is.null(Phi)) {
    stop("Either K or Phi must be specified.")
  }
  if (!is.null(K) && !is.null(Phi)) {
    warning("Both K and Phi specified. Only Phi is used.")
  }

  if (!is.null(Phi)) {

    if (!is.matrix(Phi) || !is.numeric(Phi)) {
      stop("Phi must be a numeric matrix.")
    }
    if (nrow(Phi) != nrow(X)) {
      stop("Phi must have the same numer of rows as X.")
    }
    if (nrow(Phi) <= ncol(Phi)) {
      stop("Phi must have fewer columns than rows.")
    }

    feat_rep <- list(
      type = "explicit",
      Phi = Phi
    )

  } else {

    if (!is.matrix(K) || !is.numeric(K)) {
      stop("K must be a numeric matrix.")
    }
    if (nrow(K) != ncol(K)) {
      stop("K must be square.")
    }
    if (!isSymmetric(K)) {
      stop("K must be symmetric.")
    }
    if (nrow(K) != nrow(X)) {
      stop("K must have the same number of rows and columns as X.")
    }

    feat_rep <- list(
      type = "implicit",
      K = K
    )

  }

  stopifnot(
    length(tol) == 1L,
    is.finite(tol),
    tol > 0
  )

  if (is.null(mtry)) {
    mtry <- sqrt(ncol(X))
  } else if (mtry > ncol(X)) {
    stop("Mtry cannot be larger than number of independent variables.")
  }
  if (is.null(min_node_size)) {
    min_node_size <- 5
  }

  ## Unordered factors
  if (!(unordered_factors %in% c("ignore", "partition"))) {
    stop("Unknown value for unordered_factors.")
  }

  ##  TODO: Add checks and give informative error messages for other params

  if (unordered_factors == "ignore") {
    ## Just set to ordered if "ignore"
    character.idx <- vapply(X, is.character, logical(1))
    ordered.idx <- vapply(X, is.ordered, logical(1))
    factor.idx <- vapply(X, is.factor, logical(1))
    recode.idx <- character.idx | (factor.idx & !ordered.idx)

    ## `drop = FALSE` is important when exactly one column is recoded.
    if (any(recode.idx)) {
      X[, recode.idx] <- lapply(
        X[, recode.idx, drop = FALSE],
        as.ordered
      )
    }
  }

  ## Create forest object
  forest <- Forest$new(num_trees = as.integer(num_trees),
                         mtry = as.integer(mtry),
                         min_node_size = as.integer(min_node_size),
                         replace = replace,
                         x_data = Data$new(data = X),
                         unordered_factors = unordered_factors,
                         feat_rep = feat_rep,
                         tol = tol)

  ## Grow forest
  forest$grow(num_threads = num_threads)

  ## Return forest
  return(forest)
}
