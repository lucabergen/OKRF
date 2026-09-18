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
##' @param K Gram matrix of targets (dim. \eqn{n \times n}). If possible, use Phi instead of K.
##' @param Phi Feature matrix of targets (dim. \eqn{n \times d, n > d}).
##' @param X Covariate data of class \code{data.frame}, with one row per target observation.
##' @param tol Error tolerance of approximation. Default 0.001.
##' @param num_trees Number of trees. Must be a positive integer.
##' @param mtry Number of covariates considered at each split.
##'   Default is the (rounded down) square root of the number of covariates.
##' @param min_node_size Minimum number of observations required in a node
##'   before attempting a split. Default 5.
##' @param max_depth Maximum tree depth. Default `NULL` means no explicit depth limit.
##'   The root node has depth 0.
##' @param min_leaf_size Minimum number of observations in a leaf node.
##'   Default 1.
##' @param replace Whether observations are sampled with replacement. Default `TRUE`.
##' @param sample_fraction Fraction of training observations sampled per tree. Must be in `(0, 1]`.
##'   Default is 1 for sampling with replacement and 0.632 for sampling without replacement.
##' @param unordered_factors How to handle unordered factor variables.
##'   Either "ignore" or "partition" with default "ignore".
##' @param num_threads Number of threads used for mclapply, set to 1 for debugging.
##' @examples
##' \donttest{
##' set.seed(20260916)
##'
##' ## Covariates
##' n <- 150
##'
##' X <- data.frame(x1 = runif(n, -3, 3), x2 = runif(n, -3, 3))
##'
##' ## Scalar output signal
##' y <- sin(X$x1) + 0.5 * cos(2 * X$x2) + rnorm(n)
##'
##' ## Random Fourier Features for a Gaussian output kernel
##' num_frequencies <- 40
##' ell <- 1
##'
##' omega <- rnorm(num_frequencies, mean = 0, sd = 1 / ell)
##'
##' projection <- outer(y, omega)
##'
##' Phi <- sqrt(1 / num_frequencies) * cbind(cos(projection), sin(projection))
##'
##' ## Fit an Output Kernel Random Forest
##' forest <- simpleOKRF(
##'   Phi = Phi,
##'   X = X,
##'   tol = 1e-6,
##'   num_trees = 25,
##'   mtry = 1,
##'   min_node_size = 5,
##'   replace = TRUE,
##'   unordered_factors = "ignore",
##'   num_threads = 1
##' )
##'
##' ## Predict the output representation for new covariates
##' newdata <- data.frame( x1 = c(-2, 0, 2), x2 = c(-2, 0, 2))
##'
##' Phi_prediction <- forest$predict(
##'   newdata = newdata,
##'   type = "response"
##' )
##'
##' dim(Phi_prediction)
##' }
##' @author Luca Bergen
##' @import stats
##' @export
simpleOKRF <- function(
    K = NULL,
    Phi = NULL,
    X,
    tol = 1e-3,
    scope = c("tree", "forest"),
    num_trees = 200L,
    mtry = floor(sqrt(ncol(X))),
    min_node_size = 5L,
    max_depth = NULL,
    min_leaf_size = 1L,
    replace = TRUE,
    sample_fraction = ifelse(replace, 1, 0.632),
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

  # num_trees
  if (
    length(num_trees) != 1L ||
    !is.numeric(num_trees) ||
    !is.finite(num_trees) ||
    num_trees < 1 ||
    num_trees != as.integer(num_trees)
  ) {
    stop("num_trees must be a positive integer.")
  }

  num_trees <- as.integer(num_trees)

  # mtry
  if (
    length(mtry) != 1L ||
    !is.numeric(mtry) ||
    !is.finite(mtry) ||
    mtry < 1 ||
    mtry != as.integer(mtry) ||
    mtry > ncol(X)
  ) {
    stop(
      "mtry must be an integer between 1 and the number of covariates."
    )
  }

  mtry <- as.integer(mtry)

  # min_node_size
  if (
    length(min_node_size) != 1L ||
    !is.numeric(min_node_size) ||
    !is.finite(min_node_size) ||
    min_node_size < 1 ||
    min_node_size != as.integer(min_node_size)
  ) {
    stop("min_node_size must be a positive integer.")
  }

  min_node_size <- as.integer(min_node_size)

  # max_depth
  if (is.null(max_depth)) {
    # NA_integer_ means no limitation
    max_depth <- NA_integer_

  } else if (
    length(max_depth) != 1L ||
    !is.numeric(max_depth) ||
    !is.finite(max_depth) ||
    max_depth < 0 ||
    max_depth != as.integer(max_depth)
  ) {
    stop("max_depth must be NULL or a non-negative integer.")

  } else {
    max_depth <- as.integer(max_depth)
  }

  # min_leaf_size
  if (
    length(min_leaf_size) != 1L ||
    !is.numeric(min_leaf_size) ||
    !is.finite(min_leaf_size) ||
    min_leaf_size < 1 ||
    min_leaf_size != as.integer(min_leaf_size)
  ) {
    stop("min_leaf_size must be a positive integer.")
  }

  min_leaf_size <- as.integer(min_leaf_size)

  # replace
  if (
    !is.logical(replace) ||
    length(replace) != 1L ||
    is.na(replace)
  ) {
    stop("replace must be TRUE or FALSE.")
  }

  # sample_fraction
  if (
    length(sample_fraction) != 1L ||
    !is.numeric(sample_fraction) ||
    !is.finite(sample_fraction) ||
    sample_fraction <= 0 ||
    sample_fraction > 1
  ) {
    stop("sample_fraction must be a number in (0, 1].")
  }

  # K and Phi
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
      stop("Phi must have the same number of rows as X.")
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

  # tol
  if (
    length(tol) != 1L ||
    !is.numeric(tol) ||
    !is.finite(tol) ||
    tol <= 0
  ) {
    stop("tol must be a positive finite number.")
  }

  # scope
  scope <- match.arg(scope, c("tree", "forest"))

  # num_threads
  if (
    length(num_threads) != 1L ||
    !is.numeric(num_threads) ||
    !is.finite(num_threads) ||
    num_threads < 1 ||
    num_threads != as.integer(num_threads)
  ) {
    stop("num_threads must be a positive integer.")
  }

  num_threads <- as.integer(num_threads)

  # unordered_factors
  if (!(unordered_factors %in% c("ignore", "partition"))) {
    stop("Unknown value for unordered_factors.")
  }

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
  forest <- Forest$new(
    num_trees = num_trees,
    mtry = mtry,
    min_node_size = min_node_size,
    max_depth = max_depth,
    min_leaf_size = min_leaf_size,
    replace = replace,
    sample_fraction = sample_fraction,
    x_data = Data$new(data = X),
    unordered_factors = unordered_factors,
    feat_rep = feat_rep,
    tol = tol,
    scope = scope
  )

  ## Grow forest
  forest$grow(num_threads = num_threads)

  ## Return forest
  return(forest)
}
