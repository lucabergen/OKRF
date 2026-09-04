##' A simple implementation of Output Kernel Random Forests using tree-wise
##' feature approximations based on the R package \code{simpleRF}.
##' Uses reference classes and only plain \code{R}.
##' Not optimized for computation speed.
##'
##' Unordered factor variables can be handled in different ways.
##' Use "ignore" to treat them as ordered in the order of the factor levels.
##' With "partition" all 2-partitions of the factor levels are considered for splitting.
##'
##' @title simpleOKRF
##' @param K Gram matrix of targets
##' @param Phi Feature matrix of targets
##' @param data Covariate data of class \code{data.frame}.
##' @param eps Error tolerance of approximation. Default 0.001.
##' @param num_trees Number of trees.
##' @param mtry Number of variables to possibly split at in each node.
##' @param min_node_size Minimal node size. Default 5.
##' @param replace Sample with replacement. Default TRUE.
##' @param unordered_factors How to handle unordered factor variables. Either "ignore" or "partition" with default "ignore".
##' @param num_threads Number of threads used for mclapply, set to 1 for debugging.
##' @examples
##' \donttest{
##' TODO: Include new example
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
  if (is.null(K) && is.null(Phi)) {
    stop("Either K or Phi must be specified.")
  }

  if (!is.null(K) && !is.null(Phi)) {
    warning("Both K and Phi specified. Only Phi is used.")
  }

  # TODO: Give informative error messages
  if (!is.null(Phi)) {
    stopifnot(
      is.matrix(Phi),
      nrow(Phi) > ncol(Phi)
    )
  } else {
    stopifnot(
      is.matrix(K),
      nrow(K) == ncol(K)
    )
  }
  stopifnot(
    length(tol) == 1L,
    is.finite(tol),
    tol >= 0
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
  covariate_levels <- list()

  if (unordered_factors == "ignore") {
    ## Just set to ordered if "ignore"
    character.idx <- sapply(X, is.character)
    ordered.idx <- sapply(X, is.ordered)
    factor.idx <- sapply(X, is.factor)
    recode.idx <- character.idx | (factor.idx & !ordered.idx)
    X[, recode.idx] <- lapply(X[, recode.idx], as.ordered)

    ## Save levels
    covariate_levels <- lapply(X, levels)
  }

  ## Create forest object
  forest <- ForestRegression$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry),
                                 min_node_size = as.integer(min_node_size),
                                 replace = replace, splitrule = splitrule,
                                 data = Data$new(data = model.data),
                                 formula = formula, unordered_factors = unordered_factors,
                                 covariate_levels = covariate_levels)

  ## Grow forest
  forest$grow(num_threads = num_threads)

  ## Return forest
  return(forest)
}
