test_that("invalid X inputs are rejected", {
  Phi <- matrix(seq_len(12), nrow = 6L, ncol = 2L)

  expect_error(
    simpleOKRF(Phi = Phi, X = matrix(numeric(12), nrow = 6L, ncol = 2L)),
    "X must be a data.frame"
  )

  expect_error(
    simpleOKRF(
      Phi = matrix(seq_len(2), nrow = 2L, ncol = 1L),
      X = data.frame(x = 1)
    ),
    "at least two observations"
  )

  expect_error(
    simpleOKRF(
      Phi = matrix(seq_len(6), nrow = 3L, ncol = 2L),
      X = data.frame(row.names = seq_len(3L))
    ),
    "at least one covariate"
  )

  X_with_na <- data.frame(x = c(0, 1, NA, 3))
  Phi_with_na <- matrix(seq_len(8), nrow = 4L, ncol = 2L)

  expect_error(
    simpleOKRF(Phi = Phi_with_na, X = X_with_na),
    "Missing values in split covariates"
  )
})


test_that("feature representations are validated", {
  X <- data.frame(x = seq_len(6L))
  Phi <- matrix(seq_len(12), nrow = 6L, ncol = 2L)

  expect_error(
    simpleOKRF(X = X),
    "Either K or Phi must be specified"
  )

  expect_error(
    simpleOKRF(Phi = as.vector(Phi), X = X),
    "Phi must be a numeric matrix"
  )

  expect_error(
    simpleOKRF(Phi = matrix(seq_len(10), nrow = 5L, ncol = 2L), X = X),
    "same number of rows as X"
  )

  ## The implementation requires more rows than explicit feature columns.
  X_small <- data.frame(x = seq_len(4L))
  expect_error(
    simpleOKRF(
      Phi = matrix(seq_len(16), nrow = 4L, ncol = 4L),
      X = X_small
    ),
    "fewer columns than rows"
  )

  expect_error(
    simpleOKRF(K = matrix(seq_len(12), nrow = 3L, ncol = 4L), X = X),
    "K must be square"
  )

  nonsymmetric_K <- matrix(
    c(1, 2, 0,
      2, 1, 3,
      4, 3, 1),
    nrow = 3L,
    byrow = TRUE
  )

  expect_error(
    simpleOKRF(K = nonsymmetric_K, X = data.frame(x = 1:3)),
    "K must be symmetric"
  )

  expect_error(
    simpleOKRF(K = diag(3), X = X),
    "same number of rows and columns as X"
  )

  expect_warning(
    simpleOKRF(
      K = tcrossprod(Phi),
      Phi = Phi,
      X = X,
      num_trees = 1L,
      min_node_size = 10L
    ),
    "Both K and Phi specified"
  )
})


test_that("scalar numeric parameters are validated", {
  X <- data.frame(x = seq_len(8L))
  Phi <- matrix(seq_len(16), nrow = 8L, ncol = 2L)

  expect_error(
    simpleOKRF(Phi = Phi, X = X, tol = 0),
    "tol"
  )

  expect_error(
    simpleOKRF(Phi = Phi, X = X, tol = c(1e-3, 1e-4)),
    "tol"
  )

  expect_error(
    simpleOKRF(Phi = Phi, X = X, tol = NA_real_),
    "tol"
  )

  expect_error(
    simpleOKRF(Phi = Phi, X = X, mtry = 2L),
    "mtry must be an integer between 1 and the number of covariates."
  )

  expect_error(
    simpleOKRF(
      Phi = Phi,
      X = X,
      unordered_factors = "unsupported"
    ),
    "Unknown value for unordered_factors"
  )
})


test_that("a constant covariate produces valid predictions", {
  set.seed(20260917)

  n <- 12L
  X <- data.frame(x = rep(1, n))
  Phi <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    num_trees = 3L,
    mtry = 1L,
    min_node_size = 3L,
    num_threads = 1L
  )

  weights <- forest$predict(X[1:2, , drop = FALSE], type = "weights")
  response <- forest$predict(X[1:2, , drop = FALSE], type = "response")

  expect_equal(dim(weights), c(2L, n))
  expect_equal(dim(response), c(2L, ncol(Phi)))
  expect_true(all(is.finite(weights)))
  expect_true(all(weights >= 0))
  expect_equal(rowSums(weights), rep(1, 2L), tolerance = 1e-12)
  expect_equal(response, weights %*% Phi, tolerance = 1e-12)
})


test_that("sampling without replacement works for a small odd sample size", {
  set.seed(20260917)

  n <- 9L
  X <- data.frame(
    x1 = seq_len(n),
    x2 = rev(seq_len(n))
  )
  Phi <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    num_trees = 3L,
    mtry = 1L,
    min_node_size = 2L,
    replace = FALSE,
    num_threads = 1L
  )

  weights <- forest$predict(X[1:3, , drop = FALSE], type = "weights")

  expect_equal(dim(weights), c(3L, n))
  expect_true(all(is.finite(weights)))
  expect_true(all(weights >= 0))
  expect_equal(rowSums(weights), rep(1, 3L), tolerance = 1e-12)
})


test_that("prediction validates newdata type, columns, and names", {
  set.seed(20260917)

  n <- 14L
  X <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n)
  )
  Phi <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    num_trees = 3L,
    mtry = 1L,
    min_node_size = 3L,
    num_threads = 1L
  )

  expect_error(
    forest$predict(as.matrix(X[1:2, ]), type = "weights"),
    "newdata.*data.frame"
  )

  expect_error(
    forest$predict(X[1:2, "x1", drop = FALSE], type = "weights"),
    "same number of columns"
  )

  renamed <- X[1:2, , drop = FALSE]
  names(renamed) <- c("wrong", "x2")

  expect_error(
    forest$predict(renamed, type = "weights"),
    "same column names and order"
  )

  reversed <- X[1:2, c("x2", "x1"), drop = FALSE]

  expect_error(
    forest$predict(reversed, type = "weights"),
    "same column names and order"
  )

  expect_error(
    forest$predict(X[1:2, , drop = FALSE], type = "invalid"),
    "one of"
  )
})


test_that("empty newdata returns correctly shaped empty predictions", {
  set.seed(20260917)

  n <- 12L
  X <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n)
  )
  Phi <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    num_trees = 3L,
    mtry = 1L,
    min_node_size = 3L,
    num_threads = 1L
  )

  empty_newdata <- X[integer(0), , drop = FALSE]
  weights <- forest$predict(empty_newdata, type = "weights")
  response <- forest$predict(empty_newdata, type = "response")

  expect_equal(dim(weights), c(0L, n))
  expect_equal(dim(response), c(0L, ncol(Phi)))
})


test_that("unknown factor levels are rejected after character recoding", {
  set.seed(20260917)

  n <- 15L
  X <- data.frame(
    x1 = rnorm(n),
    group = factor(
      rep(c("A", "B", "C"), length.out = n),
      levels = c("A", "B", "C")
    )
  )
  Phi <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    num_trees = 3L,
    mtry = 1L,
    min_node_size = 3L,
    unordered_factors = "ignore",
    num_threads = 1L
  )

  bad_newdata <- X[1:3, , drop = FALSE]
  bad_newdata$group <- as.character(bad_newdata$group)
  bad_newdata$group[1L] <- "D"

  expect_error(
    forest$predict(bad_newdata, type = "weights"),
    "values not present in the training data"
  )
})


test_that("a root-only tree routes every observation to the root", {

  X <- data.frame(
    x = c(0, 1, 2, 3, 4, 5)
  )

  Phi <- cbind(
    intercept = 1,
    x = X$x
  )

  forest <- simpleOKRF(
    X = X,
    Phi = Phi,
    num_trees = 1L,
    max_depth = 0L
  )

  tree <- forest$trees[[1L]]

  leaf_ids <- tree$findLeafIDs(
    predict_data = Data$new(
      data = X[1:3, , drop = FALSE]
    )
  )

  expect_equal(
    leaf_ids,
    rep(1L, 3L)
  )
})


test_that("a constant covariate routes observations to valid leaves", {

  X <- data.frame(
    x = c(0, 1, 2, 3, 4, 5)
  )

  Phi <- cbind(
    intercept = 1,
    x = X$x
  )

  X_constant <- X

  X_constant[, 1L] <- X_constant[1L, 1L]

  forest <- simpleOKRF(
    X = X_constant,
    Phi = Phi,
    num_trees = 1L
  )

  tree <- forest$trees[[1L]]

  leaf_ids <- tree$findLeafIDs(
    predict_data = Data$new(
      data = X_constant[1:3, , drop = FALSE]
    )
  )

  expect_true(
    all(leaf_ids >= 1L)
  )

  expect_true(
    all(leaf_ids <= length(tree$sampleIDs))
  )
})


test_that("every reached leaf contains training IDs", {

  X <- data.frame(
    x = c(0, 1, 2, 3, 4, 5)
  )

  Phi <- cbind(
    intercept = 1,
    x = X$x
  )

  forest <- simpleOKRF(
    X = X,
    Phi = Phi,
    num_trees = 1L
  )

  tree <- forest$trees[[1L]]

  leaf_ids <- tree$findLeafIDs(
    predict_data = Data$new(
      data = X[1:5, , drop = FALSE]
    )
  )

  train_ids <- tree$getTrainIDsByLeaf(
    leaf_ids = leaf_ids
  )

  expect_length(
    train_ids,
    5L
  )

  expect_true(
    all(lengths(train_ids) > 0L)
  )
})


test_that("child node storage is normalized", {

  X <- data.frame(
    x = c(0, 1, 2, 3, 4, 5)
  )

  Phi <- cbind(
    intercept = 1,
    x = X$x
  )

  forest <- simpleOKRF(
    X = X,
    Phi = Phi,
    num_trees = 1L
  )

  tree <- forest$trees[[1L]]

  expect_length(
    tree$child_nodeIDs,
    length(tree$sampleIDs)
  )

  invalid_child_nodes <- vapply(
    tree$child_nodeIDs,
    function(child_ids) {
      if (is.null(child_ids)) {
        return(FALSE)
      }

      length(child_ids) != 2L ||
        anyNA(child_ids) ||
        any(child_ids < 1L) ||
        any(child_ids > length(tree$sampleIDs))
    },
    logical(1)
  )

  expect_false(
    any(invalid_child_nodes)
  )
})


test_that("leaf IDs are valid integer node IDs", {

  X <- data.frame(
    x = c(0, 1, 2, 3, 4, 5)
  )

  Phi <- cbind(
    intercept = 1,
    x = X$x
  )

  forest <- simpleOKRF(
    X = X,
    Phi = Phi,
    num_trees = 1L
  )

  tree <- forest$trees[[1L]]

  predict_data <- Data$new(
    data = X[1:5, , drop = FALSE]
  )

  leaf_ids <- tree$findLeafIDs(
    predict_data = predict_data
  )

  expect_type(
    leaf_ids,
    "integer"
  )

  expect_true(
    all(leaf_ids >= 1L)
  )

  expect_true(
    all(leaf_ids <= length(tree$sampleIDs))
  )
})


