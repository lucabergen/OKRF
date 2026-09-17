test_that("explicit Phi response agrees with the weight reference", {
  set.seed(20260916)

  n <- 40L
  num_features <- 4L

  X <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n),
    x3 = sample(seq_len(5L), n, replace = TRUE)
  )

  Phi <- matrix(
    rnorm(n * num_features),
    nrow = n,
    ncol = num_features
  )

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    tol = 1e-6,
    num_trees = 8L,
    mtry = 2L,
    min_node_size = 4L,
    replace = TRUE,
    unordered_factors = "ignore",
    num_threads = 1L
  )

  newdata <- X[1:10, , drop = FALSE]

  weights <- forest$predict(
    newdata = newdata,
    type = "weights"
  )

  response <- forest$predict(
    newdata = newdata,
    type = "response"
  )

  expect_equal(
    dim(weights),
    c(nrow(newdata), n)
  )

  expect_equal(
    dim(response),
    c(nrow(newdata), num_features)
  )

  expect_equal(
    rowSums(weights),
    rep(1, nrow(newdata)),
    tolerance = 1e-10
  )

  expect_equal(
    response,
    weights %*% Phi,
    tolerance = 1e-10
  )
})


test_that("factor and character newdata give the same response", {
  set.seed(20260916)

  n <- 40L
  num_features <- 3L

  X <- data.frame(
    x1 = rnorm(n),
    x2 = factor(
      sample(c("A", "B", "C"), n, replace = TRUE),
      levels = c("A", "B", "C")
    )
  )

  Phi <- matrix(
    rnorm(n * num_features),
    nrow = n,
    ncol = num_features
  )

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    tol = 1e-6,
    num_trees = 5L,
    mtry = 1L,
    min_node_size = 4L,
    replace = TRUE,
    unordered_factors = "ignore",
    num_threads = 1L
  )

  factor_newdata <- X[1:8, , drop = FALSE]

  character_newdata <- factor_newdata
  character_newdata$x2 <- as.character(character_newdata$x2)

  factor_response <- forest$predict(
    newdata = factor_newdata,
    type = "response"
  )

  character_response <- forest$predict(
    newdata = character_newdata,
    type = "response"
  )

  expect_equal(
    character_response,
    factor_response,
    tolerance = 1e-10
  )
})


test_that("unknown factor levels are rejected", {
  set.seed(20260916)

  n <- 40L
  num_features <- 3L

  X <- data.frame(
    x1 = rnorm(n),
    x2 = factor(
      sample(c("A", "B", "C"), n, replace = TRUE),
      levels = c("A", "B", "C")
    )
  )

  Phi <- matrix(
    rnorm(n * num_features),
    nrow = n,
    ncol = num_features
  )

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    tol = 1e-6,
    num_trees = 5L,
    mtry = 1L,
    min_node_size = 4L,
    replace = TRUE,
    unordered_factors = "ignore",
    num_threads = 1L
  )

  bad_newdata <- X[1:8, , drop = FALSE]
  bad_newdata$x2 <- as.character(bad_newdata$x2)
  bad_newdata$x2[1] <- "unknown-level"

  expect_error(
    forest$predict(
      newdata = bad_newdata,
      type = "response"
    ),
    regexp = "not present in the training data"
  )
})


test_that("K representation supports weights", {
  set.seed(20260916)

  n <- 30L
  num_features <- 3L

  X <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n)
  )

  Phi <- matrix(
    rnorm(n * num_features),
    nrow = n,
    ncol = num_features
  )

  K <- tcrossprod(Phi)

  forest <- simpleOKRF(
    K = K,
    X = X,
    tol = 1e-6,
    num_trees = 4L,
    mtry = 1L,
    min_node_size = 4L,
    replace = TRUE,
    unordered_factors = "ignore",
    num_threads = 1L
  )

  newdata <- X[1:5, , drop = FALSE]

  weights <- forest$predict(
    newdata = newdata,
    type = "weights"
  )

  expect_equal(
    dim(weights),
    c(nrow(newdata), n)
  )

  expect_equal(
    rowSums(weights),
    rep(1, nrow(newdata)),
    tolerance = 1e-10
  )
})


test_that("K representation rejects explicit response prediction", {
  set.seed(20260916)

  n <- 30L
  num_features <- 3L

  X <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n)
  )

  Phi <- matrix(
    rnorm(n * num_features),
    nrow = n,
    ncol = num_features
  )

  forest <- simpleOKRF(
    K = tcrossprod(Phi),
    X = X,
    tol = 1e-6,
    num_trees = 4L,
    mtry = 1L,
    min_node_size = 4L,
    replace = TRUE,
    unordered_factors = "ignore",
    num_threads = 1L
  )

  newdata <- X[1:5, , drop = FALSE]

  expect_error(
    forest$predict(
      newdata = newdata,
      type = "response"
    ),
    regexp = "requires explicit.*Phi"
  )
})
