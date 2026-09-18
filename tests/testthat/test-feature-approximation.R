test_that("QR and pivoted Cholesky preserve the same Gram matrix", {
  set.seed(20260916)

  n <- 50L
  num_frequencies <- 20L
  num_covariates <- 3L

  ## Random Fourier feature construction.
  Z <- matrix(
    rnorm(n * num_covariates),
    nrow = n,
    ncol = num_covariates
  )

  omega <- matrix(
    rnorm(num_covariates * num_frequencies),
    nrow = num_covariates,
    ncol = num_frequencies
  )

  projection <- Z %*% omega

  Phi <- sqrt(1 / num_frequencies) * cbind(
    cos(projection),
    sin(projection)
  )

  K <- tcrossprod(Phi)

  qr_features <- simpleOKRF:::pqr_wrap(
    X = Phi,
    tol = 1e-10
  )

  cholesky_features <- simpleOKRF:::pchol_wrap(
    K = K,
    tol = 1e-10
  )

  ## Both factorizations should retain the same numerical rank.
  expect_equal(
    ncol(qr_features),
    ncol(cholesky_features)
  )

  ## The relevant invariant is the represented Gram matrix.
  expect_equal(
    tcrossprod(qr_features),
    K,
    tolerance = 1e-8
  )

  expect_equal(
    tcrossprod(cholesky_features),
    K,
    tolerance = 1e-8
  )

  expect_equal(
    tcrossprod(qr_features),
    tcrossprod(cholesky_features),
    tolerance = 1e-8
  )
})


test_that("approximateFeatures dispatches on the representation type", {
  set.seed(20260916)

  n <- 50L
  num_features <- 8L

  Phi <- matrix(
    rnorm(n * num_features),
    nrow = n,
    ncol = num_features
  )

  K <- tcrossprod(Phi)

  explicit_features <- simpleOKRF:::approximateFeatures(
    feat_rep = list(
      type = "explicit",
      Phi = Phi
    ),
    tol = 1e-10
  )

  kernel_features <- simpleOKRF:::approximateFeatures(
    feat_rep = list(
      type = "implicit",
      K = K
    ),
    tol = 1e-10
  )

  expect_equal(
    tcrossprod(explicit_features),
    K,
    tolerance = 1e-8
  )

  expect_equal(
    tcrossprod(kernel_features),
    K,
    tolerance = 1e-8
  )
})


test_that("tree-level split features match the reference approximation", {
  X <- data.frame(
    x = seq_len(10L)
  )

  Phi <- cbind(
    intercept = 1,
    x = X$x
  )

  forest <- simpleOKRF(
    X = X,
    Phi = Phi,
    num_trees = 1L,
    min_leaf_size = 2L
  )

  tree <- forest$trees[[1L]]

  reference_ids <- tree$prepared_features$reference_ids

  expected_features <- approximateFeatures(
    feat_rep = list(
      type = "explicit",
      Phi = Phi[
        reference_ids,
        ,
        drop = FALSE
      ]
    ),
    tol = forest$tol
  )

  actual_features <- tree$prepared_features$split_features[
    reference_ids,
    ,
    drop = FALSE
  ]

  expect_equal(
    actual_features,
    expected_features
  )
})


test_that("all tree training IDs are covered by the reference IDs", {
  X <- data.frame(
    x = seq_len(10L)
  )

  Phi <- cbind(
    intercept = 1,
    x = X$x
  )

  forest <- simpleOKRF(
    X = X,
    Phi = Phi,
    num_trees = 1L,
    min_leaf_size = 2L
  )

  tree <- forest$trees[[1L]]

  reference_ids <- tree$prepared_features$reference_ids

  used_ids <- sort(
    unique(
      unlist(tree$sampleIDs)
    )
  )

  expect_true(
    all(used_ids %in% reference_ids)
  )
})


test_that("every accepted split improves the split score", {

  split_score <- function(features, sample_ids) {
    node_features <- features[
      sample_ids,
      ,
      drop = FALSE
    ]

    node_sum <- colSums(node_features)

    sum(node_sum^2) / length(sample_ids)
  }

  X <- data.frame(
    x = seq_len(12L)
  )

  Phi <- cbind(
    intercept = 1,
    x = X$x
  )

  forest <- simpleOKRF(
    X = X,
    Phi = Phi,
    num_trees = 1L,
    min_leaf_size = 2L
  )

  tree <- forest$trees[[1L]]
  features <- tree$prepared_features$split_features

  split_nodes <- which(
    vapply(
      tree$child_nodeIDs,
      function(child_ids) {
        length(child_ids) == 2L
      },
      logical(1)
    )
  )

  expect_true(
    length(split_nodes) > 0L
  )

  for (node_id in split_nodes) {

    parent_ids <- tree$sampleIDs[[node_id]]
    child_ids <- tree$child_nodeIDs[[node_id]]

    left_ids <- tree$sampleIDs[[child_ids[[1L]]]]
    right_ids <- tree$sampleIDs[[child_ids[[2L]]]]

    parent_score <- split_score(
      features = features,
      sample_ids = parent_ids
    )

    child_score <- split_score(
      features = features,
      sample_ids = left_ids
    ) +
      split_score(
        features = features,
        sample_ids = right_ids
      )

    expect_gte(
      length(left_ids),
      tree$min_leaf_size
    )

    expect_gte(
      length(right_ids),
      tree$min_leaf_size
    )

    expect_gt(
      child_score,
      parent_score
    )
  }
})


test_that("response predictions equal weighted response features", {

  X <- data.frame(
      x = seq_len(10L)
    )

    Phi <- cbind(
      intercept = 1,
      x = X$x
    )

    forest <- simpleOKRF(
      X = X,
      Phi = Phi,
      num_trees = 2L,
      min_leaf_size = 2L
    )

    response_prediction <- forest$predict(
      newdata = X,
      type = "response"
    )

    weights <- forest$predict(
      newdata = X,
      type = "weights"
    )

    expected_prediction <- weights %*% Phi

    expect_equal(
      dim(response_prediction),
      dim(expected_prediction)
    )

    expect_equal(
      unname(response_prediction),
      unname(expected_prediction),
      tolerance = 1e-8
    )
  }
)


test_that("trees use tree-level split features", {
  X <- data.frame(
    x = seq_len(10L)
  )

  Phi <- cbind(
    intercept = 1,
    x = X$x
  )

  forest <- simpleOKRF(
    X = X,
    Phi = Phi,
    num_trees = 3L,
    min_leaf_size = 2L
  )

  expect_true(
    all(
      vapply(
        forest$trees,
        function(tree) {
          identical(
            tree$prepared_features$scope,
            "tree"
          )
        },
        logical(1)
      )
    )
  )
})


test_that("tree reference IDs match the bootstrap sample", {
  X <- data.frame(
    x = seq_len(10L)
  )

  Phi <- cbind(
    intercept = 1,
    x = X$x
  )

  forest <- simpleOKRF(
    X = X,
    Phi = Phi,
    num_trees = 1L,
    min_leaf_size = 2L
  )

  tree <- forest$trees[[1L]]

  expect_identical(
    tree$prepared_features$reference_ids,
    sort(unique(tree$sampleIDs[[1L]]))
  )
})
