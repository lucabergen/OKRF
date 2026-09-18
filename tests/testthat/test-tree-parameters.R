test_that("tree parameters are propagated to the forest and trees", {
  set.seed(20260917)

  n <- 20L

  X <- data.frame(
    x1 = seq_len(n),
    x2 = rev(seq_len(n))
  )

  Phi <- matrix(
    rnorm(n * 2L),
    nrow = n,
    ncol = 2L
  )

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    num_trees = 4L,
    mtry = 1L,
    min_node_size = 3L,
    max_depth = 2L,
    min_leaf_size = 4L,
    replace = FALSE,
    sample_fraction = 0.5,
    num_threads = 1L
  )

  expect_equal(forest$num_trees, 4L)
  expect_equal(forest$mtry, 1L)
  expect_equal(forest$min_node_size, 3L)
  expect_equal(forest$max_depth, 2L)
  expect_equal(forest$min_leaf_size, 4L)
  expect_false(forest$replace)
  expect_equal(forest$sample_fraction, 0.5)

  expect_length(forest$trees, 4L)

  expect_true(all(vapply(
    forest$trees,
    function(tree) {
      identical(tree$mtry, 1L) &&
        identical(tree$min_node_size, 3L) &&
        identical(tree$max_depth, 2L) &&
        identical(tree$min_leaf_size, 4L) &&
        identical(tree$sample_fraction, 0.5)
    },
    logical(1)
  )))
})


test_that("sample_fraction controls the sample size", {
  set.seed(20260917)

  n <- 20L

  X <- data.frame(
    x1 = seq_len(n),
    x2 = rnorm(n)
  )

  Phi <- matrix(
    rnorm(n * 2L),
    nrow = n,
    ncol = 2L
  )

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    num_trees = 5L,
    mtry = 1L,
    min_node_size = 2L,
    min_leaf_size = 2L,
    replace = FALSE,
    sample_fraction = 0.5,
    num_threads = 1L
  )

  expected_sample_size <- floor(0.5 * n)

  expect_true(all(vapply(
    forest$trees,
    function(tree) {
      length(tree$sampleIDs[[1L]]) == expected_sample_size
    },
    logical(1)
  )))
})


test_that("sampling without replacement produces no duplicate sample IDs", {
  set.seed(20260917)

  n <- 20L

  X <- data.frame(
    x1 = seq_len(n),
    x2 = rnorm(n)
  )

  Phi <- matrix(
    rnorm(n * 2L),
    nrow = n,
    ncol = 2L
  )

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    num_trees = 5L,
    mtry = 1L,
    min_node_size = 2L,
    min_leaf_size = 2L,
    replace = FALSE,
    sample_fraction = 0.5,
    num_threads = 1L
  )

  expect_true(all(vapply(
    forest$trees,
    function(tree) {
      sample_ids <- tree$sampleIDs[[1L]]

      length(unique(sample_ids)) == length(sample_ids)
    },
    logical(1)
  )))
})


test_that("max_depth equal to zero makes the root a leaf", {
  set.seed(20260917)

  n <- 20L

  X <- data.frame(
    x1 = seq_len(n),
    x2 = rnorm(n)
  )

  Phi <- matrix(
    rnorm(n * 2L),
    nrow = n,
    ncol = 2L
  )

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    num_trees = 4L,
    mtry = 1L,
    min_node_size = 2L,
    max_depth = 0L,
    min_leaf_size = 4L,
    replace = FALSE,
    sample_fraction = 1,
    num_threads = 1L
  )

  for (tree in forest$trees) {
    expect_length(
      tree$child_nodeIDs,
      length(tree$sampleIDs)
    )

    expect_null(
      tree$child_nodeIDs[[1L]]
    )

    expect_length(
      tree$leaf_train_ids[[1L]],
      n
    )
  }
})


test_that("all terminal nodes satisfy min_leaf_size", {
  set.seed(20260917)

  n <- 40L
  min_leaf_size <- 5L

  X <- data.frame(
    x1 = seq_len(n),
    x2 = rnorm(n)
  )

  Phi <- matrix(
    rnorm(n * 2L),
    nrow = n,
    ncol = 2L
  )

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    num_trees = 8L,
    mtry = 1L,
    min_node_size = 2L,
    min_leaf_size = min_leaf_size,
    replace = FALSE,
    sample_fraction = 1,
    num_threads = 1L
  )

  terminal_sizes <- unlist(
    lapply(
      forest$trees,
      function(tree) {
        terminal_nodes <- tree$leaf_train_ids[
          !vapply(
            tree$leaf_train_ids,
            is.null,
            logical(1)
          )
        ]

        vapply(
          terminal_nodes,
          length,
          integer(1)
        )
      }
    ),
    use.names = FALSE
  )

  expect_true(length(terminal_sizes) > 0L)
  expect_true(all(terminal_sizes >= min_leaf_size))
})


test_that("a split is only accepted when both children satisfy min_leaf_size", {
  set.seed(20260917)

  n <- 20L
  min_leaf_size <- 5L

  X <- data.frame(
    x = seq_len(n)
  )

  ## A one-dimensional signal aligned with x makes a useful
  ## improving split very likely and keeps the test deterministic.
  Phi <- matrix(
    X$x,
    nrow = n,
    ncol = 1L
  )

  forest <- simpleOKRF(
    Phi = Phi,
    X = X,
    num_trees = 1L,
    mtry = 1L,
    min_node_size = 2L,
    max_depth = 1L,
    min_leaf_size = min_leaf_size,
    replace = FALSE,
    sample_fraction = 1,
    num_threads = 1L
  )

  tree <- forest$trees[[1L]]

  ## max_depth = 1 allows the root split but makes its children terminal.
  expect_false(is.null(tree$child_nodeIDs[[1L]]))

  child_ids <- tree$child_nodeIDs[[1L]]

  child_sizes <- vapply(
    child_ids,
    function(child_id) {
      length(tree$leaf_train_ids[[child_id]])
    },
    integer(1)
  )

  expect_true(all(child_sizes >= min_leaf_size))
})


test_that(
  "a node smaller than two min_leaf_size values cannot be split",
  {
    set.seed(20260917)

    n <- 8L

    X <- data.frame(
      x = seq_len(n)
    )

    Phi <- matrix(
      X$x,
      nrow = n,
      ncol = 1L
    )

    forest <- simpleOKRF(
      Phi = Phi,
      X = X,
      num_trees = 1L,
      mtry = 1L,
      min_node_size = 2L,
      min_leaf_size = 5L,
      replace = FALSE,
      sample_fraction = 1,
      num_threads = 1L
    )

    tree <- forest$trees[[1L]]

    # Eight observations cannot produce two children
    # with at least five observations each.
    expect_length(
      tree$sampleIDs,
      1L
    )

    # Every node has an entry in child_nodeIDs.
    expect_length(
      tree$child_nodeIDs,
      length(tree$sampleIDs)
    )

    # The root is a leaf and therefore has no children.
    expect_null(
      tree$child_nodeIDs[[1L]]
    )

    expect_length(
      tree$leaf_train_ids[[1L]],
      n
    )
  }
)


test_that("invalid tree parameters are rejected", {
  X <- data.frame(
    x = seq_len(12L)
  )

  Phi <- matrix(
    rnorm(12L * 2L),
    nrow = 12L,
    ncol = 2L
  )

  expect_error(
    simpleOKRF(
      Phi = Phi,
      X = X,
      min_node_size = 0L
    ),
    "min_node_size"
  )

  expect_error(
    simpleOKRF(
      Phi = Phi,
      X = X,
      min_leaf_size = 0L
    ),
    "min_leaf_size"
  )

  expect_error(
    simpleOKRF(
      Phi = Phi,
      X = X,
      min_leaf_size = 1.5
    ),
    "min_leaf_size"
  )

  expect_error(
    simpleOKRF(
      Phi = Phi,
      X = X,
      max_depth = -1L
    ),
    "max_depth"
  )

  expect_error(
    simpleOKRF(
      Phi = Phi,
      X = X,
      max_depth = 1.5
    ),
    "max_depth"
  )

  expect_error(
    simpleOKRF(
      Phi = Phi,
      X = X,
      sample_fraction = 0
    ),
    "sample_fraction"
  )

  expect_error(
    simpleOKRF(
      Phi = Phi,
      X = X,
      sample_fraction = 1.1
    ),
    "sample_fraction"
  )

  expect_error(
    simpleOKRF(
      Phi = Phi,
      X = X,
      replace = NA
    ),
    "replace"
  )

  expect_error(
    simpleOKRF(
      Phi = Phi,
      X = X,
      num_threads = 0L
    ),
    "num_threads"
  )
})


test_that("explicit features are separated into split and response features", {
  Phi <- matrix(
    c(
      1, 0,
      0, 1,
      1, 1
    ),
    nrow = 3L,
    byrow = TRUE
  )

  prepared_features <- prepareFeatures(
    feat_rep = list(
      type = "explicit",
      Phi = Phi
    ),
    tol = 1e-8
  )

  expect_true(
    is.matrix(prepared_features$split_features)
  )

  expect_identical(
    prepared_features$response_features,
    Phi
  )

  expect_equal(
    nrow(prepared_features$split_features),
    nrow(Phi)
  )

  expect_equal(
    prepared_features$rank,
    ncol(prepared_features$split_features)
  )

  expect_identical(
    prepared_features$scope,
    "forest"
  )
})


test_that("kernel features do not create explicit response features", {
  K <- matrix(
    c(
      1, 0.2, 0.1,
      0.2, 1, 0.3,
      0.1, 0.3, 1
    ),
    nrow = 3L,
    byrow = TRUE
  )

  prepared_features <- prepareFeatures(
    feat_rep = list(
      type = "kernel",
      K = K
    ),
    tol = 1e-8
  )

  expect_true(
    is.matrix(prepared_features$split_features)
  )

  expect_equal(
    nrow(prepared_features$split_features),
    nrow(K)
  )

  expect_true(
    is.matrix(prepared_features$response_features)
  )

  expect_equal(
    nrow(prepared_features$response_features),
    0L
  )

  expect_equal(
    ncol(prepared_features$response_features),
    0L
  )
})


test_that("Tree receives prepared split and response features", {
  X <- data.frame(
    x = c(0, 1, 2, 3)
  )

  Phi <- matrix(
    c(
      1, 0,
      0, 1,
      1, 1,
      2, 1
    ),
    nrow = 4L,
    byrow = TRUE
  )

  forest <- simpleOKRF(
    X = X,
    Phi = Phi,
    num_trees = 1L
  )

  tree <- forest$trees[[1L]]

  expect_true(
    is.matrix(tree$prepared_features$split_features)
  )

  expect_identical(
    tree$prepared_features$response_features,
    Phi
  )

  expect_equal(
    nrow(tree$prepared_features$split_features),
    nrow(X)
  )
})


test_that("tree-level features use bootstrap reference IDs", {
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
    min_leaf_size = 1L
  )

  tree <- forest$trees[[1L]]

  expected_reference_ids <- sort(
    unique(
      tree$sampleIDs[[1L]]
    )
  )

  expect_identical(
    tree$prepared_features$scope,
    "tree"
  )

  expect_identical(
    tree$prepared_features$reference_ids,
    expected_reference_ids
  )

  expect_equal(
    nrow(tree$prepared_features$split_features),
    nrow(X)
  )

  expect_equal(
    tree$prepared_features$rank,
    ncol(tree$prepared_features$split_features)
  )
})


test_that("tree-level features work with an implicit kernel", {
  X <- data.frame(
    x = c(0, 1, 2, 3)
  )

  K <- outer(
    X$x,
    X$x,
    function(x, y) {
      exp(-(x - y)^2)
    }
  )

  forest <- simpleOKRF(
    X = X,
    K = K,
    num_trees = 1L,
    min_leaf_size = 1L
  )

  tree <- forest$trees[[1L]]

  expect_identical(
    tree$prepared_features$scope,
    "tree"
  )

  expect_equal(
    nrow(tree$prepared_features$split_features),
    nrow(X)
  )

  expect_equal(
    ncol(tree$prepared_features$response_features),
    0L
  )
})


test_that("tree-level reference IDs are valid", {
  X <- data.frame(
    x = seq_len(8L)
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

  reference_ids <- tree$prepared_features$reference_ids

  expect_true(
    all(reference_ids >= 1L)
  )

  expect_true(
    all(reference_ids <= nrow(X))
  )

  expect_equal(
    reference_ids,
    sort(unique(reference_ids))
  )
})


test_that("forest parameters are propagated to every tree", {
  X <- data.frame(
    x1 = seq_len(12L),
    x2 = rev(seq_len(12L))
  )

  Phi <- cbind(
    intercept = 1,
    x = X$x1
  )

  forest <- simpleOKRF(
    X = X,
    Phi = Phi,
    num_trees = 3L,
    mtry = 1L,
    min_node_size = 4L,
    min_leaf_size = 2L,
    max_depth = 2L
  )

  expect_length(
    forest$trees,
    3L
  )

  for (tree in forest$trees) {
    expect_equal(tree$mtry, forest$mtry)
    expect_equal(tree$min_node_size, forest$min_node_size)
    expect_equal(tree$min_leaf_size, forest$min_leaf_size)
    expect_equal(tree$max_depth, forest$max_depth)
    expect_equal(
      tree$sample_fraction,
      forest$sample_fraction
    )
  }
})


test_that("feature scores are computed correctly", {
  tree <- Tree$new(
    mtry = 1L,
    min_node_size = 1L,
    min_leaf_size = 1L,
    max_depth = 1L,
    unordered_factors = "ignore",
    x_data = Data$new(
      data = data.frame(x = 1:3)
    ),
    feat_rep = list(
      type = "explicit",
      Phi = diag(3L)
    ),
    tol = 1e-3,
    prepared_features = list(
      split_features = diag(3L),
      response_features = diag(3L),
      rank = 3L,
      scope = "forest",
      reference_ids = 1:3
    ),
    sample_fraction = 1
  )

  expect_equal(
    tree$computeFeatureScore(
      feature_sum = c(3, 4),
      num_observations = 2L
    ),
    12.5
  )
})
