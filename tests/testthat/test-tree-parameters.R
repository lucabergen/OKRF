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


test_that("max_depth equal to zero makes the root terminal", {
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

  expect_true(all(vapply(
    forest$trees,
    function(tree) {
      length(tree$child_nodeIDs) == 0L
    },
    logical(1)
  )))

  expect_true(all(vapply(
    forest$trees,
    function(tree) {
      length(tree$terminal_sampleIDs[[1L]]) == n
    },
    logical(1)
  )))
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
        terminal_nodes <- tree$terminal_sampleIDs[
          !vapply(
            tree$terminal_sampleIDs,
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
      length(tree$terminal_sampleIDs[[child_id]])
    },
    integer(1)
  )

  expect_true(all(child_sizes >= min_leaf_size))
})


test_that("a node smaller than two min_leaf_size values cannot be split", {
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

  ## Eight observations cannot produce two children
  ## with at least five observations each.
  expect_length(
    tree$child_nodeIDs,
    0L
  )

  expect_equal(
    length(tree$terminal_sampleIDs[[1L]]),
    n
  )
})


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
