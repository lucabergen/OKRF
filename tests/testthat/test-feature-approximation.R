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
      type = "kernel",
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
