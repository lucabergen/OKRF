
##' Compute approximate features used by the OKRF tree in splitting
approximateFeatures <- function(feat_rep, tol) {

  if (feat_rep$type == "explicit") {
    pqr_wrap(
      X = feat_rep$Phi,
      tol = tol
    )
  } else {
    pchol_wrap(
      K = feat_rep$K,
      tol = tol
    )
  }
}

##' Compute approximate features by pivoted Cholesky
##'
##' Uses LAPACKs DPSTRF
pchol_wrap <- function(K, tol) {

  C <- suppressWarnings(
    chol(K, pivot = TRUE, tol = tol)
  )

  r <- attr(C, "rank")
  p <- attr(C, "pivot")

  L <- matrix(0, nrow(K), r)

  if (r > 0) {
    # chol returns a n \times n matrix, so we use only the r relevant rows
    L[p, ] <- t(C[seq_len(r), , drop = FALSE])
  }

  L
}

##' Compute approximate features by pivoted QR used in splitting
##'
# Note that this uses DGEQP3, which does not stop at the specified tolerance, but
# computes the complete pivoted QR and then selects the lower rank corresponding
# to tol. This is only offered by DGEQP3RK, which can only be used via a Rcpp wrapper.
pqr_wrap <- function(X, tol) {

  Q <- qr(t(X), LAPACK = TRUE)
  R <- qr.R(Q)
  p <- Q$pivot

  pivots <- diag(R)^2
  rejected <- which(!is.finite(pivots) | pivots <= tol)

  r <- if (length(rejected)) {rejected[1] - 1} else {length(pivots)}

  L <- matrix(0, nrow(X), r)

  if (r > 0) {
    Rr <- R[seq_len(r), , drop = FALSE]

    signs <- sign(diag(R)[seq_len(r)])
    signs[signs == 0] <- 1
    Rr <- signs * Rr

    L[p, ] <- t(Rr)
  }

  L
}

##' Prepare features used by the OKRF tree
##'
##' `feat_rep` contains the original output representation, either `Phi`
##' or `K`. The returned object separates features used for splitting
##' from features used for terminal response predictions.
prepareFeatures <- function(feat_rep, tol, scope = "forest", reference_ids = NULL) {

  scope <- match.arg(scope, c("forest", "tree"))

  ## The current implementation is forest-level.
  ## Tree-level preparation will use reference_ids in a later step.
  if (scope == "tree") {
    stop(
      "Tree-level feature preparation is not implemented yet."
    )
  }

  split_features <- approximateFeatures(
    feat_rep = feat_rep,
    tol = tol
  )

  response_features <- if (identical(feat_rep$type, "explicit")) {
    feat_rep$Phi
  } else {
    matrix(
      numeric(0),
      nrow = 0L,
      ncol = 0L
    )
  }

  list(
    split_features = split_features,
    response_features = response_features,
    rank = as.integer(ncol(split_features)),
    scope = scope,
    reference_ids = reference_ids
  )
}
