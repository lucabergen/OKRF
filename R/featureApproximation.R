
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

##' Compute approximate features (only) for a tree-specific reference sample.
##'
##' The reference sample consists of the unique training observations
##' included in the tree's bootstrap sample.
##'
##' `split_features` has one row per training observation so that global
##' training IDs can be used for indexing. Approximate split features are
##' computed only for the reference observations. Zero rows for
##' non-reference observations are placeholders and do not affect split
##' scores because only reference observations are used during tree growth.
prepareTreeFeatures <- function(feat_rep, tol, reference_ids) {

  if (identical(feat_rep$type, "explicit")) {

    num_observations <- nrow(feat_rep$Phi)

    reference_rep <- list(
      type = "explicit",
      Phi = feat_rep$Phi[reference_ids,,drop = FALSE]
    )

  } else {

    num_observations <- nrow(feat_rep$K)

    # Create feature representation for reference sample
    reference_rep <- list(
      type = "implicit",
      K = feat_rep$K[reference_ids,reference_ids,drop = FALSE]
    )
  }


  ## Approximate features only on reference sample
  split_reference_features <- approximateFeatures(
    feat_rep = reference_rep,
    tol = tol
  )

  split_features <- matrix(
    0,
    nrow = num_observations,
    ncol = ncol(split_reference_features)
  )

  ##  split_features has nrow = n, with all-zero rows for non-reference
  ## observations and non-zero rows for reference observations

  if (ncol(split_reference_features) > 0L) {
    split_features[
      reference_ids, seq_len(ncol(split_features))
    ] <- split_reference_features
  }

  if (identical(feat_rep$type, "explicit")) {

    # Return the response features of the full training sample
    response_features <- feat_rep$Phi

  } else {

    # If K is given, no explicit feature representation is possible
    response_features <- matrix(numeric(0), nrow = 0L,ncol = 0L)

  }

  list(
    split_features = split_features,
    response_features = response_features,
    rank = as.integer(ncol(split_features)),
    scope = "tree",
    reference_ids = reference_ids
  )
}

##' Prepare features used by the OKRF tree
##'
##' `feat_rep` contains the original output representation, either `Phi`
##' or `K`. The returned object separates explicit features used for splitting
##' from explicit features used for leaf response predictions. The features used
##' for splitting are either approximated on the full sample (scope = "forest")
##' or on tree-specific reference samples only (scope = "tree"). The features
##' for response prediction are over the full sample in both cases.
##'
##' For explicit feature representations, `response_features` contains
##' the complete original feature matrix. For implicit kernel
##' representations, no explicit response feature matrix is available.
prepareFeatures <- function(feat_rep, tol, scope = "forest", reference_ids = NULL) {

  scope <- match.arg(scope, c("forest", "tree"))

  if (scope == "tree") {

    if (is.null(reference_ids) ||
        length(reference_ids) == 0L) {
      stop(
        "`reference_ids` must be provided for tree-level preparation."
      )
    }

    # Only use the unique IDs
    reference_ids <- sort(unique(as.integer(reference_ids)))

    prepareTreeFeatures(
      feat_rep = feat_rep,
      tol = tol,
      reference_ids = reference_ids
    )

  } else {

    split_features <- approximateFeatures(
      feat_rep = feat_rep,
      tol = tol
    )


    if (identical(feat_rep$type, "explicit")) {

      response_features <- feat_rep$Phi

    } else {

      # If K is given no explicit feature representation is possible
      response_features <- matrix(numeric(0), nrow = 0L, ncol = 0L)

    }

    list(
      split_features = split_features,
      response_features = response_features,
      rank = as.integer(ncol(split_features)),
      scope = scope,
      reference_ids = reference_ids
    )

  }

}
