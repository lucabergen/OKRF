
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


# Uses LAPACKs DPSTRF
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

# # Test if these return the same result
# n = 5000
# x1 <- rnorm(n)
# x2 <- rnorm(n)
# y <- 0.7*x1 + 0.2*x2^2 - 0.3*x1*x2 + rnorm(n)
# # compute RFFs
# omega <- rnorm(100, mean = 0, sd = 1)
# projection <- outer(y, omega)
# Y <- sqrt(1 / 100) * cbind(cos(projection), sin(projection))
# K <- tcrossprod(Y)
#
# ch <- pchol_wrap(K, tol = 1e-4)
# qr <- pqr_wrap(Y, tol = 1e-4)
#
# dim(ch)
# dim(qr)
# all.equal(qr, ch)
