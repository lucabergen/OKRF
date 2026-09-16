
##' @title Forest class
##' @description Output kernel random forest.
##' @importFrom parallel mclapply
##' @import methods
Forest <- setRefClass("Forest",
  fields = list(
    num_trees = "integer",
    mtry = "integer",
    min_node_size = "integer",
    unordered_factors = "character",
    feat_rep = "list",
    chol_features = "matrix",
    original_features = "matrix",
    tol = "numeric",
    approx_rank = "integer",
    x_data = "Data",
    trees = "list",
    replace = "logical",
    covariate_levels = "list"),
  methods = list(

    initialize = function(...) {

      callSuper(...)

      original_features <<- if (feat_rep$type == "explicit") {
        feat_rep$Phi
      } else {
        matrix(
          numeric(0),
          nrow = 0L,
          ncol = 0L
        )
      }

      ## Assign forest-specific Cholesky/QR features
      chol_features <- approximateFeatures(feat_rep, tol)

      ## Validate dimensions
      if (nrow(chol_features) != x_data$nrow) {
        stop(
          "`chol_features` must have one row for each observation in `x_data`."
        )
      }

      if (ncol(chol_features) != approx_rank) {
        stop("`approx_rank` must equal `ncol(chol_features)`.")
      }

      invisible(.self)
    },


    grow = function(num_threads) {

      ## Allocate empty trees
      trees <<- replicate(
        num_trees,
        Tree$new(
          mtry = mtry,
          min_node_size = min_node_size,
          unordered_factors = unordered_factors,
          x_data = x_data,
          chol_features = chol_features,
          approx_rank = approx_rank
        ),
        simplify = FALSE
      )

      ## Grow trees
      trees <<- parallel::mclapply(
        trees,
        function(tree) {
          tree$grow(replace = replace)
          tree
        },
        mc.cores = num_threads
      )

      ## Save leaf predictions for each tree
      if (ncol(original_features) > 0L) {
        for (tree in trees) {
          tree$setTerminalPredictions(original_features = original_features)
        }
      }

      invisible(.self)
    },


    predict = function(newdata, type = c("weights", "response")) {
      type <- match.arg(type)

      if (!is.data.frame(newdata)) {
        stop("`newdata` must be a data.frame.")
      }

      if (ncol(newdata) != x_data$ncol) {
        stop(
          "`newdata` must have the same number of columns as the training data."
        )
      }

      if (!identical(colnames(newdata), x_data$names)) {
        stop(
          "`newdata` must have the same column names and order as the training data."
        )
      }

      if (length(trees) == 0L) {
        stop("The forest has not been grown yet.")
      }

      predict_data <- Data$new(data = newdata)

      num_newdata  <- predict_data$nrow
      num_training <- x_data$nrow

      ## Explicit Phi response prediction
      if (type == "response") {

        if (ncol(phi_features) == 0L) {
          stop("`type = 'response'` requires explicit `Phi` features.")
        }

        forest_prediction <- matrix(
          0,
          nrow = num_newdata,
          ncol = ncol(phi_features)
        )

        for (tree in trees) {
          forest_prediction <- forest_prediction +
            tree$getTerminalPredictions(
              predict_data = predict_data
            )
        }

        forest_prediction / length(trees)

        ## Implicit weight prediction
      } else {

      forest_weights <- matrix(
        0,
        nrow = num_newdata,
        ncol = num_training
      )

      for (tree in trees) {

        terminal_sampleIDs_newdata <- tree$getTerminalSampleIDs(
          predict_data = predict_data
        )

        for (i in seq_len(num_newdata)) {

          terminal_sampleIDs <- terminal_sampleIDs_newdata[[i]]

          if (length(terminal_sampleIDs) == 0L) {
            next
          }

          sample_counts <- tabulate(
            terminal_sampleIDs,
            nbins = num_training
          )

          forest_weights[i, ] <- forest_weights[i, ] +
            sample_counts / length(terminal_sampleIDs)
        }
      }

      # Average over all trees
      forest_weights / length(trees)

      }
    },

    show = function() {
      cat("simpleOKRF Forest\n")
      cat("Number of trees:                 ", num_trees, "\n")
      cat("Sample size:                     ", x_data$nrow, "\n")
      cat("Number of independent variables: ", x_data$ncol, "\n")
      cat("Mtry:                            ", mtry, "\n")
      cat("Target node size:                ", min_node_size, "\n")
      cat("Replace                          ", replace, "\n")
      cat("Unordered factor handling        ", unordered_factors, "\n")
    },

    print = function() {
      show()
    })
)


