
##' @title Forest class
##' @description Output kernel random forest.
##' @importFrom parallel mclapply
##' @import methods
Forest <- setRefClass("Forest",
  fields = list(
    num_trees = "integer",
    mtry = "integer",
    min_node_size = "integer",
    max_depth = "integer",
    min_leaf_size = "integer",
    unordered_factors = "character",
    feat_rep = "list",
    chol_features = "matrix",
    original_features = "matrix",
    tol = "numeric",
    approx_rank = "integer",
    x_data = "Data",
    trees = "list",
    replace = "logical",
    sample_fraction = "numeric",
    covariate_levels = "list"
  ),
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
      chol_features <<- approximateFeatures(feat_rep, tol)

      approx_rank <<- as.integer(ncol(chol_features))

      ## Validate dimensions
      if (nrow(chol_features) != x_data$nrow) {
        stop(
          "`chol_features` must have one row for each observation in `x_data`."
        )
      }

      covariate_levels <<- lapply(
        x_data$data,
        levels
      )

      invisible(.self)
    },


    grow = function(num_threads) {

      ## Allocate empty trees
      trees <<- replicate(
        num_trees,
        Tree$new(
          mtry = mtry,
          min_node_size = min_node_size,
          max_depth = max_depth,
          min_leaf_size = min_leaf_size,
          unordered_factors = unordered_factors,
          x_data = x_data,
          chol_features = chol_features,
          approx_rank = approx_rank,
          sample_fraction = sample_fraction
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


    prepareNewdata = function(newdata) {

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

      if (length(covariate_levels) != x_data$ncol) {
        stop(
          "`covariate_levels` must contain one entry for each covariate."
        )
      }

      ## Work on a copy. The input data.frame is not modified by reference.
      prepared_newdata <- newdata

      for (j in seq_len(x_data$ncol)) {

        training_levels <- covariate_levels[[j]]

        ## Numeric covariates do not have factor levels and are unchanged.
        if (is.null(training_levels)) {
          next
        }

        ## Convert factor, ordered factor, or character input to character
        ## before imposing the training levels.
        values <- as.character(
          prepared_newdata[[j]]
        )

        ## Check for levels that were not present during training.
        observed_values <- unique(
          values[!is.na(values)]
        )

        unknown_values <- setdiff(
          observed_values,
          training_levels
        )

        if (length(unknown_values) > 0L) {
          stop(
            "Column `",
            x_data$names[j],
            "` contains values not present in the training data: ",
            paste(unknown_values, collapse = ", ")
          )
        }

        ## Recreate the exact ordered-factor representation used for training.
        prepared_newdata[[j]] <- ordered(
          values,
          levels = training_levels
        )
      }

      prepared_newdata
    },


    predict = function(newdata, type = c("weights", "response")) {

      type <- match.arg(type)

      if (length(trees) == 0L) {
        stop("The forest has not been grown yet.")
      }

      newdata <- prepareNewdata(
        newdata = newdata
      )

      predict_data <- Data$new(data = newdata)

      num_newdata  <- predict_data$nrow
      num_training <- x_data$nrow

      ## Explicit Phi response prediction
      if (type == "response") {

        if (ncol(original_features) == 0L) {
          stop("`type = 'response'` requires explicit `Phi` features.")
        }

        forest_prediction <- matrix(
          0,
          nrow = num_newdata,
          ncol = ncol(original_features)
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


