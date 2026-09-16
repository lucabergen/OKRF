
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
    tol = "numeric",
    approx_rank = "integer",
    x_data = "Data",
    trees = "list",
    replace = "logical",
    covariate_levels = "list"),
  methods = list(

    initialize = function(...) {

      callSuper(...)

      # TODO: Add proper feature approximation

      ## Assign forest-specific Cholesky features
      if (feat_rep$type == "explicit") {

        chol_features <<- feat_rep$Phi

        approx_rank <<- as.integer(ncol(chol_features))

      } else {
        stop("K-based feature construction is not implemented yet.")
      }

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

      trees <<- parallel::mclapply(
        trees,
        function(tree) {
          tree$grow(replace = replace)
          tree
        },
        mc.cores = num_threads
      )

      invisible(.self)
    },

    predict = function(newdata, type = c("weights", "response")) {

      type <- match.arg(type)

      ## Only weights are implemented at this stage
      if (type != "weights") {
        stop("`type = 'response'` is not implemented yet. ")
      }

      ## Validate newdata
      if (!is.data.frame(newdata)) {
        stop("`newdata` must be a data.frame.")
      }

      if (ncol(newdata) != x_data$ncol) {
        stop(paste0("`newdata` must have ", x_data$ncol, " columns."))
      }

      if (!identical(colnames(newdata), x_data$names)) {
        stop(
          "`newdata` must have same column names and order as the training data."
        )
      }

      if (length(trees) == 0L) {
        stop("The forest has not been grown yet.")
      }

      ## Wrap newdata in the Data reference class
      predict_data <- Data$new(
        data = newdata
      )

      num_newdata <- predict_data$nrow
      num_training <- x_data$nrow
      num_trees_fitted <- length(trees)

      ## Final forest weight matrix:
      ## rows = observations in newdata,
      ## columns = observations in the training data
      forest_weights <- matrix(
        0,
        nrow = num_newdata,
        ncol = num_training
      )

      ## Process each tree
      for (tree in trees) {

        ## Get the terminal-node training IDs for all new observations
        terminal_sampleIDs_newdata <- tree$getTerminalSampleIDs(
          predict_data = predict_data
        )

        ## One result must exist for every row in newdata
        if (length(terminal_sampleIDs_newdata) != num_newdata) {
          stop(
            paste0("Number of terminal-node results (",
              length(terminal_sampleIDs_newdata),
              ") does not match the number of prediction rows (",
              num_newdata,
              ")."
            )
          )
        }

        ## Convert terminal-node memberships to tree weights
        for (i in seq_along(terminal_sampleIDs_newdata)) {

          ## `i` is the row index of the current observation in `newdata`.
          terminal_sampleIDs <- terminal_sampleIDs_newdata[[i]]

          if (length(terminal_sampleIDs) == 0L) {
            stop(
              paste0("Terminal node for prediction row ",i,
                " contains no training sample IDs."
              )
            )
          }

          ## Count training observations in the terminal node.
          ## Bootstrap duplicates are intentionally retained.
          sample_counts <- tabulate(
            terminal_sampleIDs,
            nbins = num_training
          )

          ## Normalize the weights within this tree
          forest_weights[i, ] <- forest_weights[i, ] +
            sample_counts / length(terminal_sampleIDs)
        }
      }

      ## Average over all trees
      forest_weights / num_trees_fitted
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


