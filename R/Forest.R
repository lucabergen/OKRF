
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
    approx_scope = "character",
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

        approx_rank <<- as.integer(
          ncol(chol_features)
        )

      } else {

        stop(
          "K-based feature construction is not implemented yet."
        )
      }

      ## Validate dimensions
      if (nrow(chol_features) != x_data$nrow) {
        stop(
          "`chol_features` must have one row for each observation in `x_data`."
        )
      }

      if (ncol(chol_features) != approx_rank) {
        stop(
          "`approx_rank` must equal `ncol(chol_features)`."
        )
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
      # TODO: implement prediction; add checks if newdata has correct format
      stop("Prediction is not implemented yet.")
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


