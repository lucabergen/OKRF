
##' @title Forest class
##' @description Virtual class for output kernel random forest.
##' @importFrom parallel mclapply
##' @import methods
Forest <- setRefClass("Forest",
  fields = list(
    num_trees = "integer",
    mtry = "integer",
    min_node_size = "integer",
    unordered_factors = "character",
    feat_rep = "list",
    tol = "numeric",
    x_data = "Data",
    trees = "list",
    replace = "logical",
    covariate_levels = "list"),
  methods = list(

    grow = function(num_threads) {

      trees <<- replicate(
        num_trees,
        Tree$new(
          mtry = mtry,
          min_node_size = min_node_size,
          unordered_factors = unordered_factors,
          x_data = x_data,
          feat_rep = feat_rep,
          tol = tol
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
      cat("Type:                            ", treetype, "\n")
      cat("Number of trees:                 ", num_trees, "\n")
      cat("Sample size:                     ", x_data$nrow, "\n")
      cat("Number of independent variables: ", x_data$ncol, "\n")
      cat("Mtry:                            ", mtry, "\n")
      cat("Target node size:                ", min_node_size, "\n")
      cat("Replace                          ", replace, "\n")
      cat("Unordered factor handling        ", unordered_factors, "\n")
      cat("OOB prediction error:            ", predictionError(), "\n")
    },

    print = function() {
      show()
    })
)


