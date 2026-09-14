
##' @title Output kernel forest class
##' @description Subclass for output kernel forests.
##' Contains all fields and methods used special for output kernel forests.
ForestOK <- setRefClass("ForestOK",
  contains = "Forest",
  fields = list(
    feat_rep = "list",
    tol = "numeric"
  ),
  methods = list(

    grow = function(num_threads) {
      treetype <<- "OK"

      ## Create trees
      trees <<- replicate(
        num_trees,
        TreeOK$new(
          feat_rep = feat_rep,
          tol = tol
        ),
        simplify = FALSE
      )

      ## Call parent method
      callSuper(num_threads)
    },

    predict = function(newdata) {
      callSuper(newdata)
    },

    aggregatePredictions = function(predictions) {
      ## For all samples take mean over all trees
      return(rowMeans(predictions))
    },

    predictionError = function() {
      ## For each tree loop over OOB samples and get predictions
      tree_predictions <- sapply(trees, function(x) {
        oob_samples <- x$oob_sampleIDs
        result <- rep(NA, data$nrow)
        result[oob_samples] <- x$predictOOB()
        return(result)
      })

      ## Compute prediction for each sample
      predicted_response <- rowMeans(tree_predictions, na.rm = TRUE)

      ## Compare predictions with true data
      return(sum((predicted_response - data$column(1))^2, na.rm = TRUE) / data$nrow)
    })

)


