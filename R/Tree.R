
##' @title Tree class
##' @description Output kernel random forest tree.
Tree <- setRefClass("Tree",
  fields = list(
    mtry = "integer",
    min_node_size = "integer",
    unordered_factors = "character",
    x_data = "Data",
    chol_features = "matrix",
    approx_rank = "integer",
    sampleIDs = "list",
    oob_sampleIDs = "integer",
    child_nodeIDs = "list",
    split_varIDs = "integer",
    split_values = "numeric",
    split_levels_left = "list",
    terminal_sampleIDs = "list"),
  methods = list(

    initialize = function(...) {

      callSuper(...)

      if (!is.matrix(chol_features)) {
        stop(
          "`chol_features` must be a matrix."
        )
      }

      if (nrow(chol_features) != x_data$nrow) {
        stop(
          "`chol_features` must have one row for each observation in `x_data`."
        )
      }

      approx_rank <<- as.integer(
        ncol(chol_features)
      )

      invisible(.self)
    },

    grow = function(replace) {

      ## Validate Cholesky features
      if (nrow(chol_features) != x_data$nrow) {
        stop(
          "`chol_features` must have one row for each row in `x_data`."
        )
      }

      if (ncol(chol_features) != approx_rank) {
        stop(
          "`approx_rank` must equal `ncol(chol_features)`."
        )
      }


      ## Boostrap
      num_samples <- x_data$nrow

      if (replace) {
        num_bootstrap_samples <- num_samples
      } else {
        num_bootstrap_samples <- floor(0.6321 * num_samples)
      }

      bootstrap_sample <- sample.int(
        n = num_samples,
        size = num_bootstrap_samples,
        replace = replace
      )

      oob_sampleIDs <<- setdiff(
        seq_len(num_samples),
        unique(bootstrap_sample)
      )

      ## Assign bootstrap samples to root node
      sampleIDs <<- list(bootstrap_sample)

      # Call recursive splitting function on root node
      splitNode(1)

      invisible(.self)
    },

    ## Naive explicit split score, eg for unit tests
    # splitScore = function(left_ids, right_ids) {
    #
    #   left_sum <- colSums(chol_features[left_ids, , drop = FALSE])
    #   right_sum <- colSums(chol_features[right_ids, , drop = FALSE])
    #
    #   sum(left_sum^2) / length(left_ids) +
    #     sum(right_sum^2) / length(right_ids)
    # },

    splitNode = function(nodeID) {
      ## Sample possible split variables
      possible_split_varIDs <- sample.int(
        n = x_data$ncol,
        size = mtry,
        replace = FALSE
      )

      ## Split node
      split <- splitNodeInternal(
        nodeID = nodeID,
        possible_split_varIDs = possible_split_varIDs
      )

      if (!is.null(split)) {

        ## Save split information
        split_varIDs[[nodeID]] <<- split$varID
        split_values[[nodeID]] <<- split$value
        split_levels_left[[nodeID]] <<- split$values_left

        ## Create child nodes
        left_child <- length(sampleIDs) + 1
        right_child <- length(sampleIDs) + 2
        child_nodeIDs[[nodeID]] <<- c(left_child, right_child)

        ## For each sample in node, assign to left or right child
        if (length(split_levels_left[[nodeID]]) == 0) {
          ## Ordered splitting
          idx <- x_data$subset(sampleIDs[[nodeID]], split$varID) <= split$value
        } else {
          # Unordered splitting
          idx <- x_data$subset(sampleIDs[[nodeID]], split$varID) %in% split_levels_left[[nodeID]]
        }
        sampleIDs[[left_child]] <<- sampleIDs[[nodeID]][idx]
        sampleIDs[[right_child]] <<- sampleIDs[[nodeID]][!idx]

        ## Recursively call split node on child nodes
        splitNode(left_child)
        splitNode(right_child)
      } else {
        # Unlike in normal regression trees, leafs do not include predictions,
        # but observation indices (possibly with repetitions due to bootstrap)
        split_varIDs[[nodeID]] <<- NA_integer_
        makeTerminalNode(nodeID)
      }
    },

    splitNodeInternal = function(nodeID, possible_split_varIDs) {
      # TODO: Implement split function
      stop("Splitting has not been implemented yet")
    },

    makeTerminalNode = function(nodeID) {
      # Save observation indices
      terminal_sampleIDs[[nodeID]] <<- sampleIDs[[nodeID]]
    },

    )
)
