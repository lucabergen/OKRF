
##' @title Tree class
##' @description Output kernel random forest tree.
Tree <- setRefClass("Tree",
  fields = list(
    mtry = "integer",
    min_node_size = "integer",
    unordered_factors = "character",
    x_data = "Data",
    feat_rep = "list",
    tol = "numeric",
    approx_rank = "integer",
    sampleIDs = "list",
    oob_sampleIDs = "integer",
    child_nodeIDs = "list",
    split_varIDs = "integer",
    split_values = "numeric",
    split_levels_left = "list",
    terminal_sampleIDs = "list"),
  methods = list(

    grow = function(replace) {

      # Boostrap
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

      # Assign bootstrap samples to root node
      sampleIDs <<- list(bootstrap_sample)

      # Call recursive splitting function on root node
      splitNode(1)

      invisible(.self)
    },


    splitNode = function(nodeID) {
      ## Sample possible split variables
      possible_split_varIDs <- sample.int(x_data$ncol, mtry)

      ## Split node
      split <- splitNodeInternal(nodeID, possible_split_varIDs)

      if (!is.null(split)) {
        ## Assign split
        split_varIDs[[nodeID]] <<- split$varID
        split_values[[nodeID]] <<- split$value

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
