
##' @title Tree class
##' @description Output kernel random forest tree.
Tree <- setRefClass("Tree",
  fields = list(
    mtry = "integer",
    min_node_size = "integer",
    max_depth = "integer",
    min_leaf_size = "integer",
    unordered_factors = "character",
    x_data = "Data",
    prepared_features = "list",
    sampleIDs = "list",
    oob_sampleIDs = "integer",
    child_nodeIDs = "list",
    split_varIDs = "integer",
    split_values = "numeric",
    split_levels_left = "list",
    leaf_train_ids = "list",
    leaf_predictions = "matrix",
    sample_fraction = "numeric"
  ),
  methods = list(

    initialize = function(...) {

      callSuper(...)

      if (!is.list(prepared_features)) {
        stop("`prepared_features` must be a list.")
      }

      if (!is.matrix(prepared_features$split_features)) {
        stop(
          "`prepared_features$split_features` must be a matrix."
        )
      }

      if (nrow(prepared_features$split_features) != x_data$nrow) {
        stop(
          "`prepared_features$split_features` must have one row for each ",
          "observation in `x_data`."
        )
      }

      leaf_predictions <<- matrix(
        nrow = 0L,
        ncol = 0L
      )

      invisible(.self)
    },

    grow = function(replace) {

      ## Validate Cholesky features
      if (!is.list(prepared_features)) {
        stop("`prepared_features` must be a list.")
      }

      if (!is.matrix(prepared_features$split_features)) {
        stop("`prepared_features$split_features` must be a matrix.")
      }

      if (nrow(prepared_features$split_features) != x_data$nrow) {
        stop(
          "`prepared_features$split_features` must have one row for each ",
          "row in `x_data`."
        )
      }

      if (length(prepared_features$rank) != 1L ||
          !is.numeric(prepared_features$rank) ||
          !is.finite(prepared_features$rank) ||
          prepared_features$rank != ncol(prepared_features$split_features)) {
        stop(
          "`prepared_features$rank` must equal ",
          "`ncol(prepared_features$split_features)`."
        )
      }

      ## Boostrap
      num_samples <- x_data$nrow

      num_bootstrap_samples <- floor(
        sample_fraction * num_samples
      )

      if (!replace && num_bootstrap_samples > num_samples) {
        stop(
          "sample_fraction cannot produce more samples than available ",
          "when replace is FALSE."
        )
      }

      if (num_bootstrap_samples < 1L) {
        stop(
          "sample_fraction must select at least one observation."
        )
      }

      if (num_bootstrap_samples < min_leaf_size) {
        stop(
          "The sampled training data contains fewer observations ",
          "than min_leaf_size."
        )
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

    splitNode = function(nodeID, depth = 0L) {

      node_sampleIDs <- sampleIDs[[nodeID]]

      ## A terminal node must contain at least min_leaf_size observations
      if (length(node_sampleIDs) < min_leaf_size) {
        split_varIDs[[nodeID]] <<- NA_integer_
        makeTerminalNode(nodeID)
        return(invisible(NULL))
      }

      ## A binary split requires two children of at least
      ## min_leaf_size observations each
      if (length(node_sampleIDs) < 2L * min_leaf_size) {
        split_varIDs[[nodeID]] <<- NA_integer_
        makeTerminalNode(nodeID)
        return(invisible(NULL))
      }

      ## Stop at the maximum depth
      if (!is.na(max_depth) && depth >= max_depth) {
        split_varIDs[[nodeID]] <<- NA_integer_
        makeTerminalNode(nodeID)
        return(invisible(NULL))
      }

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
          ## Numeric or ordered splitting
          idx <- as.numeric(
            x_data$subset(sampleIDs[[nodeID]], split$varID)
          ) <= split$value
        } else {
          # Unordered splitting
          idx <- x_data$subset(sampleIDs[[nodeID]], split$varID) %in% split_levels_left[[nodeID]]
        }
        sampleIDs[[left_child]] <<- sampleIDs[[nodeID]][idx]
        sampleIDs[[right_child]] <<- sampleIDs[[nodeID]][!idx]

        ## Recursively call split node on child nodes
        splitNode(
          nodeID = left_child,
          depth = depth + 1L
        )

        splitNode(
          nodeID = right_child,
          depth = depth + 1L
        )

      } else {
        # Unlike in normal regression trees, leafs do not include predictions,
        # but observation indices (possibly with repetitions due to bootstrap)
        split_varIDs[[nodeID]] <<- NA_integer_
        makeTerminalNode(nodeID)
      }
    },

    findBestSplit = function(nodeID, possible_split_varIDs) {

      best_split <- list(
        score = -Inf,
        varID = NA_integer_,
        value = NA_real_,
        values_left = character(0)
      )

      for (split_varID in possible_split_varIDs) {

        data_values <- x_data$subset(sampleIDs[[nodeID]], split_varID)

        if (anyNA(data_values)) {
          stop("Missing values in split covariates are not supported yet.")
        }

        if (is.numeric(data_values) || is.ordered(data_values)) {

          best_split <- findBestSplitValueOrdered(
            nodeID = nodeID,
            split_varID = split_varID,
            best_split = best_split
          )

        } else {
          stop("Unordered factors are not implemented yet.")
        }

      }

      if (!is.finite(best_split$score)) {
        return(NULL)
      }

      best_split
    },

    findBestSplitValueOrdered = function(nodeID, split_varID, best_split) {

      node_sampleIDs <- sampleIDs[[nodeID]]

      data_values <- x_data$subset(node_sampleIDs, split_varID)

      ## Numeric and ordered variables are represented by their numeric codes
      ordered_values <- as.numeric(data_values)

      ## Sort only once
      order_idx <- order(ordered_values)

      ordered_values <- ordered_values[order_idx]

      ordered_split_features <- prepared_features$split_features[
        node_sampleIDs[order_idx],,drop = FALSE
      ]

      num_samples <- length(ordered_values)

      ## End positions of groups with equal covariate values
      group_ends <- which(
        ordered_values[-num_samples] < ordered_values[-1L]
      )

      if (length(group_ends) == 0L) {
        return(best_split)
      }

      ## Initially all observations are in the right child
      sum_left <- numeric(length = prepared_features$rank)
      sum_right <- colSums(ordered_split_features)

      n_left <- 0L
      n_right <- num_samples
      group_start <- 1L

      for (group_end in group_ends) {

        ## Move one complete value group from right to left
        group_features <- ordered_split_features[group_start:group_end,,
                                                drop = FALSE]

        group_sum <- colSums(group_features)
        group_size <- group_end - group_start + 1L

        sum_left <- sum_left + group_sum
        sum_right <- sum_right - group_sum

        n_left <- n_left + group_size
        n_right <- n_right - group_size

        ## Only admissible leaf node sizes are considered
        if (n_left >= min_leaf_size && n_right >= min_leaf_size){

          score <- sum(sum_left^2) / n_left + sum(sum_right^2) / n_right

          if (score > best_split$score) {

            split_value <- ordered_values[group_end] +
              (ordered_values[group_end + 1L] - ordered_values[group_end]) / 2

            best_split <- list(
              score = score,
              varID = as.integer(split_varID),
              value = split_value,
              values_left = character(0)
            )
          }
        }

        ## Continue with the next group
        group_start <- group_end + 1L
      }

      best_split
    },


    splitNodeInternal = function(nodeID, possible_split_varIDs) {

      node_sampleIDs <- sampleIDs[[nodeID]]

      ## Do not attempt a split below min_node_size
      if (length(node_sampleIDs) < min_node_size) {
        return(NULL)
      }

      ## A binary split gives two children with at least min_leaf_size obs. each
      if (length(node_sampleIDs) < 2L * min_leaf_size) {
        return(NULL)
      }


      ## Score of the unsplit node
      split_features <- prepared_features$split_features

      node_features <- split_features[node_sampleIDs,,drop = FALSE]

      node_sum <- colSums(node_features)

      node_score <- sum(node_sum^2) / length(node_sampleIDs)

      ## Search candidate splits
      best_split <- findBestSplit(
        nodeID = nodeID,
        possible_split_varIDs = possible_split_varIDs
      )

      if (is.null(best_split)) {
        return(NULL)
      }

      ## Only accept an improvement
      if (best_split$score <= node_score) {
        return(NULL)
      }

      best_split
    },

    # Find the ID of the leaf containing the observation
    findLeafID = function(predict_data, row_id) {
      nodeID <- 1L

      while (TRUE) {

        if (nodeID > length(child_nodeIDs)) {
          return(as.integer(nodeID))
        }

        child_node_ids <- child_nodeIDs[[nodeID]]

        # Stop when the current node is a leaf.
        if (length(child_node_ids) == 0L ||
            is.null(child_node_ids)) {
          return(as.integer(nodeID))
        }

        split_var <- split_varIDs[nodeID]
        value <- predict_data$subset(row_id, split_var)

        # Follow the left or right child for ordered splits.
        if (length(split_levels_left[[nodeID]]) == 0L) {
          if (as.numeric(value) <= split_values[nodeID]) {
            nodeID <- child_node_ids[[1L]]
          } else {
            nodeID <- child_node_ids[[2L]]
          }

          # Follow the left or right child for factor splits.
        } else {
          if (value %in% split_levels_left[[nodeID]]) {
            nodeID <- child_node_ids[[1L]]
          } else {
            nodeID <- child_node_ids[[2L]]
          }
        }
      }
    },

    # Find the ID of the leaf reached for every observation
    findLeafIDs = function(predict_data) {
      vapply(
        seq_len(predict_data$nrow),
        function(row_id) {
          findLeafID(
            predict_data = predict_data,
            row_id = row_id
          )
        },
        integer(1)
      )
    },

    ## TODO: Remove wrapper from code
    getTerminalNodeIDs = function(predict_data) {

      findLeafIDs(predict_data)

    },

    # Get the IDs of the training observation contained in a leaf
    getTrainIDsByLeaf = function(leaf_ids) {
      lapply(leaf_ids, function(leaf_id) {leaf_train_ids[[leaf_id]]})
    },

    ## TODO: Remove wrapper from code
    getTerminalSampleIDs = function(predict_data) {

      leaf_ids <- findLeafIDs(predict_data)
      getTrainIDsByLeaf(leaf_ids)

    },


    setLeafPredictions = function(response_features) {

      if (!is.matrix(response_features)) {
        stop("`response_features` must be a matrix.")
      }

      if (nrow(response_features) != x_data$nrow) {
        stop(
          "`response_features` must have one row for each observation in `x_data`."
        )
      }

      num_nodes <- length(sampleIDs)
      num_features <- ncol(response_features)

      leaf_predictions_new <- matrix(
        0,
        nrow = num_nodes,
        ncol = num_features
      )

      for (nodeID in seq_len(num_nodes)) {

        train_ids <- leaf_train_ids[[nodeID]]

        # Internal nodes do not have leaf predictions
        if (length(train_ids) > 0L) {
          leaf_predictions_new[nodeID, ] <- colMeans(
            response_features[train_ids,,drop = FALSE]
          )
        }
      }

      leaf_predictions <<- leaf_predictions_new

      invisible(.self)
    },


    # Predict the stored values for every observation.
    predictLeafValues = function(predict_data) {

      if (ncol(leaf_predictions) == 0L) {
        stop("Leaf predictions have not been initialized.")
      }

      leaf_ids <- findLeafIDs(predict_data = predict_data)
      leaf_predictions[leaf_ids,,drop = FALSE]

    },


    makeTerminalNode = function(nodeID) {
      ## Save observation indices
      leaf_train_ids[[nodeID]] <<- sampleIDs[[nodeID]]
    }

    )
)
