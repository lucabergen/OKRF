
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

      ordered_chol_features <- chol_features[node_sampleIDs[order_idx],
                                             ,drop = FALSE]

      num_samples <- length(ordered_values)

      ## End positions of groups with equal covariate values
      group_ends <- which(
        ordered_values[-num_samples] < ordered_values[-1L]
      )

      if (length(group_ends) == 0L) {
        return(best_split)
      }

      ## Initially all observations are in the right child
      sum_left <- numeric(approx_rank)
      sum_right <- colSums(ordered_chol_features)

      n_left <- 0L
      n_right <- num_samples
      group_start <- 1L

      for (group_end in group_ends) {

        ## Move one complete value group from right to left
        group_features <- ordered_chol_features[group_start:group_end,,
                                                drop = FALSE]

        group_sum <- colSums(group_features)
        group_size <- group_end - group_start + 1L

        sum_left <- sum_left + group_sum
        sum_right <- sum_right - group_sum

        n_left <- n_left + group_size
        n_right <- n_right - group_size

        ## Respect min_node_size
        if (n_left >= min_node_size && n_right >= min_node_size) {

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

      ## No further split if the node is too small
      if (length(node_sampleIDs) < 2L * min_node_size) {
        return(NULL)
      }

      ## Score of the unsplit node
      node_features <- chol_features[node_sampleIDs,,drop = FALSE]

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

    makeTerminalNode = function(nodeID) {
      # Save observation indices
      terminal_sampleIDs[[nodeID]] <<- sampleIDs[[nodeID]]
    }

    )
)
