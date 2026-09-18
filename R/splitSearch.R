
computeFeatureScore <- function(feature_sum, num_observations) {

  sum(feature_sum^2) / num_observations
}


findBestOrderedSplit <- function(ordered_values, ordered_features,
                                 min_leaf_size, split_varID, best_split) {

  num_samples <- length(ordered_values)

  # Compare each value with the next value in the sorted vector.
  # Removing the last value from the left vector and the first value
  # from the right vector creates all adjacent value pairs.
  group_ends <- which(ordered_values[-num_samples] < ordered_values[-1L])

  # No valid split exists when all values are equal. In this case, keep the
  # best split found so far.
  if (length(group_ends) == 0L) {
    return(best_split)
  }

  # Initially all observations are in the right child
  sum_left <- numeric(length = ncol(ordered_features))
  sum_right <- colSums(ordered_features)

  n_left <- 0L
  n_right <- num_samples
  group_start <- 1L

  for (group_end in group_ends) {

    # Move one complete group of equal covariate values from the right child
    # to the left child.
    group_features <- ordered_features[group_start:group_end,,drop = FALSE]

    group_sum <- colSums(group_features)
    group_size <- group_end - group_start + 1L

    sum_left <- sum_left + group_sum
    sum_right <- sum_right - group_sum

    n_left <- n_left + group_size
    n_right <- n_right - group_size

    #  Only admissible leaf node sizes are considered
    if (n_left >= min_leaf_size && n_right >= min_leaf_size) {

      score <- computeFeatureScore(
        feature_sum = sum_left,
        num_observations = n_left
      ) +
        computeFeatureScore(
          feature_sum = sum_right,
          num_observations = n_right
        )

      if (score > best_split$score) {

        # Set the split threshold halfway between the two adjacent
        # distinct covariate values.
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

    # Start the next group at the following observation
    group_start <- group_end + 1L
  }

  best_split
}
