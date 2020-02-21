#' Generate a N wands tree that has 2 features
#'
#' @description Generate a N wands tree that has 2 features.
#'
#' @param n_samples The number of samples to be generated.
#"
#' @param n_wands The number of wands to be generated.
#'
#' @param fatness How much fat from the based tree. `[0.0, 1.0]` is
#'   available value range.
#'
#' @return A generated `martix`. Rows are samples and columns are
#'   features.
#'
#' @examples
#' # Generate a 3 wands tree data that has 500 samples including a bit noise
#' # samples.
#' tree_thin <- treefit::generate_n_wands_2d_tree_expression(500, 3, 0.1)
#' plot(tree_thin)
#'
#' # Generate a 5 wands tree data that has 600 samples including many
#' # noise samples.
#' tree_fat <- treefit::generate_n_wands_2d_tree_expression(600, 5, 0.9)
#' plot(tree_fat)
#'
#' @export
generate_n_wands_2d_tree_expression <- function(n_samples, n_wands, fatness) {
  n_features <- 2
  sigma <- fatness / n_wands
  wand_lengths <- stats::rnorm(n_wands)
  tree <- matrix(,
                 nrow=n_samples,
                 ncol=n_features,
                 dimnames=list(lapply(1:n_samples,
                                      function(i) {paste0("sample", i)}),
                               lapply(1:n_features,
                                      function(i) {paste0("feature", i)})))
  for (i in 1:n_samples) {
    wand <- sample(1:n_wands, 1)
    theta <- wand / n_wands * n_features * pi
    position <- c(cos(theta), sin(theta))
    position <- position * stats::runif(1)
    position <- position + stats::rnorm(n_features, sd=sigma)
    tree[i, ] <- position
  }
  tree
}

#' Generate a tree that consists of multiple N wands tree
#'
#' @description Generate a tree that consists of multiple N wands
#'   tree. The generated tree has 2 features.
#'
#' @param n_samples_vector The vector of the number of samples to be
#'   generated. For example, `c(200, 100, 300)` means that the first
#'   tree has 200 samples, the second tree has 100 samples and the
#'   third tree has 300 samples.
#'
#' @param n_wands_vector The vector of the number of wands to be
#'   generated.  For example, `c(3, 2, 5)` means the first tree has 3
#'   wands, the second tree has 2 wands and the third tree has 5
#'   wands. The size of `n_wands_vector` must equal to the size of
#'   `n_samples_vector`.
#'
#' @param fatness How much fat from the based tree. `[0.0, 1.0]` is
#'   available value range.
#'
#' @return A generated `martix`. Rows are samples and columns are
#'   features.
#'
#' @examples
#' # Generate a 3-5-4 wands linked tree data that have
#' # 200-400-300 samples including a bit noise samples.
#' linked_tree_thin <-
#'   treefit::generate_n_wands_linked_2d_tree_expression(c(200, 400, 300),
#'                                                       c(3, 5, 4),
#'                                                       0.1)
#' plot(linked_tree_thin)
#'
#' # Generate a 4-3 wands linked tree data that have
#' # 300-200 samples including many noise samples.
#' linked_tree_fat <-
#'   treefit::generate_n_wands_linked_2d_tree_expression(c(300, 200),
#'                                                       c(4, 3),
#'                                                       0.9)
#' plot(linked_tree_fat)
#'
#' @export
generate_n_wands_linked_2d_tree_expression <- function(n_samples_vector,
                                                       n_wands_vector,
                                                       fatness) {
  n_features <- 2
  n_total_samples <- sum(n_samples_vector)
  tree <- matrix(,
                 nrow=n_total_samples,
                 ncol=n_features,
                 dimnames=list(lapply(1:n_total_samples,
                                      function(i) {paste0("sample", i)}),
                               lapply(1:n_features,
                                      function(i) {paste0("feature", i)})))
  sub_tree_offsets <- c(0.0, 0.0)
  for (i in 1:length(n_samples_vector)) {
    n_samples <- n_samples_vector[i]
    n_wands <- n_wands_vector[i]
    sub_tree <- generate_n_wands_2d_tree_expression(n_samples, n_wands, fatness)
    theta <- 2 * pi * (n_wands %/% 2 / n_wands)
    sub_tree_offsets[1] <- sub_tree_offsets[1] + -cos(theta) + 1
    sub_tree_offsets[2] <- sub_tree_offsets[2] + -sin(theta)
    sub_tree[, 1] <- sub_tree[, 1] + sub_tree_offsets[1]
    sub_tree[, 2] <- sub_tree[, 2] + sub_tree_offsets[2]
    tree <- rbind(tree, sub_tree)
  }
  tree
}
