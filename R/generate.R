#' Generate a 2-dimensional star tree data
#'
#' @description Generate a 2-dimensional star tree data that contain
#'   `n_samples` data points and fit a star tree with `n_arms` arms.
#'
#' @param n_samples The number of samples to be generated.
#"
#' @param n_arms The number of arms to be generated.
#'
#' @param fatness How fat from the based star tree. `[0.0, 1.0]` is
#'   available value range.
#'
#' @return A generated `martix`. The rows and columns correspond to
#'   samples and features.
#'
#' @examples
#' # Generate a 2-dimensional star tree data that contain 500 data points
#' # and fit a star tree with 3 arms. The generated data are a bit noisy but
#' # tree-like.
#' star.tree_like <- treefit::generate_2d_n_arms_star_data(500, 3, 0.1)
#' plot(star.tree_like)
#'
#' # Generate a 2-dimensional star tree data that contain 600 data points
#' # and fit a star tree with 5 arms. The generated data are very noisy and
#' # less tree-like.
#' star.less_tree_like <- treefit::generate_2d_n_arms_star_data(600, 5, 0.9)
#' plot(star.less_tree_like)
#'
#' @export
generate_2d_n_arms_star_data <- function(n_samples, n_arms, fatness) {
  n_features <- 2
  sigma <- fatness / n_arms
  star <- matrix(,
                 nrow=n_samples,
                 ncol=n_features,
                 dimnames=list(lapply(1:n_samples,
                                      function(i) {paste0("sample", i)}),
                               lapply(1:n_features,
                                      function(i) {paste0("feature", i)})))
  for (i in 1:n_samples) {
    arm <- sample(1:n_arms, 1)
    theta <- arm / n_arms * n_features * pi
    position <- c(cos(theta), sin(theta))
    position <- position * stats::runif(1)
    position <- position + stats::rnorm(n_features, sd=sigma)
    star[i, ] <- position
  }
  star
}

#' Generate a 2-dimensional linked star tree data
#'
#' @description Generate a 2-dimensional linked star tree data. Each
#'   star tree data contain `n_samples_vector[i]` data points and fit
#'   a star tree with `n_arms_vector[i]` arms.
#'
#' @param n_samples_vector The vector of the number of samples to be
#'   generated. For example, `c(200, 100, 300)` means that the first
#'   tree has 200 samples, the second tree has 100 samples and the
#'   third tree has 300 samples.
#'
#' @param n_arms_vector The vector of the number of arms to be
#'   generated.  For example, `c(3, 2, 5)` means the first tree fits a
#'   star tree with 3 arms, the second tree fits a star tree with 2
#'   arms and the third tree fits a star tree with 5 arms. The size of
#'   `n_arms_vector` must equal to the size of `n_samples_vector`.
#'
#' @param fatness How fat from the based tree. `[0.0, 1.0]` is
#'   available value range.
#'
#' @return A generated `martix`. The rows and columns correspond to
#'   samples and features.
#'
#' @examples
#' # Generate a 2-dimensional linked star tree data that contain
#' # 200-400-300 data points and fit a linked star tree with 3-5-4
#' # arms. The generated data are a bit noisy but tree-like.
#' linked_star.tree_like <-
#'   treefit::generate_2d_n_arms_linked_star_data(c(200, 400, 300),
#'                                                c(3, 5, 4),
#'                                                0.1)
#' plot(linked_star.tree_like)
#'
#' # Generate a 2-dimensional linked star tree data that contain
#' # 300-200 data points and fit a linked star tree with 4-3 arms.
#' # The generated data are very noisy and less tree-like.
#' linked_star.less_tree_like <-
#'   treefit::generate_2d_n_arms_linked_star_data(c(300, 200),
#'                                                c(4, 3),
#'                                                0.9)
#' plot(linked_star.less_tree_like)
#'
#' @export
generate_2d_n_arms_linked_star_data <- function(n_samples_vector,
                                                n_arms_vector,
                                                fatness) {
  n_features <- 2
  n_total_samples <- sum(n_samples_vector)
  star <- matrix(,
                 nrow=n_total_samples,
                 ncol=n_features,
                 dimnames=list(lapply(1:n_total_samples,
                                      function(i) {paste0("sample", i)}),
                               lapply(1:n_features,
                                      function(i) {paste0("feature", i)})))
  sub_star_offsets <- c(0.0, 0.0)
  for (i in 1:length(n_samples_vector)) {
    n_samples <- n_samples_vector[i]
    n_arms <- n_arms_vector[i]
    sub_star <- generate_2d_n_arms_star_data(n_samples, n_arms, fatness)
    theta <- 2 * pi * (n_arms %/% 2 / n_arms)
    sub_star_offsets[1] <- sub_star_offsets[1] + -cos(theta) + 1
    sub_star_offsets[2] <- sub_star_offsets[2] + -sin(theta)
    sub_star[, 1] <- sub_star[, 1] + sub_star_offsets[1]
    sub_star[, 2] <- sub_star[, 2] + sub_star_offsets[2]
    star <- rbind(star, sub_star)
  }
  star
}
