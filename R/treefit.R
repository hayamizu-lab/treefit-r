#' Generate perturbated counts by the Poisson distribution
#'
#' @description Generate perturbated counts from the original counts
#'   by the Poisson distribution.
#'
#' @note This is an API for advanced users. This API may be changed.
#'
#' @param counts The original counts. The rows and columns correspond
#'   to samples and features. The values are count of features.
#'
#' @param strength How much perturbated. `0.0` is weak. `1.0` is strong.
#'
#' @export
perturbate_poisson <- function(counts, strength=1.0) {
  matrix(stats::rpois(length(counts), round(counts * strength)),
         nrow=nrow(counts),
         ncol=ncol(counts),
         dimnames=dimnames(counts))
}

calculate_mst <- function(expression) {
  distance_matrix <- calculate_distance_matrix(expression)
  graph <- igraph::graph_from_adjacency_matrix(as.matrix(distance_matrix),
                                               weighted=TRUE,
                                               mode="undirected")
  weighted_mst <- igraph::mst(graph)
  igraph::delete_edge_attr(weighted_mst, "weight")
}

calculate_distance_matrix <- function(expression) {
  # Use "euclidean" method for now.
  # TODO: make method parameter
  stats::dist(expression, method="euclidean")
}

#' Generate perturbated expression by k-NN data
#'
#' @description Generate perturbated expression from the original
#'   expression based on k-NN (k-nearest neighbor) data.
#'
#' @note This is an API for advanced users. This API may be changed.
#'
#' @param expression The original expression. The rows and columns
#'   correspond to samples and features. The expression is normalized
#'   count of features.
#'
#' @param strength How much perturbated. `0.0` is weak. `1.0` is strong.
#'
#' @export
perturbate_knn <- function(expression, strength=1.0) {
  n_samples <- nrow(expression)
  n_features <- ncol(expression)
  # TODO: Improve
  k_nearest_neighbors <- round(n_samples * 0.0125)
  if (k_nearest_neighbors < 2) {
    k_nearest_neighbors <- 2
  }
  standard_deviation <- strength / sqrt(k_nearest_neighbors)
  distance_matrix <- as.matrix(calculate_distance_matrix(expression))
  perturbated_expression <- matrix(nrow=n_samples,
                                   ncol=n_features,
                                   dimnames=dimnames(expression))
  for (i in 1:n_samples) {
    sorted_indices <- order(distance_matrix[i, ])
    nearest_neighbors <-
      expression[sorted_indices[1:(k_nearest_neighbors + 1)], ]
    diffs <- apply(nearest_neighbors,
                   1,
                   function(row) {row - expression[i, ]})
    diffs <- t(diffs)
    weights <- stats::rnorm(n_features, sd=standard_deviation)
    weighted_diffs <- apply(diffs, 1, function(diff) {
      diff * weights
    })
    weighted_diffs <- t(weighted_diffs)
    perturbated_expression[i, ] <- expression[i, ] + colSums(weighted_diffs)
  }
  perturbated_expression
}

calculate_low_dimension_laplacian_eigenvectors <- function(mst, p) {
  e <- eigen(igraph::laplacian_matrix(mst))
  n_target_vectors <- nrow(e$vectors)
  ## Remove zero eigenvalues
  while (n_target_vectors > 0 && all.equal(e$values[n_target_vectors], 0) == TRUE) {
    n_target_vectors <- n_target_vectors - 1
  }
  low_dimension_values <-
    e$values[n_target_vectors:(n_target_vectors - p + 1)]
  low_dimension_vectors <-
    e$vectors[, n_target_vectors:(n_target_vectors - p + 1)]
  if (anyDuplicated(low_dimension_values)) {
    message(paste("low dimension eigenvalues have no duplicated values: ",
                  toString(low_dimension_values),
                  "\n",
                  sep=""))
    # TODO: We need to replace only eigenvalues that has the same eigenvalue.
    low_dimension_vectors = pracma::gramSchmidt(low_dimension_vectors)$Q
  }
  apply(low_dimension_vectors,
        2,
        function(vector) {vector / norm(vector, type="2")})
}

calculate_canonical_correlation <- function(u, v) {
  uTv <- t(u) %*% v
  svd(uTv)$d
}

calculate_grassmann_distance_max_cca <- function(canonical_correlation) {
  max_cos_theta <- sort(canonical_correlation, decreasing=TRUE)[1]
  sqrt(max(0, 1 - max_cos_theta ** 2))
}

calculate_grassmann_distance_rms_cca <- function(canonical_correlation) {
  n_features <- length(canonical_correlation)
  (n_features - sum(canonical_correlation ** 2)) / n_features
}

is_seurat <- function(object) {
  inherits(object, "Seurat")
}

do_normalize <- function(original, target, normalize, verbose) {
  if (is.null(normalize)) {
    if (is_seurat(original)) {
      assay_name <- target
      original <- Seurat::NormalizeData(original,
                                        assay_name,
                                        verbose=verbose)
      original <- Seurat::FindVariableFeatures(original,
                                               assay_name,
                                               verbose=verbose)
      original <- Seurat::ScaleData(original,
                                    assay=assay_name,
                                    verbose=verbose)
      original
    } else {
      # TODO: Default normalization
      target
    }
  } else if (normalize) {
    if (is_seurat(original)) {
      Seurat::DefaultAssay(original) <- target
      normalize(original)
    } else {
      normalize(target)
    }
  } else {
    target
  }
}

seurat_pca_reduction_name <- function(base_assay_name) {
  paste0(base_assay_name, "_pca")
}

do_reduce_dimension <- function(original, target, reduce_dimension, verbose) {
  if (is.null(reduce_dimension) || is.numeric(reduce_dimension)) {
    if (is.numeric(reduce_dimension)) {
      n_dimensions <- reduce_dimension
    } else {
      n_dimensions <- NULL
    }
    if (is_seurat(original)) {
      assay_name <- target
      pca_reduction_name <- seurat_pca_reduction_name(assay_name)
      if (is.null(n_dimensions)) {
        n_dimensions <- 50 # The default value of Seurat::RunPCA
      }
      original <- Seurat::RunPCA(original,
                                 assay_name,
                                 npcs=n_dimensions,
                                 reduction.name=pca_reduction_name,
                                 reduction.key="TREEFITPC_",
                                 verbose=verbose)
      original
    } else {
      target <- stats::prcomp(target, rank.=n_dimensions)$x
    }
  } else if (reduce_dimension) {
    if (is_seurat(original)) {
      Seurat::DefaultAssay(original) <- target
      reduce_dimension(original)
    } else {
      reduce_dimension(target)
    }
  } else {
    target
  }
}

do_build_tree <- function(original, target, build_tree, verbose) {
  if (is.null(build_tree)) {
    if (is_seurat(original)) {
      assay_name <- target
      pca_reduction_name <- seurat_pca_reduction_name(assay_name)
      if (pca_reduction_name %in% names(original@reductions)) {
        if (verbose) {
          message("use PCA")
        }
        data <- Seurat::Embeddings(original, pca_reduction_name)
      } else {
        if (verbose) {
          message("use scaled data")
        }
        data <- Seurat::GetAssayData(original, "scale.data", assay_name)
        data <- t(data)
      }
      calculate_mst(data)
    } else {
      calculate_mst(target)
    }
  } else {
    if (is_seurat(original)) {
      Seurat::DefaultAssay(original) <- target
      build_tree(original)
    } else {
      build_tree(target)
    }
  }
}

calculate_eigenvectors_list <- function(original,
                                        perturbations,
                                        normalize,
                                        reduce_dimension,
                                        build_tree,
                                        max_p,
                                        verbose,
                                        n_perturbations) {
  poisson_strength <- 1.0
  # TODO: Improve
  knn_strength <- 0.2 * (500 / 200) ** 0.5
  targets <- c()

  if ("poisson" %in% perturbations || is.null(perturbations)) {
    if (verbose) {
      message("Poisson perturbation")
    }
    if (is_seurat(original)) {
      default_assay_name <- Seurat::DefaultAssay(original)
      counts <- t(as.matrix(Seurat::GetAssayData(original, "counts")))
      for (i in 1:n_perturbations) {
        assay_name <- paste0(default_assay_name, "Perturbed", i)
        perturbated_counts <- perturbate_poisson(counts, poisson_strength)
        original[[assay_name]] <-
          Seurat::CreateAssayObject(counts=t(perturbated_counts))
        targets <- c(targets, list(assay_name))
      }
    } else if (!is.null(original$counts)) {
      for (i in 1:n_perturbations) {
        perturbated_counts <- perturbate_poisson(original$counts,
                                                 poisson_strength)
        targets <- c(targets, list(perturbated_counts))
      }
    } else if (!is.null(perturbations)) {
      stop(paste("no count data:", original))
    }

    if (!is.null(targets)) {
      if (verbose) {
        message("Normalization")
      }
      if (is_seurat(original)) {
        for (assay_name in targets) {
          original <- do_normalize(original, assay_name, normalize, verbose)
        }
      } else {
        targets <- lapply(targets,
                          function(target) {
                            do_normalize(original, target, normalize, verbose)
                          })
      }
    }
  }

  if ("knn" %in% perturbations || is.null(perturbations)) {
    if (verbose) {
      message("k-NN perturbation")
    }
    if (is.null(targets)) {
      if (is.null(original$expression)) {
        stop(paste("no expression data:", original))
      }

      for (i in 1:n_perturbations) {
        perturbated_expression <-
          perturbate_knn(original$expression, knn_strength)
        targets <- c(targets, list(perturbated_expression))
      }
    } else {
      if (is_seurat(original)) {
        for (assay_name in targets) {
          expression <- Seurat::GetAssayData(original, "data", assay_name)
          expression <- t(as.matrix(expression))
          original <-
            Seurat::SetAssayData(original,
                                 "data",
                                 t(perturbate_knn(expression, knn_strength)),
                                 assay_name)
        }
      } else {
        targets <- lapply(targets,
                          function(expression) {
                            perturbate_knn(expression, knn_strength)
                          })
      }
    }
  }

  lapply(targets,
         function(target) {
           if (verbose) {
             message("Dimensionality reduction")
           }
           if (is_seurat(original)) {
             original <- do_reduce_dimension(original,
                                             target,
                                             reduce_dimension,
                                             verbose)
           } else {
             target <- do_reduce_dimension(original,
                                           target,
                                           reduce_dimension,
                                           verbose)
           }
           if (verbose) {
             message("Build tree")
           }
           tree <- do_build_tree(original, target, build_tree, verbose)
           if (verbose) {
             message("Calculate laplacian eigenvectors")
           }
           calculate_low_dimension_laplacian_eigenvectors(tree, max_p)
         })
}

#' Estimate the goodness-of-fit between tree models and data
#'
#' @description Estimate the goodness-of-fit between tree models and
#'   data.
#'
#' @param target The target data to be estimated. It must be one of them:
#'
#'   * `list(counts=COUNTS, expression=EXPRESSION)`: You must specify
#'     at least one of `COUNTS` and `EXPRESSION`.  They are
#'     `matrix`. The rows and columns correspond to samples such
#'     cells and features such as genes. `COUNTS`'s value is count
#'     data such as the number of genes expressed.  `EXPRESSION`'s
#'     value is normalized count data.
#'   * `Seurat` object
#'
#' @param name The name of `target` as string.
#'
#' @param perturbations How to perturbate the target data.
#'
#'   If this is `NULL`, all available perturbation methods are used.
#'
#'   You can specify used perturbation methods as `list`. Here are
#'   available methods:
#'
#"     * `"poisson"`: A perturbation method for counts data.
#"     * `"knn"`: A perturbation method for expression data.
#'
#' @param normalize How to normalize counts data.
#'
#'   If this is `NULL`, the default normalization is applied.
#'
#'   You can specify a function that normalizes counts data.
#'
#' @param reduce_dimension How to reduce dimension of expression data.
#'
#'   If this is `NULL`, the default dimensionality reduction is applied.
#'
#'   You can specify a function that reduces dimension of expression data.
#'
#' @param build_tree How to build a tree of expression data.
#'
#'   If this is `NULL`, MST is built.
#'
#'   You can specify a function that builds tree of expression data.
#'
#' @param max_p How many low dimension Laplacian eigenvectors are used.
#'
#'   The default is 20.
#'
#' @param verbose Show messages if `TRUE`, don't show any messages otherwise.
#'
#'   The default is `FALSE`
#'
#' @param n_perturbations How many times to perturb.
#'
#'   The default is 20.
#'
#' @return An estimated result as a `treefit` object. It has the
#'   following elements:
#'
#'   * `max_cca_distance`: The result of max canonical correlation
#'     analysis distance as `data.frame`.
#'   * `rms_cca_distance`: The result of root mean square canonical
#'     correlation analysis distance as `data.frame`.
#'   * `n_principal_paths_candidates`: The candidates of the number of
#'     principal paths.
#'
#'   `data.frame` of `max_cca_distance` and `rms_cca_distance` has the
#'   same structure. They have the following columns:
#'
#'   * `p`: Dimensionality of the feature space of tree structures.
#'   * `mean`: The mean of the target distance values.
#'   * `standard_deviation`: The standard deviation of the target
#'     distance values.
#'
#' @examples
#' # Generate a star tree data that have normalized expression values
#' # not count data.
#' star <- treefit::generate_2d_n_arms_star_data(300, 3, 0.1)
#' # Estimate tree-likeness of the tree data.
#' fit <- treefit::treefit(list(expression=star))
#'
#' @export
treefit <- function(target,
                    name=NULL,
                    perturbations=NULL,
                    normalize=NULL,
                    reduce_dimension=NULL,
                    build_tree=NULL,
                    max_p=20,
                    verbose=FALSE,
                    n_perturbations=20) {
  eigenvectors_list <- calculate_eigenvectors_list(target,
                                                   perturbations,
                                                   normalize,
                                                   reduce_dimension,
                                                   build_tree,
                                                   max_p,
                                                   verbose,
                                                   n_perturbations)
  if (is.null(name)) {
    name <- deparse(substitute(target))
  }
  ps <- c()
  max_cca_distance.means <- c()
  max_cca_distance.standard_deviations <- c()
  rms_cca_distance.means <- c()
  rms_cca_distance.standard_deviations <- c()
  for (p in 1:max_p) {
    eigenvectors_pairs <- utils::combn(1:length(eigenvectors_list), 2)
    n_eigenvectors_pairs <- ncol(eigenvectors_pairs)
    ps <- c(ps, p)
    max_cca_distance.values <- c()
    rms_cca_distance.values <- c()
    for (i in 1:n_eigenvectors_pairs) {
      u <- eigenvectors_list[[eigenvectors_pairs[1, i]]][, 1:p]
      v <- eigenvectors_list[[eigenvectors_pairs[2, i]]][, 1:p]
      canonical_correlation <- calculate_canonical_correlation(u, v)
      max_cca_distance.values <-
        c(max_cca_distance.values,
          calculate_grassmann_distance_max_cca(canonical_correlation))
      rms_cca_distance.values <-
        c(rms_cca_distance.values,
          calculate_grassmann_distance_rms_cca(canonical_correlation))
    }
    max_cca_distance.means <-
      c(max_cca_distance.means,
        mean(max_cca_distance.values))
    max_cca_distance.standard_deviations <-
      c(max_cca_distance.standard_deviations,
        stats::sd(max_cca_distance.values))
    rms_cca_distance.means <-
      c(rms_cca_distance.means,
        mean(rms_cca_distance.values))
    rms_cca_distance.standard_deviations <-
      c(rms_cca_distance.standard_deviations,
        stats::sd(rms_cca_distance.values))
  }
  n_principal_paths_candidates <- c()
  for (p in 1:(max_p - 1)) {
    if (p == 1) {
      rms_cca_distance_mean_before <- Inf
    } else {
      rms_cca_distance_mean_before <- rms_cca_distance.means[p - 1]
    }
    rms_cca_distance_mean <- rms_cca_distance.means[p]
    rms_cca_distance_mean_after <- rms_cca_distance.means[p + 1]
    if (rms_cca_distance_mean_before > rms_cca_distance_mean &&
          rms_cca_distance_mean < rms_cca_distance_mean_after) {
      n_principal_paths_candidates <-
        c(n_principal_paths_candidates,
          p + 1)
    }
  }

  attributes <- list(
    name=name,
    max_cca_distance=
      data.frame(p=ps,
                 mean=max_cca_distance.means,
                 standard_deviation=max_cca_distance.standard_deviations),
    rms_cca_distance=
      data.frame(p=ps,
                 mean=rms_cca_distance.means,
                 standard_deviation=rms_cca_distance.standard_deviations),
    n_principal_paths_candidates=n_principal_paths_candidates
  )
  structure(attributes, class=c("treefit"))
}

#' Plot estimated results
#'
#' @description Plot estimate results to get insight.
#'
#' @param x The estimated result by `treefit()` to be visualized.
#'
#' @param ... The more estimated results to be visualized together or
#'   other graphical parameters.
#'
#' @method plot treefit
#'
#' @examples
#' # Generate a tree data.
#' tree <- treefit::generate_2d_n_arms_star_data(200, 3, 0.1)
#' # Estimate the goodness-of-fit between tree models and the tree data.
#' fit <- treefit::treefit(list(expression=tree), "tree")
#' # Visualize the estimated result.
#' plot(fit)
#'
#' # You can mix multiple estimated results by adding "name" column.
#' tree2 <- treefit::generate_2d_n_arms_star_data(200, 3, 0.9)
#' fit2 <- treefit::treefit(list(expression=tree2), "tree2")
#' plot(fit, fit2)
#'
#' @import patchwork
#'
#' @export
plot.treefit <- function(x, ...) {
  fits <- list(x, ...)
  if (length(fits) == 1) {
    all_max_cca_distances <- fits[[1]]$max_cca_distance
    all_rms_cca_distances <- fits[[1]]$rms_cca_distance
  } else {
    all_max_cca_distances <- cbind(fits[[1]]$max_cca_distance,
                                   name=fits[[1]]$name)
    all_rms_cca_distances <- cbind(fits[[1]]$rms_cca_distance,
                                   name=fits[[1]]$name)
    for (fit in fits[2:length(fits)]) {
      all_max_cca_distances <- rbind(all_max_cca_distances,
                                     cbind(fit$max_cca_distance,
                                           name=fit$name))
      all_rms_cca_distances <- rbind(all_rms_cca_distances,
                                     cbind(fit$rms_cca_distance,
                                           name=fit$name))
    }
  }
  plot_data_frame <- function(title, value_label, data_frame) {
    if ("name" %in% names(data_frame)) {
      group <- "name"
      color <- "name"
    } else {
      group <- NULL
      color <- NULL
    }
    ggplot2::ggplot(data_frame) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5),
                     legend.position="top",
                     legend.title=ggplot2::element_blank()) +
      ggplot2::geom_line(
        ggplot2::aes_string("p",
                            "mean",
                            group=group,
                            color=color)) +
      ggplot2::geom_pointrange(
        ggplot2::aes_string("p",
                            "mean",
                            group=group,
                            color=color,
                            ymin="mean - standard_deviation",
                            ymax="mean + standard_deviation"),
        shape=1) +
      ggplot2::scale_x_continuous(breaks=seq(1, max(data_frame$p))) +
      # ggplot2::coord_cartesian(ylim=c(0, 2)) +
      ggplot2::labs(x="p: Dimensionality of the feature space of trees",
                    y=paste(value_label, "(mean and SD)"))
  }
  plot_data_frame("Analysis of the structural instability\nof the estimated trees",
                  "max_cca_distance",
                  all_max_cca_distances) |
    plot_data_frame("Prediction for\nthe number of principal paths",
                    "rms_cca_distance",
                    all_rms_cca_distances)
}
