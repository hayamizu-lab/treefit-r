# Copyright (C) 2020  Momoko Hayamizu <hayamizu@ism.ac.jp>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

test_that("calculate_distance_matrix", {
  values <- matrix(c(1, 0,
                     0, 1,
                     0, 0),
                   ncol=2,
                   byrow=TRUE)
  distance_matrix <- matrix(c(0,       sqrt(2), 1,
                              sqrt(2), 0,       1,
                              1,       1,       0),
                            ncol=3,
                            byrow=TRUE)
  expect_equivalent(as.matrix(calculate_distance_matrix(values)),
                    distance_matrix)
})

test_that("perturbate_knn", {
  expression <- generate_2d_n_arms_star_data(100, 3, 0.1)
  n_perturbations <- 5
  strength <- 0.2
  min_diff_mean <- 0.001
  max_diff_mean <- 0.02
  min_diff_variance <- 0.00001
  max_diff_variance <- 0.001
  perturbated_expression_list <-
    lapply(1:n_perturbations,
           function(i) {
             perturbate_knn(expression, strength)
           })
  original_distance_matrix <- calculate_distance_matrix(expression)
  perturbated_distance_matrix_list <-
    lapply(perturbated_expression_list,
           calculate_distance_matrix)
  diffs <-
    lapply(perturbated_distance_matrix_list,
           function(perturbated_distance_matrix) {
             abs(original_distance_matrix - perturbated_distance_matrix)
           })
  diff_means <- lapply(diffs, mean)
  diff_variances <- lapply(diffs, var)
  expect_equal(diff_means > min_diff_mean,
               rep(TRUE, n_perturbations))
  expect_equal(diff_means < max_diff_mean,
               rep(TRUE, n_perturbations))
  expect_equal(diff_variances > min_diff_variance,
               rep(TRUE, n_perturbations))
  expect_equal(diff_variances < max_diff_variance,
               rep(TRUE, n_perturbations))
})

test_that("calculate_mst", {
  values <- matrix(c(1, 0,
                     0, 3,
                     0, 0),
                   ncol=2,
                   byrow=TRUE)
  mst <- calculate_mst(values)
  expect_equal(igraph::edge_attr_names(mst),
               character(0))
  expect_equivalent(igraph::as_adjacency_matrix(mst, sparse=FALSE),
                    matrix(c(0, 0, 1,
                             0, 0, 1,
                             1, 1, 0),
                           ncol=3,
                           byrow=TRUE))
})

test_that("treefit: 2 arms", {
  star <- generate_2d_n_arms_star_data(200, 2, 0.1)
  fit <- treefit(list(expression=star),
                 "tree-like")
  expect_equal(list("tree-like",
                    2),
               list(fit$name,
                    fit$n_principal_paths_candidates[1]))
})

test_that("treefit: 3 arms", {
  star <- generate_2d_n_arms_star_data(200, 3, 0.1)
  fit <- treefit(list(expression=star),
                 "tree-like")
  expect_equal(list("tree-like",
                    3),
               list(fit$name,
                    fit$n_principal_paths_candidates[1]))
})
