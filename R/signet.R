##
## SIGNET: Signed Network Spectral Clustering
## Copyright (C) 2026 Zongyue Teng, Dandan Liu and Panpan Zhang
## Zongyue Teng <zongyue.teng@vanderbilt.edu>
##
## This file is part of the R package SIGNET.
##
## The R package SIGNET is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package SIGNET is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

#' @importFrom rARPACK eigs
#' @importFrom stats kmeans
#' @importFrom cluster silhouette
#' @importFrom Rdpack reprompt
NULL

##  signet: Signed Network Spectral Clustering

#' Signed Network Spectral Clustering
#'
#' Function short description
#'
#' @param adj_mat is an adjacency matrix.
#' @param K is the pre-specified number of clusters.
#' @param Kmin is the minimum candidate value for the number of clusters if \code{K}
#' is unknown. The default value is \code{2}.
#' @param Kmax is the maximum candidate value for the number of clusters if \code{K}
#' is unknown.
#' @param method allows the users to select the method for selecting the number of
#' clusters among \code{eigen-gap} (eigen-gap heuristic), \code{sil-max} (maximum
#' silhouette score) and \code{sil-gap} (silhouette-based gap criterion).
#' @param dist_fun allows the users to specify their own distance function used to
#' calculate the silhouette score. The default function is \code{1 - adj_mat}.
#'
#'@return A list containing:
#'\describe{
#'\item{method}{The method used to select the number of clusters.}
#'\item{K_sel}{The selected number of clusters.}
#'\item{cluster_res}{Cluster assignment under \code{K_sel}.}
#'}
#'
#' @references
#' \insertCite{kunegis2010spectral,teng2026signet}{SIGNET}.
#'
#'@example inst/examples/ex-signet.R
#'
#' @keywords internal
#'
#' @export

signet <- function(adj_mat,
                   K = NULL,
                   Kmax = NULL,
                   Kmin = 2,
                   method = c("eigen-gap", "sil-max", "sil-gap"),
                   dist_fun = NULL) {
  # Sanity checks
  if (!isSymmetric(adj_mat)) {
    stop("The adjacency matrix is not symmetric!")
  }
  if (!is.null(K)) {
    if (K != as.integer(K) | K < 2 | K >= dim(adj_mat)[1]) {
      stop("Invalid input of K!")
    }
    if (!is.null(Kmax)) {
      stop("Conflicting specification: Please specify one of K and Kmax only!")
    }
  } else {
    if (is.null(Kmax)) {
      stop("No specification of K or Kmax!")
    } else {
      if (Kmax != as.integer(Kmax) | Kmax < 2 | Kmax >= dim(adj_mat)[2]) {
        stop("Invalid input of Kmax!")
      }
      if (Kmin != as.integer(Kmin) |
          Kmin < 2 | Kmin >= dim(adj_mat)[2]) {
        stop("Invalid input of Kmin!")
      }
      if (Kmax <= Kmin) {
        stop("Kmax must be greater than Kmin!")
      }
    }
    method <- match.arg(method)
  }

  diag(adj_mat) <- 0
  num_nodes <- dim(adj_mat)[1] ## number of nodes

  D_abs <- diag(rowSums(abs(adj_mat)))
  L <- D_abs - adj_mat ## Laplacian
  D_abs_inv_sqrt <- diag(1 / sqrt(diag(D_abs)))
  L_norm <- D_abs_inv_sqrt %*% L %*% D_abs_inv_sqrt ## normalized Laplacian

  L_norm_eigen <- suppressWarnings(rARPACK::eigs(L_norm, num_nodes, which = "LM"))
  L_norm_eigenvalues <- L_norm_eigen$values
  L_norm_eigenvectors <- L_norm_eigen$vectors

  # Known K: Standard spectral clustering
  if (!is.null(K)) {
    if (L_norm_eigenvalues[num_nodes] < 1e-10) {
      selected_eigenvectors <- L_norm_eigenvectors[, (num_nodes - K):(num_nodes - 1)]
    } else {
      selected_eigenvectors <- L_norm_eigenvectors[, (num_nodes - K + 1):num_nodes]
    }
    selected_eigenvectors <- as.matrix(selected_eigenvectors)
    selected_eigenvectors_standard <- selected_eigenvectors / apply(selected_eigenvectors, 1, function(x) {
      sqrt(sum(x^2))
    })

    kmeans_cluster <- stats::kmeans(selected_eigenvectors_standard,
                             centers = K,
                             nstart = num_nodes)
    cluster_res <- kmeans_cluster$cluster
    K_sel <- K
    method <- "standard"
  } else {
    # Unknown K
    K_candidate <- seq(Kmin, Kmax, by = 1)

    # Eigen-gap approach
    if (method == "eigen-gap") {
      ranked_eigenvalues <- sort(L_norm_eigenvalues)
      gap <- diff(ranked_eigenvalues)
      gap_candidate <- gap[K_candidate]
      num_cluster <- K_candidate[which.max(gap_candidate)]

      if (L_norm_eigenvalues[num_nodes] < 1e-10) {
        selected_eigenvectors <- L_norm_eigenvectors[, (num_nodes - num_cluster):(num_nodes - 1)]
      } else {
        selected_eigenvectors <- L_norm_eigenvectors[, (num_nodes - num_cluster + 1):num_nodes]
      }
      selected_eigenvectors <- as.matrix(selected_eigenvectors)
      selected_eigenvectors_standard <- selected_eigenvectors / apply(selected_eigenvectors, 1, function(x) {
        sqrt(sum(x^2))
      })

      kmeans_cluster <- stats::kmeans(selected_eigenvectors_standard,
                               centers = num_cluster,
                               nstart = num_nodes)
      cluster_res <- kmeans_cluster$cluster
      K_sel <- num_cluster

    } else {
      # Silhouette-based approaches
      if (is.null(dist_fun)) {
        dist_mat <- 1 - adj_mat ## adjacency-based distance (default)
      } else {
        if (!is.function(dist_fun)) {
          stop("`dist_fun` must be a function or NULL.")
        }
        dist_mat <- dist_fun(adj_mat)
      }
      if (!is.matrix(dist_mat) || any(dim(dist_mat) != dim(adj_mat))) {
        stop("`dist_fun(adj_mat)` must return a matrix with the same dimensions as `adj_mat`!")
      }
      diag(dist_mat) <- 0

      res <- lapply(K_candidate, function(k) {
        if (L_norm_eigenvalues[num_nodes] < 1e-10) {
          selected_eigenvectors <- L_norm_eigenvectors[, (num_nodes - k):(num_nodes - 1)]
        } else {
          selected_eigenvectors <- L_norm_eigenvectors[, (num_nodes - k + 1):num_nodes]
        }
        selected_eigenvectors <- as.matrix(selected_eigenvectors)
        selected_eigenvectors_standard <- selected_eigenvectors / apply(selected_eigenvectors, 1, function(x) {
          sqrt(sum(x^2))
        })

        kmeans_cluster <- stats::kmeans(selected_eigenvectors_standard,
                                 centers = k,
                                 nstart = num_nodes)
        cluster_res <- kmeans_cluster$cluster

        sil <- cluster::silhouette(cluster_res, dist = dist_mat)

        list(cluster_res = cluster_res,
             sil = mean(sil[, 3]),
             K = k)
      })

      sil <- unlist(lapply(res, function(x) {
        x[["sil"]]
      }))

      if (method == "sil-max") {
        ix <- which.max(sil)
        cluster_res <- res[[ix]]$cluster_res
        K_sel <- res[[ix]]$K
      }

      if (method == "sil-gap") {
        diff <- diff(sil)
        ix_neg <- which(diff < 0)

        if (length(ix_neg) == 0) {
          stop(
            "No negative change in average silhouette score was detected over the candidate K range. Try a greater Kmax or use a different method!"
          )
        }

        ix <- ix_neg[which.min(diff[ix_neg])]
        cluster_res <- res[[ix]]$cluster_res
        K_sel <- res[[ix]]$K
      }
    }
  }

  return(list(
    method = method,
    K = K_sel,
    cluster_res = cluster_res
  ))
}
