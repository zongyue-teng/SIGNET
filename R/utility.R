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

##  data_to_adj: Construct adjacency matrix from data
#'
#' Function short description
#'
#' @param data is the data matrix.
#' @param adj_fun allows the users to specify their own function used to
#' construct the adjacency matrix. The default function is Pearson correlation \code{cor(data, use = "p")}.
#'
#' @return The adjacency matrix.
#'
#' @keywords adjacency
#'
#' @export

data_to_adj <- function(data, adj_fun = NULL){

  # Sanity checks
  if (!is.matrix(data) & !is.data.frame(data)) {
    stop("`data` must be a matrix or data frame!")
  }

  if(sum(sapply(data, is.numeric)) != ncol(data)){
    stop("All columns in `data` must be numeric or integer to compute the adjacency matrix!")
  }

  if (is.null(adj_fun)) {
    adj_mat <- cor(data, use = "pairwise.complete.obs")
  } else {
    if (!is.function(adj_fun)) {
      stop("`adj_fun` must be a function or NULL.")
    }
    adj_mat <- adj_fun(data)
  }


  diag(adj_mat) <- 0  ## no self-loop
  return(adj_mat)
}



##  eval_res: Evaluate clustering results
#'
#' Function short description
#'
#' @importFrom mclust adjustedRandIndex
#'
#' @param memberships is a list of membership vectors
#' @param eval_fun allows the users to specify their own function used to
#' evaluate clustering results. The default function is adjusted Rand index \code{mclust::adjustedRandIndex()}.
#'
#' @return A matrix containing the pairwise evaluation results between membership vectors.
#'
#' @references
#' \insertCite{hubert1985comparing}{SIGNET}.
#'
#' @keywords evaluation
#'
#' @export

eval_res <- function(memberships, eval_fun = NULL){

  # Sanity checks
  if (is.matrix(memberships) | is.data.frame(memberships)) {
    memberships <- as.list(as.data.frame(memberships))
  }

  if (!is.list(memberships)) {
    stop("`memberships` must be a list, matrix, or data frame of membership vectors!")
  }

  if(length(unique(sapply(memberships, length))) != 1){
    stop("Membership vectors need to be of the same length!")
  }

  num <- length(memberships)
  if(num < 2){
    stop("At least two membership vectors are required!")
  }

  if (!is.null(eval_fun) & !is.function(eval_fun)) {
    stop("`eval_fun` must be a function or NULL.")
  }

  eval_mat <- matrix(NA, nrow = num, ncol = num)

  for(i in 1:nrow(eval_mat)){
    for(j in i:nrow(eval_mat)){

      if(is.null(eval_fun)){
        eval_mat[j,i] <- eval_mat[i,j] <- mclust::adjustedRandIndex(memberships[[i]], memberships[[j]])
      } else{
        eval_mat[j,i] <- eval_mat[i,j] <- eval_fun(memberships[[i]], memberships[[j]])
      }

    }
  }
  diag(eval_mat) <- 1
  name <- names(memberships)
  if(is.null(name)){
    name <- paste0("V", 1:num)
  }
  colnames(eval_mat) <- name
  rownames(eval_mat) <- name

  return(eval_mat)
}

