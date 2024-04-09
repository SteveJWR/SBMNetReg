#' Degree
#'
#' @param data data frame
#' @param G Adjacency matrix
#'
#' @return degree
#' @export
degree <- function(data, G) {
  d <- Matrix::colSums(G)
  return(d)
}

#' Degree Ratio
#'
#' @param data data frame
#' @param G Adjacency matrix
#'
#' @return degree_ratio
#' @export
degree_ratio <- function(data, G) {

  d <- Matrix::colSums(G)
  mean_d = mean(d)
  return(d/mean_d)
}


#' Number of Treated Neighbors
#'
#' @param data data frame
#' @param G Adjacency matrix
#' @param treatment_column # which column name indicates the treatment
#'
#' @return Number of treated neighbors
#' @export
#'
treated_neighbors <- function(data, G, treatment_column = 'A') {
  return(as.numeric(G %*% data[[treatment_column]]))
}


#' Fraction of Treated Neighbors
#'
#' @param data data frame
#' @param G Adjacency matrix
#' @param treatment_column # which column name indicates the treatment
#'
#' @return fraction of treated neighbors
#' @export
#'
frac_treated_neighbors <- function(data, G, treatment_column = 'A') {

  treated_neighbors <- treated_neighbors(data, G, treatment_column = 'A')
  deg <- degree(data, G)
  frac_treat = treated_neighbors / deg
  frac_treat[is.na(frac_treat)] <- 0 # replaces undefined with 0
  return(frac_treat)
}



#' Power of treated neighbors
#'
#' @param data data frame
#' @param G Adjacency matrix
#' @param treatment_column Which column name indicates the treatment
#' @param power Which power to raise the adjacency matrix to
#'
#' @return number of treated neighbours in the power graph
#' @export
#'
power_treated_neighbors <- function(data, G, treatment_column = 'A', power = 2) {
  G_power = as.matrix(G) %^% power
  power_treated_neighbors <- as.numeric(G_power %*% data[[treatment_column]])
  return(power_treated_neighbors)
}




#' Network Features
#'
#' @param data base data frame
#' @param G Adjacency matrix
#' @param feat_list List of network feature functions
#'
#' @return data frame with network features
#' @export
network_features <- function(data, G, feat_list = list('degree'= degree)) {
  num_F = length(feat_list)
  feature_names = names(feat_list)
  data_tmp = data
  for(i in seq(num_F)){
    feat = feat_list[[i]]
    data_tmp[[feature_names[i]]] = feat(data, G)
  }
  return(data_tmp)
}





# graph features.

# All functions here should satisfy the formatting:
# function_name <- function(data, G, ...):
# Function: Get network features






