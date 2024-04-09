# Utility functions

# TODO: Document and Test these, perhaps some of them can be hidden.

#' Compute Synthetic ARD from a list of discrete traits
#'
#' @param traits Vector or matrix of traits
#' @param G Adjacency Matrix
#'
#' @return ARD responses
#' @export
computeARD <- function (traits, G)
{
  if(is.matrix(traits)){
    n = nrow(G)
    Tr = ncol(traits)
    X = matrix(nrow = n, ncol = Tr)
    for (k in seq(Tr)) {
      ind = 1 * (traits[,k] == 1)
      X[, k] = as.numeric(G %*% ind)
    }
  } else {
    n = nrow(G)
    Tr = max(traits)
    X = matrix(nrow = n, ncol = Tr)
    for (k in seq(Tr)) {
      ind = 1 * (traits == k)
      X[, k] = as.numeric(G %*% ind)
    }
  }
  return(X)
}


#' Computes K-means clustering for ARD responses.
#'
#' @param X ARD responses, should be normalized
#' @param K Number of clusters
#' @param alg Algorithm Used
#'
#' @return K-means clustering
#' @export
clusterARDKmeans <- function(X, K, alg = 'Lloyd'){
  res = stats::kmeans(X, centers = K, algorithm = alg)
  return(res)
}

#' Computes Hierarchical Clustering of the ARD responses
#'
#' @param X ARD responses, should be normalized
#' @param K Number of clusters
#' @param alg Algorithm used For merging
#'
#' @return Hierarchical Clustering
#' @export
clusterARDHclust <- function(X, K, alg = 'complete'){
  D = stats::dist(X)
  tree = hclust(D, method = alg)
  res = cutree(tree, k = K)
  return(res)
}



#' Label swapping to align an estimated labels with a reference set
#'
#' @param Z.ref Reference labels
#' @param Z.hat Estimated labels
#'
#' @return Aligned Estimated labels
#' @export
labelSwitching <- function(Z.ref,Z.hat){
  K = max(Z.hat)
  switch.model <- label.switching::ecr(Z.hat,matrix(Z.ref, nrow = 1), K)
  return(switch.model$permutations[Z.hat])
}


#' Label swapping to align an estimated labels with a reference set
#'
#' @param Z.ref Reference labels
#' @param Z.hat Estimated labels
#'
#' @return permutations
#' @export
labelSwitchingPermutations <- function(Z.ref,Z.hat){
  K = max(Z.hat)
  switch.model <- label.switching::ecr(Z.hat,matrix(Z.ref, nrow = 1), K)
  return(switch.model$permutations)
}



#' Simulate Draw of an SBM
#'
#' @param n number of nodes to simulate
#' @param P Probability matrix
#' @param PI Proportions of each group
#' @param Z Group membership
#'
#' @return Adjacency matrix and group membership
#' @export
#'
generateSBM <- function(n,P,PI,Z){
  if(missing(Z)){
    if(nrow(P) != ncol(P)){
      stop("P must be square")
    }
    if(sum(abs(P  - t(P))) > 0.01){
      stop("P must be symmetric")
    }
    if(length(PI) != nrow(P)){
      stop("proportions and P must be of the same dimension")
    }
    n.groups = round(n*PI/(sum(PI)))
    K = length(PI)
    while(sum(n.groups) < n){
      i = sample(seq(K), size = 1)
      n.groups[i] = n.groups[i] + 1
    }
    while(sum(n.groups) > n){
      i = sample(seq(K), size = 1)
      n.groups[i] = n.groups[i] - 1
    }
    K = length(PI)
    groups = rep(seq(K),times = n.groups)
    g = igraph::sample_sbm(n,pref.matrix = P, block.sizes = n.groups)
    G = igraph::as_adjacency_matrix(g)

  } else {
    if(missing(n)){
      n = length(Z)
    }
    K = max(Z)
    n.groups <- rep(NA,K)
    for(i in seq(K)){
      n.groups[i] = sum(Z == i)
    }

    # permute the ordering at the end so that it will agree with Z
    g = igraph::sample_sbm(n = sum(n.groups),pref.matrix = P, block.sizes = n.groups)
    G = igraph::as_adjacency_matrix(g)

    Z.original = rep(seq(K),times = n.groups)
    perm = Matrix::invPerm(order(Z))
    # marks the end of the groups

    # permute the appropriate Z
    G <- G[perm,perm]
    groups = Z

  }
  return(list('G' = G, "Z" = Z))
}



#' Generate Set Of SBM Draws
#'
#' @param SBM Stochastic Block Model
#' @param L Number of draws
#'
#' @return List of Adjacency Matrices
#' @export
Generate_G_set <- function(SBM, L = 1000){
  if(!(all(c('Z', 'P' ) %in% names(SBM)))){
    L_per_SBM = ceiling(L/length(SBM))
    G_set = list()
    B = length(SBM)
    for(b in seq(B)){
      P = SBM[[b]]$P
      Z = SBM[[b]]$Z
      for(l in seq(L_per_SBM)) {
        g.sim <- generateSBM(n,P,Z = Z)
        G_set[[l + (b-1)*L_per_SBM]] = g.sim$G
      }
    }
  } else {
    P = SBM$P
    Z = SBM$Z
    G_set = list()
    for(l in seq(L)) {
      g.sim <- generateSBM(n,P,Z = Z)
      G_set[[l]] = g.sim$G
    }
  }
  return(G_set)
}





#' Use Leiden Clustering to compute the clusters
#'
#' @param G Adjacency matrix
#' @param K Number of clusters
#' @param adj_mode Only 'undirected' is supported
#' @param max_resolution_parameter Maximum resolution parameter (tuning parameter for finding clusters)
#' @param n_iterations Number of iterations for Leiden algorithm
#'
#' @return Cluster membership
#' @export
#'
cluster_leiden_K <- function(G, K = 4, adj_mode = "undirected", max_resolution_parameter = 10, n_iterations = 5){
  g <- igraph::graph_from_adjacency_matrix(G, mode = adj_mode)
  tmp_cluster = igraph::cluster_leiden(g, objective_function = "modularity",
                               n_iterations = n_iterations, resolution_parameter = max_resolution_parameter)
  resolution_parameter_upper = max_resolution_parameter
  resolution_parameter_lower = 0
  iter = 0
  if(tmp_cluster$nb_clusters == K){
    return(tmp$membership) # we are done
  }
  if(tmp_cluster$nb_clusters <  K) {
    stop("Not enough clusters")
  } else {
    while(tmp_cluster$nb_clusters != K){
      iter = iter + 1
      resolution_parameter_tmp = (resolution_parameter_upper + resolution_parameter_lower)/2
      tmp_cluster = igraph::cluster_leiden(g, objective_function = "modularity",
                                   n_iterations = n_iterations, resolution_parameter = resolution_parameter_tmp)
      if(tmp_cluster$nb_clusters > K){
        resolution_parameter_upper = resolution_parameter_tmp
      } else {
        resolution_parameter_lower = resolution_parameter_tmp
      }
      if( iter > 100){
        stop("Too many iterations")
      }

    }
    return(tmp_cluster$membership)
  }
}



#' Logistic Function
#' @export
logistic <- function(x){
  1/(1 + exp(-x))
}



#' Sample From Probability Matrix
#'
#' @param probabilities
#'
#' @return Sampled clusters
#' @export
sample_clusters_from_prob_matrix <- function(probabilities) {
  draws <- apply(probabilities, 1, function(row_probs) {
    sample.int(length(row_probs), size = 1, prob = row_probs)
  })
  return(draws)
}


#' Importing mclustBIC function
#' @export
mclustBIC <- mclust::mclustBIC
