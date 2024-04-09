


# Function to estimate the SBM using the subgraph estimator, applies to RDS as well
#' Estimate SBM using subgraph and cluster labels
#'
#' @param Z cluster memberships
#' @param G adjacency matrix
#'
#' @return Stochastic block model (P and Z)
#' @export
SBM_estimate_subgraph <- function(Z,G){
  K = max(unique(Z))
  counts = as.numeric(table(Z))
  P_hat = matrix(0, K, K)
  for(k in seq(K)){
    for(j in seq(K)){
      if(k > j){
        id_k = which(Z == k)
        id_j = which(Z == j)
        P_hat[k,j] = sum(G[id_k,id_j]) / (counts[k] * counts[j])
      }
      else if(k == j){
        id_k = which(Z == k)
        id_j = which(Z == j)
        P_hat[k,j] = sum(G[id_k,id_j]) / (counts[k] * counts[j] - counts[k])
      }
    }
  }
  P_hat = P_hat + t(P_hat)
  diag(P_hat) = diag(P_hat) / 2
  SBM = list()
  SBM$Z = Z
  SBM$P = P_hat
  return(SBM)
}








#' Estimate ARD blocks from trait responses.
#'
#' @param traits Self traits
#' @param X ARD responses
#' @param K Number of clusters
#' @param alg Which algorithm will be used for agglomerative clustering
#' @param Z.ref Reference cluster labels
#' @param constrained_opt Use constrained optimization
#' @param solver Solver to use in CVX if using constrained optimization
#' @param cluster_traits_quantile Quantile for total trait threshold
#'
#' @return Stochastic block model (P and Z)
#' @export
SBM_estimate_ARD <- function (traits, X, K,
                              alg = "ward.D", Z.ref = NULL,
                              constrained_opt = T,solver = "ECOS",
                              cluster_traits_quantile = 0.75)
{
  if (is.matrix(traits)) {
    if (ncol(traits) != ncol(X)) {
      errorCondition("traits vector must align with ARD vector")
    }
  }
  num_T = ncol(X)
  X.norm = X
  nt_vec = numeric(num_T)
  for (t in seq(num_T)) {
    if (is.matrix(traits)) {
      nt = sum(traits[, t] == 1)
    }
    else {
      nt = sum(traits == t)
    }
    X.norm[, t] = X[, t]/nt
    nt_vec[t] = nt
  }
  trait_thresh = as.numeric(quantile(nt_vec, cluster_traits_quantile))
  trait_idx = which(nt_vec >= trait_thresh) # pick a threshold for clustering

  X.norm[is.nan(X.norm)] = 0
  idx_non_zero = which(rowSums(X) > 0)
  idx_zero = which(rowSums(X) == 0)
  if(length(idx_zero) > 0){
    Z.hat.nonzero = clusterARDHclust(X.norm[idx_non_zero,trait_idx], K - 1, alg = alg)
    Z.hat = rep(K, nrow(X))
    Z.hat[idx_non_zero] = Z.hat.nonzero
  } else {
    Z.hat = clusterARDHclust(X.norm[,trait_idx], K, alg = alg)
  }

  if (!is.null(Z.ref)) {
    Z.hat = labelSwitching(Z.ref, Z.hat)
  }
  SBM = list()
  SBM$Z <- Z.hat
  SBM$P <- ARD_SBM_regression(traits, X, Z.hat, constrained_opt = constrained_opt,
                              solver = solver)
  SBM$trait_threshold = trait_thresh
  return(SBM)
}



#' Regression Estimationg for Estimating the SBM given a set of cluster memberships
#'
#' @param traits Traits vector or matrix
#' @param X ARD responses
#' @param Z Cluster memberships
#' @param constrained_opt Use constrained optimization
#' @param solver Solver to use in CVX if using constrained optimization
#'
#' @return Estimated SBM P matrix
#' @export
ARD_SBM_regression <- function(traits,X, Z, constrained_opt = T, solver = 'ECOS'){
  K = max(Z)
  if(is.matrix(traits)){
    if(ncol(traits) != ncol(X)){
      errorCondition('traits vector must align with ARD vector')
    }
  }
  num_T = ncol(X)
  # Probability of connecting to a cluster given a trait
  Pkt = matrix(NA, nrow = K, ncol = num_T)
  # Probability of a trait given a cluster
  Omega = matrix(NA, nrow = K, ncol = num_T)
  # weights for the number of counts for each feature combination
  Weights = matrix(0, nrow = K, ncol = num_T)

  for (k in seq(K)) {
    idx.k = which(Z == k)
    n_k = length(idx.k)
    for (t in seq(num_T)) {
      if(is.matrix(traits)){
        idx.t = which(traits[,t] == 1)
      } else {
        idx.t = which(traits == t)
      }

      n_t = length(idx.t)
      if(is.matrix(traits)){
        n_kt = sum(Z == k & traits[,t] == 1)
      } else {
        n_kt = sum(Z == k & traits == t)
      }

      Omega[k, t] = n_kt/n_t
      X.sub = X[idx.k, t]
      num_counts = sum(X.sub)
      Weights[k, t] = Weights[k, t] + n_k
      Pkt[k, t] = num_counts/(n_k * n_t)
    }
  }

  if(constrained_opt){
    # Define the SBM cross block probabilities
    P_mat <- CVXR::Variable(K,K)
    # Define the objective function (least squares)
    objective <- CVXR::Minimize(sum(((P_mat %*% Omega - Pkt)^2) * Weights))
    # Define the constraints
    constraints <- list(P_mat >= 0,
                        P_mat <= 1,
                        P_mat == t(P_mat))  # Example constraint: coefficients must be non-negative
    # Formulate the problem
    problem <- CVXR::Problem(objective, constraints)
    # Solve the problem
    result <- CVXR::solve(problem, solver = solver)
    # Extract the solution
    P_mat_solution <- result$getValue(P_mat)
    # trim the solution
    P.hat = P_mat_solution
  } else {
    P.hat = solve(Omega %*% t(Omega)) %*% Omega %*%  t(Pkt)
  }
  # Correction For Slight Errors in the Matrix
  P.hat = P.hat + t(P.hat)
  P.hat = P.hat/2
  P.hat[P.hat <= 0] = 0
  P.hat[P.hat >= 1] = 1
  return(P.hat)
}


#' Estimate draws of the ARD model from trait responses.
#'
#' @param traits Self traits
#' @param X ARD responses
#' @param K Number of clusters
#' @param L Number of Bootstrap iterations
#' @param alg Algorithm Used for Agglomerative Clustering
#' @param cluster_traits_quantile Quantile for total trait threshold
#'
#' @return List of Draws from the Stochastic Blockmodel
#' @export
SBM_uncertainty_ARD <- function(traits, X, K,
                                L = 1000, alg = 'ward.D',
                                cluster_traits_quantile = 0.75){
  base_SBM <- SBM_estimate_ARD(traits, X, K,
                               alg = alg, constrained_opt = F,
                               cluster_traits_quantile = cluster_traits_quantile)

  if(is.matrix(traits)){
    if(ncol(traits) != ncol(X)){
      errorCondition('traits vector must align with ARD vector')
    }
  }
  num_T = ncol(X)
  X.norm = X # placeholder
  nt_vec = numeric(num_T)
  for(t in seq(num_T)){
    if(is.matrix(traits)){
      nt = sum(traits[,t] == 1)
    } else {
      nt = sum(traits == t)
    }

    X.norm[,t] = X[,t]/nt
    nt_vec[t] = nt
  }
  trait_thresh = as.numeric(quantile(nt_vec, cluster_traits_quantile))
  trait_idx = which(nt_vec >= trait_thresh)
  X.norm[is.nan(X.norm)] = 0

  gmm_fit <- mclust::Mclust(X.norm[,trait_idx], G = K)
  class.probs <- gmm_fit$z
  n = nrow(X)
  SBM_list <- list()
  for (l in seq(L)) {
    Z.hat.boot = sample_clusters_from_prob_matrix(class.probs)
    P.hat.boot = matrix(NA, nrow = K, ncol = K)
    for (k1 in seq(K)) {
      idx.k1 = which(Z.hat.boot == k1)
      n_k1 = length(idx.k1)
      for (k2 in seq(K)) {
        idx.k2 = which(Z.hat.boot == k2)
        n_k2 = length(idx.k2)
        P.hat.boot[k1, k2] = rbinom(1, n_k1 * n_k2, base_SBM$P[k1,
                                                               k2])/(n_k1 * n_k2)
      }
    }
    P.hat.boot = P.hat.boot + t(P.hat.boot)
    diag(P.hat.boot) = diag(P.hat.boot)/2
    P.hat.boot[P.hat.boot <= 0] = 0
    P.hat.boot[P.hat.boot >= 1] = 1
    SBM_tmp <- list()
    SBM_tmp$Z <- Z.hat.boot
    SBM_tmp$P <- P.hat.boot
    SBM_list[[l]] <- SBM_tmp
  }
  return(SBM_list)
}


#' Estimate draws of the ARD model from trait responses.
#'
#' @param traits Self traits
#' @param X ARD responses
#' @param K Number of clusters
#' @param L Number of Bootstrap iterations
#' @param alg Algorithm Used for Agglomerative Clustering
#' @param frac_observed Correction if there is a fraction of the observed data
#' @param verbose Print progress
#' @param cluster_traits_quantile Quantile for total trait threshold
#'
#' @return Bootstrap draws of the Stochastic Blockmodel
#' @export
SBM_bootstrap_ARD <- function (traits, X, K,
                               L = 1000, alg = "ward.D",
                               frac_observed = 1, verbose = T,
                               cluster_traits_quantile = 0.75)
{
  base_SBM <- SBM_estimate_ARD(traits, X, K,
                               alg = alg, constrained_opt = F,
                               cluster_traits_quantile = cluster_traits_quantile)
  Z.base = base_SBM$Z
  G_list = Generate_G_set(base_SBM, L = L)
  SBM_list <- list()
  for (l in seq(L)) {
    if (verbose & (l%%100 == 0)) {
      cat(paste("Bootstrap ", l, " of ", L), end = "\r")
    }
    X.boot = computeARD(traits, G_list[[l]])
    if (frac_observed < 1) {
      X.boot = X.boot/frac_observed
    }
    SBM_list[[l]] <- tryCatch({
      SBM_estimate_ARD(traits, X.boot, K, alg = alg,
                       Z.ref = Z.base, constrained_opt = F,
                       cluster_traits_quantile = cluster_traits_quantile)
    }, error = function(err) {
      NULL
    })
  }
  SBM_list <- Filter(Negate(is.null), SBM_list)
  return(SBM_list)
}
