#' Sample A Binomial Treatment
#'
#' @param n Number of samples
#' @param p Probability of Treatment Assignment
#'
#' @return Treatment Assignments
#' @export
binomialTreatment <- function(n,p){
  A <- rbinom(n,size = 1, prob = p)
  return(A)
}


#' Cluster Treatment Allocation
#'
#' @param Z Cluster Membership
#' @param p.treat Probability of Fully Treated Clusters
#'
#' @return Treatment Assignments
#' @export
clusterTreatment <- function(Z, p.treat = 1/2){
  unique.clusters <- unique(Z)
  K <- length(unique.clusters)
  K.treat <- round(K*p.treat)
  treatments <- c(rep(1,K.treat),rep(0,K - K.treat))
  A.block <- sample(treatments, replace = F)
  A = rep(NA, length(Z))
  for(i in seq(K)){
    z.true = unique.clusters[i]
    z.idx = which(Z == z.true)
    A[z.idx] = A.block[i]
  }
  return(A)
}


#' Cluster Randomization
#'
#' @param Z Cluster Membership
#' @param levels Saturation levels for each of the clusters assigned in Z
#'
#' @return Treatment Assignments
#' @export
saturationRandomizationTreatment <- function(Z, levels){
  K = max(Z)
  # Z must be an integer between 1 and K
  if(length(levels) != K){
    stop("saturation levels must be the same as the number of clusters")
  }
  n = length(Z)
  A = rep(0,n)
  for(i in seq(K)){
    idx.i = which(Z == i)
    ni = length(idx.i)
    ni.treat = min(round(ni*levels[i]), ni) # prevents sampling more than the number of units in the cluster
    idx.i.treat <- sample(idx.i, ni.treat, replace = F)
    A[idx.i.treat] = 1
  }
  return(A)
}


