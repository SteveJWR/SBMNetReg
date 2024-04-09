# data generation linear model
# we can generate the features before hand from G

# noise, Normal, input is a dataset,
# formula, coef, G,  data, features, sigma, covar = NULL if

#' Simulate A Linear Model From The Features
#'
#' @param data Input Data
#' @param fmla Formula for the model
#' @param coef Coefficients
#' @param sd Standard Deviation for iid data
#' @param noise_covar Covariance matrix for the noise for non-iid
#'
#' @return List with the outcomes, the features and the noise covariance
#' @export
simLinearModel <- function(data, fmla, coef, sd = 1, noise_covar = NULL){
  # generate the outcomes
  n = nrow(data)
  p = ncol(data)
  X = model.matrix(fmla, data)
  if(is.null(noise_covar)){
    eps = rnorm(n, 0, sd)
  } else {
    eps = MASS::mvrnorm( mu = rep(0,n), Sigma = cov_matrix)
  }

  y = X %*% coef + eps
  return(list('y' = y, 'X' = X, 'noise_covar' = noise_covar))
}


#' Simulate a simple logistic model
#'
#' @param data Input Data
#' @param fmla Formula for the model
#' @param coef Coefficients
#'
#' @return List with the outcomes, the features and the noise covariance
#' @export
simLogisticModel <- function(data, fmla, coef){
  # generate the outcomes
  n = nrow(data)
  p = ncol(data)
  X = model.matrix(fmla, data)
  p = logistic(as.numeric(X %*% coef))
  y = rbinom(n, 1, p)
  return(list('y' = y, 'X' = X))
}





# simHearingModel <- function(){
#
# }

# TODO: Add details
#' @export
simSimpleContagion <- function(A,G,T_steps = 3, q = 0.05, single_infection_time = T){
  infected = A
  infected_last_round = infected
  n = length(A)
  for(step in seq(T_steps)){

    # find the infectious nodes
    if(single_infection_time){
      infectious_neighbors = which(infected_last_round == 1)
    } else {
      infectious_neighbors = which(infected == 1)
    }

    if(length(infectious_neighbors) == 1){
      infect_degree =  G[infectious_neighbors,]
    } else {
      infect_degree =  Matrix::colSums(G[infectious_neighbors,])
    }

    # infection probability based on number of infected neighbours
    infection_probs = 1 - (1 - q)**(infect_degree)

    # infection occurs
    neighbors = which(runif(n) < infection_probs & !infected)

    infected_last_round = numeric(n)
    infected_last_round[neighbors] = 1

    # find the neighbors that are newly infected
    infected[neighbors] = 1
  }
  return(infected)
}


# Similar structure of the Beaman paper
# Complex Contagion

# TODO: Add details
#' @export
simComplexContagion <- function (A, G, thresholds, T_steps = 4)
{
  infected = A
  for (step in seq(T_steps)) {
    treated_neighbors = which(infected == 1)
    if(length(treated_neighbors) > 1 ){
      infected_neighbors = which(Matrix::colSums(G[treated_neighbors,
      ]) > thresholds)
    } else {
      # Handle the singleton case
      infected_neighbors = which(G[treated_neighbors,
      ] > thresholds)
    }
    infected[infected_neighbors] = 1
  }
  return(infected)
}

# Beaman truncated normal adoption threshold
#' @export
BeamanAdoptionThreshold <- function(n, mean = 2, sd = 0.5){
  # generate the adoption threshold
  lower = 0
  return(rtruncnorm(n, a = lower, b = Inf, mean = mean, sd = sd))
}




