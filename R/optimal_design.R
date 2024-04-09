

# Defining a function to calculate the variance of the linear combination of the regression features
#' Variance Function Generator
#'
#' @param A Treatment Assignment
#' @param phi Coefficients for the linear combination
#' @param formula Formula for the regression model
#' @param data Observed data frame
#' @param network_feature_list List of network features to include in the model
#' @param G_set Set of samples of the graph
#' @param Sigma Correlation matrix of the additive noise
#' @param verbose Print progress
#'
#' @return Variance of the linear combination of coefficients
#' @export
linear_variance_function <- function(A, phi, formula, data, network_feature_list, G_set, Sigma, verbose = T){
  if(length(G_set) <= 1){
    stop("G_set must have more than one sample")
  }
  # default setting of the name to the treatment assignment
  data_copy = data
  data_copy$A = A

  B = length(G_set)
  for(b in seq(B)){
    if(( b %% 25 == 0 )& verbose){
      cat(paste('Simulation', b, '/', B), end = '\r')
    }
    data_samp = network_features(data_copy, as.matrix(G_set[[b]]),
                                 feat_list = network_feature_list)
    if(b == 1){
      data_mean = model.matrix(formula, data_samp)/B
    } else {
      data_mean = data_mean + model.matrix(formula, data_samp)/B
    }
  }

  F_mat = as.matrix(data_mean)
  Sig_F = (t(F_mat) %*% F_mat)

  # avoiding errors from non-invertible matrices.
  if(min(eigen(Sig_F)$values) < 10**(-7)){
    V_a = 10^7
  } else {
    Sig_F_inv = solve(Sig_F)
    V_a = as.numeric(t(phi) %*% Sig_F_inv %*% t(F_mat) %*% Sigma %*% F_mat%*% Sig_F_inv %*% phi)
  }
  V_a = min(c(V_a, 10^7))
  return(V_a)
}

#' Function Factory for Bayesian Optimization of the linear variance function
#'
#' @param phi Coefficients for the linear combination
#' @param formula Formula for the regression model
#' @param data Observed data frame
#' @param network_feature_list List of network features to include in the model
#' @param G_set Set of samples of the graph
#' @param Sigma Correlation matrix of the additive noise
#'
#' @return Variance where one can compute this using a function of group saturations tau
#' @export
linear_variance_function_factory <- function(phi, formula, data, network_feature_list, G_set, clusters, Sigma){
  # force the input parameters so the Bayesian Optimization can be called
  force(phi)
  force(G_set)
  force(data)
  force(clusters)
  force(network_feature_list)
  force(Sigma)

  var_func <- function(tau){
    A = saturationRandomizationTreatment(clusters, tau)
    return( linear_variance_function(A,phi, formula, data, network_feature_list, G_set, Sigma, verbose = F))
  }
  return(var_func)
}





