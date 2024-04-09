# estimators

# SBM_lm


#' Generates Draws from the SBM and estimates a linear model using the averaged coefficients
#'
#' @param formula Regression formula
#' @param data Dataset
#' @param SBM Model Parameters for the SBM
#' @param network_feature_list List of network features to include in the model
#' @param B Number of simulations
#' @param verbose Print progress
#'
#' @return List with the model and the mean data
#' @export
SBM_lm <- function(formula, data, SBM, network_feature_list = list('degree'= degree, 'frac.treated' =  frac_treated_neighbors), B = 1000, verbose = T){
  # simulate draws of the ARD model

  # Generate Network Draws Of Features
  for(b in seq(B)){
    if(verbose){
      if( b %% 25 == 0 ){
        cat(paste('Simulation', b, '/', B), end = '\r')
      }
    }

    g.sim <- generateSBM(P = SBM$P,Z = SBM$Z)
    data_samp = network_features(data, g.sim$G, feat_list = network_feature_list)
    if(b == 1){
      data_mean = model.matrix(formula, data_samp)/B
    } else {
      data_mean = data_mean + model.matrix(formula, data_samp)/B
    }
  }
  data_mean = as.data.frame(data_mean)
  data_mean$Y = data$Y
  model = lm(formula, data_mean)

  return(list('model' = model, 'mean_data' = data_mean))
}


#' GLM Fit from sample
#'
#' @param formula Regression formula
#' @param data Dataset
#' @param family GLM family
#' @param SBM Model Parameters for the SBM
#' @param network_feature_list List of network features to include in the model
#' @param B Number of Feature Samples
#' @param EM_threshold Tuning threshold for EM algorithm
#' @param max_iter Maximum number of iterations for the EM algorithm
#'
#' @return Glm model for binary data
#' @export
#'
SBM_glm <- function(formula, data, family = quasibinomial(), SBM, network_feature_list = list('degree'= degree, 'frac.treated' =  frac_treated_neighbors), B = 1000, EM_threshold = 1e-6, max_iter =100){
  # Currently this is only guaranteed to work for binary outcomes.

  # simulate draws of the ARD model
  # in this case, we draw the network features

  # Generate Network Draws Of Features
  data_list = list()
  for(b in seq(B)){
    if( b %% 25 == 0 ){
      cat(paste('Simulation', b, '/', B), end = '\r')
    }
    g.sim <- generateSBM(P = SBM$P,Z = SBM$Z)
    data_samp = network_features(data, g.sim$G, feat_list = network_feature_list)
    if(b == 1){
      data_long = model.matrix(formula, data_samp)
      data_list[[b]] = data_long
    } else {
      data_tmp = model.matrix(formula, data_samp)
      data_list[[b]] = data_tmp
    }
  }

  data_long = as.data.frame(do.call("rbind", data_list))

  rownames(data_long) = NULL
  data_long$Y = rep(data$Y, B)
  data_long$sim_id = rep(seq(B), each = nrow(data))
  data_long$obs_id = rep(seq(nrow(data)), times = B)
  model = glm(formula, data_long, family = family)
  old_coefs = model$coefficients


  abs_diff = EM_threshold + 1
  iter = 0
  while ((abs_diff > EM_threshold) & (iter < max_iter)){
    iter = iter + 1
    cat(paste('EM Step: ', iter), end = '\r')
    # E-step
    # compute the probabilities (weights)
    PY1 = as.numeric(predict(model, newdata = data_long, type = 'response'))
    data_long$weights = PY1*(data_long$Y) + (1 - PY1)*(1 - data_long$Y)
    data_long = as.data.frame(data_long %>%
                                dplyr::group_by(obs_id) %>%
                                dplyr::mutate(weights = weights/sum(weights)))
    weights = data_long$weights
    # M-step
    # update the model
    model = glm(formula, data_long, family = family, weights = weights)
    # check the convergence when the coefficients stop changing
    abs_diff = sum(abs(model$coefficients - old_coefs))
    old_coefs = model$coefficients
  }

  # TODO: We can add additional pieces for the variance
  return(list('model' = model, 'long_data' = data_long))
}
# currently only implemented for binary data.



